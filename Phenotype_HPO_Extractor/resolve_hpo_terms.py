#!/usr/bin/env python3
"""Resolve HPO terms from either a disease name or a phenotype paragraph.

usage: resolve_hpo_terms.py [-h] --input INPUT [--outdir OUTDIR]


Mode
----
The script auto-detects whether the input looks like a disease label
or a phenotype paragraph.

Data sources
------------
- data/hpo/phenotype.hpoa
- data/hpo/hp.obo
- data/orphanet/en_product4.xml
- data/mondo/mondo.obo

Outputs
-------
- TSV with matched rows
- JSON summary with resolution metadata
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from dataclasses import dataclass, field
from difflib import SequenceMatcher
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple


ROOT = Path(__file__).resolve().parent
HPOA_PATH = ROOT / "data" / "hpo" / "phenotype.hpoa"
HPO_OBO_PATH = ROOT / "data" / "hpo" / "hp.obo"
ORPHA_XML_PATH = ROOT / "data" / "orphanet" / "en_product4.xml"
MONDO_OBO_PATH = ROOT / "data" / "mondo" / "mondo.obo"


STOPWORDS = {
    "and",
    "are",
    "as",
    "an",
    "be",
    "by",
    "for",
    "from",
    "in",
    "is",
    "of",
    "or",
    "on",
    "presented",
    "patient",
    "patients",
    "the",
    "to",
    "with",
    "without",
    "disease",
    "syndrome",
    "disorder",
    "condition",
    "type",
    "forms",
    "form",
    "familial",
    "congenital",
    "acquired",
}

GENERIC_SINGLE_WORD_BLACKLIST = {
    "affected",
    "anterior",
    "bilateral",
    "distal",
    "focal",
    "generalized",
    "left",
    "lower",
    "mild",
    "moderate",
    "proximal",
    "right",
    "severe",
    "upper",
}


def normalize_text(text: str) -> str:
    text = text.lower()
    text = text.replace("-", " ")
    text = re.sub(r"[^a-z0-9:]+", " ", text)
    tokens = [tok for tok in text.split() if tok and tok not in STOPWORDS]
    return " ".join(tokens)


def tokenize(text: str) -> List[str]:
    norm = normalize_text(text)
    return norm.split() if norm else []


def looks_like_id(text: str) -> bool:
    return bool(re.search(r"\b(?:ORPHA:\d+|OMIM:\d+|MONDO:\d+|HP:\d+)\b", text, re.I))


def read_single_value(path: Path) -> str:
    content = path.read_text(encoding="utf-8", errors="ignore").strip()
    return content


def collect_input_chunks(text: str) -> List[str]:
    """Return candidate chunks from a possibly multi-line disease paragraph."""
    raw_chunks: List[str] = [text.strip()]
    raw_chunks.extend(line.strip() for line in text.splitlines() if line.strip())
    raw_chunks.extend(
        chunk.strip()
        for chunk in re.split(r"[.;\n]+", text)
        if chunk.strip()
    )

    # Keep ordering, drop duplicates.
    seen: Set[str] = set()
    chunks: List[str] = []
    for chunk in raw_chunks:
        norm = re.sub(r"\s+", " ", chunk).strip()
        if not norm or norm in seen:
            continue
        seen.add(norm)
        chunks.append(norm)
    return chunks


@dataclass
class HpoTerm:
    hpo_id: str
    name: str
    synonyms: List[str] = field(default_factory=list)


@dataclass
class DiseaseRecord:
    disease_id: str
    name: str
    source: str
    synonyms: Set[str] = field(default_factory=set)
    hpo_ids: Set[str] = field(default_factory=set)


def parse_hpoa(path: Path) -> Tuple[Dict[str, DiseaseRecord], Dict[str, Set[str]]]:
    diseases: Dict[str, DiseaseRecord] = {}
    disease_to_hpos: Dict[str, Set[str]] = defaultdict(set)
    with path.open(encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            disease_id = cols[0].strip()
            disease_name = cols[1].strip()
            hpo_id = cols[3].strip()
            disease_to_hpos[disease_id].add(hpo_id)
            if disease_id not in diseases:
                diseases[disease_id] = DiseaseRecord(disease_id=disease_id, name=disease_name, source="HPOA")
    for disease_id, hpos in disease_to_hpos.items():
        diseases[disease_id].hpo_ids = set(hpos)
    return diseases, disease_to_hpos


def parse_hpo_obo(path: Path) -> Dict[str, HpoTerm]:
    terms: Dict[str, HpoTerm] = {}
    current: Dict[str, object] = {}
    in_term = False

    def flush() -> None:
        nonlocal current, in_term
        if not in_term:
            return
        hpo_id = current.get("id")
        name = current.get("name")
        if hpo_id and name:
            terms[hpo_id] = HpoTerm(
                hpo_id=hpo_id,
                name=name,
                synonyms=list(current.get("synonyms", [])),
            )
        current = {}
        in_term = False

    with path.open(encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line == "[Term]":
                flush()
                in_term = True
                continue
            if not in_term:
                continue
            if not line.strip():
                continue
            if line.startswith("id: "):
                current["id"] = line[4:].strip()
            elif line.startswith("name: "):
                current["name"] = line[6:].strip()
            elif line.startswith("synonym: "):
                m = re.match(r'synonym:\s+"(.+?)"', line)
                if m:
                    current.setdefault("synonyms", []).append(m.group(1))
    flush()
    return terms


def parse_orphanet_xml(path: Path) -> Dict[str, DiseaseRecord]:
    records: Dict[str, DiseaseRecord] = {}
    context = ET.iterparse(path, events=("end",))
    for _, elem in context:
        if elem.tag != "Disorder":
            continue
        orpha_code = elem.findtext("OrphaCode", default="").strip()
        name = elem.findtext("Name", default="").strip()
        if not orpha_code or not name:
            elem.clear()
            continue
        disease_id = f"ORPHA:{orpha_code}"
        rec = records.get(disease_id)
        if rec is None:
            rec = DiseaseRecord(disease_id=disease_id, name=name, source="ORPHANET")
            records[disease_id] = rec
        else:
            if not rec.name:
                rec.name = name
        # Collect HPO IDs from the association list if present.
        for assoc in elem.findall("./HPODisorderAssociationList/HPODisorderAssociation"):
            hpo_id = assoc.findtext("./HPO/HPOId", default="").strip()
            if hpo_id:
                rec.hpo_ids.add(hpo_id)
        elem.clear()
    return records


def parse_mondo_obo(path: Path) -> Dict[str, DiseaseRecord]:
    records: Dict[str, DiseaseRecord] = {}
    current: Dict[str, object] = {}
    in_term = False

    def flush() -> None:
        nonlocal current, in_term
        if not in_term:
            return
        node_id = current.get("id")
        name = current.get("name")
        if not node_id or not name:
            current = {}
            in_term = False
            return

        xrefs = current.get("xrefs", [])
        if not xrefs:
            current = {}
            in_term = False
            return

        # Keep only disease-relevant terms with Orphanet/OMIM-style cross refs.
        disease_xrefs = [xref for xref in xrefs if xref.startswith("Orphanet:") or xref.startswith("OMIM:")]
        if not disease_xrefs:
            current = {}
            in_term = False
            return

        rec = DiseaseRecord(disease_id=node_id, name=name, source="MONDO")
        rec.synonyms.update(current.get("synonyms", []))
        for xref in disease_xrefs:
            rec.synonyms.add(xref)
        records[node_id] = rec
        current = {}
        in_term = False

    with path.open(encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line == "[Term]":
                flush()
                in_term = True
                continue
            if not in_term:
                continue
            if not line.strip():
                continue
            if line.startswith("id: "):
                current["id"] = line[4:].strip()
            elif line.startswith("name: "):
                current["name"] = line[6:].strip()
            elif line.startswith("synonym: "):
                m = re.match(r'synonym:\s+"(.+?)"', line)
                if m:
                    current.setdefault("synonyms", []).append(m.group(1))
            elif line.startswith("xref: "):
                xref = line[6:].split(" ", 1)[0].strip()
                current.setdefault("xrefs", []).append(xref)
    flush()
    return records


def build_disease_index(
    hpoa_records: Dict[str, DiseaseRecord],
    orphanet_records: Dict[str, DiseaseRecord],
    mondo_records: Dict[str, DiseaseRecord],
) -> Dict[str, DiseaseRecord]:
    index: Dict[str, DiseaseRecord] = {}

    def merge_record(src: DiseaseRecord) -> None:
        rec = index.get(src.disease_id)
        if rec is None:
            index[src.disease_id] = DiseaseRecord(
                disease_id=src.disease_id,
                name=src.name,
                source=src.source,
                synonyms=set(src.synonyms),
                hpo_ids=set(src.hpo_ids),
            )
            return
        if src.name and not rec.name:
            rec.name = src.name
        rec.synonyms.update(src.synonyms)
        rec.hpo_ids.update(src.hpo_ids)
        if rec.source != src.source:
            rec.source = f"{rec.source};{src.source}"

    for rec in hpoa_records.values():
        merge_record(rec)
    for rec in orphanet_records.values():
        merge_record(rec)
    for rec in mondo_records.values():
        merge_record(rec)

    # Add disease name synonyms from HPOA itself.
    for rec in hpoa_records.values():
        if rec.name:
            index[rec.disease_id].synonyms.add(rec.name)

    return index


def build_hpo_phrase_index(hpo_terms: Dict[str, HpoTerm]) -> Tuple[Dict[str, Set[str]], Dict[str, Dict[str, object]]]:
    token_index: Dict[str, Set[str]] = defaultdict(set)
    phrase_meta: Dict[str, Dict[str, object]] = {}

    def add_phrase(hpo_id: str, phrase: str, source: str) -> None:
        norm = normalize_text(phrase)
        if len(norm) < 3:
            return
        key = norm
        meta = phrase_meta.get(key)
        if meta is None:
            phrase_meta[key] = {
                "hpo_ids": {hpo_id},
                "phrase": phrase,
                "source": {source},
                "tokens": key.split(),
            }
        else:
            meta["hpo_ids"].add(hpo_id)
            meta["source"].add(source)
        for tok in set(key.split()):
            if tok not in STOPWORDS and len(tok) > 2:
                token_index[tok].add(key)

    for hpo_id, term in hpo_terms.items():
        add_phrase(hpo_id, term.name, "label")
        for syn in term.synonyms:
            add_phrase(hpo_id, syn, "synonym")

    return token_index, phrase_meta


def find_exact_phrase_spans(text_tokens: Sequence[str], phrase_tokens: Sequence[str]) -> List[Tuple[int, int]]:
    if not phrase_tokens or len(phrase_tokens) > len(text_tokens):
        return []
    spans: List[Tuple[int, int]] = []
    window = len(phrase_tokens)
    for i in range(len(text_tokens) - window + 1):
        if list(text_tokens[i : i + window]) == list(phrase_tokens):
            spans.append((i, i + window))
    return spans


def disease_lookup_candidates(query: str, disease_index: Dict[str, DiseaseRecord]) -> List[Tuple[float, DiseaseRecord, str]]:
    qnorm = normalize_text(query)
    qtokens = set(qnorm.split())
    scored: List[Tuple[float, DiseaseRecord, str]] = []

    for rec in disease_index.values():
        variants: List[str] = [rec.name, rec.disease_id]
        variants.extend(sorted(rec.synonyms))
        best = 0.0
        best_label = rec.name
        for variant in variants:
            vnorm = normalize_text(variant)
            if not vnorm:
                continue
            if qnorm == vnorm:
                score = 1.0
            elif qnorm in vnorm or vnorm in qnorm:
                score = 0.95
            else:
                vtokens = set(vnorm.split())
                inter = len(qtokens & vtokens)
                union = len(qtokens | vtokens) or 1
                jaccard = inter / union
                seq = SequenceMatcher(None, qnorm, vnorm).ratio()
                score = max(jaccard, seq * 0.9)
            if score > best:
                best = score
                best_label = variant
        if best > 0:
            scored.append((best, rec, best_label))

    scored.sort(key=lambda x: (-x[0], x[1].name.lower()))
    return scored


def resolve_disease_from_text(text: str, disease_index: Dict[str, DiseaseRecord]) -> Tuple[Optional[DiseaseRecord], float, List[Tuple[float, DiseaseRecord, str]], str]:
    """Resolve a disease from a multi-line / multi-clause text blob.

    The function tries the full text first, then each line and sentence-like
    chunk. This makes disease-name inputs resilient to paragraphs and OCR-like
    formatting where the disease name is embedded in a longer block.
    """
    chunks = collect_input_chunks(text)
    best_rec: Optional[DiseaseRecord] = None
    best_score = 0.0
    best_candidates: List[Tuple[float, DiseaseRecord, str]] = []
    best_chunk = ""

    for chunk in chunks:
        id_match = re.search(r"\b(ORPHA:\d+|OMIM:\d+|MONDO:\d+)\b", chunk, re.I)
        if id_match:
            did = id_match.group(1).upper()
            for rec in disease_index.values():
                if rec.disease_id.upper() == did:
                    return rec, 1.0, [(1.0, rec, rec.name)], chunk

        candidates = disease_lookup_candidates(chunk, disease_index)
        if not candidates:
            continue
        top_score, top_rec, top_label = candidates[0]
        if top_score > best_score:
            best_score = top_score
            best_rec = top_rec
            best_candidates = candidates[:10]
            best_chunk = chunk

            # If we have an exact or near-exact hit, stop early.
            if top_score >= 0.99:
                break

    return best_rec, best_score, best_candidates, best_chunk


def extract_hpo_terms_from_text(text: str, token_index: Dict[str, Set[str]], phrase_meta: Dict[str, Dict[str, object]]) -> List[Dict[str, object]]:
    text_norm = normalize_text(text)
    text_tokens = text_norm.split()

    candidate_phrases: Set[str] = set()
    for tok in set(text_tokens):
        candidate_phrases.update(token_index.get(tok, set()))

    # Prefer longer clinical phrases first so we keep "proximal muscle weakness"
    # rather than only the shorter nested phrase "muscle weakness".
    ordered_candidates = sorted(
        candidate_phrases,
        key=lambda k: (-(len(phrase_meta[k]["tokens"]) if k in phrase_meta else 0), k),
    )

    accepted_spans: List[Tuple[int, int]] = []
    seen_hpo: Dict[str, Dict[str, object]] = {}

    for phrase_key in ordered_candidates:
        meta = phrase_meta.get(phrase_key)
        if not meta:
            continue
        phrase_tokens = meta["tokens"]
        if len(phrase_tokens) == 1 and phrase_tokens[0] in GENERIC_SINGLE_WORD_BLACKLIST:
            continue
        spans = find_exact_phrase_spans(text_tokens, phrase_tokens)
        if not spans:
            continue

        matched_source = ",".join(sorted(meta["source"]))
        phrase_len = len(phrase_tokens)
        base_score = 0.92 if "label" in meta["source"] else 0.88
        score = min(1.0, base_score + min(0.08, phrase_len * 0.01))

        for start, end in spans:
            if any(not (end <= a or start >= b) for a, b in accepted_spans):
                continue
            accepted_spans.append((start, end))
            for hpo_id in meta["hpo_ids"]:
                prev = seen_hpo.get(hpo_id)
                if prev and prev["match_score"] >= score:
                    continue
                seen_hpo[hpo_id] = {
                    "hpo_id": hpo_id,
                    "matched_phrase": meta["phrase"],
                    "normalized_phrase": phrase_key,
                    "match_score": round(score, 4),
                    "source": matched_source,
                    "span_start": start,
                    "span_end": end,
                }

    matches = list(seen_hpo.values())
    matches.sort(key=lambda x: (-x["match_score"], x["span_start"], x["hpo_id"]))
    return matches


def load_hpo_name_map(hpo_terms: Dict[str, HpoTerm]) -> Dict[str, str]:
    return {hpo_id: term.name for hpo_id, term in hpo_terms.items()}


def write_tsv(path: Path, rows: Sequence[Dict[str, object]], fieldnames: Sequence[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def write_hpo_ids_only(path: Path, hpo_ids: Sequence[str]) -> None:
    with path.open("w", encoding="utf-8") as fh:
        for hpo_id in hpo_ids:
            fh.write(f"{hpo_id}\n")


def score_germline_matches(matches: Sequence[Dict[str, object]]) -> float:
    return round(sum(float(m.get("match_score", 0.0)) for m in matches[:10]), 4)


def infer_mode(
    input_text: str,
    disease_index: Dict[str, DiseaseRecord],
    token_index: Dict[str, Set[str]],
    phrase_meta: Dict[str, Dict[str, object]],
) -> Tuple[str, Dict[str, object]]:
    resolved, disease_score, candidates, matched_chunk = resolve_disease_from_text(input_text, disease_index)
    germline_matches = extract_hpo_terms_from_text(input_text, token_index, phrase_meta)
    germline_score = score_germline_matches(germline_matches)

    if resolved and disease_score >= 0.98 and len(germline_matches) <= 2:
        inferred_mode = "somatic"
    elif len(germline_matches) >= 3 and germline_score >= 2.0:
        inferred_mode = "germline"
    elif resolved and disease_score >= 0.90 and germline_score < 1.8:
        inferred_mode = "somatic"
    else:
        inferred_mode = "germline"

    details = {
        "resolved": resolved,
        "disease_score": round(disease_score, 4),
        "disease_candidates": candidates[:10],
        "matched_chunk": matched_chunk,
        "germline_matches": germline_matches,
        "germline_score": germline_score,
    }
    return inferred_mode, details


def run_somatic(
    input_text: str,
    disease_index: Dict[str, DiseaseRecord],
    hpo_name_map: Dict[str, str],
    outdir: Path,
    detection: Dict[str, object],
) -> Dict[str, object]:
    resolved = detection["resolved"]
    score = float(detection["disease_score"])
    candidates = detection["disease_candidates"]
    matched_chunk = detection["matched_chunk"]
    result_rows: List[Dict[str, object]] = []
    summary: Dict[str, object] = {
        "mode": "somatic",
        "inferred_mode": "somatic",
        "input": input_text.strip(),
        "resolved": None,
        "confidence": round(score, 4),
        "matched_chunk": matched_chunk,
        "disease_score": round(score, 4),
        "germline_score": detection["germline_score"],
        "top_candidates": [],
        "matched_hpo_count": 0,
    }

    if resolved:
        hpo_ids = sorted(resolved.hpo_ids)
        for hpo_id in hpo_ids:
            result_rows.append(
                {
                    "mode": "somatic",
                    "query": input_text.strip(),
                    "resolved_id": resolved.disease_id,
                    "resolved_name": resolved.name,
                    "resolved_source": resolved.source,
                    "confidence": round(score, 4),
                    "hpo_id": hpo_id,
                    "hpo_name": hpo_name_map.get(hpo_id, ""),
                    "match_type": "disease_annotation",
                }
            )
        summary["resolved"] = {
            "disease_id": resolved.disease_id,
            "disease_name": resolved.name,
            "source": resolved.source,
        }
        summary["matched_hpo_count"] = len(hpo_ids)
        summary["hpo_ids"] = hpo_ids

    summary["top_candidates"] = [
        {
            "disease_id": rec.disease_id,
            "disease_name": rec.name,
            "score": round(sc, 4),
            "matched_variant": label,
        }
        for sc, rec, label in candidates[:10]
    ]

    out_tsv = outdir / "resolution_matches.tsv"
    fieldnames = [
        "mode",
        "query",
        "resolved_id",
        "resolved_name",
        "resolved_source",
        "confidence",
        "hpo_id",
        "hpo_name",
        "match_type",
    ]
    write_tsv(out_tsv, result_rows, fieldnames)

    out_json = outdir / "resolution_summary.json"
    out_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    out_ids = outdir / "hpo_ids.txt"
    write_hpo_ids_only(out_ids, sorted({row["hpo_id"] for row in result_rows}))
    return summary


def run_germline(
    input_text: str,
    token_index: Dict[str, Set[str]],
    phrase_meta: Dict[str, Dict[str, object]],
    hpo_name_map: Dict[str, str],
    outdir: Path,
    detection: Dict[str, object],
) -> Dict[str, object]:
    matches = detection["germline_matches"]
    rows: List[Dict[str, object]] = []
    for m in matches:
        rows.append(
            {
                "mode": "germline",
                "query": input_text.strip(),
                "resolved_id": "",
                "resolved_name": "",
                "resolved_source": "",
                "confidence": "",
                "hpo_id": m["hpo_id"],
                "hpo_name": hpo_name_map.get(m["hpo_id"], ""),
                "match_type": "phenotype_extraction",
                "matched_phrase": m["matched_phrase"],
                "normalized_phrase": m["normalized_phrase"],
                "match_score": m["match_score"],
                "source": m["source"],
            }
        )

    summary = {
        "mode": "germline",
        "inferred_mode": "germline",
        "input": input_text.strip(),
        "matched_hpo_count": len(matches),
        "top_matches": matches[:20],
        "hpo_ids": [m["hpo_id"] for m in matches],
        "disease_score": detection["disease_score"],
        "germline_score": detection["germline_score"],
    }

    out_tsv = outdir / "resolution_matches.tsv"
    fieldnames = [
        "mode",
        "query",
        "resolved_id",
        "resolved_name",
        "resolved_source",
        "confidence",
        "hpo_id",
        "hpo_name",
        "match_type",
        "matched_phrase",
        "normalized_phrase",
        "match_score",
        "source",
    ]
    write_tsv(out_tsv, rows, fieldnames)

    out_json = outdir / "resolution_summary.json"
    out_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    out_ids = outdir / "hpo_ids.txt"
    write_hpo_ids_only(out_ids, sorted({m["hpo_id"] for m in matches}))
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Resolve HPO terms from disease names or phenotype descriptions.")
    parser.add_argument("--input", required=True, help="Input text file")
    parser.add_argument("--outdir", default="hpo_resolution", help="Output directory")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not HPOA_PATH.exists():
        raise FileNotFoundError(f"Missing file: {HPOA_PATH}")
    if not HPO_OBO_PATH.exists():
        raise FileNotFoundError(f"Missing file: {HPO_OBO_PATH}")
    if not ORPHA_XML_PATH.exists():
        raise FileNotFoundError(f"Missing file: {ORPHA_XML_PATH}")
    if not MONDO_OBO_PATH.exists():
        raise FileNotFoundError(f"Missing file: {MONDO_OBO_PATH}")

    hpoa_records, _ = parse_hpoa(HPOA_PATH)
    hpo_terms = parse_hpo_obo(HPO_OBO_PATH)
    orphanet_records = parse_orphanet_xml(ORPHA_XML_PATH)
    mondo_records = parse_mondo_obo(MONDO_OBO_PATH)
    disease_index = build_disease_index(hpoa_records, orphanet_records, mondo_records)
    token_index, phrase_meta = build_hpo_phrase_index(hpo_terms)
    hpo_name_map = load_hpo_name_map(hpo_terms)

    input_text = read_single_value(Path(args.input))
    inferred_mode, detection = infer_mode(input_text, disease_index, token_index, phrase_meta)

    if inferred_mode == "somatic":
        summary = run_somatic(input_text, disease_index, hpo_name_map, outdir, detection)
    else:
        summary = run_germline(input_text, token_index, phrase_meta, hpo_name_map, outdir, detection)

    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
