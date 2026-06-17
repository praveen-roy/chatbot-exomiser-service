#!/usr/bin/env python3
"""
run_exomiser.py — Command-line wrapper for Exomiser 15+

Usage:
    python run_exomiser.py --vcf sample.vcf.gz --hpo hpo_ids.txt [options]

Minimal example:
    python run_exomiser.py \
        --vcf /path/to/sample.vcf.gz \
        --hpo /path/to/hpo_ids.txt

Full example:
    python run_exomiser.py \
        --vcf /path/to/sample.vcf.gz \
        --hpo /path/to/hpo_ids.txt \
        --assembly hg38 \
        --mode PASS_ONLY \
        --inheritance AUTOSOMAL_DOMINANT AUTOSOMAL_RECESSIVE_HOM_ALT AUTOSOMAL_RECESSIVE_COMP_HET X_DOMINANT X_RECESSIVE_HOM_ALT X_RECESSIVE_COMP_HET MITOCHONDRIAL\
        --prioritiser hiphive \
        --output-formats HTML JSON TSV_GENE TSV_VARIANT VCF \
        --out-dir ./results \
        --out-prefix my_sample \
        --work-dir ~/exomiser_runs \
        --jar /home/ubuntu/exomiser/exomiser-cli-15.0.0/exomiser-cli-15.0.0.jar \
        --data-dir /home/ubuntu/exomiser/exomiser-cli-15.0.0/data \
        --data-version 2512 \
        --transcript-source ensembl \
        --dry-run

HPO .txt file format (one ID per line, comments with # ignored):
    HP:0001250
    HP:0001263
    # This is a comment
    0000252        <- HP: prefix is optional
"""

import argparse
import os
import sys
import subprocess
import glob
import yaml
from pathlib import Path
from datetime import datetime

# ── Fixed analysis defaults ───────────────────────────────────────────────────
FREQ_SOURCES = [
    "UK10K",
    "GNOMAD_E_AFR", "GNOMAD_E_AMR", "GNOMAD_E_EAS", "GNOMAD_E_NFE", "GNOMAD_E_SAS",
    "GNOMAD_G_AFR", "GNOMAD_G_AMR", "GNOMAD_G_EAS", "GNOMAD_G_NFE", "GNOMAD_G_SAS",
]
PATH_SOURCES = ["REVEL", "MVP", "ALPHA_MISSENSE", "SPLICE_AI"]
REMOVE_EFFECTS = [
    "FIVE_PRIME_UTR_EXON_VARIANT", "FIVE_PRIME_UTR_INTRON_VARIANT",
    "THREE_PRIME_UTR_EXON_VARIANT", "THREE_PRIME_UTR_INTRON_VARIANT",
    "NON_CODING_TRANSCRIPT_EXON_VARIANT", "NON_CODING_TRANSCRIPT_INTRON_VARIANT",
    "CODING_TRANSCRIPT_INTRON_VARIANT",
    "UPSTREAM_GENE_VARIANT", "DOWNSTREAM_GENE_VARIANT",
    "INTERGENIC_VARIANT", "REGULATORY_REGION_VARIANT",
]
MAF = {
    "AUTOSOMAL_DOMINANT":           0.1,
    "AUTOSOMAL_RECESSIVE_HOM_ALT":  0.1,
    "AUTOSOMAL_RECESSIVE_COMP_HET": 2.0,
    "X_DOMINANT":                   0.1,
    "X_RECESSIVE_HOM_ALT":          0.1,
    "X_RECESSIVE_COMP_HET":         2.0,
    "MITOCHONDRIAL":                0.2,
}
ALL_INHERITANCE_MODES = list(MAF.keys())
DEFAULT_INHERITANCE_MODES = [
    "AUTOSOMAL_DOMINANT",
    "AUTOSOMAL_RECESSIVE_HOM_ALT",
    "AUTOSOMAL_RECESSIVE_COMP_HET",
    "X_DOMINANT",
    "X_RECESSIVE_HOM_ALT",
    "X_RECESSIVE_COMP_HET",
]

# ── Defaults ──────────────────────────────────────────────────────────────────
DEFAULT_JAR       = "/home/ubuntu/exomiser/exomiser-cli-15.0.0/exomiser-cli-15.0.0.jar"
DEFAULT_DATA_DIR  = "/home/ubuntu/exomiser/exomiser-cli-15.0.0/data"
DEFAULT_DATA_VER  = "2512"
DEFAULT_ASSEMBLY  = "hg38"
DEFAULT_TRANSCRIPT= "ensembl"
DEFAULT_WORK_DIR  = os.path.expanduser("/home/ubuntu/exomiser/exomiser-cli-15.0.0")
DEFAULT_OUT_DIR   = os.path.expanduser("/home/ubuntu/exomiser/exomiser-cli-15.0.0/results")
DEFAULT_OUT_FMT   = ["HTML", "JSON", "TSV_GENE", "TSV_VARIANT", "VCF"]


# ─────────────────────────────────────────────────────────────────────────────
def parse_hpo_file(hpo_file: str) -> list[str]:
    """
    Parse a .txt file of HPO IDs. Supports:
      - One ID per line
      - Lines starting with # are treated as comments and skipped
      - HP: prefix is optional (will be added automatically)
      - Blank lines are skipped
    """
    hpo_ids = []
    skipped = []
    path = Path(hpo_file)

    if not path.exists():
        print(f"[ERROR] HPO file not found: {hpo_file}", file=sys.stderr)
        sys.exit(1)

    with open(path) as f:
        for lineno, line in enumerate(f, 1):
            raw = line.strip()
            if not raw or raw.startswith("#"):
                continue
            # Normalise: strip surrounding quotes, force uppercase
            raw = raw.strip("'\"").upper()
            if not raw.startswith("HP:"):
                raw = "HP:" + raw
            # Validate format: HP: + exactly 7 digits
            if len(raw) == 10 and raw[3:].isdigit():
                hpo_ids.append(raw)
            else:
                skipped.append((lineno, line.strip()))

    if skipped:
        print(f"[WARN] Skipped {len(skipped)} invalid line(s) in HPO file:")
        for lineno, val in skipped:
            print(f"       Line {lineno}: '{val}'")

    if not hpo_ids:
        print("[ERROR] No valid HPO IDs found in file.", file=sys.stderr)
        sys.exit(1)

    return hpo_ids


def build_yaml(args, vcf_abs: str, hpo_ids: list[str],
               results_dir: str, ped_abs: str | None) -> dict:
    """Build the Exomiser analysis YAML dict."""

    inh_modes = {m: MAF[m] for m in args.inheritance}

    steps = [
        {"geneBlacklistFilter":  {}},
        {"failedVariantFilter":  {}},
        {"alleleBalanceFilter":  {}},
        {"variantEffectFilter":  {"remove": REMOVE_EFFECTS}},
        {"frequencyFilter":      {"maxFrequency": 2.0}},
        {"pathogenicityFilter":  {"keepNonPathogenic": True}},
        {"inheritanceFilter":    {}},
        {"omimPrioritiser":      {}},
    ]
    p = args.prioritiser.lower()
    if p == "hiphive": steps.append({"hiPhivePrioritiser": {}})
    elif p == "phive": steps.append({"phivePrioritiser":   {}})
    elif p == "phenix": steps.append({"phenixPrioritiser": {}})

    analysis_section = {
        "genomeAssembly":       args.assembly,
        "vcf":                  vcf_abs,
        "hpoIds":               hpo_ids,
        "inheritanceModes":     inh_modes,
        "analysisMode":         args.mode,
        "frequencySources":     FREQ_SOURCES,
        "pathogenicitySources": PATH_SOURCES,
        "steps":                steps,
    }
    if ped_abs:
        analysis_section["ped"] = ped_abs
    if args.proband:
        analysis_section["proband"] = args.proband

    return {
        "analysis": analysis_section,
        "outputOptions": {
            "outputContributingVariantsOnly": False,
            "numGenes":        0,
            "outputDirectory": results_dir,
            "outputFileName":  args.out_prefix,
            "outputFormats":   args.output_formats,
        },
    }


def resolve_jar(jar_pattern: str) -> str:
    matches = glob.glob(jar_pattern)
    if matches:
        return matches[0]
    if os.path.isfile(jar_pattern):
        return jar_pattern
    print(f"[ERROR] JAR not found: {jar_pattern}", file=sys.stderr)
    sys.exit(1)


def log(msg: str):
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}")


# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Command-line Exomiser 15 runner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Required
    req = parser.add_argument_group("required arguments")
    req.add_argument("--vcf",  required=True, help="Input VCF file (.vcf or .vcf.gz)")
    req.add_argument("--hpo",  required=True, help="HPO ID list (.txt, one ID per line)")

    # Exomiser paths
    paths = parser.add_argument_group("Exomiser paths")
    paths.add_argument("--jar",              default=DEFAULT_JAR,
                       help=f"Path to Exomiser JAR (default: {DEFAULT_JAR})")
    paths.add_argument("--data-dir",         default=DEFAULT_DATA_DIR,
                       dest="data_dir",
                       help=f"Exomiser data directory (default: {DEFAULT_DATA_DIR})")
    paths.add_argument("--data-version",     default=DEFAULT_DATA_VER,
                       dest="data_version",
                       help=f"Data bundle version prefix (default: {DEFAULT_DATA_VER})")
    paths.add_argument("--transcript-source",default=DEFAULT_TRANSCRIPT,
                       dest="transcript_source",
                       choices=["ensembl", "refseq", "ucsc"],
                       help=f"Transcript source (default: {DEFAULT_TRANSCRIPT})")
    paths.add_argument("--work-dir",         default=DEFAULT_WORK_DIR,
                       dest="work_dir",
                       help=f"Directory where YAML is written (default: {DEFAULT_WORK_DIR})")
    paths.add_argument("--out-dir",          default=DEFAULT_OUT_DIR,
                       dest="out_dir",
                       help=f"Results output directory (default: {DEFAULT_OUT_DIR})")
    paths.add_argument("--out-prefix",       default="exomiser_results",
                       dest="out_prefix",
                       help="Output file name prefix (default: exomiser_results)")

    # Optional input
    opt = parser.add_argument_group("optional inputs")
    opt.add_argument("--ped",      default=None, help="PED file for family/trio analysis")
    opt.add_argument("--proband",  default=None, help="Proband sample ID")

    # Analysis settings
    ana = parser.add_argument_group("analysis settings")
    ana.add_argument("--assembly",  default=DEFAULT_ASSEMBLY,
                     choices=["hg38", "hg19"],
                     help=f"Genome assembly (default: {DEFAULT_ASSEMBLY})")
    ana.add_argument("--mode",      default="PASS_ONLY",
                     choices=["PASS_ONLY", "FULL"],
                     help="Analysis mode (default: PASS_ONLY)")
    ana.add_argument("--inheritance", nargs="+",
                     default=DEFAULT_INHERITANCE_MODES,
                     choices=ALL_INHERITANCE_MODES,
                     metavar="MODE",
                     help=(
                         "Inheritance modes to analyse. Default: AD, AR_HOM, AR_CH, XD, XR_HOM, XR_CH.\n"
                         f"Available: {', '.join(ALL_INHERITANCE_MODES)}"
                     ))
    ana.add_argument("--prioritiser", default="hiphive",
                     choices=["hiphive", "phive", "phenix", "none"],
                     help="Prioritiser to use (default: hiphive)")
    ana.add_argument("--output-formats", nargs="+",
                     default=DEFAULT_OUT_FMT,
                     dest="output_formats",
                     choices=["HTML", "JSON", "TSV_GENE", "TSV_VARIANT", "VCF", "PARQUET"],
                     metavar="FORMAT",
                     help=f"Output formats (default: {' '.join(DEFAULT_OUT_FMT)})")

    # JVM memory
    jvm = parser.add_argument_group("JVM settings")
    jvm.add_argument("--xms", default="4g", help="JVM initial heap size (default: 4g)")
    jvm.add_argument("--xmx", default="8g", help="JVM max heap size (default: 8g)")

    # Misc
    parser.add_argument("--dry-run", action="store_true",
                        help="Generate YAML only — do not run Exomiser")
    parser.add_argument("--yaml-only", action="store_true",
                        help="Write analysis.yml and print the run command, then exit")

    args = parser.parse_args()

    # ── Resolve paths ─────────────────────────────────────────────────────────
    vcf_abs  = os.path.abspath(os.path.expanduser(args.vcf))
    work_dir = os.path.abspath(os.path.expanduser(args.work_dir))
    out_dir  = os.path.abspath(os.path.expanduser(args.out_dir))

    if not os.path.isfile(vcf_abs):
        print(f"[ERROR] VCF not found: {vcf_abs}", file=sys.stderr)
        sys.exit(1)

    ped_abs = None
    if args.ped:
        ped_abs = os.path.abspath(os.path.expanduser(args.ped))
        if not os.path.isfile(ped_abs):
            print(f"[ERROR] PED not found: {ped_abs}", file=sys.stderr)
            sys.exit(1)

    # ── Parse HPO file ────────────────────────────────────────────────────────
    log(f"Reading HPO IDs from: {args.hpo}")
    hpo_ids = parse_hpo_file(args.hpo)
    log(f"Loaded {len(hpo_ids)} HPO ID(s): {', '.join(hpo_ids)}")

    # ── Create directories ────────────────────────────────────────────────────
    os.makedirs(work_dir, exist_ok=True)
    os.makedirs(out_dir,  exist_ok=True)

    # ── Build & write YAML ────────────────────────────────────────────────────
    analysis_yaml = build_yaml(args, vcf_abs, hpo_ids, out_dir, ped_abs)
    yaml_path = os.path.join(work_dir, "analysis.yml")

    with open(yaml_path, "w") as fh:
        yaml.dump(analysis_yaml, fh, default_flow_style=False,
                  sort_keys=False, allow_unicode=True)

    log(f"YAML written → {yaml_path}")

    # ── Build run command ─────────────────────────────────────────────────────
    jar = resolve_jar(args.jar)

    cmd = [
        "java",
        f"-Xms{args.xms}",
        f"-Xmx{args.xmx}",
        f"-Dexomiser.data-directory={args.data_dir}",
        f"-Dexomiser.{args.assembly}.data-version={args.data_version}",
        f"-Dexomiser.{args.assembly}.transcript-source={args.transcript_source}",
        "-jar", jar,
        "analyse",
        "--analysis", yaml_path,
    ]

    log("Command:\n  " + " \\\n  ".join(cmd))

    if args.dry_run or args.yaml_only:
        log("Dry-run mode — skipping execution.")
        sys.exit(0)

    # ── Execute ───────────────────────────────────────────────────────────────
    log("Starting Exomiser…")
    print("-" * 60)
    try:
        proc = subprocess.run(cmd, text=True)
    except FileNotFoundError:
        print(f"\n[ERROR] java not found. Is Java 11+ installed?", file=sys.stderr)
        sys.exit(1)
    print("-" * 60)

    if proc.returncode == 0:
        log("✅  Exomiser finished successfully.")
        log(f"Results → {out_dir}/")
        result_files = sorted(Path(out_dir).glob(f"{args.out_prefix}*"))
        if result_files:
            log("Output files:")
            for f in result_files:
                print(f"   {f}")
    else:
        log(f"❌  Exomiser exited with code {proc.returncode}")
        sys.exit(proc.returncode)


if __name__ == "__main__":
    main()
