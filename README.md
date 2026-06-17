# chatbot-exomiser-service

Exomiser worker scripts for Geneie (Germline filter_3).

## Layout
- Phenotype_HPO_Extractor/resolve_hpo_terms.py — phenotype → HPO IDs
- exomiser-cli-15.0.0/run_exomiser.py — Exomiser CLI wrapper

## Manual test (on this host)
# Step 1: HPO
cd Phenotype_HPO_Extractor
python3 resolve_hpo_terms.py --input input_dis.txt

# Step 2: Exomiser (hg38 example)
cd ../exomiser-cli-15.0.0
python3 run_exomiser.py --vcf <path-to-vcf> --hpo ../Phenotype_HPO_Extractor/hpo_resolution/hpo_ids.txt --inheritance AUTOSOMAL_DOMINANT --assembly hg38
