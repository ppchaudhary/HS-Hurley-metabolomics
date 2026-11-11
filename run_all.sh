#!/usr/bin/env bash
# Reproduce all analyses from the raw feature tables in data/. Run from the repository root.
set -euo pipefail
mkdir -p results
Rscript scripts/01_serum_main.R
Rscript scripts/02_skin_main.R
Rscript scripts/03_skin_followup.R
python3 scripts/04_pathway_mummichog.py --sample Skin  --tables results/skin/tables  --meta data/Metadata_skin.csv  --out results/pathway_skin
python3 scripts/04_pathway_mummichog.py --sample Serum --tables results/serum/tables --meta data/Metadata_serum.csv --out results/pathway_serum
python3 scripts/05_compare_skin_serum.py results/skin/tables/T04_trend_with_Hurley_stage.csv \
        results/serum/tables/T04_trend_with_Hurley_stage.csv results/pathway_skin results/pathway_serum results/comparison
echo "All analyses complete: see results/"
