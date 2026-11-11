# Metabolomic survey of skin and serum across Hurley stages of hidradenitis suppurativa

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

Analysis code for the research letter *"A descriptive metabolomic survey of skin and serum across Hurley stages of hidradenitis suppurativa"* (Chaudhary PP, et al.; submitted to the *Journal of Investigative Dermatology*).

Running the pipeline on the feature tables in `data/` reproduces every number, figure and supplementary table in the letter.

## Contents

| Path | Purpose |
|---|---|
| `scripts/01_serum_main.R` | Serum analysis (45 patients): QC, PQN normalization, PCA/NMDS/PERMANOVA, stage trend (limma, Spearman, Jonckheere–Terpstra), pairwise volcano plots, heatmaps, lipid classes, Random Forest, sPLS-DA, sensitivity analyses |
| `scripts/02_skin_main.R` | Skin analysis (15 patients), same pipeline and quality-control rule |
| `scripts/03_skin_followup.R` | Skin follow-up: background-ion removal (mass defect), isotope/adduct/split-peak merging, robust candidates, lipid ratios, chain length and saturation |
| `scripts/04_pathway_mummichog.py` | Pathway enrichment with mummichog (as in MetaboAnalyst Functional Analysis) and lipid-class enrichment by label permutation; run once per sample type |
| `scripts/05_compare_skin_serum.py` | Skin–serum feature matching (5 ppm) and side-by-side pathway and lipid-class comparison |
| `run_all.sh` | Runs all steps in order |
| `data/` | Input feature tables and metadata (see `data/README.md`) |

## Requirements

- R ≥ 4.5 with limma, vegan, clinfun, randomForest, Boruta, mixOmics, ComplexHeatmap, EnhancedVolcano, circlize, ggplot2, ggrepel, patchwork, dplyr, tidyr, scales, ragg and systemfonts. Missing packages are installed automatically from CRAN/Bioconductor on first run; the versions used are listed in `R_packages.txt`.
- Python ≥ 3.9: `pip install -r requirements.txt` (this installs the `mummichog` command, which must be on your PATH).
- A bash shell (macOS/Linux; on Windows use WSL or Git Bash).

## Usage

```bash
git clone https://github.com/ppchaudhary/HS-Hurley-metabolomics.git
cd HS-Hurley-metabolomics
pip install -r requirements.txt
./run_all.sh            # about 40 minutes; most of the time is Random Forest permutations
```

Steps can also be run individually from the repository root, in the order listed in `run_all.sh`.

Results are written to `results/`: `serum`, `skin`, `skin_followup`, `pathway_skin` and `pathway_serum` each contain `figures/` (PDF, PNG and 600 dpi TIFF) and `tables/` (CSV); `comparison` contains CSV tables.

## Key settings

Thresholds are set at the top of each script:

- detection filter: ≥ 80% of samples in at least one stage;
- significance: FDR < 0.05 (Benjamini–Hochberg);
- volcano cut-offs: p < 0.05 and |log2 FC| > 1;
- random seed: 2026.

Scripts 01 and 02 flag samples automatically: low correlation (median correlation < cohort median – 3 MAD) and low signal (total ion signal < ¼ of the median). Script 03 defines robust candidates by excluding the two lowest-signal skin samples (9-MT, 10-MT).

## Reproducibility

The full pipeline was rerun from the input files in a clean directory, and all reported results were reproduced (see `REPRODUCIBILITY.md`). Most stochastic steps use a fixed seed. The mummichog permutations and the PERMANOVA tests within the sensitivity summaries are not seeded, so their p-values can differ in the third decimal place between runs, without changing any conclusion. Small differences can also occur with other package versions.

## Data

The input feature tables are not included in this repository. They are available from the corresponding author on reasonable request. Place them in `data/` as described in `data/README.md` before running the pipeline.

## Citation

If you use this code, please cite the letter (citation to be added on publication) and the archived code (DOI above). Please also cite the software listed in the Supplementary Materials and Methods.

## Contact

Prem Prashant Chaudhary, NIAID, National Institutes of Health — [chaudharyp2@nih.gov]

## License

MIT (see `LICENSE`).