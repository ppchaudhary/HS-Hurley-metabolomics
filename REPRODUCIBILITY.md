# Reproducibility check

The complete pipeline (`run_all.sh`) was rerun on 10 November 2026 in a clean directory containing only the four input files in `data/`, using R 4.5.3 and Python 3 with mummichog 2.7.0. It ran without errors (about 40 minutes), and its outputs were compared with the results reported in the letter.

| Output | Result of rerun |
|---|---|
| Serum trend statistics (1,997 features; 97 at FDR < 0.05) | Identical |
| Serum pairwise contrasts, PERMANOVA, Random Forest importance | Identical |
| Serum Random Forest: 79% accuracy, pseudo-R² 0.60, permutation p = 0.005 | Identical |
| Skin trend statistics (3,010 features; 0 at FDR < 0.05) | Identical |
| Skin follow-up: 1,831 compounds, 23 robust candidates, chain/saturation trends | Identical |
| Skin–serum matching (232 features; 64% same direction; ρ = 0.21) | Identical |
| Lipid-class enrichment (10,000 permutations, seeded) | Identical |
| PERMANOVA p-values within the sensitivity summaries (T10) | Differ by ≤ 0.007 (unseeded permutations) |
| mummichog pathway p-values | Differ by ≤ 0.003 (mummichog permutations cannot be seeded); same pathways and hits |

All conclusions are unchanged. Exact values can differ slightly with other package versions; the versions used are listed in `R_packages.txt` and `requirements.txt`.
