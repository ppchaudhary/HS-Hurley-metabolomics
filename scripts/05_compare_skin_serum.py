"""
Skin vs serum comparison (HS, Hurley stages).
Inputs : trend tables from the skin and serum R analyses (T04), and the
         pathway folders made by HS_pathway_mummichog.py.
Outputs: Comparison_Skin_vs_Serum/ matched_features_5ppm.csv,
         pathway_comparison.csv, lipid_class_enrichment_comparison.csv
Usage  : python HS_compare_skin_serum.py SKIN_T04.csv SERUM_T04.csv Pathway_Skin Pathway_Serum
"""
import sys, os, numpy as np, pandas as pd
from scipy.stats import spearmanr
sk_f, se_f, pk, ps = sys.argv[1:5]
OUT = sys.argv[5] if len(sys.argv) > 5 else "Comparison_Skin_vs_Serum"
os.makedirs(OUT, exist_ok=True)
sk, se = pd.read_csv(sk_f), pd.read_csv(se_f)
sem = se.sort_values("mz").reset_index(drop=True); arr = sem.mz.values
rows = []
for _, r in sk.iterrows():                       # nearest serum feature within 5 ppm
    tol = r.mz * 5e-6
    lo, hi = np.searchsorted(arr, r.mz - tol), np.searchsorted(arr, r.mz + tol, "right")
    if hi > lo:
        s = sem.iloc[lo + np.argmin(abs(arr[lo:hi] - r.mz))]
        rows.append(dict(mz_skin=r.mz, mz_serum=s.mz, Label_skin=r.Label, Label_serum=s.Label,
                         Class=r.Class, fc_skin=r.log2FC_per_stage, p_skin=r.P_limma,
                         fc_serum=s.log2FC_per_stage, p_serum=s.P_limma, fdr_serum=s.FDR_limma))
m = pd.DataFrame(rows); m.to_csv(f"{OUT}/matched_features_5ppm.csv", index=False)
rho, p = spearmanr(m.fc_skin, m.fc_serum)
print(f"matched {len(m)}; same direction {np.mean(np.sign(m.fc_skin)==np.sign(m.fc_serum)):.2f}; rho {rho:.2f} p {p:.3g}")
# pathways side by side
recs = []
for lab, d in (("Skin", pk), ("Serum", ps)):
    t = pd.read_csv(f"{d}/tables/pathways_all_comparisons.csv")
    t = t[t.overlap_size > 0]
    for _, r in t.iterrows():
        recs.append((r.pathway, lab, r.comparison.split(" (")[0],
                     f"{int(r.overlap_size)}/{int(r.pathway_size)}, p={r['p-value']:.3g}"))
d = pd.DataFrame(recs, columns=["Pathway", "Sample", "Comparison", "cell"])
w = d.pivot_table(index="Pathway", columns=["Sample", "Comparison"], values="cell", aggfunc="first").fillna("-")
w.columns = [f"{a}: {b}" for a, b in w.columns]
w["In both sample types"] = w.index.isin(set(d[d.Sample == "Skin"].Pathway) & set(d[d.Sample == "Serum"].Pathway))
w.reset_index().to_csv(f"{OUT}/pathway_comparison.csv", index=False)
# lipid-class enrichment side by side
a = pd.read_csv(f"{pk}/tables/lipid_class_enrichment.csv"); b = pd.read_csv(f"{ps}/tables/lipid_class_enrichment.csv")
c = a[["Class", "n_features", "Shift_vs_other_features", "P_permutation"]].merge(
    b[["Class", "n_features", "Shift_vs_other_features", "P_permutation"]], on="Class", how="outer", suffixes=("_skin", "_serum"))
for s in ("skin", "serum"):
    c[f"Direction_{s}"] = c[f"Shift_vs_other_features_{s}"].apply(lambda v: "" if pd.isna(v) else ("Up" if v > 0 else "Down"))
c["Same_direction"] = (c.Direction_skin == c.Direction_serum) & (c.Direction_skin != "")
c.to_csv(f"{OUT}/lipid_class_enrichment_comparison.csv", index=False)
print("Done:", OUT)
