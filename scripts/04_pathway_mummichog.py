"""
HS metabolomics - pathway analysis with mummichog (the algorithm behind
MetaboAnalyst "Functional Analysis" / MS Peaks to Pathways), for SKIN or SERUM.

Uses the result tables of the main R analysis (T02, T04, T05) and builds one
self-contained folder per sample type:

  Pathway_<Sample>/
    inputs/          MetaboAnalyst-ready peak lists (m.z, p.value, t.score)
    mummichog_runs/  raw mummichog output per comparison (incl. its HTML report)
    tables/          pathway results, matched compounds, lipid-class enrichment
    figures/         MetaboAnalyst-style bubble plots, lipid-class enrichment
    scripts/         a copy of this script

Part A  mummichog: m/z -> human metabolic network (MFN), pathway
        over-representation among significant features (p < 0.05)
Part B  Lipid-class enrichment along the stage trend: competitive test
        (class vs other features), significance from 10,000 permutations of
        patients' stage labels (respects correlation between features)

Requires: pip install mummichog pandas numpy matplotlib
Usage:
  python HS_pathway_mummichog.py --sample Skin  --tables HS_skin_results_autoQC/tables --meta Metadata.csv
  python HS_pathway_mummichog.py --sample Serum --tables HS_serum_results/tables      --meta Metadata_serum.csv
"""
import os, glob, shutil, subprocess, argparse, textwrap, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--sample", required=True)            # label, e.g. Skin / Serum
ap.add_argument("--tables", required=True)            # folder with T02, T04, T05
ap.add_argument("--meta", required=True)              # Sample, Hurley Stage
ap.add_argument("--out", default=None)
ap.add_argument("--mode", default="positive")         # ion mode
ap.add_argument("--ppm", type=float, default=5)
ap.add_argument("--pcut", type=float, default=0.05)
ap.add_argument("--nperm_mcg", type=int, default=1000)
ap.add_argument("--nperm_class", type=int, default=10000)
ap.add_argument("--replot", action="store_true")   # reuse existing mummichog runs
args = ap.parse_args()

OUT = args.out or f"Pathway_{args.sample}"
for sub in ("inputs", "mummichog_runs", "tables", "figures", "scripts"):
    os.makedirs(os.path.join(OUT, sub), exist_ok=True)
shutil.copy(__file__, os.path.join(OUT, "scripts", os.path.basename(__file__)))
SEED, MIN_CLASS = 2026, 5
plt.rcParams.update({"font.family": "sans-serif", "font.size": 9.5})

meta = pd.read_csv(args.meta, encoding="utf-8-sig")
meta.columns = ["Sample", "Stage"]
n_per = meta["Stage"].value_counts().sort_index().to_dict()
n_all = len(meta)

trend = pd.read_csv(f"{args.tables}/T04_trend_with_Hurley_stage.csv")
COMP = {"trend": (f"Stage trend (Hurley 1 \u2192 2 \u2192 3, n = {n_all})",
                  trend["mz"], trend["P_limma"], trend["t"], trend["ID"])}
for key, (a, b) in {"H3vsH2": (3, 2), "H3vsH1": (3, 1), "H2vsH1": (2, 1)}.items():
    f = f"{args.tables}/T05_limma_S{a}vsS{b}.csv"
    if os.path.exists(f):
        d = pd.read_csv(f)
        COMP[key] = (f"Hurley {a} vs Hurley {b} (n = {n_per.get(a, 0)} vs {n_per.get(b, 0)})",
                     d["mz"], d["P.Value"], d["t"], d["ID"])

# ---------------- Part A: mummichog ------------------------------------------
def run_mummichog(key, title, mz, p, t, ids):
    pd.DataFrame({"m.z": mz, "p.value": p, "t.score": t}).to_csv(
        f"{OUT}/inputs/{key}_MetaboAnalyst_input.txt", sep="\t", index=False)
    wd = os.path.join(OUT, "mummichog_runs", key)
    if not (args.replot and glob.glob(f"{wd}/*{key}")):
      if os.path.isdir(wd):
          shutil.rmtree(wd)
      os.makedirs(wd)
      pd.DataFrame({"m/z": mz, "retention_time": 100, "p-value": p, "t-score": t,
                    "custom_id": ids}).to_csv(f"{wd}/input.txt", sep="\t", index=False)
      # mummichog's HTML step can crash when nothing is significant; tables are
      # written before that, so read them directly
      subprocess.run(["mummichog", "-f", "input.txt", "-o", key, "-m", args.mode,
                    "-u", str(args.ppm), "-c", str(args.pcut), "-p", str(args.nperm_mcg)],
                   cwd=wd, check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    rd = sorted(glob.glob(f"{wd}/*{key}"))[-1] + "/tables"
    pw = pd.read_csv(f"{rd}/mcg_pathwayanalysis_{key}.tsv", sep="\t")
    pw.columns = [c.strip() for c in pw.columns]
    pw = pw[pw["pathway_size"] > 0].copy()
    ui = pd.read_csv(f"{rd}/userInput_to_EmpiricalCompounds.tsv", sep="\t")
    sig = ui[ui["p_value"] < args.pcut]
    n_ref = len(pd.read_csv(f"{rd}/ListOfEmpiricalCompounds.tsv", sep="\t"))
    pw["expected"] = pw["pathway_size"] * sig["EID"].nunique() / max(n_ref, 1)
    pw["enrichment_factor"] = pw["overlap_size"] / pw["expected"].replace(0, np.nan)
    pw.insert(0, "comparison", title)
    pw.sort_values("p-value").to_csv(f"{OUT}/tables/pathways_{key}.csv", index=False)
    sig[["m/z", "p_value", "statistic", "EID", "str_row_ion", "compound_names"]].to_csv(
        f"{OUT}/tables/matched_significant_features_{key}.csv", index=False)
    info = dict(comparison=title, features=len(mz), significant=int((p < args.pcut).sum()),
                matched=int(ui["m/z"].nunique()), significant_matched=int(sig["m/z"].nunique()),
                pathways_with_hit=int((pw["overlap_size"] > 0).sum()),
                pathways_p05=int(((pw["overlap_size"] > 0) & (pw["p-value"] < 0.05)).sum()))
    print(json.dumps(info))
    return pw, info

def bubble_plot(ax, pw, title, info):
    hit = pw[pw["overlap_size"] > 0].copy()
    n_non = int((pw["overlap_size"] == 0).sum())
    ymax = max(5, (-np.log10(hit["p-value"])).max() * 1.15 if len(hit) else 5)
    xmax = max(1, hit["enrichment_factor"].max() * 1.9 if len(hit) else 1)
    ax.axhline(-np.log10(0.05), ls="--", lw=0.6, c="grey")
    ax.text(xmax, -np.log10(0.05), " p = 0.05", va="bottom", ha="right", fontsize=7, color="grey")
    if len(hit):
        sc = ax.scatter(hit["enrichment_factor"], -np.log10(hit["p-value"]),
                        s=150 + 120 * hit["overlap_size"], c=-np.log10(hit["p-value"]),
                        cmap="YlOrRd", vmin=0, vmax=ymax, edgecolor="black", linewidth=0.5, zorder=3)
        # pathways sharing the same point get one combined label
        lab = hit[hit["p-value"] < 0.05].copy()
        lab["key"] = list(zip(lab["enrichment_factor"].round(4), lab["p-value"].round(6)))
        for i, (k, g) in enumerate(lab.groupby("key", sort=False)):
            txt = "\n".join(f"{textwrap.shorten(r['pathway'], 40)} "
                            f"({int(r['overlap_size'])}/{int(r['pathway_size'])})"
                            for _, r in g.iterrows())
            ax.annotate(txt, (k[0], -np.log10(k[1])), xytext=(14, 10 + 12 * (i % 2)),
                        textcoords="offset points", ha="left", va="bottom", fontsize=6.8,
                        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="grey", lw=0.4, alpha=0.9),
                        arrowprops=dict(arrowstyle="-", lw=0.4, color="grey"))
        cb = plt.colorbar(sc, ax=ax, fraction=0.04, pad=0.03)
        cb.set_label("-log10(p)")
    else:
        ax.text(0.5, 0.55, "No pathway contained a\nsignificant matched compound",
                transform=ax.transAxes, ha="center", va="center", fontsize=10, color="grey")
    ax.set_xlim(-0.05 * xmax, xmax); ax.set_ylim(0, ymax)
    ax.set_xlabel("Enrichment factor (hits / expected)")
    ax.set_ylabel("-log10(p)  (mummichog, gamma-adjusted)")
    ax.set_title(title, fontweight="bold", loc="left", fontsize=10)
    ax.text(0.98, 0.97, f"{info['matched']} of {info['features']} features matched a compound\n"
            f"{info['significant_matched']} of {info['significant']} significant features matched\n"
            f"{n_non} further pathways: no significant hit",
            transform=ax.transAxes, ha="right", va="top", fontsize=7.2, color="grey")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)

res = {k: run_mummichog(k, *v) for k, v in COMP.items()}
pd.concat([r[0] for r in res.values()]).to_csv(f"{OUT}/tables/pathways_all_comparisons.csv", index=False)
pd.DataFrame([r[1] for r in res.values()]).to_csv(f"{OUT}/tables/summary_by_comparison.csv", index=False)

n = len(res); ncol = 2; nrow = int(np.ceil(n / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(13, 5.2 * nrow))
axes = np.atleast_1d(axes).ravel()
for ax, (k, (pw, info)) in zip(axes, res.items()):
    bubble_plot(ax, pw, COMP[k][0], info)
for ax in axes[n:]:
    ax.axis("off")
fig.suptitle(f"{args.sample}: pathway enrichment (mummichog, human MFN, {args.mode} mode, "
             f"{args.ppm:g} ppm, p < {args.pcut})", fontweight="bold", x=0.01, ha="left")
fig.tight_layout()
for ext in ("pdf", "png", "tiff"):
    fig.savefig(f"{OUT}/figures/pathway_bubble_mummichog.{ext}", dpi=600)
plt.close(fig)

# ---------------- Part B: lipid-class enrichment ------------------------------
d = pd.read_csv(f"{args.tables}/T02_log2_PQN_data.csv")
samples = [s for s in meta["Sample"] if s in d.columns]
y = meta.set_index("Sample").loc[samples, "Stage"].values.astype(float)
X = d[samples].values
Xc = X - X.mean(1, keepdims=True)
Xs = Xc / np.sqrt((Xc ** 2).sum(1, keepdims=True))
cls = d["Class"].astype(str).copy()
ether = d["Name"].astype(str).str.match(r"^(PC|PE|PS)[ (]?\(?[OP]-")
cls[ether] = cls[ether] + "-O"
classes = [c for c, k in cls.value_counts().items()
           if k >= MIN_CLASS and c not in ("Unannotated", "Other annotated", "nan")]
masks = {c: (cls == c).values for c in classes}

def class_stats(yv):
    yc = yv - yv.mean()
    r = Xs @ (yc / np.sqrt((yc ** 2).sum()))
    return {c: r[m].mean() - r[~m].mean() for c, m in masks.items()}, r

obs, r = class_stats(y)
rng = np.random.default_rng(SEED)
null = {c: np.empty(args.nperm_class) for c in classes}
for b in range(args.nperm_class):
    s, _ = class_stats(rng.permutation(y))
    for c in classes:
        null[c][b] = s[c]
ce = pd.DataFrame([{"Class": c, "n_features": int(masks[c].sum()),
                    "Mean_r_with_stage": r[masks[c]].mean(),
                    "Shift_vs_other_features": obs[c],
                    "P_permutation": (np.sum(np.abs(null[c]) >= abs(obs[c])) + 1) / (args.nperm_class + 1)}
                   for c in classes]).sort_values("P_permutation").reset_index(drop=True)
m = len(ce)
ce["FDR"] = np.minimum(1, np.minimum.accumulate((ce["P_permutation"] * m / np.arange(1, m + 1))[::-1])[::-1])
ce.to_csv(f"{OUT}/tables/lipid_class_enrichment.csv", index=False)
print(ce.round(4).to_string(index=False))

ce_p = ce.sort_values("Shift_vs_other_features")
fig, ax = plt.subplots(figsize=(6.4, 0.32 * len(ce_p) + 1.4))
ax.axvline(0, color="grey", lw=0.6)
ax.scatter(ce_p["Shift_vs_other_features"], range(len(ce_p)), s=20 + 5 * ce_p["n_features"],
           c=["#B2182B" if v > 0 else "#2166AC" for v in ce_p["Shift_vs_other_features"]],
           edgecolor="black", linewidth=0.4, zorder=3)
ax.set_yticks(range(len(ce_p)))
ax.set_yticklabels([f"{c} (n = {k})  p = {p:.3f}" for c, k, p in
                    zip(ce_p["Class"], ce_p["n_features"], ce_p["P_permutation"])])
ax.set_xlabel("Shift in correlation with Hurley stage\n(class mean r minus all other features)")
ax.set_title(f"{args.sample}: lipid-class enrichment along Hurley stage", fontweight="bold", loc="left")
for s_ in ("top", "right"):
    ax.spines[s_].set_visible(False)
fig.tight_layout()
for ext in ("pdf", "png", "tiff"):
    fig.savefig(f"{OUT}/figures/lipid_class_enrichment.{ext}", dpi=600)
print("Done:", OUT)
