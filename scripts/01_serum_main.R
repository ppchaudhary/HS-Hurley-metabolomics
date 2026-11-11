###############################################################################
#  SERUM metabolomics in Hidradenitis Suppurativa (HS) across Hurley stages
#  ---------------------------------------------------------------------------
#  Input : Serum_data_HS.csv (col1 = m/z, col2 = Name, then one column/sample)
#          Metadata_serum.csv (Sample, Hurley Stage)
#  Output: ./HS_serum_results/  figures (PDF + 600 dpi PNG + TIFF) and CSV tables
#
#  Workflow
#   1  Import, unique feature labels (Name > Name (m/z) if duplicated > m/z)
#   2  QC of raw data (TIC, sample correlation, raw PCA)
#   3  Filtering, half-minimum imputation, PQN normalisation, log2, Pareto
#   4  Ordination: PCA, NMDS (Bray-Curtis) + PERMANOVA + dispersion test
#   5  Ordinal trend with Hurley stage (limma + Spearman + Jonckheere-Terpstra)
#   6  Pairwise limma contrasts + EnhancedVolcano
#   7  Heatmaps (top trend features; stage-mean heatmap)
#   8  Trend box plots of top metabolites
#   9  Lipid-class level analysis
#  10  Random Forest (Stage 2 vs 3 classification; stage regression)
#  11  sPLS-DA (mixOmics)
#  12  Sensitivity analyses (without QC-flagged samples, chosen automatically)
###############################################################################

## ======================= 0. USER SETTINGS ==================================
DATA_FILE <- "data/Serum_data_HS.csv"
META_FILE <- "data/Metadata_serum.csv"
OUT_DIR   <- "results/serum"

DETECT_FRAC <- 0.80      # feature kept if >0 in >= 80% of samples of >= 1 stage
FC_CUT      <- 1         # |log2FC| threshold for volcano
P_CUT       <- 0.05      # p threshold for volcano / trend calls
TOP_N_HEAT  <- 50        # features in heatmap
TOP_N_BOX   <- 12        # features in trend box plots
TOP_N_RF    <- 25        # features in RF importance plots
N_PERM_RF   <- 200       # label permutations for RF significance
RUN_BORUTA  <- TRUE      # Boruta feature selection (slower)
# QC flags for sensitivity analyses. "auto" = chosen from the data:
#   SUSPECT_SAMPLES: median Spearman correlation to other samples < median - 3*MAD
#   LOW_SIGNAL     : total ion signal < 1/4 of the median
# or give sample names, e.g. c("51-HS", "11-HS")
SUSPECT_SAMPLES <- "auto"
LOW_SIGNAL      <- "auto"
COLLAPSE_DUP_MZ <- TRUE  # average rows sharing the same m/z (near-identical duplicate peaks)
FONT <- "Arial"
SEED <- 2026

## ======================= 0b. PACKAGES ======================================
cran <- c("ggplot2", "ggrepel", "patchwork", "dplyr", "tidyr", "vegan",
          "randomForest", "Boruta", "clinfun", "circlize", "ragg",
          "systemfonts", "scales")
bioc <- c("limma", "EnhancedVolcano", "ComplexHeatmap", "mixOmics")
miss_c <- cran[!cran %in% rownames(installed.packages())]
if (length(miss_c)) install.packages(miss_c, repos = "https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
miss_b <- bioc[!bioc %in% rownames(installed.packages())]
if (length(miss_b)) BiocManager::install(miss_b, update = FALSE)

suppressPackageStartupMessages({
  library(ggplot2); library(ggrepel); library(patchwork); library(dplyr)
  library(tidyr);   library(vegan);   library(limma);     library(EnhancedVolcano)
  library(ComplexHeatmap); library(circlize); library(randomForest)
  library(mixOmics); library(clinfun); library(scales)
})

set.seed(SEED)
dir.create(OUT_DIR, showWarnings = FALSE)
FIG <- file.path(OUT_DIR, "figures"); dir.create(FIG, showWarnings = FALSE)
TAB <- file.path(OUT_DIR, "tables");  dir.create(TAB, showWarnings = FALSE)

# font fallback (Arial may not be installed on Linux servers)
if (!FONT %in% systemfonts::system_fonts()$family) FONT <- "sans"

## ======================= 0c. STYLE =========================================
# Sequential, colour-blind-safe: light -> dark with increasing severity
STAGE_COL <- c("Hurley 1" = "#FDAE6B", "Hurley 2" = "#E6550D", "Hurley 3" = "#7F2704")
DIR_COL   <- c("Up with stage" = "#B2182B", "Down with stage" = "#2166AC")

theme_pub <- function(base = 11) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(axis.text        = element_text(colour = "black"),
          axis.line        = element_line(linewidth = 0.5),
          axis.ticks       = element_line(linewidth = 0.5, colour = "black"),
          legend.title     = element_text(face = "bold"),
          legend.key.size  = unit(0.45, "cm"),
          plot.title       = element_text(face = "bold", size = base + 2),
          plot.subtitle    = element_text(size = base - 1, colour = "grey30"),
          strip.background = element_blank(),
          strip.text       = element_text(face = "bold"))
}
theme_set(theme_pub())

save_gg <- function(p, name, w = 7, h = 5) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h,
         device = cairo_pdf)
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h,
         dpi = 600, device = ragg::agg_png, bg = "white")
  ggsave(file.path(FIG, paste0(name, ".tiff")), p, width = w, height = h,
         dpi = 600, device = ragg::agg_tiff, compression = "lzw", bg = "white")
}
save_ht <- function(ht, name, w = 8, h = 10) {
  cairo_pdf(file.path(FIG, paste0(name, ".pdf")), width = w, height = h,
            family = FONT)
  draw(ht, merge_legend = TRUE); dev.off()
  ragg::agg_png(file.path(FIG, paste0(name, ".png")), width = w, height = h,
                units = "in", res = 600, background = "white")
  draw(ht, merge_legend = TRUE); dev.off()
  ragg::agg_tiff(file.path(FIG, paste0(name, ".tiff")), width = w, height = h,
                 units = "in", res = 600, background = "white", compression = "lzw")
  draw(ht, merge_legend = TRUE); dev.off()
}
trunc_lab <- function(x, n = 45)
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "\u2026"), x)
fmt_p <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
                            formatC(p, format = "f", digits = 3))

## ======================= 1. IMPORT ==========================================
raw  <- read.csv(DATA_FILE, check.names = FALSE, fileEncoding = "UTF-8-BOM",
                 stringsAsFactors = FALSE)
meta <- read.csv(META_FILE, check.names = FALSE, fileEncoding = "UTF-8-BOM",
                 stringsAsFactors = FALSE)
colnames(meta) <- c("Sample", "Stage_num")
meta$Sample    <- trimws(meta$Sample)
meta$Stage     <- factor(paste("Hurley", meta$Stage_num),
                         levels = c("Hurley 1", "Hurley 2", "Hurley 3"))
meta$num       <- suppressWarnings(as.numeric(sub("-.*", "", meta$Sample)))
if (anyNA(meta$num)) meta$num <- seq_len(nrow(meta))
meta           <- meta[order(meta$Stage_num, meta$num), ]
rownames(meta) <- meta$Sample

stopifnot(all(meta$Sample %in% colnames(raw)))
mz   <- as.numeric(raw[[1]])
nm   <- trimws(as.character(raw[[2]])); nm[nm == "" | nm == "NA"] <- NA
X    <- as.matrix(raw[, meta$Sample]); storage.mode(X) <- "numeric"

# --- collapse duplicated m/z rows (the export contains ~780 m/z values listed
#     twice with near-identical intensities; keeping both double-counts them)
if (COLLAPSE_DUP_MZ && any(duplicated(mz))) {
  key <- sprintf("%.4f", mz)
  cat(sprintf("Collapsing %d duplicated m/z rows by averaging\n", sum(duplicated(key))))
  grp <- factor(key, levels = unique(key))
  X   <- rowsum(X, grp, reorder = FALSE) / as.vector(table(grp)[levels(grp)])
  nm  <- tapply(nm, grp, function(v) {
    v <- unique(v[!is.na(v)]); if (length(v)) paste(v, collapse = " / ") else NA })[levels(grp)]
  mz  <- as.numeric(levels(grp))
  nm  <- as.character(nm)
}

# --- feature labels: Name if unique, "Name (m/z)" if duplicated, else "m/z x"
dup_nm <- !is.na(nm) & (duplicated(nm) | duplicated(nm, fromLast = TRUE))
lab <- ifelse(is.na(nm), paste0("m/z ", sprintf("%.4f", mz)),
              ifelse(dup_nm, paste0(nm, " (", sprintf("%.4f", mz), ")"), nm))
lab <- make.unique(lab, sep = " _")

# --- lipid class from annotation
lipid_class <- function(n) {
  if (is.na(n)) return("Unannotated")
  pats <- c("PE-Cer" = "^PE-Cer", "HexCer" = "^HexCer", "LPC" = "^(LPC|LysoPC)",
            "LPE" = "^(LPE|LysoPE)", "PC" = "^PC[ (]", "PE" = "^PE[ (]",
            "PS" = "^PS[ (]", "PG" = "^PG[ (]", "PI" = "^PI[ (]", "PA" = "^PA[ (]",
            "SM" = "^SM[ (]", "Cer" = "^Cer[ (]", "TG" = "^TG[ (]",
            "DG" = "^DG[ (]", "MG" = "^MG[ (]", "CE" = "^CE[ (]",
            "FA" = "^FA[ (]", "CAR" = "^CAR[ (]",
            "Peptide" = "Tryptic Peptide")
  hit <- names(pats)[vapply(pats, grepl, logical(1), x = n)]
  if (length(hit)) hit[1] else "Other annotated"
}
annot <- data.frame(ID = sprintf("F%04d", seq_along(mz)), mz = mz, Name = nm,
                    Label = lab, Class = vapply(nm, lipid_class, character(1)),
                    stringsAsFactors = FALSE)
rownames(annot) <- annot$ID
rownames(X) <- annot$ID
cat(sprintf("Features: %d | Named: %d | Samples: %d\n", nrow(X), sum(!is.na(nm)), ncol(X)))
print(table(meta$Stage))

## ======================= 2. RAW QC ==========================================
qc <- data.frame(Sample = colnames(X), TIC = colSums(X),
                 Zeros = colSums(X == 0), Stage = meta$Stage)
qc$Sample <- factor(qc$Sample, levels = meta$Sample)
p_tic <- ggplot(qc, aes(Sample, TIC / 1e6, fill = Stage)) +
  geom_col(width = 0.75, colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = STAGE_COL, limits = names(STAGE_COL)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = expression("Total ion signal ("*10^6*")"),
       title = "Total ion signal per sample (raw)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_gg(p_tic, "01a_QC_total_ion_signal", 7, 4)

# sample-sample Spearman correlation on log raw data
cor_raw <- cor(log2(X + 1), method = "spearman")
ha_col <- HeatmapAnnotation(Stage = meta$Stage, col = list(Stage = STAGE_COL),
                            annotation_name_gp = gpar(fontsize = 9, fontface = "bold"))
ht_cor <- Heatmap(cor_raw, name = "Spearman\nrho",
                  col = colorRamp2(c(0.3, 0.65, 1), c("#2166AC", "white", "#B2182B")),
                  top_annotation = ha_col, cluster_rows = TRUE, cluster_columns = TRUE,
                  row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 9),
                  cell_fun = if (ncol(X) <= 20) function(j, i, x, y, w, h, fill)
                    grid.text(sprintf("%.2f", cor_raw[i, j]), x, y,
                              gp = gpar(fontsize = 6)) else NULL,
                  column_title = "Sample correlation (raw, log2)",
                  column_title_gp = gpar(fontface = "bold"))
save_ht(ht_cor, "01b_QC_sample_correlation", max(8, ncol(X) * 0.18 + 2), max(7, ncol(X) * 0.18 + 1))

# ---- automatic QC flags --------------------------------------------------
cc <- cor_raw; diag(cc) <- NA
med_cor <- apply(cc, 2, median, na.rm = TRUE)
if (identical(SUSPECT_SAMPLES, "auto"))
  SUSPECT_SAMPLES <- names(which(med_cor < median(med_cor) - 3 * mad(med_cor)))
if (identical(LOW_SIGNAL, "auto"))
  LOW_SIGNAL <- as.character(qc$Sample[qc$TIC < median(qc$TIC) / 4])
qc$Median_correlation <- med_cor[as.character(qc$Sample)]
qc$Flag <- ifelse(qc$Sample %in% SUSPECT_SAMPLES & qc$Sample %in% LOW_SIGNAL, "Low correlation + low signal",
           ifelse(qc$Sample %in% SUSPECT_SAMPLES, "Low correlation",
           ifelse(qc$Sample %in% LOW_SIGNAL, "Low signal", "")))
write.csv(qc, file.path(TAB, "T00_QC_per_sample.csv"), row.names = FALSE)
cat("QC flags - low correlation:", paste(SUSPECT_SAMPLES, collapse = ", "),
    "| low signal:", paste(LOW_SIGNAL, collapse = ", "), "\n")

## ======================= 3. PRE-PROCESSING =================================
# 3a. detection filter
det <- sapply(levels(meta$Stage), function(s)
  rowMeans(X[, meta$Sample[meta$Stage == s], drop = FALSE] > 0))
keep <- apply(det, 1, max) >= DETECT_FRAC
Xf <- X[keep, ]
cat(sprintf("Features after detection filter: %d\n", nrow(Xf)))

# 3b. half-minimum imputation of zeros (per feature)
Xi <- t(apply(Xf, 1, function(v) {
  if (any(v == 0)) v[v == 0] <- min(v[v > 0]) / 2; v }))

# 3c. PQN normalisation (Dieterle et al. 2006)
pqn <- function(M) {
  Mt   <- sweep(M, 2, colSums(M), "/") * median(colSums(M))   # integral norm
  ref  <- apply(Mt, 1, median)                                  # reference
  dil  <- apply(Mt / ref, 2, median)                            # dilution factors
  list(norm = sweep(Mt, 2, dil, "/"), factor = dil)
}
pq   <- pqn(Xi)
Xn   <- pq$norm                         # PQN-normalised intensities
L    <- log2(Xn)                        # log2 PQN (statistics)
pareto <- function(M) {                 # features in rows
  cen <- M - rowMeans(M); cen / sqrt(apply(M, 1, sd)) }
P    <- pareto(L)                       # log2 + Pareto (multivariate)

write.csv(data.frame(annot[rownames(Xn), ], Xn, check.names = FALSE),
          file.path(TAB, "T01_PQN_normalised_intensities.csv"), row.names = FALSE)
write.csv(data.frame(annot[rownames(L), ], L, check.names = FALSE),
          file.path(TAB, "T02_log2_PQN_data.csv"), row.names = FALSE)

# before / after boxplots
long_box <- function(M, what) {
  data.frame(M, check.names = FALSE) |>
    tidyr::pivot_longer(everything(), names_to = "Sample", values_to = "v") |>
    mutate(Sample = factor(Sample, levels = meta$Sample),
           Stage = meta[as.character(Sample), "Stage"], Data = what)
}
bx <- rbind(long_box(log2(Xi), "Raw (log2)"), long_box(L, "PQN normalised (log2)"))
bx$Data <- factor(bx$Data, levels = c("Raw (log2)", "PQN normalised (log2)"))
p_norm <- ggplot(bx, aes(Sample, v, fill = Stage)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.3, linewidth = 0.3) +
  facet_wrap(~Data, ncol = 1) + scale_fill_manual(values = STAGE_COL, limits = names(STAGE_COL)) +
  labs(x = NULL, y = expression(log[2]~"intensity"),
       title = "Intensity distributions before and after normalisation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_gg(p_norm, "02_normalisation_boxplots", 7.5, 6.5)

## ======================= 4. ORDINATION =====================================
nmds_plot <- function(Xmat, md, title) {
  # Bray-Curtis on log-transformed PQN intensities (limits dominance of the
  # few very abundant ions)
  Dm <- vegdist(t(log2(Xmat + 1)), method = "bray")
  set.seed(SEED)
  # "hybrid" NMDS (linear + monotone): the default global NMDS degenerates
  # (stress ~ 0) here because a few divergent samples dominate the rank order
  fit <- metaMDS(Dm, k = 2, trymax = 200, trace = FALSE, model = "hybrid")
  pf  <- adonis2(Dm ~ Stage, data = md, permutations = 9999)
  pt  <- adonis2(Dm ~ Stage_num, data = md, permutations = 9999)
  bdp <- permutest(betadisper(Dm, md$Stage), permutations = 9999)
  df  <- data.frame(scores(fit, display = "sites"), md)
  ok  <- names(which(table(md$Stage) >= 3))
  sub <- sprintf(paste0("Hybrid NMDS stress = %.3f | PERMANOVA stage R\u00b2 = %.2f, p = %s",
                        "\nTrend R\u00b2 = %.2f, p = %s | Dispersion p = %s"),
                 fit$stress, pf$R2[1], fmt_p(pf$`Pr(>F)`[1]),
                 pt$R2[1], fmt_p(pt$`Pr(>F)`[1]), fmt_p(bdp$tab$`Pr(>F)`[1]))
  p <- ggplot(df, aes(NMDS1, NMDS2, fill = Stage, colour = Stage)) +
    stat_ellipse(data = subset(df, Stage %in% ok), geom = "polygon",
                 alpha = 0.12, level = 0.95, type = "norm", linewidth = 0.4,
                 show.legend = FALSE) +
    geom_point(size = 3.2, shape = 21, colour = "black", stroke = 0.4) +
    geom_text_repel(aes(label = Sample), size = 3, family = FONT,
                    colour = "black", seed = SEED, max.overlaps = 30) +
    scale_fill_manual(values = STAGE_COL, limits = names(STAGE_COL), drop = FALSE) +
    scale_colour_manual(values = STAGE_COL, limits = names(STAGE_COL)) +
    labs(title = title, subtitle = sub) +
    theme(plot.subtitle = element_text(size = 8.5))
  list(plot = p, D = Dm, fit = fit, pf = pf, pt = pt, bd = bdp)
}

ell_ok <- names(which(table(meta$Stage) >= 3))   # ellipse only if n >= 3

# 4a. PCA
pca <- prcomp(t(P), center = FALSE, scale. = FALSE)
ve  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pc  <- data.frame(pca$x[, 1:2], meta)
p_pca <- ggplot(pc, aes(PC1, PC2, colour = Stage, fill = Stage)) +
  stat_ellipse(data = subset(pc, Stage %in% ell_ok), geom = "polygon",
               alpha = 0.12, level = 0.95, type = "norm", linewidth = 0.4, show.legend = FALSE) +
  geom_point(size = 3.2, shape = 21, colour = "black", stroke = 0.4) +
  geom_text_repel(aes(label = Sample), size = 3, family = FONT,
                  colour = "black", seed = SEED) +
  scale_colour_manual(values = STAGE_COL, limits = names(STAGE_COL)) + scale_fill_manual(values = STAGE_COL, limits = names(STAGE_COL)) +
  labs(x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)"),
       title = "PCA (log2 PQN, Pareto scaled)")
save_gg(p_pca, "03a_PCA", 6.5, 5)

# 4b. NMDS (Bray-Curtis) - all samples
nm_all <- nmds_plot(Xn, meta, "NMDS (Bray-Curtis), all samples")
save_gg(nm_all$plot, "03b_NMDS_BrayCurtis", 6.8, 5.5)
nmds <- nm_all$fit; perm_f <- nm_all$pf; perm_t <- nm_all$pt; bd <- nm_all$bd

# 4c. NMDS without the low-correlation (technically divergent) samples; re-normalised
keep_s <- setdiff(meta$Sample, SUSPECT_SAMPLES)
nm_sub <- nmds_plot(pqn(Xi[, keep_s])$norm, meta[keep_s, ],
                    sprintf("NMDS without low-correlation samples (n = %d)", length(SUSPECT_SAMPLES)))
save_gg(nm_all$plot + nm_sub$plot + plot_layout(guides = "collect"),
        "03c_NMDS_all_vs_without_flagged", 13, 5.5)

ord_stats <- data.frame(
  Test = c("PERMANOVA stage (factor)", "PERMANOVA stage (ordinal trend)",
           "Beta-dispersion (PERMDISP)", "NMDS stress"),
  R2 = c(perm_f$R2[1], perm_t$R2[1], NA, NA),
  F  = c(perm_f$F[1], perm_t$F[1], bd$tab$F[1], NA),
  P  = c(perm_f$`Pr(>F)`[1], perm_t$`Pr(>F)`[1], bd$tab$`Pr(>F)`[1], NA),
  Value = c(NA, NA, NA, nmds$stress))
write.csv(ord_stats, file.path(TAB, "T03_ordination_PERMANOVA.csv"), row.names = FALSE)

## ======================= 5. ORDINAL TREND ANALYSIS =========================
run_trend <- function(Lm, md) {
  des <- model.matrix(~ md$Stage_num)
  fit <- eBayes(lmFit(Lm, des))
  tt  <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  sp  <- t(apply(Lm, 1, function(v) {
    ct <- suppressWarnings(cor.test(v, md$Stage_num, method = "spearman", exact = FALSE))
    c(ct$estimate, ct$p.value) }))
  jt  <- apply(Lm, 1, function(v)
    suppressWarnings(jonckheere.test(v, md$Stage_num, alternative = "two.sided")$p.value))
  data.frame(annot[rownames(Lm), ],
             log2FC_per_stage = tt$logFC, AveExpr = tt$AveExpr, t = tt$t,
             P_limma = tt$P.Value, FDR_limma = p.adjust(tt$P.Value, "BH"),
             Spearman_rho = sp[, 1], P_Spearman = sp[, 2],
             FDR_Spearman = p.adjust(sp[, 2], "BH"),
             P_JT = jt, FDR_JT = p.adjust(jt, "BH"),
             Direction = ifelse(tt$logFC > 0, "Up with stage", "Down with stage"),
             stringsAsFactors = FALSE)
}
trend <- run_trend(L, meta)
trend <- trend[order(trend$P_limma), ]
trend$Consensus <- trend$P_limma < P_CUT & trend$P_Spearman < P_CUT & trend$P_JT < P_CUT
write.csv(trend, file.path(TAB, "T04_trend_with_Hurley_stage.csv"), row.names = FALSE)
cat(sprintf("Trend: limma p<%.2f = %d | FDR<0.05 = %d | consensus (3 tests) = %d\n",
            P_CUT, sum(trend$P_limma < P_CUT), sum(trend$FDR_limma < 0.05),
            sum(trend$Consensus)))

# trend "volcano": slope per stage vs significance
top_tr <- head(trend, 20)
p_trv <- ggplot(trend, aes(log2FC_per_stage, -log10(P_limma))) +
  geom_point(aes(colour = ifelse(P_limma < P_CUT, Direction, "NS")),
             size = 1.2, alpha = 0.75) +
  scale_colour_manual(values = c(DIR_COL, NS = "grey75"), name = NULL) +
  geom_hline(yintercept = -log10(P_CUT), linetype = 2, linewidth = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_text_repel(data = top_tr, aes(label = trunc_lab(Label, 35)), size = 2.6,
                  family = FONT, max.overlaps = 50, segment.size = 0.25,
                  min.segment.length = 0, seed = SEED) +
  labs(x = expression(log[2]~"fold change per Hurley stage"),
       y = expression(-log[10]~italic(P)~"(limma trend)"),
       title = "Metabolites associated with increasing Hurley stage",
       subtitle = sprintf("%d features p < %.2f; %d FDR < 0.05",
                          sum(trend$P_limma < P_CUT), P_CUT,
                          sum(trend$FDR_limma < 0.05)))
save_gg(p_trv, "04_trend_volcano_stage", 7.5, 6)

## ======================= 6. PAIRWISE CONTRASTS + ENHANCEDVOLCANO ============
des_g <- model.matrix(~ 0 + Stage, data = meta)
colnames(des_g) <- c("S1", "S2", "S3")
fit_g <- lmFit(L, des_g)
cm <- makeContrasts(S3vsS2 = S3 - S2, S2vsS1 = S2 - S1, S3vsS1 = S3 - S1,
                    levels = des_g)
fit_c <- eBayes(contrasts.fit(fit_g, cm))

volc <- function(coef, title, note = NULL) {
  tt <- topTable(fit_c, coef = coef, number = Inf, sort.by = "none")
  tt <- data.frame(annot[rownames(tt), ], tt)
  tt$FDR <- p.adjust(tt$P.Value, "BH")
  write.csv(tt[order(tt$P.Value), ],
            file.path(TAB, paste0("T05_limma_", coef, ".csv")), row.names = FALSE)
  sig <- tt[tt$P.Value < P_CUT & abs(tt$logFC) > FC_CUT, ]
  up  <- head(sig[sig$logFC > 0, ][order(sig$P.Value[sig$logFC > 0]), ], 10)
  dn  <- head(sig[sig$logFC < 0, ][order(sig$P.Value[sig$logFC < 0]), ], 10)
  sel_ids <- c(up$ID, dn$ID)
  tt$Plab <- ifelse(tt$ID %in% sel_ids, trunc_lab(tt$Label, 35), "")
  nsig <- c(up = sum(sig$logFC > 0), down = sum(sig$logFC < 0),
            fdr = sum(tt$FDR < 0.05 & abs(tt$logFC) > FC_CUT))
  EnhancedVolcano(tt, lab = tt$Plab, x = "logFC", y = "P.Value",
                  selectLab = tt$Plab[tt$ID %in% sel_ids],
                  ylim = c(0, max(-log10(tt$P.Value)) * 1.15),
                  pCutoff = P_CUT, FCcutoff = FC_CUT,
                  pointSize = 1.6, labSize = 2.8, labFace = "plain",
                  col = c("grey75", "#7FA7C9", "#9E9AC8", "#B2182B"),
                  colAlpha = 0.8, drawConnectors = TRUE, widthConnectors = 0.3,
                  colConnectors = "grey40", max.overlaps = 60,
                  boxedLabels = FALSE, legendPosition = "bottom",
                  legendLabSize = 9, legendIconSize = 3,
                  legendLabels = c("NS", expression(log[2]~FC),
                                   "p-value", expression(p~and~log[2]~FC)),
                  title = title,
                  subtitle = paste0(sprintf("Up: %d | Down: %d (p < %.2f, |log2FC| > %g) | FDR < 0.05: %d",
                                            nsig["up"], nsig["down"], P_CUT, FC_CUT, nsig["fdr"]),
                                    if (!is.null(note)) paste0("\n", note) else ""),
                  caption = paste0("limma moderated t-test; total features = ", nrow(tt)),
                  gridlines.major = FALSE, gridlines.minor = FALSE,
                  axisLabSize = 12, titleLabSize = 13, subtitleLabSize = 9,
                  captionLabSize = 8) +
    theme(text = element_text(family = FONT))
}
v32 <- volc("S3vsS2", "Hurley 3 vs Hurley 2")
n_h1 <- sum(meta$Stage == "Hurley 1")
note_h1 <- if (n_h1 < 5) sprintf("Exploratory: Hurley 1 n = %d", n_h1) else NULL
v21 <- volc("S2vsS1", "Hurley 2 vs Hurley 1", note_h1)
v31 <- volc("S3vsS1", "Hurley 3 vs Hurley 1", note_h1)
save_gg(v32, "05a_EnhancedVolcano_H3_vs_H2", 7.5, 7.5)
save_gg(v21, "05b_EnhancedVolcano_H2_vs_H1", 7.5, 7.5)
save_gg(v31, "05c_EnhancedVolcano_H3_vs_H1", 7.5, 7.5)

## ======================= 7. HEATMAPS ========================================
top_ids <- head(trend$ID, TOP_N_HEAT)
Z <- t(scale(t(L[top_ids, ])))
rownames(Z) <- trunc_lab(annot[top_ids, "Label"], 50)
cls <- annot[top_ids, "Class"]
cls_lv <- sort(unique(cls))
cls_col <- setNames(colorRampPalette(c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                                       "#66A61E", "#E6AB02", "#A6761D", "#666666"))(length(cls_lv)),
                    cls_lv)
if ("Unannotated" %in% cls_lv) cls_col["Unannotated"] <- "grey85"
ra <- rowAnnotation(
  Trend = trend[match(top_ids, trend$ID), "Direction"],
  Class = cls,
  `-log10 P` = anno_barplot(-log10(trend[match(top_ids, trend$ID), "P_limma"]),
                            width = unit(1.5, "cm"),
                            gp = gpar(fill = "grey50", col = NA)),
  col = list(Trend = DIR_COL, Class = cls_col),
  annotation_name_gp = gpar(fontsize = 8, fontface = "bold"))
z_col <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
ht_top <- Heatmap(Z, name = "Z-score", col = z_col,
                  top_annotation = ha_col, right_annotation = ra,
                  column_split = meta$Stage, cluster_column_slices = FALSE,
                  column_title_gp = gpar(fontface = "bold", fontsize = 10),
                  cluster_rows = TRUE, cluster_columns = TRUE,
                  clustering_distance_rows = "pearson",
                  row_names_gp = gpar(fontsize = 6.5),
                  column_names_gp = gpar(fontsize = 8),
                  row_names_max_width = unit(8, "cm"),
                  border = TRUE, rect_gp = gpar(col = "white", lwd = 0.3))
save_ht(ht_top, "06a_heatmap_top_trend_features", 11.5, 11)

# stage-mean heatmap
M <- sapply(levels(meta$Stage), function(s)
  rowMeans(L[top_ids, meta$Stage == s, drop = FALSE]))
Zm <- t(scale(t(M))); rownames(Zm) <- rownames(Z)
ha_m <- HeatmapAnnotation(Stage = levels(meta$Stage),
                          col = list(Stage = STAGE_COL), show_legend = FALSE,
                          annotation_name_gp = gpar(fontsize = 8, fontface = "bold"))
ord  <- order(trend[match(top_ids, trend$ID), "log2FC_per_stage"])
ht_m <- Heatmap(Zm[ord, ], name = "Z-score\n(stage mean)", col = z_col,
                top_annotation = ha_m, cluster_rows = FALSE, cluster_columns = FALSE,
                right_annotation = rowAnnotation(
                  Trend = trend[match(top_ids, trend$ID), "Direction"][ord],
                  col = list(Trend = DIR_COL), annotation_name_gp = gpar(fontsize = 8)),
                row_names_gp = gpar(fontsize = 6.5), column_names_gp = gpar(fontsize = 9),
                column_names_rot = 0, column_names_centered = TRUE,
                row_names_max_width = unit(8, "cm"), border = TRUE,
                rect_gp = gpar(col = "white", lwd = 0.5),
                column_title = "Stage-mean profile of top trend metabolites",
                column_title_gp = gpar(fontface = "bold"))
save_ht(ht_m, "06b_heatmap_stage_means", 7.5, 11)

## ======================= 8. TREND BOX PLOTS ================================
box_ids <- head(trend$ID, TOP_N_BOX)
bd_df <- data.frame(t(L[box_ids, ]), check.names = FALSE) |>
  mutate(Sample = rownames(t(L[box_ids, ])), Stage = meta$Stage) |>
  pivot_longer(all_of(box_ids), names_to = "ID", values_to = "log2int") |>
  left_join(trend[, c("ID", "Label", "P_limma", "FDR_limma")], by = "ID") |>
  mutate(Facet = paste0(trunc_lab(Label, 32), "\np = ", fmt_p(P_limma),
                        ", FDR = ", fmt_p(FDR_limma)),
         Facet = factor(Facet, levels = unique(Facet[match(box_ids, ID)])))
p_box <- ggplot(bd_df, aes(Stage, log2int, fill = Stage)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85, linewidth = 0.3) +
  geom_jitter(width = 0.12, size = 1.4, shape = 21, colour = "black",
              stroke = 0.3, show.legend = FALSE) +
  stat_summary(aes(group = 1), fun = median, geom = "line", colour = "black",
               linewidth = 0.4, linetype = 2) +
  facet_wrap(~Facet, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = STAGE_COL, limits = names(STAGE_COL)) +
  scale_x_discrete(labels = c("1", "2", "3")) +
  labs(x = "Hurley stage", y = expression(log[2]~"PQN intensity"),
       title = "Top metabolites associated with Hurley stage") +
  theme(strip.text = element_text(size = 7.5), legend.position = "none")
save_gg(p_box, "07_trend_boxplots_top_metabolites", 10, 2.6 * ceiling(TOP_N_BOX / 4))

## ======================= 9. LIPID-CLASS LEVEL ==============================
cls_all <- annot[rownames(Xn), "Class"]
use_cls <- names(which(table(cls_all[!cls_all %in% c("Unannotated", "Other annotated")]) >= 3))
if (length(use_cls)) {
  CS <- t(sapply(use_cls, function(k) colSums(Xn[cls_all == k, , drop = FALSE])))
  CL <- log2(CS)
  fit_cl <- eBayes(lmFit(CL, model.matrix(~ meta$Stage_num)))
  cl_tab <- topTable(fit_cl, coef = 2, number = Inf, sort.by = "P")
  cl_tab$Class <- rownames(cl_tab)
  cl_tab$n_features <- as.integer(table(cls_all)[cl_tab$Class])
  cl_tab$Spearman_rho <- apply(CL[cl_tab$Class, , drop = FALSE], 1,
                               function(v) cor(v, meta$Stage_num, method = "spearman"))
  write.csv(cl_tab, file.path(TAB, "T06_lipid_class_trend.csv"), row.names = FALSE)
  cl_df <- data.frame(t(CL), check.names = FALSE) |> mutate(Stage = meta$Stage) |>
    pivot_longer(-Stage, names_to = "Class", values_to = "v") |>
    left_join(cl_tab[, c("Class", "P.Value", "n_features")], by = "Class") |>
    mutate(Facet = paste0(Class, " (n = ", n_features, ")\np = ", fmt_p(P.Value)),
           Facet = factor(Facet, levels = unique(Facet[order(P.Value)])))
  p_cls <- ggplot(cl_df, aes(Stage, v, fill = Stage)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85, linewidth = 0.3) +
    geom_jitter(width = 0.12, size = 1.3, shape = 21, colour = "black", stroke = 0.3) +
    facet_wrap(~Facet, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = STAGE_COL, limits = names(STAGE_COL)) + scale_x_discrete(labels = c("1", "2", "3")) +
    labs(x = "Hurley stage", y = expression(log[2]~"summed class intensity"),
         title = "Lipid / compound class totals across Hurley stages") +
    theme(legend.position = "none", strip.text = element_text(size = 8))
  save_gg(p_cls, "08_lipid_class_trends", 9, 2.4 * ceiling(length(use_cls) / 4) + 0.6)
}

## ======================= 10. RANDOM FOREST =================================
Xrf <- t(P)                                   # samples x features
colnames(Xrf) <- rownames(P)
dir_of <- setNames(trend$Direction, trend$ID)
top50  <- head(trend$ID, TOP_N_HEAT)

imp_plot <- function(imp, metric, title, xl) {
  d <- data.frame(ID = rownames(imp), Imp = imp[, metric]) |>
    arrange(desc(Imp)) |> head(TOP_N_RF) |>
    mutate(Label = trunc_lab(annot[ID, "Label"], 40), Direction = dir_of[ID],
           InTrend = ifelse(ID %in% top50, "Yes", "No"),
           Label = factor(Label, levels = rev(Label)))
  ggplot(d, aes(Imp, Label)) +
    geom_segment(aes(x = 0, xend = Imp, yend = Label), colour = "grey70", linewidth = 0.4) +
    geom_point(aes(fill = Direction, shape = InTrend), size = 3, colour = "black", stroke = 0.4) +
    scale_fill_manual(values = DIR_COL) +
    scale_shape_manual(values = c(Yes = 21, No = 24), name = "In top-50\ntrend list") +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    labs(x = xl, y = NULL, title = title) +
    theme(axis.text.y = element_text(size = 8))
}

# 10a. Classification: Hurley 2 vs Hurley 3
idx <- meta$Stage %in% c("Hurley 2", "Hurley 3")
y23 <- droplevels(meta$Stage[idx]); X23 <- Xrf[idx, ]
set.seed(SEED)
rf_c <- randomForest(X23, y23, ntree = 1000, importance = TRUE)
oob_acc <- 1 - tail(rf_c$err.rate[, "OOB"], 1)
loo_pred <- sapply(seq_len(nrow(X23)), function(i) {
  set.seed(SEED + i)
  m <- randomForest(X23[-i, ], y23[-i], ntree = 1000)
  as.character(predict(m, X23[i, , drop = FALSE])) })
loo_acc <- mean(loo_pred == as.character(y23))
set.seed(SEED)
null_acc <- replicate(N_PERM_RF, {
  m <- randomForest(X23, sample(y23), ntree = 500)
  1 - tail(m$err.rate[, "OOB"], 1) })
p_perm_c <- (sum(null_acc >= oob_acc) + 1) / (N_PERM_RF + 1)
cat(sprintf("RF H2 vs H3: OOB acc = %.2f | LOOCV acc = %.2f | perm p = %s\n",
            oob_acc, loo_acc, fmt_p(p_perm_c)))
p_rfc <- imp_plot(importance(rf_c), "MeanDecreaseAccuracy",
                  "Random Forest: Hurley 2 vs Hurley 3", "Mean decrease in accuracy") +
  labs(subtitle = sprintf("OOB accuracy = %.2f | LOOCV accuracy = %.2f | permutation p = %s",
                          oob_acc, loo_acc, fmt_p(p_perm_c)))

# 10b. Regression on ordinal stage (all samples)
set.seed(SEED)
rf_r <- suppressWarnings(randomForest(Xrf, meta$Stage_num, ntree = 1000, importance = TRUE))
r2_obs <- tail(rf_r$rsq, 1)
set.seed(SEED)
null_r2 <- replicate(N_PERM_RF, tail(suppressWarnings(randomForest(Xrf, sample(meta$Stage_num),
                                                   ntree = 500))$rsq, 1))
p_perm_r <- (sum(null_r2 >= r2_obs) + 1) / (N_PERM_RF + 1)
cat(sprintf("RF regression: OOB pseudo-R2 = %.2f | perm p = %s\n", r2_obs, fmt_p(p_perm_r)))
p_rfr <- imp_plot(importance(rf_r), "%IncMSE",
                  "Random Forest regression on Hurley stage", "% increase in MSE") +
  labs(subtitle = sprintf("OOB pseudo-R\u00b2 = %.2f | permutation p = %s",
                          r2_obs, fmt_p(p_perm_r)))
save_gg(p_rfc + p_rfr + plot_layout(guides = "collect") &
          theme(plot.subtitle = element_text(size = 8)),
        "09a_RandomForest_importance", 14, 7)

perm_df <- rbind(data.frame(Model = "Classification (H2 vs H3): OOB accuracy",
                            v = null_acc, obs = oob_acc),
                 data.frame(Model = "Regression (stage): OOB pseudo-R\u00b2",
                            v = null_r2, obs = r2_obs))
p_perm <- ggplot(perm_df, aes(v)) +
  geom_histogram(bins = 25, fill = "grey75", colour = "white") +
  geom_vline(aes(xintercept = obs), colour = "#B2182B", linewidth = 0.8) +
  facet_wrap(~Model, scales = "free") +
  labs(x = "Statistic under permuted labels", y = "Count",
       title = "Random Forest permutation tests (red = observed)")
save_gg(p_perm, "09b_RandomForest_permutation", 9, 3.8)

rf_imp <- data.frame(annot[colnames(Xrf), ],
                     MDA_H2vsH3 = importance(rf_c)[, "MeanDecreaseAccuracy"],
                     MDG_H2vsH3 = importance(rf_c)[, "MeanDecreaseGini"],
                     IncMSE_stage = importance(rf_r)[, "%IncMSE"],
                     In_top50_trend = colnames(Xrf) %in% top50)
rf_imp <- rf_imp[order(-rf_imp$MDA_H2vsH3), ]

# 10c. Boruta (optional)
if (RUN_BORUTA) {
  set.seed(SEED)
  bor <- suppressWarnings(Boruta::Boruta(X23, y23, maxRuns = 300, doTrace = 0))
  rf_imp$Boruta_H2vsH3 <- as.character(bor$finalDecision[rownames(rf_imp)])
  cat("Boruta decisions:\n"); print(table(bor$finalDecision))
}
write.csv(rf_imp, file.path(TAB, "T07_RandomForest_importance.csv"), row.names = FALSE)

## ======================= 11. sPLS-DA =======================================
set.seed(SEED)
spl <- splsda(t(L), meta$Stage, ncomp = 2, keepX = c(30, 30), scale = TRUE)
perf_spl <- perf(spl, validation = "loo", progressBar = FALSE)
ber <- perf_spl$error.rate$BER[, "max.dist"]
sc  <- data.frame(spl$variates$X[, 1:2], meta)
colnames(sc)[1:2] <- c("C1", "C2")
ev  <- round(100 * spl$prop_expl_var$X[1:2], 1)
p_spl <- ggplot(sc, aes(C1, C2, fill = Stage, colour = Stage)) +
  stat_ellipse(data = subset(sc, Stage %in% ell_ok), geom = "polygon",
               alpha = 0.12, level = 0.95, type = "norm", linewidth = 0.4, show.legend = FALSE) +
  geom_point(size = 3.2, shape = 21, colour = "black", stroke = 0.4) +
  geom_text_repel(aes(label = Sample), size = 3, family = FONT,
                  colour = "black", seed = SEED) +
  scale_fill_manual(values = STAGE_COL, limits = names(STAGE_COL)) + scale_colour_manual(values = STAGE_COL, limits = names(STAGE_COL)) +
  labs(x = paste0("Component 1 (", ev[1], "%)"), y = paste0("Component 2 (", ev[2], "%)"),
       title = "sPLS-DA",
       subtitle = sprintf("LOOCV balanced error rate: comp 1 = %.2f, comp 2 = %.2f",
                          ber[1], ber[2]))
ld <- selectVar(spl, comp = 1)$value
ld <- head(ld[order(-abs(ld$value.var)), , drop = FALSE], 20)
ld_df <- data.frame(ID = rownames(ld), Loading = ld$value.var) |>
  mutate(Label = trunc_lab(annot[ID, "Label"], 40),
         Highest = sapply(ID, function(i)
           levels(meta$Stage)[which.max(tapply(L[i, ], meta$Stage, median))]),
         Label = factor(Label, levels = Label[order(Loading)]))
p_ld <- ggplot(ld_df, aes(Loading, Label, fill = Highest)) +
  geom_col(colour = "black", linewidth = 0.25, width = 0.75) +
  scale_fill_manual(values = STAGE_COL, name = "Highest\nmedian in") +
  labs(x = "Loading (component 1)", y = NULL, title = "Top sPLS-DA loadings") +
  theme(axis.text.y = element_text(size = 8))
save_gg(p_spl + p_ld + plot_layout(widths = c(1, 1.1)), "10_sPLSDA", 13, 5.5)
write.csv(data.frame(annot[ld_df$ID, ], Loading_comp1 = ld_df$Loading),
          file.path(TAB, "T08_sPLSDA_loadings_comp1.csv"), row.names = FALSE)

## ======================= 12. SENSITIVITY ANALYSES ==========================
both <- unique(c(SUSPECT_SAMPLES, LOW_SIGNAL))
sens_sets <- list(setdiff(meta$Sample, SUSPECT_SAMPLES),
                  setdiff(meta$Sample, LOW_SIGNAL),
                  setdiff(meta$Sample, both))
names(sens_sets) <- c(sprintf("Without low-correlation (n = %d)", length(SUSPECT_SAMPLES)),
                      sprintf("Without low-signal (n = %d)", length(LOW_SIGNAL)),
                      sprintf("Without all flagged (n = %d)", length(both)))
sens_sets <- sens_sets[c(length(SUSPECT_SAMPLES), length(LOW_SIGNAL), length(both)) > 0]
sens_sets <- sens_sets[!duplicated(sens_sets)]
writeLines(c(paste("Low-correlation:", paste(SUSPECT_SAMPLES, collapse = ", ")),
             paste("Low-signal:", paste(LOW_SIGNAL, collapse = ", "))),
           file.path(TAB, "T10b_flagged_samples.txt"))
sens_res <- lapply(names(sens_sets), function(nm_s) {
  s  <- sens_sets[[nm_s]]
  md <- meta[s, ]
  # renormalise within the subset so excluded samples do not shape PQN reference
  Ls <- log2(pqn(Xi[, s])$norm)
  tr <- run_trend(Ls, md)
  Ds <- vegdist(t(log2(pqn(Xi[, s])$norm + 1)), "bray")
  pf <- adonis2(Ds ~ Stage_num, data = md, permutations = 9999)
  list(name = nm_s, trend = tr, R2 = pf$R2[1], P = pf$`Pr(>F)`[1])
})
sens_tab <- trend[, c("ID", "Label", "log2FC_per_stage", "P_limma", "FDR_limma")]
for (r in sens_res) {
  m <- r$trend[match(sens_tab$ID, r$trend$ID), ]
  sens_tab[[paste0("log2FC_", r$name)]] <- m$log2FC_per_stage
  sens_tab[[paste0("P_", r$name)]]      <- m$P_limma
}
write.csv(sens_tab, file.path(TAB, "T09_sensitivity_trend.csv"), row.names = FALSE)
sens_sum <- data.frame(
  Analysis = c("All samples", sapply(sens_res, `[[`, "name")),
  n = c(nrow(meta), sapply(sens_sets, length)),
  PERMANOVA_trend_R2 = c(perm_t$R2[1], sapply(sens_res, `[[`, "R2")),
  PERMANOVA_trend_P  = c(perm_t$`Pr(>F)`[1], sapply(sens_res, `[[`, "P")),
  Top50_still_p05_same_direction = c(TOP_N_HEAT, sapply(sens_res, function(r) {
    m <- r$trend[match(top50, r$trend$ID), ]
    o <- trend[match(top50, trend$ID), ]
    sum(m$P_limma < P_CUT & sign(m$log2FC_per_stage) == sign(o$log2FC_per_stage)) })))
write.csv(sens_sum, file.path(TAB, "T10_sensitivity_summary.csv"), row.names = FALSE)
print(sens_sum)

sc_df <- do.call(rbind, lapply(sens_res, function(r) {
  m <- r$trend[match(trend$ID, r$trend$ID), ]
  data.frame(Analysis = r$name, main = trend$log2FC_per_stage,
             sens = m$log2FC_per_stage, Top = trend$ID %in% top50)
}))
sc_df$Analysis <- factor(sc_df$Analysis, levels = names(sens_sets))
p_sens <- ggplot(sc_df, aes(main, sens)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50", linewidth = 0.3) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_point(data = subset(sc_df, !Top), colour = "grey75", size = 0.8, alpha = 0.6) +
  geom_point(data = subset(sc_df, Top), colour = "#B2182B", size = 1.6) +
  facet_wrap(~Analysis, nrow = 1) +
  labs(x = expression(log[2]~"FC per stage (all samples)"),
       y = expression(log[2]~"FC per stage (subset)"),
       title = "Sensitivity of stage-trend effects to sample exclusion",
       subtitle = sprintf("Red = top %d trend features in the full analysis", TOP_N_HEAT))
save_gg(p_sens, "11_sensitivity_trend_effects", 11, 4.2)

## ======================= 13. SESSION ========================================
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
cat("\nDone. Figures in", FIG, "| tables in", TAB, "\n")
