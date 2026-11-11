###############################################################################
#  HS skin metabolomics - FOLLOW-UP ANALYSIS
#  ---------------------------------------------------------------------------
#  Step 1  Remove likely salt-cluster / background ions (negative mass defect)
#  Step 2  Group isotopes and adducts (13C, Na/H, K/H, K/Na) into one feature
#          -> re-run the stage-trend analysis on the cleaned feature set
#  Step 4  Lipid class ratios (e.g. Cer/SM, PC/PE, ether/diacyl) vs stage
#  Step 5  Chain length and saturation (weighted mean C, DB; SFA/MUFA/PUFA)
#
#  Input : Skin_data_HS.csv, Metadata.csv (same as main script)
#  Output: ./HS_results_followup/  figures 12-19 (PDF/PNG/TIFF) + tables T11-T19
###############################################################################

## ======================= 0. SETTINGS =======================================
DATA_FILE <- "data/Skin_data_HS.csv"
META_FILE <- "data/Metadata_skin.csv"
OUT_DIR   <- "results/skin_followup"

DETECT_FRAC <- 0.80
P_CUT       <- 0.05
FC_CUT      <- 1
TOP_N_HEAT  <- 50
TOP_N_BOX   <- 12
LOW_SIGNAL  <- c("9-MT", "10-MT")
SUSPECT     <- c("1-MT", "2-MT")

# Step 1: mass-defect filter
# Fractional mass >= MD_NEG_START is treated as a NEGATIVE mass defect.
# Organic (CHNOPS) metabolites below ~m/z 950 cannot have a negative mass
# defect; such ions are typically salt clusters, matrix/calibrant ions.
MD_NEG_START <- 0.85
MD_MAX_MZ    <- 950     # above this, H-rich lipids (e.g. TG) can wrap past .99

# Step 2: isotope / adduct grouping
GROUP_PPM <- 5          # m/z tolerance
GROUP_R   <- 0.90       # min Pearson r (log2 intensities) between partners
DELTAS <- c("Split peak (same ion)" = 0,          # peaks picked twice, < 5 ppm apart
            "13C isotope (M+1)"   = 1.003355,
            "13C2 isotope (M+2)"  = 2.006710,
            "[M+Na] vs [M+H]"     = 21.981944,
            "[M+K] vs [M+H]"      = 37.955882,
            "[M+K] vs [M+Na]"     = 15.973938)
NULL_DELTAS <- c(4.4, 6.7, 9.3, 12.1, 14.6, 25.3, 31.7, 44.2)  # nonsense shifts
FONT <- "Arial"
SEED <- 2026
CHAIN_MIN_SPECIES <- 8   # min parsed species per class for chain/saturation metrics

## ======================= 0b. PACKAGES & STYLE ==============================
cran <- c("ggplot2", "ggrepel", "patchwork", "dplyr", "tidyr", "clinfun",
          "circlize", "ragg", "systemfonts", "scales")
bioc <- c("limma", "EnhancedVolcano", "ComplexHeatmap")
miss <- cran[!cran %in% rownames(installed.packages())]
if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
miss <- bioc[!bioc %in% rownames(installed.packages())]
if (length(miss)) BiocManager::install(miss, update = FALSE)
suppressPackageStartupMessages({
  library(ggplot2); library(ggrepel); library(patchwork); library(dplyr)
  library(tidyr); library(limma); library(EnhancedVolcano)
  library(ComplexHeatmap); library(circlize); library(clinfun); library(scales)
})
set.seed(SEED)
FIG <- file.path(OUT_DIR, "figures"); TAB <- file.path(OUT_DIR, "tables")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
if (!FONT %in% systemfonts::system_fonts()$family) FONT <- "sans"

STAGE_COL <- c("Hurley 1" = "#FDAE6B", "Hurley 2" = "#E6550D", "Hurley 3" = "#7F2704")
DIR_COL   <- c("Up with stage" = "#B2182B", "Down with stage" = "#2166AC")
theme_pub <- function(base = 11) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(axis.text = element_text(colour = "black"),
          axis.line = element_line(linewidth = 0.5),
          axis.ticks = element_line(linewidth = 0.5, colour = "black"),
          legend.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", size = base + 2),
          plot.subtitle = element_text(size = base - 1, colour = "grey30"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"))
}
theme_set(theme_pub())
save_gg <- function(p, name, w = 7, h = 5) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 600,
         device = ragg::agg_png, bg = "white")
  ggsave(file.path(FIG, paste0(name, ".tiff")), p, width = w, height = h, dpi = 600,
         device = ragg::agg_tiff, compression = "lzw", bg = "white")
}
save_ht <- function(ht, name, w = 8, h = 10) {
  cairo_pdf(file.path(FIG, paste0(name, ".pdf")), width = w, height = h, family = FONT)
  draw(ht, merge_legend = TRUE); dev.off()
  ragg::agg_png(file.path(FIG, paste0(name, ".png")), width = w, height = h,
                units = "in", res = 600, background = "white")
  draw(ht, merge_legend = TRUE); dev.off()
  ragg::agg_tiff(file.path(FIG, paste0(name, ".tiff")), width = w, height = h,
                 units = "in", res = 600, background = "white", compression = "lzw")
  draw(ht, merge_legend = TRUE); dev.off()
}
trunc_lab <- function(x, n = 45) ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "\u2026"), x)
fmt_p <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
                            formatC(p, format = "f", digits = 3))

## ======================= 1. IMPORT (same as main script) ====================
raw  <- read.csv(DATA_FILE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
meta <- read.csv(META_FILE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
colnames(meta) <- c("Sample", "Stage_num")
meta$Sample <- trimws(meta$Sample)
meta$Stage  <- factor(paste("Hurley", meta$Stage_num), levels = names(STAGE_COL))
meta <- meta[order(meta$Stage_num, as.numeric(sub("-MT", "", meta$Sample))), ]
rownames(meta) <- meta$Sample

mz <- as.numeric(raw[[1]])
nm <- trimws(as.character(raw[[2]])); nm[nm == "" | nm == "NA"] <- NA
X  <- as.matrix(raw[, meta$Sample]); storage.mode(X) <- "numeric"
# collapse duplicated m/z rows (as in main script)
key <- sprintf("%.4f", mz); grp <- factor(key, levels = unique(key))
X   <- rowsum(X, grp, reorder = FALSE) / as.vector(table(grp)[levels(grp)])
nm  <- as.character(tapply(nm, grp, function(v) {
  v <- unique(v[!is.na(v)]); if (length(v)) paste(v, collapse = " / ") else NA })[levels(grp)])
mz  <- as.numeric(levels(grp))
ids <- sprintf("F%04d", seq_along(mz))
rownames(X) <- ids
feat <- data.frame(ID = ids, mz = mz, Name = nm, MeanRaw = rowMeans(X),
                   stringsAsFactors = FALSE, row.names = ids)
cat(sprintf("Unique m/z features: %d\n", nrow(X)))

## ======================= helpers: preprocessing & trend =====================
pqn <- function(M) {
  Mt  <- sweep(M, 2, colSums(M), "/") * median(colSums(M))
  ref <- apply(Mt, 1, median)
  sweep(Mt, 2, apply(Mt / ref, 2, median), "/")
}
preprocess <- function(M, md) {
  det  <- sapply(levels(md$Stage), function(s)
    rowMeans(M[, md$Sample[md$Stage == s], drop = FALSE] > 0))
  M    <- M[apply(det, 1, max) >= DETECT_FRAC, , drop = FALSE]
  M    <- t(apply(M, 1, function(v) { if (any(v == 0)) v[v == 0] <- min(v[v > 0]) / 2; v }))
  pqn(M)
}
run_trend <- function(Lm, md) {
  fit <- eBayes(lmFit(Lm, model.matrix(~ md$Stage_num)))
  tt  <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  sp  <- t(apply(Lm, 1, function(v) {
    ct <- suppressWarnings(cor.test(v, md$Stage_num, method = "spearman", exact = FALSE))
    c(ct$estimate, ct$p.value) }))
  jt  <- apply(Lm, 1, function(v)
    suppressWarnings(jonckheere.test(v, md$Stage_num, alternative = "two.sided")$p.value))
  data.frame(ID = rownames(Lm), log2FC_per_stage = tt$logFC, t = tt$t,
             P_limma = tt$P.Value, FDR_limma = p.adjust(tt$P.Value, "BH"),
             Spearman_rho = sp[, 1], P_Spearman = sp[, 2], P_JT = jt,
             Direction = ifelse(tt$logFC > 0, "Up with stage", "Down with stage"),
             stringsAsFactors = FALSE)
}

## ======================= STEP 1. MASS-DEFECT FILTER =========================
frac <- feat$mz - floor(feat$mz)
feat$MassDefect <- ifelse(frac >= MD_NEG_START, frac - 1, frac)
feat$MD_flag <- feat$MassDefect < 0 & feat$mz < MD_MAX_MZ
cat(sprintf("Step 1: %d features with negative mass defect below m/z %d removed\n",
            sum(feat$MD_flag), MD_MAX_MZ))
write.csv(feat[feat$MD_flag, c("ID", "mz", "Name", "MassDefect", "MeanRaw")],
          file.path(TAB, "T11_removed_negative_mass_defect.csv"), row.names = FALSE)

md_df <- feat |>
  mutate(Status = case_when(MD_flag ~ "Removed (negative mass defect)",
                            MassDefect < 0 ~ "Kept (m/z > 950, ambiguous)",
                            TRUE ~ "Kept"))
lab_md <- md_df |> filter(MD_flag, !is.na(Name)) |> arrange(desc(MeanRaw)) |> head(8)
p_md <- ggplot(md_df, aes(mz, MassDefect, colour = Status)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_point(aes(size = log10(MeanRaw)), alpha = 0.6) +
  geom_text_repel(data = lab_md, aes(label = trunc_lab(Name, 30)), size = 2.7,
                  colour = "black", family = FONT, seed = SEED, max.overlaps = 30) +
  scale_colour_manual(values = c("Kept" = "grey60",
                                 "Kept (m/z > 950, ambiguous)" = "#7FA7C9",
                                 "Removed (negative mass defect)" = "#B2182B"),
                      name = NULL) +
  scale_size_continuous(range = c(0.3, 2.5), name = "log10 mean\nintensity") +
  labs(x = "m/z", y = "Mass defect (Da)",
       title = "Step 1: mass-defect filter",
       subtitle = sprintf("%d of %d features removed (negative mass defect, m/z < %d)",
                          sum(feat$MD_flag), nrow(feat), MD_MAX_MZ)) +
  theme(legend.position = "right")
save_gg(p_md, "12_mass_defect_filter", 8, 5)

keep1 <- feat$ID[!feat$MD_flag]

## ======================= STEP 2. ISOTOPE / ADDUCT GROUPING ==================
mz1 <- feat[keep1, "mz"]; o <- order(mz1); ids1 <- keep1[o]; mz1 <- mz1[o]
L1  <- log2(X[ids1, ] + 1)
R1  <- cor(t(L1))
mean1 <- feat[ids1, "MeanRaw"]

find_links <- function(delta, iso = FALSE) {
  res <- list()
  for (i in seq_along(mz1)) {
    tgt <- mz1[i] + delta; tol <- tgt * GROUP_PPM * 1e-6
    lo <- findInterval(tgt - tol, mz1) + 1; hi <- findInterval(tgt + tol, mz1)
    if (hi < lo) next
    for (j in lo:hi) {
      if (delta == 0 && j <= i) next          # split peaks: count each pair once
      if (R1[i, j] >= GROUP_R && (!iso || mean1[j] < mean1[i]))
        res[[length(res) + 1]] <- c(i, j)
    }
  }
  if (length(res)) do.call(rbind, res) else matrix(integer(0), ncol = 2)
}
links <- lapply(names(DELTAS), function(k) find_links(DELTAS[[k]], iso = grepl("isotope", k)))
names(links) <- names(DELTAS)
null_n <- sapply(NULL_DELTAS, function(d) nrow(find_links(d)))
link_tab <- data.frame(Type = names(DELTAS), Links = sapply(links, nrow),
                       Expected_by_chance = mean(null_n))
link_tab$Est_false_rate <- pmin(1, link_tab$Expected_by_chance / pmax(link_tab$Links, 1))
print(link_tab)

# union-find to build groups
parent <- seq_along(mz1)
findp <- function(x) { while (parent[x] != x) x <- parent[x]; x }
for (k in names(links)) for (r in seq_len(nrow(links[[k]]))) {
  a <- findp(links[[k]][r, 1]); b <- findp(links[[k]][r, 2])
  if (a != b) parent[b] <- a
}
comp <- sapply(seq_along(mz1), findp)
grp_df <- data.frame(ID = ids1, mz = mz1, Name = feat[ids1, "Name"],
                     MeanRaw = mean1, Group = comp, stringsAsFactors = FALSE)
grp_df <- grp_df |> group_by(Group) |>
  mutate(GroupSize = n(), IsRep = row_number() == which.max(MeanRaw),
         GroupName = { nn <- Name[order(-MeanRaw)]; nn <- nn[!is.na(nn)]
                       if (length(nn)) nn[1] else NA_character_ },
         Members = paste(sprintf("%.4f", sort(mz)), collapse = "; ")) |>
  ungroup()
reps <- grp_df |> filter(IsRep)
cat(sprintf("Step 2: %d features -> %d compounds (%d groups of 2+)\n",
            nrow(grp_df), nrow(reps), sum(reps$GroupSize > 1)))
write.csv(grp_df |> filter(GroupSize > 1) |> arrange(Group, mz),
          file.path(TAB, "T12_isotope_adduct_groups.csv"), row.names = FALSE)

p_link <- link_tab |>
  pivot_longer(c(Links, Expected_by_chance), names_to = "What", values_to = "n") |>
  mutate(What = recode(What, Links = "Observed links",
                       Expected_by_chance = "Expected by chance"),
         Type = factor(Type, levels = names(DELTAS))) |>
  ggplot(aes(Type, n, fill = What)) +
  geom_col(position = position_dodge(0.75), width = 0.7, colour = "black", linewidth = 0.25) +
  scale_fill_manual(values = c("Observed links" = "#E6550D", "Expected by chance" = "grey75"),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Feature pairs linked",
       title = "Step 2: isotope and adduct grouping",
       subtitle = sprintf("%d ppm, r \u2265 %.2f | %d features merged into %d compounds",
                          GROUP_PPM, GROUP_R, nrow(grp_df), nrow(reps))) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top")
save_gg(p_link, "13_isotope_adduct_grouping", 7, 5)

## ======================= cleaned feature set & annotation ===================
lipid_class <- function(n) {
  if (is.na(n) || grepl(" / ", n)) return(NA_character_)
  if (grepl("phosphocholine", n, ignore.case = TRUE)) return("PC")
  if (grepl("phosphoethanolamine", n, ignore.case = TRUE)) return("PE")
  pats <- c("PE-Cer" = "^PE-Cer", "HexCer" = "^HexCer", "CerP" = "^CerP",
            "Cer" = "^Cer[ (]", "SM" = "^SM[ (]", "LPC" = "^(LPC|LysoPC)",
            "LPE" = "^(LPE|LysoPE)", "PC" = "^PC[ (]", "PE" = "^PE[ (]",
            "PS" = "^PS[ (]", "PI" = "^PI[ (]", "PG" = "^PG[ (]", "PA" = "^PA[ (]",
            "DG" = "^DG[ (]", "TG" = "^TG[ (]", "MG" = "^MG[ (]", "CE" = "^CE[ (]")
  hit <- names(pats)[vapply(pats, grepl, logical(1), x = n)]
  if (!length(hit)) return(NA_character_)
  cl <- hit[1]
  if (cl %in% c("PC", "PE", "PS", "DG") && grepl("^[A-Z]+ ?\\(?[OP]-", n)) cl <- paste0(cl, "-O")
  cl
}
ACYL <- c(lauroyl = "12:0", myristoyl = "14:0", pentadecanoyl = "15:0",
          palmitoyl = "16:0", palmitoleoyl = "16:1", heptadecanoyl = "17:0",
          stearoyl = "18:0", oleoyl = "18:1", linoleoyl = "18:2",
          arachidonoyl = "20:4", docosahexaenoyl = "22:6")
composition <- function(n) {             # returns c(C, DB) summed over chains
  if (is.na(n) || grepl(" / ", n)) return(c(NA, NA))
  low <- tolower(n)
  if (grepl("sn-glycero", low)) {        # systematic names
    m <- regmatches(low, gregexpr(paste(names(ACYL), collapse = "|"), low))[[1]]
    if (!length(m)) return(c(NA, NA))
    if (grepl("^1,2-di", low)) m <- rep(m[1], 2)
    if (length(m) != 2) return(c(NA, NA))
    cd <- do.call(rbind, strsplit(ACYL[m], ":"))
    return(colSums(matrix(as.numeric(cd), ncol = 2)))
  }
  m <- regmatches(n, gregexpr("[0-9]+:[0-9]+", n))[[1]]
  if (!length(m)) return(c(NA, NA))
  cd <- matrix(as.numeric(unlist(strsplit(m, ":"))), ncol = 2, byrow = TRUE)
  colSums(cd)
}

Xc <- X[reps$ID, ]
ann <- data.frame(ID = reps$ID, mz = reps$mz, Name = reps$GroupName,
                  GroupSize = reps$GroupSize, stringsAsFactors = FALSE, row.names = reps$ID)
dupn <- !is.na(ann$Name) & (duplicated(ann$Name) | duplicated(ann$Name, fromLast = TRUE))
ann$Label <- ifelse(is.na(ann$Name), paste0("m/z ", sprintf("%.4f", ann$mz)),
                    ifelse(dupn, paste0(ann$Name, " (", sprintf("%.4f", ann$mz), ")"), ann$Name))
ann$Label <- make.unique(ann$Label, sep = " _")
ann$Class <- vapply(ann$Name, lipid_class, character(1))
comp_m <- t(vapply(ann$Name, composition, numeric(2)))
ann$Carbons <- comp_m[, 1]; ann$DB <- comp_m[, 2]

## ======================= RE-RUN TREND ON CLEANED SET ========================
# reference: original (collapsed, uncleaned) set for comparison
Lo <- log2(preprocess(X, meta))
tr_orig <- run_trend(Lo, meta)

Xn <- preprocess(Xc, meta); L <- log2(Xn)
trend <- run_trend(L, meta)
trend <- cbind(ann[trend$ID, c("mz", "Name", "Label", "Class", "GroupSize")], trend)
# robustness: re-normalise without the low-signal samples
s_rb  <- setdiff(meta$Sample, LOW_SIGNAL)
tr_rb <- run_trend(log2(preprocess(Xc[, s_rb], meta[s_rb, ])), meta[s_rb, ])
m     <- tr_rb[match(trend$ID, tr_rb$ID), ]
trend$P_without_lowsignal      <- m$P_limma
trend$log2FC_without_lowsignal <- m$log2FC_per_stage
trend$Robust <- trend$P_limma < P_CUT & !is.na(m$P_limma) & m$P_limma < P_CUT &
                sign(m$log2FC_per_stage) == sign(trend$log2FC_per_stage)
s_bt  <- setdiff(meta$Sample, SUSPECT)
tr_bt <- run_trend(log2(preprocess(Xc[, s_bt], meta[s_bt, ])), meta[s_bt, ])
m2    <- tr_bt[match(trend$ID, tr_bt$ID), ]
trend$P_without_1MT_2MT <- m2$P_limma
trend <- trend[order(trend$P_limma), ]
write.csv(trend, file.path(TAB, "T13_trend_cleaned_features.csv"), row.names = FALSE)
write.csv(trend[trend$Robust, ], file.path(TAB, "T14_robust_candidates.csv"), row.names = FALSE)

cmp <- data.frame(
  Metric = c("Features tested", "Trend p < 0.05", "  up with stage", "  down with stage",
             "Trend FDR < 0.05", "Robust to removing 9-MT, 10-MT (p < 0.05, same direction)"),
  Before_cleaning = c(nrow(tr_orig), sum(tr_orig$P_limma < P_CUT),
                      sum(tr_orig$P_limma < P_CUT & tr_orig$log2FC_per_stage > 0),
                      sum(tr_orig$P_limma < P_CUT & tr_orig$log2FC_per_stage < 0),
                      sum(tr_orig$FDR_limma < 0.05), NA),
  After_cleaning = c(nrow(trend), sum(trend$P_limma < P_CUT),
                     sum(trend$P_limma < P_CUT & trend$log2FC_per_stage > 0),
                     sum(trend$P_limma < P_CUT & trend$log2FC_per_stage < 0),
                     sum(trend$FDR_limma < 0.05), sum(trend$Robust)))
write.csv(cmp, file.path(TAB, "T15_before_after_cleaning.csv"), row.names = FALSE)
print(cmp)

# Figure 14: trend volcano (cleaned)
trend$Col <- ifelse(trend$Robust, paste(trend$Direction, "(robust)"),
                    ifelse(trend$P_limma < P_CUT, trend$Direction, "NS"))
COL14 <- c("Up with stage (robust)" = "#B2182B", "Down with stage (robust)" = "#2166AC",
           "Up with stage" = "#F4A582", "Down with stage" = "#92C5DE", "NS" = "grey80")
top_lab <- head(trend[trend$Robust, ], 20)
p_tr <- ggplot(trend, aes(log2FC_per_stage, -log10(P_limma), colour = Col)) +
  geom_point(size = 1.3, alpha = 0.85) +
  geom_hline(yintercept = -log10(P_CUT), linetype = 2, linewidth = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_text_repel(data = top_lab, aes(label = trunc_lab(Label, 35)), size = 2.6,
                  colour = "black", family = FONT, max.overlaps = 50,
                  segment.size = 0.25, min.segment.length = 0, seed = SEED) +
  scale_colour_manual(values = COL14, breaks = names(COL14), name = NULL) +
  labs(x = expression(log[2]~"fold change per Hurley stage"),
       y = expression(-log[10]~italic(P)~"(limma trend)"),
       title = "Stage trend after cleaning (Steps 1-2)",
       subtitle = sprintf("%d compounds | %d p < 0.05 (%d up, %d down) | %d robust | %d FDR < 0.05",
                          nrow(trend), sum(trend$P_limma < P_CUT),
                          sum(trend$P_limma < P_CUT & trend$log2FC_per_stage > 0),
                          sum(trend$P_limma < P_CUT & trend$log2FC_per_stage < 0),
                          sum(trend$Robust), sum(trend$FDR_limma < 0.05)))
save_gg(p_tr, "14_trend_volcano_cleaned", 8, 6)

# Figure 15: EnhancedVolcano H3 vs H2 (cleaned)
des <- model.matrix(~ 0 + Stage, data = meta); colnames(des) <- c("S1", "S2", "S3")
fit <- eBayes(contrasts.fit(lmFit(L, des), makeContrasts(S3 - S2, levels = des)))
tt  <- topTable(fit, coef = 1, number = Inf, sort.by = "none")
tt  <- cbind(ann[rownames(tt), c("mz", "Name", "Label", "Class")], tt)
tt$FDR <- p.adjust(tt$P.Value, "BH")
write.csv(tt[order(tt$P.Value), ], file.path(TAB, "T16_limma_H3vsH2_cleaned.csv"),
          row.names = FALSE)
sig <- tt[tt$P.Value < P_CUT & abs(tt$logFC) > FC_CUT, ]
sel <- c(rownames(head(sig[sig$logFC > 0, ][order(sig$P.Value[sig$logFC > 0]), ], 10)),
         rownames(head(sig[sig$logFC < 0, ][order(sig$P.Value[sig$logFC < 0]), ], 10)))
tt$Plab <- ifelse(rownames(tt) %in% sel, trunc_lab(tt$Label, 35), "")
v <- EnhancedVolcano(tt, lab = tt$Plab, x = "logFC", y = "P.Value",
                     selectLab = tt$Plab[rownames(tt) %in% sel],
                     ylim = c(0, max(-log10(tt$P.Value)) * 1.15),
                     pCutoff = P_CUT, FCcutoff = FC_CUT, pointSize = 1.6, labSize = 2.8,
                     col = c("grey75", "#7FA7C9", "#9E9AC8", "#B2182B"), colAlpha = 0.8,
                     drawConnectors = TRUE, widthConnectors = 0.3, colConnectors = "grey40",
                     max.overlaps = 60, legendPosition = "bottom", legendLabSize = 9,
                     legendLabels = c("NS", expression(log[2]~FC), "p-value",
                                      expression(p~and~log[2]~FC)),
                     title = "Hurley 3 vs Hurley 2 (cleaned features)",
                     subtitle = sprintf("Up: %d | Down: %d (p < %.2f, |log2FC| > %g) | FDR < 0.05: %d",
                                        sum(sig$logFC > 0), sum(sig$logFC < 0), P_CUT, FC_CUT,
                                        sum(tt$FDR < 0.05 & abs(tt$logFC) > FC_CUT)),
                     caption = paste0("limma moderated t-test; ", nrow(tt), " compounds"),
                     gridlines.major = FALSE, gridlines.minor = FALSE,
                     axisLabSize = 12, titleLabSize = 13, subtitleLabSize = 9,
                     captionLabSize = 8) +
  theme(text = element_text(family = FONT))
save_gg(v, "15_EnhancedVolcano_H3_vs_H2_cleaned", 7.5, 7.5)

# Figure 16: heatmap of top cleaned trend compounds
top_ids <- head(trend$ID, TOP_N_HEAT)
Z <- t(scale(t(L[top_ids, ]))); rownames(Z) <- trunc_lab(ann[top_ids, "Label"], 50)
ha <- HeatmapAnnotation(Stage = meta$Stage, col = list(Stage = STAGE_COL),
                        annotation_name_gp = gpar(fontsize = 9, fontface = "bold"))
tsub <- trend[match(top_ids, trend$ID), ]
ra <- rowAnnotation(Trend = tsub$Direction,
                    Robust = ifelse(tsub$Robust, "Yes", "No"),
                    col = list(Trend = DIR_COL, Robust = c(Yes = "black", No = "grey90")),
                    annotation_name_gp = gpar(fontsize = 8, fontface = "bold"))
ht <- Heatmap(Z, name = "Z-score", col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
              top_annotation = ha, right_annotation = ra,
              column_split = meta$Stage, cluster_column_slices = FALSE,
              column_title_gp = gpar(fontface = "bold", fontsize = 10),
              clustering_distance_rows = "pearson", row_names_gp = gpar(fontsize = 6.5),
              column_names_gp = gpar(fontsize = 8), row_names_max_width = unit(8, "cm"),
              border = TRUE, rect_gp = gpar(col = "white", lwd = 0.3))
save_ht(ht, "16_heatmap_top_trend_cleaned", 11, 11)

# Figure 17: box plots of top robust candidates
bx_ids <- trend$ID[trend$Robust]          # all robust candidates
if (length(bx_ids)) {
  bd <- data.frame(t(L[bx_ids, , drop = FALSE]), check.names = FALSE) |>
    mutate(Stage = meta$Stage) |>
    pivot_longer(all_of(bx_ids), names_to = "ID", values_to = "v") |>
    left_join(trend[, c("ID", "Label", "P_limma", "P_without_lowsignal")], by = "ID") |>
    mutate(Facet = paste0(trunc_lab(Label, 32), "\np = ", fmt_p(P_limma),
                          " (w/o 9,10-MT: ", fmt_p(P_without_lowsignal), ")"),
           Facet = factor(Facet, levels = unique(Facet[match(bx_ids, ID)])))
  p_bx <- ggplot(bd, aes(Stage, v, fill = Stage)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85, linewidth = 0.3) +
    geom_jitter(width = 0.12, size = 1.4, shape = 21, colour = "black", stroke = 0.3) +
    stat_summary(aes(group = 1), fun = median, geom = "line", linetype = 2, linewidth = 0.4) +
    facet_wrap(~Facet, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = STAGE_COL) + scale_x_discrete(labels = c("1", "2", "3")) +
    labs(x = "Hurley stage", y = expression(log[2]~"PQN intensity"),
         title = sprintf("All %d robust stage-associated compounds", length(bx_ids))) +
    theme(legend.position = "none", strip.text = element_text(size = 7.5))
  save_gg(p_bx, "17_robust_candidates_boxplots", 12, 2.5 * ceiling(length(bx_ids) / 5) + 0.6)
}

# Figure 17b: heatmap of all robust candidates
if (length(bx_ids) > 1) {
  rb <- trend[match(bx_ids, trend$ID), ]
  rb <- rb[order(rb$Direction, -abs(rb$log2FC_per_stage)), ]
  Zr <- t(scale(t(L[rb$ID, ]))); rownames(Zr) <- trunc_lab(rb$Label, 40)
  ht_r <- Heatmap(Zr, name = "Z-score",
                  col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
                  top_annotation = ha, column_split = meta$Stage,
                  cluster_column_slices = FALSE, cluster_rows = FALSE,
                  row_split = factor(rb$Direction, levels = names(DIR_COL)),
                  row_title_gp = gpar(fontsize = 9, fontface = "bold"),
                  left_annotation = rowAnnotation(Trend = rb$Direction,
                                                  col = list(Trend = DIR_COL), show_legend = FALSE,
                                                  show_annotation_name = FALSE),
                  right_annotation = rowAnnotation(
                    `-log10 P` = anno_barplot(-log10(rb$P_limma), width = unit(1.4, "cm"),
                                              gp = gpar(fill = "grey50", col = NA)),
                    annotation_name_gp = gpar(fontsize = 8, fontface = "bold")),
                  column_title_gp = gpar(fontface = "bold", fontsize = 10),
                  row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
                  border = TRUE, rect_gp = gpar(col = "white", lwd = 0.4))
  save_ht(ht_r, "17b_heatmap_robust_candidates", 8.5, 0.22 * length(bx_ids) + 3)
}

## ======================= STEP 4. LIPID CLASS RATIOS =========================
fit_metric <- function(v, md) {
  f  <- summary(lm(v ~ md$Stage_num))$coefficients
  sp <- suppressWarnings(cor.test(v, md$Stage_num, method = "spearman", exact = FALSE))
  c(Slope = f[2, 1], P_lm = f[2, 4], Spearman_rho = unname(sp$estimate), P_Spearman = sp$p.value)
}
cls  <- ann[rownames(Xn), "Class"]
ncls <- table(cls)
csum <- function(k) colSums(Xn[!is.na(cls) & cls == k, , drop = FALSE])
RATIOS <- list(c("Cer", "SM"), c("HexCer", "Cer"), c("PC", "PE"), c("PE-O", "PE"),
               c("PC-O", "PC"), c("DG", "TG"), c("PS", "PE"), c("LPC", "PC"), c("LPE", "PE"))
RATIOS <- Filter(function(r) all(r %in% names(ncls)) && all(ncls[r] >= 2), RATIOS)
rat_mat <- t(sapply(RATIOS, function(r) log2(csum(r[1]) / csum(r[2]))))
rownames(rat_mat) <- sapply(RATIOS, function(r)
  sprintf("%s/%s (n = %d/%d)", r[1], r[2], ncls[r[1]], ncls[r[2]]))
rat_res <- data.frame(Ratio = rownames(rat_mat),
                      t(apply(rat_mat, 1, fit_metric, md = meta)), check.names = FALSE)
rat_res$FDR <- p.adjust(rat_res$P_lm, "BH")
rat_res <- rat_res[order(rat_res$P_lm), ]
write.csv(rat_res, file.path(TAB, "T17_lipid_ratios_trend.csv"), row.names = FALSE)
print(rat_res)

rd <- data.frame(t(rat_mat), check.names = FALSE) |> mutate(Stage = meta$Stage) |>
  pivot_longer(-Stage, names_to = "Ratio", values_to = "v") |>
  left_join(rat_res[, c("Ratio", "P_lm", "FDR")], by = "Ratio") |>
  mutate(Facet = paste0(Ratio, "\np = ", fmt_p(P_lm), ", FDR = ", fmt_p(FDR)),
         Facet = factor(Facet, levels = unique(Facet[order(P_lm)])))
p_rat <- ggplot(rd, aes(Stage, v, fill = Stage)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85, linewidth = 0.3) +
  geom_jitter(width = 0.12, size = 1.4, shape = 21, colour = "black", stroke = 0.3) +
  facet_wrap(~Facet, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = STAGE_COL) + scale_x_discrete(labels = c("1", "2", "3")) +
  labs(x = "Hurley stage", y = expression(log[2]~"ratio of summed class intensity"),
       title = "Step 4: lipid class ratios across Hurley stages") +
  theme(legend.position = "none", strip.text = element_text(size = 8))
save_gg(p_rat, "18_lipid_ratios", 8.5, 2.7 * ceiling(nrow(rat_mat) / 3) + 0.5)

## ======================= STEP 5. CHAIN LENGTH & SATURATION ==================
has_comp <- !is.na(cls) & !is.na(ann[rownames(Xn), "Carbons"])
comp_ann <- ann[rownames(Xn), ][has_comp, ]
comp_ann$Class2 <- sub("-O$", "", comp_ann$Class)           # pool ether with parent class
write.csv(comp_ann[, c("ID", "mz", "Name", "Class", "Carbons", "DB")],
          file.path(TAB, "T18_parsed_lipid_composition.csv"), row.names = FALSE)
use_cl <- names(which(table(comp_ann$Class2) >= CHAIN_MIN_SPECIES))
wmean <- function(ids, w) colSums(Xn[ids, , drop = FALSE] * w) / colSums(Xn[ids, , drop = FALSE])
ch_list <- list()
for (k in use_cl) {
  a <- comp_ann[comp_ann$Class2 == k, ]
  ch_list[[paste0(k, ": mean carbons")]]      <- wmean(a$ID, a$Carbons)
  ch_list[[paste0(k, ": mean double bonds")]] <- wmean(a$ID, a$DB)
}
gpl <- comp_ann[comp_ann$Class2 %in% c("PC", "PE", "PS", "PI", "PG", "PA"), ]
if (nrow(gpl) >= 5) {
  tot <- colSums(Xn[gpl$ID, , drop = FALSE])
  ch_list[["Glycerophospholipids: % saturated (0 DB)"]] <-
    100 * colSums(Xn[gpl$ID[gpl$DB == 0], , drop = FALSE]) / tot
  ch_list[["Glycerophospholipids: % mono-unsaturated (1 DB)"]] <-
    100 * colSums(Xn[gpl$ID[gpl$DB == 1], , drop = FALSE]) / tot
  ch_list[["Glycerophospholipids: % poly-unsaturated (>=2 DB)"]] <-
    100 * colSums(Xn[gpl$ID[gpl$DB >= 2], , drop = FALSE]) / tot
}
ch_mat <- do.call(rbind, ch_list)
n_sp <- sapply(rownames(ch_mat), function(r) {
  k <- sub(":.*", "", r)
  if (k == "Glycerophospholipids") nrow(gpl) else sum(comp_ann$Class2 == k) })
ch_res <- data.frame(Metric = rownames(ch_mat), n_species = n_sp,
                     t(apply(ch_mat, 1, fit_metric, md = meta)), check.names = FALSE)
ch_res$FDR <- p.adjust(ch_res$P_lm, "BH")
ch_res <- ch_res[order(ch_res$P_lm), ]
write.csv(ch_res, file.path(TAB, "T19_chain_saturation_trend.csv"), row.names = FALSE)
print(ch_res)

cd <- data.frame(t(ch_mat), check.names = FALSE) |> mutate(Stage = meta$Stage) |>
  pivot_longer(-Stage, names_to = "Metric", values_to = "v") |>
  left_join(ch_res[, c("Metric", "n_species", "P_lm", "FDR")], by = "Metric") |>
  mutate(Facet = paste0(Metric, " (n = ", n_species, ")\np = ", fmt_p(P_lm),
                        ", FDR = ", fmt_p(FDR)),
         Facet = factor(Facet, levels = unique(Facet[order(P_lm)])))
p_ch <- ggplot(cd, aes(Stage, v, fill = Stage)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85, linewidth = 0.3) +
  geom_jitter(width = 0.12, size = 1.3, shape = 21, colour = "black", stroke = 0.3) +
  facet_wrap(~Facet, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = STAGE_COL) + scale_x_discrete(labels = c("1", "2", "3")) +
  labs(x = "Hurley stage", y = "Intensity-weighted value",
       title = "Step 5: lipid chain length and saturation across Hurley stages") +
  theme(legend.position = "none", strip.text = element_text(size = 7.5))
save_gg(p_ch, "19_chain_length_saturation", 11, 2.6 * ceiling(nrow(ch_mat) / 4) + 0.6)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
cat("\nDone. Figures in", FIG, "| tables in", TAB, "\n")
