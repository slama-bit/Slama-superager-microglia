# =============================================================================
# WP0 — Artifact Firewall & Microglia Object Construction  (v2.1)
# =============================================================================
# RECONSTRUCTED SCRIPT — NOT THE EXECUTED CODE
# Rebuilt from analysis session logs: the WP0_clean_v2 notebook source plus the
# v2.1 identity-panel correction (C3 removed, HEXB added) that was applied
# before the final full re-run. The original executed notebook
# (WP0_clean_v2_executed.ipynb) was lost with the compute instance.
# Re-running this script is NOT guaranteed to reproduce the frozen object
# MD5 (83cefa0c09ece745eb4dc46be79304f7; 10,676 nuclei, 38 donors).
# =============================================================================
#
# 
# **Dataset:** 10x Multiome hippocampal, GSE268609 (Disouky et al. Nature 2026)  
# **Input:** 153,530 nuclei, 39 donors, 5 groups (AD, MCI->PCI, HA, YA, SA)  
# **Output:** Frozen microglia object + QC appendix  
# 
# **Key references:**
# - Disouky et al. Nature 2026 — parent study
# - Sun et al. Cell 2023 — contamination removal (marker dotplot approach)
# - Marsh et al. Nature Neuroscience 2022 — module score contamination + group-specific state protection
# - Lun et al. F1000Research 2016 — MAD-based QC
# - Luecken & Theis MSB 2019 — per-sample adaptive thresholds
# - Wolock et al. Cell Systems 2019 — Scrublet doublet removal

# -----------------------------------------------------------------------------
# ## 0. Configuration & Setup

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(ggrepel)
  library(tidyr)
  library(tibble)
  library(GenomeInfoDb)
  library(AnnotationHub)
  library(readxl)
  library(tools)
})

# ── Paths and parameters ──────────────────────────────────────────────────────
CONFIG <- list(
  input_rds         = "/home/jupyter/seurat_object/GSE268609_seurat_object.rds",
  fragment_file     = "/home/jupyter/fragments/mg_fragments.tsv.gz",
  metadata_xlsx     = "/home/jupyter/WP0_Sal/plots/4_5994754183109025075.xlsx",
  macs2_path        = "/opt/conda/bin/macs2",
  output_dir        = "/home/jupyter/WP0_clean",
  plot_dir          = "/home/jupyter/WP0_clean/plots",
  mad_multiplier    = 3,
  fallback_n        = 50,
  doublet_thresh    = 0.25,
  yield_min         = 20,
  atac_failed_donor = "20",
  seed              = 42
)

# Consistent colour scheme across all plots
GROUP_COLORS <- c(
  "YA"  = "#2166ac",
  "HA"  = "#92c5de",
  "PCI" = "#f4a582",
  "AD"  = "#d6604d",
  "SA"  = "#4dac26"
)

dir.create(CONFIG$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(CONFIG$plot_dir,   recursive = TRUE, showWarnings = FALSE)
set.seed(CONFIG$seed)

MAD_K      <- CONFIG$mad_multiplier
FALLBACK_N <- CONFIG$fallback_n

# ── Tracking infrastructure ───────────────────────────────────────────────────
exclusion_log  <- list()
removal_audits <- list()

removal_audit <- function(before, after, step_name) {
  b <- as.data.frame(table(before$Group)) |> setNames(c("Group","before"))
  a <- as.data.frame(table(after$Group))  |> setNames(c("Group","after"))
  d <- merge(b, a, by = "Group", all.x = TRUE)
  d$after[is.na(d$after)] <- 0
  d$removed     <- d$before - d$after
  d$pct_removed <- round(d$removed / d$before * 100, 2)
  d$step        <- step_name
  mu <- mean(d$pct_removed); sg <- sd(d$pct_removed)
  d$bias_flag   <- d$pct_removed > (mu + 2 * sg)
  d
}

cat("Setup complete. Output:", CONFIG$output_dir, "\n")

# -----------------------------------------------------------------------------
# ## Step 0.0 — Extract Microglia + Rename MCI → PCI

full_obj <- readRDS(CONFIG$input_rds)
DefaultAssay(full_obj) <- "RNA"

# Rename MCI -> PCI for consistency with analysis plan
full_obj@meta.data$Group <- as.character(full_obj@meta.data$Group)
full_obj@meta.data$Group[full_obj@meta.data$Group == "MCI"] <- "PCI"
full_obj@meta.data$Group <- factor(full_obj@meta.data$Group,
                                    levels = c("YA","HA","PCI","AD","SA"))
full_obj@meta.data$donor_id <- as.character(full_obj@meta.data$orig.ident)

cat("Group levels after rename:\n")
print(table(full_obj@meta.data$Group))

mg0 <- subset(full_obj, subset = Cluster == "Microglia")
rm(full_obj); gc()

stopifnot(ncol(mg0) == 12165)
stopifnot(length(unique(mg0$orig.ident)) == 39)

cat("\nMicroglia extracted:", ncol(mg0), "\n")
cat("Donors:", length(unique(mg0$orig.ident)), "\n")
cat("Per group:\n"); print(table(mg0$Group))

write.csv(
  data.frame(cell_id = colnames(mg0),
             donor   = mg0$orig.ident,
             group   = mg0$Group),
  file.path(CONFIG$output_dir, "cell_registry_raw.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# ## Step 0.0b — Low-Yield Donor Removal
# 
# Remove donors with fewer than 20 microglia nuclei before QC filtering.  
# These donors cannot contribute reliable pseudobulk estimates in WP3/WP4.

yield_before <- as.data.frame(mg0@meta.data) |>
  dplyr::mutate(orig.ident = as.character(orig.ident)) |>
  dplyr::group_by(orig.ident, Group) |>
  dplyr::summarise(n_nuclei = dplyr::n(), .groups = "drop") |>
  dplyr::arrange(n_nuclei)

cat("Yield check — bottom 5 donors:\n")
print(head(yield_before, 5))

low_yield_donors <- yield_before$orig.ident[
  yield_before$n_nuclei < CONFIG$yield_min]

cat("\nDonors below yield threshold (<", CONFIG$yield_min, "cells):\n")
print(yield_before[yield_before$orig.ident %in% low_yield_donors, ])

if (length(low_yield_donors) > 0) {
  exclusion_log[["step0.0b_yield"]] <- as.data.frame(mg0@meta.data) |>
    dplyr::filter(as.character(orig.ident) %in% low_yield_donors) |>
    dplyr::mutate(orig.ident = as.character(orig.ident)) |>
    dplyr::group_by(orig.ident, Group) |>
    dplyr::summarise(n_removed = dplyr::n(), step = "step0.0b",
                     reason = paste0("fewer than ", CONFIG$yield_min, " microglia"),
                     .groups = "drop")

  mg0 <- subset(mg0, subset = !(as.character(orig.ident) %in% low_yield_donors))
  cat("Removed:", sum(exclusion_log[["step0.0b_yield"]]$n_removed),
      "cells from", length(low_yield_donors), "donor(s)\n")
}

cat("Nuclei remaining:", ncol(mg0), "\n")
cat("Donors remaining:", length(unique(mg0$orig.ident)), "\n")

# -----------------------------------------------------------------------------
# ## Step 0.1 — RNA QC Visualization (pre-filter)

# Yield per donor barplot
p_yield <- ggplot(yield_before,
  aes(x = reorder(orig.ident, n_nuclei), y = n_nuclei, fill = Group)) +
  geom_col() +
  geom_hline(yintercept = CONFIG$yield_min,
             linetype = "dashed", colour = "red", linewidth = 0.8) +
  scale_fill_manual(values = GROUP_COLORS) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
  labs(title = "Microglia nuclei per donor (pre-QC)",
       subtitle = paste0("Red line = yield threshold (", CONFIG$yield_min, ")"),
       x = "Donor", y = "Nuclei")

# RNA QC violins with boxplot overlay
p_rna_group <- (
  ggplot(as.data.frame(mg0@meta.data),
    aes(x = Group, y = nCount_RNA, fill = Group)) +
  geom_violin(scale="width", alpha=0.8) +
  geom_boxplot(width=0.1, fill="white", outlier.size=0.3) +
  scale_fill_manual(values=GROUP_COLORS) +
  scale_y_log10() +
  theme_classic() + theme(legend.position="none") +
  labs(title="RNA counts", y="nCount_RNA (log10)", x=NULL)
) | (
  ggplot(as.data.frame(mg0@meta.data),
    aes(x = Group, y = nFeature_RNA, fill = Group)) +
  geom_violin(scale="width", alpha=0.8) +
  geom_boxplot(width=0.1, fill="white", outlier.size=0.3) +
  scale_fill_manual(values=GROUP_COLORS) +
  scale_y_log10() +
  theme_classic() + theme(legend.position="none") +
  labs(title="RNA features", y="nFeature_RNA (log10)", x=NULL)
) | (
  ggplot(as.data.frame(mg0@meta.data),
    aes(x = Group, y = percent.MT, fill = Group)) +
  geom_violin(scale="width", alpha=0.8) +
  geom_boxplot(width=0.1, fill="white", outlier.size=0.3) +
  scale_fill_manual(values=GROUP_COLORS) +
  theme_classic() + theme(legend.position="none") +
  labs(title="MT%", y="percent.MT", x=NULL)
)

ggsave(file.path(CONFIG$plot_dir, "fig_QC0_yield_per_donor.pdf"), p_yield, w=14, h=5)
ggsave(file.path(CONFIG$plot_dir, "fig_QC1_RNA_QC_group.pdf"), p_rna_group, w=14, h=5)


# -----------------------------------------------------------------------------
# ## Step 0.2 — ATAC QC + TSS Enrichment + Nucleosome Signal

frag_obj <- CreateFragmentObject(
  path               = CONFIG$fragment_file,
  cells              = colnames(mg0),
  validate.fragments = FALSE
)

ah    <- AnnotationHub()
ensdb <- ah[["AH116291"]]
suppressWarnings(seqlevelsStyle(ensdb) <- "UCSC")
annotations <- suppressWarnings(GetGRangesFromEnsDb(ensdb = ensdb))
genome(annotations) <- "hg38"

mg0[["Peaks"]] <- CreateChromatinAssay(
  counts      = GetAssayData(mg0, assay="Peaks", layer="counts"),
  genome      = "hg38",
  fragments   = list(frag_obj),
  annotation  = annotations
)

DefaultAssay(mg0) <- "Peaks"
mg0 <- TSSEnrichment(mg0)
mg0 <- NucleosomeSignal(mg0)
DefaultAssay(mg0) <- "RNA"

cat("TSS range:",
    round(min(mg0$TSS.enrichment, na.rm=TRUE), 2), "-",
    round(max(mg0$TSS.enrichment, na.rm=TRUE), 2), "\n")

mg0$atac_failed_donor <- as.character(mg0$orig.ident) == CONFIG$atac_failed_donor
cat("ATAC-failed cells (donor", CONFIG$atac_failed_donor, "):",
    sum(mg0$atac_failed_donor), "\n")

writeLines(
  c(paste0("Donor ", CONFIG$atac_failed_donor, " (SA): TSS=0. Flagged as ATAC-failed."),
    "Retained for RNA analyses. Excluded from peak calling."),
  file.path(CONFIG$output_dir, "donor_atac_flags.txt")
)

p_tss <- ggplot(as.data.frame(mg0@meta.data),
  aes(x=Group, y=TSS.enrichment, fill=Group)) +
  geom_violin(scale="width", alpha=0.8) +
  geom_boxplot(width=0.1, fill="white", outlier.size=0.3) +
  geom_hline(yintercept=2, linetype="dashed", colour="red") +
  scale_fill_manual(values=GROUP_COLORS) +
  theme_classic() + theme(legend.position="none") +
  labs(title="TSS enrichment by group", x=NULL)

p_nuc <- ggplot(as.data.frame(mg0@meta.data),
  aes(x=Group, y=nucleosome_signal, fill=Group)) +
  geom_violin(scale="width", alpha=0.8) +
  geom_boxplot(width=0.1, fill="white", outlier.size=0.3) +
  scale_fill_manual(values=GROUP_COLORS) +
  theme_classic() + theme(legend.position="none") +
  labs(title="Nucleosome signal by group", x=NULL)

ggsave(file.path(CONFIG$plot_dir, "fig_QC2_ATAC_QC.pdf"), p_tss|p_nuc, w=10, h=5)

# -----------------------------------------------------------------------------
# ## Step 0.3 — Doublet Removal
# 
# Binary flags non-informative (all False). Score-based at 95th percentile (doublet_score > 0.25).  
# **Reference:** Wolock et al. 2019 Cell Systems

DefaultAssay(mg0) <- "RNA"

cat("doublet_score 95th percentile:",
    round(quantile(mg0$doublet_score, 0.95), 3), "\n")

mg0$doublet_flag <- mg0$doublet_score > CONFIG$doublet_thresh

p_db <- ggplot(as.data.frame(mg0@meta.data),
  aes(x=doublet_score, fill=Group)) +
  geom_histogram(bins=100, alpha=0.7) +
  geom_vline(xintercept=CONFIG$doublet_thresh,
             linetype="dashed", colour="red") +
  scale_fill_manual(values=GROUP_COLORS) +
  facet_wrap(~Group, scales="free_y") +
  theme_classic() +
  labs(title="Doublet score distribution",
       subtitle=paste0("Threshold = ", CONFIG$doublet_thresh))

ggsave(file.path(CONFIG$plot_dir, "fig_QC3_doublet_scores.pdf"), p_db, w=12, h=8)

cat("Doublets flagged:", sum(mg0$doublet_flag), "\n")
cat("Per group:\n"); print(table(mg0$Group[mg0$doublet_flag]))

exclusion_log[["step0.3_doublets"]] <- as.data.frame(mg0@meta.data) |>
  dplyr::filter(doublet_flag) |>
  dplyr::mutate(orig.ident = as.character(orig.ident)) |>
  dplyr::group_by(orig.ident, Group) |>
  dplyr::summarise(n_removed = dplyr::n(), step="step0.3",
                   reason="doublet_score>0.25", .groups="drop")

mg1 <- subset(mg0, subset = doublet_flag == FALSE)
cat("After doublet removal:", ncol(mg1), "\n")

removal_audits[["step0.3"]] <- removal_audit(mg0, mg1, "step0.3")
write.csv(removal_audits[["step0.3"]],
          file.path(CONFIG$output_dir, "removal_audit_step0.3.csv"), row.names=FALSE)
print(removal_audits[["step0.3"]])

# -----------------------------------------------------------------------------
# ## Step 0.4 — RNA MAD Filtering (Donor-Aware)
# 
# 3-MAD on log10(nCount_RNA) and log10(nFeature_RNA) within each donor.  
# Donors < 50 cells use global fallback.  
# **MT filter not applied:** MAD ceiling 2.2% removes 48% PCI — fatal bias.  
# Parent paper used fixed 10% global ceiling applied at CellRanger level.

DefaultAssay(mg1) <- "RNA"
mg1 <- NormalizeData(mg1, verbose=FALSE)

donors_list <- unique(as.character(mg1$orig.ident))

global_lc <- log10(mg1$nCount_RNA   + 1)
global_lf <- log10(mg1$nFeature_RNA + 1)
global_rna <- list(
  count_lo   = 10^(median(global_lc) - MAD_K * mad(global_lc)) - 1,
  count_hi   = 10^(median(global_lc) + MAD_K * mad(global_lc)) - 1,
  feature_lo = 10^(median(global_lf) - MAD_K * mad(global_lf)) - 1
)

donor_rna_thresh <- lapply(setNames(donors_list, donors_list), function(d) {
  idx <- which(as.character(mg1$orig.ident) == d)
  n   <- length(idx)
  if (n < FALLBACK_N) return(c(list(donor=d, n=n, method="global_fallback"), global_rna))
  lc <- log10(mg1$nCount_RNA[idx]   + 1)
  lf <- log10(mg1$nFeature_RNA[idx] + 1)
  list(donor=d, n=n, method="donor_MAD",
       count_lo   = 10^(median(lc) - MAD_K * mad(lc)) - 1,
       count_hi   = 10^(median(lc) + MAD_K * mad(lc)) - 1,
       feature_lo = 10^(median(lf) - MAD_K * mad(lf)) - 1)
})

meta1 <- as.data.frame(mg1@meta.data)
meta1$orig.ident <- as.character(meta1$orig.ident)
meta1$rna_pass <- mapply(function(d, cnt, feat) {
  t <- donor_rna_thresh[[d]]
  cnt >= t$count_lo & cnt <= t$count_hi & feat >= t$feature_lo
}, d=meta1$orig.ident, cnt=meta1$nCount_RNA, feat=meta1$nFeature_RNA)

mg1 <- AddMetaData(mg1, data.frame(rna_pass=meta1$rna_pass, row.names=rownames(meta1)))

exclusion_log[["step0.4_RNA"]] <- meta1 |>
  dplyr::filter(!rna_pass) |>
  dplyr::mutate(orig.ident = as.character(orig.ident)) |>
  dplyr::group_by(orig.ident, Group) |>
  dplyr::summarise(n_removed=dplyr::n(), step="step0.4",
                   reason="RNA_MAD_outlier", .groups="drop")

mg2 <- subset(mg1, subset = rna_pass == TRUE)
cat("After RNA MAD filtering:", ncol(mg2), "\n")

removal_audits[["step0.4"]] <- removal_audit(mg1, mg2, "step0.4")
write.csv(removal_audits[["step0.4"]],
          file.path(CONFIG$output_dir, "removal_audit_step0.4.csv"), row.names=FALSE)
print(removal_audits[["step0.4"]])

writeLines(
  c("MT filter not applied at Step 0.4.",
    "MAD ceiling = 2.2% would remove 48% of PCI nuclei — fatal group bias.",
    "Parent paper (Disouky et al. 2026) used fixed 10% ceiling applied globally."),
  file.path(CONFIG$output_dir, "MT_filter_decision.txt")
)

# -----------------------------------------------------------------------------
# ## Step 0.5 — ATAC MAD Filtering (Donor-Aware)

global_lcp <- log10(mg2$nCount_Peaks  + 1)
global_lfp <- log10(mg2$nFeature_Peaks + 1)
global_atac <- list(
  count_lo = 10^(median(global_lcp) - MAD_K * mad(global_lcp)) - 1,
  count_hi = 10^(median(global_lcp) + MAD_K * mad(global_lcp)) - 1,
  feat_lo  = 10^(median(global_lfp) - MAD_K * mad(global_lfp)) - 1
)

donor_atac_thresh <- lapply(setNames(donors_list, donors_list), function(d) {
  idx <- which(as.character(mg2$orig.ident) == d)
  n   <- length(idx)
  if (n < FALLBACK_N) return(c(list(donor=d, n=n, method="global_fallback"), global_atac))
  lcp <- log10(mg2$nCount_Peaks[idx]  + 1)
  lfp <- log10(mg2$nFeature_Peaks[idx] + 1)
  list(donor=d, n=n, method="donor_MAD",
       count_lo = 10^(median(lcp) - MAD_K * mad(lcp)) - 1,
       count_hi = 10^(median(lcp) + MAD_K * mad(lcp)) - 1,
       feat_lo  = 10^(median(lfp) - MAD_K * mad(lfp)) - 1)
})

meta2 <- as.data.frame(mg2@meta.data)
meta2$orig.ident <- as.character(meta2$orig.ident)
meta2$atac_pass <- mapply(function(d, pcnt, pfeat) {
  t <- donor_atac_thresh[[d]]
  pcnt >= t$count_lo & pcnt <= t$count_hi & pfeat >= t$feat_lo
}, d=meta2$orig.ident, pcnt=meta2$nCount_Peaks, pfeat=meta2$nFeature_Peaks)

mg2 <- AddMetaData(mg2, data.frame(atac_pass=meta2$atac_pass, row.names=rownames(meta2)))

exclusion_log[["step0.5_ATAC"]] <- meta2 |>
  dplyr::filter(!atac_pass) |>
  dplyr::mutate(orig.ident = as.character(orig.ident)) |>
  dplyr::group_by(orig.ident, Group) |>
  dplyr::summarise(n_removed=dplyr::n(), step="step0.5",
                   reason="ATAC_MAD_outlier", .groups="drop")

mg3 <- subset(mg2, subset = atac_pass == TRUE)
cat("After ATAC MAD filtering:", ncol(mg3), "\n")

removal_audits[["step0.5"]] <- removal_audit(mg2, mg3, "step0.5")
write.csv(removal_audits[["step0.5"]],
          file.path(CONFIG$output_dir, "removal_audit_step0.5.csv"), row.names=FALSE)
print(removal_audits[["step0.5"]])

# -----------------------------------------------------------------------------
# ## Step 0.6 — Per-Cell Identity Filter (Sun et al. 2023 approach)
# 
# Compute detection-fraction scores for microglial identity and contamination markers.  
# Conservative rule: remove only nuclei with **zero** microglial markers AND **2+** contamination markers.  
# Store scores as metadata covariates for downstream models.  
# **Reference:** Sun et al. Cell 2023; Marsh et al. Nature Neuroscience 2022

DefaultAssay(mg3) <- "RNA"
counts <- GetAssayData(mg3, layer="counts")

# Microglia identity markers
# Identity panel v2.1. C3 (ENSG00000125730) removed: 0% detection in nuclei and
# expressed by reactive astrocytes. HEXB added (12.9% detection).
mg_markers <- c(
  "ENSG00000183160",  # TMEM119
  "ENSG00000173369",  # C1QB
  "ENSG00000106605",  # P2RY12
  "ENSG00000167397",  # CX3CR1
  "ENSG00000184007",  # CSF1R
  "ENSG00000049860"   # HEXB
)

# Contamination markers
contam_markers <- c(
  "ENSG00000171885",  # AQP4   — astrocyte
  "ENSG00000131095",  # GFAP   — astrocyte
  "ENSG00000132341",  # SNAP25 — neuron
  "ENSG00000108700",  # SYT1   — neuron
  "ENSG00000123560",  # PLP1   — oligodendrocyte
  "ENSG00000168314",  # MOBP   — oligodendrocyte
  "ENSG00000102755",  # FLT1   — endothelial
  "ENSG00000135083",  # CD163  — BAM
  "ENSG00000143878"   # MRC1   — BAM
)

mg_present     <- mg_markers[mg_markers %in% rownames(counts)]
contam_present <- contam_markers[contam_markers %in% rownames(counts)]

cat("Microglia markers found:", length(mg_present), "of", length(mg_markers), "\n")
cat("Contam markers found:   ", length(contam_present), "of", length(contam_markers), "\n")

# Detection fraction scores (0-1)
mg_score     <- colMeans(counts[mg_present,     ] > 0)
contam_score <- colMeans(counts[contam_present, ] > 0)

# Store as metadata
mg3$mg_identity_score <- mg_score
mg3$contam_score      <- contam_score

# Conservative rule: remove ONLY if zero MG markers AND 2+ contam markers
zero_mg      <- mg_score == 0
high_contam  <- contam_score >= (2 / length(contam_present))
clear_contam <- zero_mg & high_contam

cat("\nNuclei with zero MG markers:           ", sum(zero_mg), "\n")
cat("Of those, also have 2+ contam markers: ", sum(clear_contam), "\n")
cat(sprintf("Removal rate: %.2f%%\n", mean(clear_contam)*100))

cat("\nPer group:\n")
print(tapply(clear_contam, mg3$Group, function(x)
  sprintf("%d removed / %d (%.1f%%)", sum(x), length(x), mean(x)*100)))

exclusion_log[["step0.6_identity"]] <- as.data.frame(mg3@meta.data)[clear_contam, ] |>
  dplyr::mutate(orig.ident = as.character(orig.ident)) |>
  dplyr::group_by(orig.ident, Group) |>
  dplyr::summarise(n_removed=dplyr::n(), step="step0.6",
                   reason="zero MG markers AND 2+ contam markers", .groups="drop")

mg4 <- subset(mg3, cells = colnames(mg3)[!clear_contam])
cat("\nAfter identity filter:", ncol(mg4), "\n")

removal_audits[["step0.6"]] <- removal_audit(mg3, mg4, "step0.6")
write.csv(removal_audits[["step0.6"]],
          file.path(CONFIG$output_dir, "removal_audit_step0.6.csv"), row.names=FALSE)
print(removal_audits[["step0.6"]])

# Identity score violin
p_id <- ggplot(as.data.frame(mg4@meta.data),
  aes(x=Group, y=mg_identity_score, fill=Group)) +
  geom_violin(scale="width", alpha=0.8) +
  geom_boxplot(width=0.1, fill="white", outlier.size=0.3) +
  scale_fill_manual(values=GROUP_COLORS) +
  theme_classic() + theme(legend.position="none") +
  labs(title="Microglial identity score by group",
       subtitle="Fraction of MG markers detected per nucleus",
       x=NULL, y="MG identity score")

ggsave(file.path(CONFIG$plot_dir, "fig_QC6_identity_score.pdf"), p_id, w=7, h=5)

# -----------------------------------------------------------------------------
# ## Step 0.7 — Canonical Marker Verification
# 
# Verify microglial identity using canonical markers. **No filtering applied.**  
# Detection rates shown as bar chart — more informative than violin for sparse snRNA-seq markers.

DefaultAssay(mg4) <- "RNA"

canonical_markers <- c(
  "ENSG00000106605",  # P2RY12
  "ENSG00000167397",  # CX3CR1
  "ENSG00000183160",  # TMEM119
  "ENSG00000184007",  # CSF1R
  "ENSG00000049860"   # HEXB
)
markers_present <- canonical_markers[canonical_markers %in% rownames(mg4)]
symbols         <- c("P2RY12","CX3CR1","TMEM119","CSF1R","HEXB")[canonical_markers %in% rownames(mg4)]

counts_mat <- GetAssayData(mg4, layer="counts")[markers_present, , drop=FALSE]
detection_df <- data.frame(
  ensembl_id   = markers_present,
  symbol       = symbols,
  pct_detected = round(rowMeans(counts_mat > 0) * 100, 1)
)
print(detection_df)

write.csv(detection_df,
          file.path(CONFIG$output_dir, "canonical_marker_detection.csv"),
          row.names=FALSE)

detection_df$symbol <- factor(detection_df$symbol,
                               levels=c("P2RY12","CX3CR1","TMEM119","CSF1R","HEXB"))

p_detection <- ggplot(detection_df,
  aes(x=symbol, y=pct_detected, fill=symbol)) +
  geom_col(width=0.6) +
  geom_text(aes(label=paste0(pct_detected,"%")), vjust=-0.5, size=4) +
  scale_fill_brewer(palette="Set2") +
  theme_classic() + theme(legend.position="none") +
  ylim(0,35) +
  labs(title="Canonical microglial marker detection rates",
       subtitle="Post-QC — low detection expected in snRNA-seq",
       x="Marker", y="% nuclei detected")

ggsave(file.path(CONFIG$plot_dir, "fig_QC7_canonical_markers.pdf"), p_detection, w=8, h=5)

# -----------------------------------------------------------------------------
# ## Step 0.8 — Final Bias Audit

mg5 <- mg4

retention_final <- data.frame(
  Group   = c("AD","PCI","HA","YA","SA"),
  n_start = c(3900, 1618, 3606, 1124, 1917),
  n_final = as.integer(table(mg5$Group)[c("AD","PCI","HA","YA","SA")])
) |>
  dplyr::mutate(
    n_removed    = n_start - n_final,
    pct_retained = round(n_final / n_start * 100, 1)
  )

print(retention_final)
write.csv(retention_final,
          file.path(CONFIG$output_dir, "final_group_retention.csv"), row.names=FALSE)

all_steps <- data.frame(
  step = c("Raw","Post-yield","Post-doublet",
           "Post-RNA MAD","Post-ATAC MAD","Post-identity"),
  n    = c(12165,
           ncol(mg0),
           sum(!mg0$doublet_flag),
           sum(meta1$rna_pass),
           sum(meta2$atac_pass),
           ncol(mg5))
) |> dplyr::mutate(removed = dplyr::lag(n, default=12165) - n)

p_waterfall <- ggplot(all_steps,
  aes(x=reorder(step, seq_len(nrow(all_steps))), y=n)) +
  geom_col(fill="steelblue", alpha=0.8) +
  geom_text(aes(label=paste0("-",removed)), vjust=-0.3, size=3.5) +
  theme_classic() +
  theme(axis.text.x=element_text(angle=30, hjust=1)) +
  labs(title="WP0 nuclei waterfall", x="Step", y="Nuclei retained")

ggsave(file.path(CONFIG$plot_dir, "fig_QC8_waterfall.pdf"), p_waterfall, w=10, h=5)

# -----------------------------------------------------------------------------
# ## Step 0.9 — Pooled Peak Re-calling
# 
# Re-call peaks on clean microglia (pooled, no group stratification).  
# Donor 20 excluded (TSS=0).  
# **Reference:** Disouky et al. 2026; Signac CallPeaks vignette

frag5 <- CreateFragmentObject(
  path=CONFIG$fragment_file, cells=colnames(mg5), validate.fragments=FALSE)

mg5[["Peaks"]] <- CreateChromatinAssay(
  counts=GetAssayData(mg5, assay="Peaks", layer="counts"),
  genome="hg38", fragments=list(frag5), annotation=annotations)

mg5_sub <- subset(mg5, subset = orig.ident != CONFIG$atac_failed_donor)
frag_sub <- CreateFragmentObject(
  path=CONFIG$fragment_file, cells=colnames(mg5_sub), validate.fragments=FALSE)
mg5_sub[["Peaks"]] <- CreateChromatinAssay(
  counts=GetAssayData(mg5_sub, assay="Peaks", layer="counts"),
  genome="hg38", fragments=list(frag_sub), annotation=annotations)

cat("Calling peaks on", ncol(mg5_sub), "cells (donor",
    CONFIG$atac_failed_donor, "excluded)...\n")

dir.create(file.path(CONFIG$output_dir, "macs2_peaks"),
           recursive=TRUE, showWarnings=FALSE)

peaks_mg <- CallPeaks(
  object=mg5_sub[["Peaks"]], group.by=NULL,
  macs2.path=CONFIG$macs2_path,
  outdir=file.path(CONFIG$output_dir, "macs2_peaks/"),
  verbose=TRUE
)

peaks_mg <- keepStandardChromosomes(peaks_mg, pruning.mode="coarse")
saveRDS(peaks_mg, file.path(CONFIG$output_dir, "WP0_microglia_peaks.rds"))

cat("Peaks called:", length(peaks_mg), "\n")
cat("Median width:", median(width(peaks_mg)), "bp\n")

write.csv(
  data.frame(n_peaks=length(peaks_mg), median_width_bp=median(width(peaks_mg)),
             n_cells=ncol(mg5_sub),
             donor_excluded=paste0("Donor ", CONFIG$atac_failed_donor),
             method="MACS2 pooled, no group stratification"),
  file.path(CONFIG$output_dir, "peak_recall_summary.csv"), row.names=FALSE
)

# -----------------------------------------------------------------------------
# ## Step 0.10 — Rebuild ATAC Count Matrix with New Peaks
# 
# Quantify fragment counts within the 82,250 microglia-specific peaks.  
# Replaces inherited global peak counts with microglia-relevant counts.

peaks_mg <- readRDS(file.path(CONFIG$output_dir, "WP0_microglia_peaks.rds"))

frag_final <- CreateFragmentObject(
  path=CONFIG$fragment_file, cells=colnames(mg5), validate.fragments=FALSE)

cat("Building count matrix — takes 10-20 min...\n")
macs2_counts <- FeatureMatrix(
  fragments = list(frag_final),
  features  = peaks_mg,
  cells     = colnames(mg5)
)

cat("Count matrix built:", nrow(macs2_counts), "peaks x",
    ncol(macs2_counts), "cells\n")

mg5[["peaks_mg"]] <- CreateChromatinAssay(
  counts    = macs2_counts,
  ranges    = peaks_mg,
  genome    = "hg38",
  fragments = list(frag_final)
)

cat("Assay peaks_mg added\n")
cat("Assays in object:", paste(Assays(mg5), collapse=", "), "\n")

# -----------------------------------------------------------------------------
# ## Step 0.11 — Donor QC Summary + donor_QC_PCs
# 
# Per-donor metrics: RNA counts, features, ATAC fragments, TSS enrichment,  
# nucleosome signal, MG identity score, doublet rate, microglia yield.  
# **FRIP:** Not computable — single merged fragment file.

DefaultAssay(mg5) <- "RNA"

raw_per_donor <- as.data.frame(table(mg0$orig.ident)) |>
  setNames(c("orig.ident","n_raw")) |>
  dplyr::mutate(orig.ident=as.character(orig.ident))

doublet_removed <- dplyr::bind_rows(
  lapply(exclusion_log, function(x) {
    x$orig.ident <- as.character(x$orig.ident); x
  })
) |>
  dplyr::filter(step=="step0.3") |>
  dplyr::group_by(orig.ident) |>
  dplyr::summarise(n_doublets_removed=sum(n_removed), .groups="drop") |>
  dplyr::mutate(orig.ident=as.character(orig.ident))

donor_qc <- as.data.frame(mg5@meta.data) |>
  dplyr::mutate(orig.ident=as.character(orig.ident)) |>
  dplyr::group_by(orig.ident, Group) |>
  dplyr::summarise(
    median_nCount_RNA     = round(median(nCount_RNA), 1),
    median_nFeature_RNA   = round(median(nFeature_RNA), 1),
    median_percent_MT     = round(median(percent.MT), 3),
    median_nCount_Peaks   = round(median(nCount_Peaks), 1),
    median_nFeature_Peaks = round(median(nFeature_Peaks), 1),
    median_TSS            = round(median(TSS.enrichment, na.rm=TRUE), 2),
    pct_TSS_below_2       = round(mean(TSS.enrichment < 2, na.rm=TRUE)*100, 1),
    median_nucleosome     = round(median(nucleosome_signal, na.rm=TRUE), 2),
    median_mg_score       = round(median(mg_identity_score, na.rm=TRUE), 3),
    n_final               = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::left_join(raw_per_donor, by="orig.ident") |>
  dplyr::mutate(pct_retained=round(n_final/n_raw*100,1)) |>
  dplyr::left_join(doublet_removed, by="orig.ident") |>
  dplyr::mutate(
    n_doublets_removed = tidyr::replace_na(n_doublets_removed, 0),
    doublet_rate_pct   = round(n_doublets_removed/n_raw*100, 2),
    FRIP               = NA_real_,
    FRIP_note          = "Not computable - single merged fragment file"
  )

cat("Donor QC summary:", nrow(donor_qc), "donors\n")
print(head(donor_qc, 5))

# PCA on QC metrics
qc_numeric <- donor_qc |>
  dplyr::select(orig.ident, median_nCount_RNA, median_nFeature_RNA,
                median_percent_MT, median_nCount_Peaks, median_nFeature_Peaks,
                median_TSS, median_nucleosome, median_mg_score,
                n_final, pct_retained, doublet_rate_pct) |>
  tidyr::drop_na() |>
  tibble::column_to_rownames("orig.ident")

pca_res <- prcomp(scale(as.matrix(qc_numeric)), center=FALSE, scale.=FALSE)
var_exp <- round(pca_res$sdev^2 / sum(pca_res$sdev^2) * 100, 1)
cat("Variance explained:", paste0("PC",1:5,": ",var_exp[1:5],"%", collapse=" | "), "\n")

donor_QC_PCs <- as.data.frame(pca_res$x[, 1:5])
colnames(donor_QC_PCs) <- paste0("donor_QC_PC", 1:5)
donor_QC_PCs$orig.ident <- rownames(donor_QC_PCs)

write.csv(donor_QC_PCs,
          file.path(CONFIG$output_dir, "donor_QC_PCs.csv"), row.names=FALSE)

cell_pc <- data.frame(orig.ident=as.character(mg5$orig.ident), row.names=colnames(mg5)) |>
  dplyr::left_join(donor_QC_PCs, by="orig.ident") |>
  dplyr::select(-orig.ident)
mg5 <- AddMetaData(mg5, metadata=cell_pc)

pca_plot_df <- dplyr::left_join(donor_QC_PCs,
  donor_qc[,c("orig.ident","Group")], by="orig.ident")

p_pca <- ggplot(pca_plot_df,
  aes(donor_QC_PC1, donor_QC_PC2, colour=Group, label=orig.ident)) +
  geom_point(size=3) + geom_text_repel(size=2.5, max.overlaps=20) +
  scale_colour_manual(values=GROUP_COLORS) + theme_classic() +
  labs(title="Donor QC PCA",
       subtitle=paste0("PC1: ",var_exp[1],"% | PC2: ",var_exp[2],"%"))

ggsave(file.path(CONFIG$plot_dir, "fig_QC11_donor_QC_PCA.pdf"), p_pca, w=8, h=7)

# -----------------------------------------------------------------------------
# ## Step 0.12 — Clinical Metadata Join

meta_clean <- read_excel(CONFIG$metadata_xlsx, sheet=1) |>
  dplyr::select(
    orig.ident       = `New ID`,
    age_at_death     = `Age at Death`,
    sex              = `Sex`,
    PMI              = `PMI`,
    cognitive_status = `Cognitive Status`,
    AD_neuropath     = `Overall AD neuropathological Change  0=Not  1=Low  2=Intermediate  3=High`,
    BRAAK            = `BRAAK`,
    c_score          = `C score`,
    thal             = `Thal  `,
    CAA_score        = `Overall CAA Score`,
    APOE             = `APOE Genotype`
  ) |>
  dplyr::mutate(
    orig.ident = as.character(orig.ident),
    sex        = tolower(trimws(sex)),
    APOE       = dplyr::na_if(trimws(APOE), "Missing/ unknown"),
    APOE       = dplyr::na_if(APOE, "Missing/Unk")
  )

cat("Clinical metadata loaded:", nrow(meta_clean), "donors\n")

donor_qc <- dplyr::left_join(donor_qc, meta_clean, by="orig.ident")

write.csv(donor_qc,
          file.path(CONFIG$output_dir, "donor_QC_summary.csv"), row.names=FALSE)

cat("\nMissing values per clinical column:\n")
for (col in c("age_at_death","sex","PMI","BRAAK","thal",
              "c_score","AD_neuropath","CAA_score","APOE")) {
  if (col %in% colnames(donor_qc))
    cat(sprintf("  %-15s : %d missing\n", col, sum(is.na(donor_qc[[col]]))))
}

cat("\nDonor QC summary saved:", nrow(donor_qc), "donors,",
    ncol(donor_qc), "columns\n")

# -----------------------------------------------------------------------------
# ## Step 0.13 — Freeze Object

mg5$WP0_passed  <- TRUE
mg5$WP0_version <- "WP0_v2.1"
mg5$WP0_date    <- as.character(Sys.Date())
DefaultAssay(mg5) <- "RNA"

cat("=== WP0 FROZEN OBJECT ===\n")
cat("Nuclei retained:", ncol(mg5), "\n")
cat("Nuclei removed:",  12165 - ncol(mg5), "\n")
cat("Pct retained:",    round(ncol(mg5)/12165*100,1), "%\n")
cat("Donors:",          length(unique(mg5$orig.ident)), "\n")
cat("Groups:\n");       print(table(mg5$Group))
cat("Assays:",          paste(Assays(mg5), collapse=", "), "\n")
cat("donor_QC_PCs:",    all(paste0("donor_QC_PC",1:5) %in% colnames(mg5@meta.data)), "\n")
cat("MG identity score:","mg_identity_score" %in% colnames(mg5@meta.data), "\n")

frozen_path <- file.path(CONFIG$output_dir, "microglia_WP0_frozen.rds")
saveRDS(mg5, frozen_path)

checksum <- tools::md5sum(frozen_path)
write(paste0(checksum,"  ",frozen_path), paste0(frozen_path, ".md5"))
Sys.chmod(frozen_path, mode="0444")

cat("\nFrozen:", frozen_path, "\n")
cat("MD5:",     checksum, "\n")

# -----------------------------------------------------------------------------
# ## Step 0.14 — QC Appendix Assembly

all_exclusions <- dplyr::bind_rows(
  lapply(exclusion_log, function(x) { x$orig.ident <- as.character(x$orig.ident); x })
)
write.csv(all_exclusions,
          file.path(CONFIG$output_dir, "FINAL_exclusion_log.csv"), row.names=FALSE)

writeLines(c(
  "# WP0 QC APPENDIX v2.1",
  paste0("Generated  : ", Sys.time()),
  paste0("Frozen RDS : ", frozen_path),
  paste0("MD5        : ", checksum),
  "",
  "## KEY DECISIONS",
  "Yield filter: donors < 20 microglia removed (donor 32, 4 cells)",
  "MT filter not applied: MAD ceiling 2.2% would remove 48% PCI",
  "Doublet removal: score-based (binary flags non-informative)",
  "Identity filter: zero MG markers AND 2+ contam markers (Sun et al. 2023)",
  "Identity panel v2.1: TMEM119, C1QB, P2RY12, CX3CR1, CSF1R, HEXB (C3 removed)",
  "Peak calling pooled: prevents group-biased peak set",
  "Donor 20 excluded from peak calling: TSS=0",
  "FeatureMatrix: ATAC counts rebuilt against microglia-specific peaks",
  "",
  "## LIMITATIONS",
  "FRIP not computable: single merged fragment file",
  "Ambient RNA correction not applied: parent paper did not apply it",
  paste0("\nTotal retained: ", ncol(mg5), "/12165 (",
         round(ncol(mg5)/12165*100,1), "%)")  
), file.path(CONFIG$output_dir, "WP0_appendix_manifest.txt"))

cat("=== WP0 v2.1 COMPLETE ===\n")
cat("Output:", CONFIG$output_dir, "\n")
cat("Total nuclei retained:", ncol(mg5), "/ 12165 (",
    round(ncol(mg5)/12165*100,1), "%)\n")
