# =============================================================================
# WP1 — RNA-first state discovery, contamination exclusion, annotation,
#       multimodal (WNN) validation and Milo graph construction
# =============================================================================
# RECONSTRUCTED SCRIPT — NOT THE EXECUTED CODE
# Rebuilt from analysis session logs (the Sal_Wp1 notebook pipeline, June
# 2026), including the "Cell 16" fix: Harmony on donor_id rather than
# orig.ident. The executed notebook was lost with the compute instance.
# Re-running is NOT guaranteed to reproduce the frozen WP1 object.
#
# Values recorded from the executed run are given as `# executed run:` comments.
# Lines marked `# TODO` were not recoverable from the logs. Verify them before
# trusting any output.
#
# Not reconstructed here:
#   - scVI sensitivity integration: the Python training script and latent
#     matrix were lost.
#   - The interim edgeR cluster-abundance screen: it is superseded by WP2.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)        # executed run: v5.5.0
  library(Signac)        # executed run: v1.17.1
  library(harmony)
  library(UCell)         # executed run: v2.14.0 (no `layer` argument)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(cluster)       # silhouette
  library(clustree)
  library(randomForest)
  library(mclust)        # adjustedRandIndex
  library(miloR)
  library(SingleCellExperiment)
})

# -----------------------------------------------------------------------------
# 0. Configuration
# -----------------------------------------------------------------------------
CONFIG <- list(
  wp0_rds        = "/home/jupyter/WP0_clean/microglia_WP0_frozen.rds",
  wp0_md5        = "83cefa0c09ece745eb4dc46be79304f7",
  gene_map_csv   = "/home/jupyter/WP1_Sall/wp_recluster/gene_map.csv",  # ensembl_gene_id, hgnc_symbol
  sun2023_csv    = "/home/jupyter/WP1_Sall/sun2023_state_markers.csv",  # columns: microgliaState, gene
  out_dir        = "/home/jupyter/WP1_Sall",
  contam_max     = 0.4,          # WP0 contam_score ceiling
  prelim_dims    = 1:20,         # screening clustering only
  prelim_res     = 0.3,
  n_hvg          = 3000,
  n_pcs_compute  = 50,
  N_PCS          = 25,           # executed run: PC25 = 78.5% cumulative variance
  k_param        = 20,
  res_grid       = c(0.2, 0.4, 0.5, 0.6, 0.8),
  final_res      = 0.4,          # executed run: 14 clusters
  rf_iter        = 20,
  lsi_dims       = 2:30,         # LSI1 excluded (depth-correlated)
  milo_prop      = 0.1,
  seed           = 42
)
set.seed(CONFIG$seed)
for (d in c("objects", "plots", "tables")) {
  dir.create(file.path(CONFIG$out_dir, d), recursive = TRUE, showWarnings = FALSE)
}
p_out <- function(...) file.path(CONFIG$out_dir, ...)

gene_map <- read.csv(CONFIG$gene_map_csv, stringsAsFactors = FALSE)
to_ensembl <- function(symbols, obj_features) {
  ids <- gene_map$ensembl_gene_id[gene_map$hgnc_symbol %in% symbols]
  intersect(ids, obj_features)
}

# -----------------------------------------------------------------------------
# 1. Load WP0 frozen object, add biobank, apply contamination-score ceiling
# -----------------------------------------------------------------------------
stopifnot(unname(tools::md5sum(CONFIG$wp0_rds)) == CONFIG$wp0_md5)
mg <- readRDS(CONFIG$wp0_rds)
mg$donor_id <- as.character(mg$orig.ident)

# Biobank is perfectly collinear with SA. It is stored for WP2-WP8 models ONLY.
# It is never regressed (SCTransform) or integrated (Harmony).
mg$biobank <- ifelse(mg$Group == "SA", "NU_ADRC", "UW_BRaIN")

mg <- subset(mg, subset = contam_score <= CONFIG$contam_max)
cat("After contam_score <=", CONFIG$contam_max, ":", ncol(mg), "cells\n")
# executed run: 10,488 cells

# -----------------------------------------------------------------------------
# 2. Screening normalization + preliminary clusters (contamination audit only)
# -----------------------------------------------------------------------------
mg <- SCTransform(mg, vst.flavor = "v2", variable.features.n = CONFIG$n_hvg,
                  vars.to.regress = "percent.MT", assay = "RNA",
                  new.assay.name = "SCT", conserve.memory = TRUE, verbose = FALSE)
mg <- RunPCA(mg, assay = "SCT", npcs = CONFIG$n_pcs_compute,
             reduction.name = "pca_sct", verbose = FALSE)
mg <- FindNeighbors(mg, dims = CONFIG$prelim_dims, reduction = "pca_sct", verbose = FALSE)
mg <- FindClusters(mg, resolution = CONFIG$prelim_res, random.seed = CONFIG$seed,
                   verbose = FALSE)
mg$prelim_cluster <- Idents(mg)
cat("Preliminary clusters:", nlevels(mg$prelim_cluster), "\n")
# executed run: 12 preliminary clusters

# -----------------------------------------------------------------------------
# 3. Contamination module scoring (UCell)
# Module audit decisions:
#   - FLT1 dropped (monocytes/macrophages express it).
#   - GJA1 and GFAP dropped (broad / reactive).
#   - CSPG4 dropped (pericytes). OCLN dropped (epithelial).
#   - MRC1 dropped (activated microglia in AD).
#   - COLEC12 dropped (no human snRNA-seq validation).
# -----------------------------------------------------------------------------
contamination_modules <- list(
  Neuronal    = c("SNAP25", "SYT1", "RBFOX3", "STMN2"),
  Astrocyte   = c("AQP4", "ALDH1L1", "SLC1A2", "SLC1A3"),
  Oligodendro = c("MOBP", "MOG", "MAG"),
  OPC         = c("PDGFRA", "OLIG1", "OLIG2"),
  Endothelial = c("CLDN5", "VWF", "SLC2A1", "PECAM1"),
  BAM         = c("FOLR2", "F13A1", "LYVE1"),
  Tcell       = c("CD3E", "CD3D", "TRAC"),
  Microglial  = c("P2RY12", "CX3CR1", "TMEM119", "CSF1R", "HEXB")
)
contam_ens <- lapply(contamination_modules, to_ensembl, obj_features = rownames(mg[["RNA"]]))
print(sapply(contam_ens, length))

mg <- AddModuleScore_UCell(mg, features = contam_ens, assay = "RNA", name = "")

# -----------------------------------------------------------------------------
# 4. HuMicA three-criterion exclusion (applied per preliminary cluster)
# A cluster is excluded only if ALL of these hold:
#   (1) its microglial score is below the dataset median;
#   (2) any contamination/microglial ratio is >= 2;
#   (3) FindMarkers top genes confirm it (manual review below).
# -----------------------------------------------------------------------------
humica_audit <- function(obj, cluster_col) {
  obj@meta.data %>%
    group_by(.data[[cluster_col]]) %>%
    summarise(n_cells = n(),
              across(c(Neuronal, Astrocyte, Oligodendro, OPC, Endothelial, BAM, Microglial),
                     mean), .groups = "drop") %>%
    mutate(mg_below_median = Microglial < median(Microglial),
           neuronal_ratio  = Neuronal    / (Microglial + 1e-6),
           astrocyte_ratio = Astrocyte   / (Microglial + 1e-6),
           oligo_ratio     = Oligodendro / (Microglial + 1e-6),
           BAM_ratio       = BAM         / (Microglial + 1e-6),
           ratio_flag      = neuronal_ratio >= 2 | astrocyte_ratio >= 2 |
                             oligo_ratio >= 2 | BAM_ratio >= 2,
           exclude_flag    = mg_below_median & ratio_flag)
}

audit_prelim <- humica_audit(mg, "prelim_cluster")
write.csv(audit_prelim, p_out("tables", "humica_audit_prelim.csv"), row.names = FALSE)
flagged <- as.character(audit_prelim$prelim_cluster[audit_prelim$exclude_flag])
cat("Flagged by criteria 1+2:", flagged, "\n")
# executed run: clusters 5 (oligodendrocyte, ratio 2.79) and 11 (BAM, ratio 5.50)

# Criterion 3: confirm with marker genes before excluding
Idents(mg) <- "prelim_cluster"
for (cl in flagged) {
  mk <- FindMarkers(mg, ident.1 = cl, assay = "RNA", only.pos = TRUE,
                    logfc.threshold = 0.25, max.cells.per.ident = 500, verbose = FALSE)
  top <- head(rownames(mk), 20)
  cat("Cluster", cl, "top markers:",
      gene_map$hgnc_symbol[match(top, gene_map$ensembl_gene_id)], "\n")
}
# TODO: set after manual review of the markers above.
exclude_clusters <- flagged

mg <- subset(mg, subset = !(prelim_cluster %in% exclude_clusters))
cat("After contamination exclusion:", ncol(mg), "cells\n")
print(table(mg$Group))
# executed run: 9,311 cells (YA 629, HA 2,822, PCI 1,107, AD 3,218, SA 1,535)
saveRDS(mg, p_out("objects", "WP1_post_contamination_exclusion.rds"))

# -----------------------------------------------------------------------------
# 5. Final normalization + PCA on the clean object
# -----------------------------------------------------------------------------
mg <- SCTransform(mg, vst.flavor = "v2", variable.features.n = CONFIG$n_hvg,
                  vars.to.regress = "percent.MT", assay = "RNA",
                  new.assay.name = "SCT", conserve.memory = TRUE, verbose = FALSE)
mg <- RunPCA(mg, assay = "SCT", npcs = CONFIG$n_pcs_compute,
             reduction.name = "pca_sct", verbose = FALSE)
var_exp <- mg[["pca_sct"]]@stdev^2
cat("Cumulative variance at PC", CONFIG$N_PCS, ":",
    round(100 * sum(var_exp[1:CONFIG$N_PCS]) / sum(var_exp), 1), "% (of 50 computed PCs)\n")

# -----------------------------------------------------------------------------
# 6. Harmony on donor only (Cell 16 fix)
# -----------------------------------------------------------------------------
mg <- RunHarmony(mg, group.by.vars = "donor_id", reduction = "pca_sct",
                 reduction.save = "harmony", dims.use = 1:CONFIG$N_PCS, verbose = FALSE)

# -----------------------------------------------------------------------------
# 7. UMAPs (uncorrected vs Harmony)
# -----------------------------------------------------------------------------
mg <- RunUMAP(mg, reduction = "pca_sct", dims = 1:CONFIG$N_PCS,
              reduction.name = "umap_uncorrected", seed.use = CONFIG$seed, verbose = FALSE)
mg <- RunUMAP(mg, reduction = "harmony", dims = 1:CONFIG$N_PCS,
              reduction.name = "umap_harmony", seed.use = CONFIG$seed, verbose = FALSE)
p_int <- DimPlot(mg, reduction = "umap_uncorrected", group.by = "donor_id") + NoLegend() +
         ggtitle("Uncorrected — donor") |
         DimPlot(mg, reduction = "umap_harmony", group.by = "donor_id") + NoLegend() +
         ggtitle("Harmony — donor")
ggsave(p_out("plots", "umap_uncorrected_vs_harmony.png"), p_int, width = 12, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# 8. SNN graph + Leiden resolution grid
# -----------------------------------------------------------------------------
mg <- FindNeighbors(mg, reduction = "harmony", dims = 1:CONFIG$N_PCS,
                    k.param = CONFIG$k_param, verbose = FALSE)
harmony_dist <- dist(Embeddings(mg, "harmony")[, 1:CONFIG$N_PCS])
res_summary <- lapply(CONFIG$res_grid, function(r) {
  mg <<- FindClusters(mg, resolution = r, algorithm = 4,
                      random.seed = CONFIG$seed, verbose = FALSE)
  col <- paste0("SCT_snn_res.", r)
  sil <- silhouette(as.integer(mg@meta.data[[col]]), harmony_dist)
  data.frame(resolution = r, n_clusters = nlevels(mg@meta.data[[col]]),
             silhouette = round(mean(sil[, 3]), 4))
})
res_summary <- bind_rows(res_summary)
print(res_summary)
write.csv(res_summary, p_out("tables", "resolution_grid.csv"), row.names = FALSE)
# executed run: 0.2=11, 0.4=14, 0.5=15, 0.6=16, 0.8=19 clusters

ggsave(p_out("plots", "clustree.png"), clustree(mg, prefix = "SCT_snn_res."),
       width = 9, height = 10, dpi = 300)

# -----------------------------------------------------------------------------
# 9. Random Forest stability (Olah et al. 2020 criterion: > 75%)
# -----------------------------------------------------------------------------
harmony_df <- as.data.frame(Embeddings(mg, "harmony")[, 1:CONFIG$N_PCS])
rf_stability <- function(cluster_col, n_iter = CONFIG$rf_iter) {
  acc <- sapply(seq_len(n_iter), function(i) {
    set.seed(i)
    tr <- sample(nrow(harmony_df), floor(0.8 * nrow(harmony_df)))
    te <- setdiff(seq_len(nrow(harmony_df)), tr)
    train <- harmony_df[tr, ]
    train$cluster <- factor(mg@meta.data[[cluster_col]][tr])
    fit <- randomForest(cluster ~ ., data = train, ntree = 100)
    mean(predict(fit, harmony_df[te, ]) == factor(mg@meta.data[[cluster_col]][te]))
  })
  data.frame(resolution = cluster_col, mean_acc = round(100 * mean(acc), 1),
             min_acc = round(100 * min(acc), 1), pass = mean(acc) > 0.75)
}
rf_res <- bind_rows(rf_stability("SCT_snn_res.0.4"), rf_stability("SCT_snn_res.0.5"))
print(rf_res)
write.csv(rf_res, p_out("tables", "rf_stability.csv"), row.names = FALSE)
# executed run: res 0.4 = 90.2% (min 88.4%), res 0.5 = 89.9% (min 88.3%), both PASS

# -----------------------------------------------------------------------------
# 10. Select final resolution and re-audit contamination at that resolution
# -----------------------------------------------------------------------------
final_col <- paste0("SCT_snn_res.", CONFIG$final_res)
mg$leiden_final <- mg@meta.data[[final_col]]
Idents(mg) <- "leiden_final"

audit_final <- humica_audit(mg, "leiden_final")
write.csv(audit_final, p_out("tables", "humica_audit_final_res.csv"), row.names = FALSE)
cat("Clusters flagged at final resolution:",
    as.character(audit_final$leiden_final[audit_final$exclude_flag]), "\n")
# executed run: no cluster met the full rule at res 0.4. The suspicious res-0.6
# clusters 14/15 had normal microglial marker expression and were retained.

# -----------------------------------------------------------------------------
# 11. Sun et al. 2023 state annotation (z-scored relative enrichment)
# Twelve states: MG0-MG8 and MG10-MG12. MG9 does not exist in the paper.
# Raw UCell maxima are dominated by MG0 in every cluster, so annotation
# uses the per-module z-score across clusters.
# -----------------------------------------------------------------------------
sun_tbl <- read.csv(CONFIG$sun2023_csv, stringsAsFactors = FALSE)
sun_names <- c(MG0 = "MG0_Homeostatic", MG1 = "MG1_Neuronal_Surveillance",
               MG2 = "MG2_Inflammatory_I", MG3 = "MG3_Ribosome_Biogenesis",
               MG4 = "MG4_Lipid_Processing", MG5 = "MG5_Phagocytic",
               MG6 = "MG6_Stress_Signature", MG7 = "MG7_Glycolytic",
               MG8 = "MG8_Inflammatory_II", MG10 = "MG10_Inflammatory_III",
               MG11 = "MG11_Antiviral", MG12 = "MG12_Cycling")
sun_modules <- split(sun_tbl$gene, sun_tbl$microgliaState)
sun_modules <- sun_modules[names(sun_names)]
names(sun_modules) <- sun_names
sun_ens <- lapply(sun_modules, to_ensembl, obj_features = rownames(mg[["RNA"]]))
print(sapply(sun_ens, length))   # MT- genes in MG1 may be absent in nuclei

mg <- AddModuleScore_UCell(mg, features = sun_ens, assay = "RNA", name = "")
sun_cols <- names(sun_ens)

sun_cluster <- mg@meta.data %>%
  group_by(leiden_final) %>%
  summarise(n_cells = n(), across(all_of(sun_cols), mean), .groups = "drop")
score_z <- scale(as.matrix(sun_cluster[, sun_cols]))
rownames(score_z) <- paste0("C", sun_cluster$leiden_final)

label_df <- data.frame(
  cluster       = sun_cluster$leiden_final,
  n_cells       = sun_cluster$n_cells,
  top_state     = apply(score_z, 1, function(x) names(which.max(x))),
  top_z         = round(apply(score_z, 1, max), 2),
  second_state  = apply(score_z, 1, function(x) names(sort(x, decreasing = TRUE))[2]),
  second_z      = round(apply(score_z, 1, function(x) sort(x, decreasing = TRUE)[2]), 2)
)
print(label_df)
write.csv(label_df, p_out("tables", "sun2023_zscore_annotation.csv"), row.names = FALSE)

png(p_out("plots", "sun2023_zscore_heatmap.png"), width = 2400, height = 1800, res = 300)
pheatmap(t(score_z), color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         display_numbers = TRUE, number_format = "%.1f",
         main = "Sun et al. 2023 — z-scored enrichment per cluster")
dev.off()

# Interim labels only. Biological names are assigned after WP2/WP3 support.
mg$interim_label <- setNames(label_df$top_state, as.character(label_df$cluster))[
  as.character(mg$leiden_final)]

# -----------------------------------------------------------------------------
# 12. WNN multimodal validation (secondary substrate, not the source of truth)
# Donor 20 (failed ATAC) carries no usable ATAC signal; check its weights below.
# -----------------------------------------------------------------------------
DefaultAssay(mg) <- "peaks_mg"
mg <- RunTFIDF(mg, verbose = FALSE)
mg <- FindTopFeatures(mg, min.cutoff = "q75")
mg <- RunSVD(mg, verbose = FALSE)
cat("LSI1 vs depth correlation:",
    round(cor(Embeddings(mg, "lsi")[, 1], log10(mg$nCount_peaks_mg + 1)), 3), "\n")
# executed run: |r| ~ 1.0, so LSI1 is excluded

mg <- RunHarmony(mg, group.by.vars = "donor_id", reduction = "lsi",
                 dims.use = CONFIG$lsi_dims, reduction.save = "harmony_atac",
                 project.dim = FALSE, verbose = FALSE)

DefaultAssay(mg) <- "SCT"
mg <- FindMultiModalNeighbors(
  mg,
  reduction.list       = list("harmony", "harmony_atac"),
  dims.list            = list(1:CONFIG$N_PCS, seq_along(CONFIG$lsi_dims)),
  modality.weight.name = c("RNA.weight", "ATAC.weight"),
  verbose              = FALSE
)
print(tapply(mg$ATAC.weight, mg$atac_failed_donor, summary))

mg <- RunUMAP(mg, nn.name = "weighted.nn", reduction.name = "umap_wnn",
              reduction.key = "wnnUMAP_", seed.use = CONFIG$seed, verbose = FALSE)
mg <- FindClusters(mg, graph.name = "wsnn", resolution = CONFIG$final_res,
                   algorithm = 4, random.seed = CONFIG$seed, verbose = FALSE)
mg$wnn_cluster <- Idents(mg)
Idents(mg) <- "leiden_final"
cat("WNN vs RNA-only ARI:",
    round(adjustedRandIndex(mg$wnn_cluster, mg$leiden_final), 3), "\n")

# -----------------------------------------------------------------------------
# 13. Donor and SA representation per cluster
# -----------------------------------------------------------------------------
rep_tbl <- mg@meta.data %>%
  group_by(leiden_final) %>%
  summarise(n_cells     = n(),
            n_donors    = n_distinct(donor_id),
            n_SA_donors = n_distinct(donor_id[Group == "SA"]),
            max_donor_frac = max(table(donor_id)) / n(),
            .groups = "drop")
print(rep_tbl)
write.csv(rep_tbl, p_out("tables", "cluster_donor_representation.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 14. Compartment subclustering (homeostatic and DAM / lipid compartments)
# TODO: confirm the compartments and the final fine-state count from the
# executed run. The May run used res = 0.2 within each compartment and
# logfc.threshold = 0.1 for markers.
# -----------------------------------------------------------------------------
subcluster_compartment <- function(obj, cells, prefix, res = 0.2) {
  sub <- subset(obj, cells = cells)
  sub <- FindNeighbors(sub, reduction = "harmony", dims = 1:CONFIG$N_PCS,
                       k.param = CONFIG$k_param, verbose = FALSE)
  sub <- FindClusters(sub, resolution = res, algorithm = 4,
                      random.seed = CONFIG$seed, verbose = FALSE)
  setNames(paste0(prefix, "-", Idents(sub)), colnames(sub))
}
homeo_clusters <- character(0)  # TODO: clusters whose top_state is MG0_Homeostatic
dam_clusters   <- character(0)  # TODO: lipid / DAM compartment clusters (e.g. MG4)
mg$wp1_fine_label <- paste0("C", mg$leiden_final)
if (length(homeo_clusters)) {
  lab <- subcluster_compartment(mg, colnames(mg)[mg$leiden_final %in% homeo_clusters], "HOMEO")
  mg$wp1_fine_label[match(names(lab), colnames(mg))] <- lab
}
if (length(dam_clusters)) {
  lab <- subcluster_compartment(mg, colnames(mg)[mg$leiden_final %in% dam_clusters], "DAM")
  mg$wp1_fine_label[match(names(lab), colnames(mg))] <- lab
}

# -----------------------------------------------------------------------------
# 15. Milo neighbourhood graph (testing is WP2)
# -----------------------------------------------------------------------------
sce <- as.SingleCellExperiment(mg, assay = "SCT")
reducedDim(sce, "HARMONY") <- Embeddings(mg, "harmony")[, 1:CONFIG$N_PCS]
milo <- Milo(sce)
milo <- buildGraph(milo, k = CONFIG$k_param, d = CONFIG$N_PCS, reduced.dim = "HARMONY")
milo <- makeNhoods(milo, prop = CONFIG$milo_prop, k = CONFIG$k_param,
                   d = CONFIG$N_PCS, refined = TRUE, reduced_dims = "HARMONY")
milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), samples = "donor_id")
nh <- colSums(nhoods(milo))
cat("Neighbourhoods:", ncol(nhoods(milo)), "| median size:", median(nh), "\n")
# annotateNhoods() needs da.res from testNhoods(), so it is deferred to WP2.
saveRDS(milo, p_out("objects", "milo_wp1.rds"))

# -----------------------------------------------------------------------------
# 16. Freeze
# -----------------------------------------------------------------------------
mg@misc$wp1_frozen <- list(
  date = as.character(Sys.Date()), input_md5 = CONFIG$wp0_md5,
  n_cells = ncol(mg), n_donors = length(unique(mg$donor_id)),
  sct_vars_to_regress = "percent.MT", harmony_group_by = "donor_id",
  N_PCS = CONFIG$N_PCS, k = CONFIG$k_param, final_res = CONFIG$final_res,
  n_clusters = nlevels(mg$leiden_final), excluded_prelim_clusters = exclude_clusters,
  milo_nhoods = ncol(nhoods(milo)), milo_median_nhood = median(nh),
  note = "Reconstructed pipeline; see header."
)
DefaultAssay(mg) <- "SCT"
frozen <- p_out("objects", "WP1_FROZEN.rds")
saveRDS(mg, frozen)
md5 <- unname(tools::md5sum(frozen))
writeLines(paste0(md5, "  WP1_FROZEN.rds"), paste0(frozen, ".md5"))
Sys.chmod(frozen, mode = "0444")
cat("WP1 frozen:", ncol(mg), "cells |", nlevels(mg$leiden_final), "clusters | MD5", md5, "\n")
writeLines(capture.output(sessionInfo()), p_out("tables", "sessionInfo_WP1.txt"))
