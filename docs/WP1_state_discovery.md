# WP1 — RNA-First State Discovery, Annotation & Validation

**Status:** ✅ Complete, frozen
**Input:** WP0 frozen object (MD5 `83cefa0c09ece745eb4dc46be79304f7`)
**Output:** `WP1_FROZEN.rds`, 9,311 nuclei, 38 donors. MD5 not recovered.
**Script:** `scripts/WP1_state_discovery.R` (reconstructed; the executed notebook was lost)

> Some values could not be recovered from the logs: the final MD5, the subclustering compartments and the fine-state count. They are marked `# TODO` in the script.

## Pipeline

0. **Pre-filter and contamination exclusion.**
   - Kept nuclei with WP0 `contam_score <= 0.4`, leaving 10,488.
   - Made a screening SCTransform/PCA and preliminary clusters (res 0.3, 12 clusters).
   - Scored UCell contamination modules.
   - Applied the HuMicA three-criterion rule. It excluded cluster 5 (oligodendrocyte, ratio 2.79) and cluster 11 (BAM, ratio 5.50).
   - **Result: 9,311 nuclei** (YA 629, HA 2,822, PCI 1,107, AD 3,218, SA 1,535).
1. **Normalization.** SCTransform v2 with 3,000 HVGs.
   - `vars.to.regress = "percent.MT"` only.
   - Biobank is not regressed. See the decision log.
2. **Dimensionality reduction.** PCA on 50 PCs, with `N_PCS = 25` used downstream (78.5% cumulative variance).
3. **Integration.** Harmony with `group.by.vars = "donor_id"` only.
   - Diagnosis and biobank are not used for integration.
4. **Clustering.** Built a kNN graph (k = 20) and ran Leiden over a resolution grid of 0.2–0.8.
   - Selected resolution 0.4, giving 14 clusters. The grid gave 0.2 = 11, 0.4 = 14, 0.5 = 15, 0.6 = 16 and 0.8 = 19 clusters.
   - Selection criteria: clustree (Zappia & Oshlack 2018), silhouette, and Random Forest stability (Olah 2020 criterion, > 75%): res 0.4 = 90.2% and res 0.5 = 89.9%.
5. **Contamination re-audit at res 0.4.** No cluster met the full rule. The contamination module lists were audited against the literature:
   - FLT1 was dropped from the endothelial module, because macrophages express it.
   - COLEC12 was dropped, because it has no validation in human snRNA-seq.
   - Suspicious res-0.6 clusters 14 and 15 were retained, because they had normal microglial marker expression.
6. **Annotation.** Ran UCell scoring against the Sun et al. 2023 human microglia state atlas.
   - States were assigned using z-scored relative enrichment, not raw UCell scores.
   - MG9 does not exist in the published atlas.
7. **Subclustering.** Subclustered the homeostatic and DAM compartments to produce fine-grained state labels (`wp1_fine_label`). The final count was not recovered.
8. **Multimodal validation.** Built a WNN graph from RNA and ATAC (LSI with dim 1 excluded, since it correlates with depth).
   - ATAC Harmony used `donor_id` on LSI dims 2–30.
   - Donor 20 has no usable ATAC signal; its ATAC weights are checked.
9. **Sensitivity checks.**
   - scVI integration and unintegrated PCA were compared by ARI against the Harmony clustering. This was done in an earlier run; it is not in the reconstructed script, because the scVI latent matrix was lost.
   - Moderate disagreement between methods is expected (Antonsson & Melsted, *Genome Res* 2025).
10. **Milo graph.** Built the neighbourhood graph in preparation for WP2.

## Interim observations (not inferential claims)

- Abundance screen with edgeR (donor-level, global BH): only MG-NEU depletion in AD survived FDR < 0.05. **No SA contrast survived.**
- Candidate states enriched in SA and HA (for example a DAM-SPP1 sub-state) are noted for testing in WP2 and WP3. They are not established findings.

## Engineering notes

- SCTransform hung under Papermill with the `future` backend. Fixed with `conserve.memory = TRUE` and `BiocParallel::SerialParam()`.
- An edgeR group factor bug (`MCI` vs `PCI`) dropped 6 donors and was fixed.
- Seurat v5 `Assay5` requires `LayerData()` syntax throughout.
- Gene rownames are Ensembl IDs. Symbol map: `gene_map.csv` (`ensembl_gene_id`, `hgnc_symbol`).
