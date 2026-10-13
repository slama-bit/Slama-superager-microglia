# WP1 — RNA-First State Discovery, Annotation & Validation

**Status:** ✅ Complete, frozen
**Input:** WP0 frozen object (MD5 `83cefa0c09ece745eb4dc46be79304f7`)
**Output:** `WP1_FROZEN.rds`, MD5 [CONFIRM]

> Values marked [CONFIRM] are ones that changed between runs. Fill them in from the final frozen object before publishing.

## Pipeline

1. **Normalization.** SCTransform v2 with 3,000 HVGs.
   - `vars.to.regress = "percent.MT"` only.
   - Biobank is not regressed. See the decision log.
2. **Dimensionality reduction.** PCA, using `N_PCS = 25` downstream. [CONFIRM]
3. **Integration.** Harmony with `group.by.vars = "donor_id"` only.
   - Diagnosis and biobank are not used for integration.
4. **Clustering.** Built a kNN graph (k = 20) and ran Leiden over a resolution grid of 0.2–0.8.
   - Selected resolution 0.4, giving [CONFIRM: 14] clusters.
   - Selection criteria: clustree (Zappia & Oshlack 2018), silhouette, and Random Forest stability (90.2%; Olah 2020 criterion).
5. **Contamination audit.** Applied the HuMicA three-criterion exclusion rule. The contamination module lists were audited against the literature:
   - FLT1 was dropped from the endothelial module, because macrophages express it.
   - COLEC12 was dropped, because it has no validation in human snRNA-seq.
   - Cluster 3 was retained, because co-expression analysis showed it is microglial.
6. **Annotation.** Ran UCell scoring against the Sun et al. 2023 human microglia state atlas.
   - States were assigned using z-scored relative enrichment, not raw UCell scores.
   - MG9 does not exist in the published atlas.
7. **Subclustering.** Subclustered the homeostatic and DAM compartments to produce fine-grained state labels (`wp1_fine_label`). [CONFIRM final count]
8. **Multimodal validation.** Built a WNN graph from RNA and ATAC (LSI with dim 1 excluded, since it correlates with depth).
   - Donor 20 was handled as NA in the ATAC modality.
9. **Sensitivity checks.**
   - scVI integration and unintegrated PCA were compared by ARI against the Harmony clustering.
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
