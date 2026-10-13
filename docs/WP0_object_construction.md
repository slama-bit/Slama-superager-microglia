# WP0 — Artifact Firewall & Microglia Object Construction

**Status:** ✅ Complete, frozen
**Notebook:** `notebooks/WP0_clean_v2_executed.ipynb`

## Output

| Field | Value |
|---|---|
| File | `microglia_WP0_frozen.rds` (884 MB, read-only) |
| MD5 | `83cefa0c09ece745eb4dc46be79304f7` |
| Version | WP0_v2.1 |
| Nuclei | 10,676 |
| Donors | 38 |
| Assays | `RNA`, `Peaks` (deposited), `peaks_mg` (re-called) |
| Metadata added | `donor_QC_PC1`–`PC5`, clinical metadata, `atac_failed_donor` |

## Pipeline

1. **Extraction.** Took nuclei with `Cluster == "Microglia"` from the deposited object (12,165 nuclei). Corrected an earlier gate that wrongly included OPCs. Screened for BAM markers.
2. **Relabeling.** Renamed `Group`: `MCI` → `PCI`.
3. **Low-yield donor removal.** Removed donor 32, which had fewer than 20 microglia nuclei.
4. **QC.** Applied within-donor 3-MAD filters on RNA and ATAC metrics.
   - No MT filter was applied. A MAD ceiling of 2.2% would have removed 48% of PCI nuclei, creating a fatal group bias.
5. **Doublets.** Removed nuclei with `doublet_score > 0.25`.
6. **Identity filter (per cell).** Used detection fractions, following the Sun et al. 2023 approach. A nucleus was removed only if it had zero microglial markers **and** two or more contamination markers.
   - Identity panel, using Ensembl IDs: TMEM119, C1QB, P2RY12, CX3CR1, CSF1R, HEXB.
   - C3 was removed from the panel because it had 0% detection and is also expressed by astrocytes.
7. **Peak re-calling.** Ran MACS2 pooled across all microglia (`group.by = NULL`), giving **82,250 peaks**.
   - Donor 20 (SA) was excluded because its ATAC library failed (TSS = 0). This donor is retained for RNA.
8. **Peak matrix.** Rebuilt the `peaks_mg` assay with `FeatureMatrix` against the new peak set.
9. **Donor QC PCA.** Computed PC1–5 on donor-level QC summaries, for use as covariates in downstream models.
10. **Freeze.** Saved the RDS, set it read-only (`chmod 444`), and recorded the MD5.

## Per-group yield

| YA | HA | PCI | AD | SA |
|---|---|---|---|---|
| 902 | 3,172 | 1,363 | 3,537 | 1,702 |

## Known limitations

- FRiP is not computable, because only a single merged fragment file is available.
- Donor 20 has RNA only, so it is excluded from all ATAC analyses.
- Ambient RNA correction (DecontX) was planned. [CONFIRM: applied or dropped]
