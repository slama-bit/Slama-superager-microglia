# Roadmap — WP2 to WP8

| WP | Method | Depends on | Notes |
|---|---|---|---|
| **WP2** | Differential abundance with Milo (neighbourhoods) and propeller (state proportions); primary contrasts C1–C4, SA vs HA primary | WP1 freeze | Biobank and donor QC PCs as covariates. Resolve the open SA/biobank decision first |
| **WP3** | Pseudobulk DE per state with edgeR (Model A formula) | WP1 | Minimum nuclei per donor per state; donor 20 flagged |
| **WP4** | ATAC pseudobulk differential accessibility on the `peaks_mg` assay | WP1 | Excludes donor 20 |
| **WP5** | SCENIC+ regulatory networks | WP4 stable | Python environment |
| **WP6** | Trajectory with PHATE / pseudotime | WP1, WP3 | |
| **WP7** | CellChat / LIANA | Full parent atlas | Needs non-microglial cell types from the parent dataset |
| **WP8** | Sensitivity: donor jackknife, old-only models, negative controls | Runs alongside WP3–WP6 | |

## Evidence tiers for SA claims

1. **Tier 1:** a discrete SA-specific state.
2. **Tier 2:** an SA program distributed across shared states.
3. **Tier 3:** a preserved-young / anti-AD configuration in SA.
