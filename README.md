# SuperAger Hippocampal Microglia — Single-Nucleus Multiome Analysis

Independent re-analysis of hippocampal microglia across cognitive aging and Alzheimer's disease using the 10x Multiome (snRNA-seq + snATAC-seq) dataset **GSE268609** (Disouky et al., *Nature* 2026, DOI: [10.1038/s41586-026-10169-4](https://doi.org/10.1038/s41586-026-10169-4)).

## Scientific question

Do SuperAgers (SA) harbor:

- **Outcome A:** a discrete SA-specific microglial state,
- **Outcome B:** an SA-specific transcriptional program within shared states, or
- **Outcome C:** a preserved-young / anti-AD microglial configuration?

The primary contrast is **SA vs HA** (healthy aged). The **donor** is the unit of inference throughout.

## Cohort

| Group | Description | Nuclei (WP0 frozen) |
|---|---|---|
| YA | Young adults | 902 |
| HA | Healthy aged | 3,172 |
| PCI | Pre-clinical impairment (deposited as `MCI`, renamed) | 1,363 |
| AD | Alzheimer's disease | 3,537 |
| SA | SuperAgers | 1,702 |
| **Total** | 38 donors | **10,676** |

## Status

| WP | Scope | Status |
|---|---|---|
| WP0 | Artifact firewall and microglia object construction | ✅ Complete, frozen |
| WP1 | RNA-first state discovery, annotation, multimodal validation | ✅ Complete, frozen |
| WP2 | Differential abundance (Milo + propeller) | ⏳ Next |
| WP3 | Pseudobulk DE per state (edgeR) | Pending |
| WP4 | ATAC pseudobulk differential accessibility | Pending |
| WP5 | Regulatory networks (SCENIC+) | Pending |
| WP6 | Trajectory (PHATE / pseudotime) | Pending |
| WP7 | Cell–cell communication (CellChat / LIANA) | Pending |
| WP8 | Sensitivity package (donor jackknife, old-only models, negative controls) | Pending |

## Repository layout

```
.
├── README.md
├── docs/
│   ├── WP0_object_construction.md
│   ├── WP1_state_discovery.md
│   ├── decisions_log.md
│   └── roadmap.md
├── notebooks/        # executed .ipynb per WP (papermill outputs)
├── scripts/          # reconstructed WP0 / WP1 pipelines (see Code status)
├── results/
│   ├── tables/       # small CSV outputs (QC tables, marker lists, abundance tests)
│   └── plots/        # PNG / PDF figures
└── env/              # sessionInfo(), conda env export
```

## Code status — read before using

> **The scripts in `scripts/` are reconstructed, not the executed code.** The compute instance holding the original notebooks was lost. The scripts were rebuilt from the analysis session logs. They follow the final locked pipeline, but re-running them is **not guaranteed to reproduce the frozen objects or their MD5s**. Values from the executed run are kept as `# executed run:` comments. Steps that could not be recovered are marked `# TODO`.

| Script | Source of reconstruction |
|---|---|
| `scripts/WP0_object_construction.R` | Full WP0_clean_v2 notebook source, plus the v2.1 identity-panel fix |
| `scripts/WP1_state_discovery.R` | Final WP1 pipeline steps, plus the Harmony `donor_id` fix |

Large objects (`.rds`, fragment files, peak matrices) are **not** versioned. See `.gitignore`. Each frozen object is identified by its MD5 checksum in the docs.

## Environment

- Vertex AI JupyterLab (Google Cloud, Debian 11), 64 GB RAM
- R 4.5.3, Seurat v5.5.0, Signac 1.17.1, Harmony, UCell, edgeR, miloR
- Python: scVI (sensitivity integration), SCENIC+ (planned WP5)
- MACS2, EnsDb.Hsapiens.v86 (hg38)
- Notebook execution via Papermill

## Key caveat

All SA donors come from one biobank (NU_ADRC); all other groups come from UW_BRaIN. **SA status and biobank are perfectly collinear.** Biobank is therefore never regressed during normalization or integration (that would remove SA biology). Any SA-specific effect in WP2–WP8 cannot be separated from site effects and is interpreted accordingly. See [`docs/decisions_log.md`](docs/decisions_log.md).

## Documentation

- [WP0 — object construction](docs/WP0_object_construction.md)
- [WP1 — state discovery](docs/WP1_state_discovery.md)
- [Decision log](docs/decisions_log.md)
- [Roadmap WP2–WP8](docs/roadmap.md)

## Citation

If you use this analysis, cite the original dataset: Disouky et al., *Nature* (2026), GEO accession GSE268609.

## Author

Islam Elkhosht
