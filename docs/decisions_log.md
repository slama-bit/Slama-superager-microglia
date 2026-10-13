# Methodological Decisions Log

Locked decisions, each with its rationale. Each was challenged against primary literature or published code before being accepted.

| # | Decision | Rationale |
|---|---|---|
| D1 | Donor is the unit of inference | Treating cells as replicates inflates false positives; use pseudobulk or donor-level models |
| D2 | `MCI` relabelled to `PCI` everywhere | Deposited label differs from the study nomenclature. A mismatch previously dropped 6 donors from edgeR |
| D3 | No MT% filter in WP0 | A MAD ceiling of 2.2% would have removed 48% of PCI nuclei, creating group bias |
| D4 | Identity filter applied per cell using detection fractions, not by cluster | Cluster numbers are non-deterministic across runs. AddModuleScore outputs are centered at zero and do not work with ratio rules |
| D5 | C3 removed from the identity panel; HEXB added | C3 had 0% detection and is expressed by reactive astrocytes |
| D6 | Pooled MACS2 peak call (`group.by = NULL`) | Peak sets stratified by group would bias WP4 differential accessibility |
| D7 | Donor 20 excluded from ATAC only | Failed ATAC library (TSS = 0); its RNA is valid |
| D8 | SCTransform `vars.to.regress = "percent.MT"` only | Regressing biobank would remove SA signal, because SA and NU_ADRC are perfectly collinear |
| D9 | Harmony corrects on `donor_id` only | Same collinearity. Diagnosis must never enter integration |
| D10 | Biobank is a covariate only in WP2–WP8 models | It is the only place it can appear. SA effects there remain non-identifiable from site effects |
| D11 | Sun 2023 states assigned by z-scored UCell enrichment | Raw UCell scores are not comparable across modules |
| D12 | FLT1 excluded from the endothelial contamination module | Documented expression in macrophages |
| D13 | WNN excludes LSI dim 1 | Correlates with sequencing depth (r ≈ 1.0) |

## Open decision (before WP2)

**How SA is handled given the biobank confound.**

- **Option (a):** keep SA contrasts and report them as confounded or exploratory.
- **Option (b):** make contrasts within UW_BRaIN (YA/HA/PCI/AD) primary, and treat SA as descriptive.
