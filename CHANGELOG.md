# Changelog

Version numbers match `VERSION` in `app.R`.

## Unreleased

Changes since the first public release. `VERSION` remains `0.3.4`.

### Reference data (2026-10-13)

`data/Final_Taxa_Emergence_Timeline_with_Cellular_Organisms.csv`:

- Cellular organisms: 3500 → 4200 MYA.
- Taxon spellings corrected: `Haplorhini` → `Haplorrhini`, `Homoninae` → `Homininae`.

These changes affect only the MYA shown in stratum labels. They do not change stratum assignments or any enrichment calculation (Fisher's exact tests, FDR correction, GO enrichment). The spelling corrections have no runtime effect because `app.R` already mapped the old spellings to the correct taxon names.

### Documentation

- Added a citation for the associated manuscript, a README screenshot, and an `renv.lock` recording package versions.
- Clarified isoform handling, the GO background, and the GO p-values shown in the plots.
- Moved version history to this file; fixed the live-app link and repository layout.

## 0.3.4 — first public release (2026-06-30)

- Public deployment-ready. Google Sheets read via `gs4_deauth()`; cached RDS fallback.

## 0.3.3

- Replaced fgsea (ranked GSEA, not appropriate for ORA) with topGO Fisher.
- GO p-value sliders now actually filter the plot.
- MYA annotations joined by taxon name rather than row order.
- Various small bug fixes.

## 0.3.2a

- Last pre-fix release.
