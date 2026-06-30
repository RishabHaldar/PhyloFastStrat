# PhyloFastStrat

An interactive Shiny app for exploring the evolutionary origins of disease-associated genes. Maps any set of human genes onto a phylostratigraphic timeline, tests for enrichment in each phylostratum, and runs GO biological process over-representation analysis at each evolutionary stratum.

**Live app:** _[paste your shinyapps.io URL here once deployed]_

---

## For users

You don't need to install anything. Open the link above in any web browser — Chrome, Safari, Firefox, Edge all work. The app also runs on phones and tablets, though the plots are easier to read on a laptop.

**First load is slow.** The app sleeps when nobody's using it, so the first person to open it after a quiet stretch waits ~10–20 seconds while it wakes up. Subsequent loads are instant. If you see a blank screen, give it a moment.

### How to use it

1. **Pick a gene list.** The default mode ("Preset Gene Lists") gives you a dropdown of curated disease gene sets. Flip the "Run Mode" toggle to "Custom Gene Lists" to either paste in your own genes or upload a one-column CSV.
2. **Click "Run analysis on selected genes."** The four tabs along the top will populate:
   - *Origin of genes* — a line plot showing how many genes from your list emerged at each evolutionary stratum, compared to the rest of the human genome. Click a point on the plot to see exactly which genes emerged at that stratum.
   - *Rate comparison* — odds-ratio enrichment of your gene list per stratum. Asterisks mark FDR-adjusted Fisher p-values < 0.05.
   - *Results table* — sortable, exportable table of every gene in your list and the stratum where it appeared. Use the buttons at the top to download as CSV or Excel.
   - *Functional enrichment plots* — GO biological process over-representation in your gene list, split by stratum. Click a stratum's column to drill into its top terms.
3. **Download anything you want to keep.** The Results table has CSV/Excel/copy buttons; the GO plot has a "Download GO Table" button.

### If it breaks

Refresh the page (Cmd/Ctrl + R). 90% of the time that fixes it. If it still doesn't work, [open an issue on GitHub](../../issues) or contact the maintainer.

---

## What PhyloFastStrat does

The pipeline behind the app, briefly:

- Every human protein is assigned a **phylostratum** — the most recent common ancestor in which a homolog can be detected — using [`phylostratr`](https://github.com/arendsee/phylostratr) with BLAST against a panel of representative species. Detection uses an adjusted p-value cutoff of 0.001.
- For a given user gene list, the app counts how many genes fall in each stratum and tests enrichment per stratum with **Fisher's exact test** on a 2×2 contingency table (disease vs. complement × in-stratum vs. out-of-stratum), with Benjamini-Hochberg FDR correction across strata.
- GO biological process **over-representation** is run separately within each stratum using `topGO`'s classic Fisher algorithm against the full GO-annotated human gene universe. Fold enrichment is reported as observed / expected gene counts.

---

## For developers

### Run locally

```r
# Install dependencies (one-time)
install.packages("BiocManager")
BiocManager::install(c(
  "shiny", "shinyjs", "shinycssloaders", "shinyWidgets", "shinythemes",
  "readxl", "DT", "tidyverse", "ggtext", "ggrepel", "ggpubr",
  "googlesheets4",
  "org.Hs.eg.db", "topGO", "GO.db", "AnnotationDbi", "STRINGdb",
  "phylostratr"
))

# Clone and run
# git clone https://github.com/YOUR_USERNAME/PhyloFastStrat.git
# cd PhyloFastStrat
# Open app.R in RStudio and click "Run App"
```

### Repository layout

```
PhyloFastStrat/
├── app.R                              # Shiny app
├── data/
│   ├── ph_0.001.rds                   # Phylostratigraphy result table
│   ├── uniprot_gene_conv_df_updated.rds  # UniProt ↔ gene symbol map
│   ├── Final_Taxa_Emergence_Timeline_with_Cellular_Organisms.csv  # MYA per stratum
│   └── MasterGeneLists_fallback.rds   # (optional) cached gene list snapshot
├── deploy.R                           # One-line redeploy helper
├── README.md
└── .gitignore
```

The disease gene lists themselves are **not stored in the repo**. They live in a public Google Sheet, and the app reads them at startup via `googlesheets4::gs4_deauth()`. To use your own gene lists, either edit the master sheet (if you have access) or change `MASTER_SHEET_URL` near the top of `app.R` to point at your own sheet (which must be set to "Anyone with the link: Viewer").

### Deploy your own instance to shinyapps.io

1. Make a free account at https://www.shinyapps.io
2. Get your account token from Account → Tokens, and run the `rsconnect::setAccountInfo(...)` snippet it gives you.
3. From R, with your working directory set to the repo root:
   ```r
   source("deploy.R")
   ```
   The first deploy installs all Bioconductor dependencies on the shinyapps.io server and can take 20–30 minutes. Subsequent deploys are usually under 2 minutes.

### Updating the gene list

The master gene list lives at the Google Sheet URL hardcoded in `app.R` (`MASTER_SHEET_URL`). Edit the sheet, and the app picks up the changes the next time it starts. (On shinyapps.io free tier the app sleeps every ~5 min, so changes propagate fast.)

To keep a fallback snapshot in the repo (so the app still works if the sheet ever goes down), run interactively:

```r
library(googlesheets4)
gs4_deauth()
sheet <- read_sheet("YOUR_SHEET_URL")
saveRDS(sheet, "data/MasterGeneLists_fallback.rds")
```

Then commit and redeploy.

---

## Version history

- **0.3.4** — Public deployment-ready. Google Sheets read via `gs4_deauth()`; cached RDS fallback.
- **0.3.3** — Replaced fgsea (ranked GSEA, not appropriate for ORA) with topGO Fisher. GO p-value sliders now actually filter the plot. MYA annotations joined by taxon name rather than row order. Various small bug fixes.
- **0.3.2a** — Last pre-fix release.

## Citation

If you use PhyloFastStrat in published work, please cite [_paper / preprint citation pending_]. The underlying phylostratigraphy method is from Domazet-Lošo & Tautz (2008) and uses Arendsee's `phylostratr` R package.

## License

MIT — see [LICENSE](LICENSE).
