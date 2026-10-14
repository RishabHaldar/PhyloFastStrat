# PhyloFastStrat

PhyloFastStrat is an R Shiny application for analyzing the evolutionary origins of human gene sets through phylostratigraphy, statistical enrichment, and Gene Ontology analysis.

[Launch the application](https://rishabhaldar.shinyapps.io/PhyloFastStrat/)

![PhyloFastStrat app showing gene emergence across phylostrata](images/app_screenshot.png)

## Features

- Visualize the distribution of genes across evolutionary phylostrata.
- Evaluate enrichment within each phylostratum using Fisher's exact test and Benjamini–Hochberg FDR correction.
- Perform Gene Ontology (GO) biological-process enrichment analysis.
- Explore curated disease-associated gene sets or upload custom gene lists.
- Download gene assignments and enrichment results.

## Methods

PhyloFastStrat uses precomputed phylostratigraphic assignments generated with [`phylostratr`](https://github.com/arendsee/phylostratr). Evolutionary ages are displayed using approximate divergence-time estimates.

For each phylostratum, Fisher's exact test evaluates enrichment of the selected gene set against the remaining human gene background. P-values are adjusted across phylostrata using the Benjamini–Hochberg procedure. By default, the per-stratum test counts UniProt entries; selecting **Exclude alternative isoforms** enables gene-level counting.

The application retains alternative UniProt entries by default. Selecting **Exclude alternative isoforms** assigns each mapped gene a single phylostratum based on its earliest-emerging entry. This setting was used for the gene-level analyses in the associated atherosclerosis study.

Gene Ontology biological-process enrichment is performed using `topGO` with its classic Fisher algorithm and the annotated human gene background. Fold enrichment is calculated as the ratio of observed to expected gene counts. The main GO visualization displays unadjusted p-values; adjusted results are available in the detailed output.

## Data and inputs

PhyloFastStrat uses precomputed evolutionary assignments and a UniProt-to-gene-symbol mapping table. Curated disease gene lists are retrieved from a public Google Sheet when the application starts.

Custom gene lists can be supplied through the application. Gene-symbol matching is case-sensitive and does not resolve aliases. Symbols containing `RNU` are excluded during input processing.

## Installation

The repository includes an `renv.lock` file created with R 4.4.2 and Bioconductor 3.20.

```bash
git clone https://github.com/RishabHaldar/PhyloFastStrat.git
cd PhyloFastStrat
```

From the repository directory, start R and run:

```r
install.packages("renv")
renv::restore()
shiny::runApp(".")
```

If `renv::restore()` asks whether to activate the project, choose to activate it. Activation may restart R or RStudio; after the restart, run `shiny::runApp(".")`.

The application requires access to the configured Google Sheet. An offline fallback gene-list file is not distributed with this repository.

## Version history

See [CHANGELOG.md](CHANGELOG.md).

## Citation

Natterson-Horowitz B, Haldar R, Sawyer D. Evolutionary Origins of Atherosclerosis: Ancient Conserved Biology Underlies Disease Vulnerability and Constrains Therapeutic Modification. Manuscript under review, 2026.

## License

MIT License. See [LICENSE](LICENSE).
