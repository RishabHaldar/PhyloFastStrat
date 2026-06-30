# PhyloFastStrat
VERSION <- "0.3.4"

# ---------------------------------------------------------------------------
# Changes vs 0.3.3:
#   - MasterGeneLists is read from Google Sheets at app launch using
#     gs4_deauth() (no OAuth). REQUIRES the sheet to be set to
#     "Anyone with the link can view" in Google Sheets sharing settings.
#   - Sheet read wrapped in tryCatch with a cached RDS fallback so the app
#     still starts if the sheet is temporarily unreachable.
#   - All file paths now relative; ready for shinyapps.io deployment.
# ---------------------------------------------------------------------------

library(BiocManager)
options(repos = BiocManager::repositories())
library(shiny)
library(shinyjs)
library(shinycssloaders)
library(shinyWidgets)
library(readxl)
library(org.Hs.eg.db)
library(phylostratr)
library(tidyverse)
library(topGO)
library(GO.db)
library(AnnotationDbi)
library(STRINGdb)
library(DT)
library(shinythemes)
library(ggtext)
library(googlesheets4)

##### Configuration #####
MASTER_SHEET_URL <- "https://docs.google.com/spreadsheets/d/1oEktKQNY8kaEFpImEW-wnoBz4dqLa3eWgrji23AFU4U/edit?gid=1728211307#gid=1728211307"

##### Load gene lists from public Google Sheet #####
# This requires the sheet to be set to "Anyone with the link: Viewer" in
# Google Sheets. If the sheet is private, the read will fail and the app
# will fall back to a cached RDS at data/MasterGeneLists_fallback.rds.
MasterGeneLists <- tryCatch({
  googlesheets4::gs4_deauth()
  sheet <- googlesheets4::read_sheet(MASTER_SHEET_URL)
  message("Loaded MasterGeneLists from Google Sheet (", ncol(sheet), " lists, ",
          nrow(sheet), " rows).")
  sheet
}, error = function(e) {
  message("Could not read Google Sheet: ", conditionMessage(e))
  fallback_path <- "data/MasterGeneLists_fallback.rds"
  if (file.exists(fallback_path)) {
    message("Loading fallback gene lists from ", fallback_path)
    readRDS(fallback_path)
  } else {
    stop("Could not load gene lists from Google Sheet, and no fallback ",
         "RDS found at ", fallback_path, ". Make sure the sheet is set to ",
         "'Anyone with the link: Viewer', or place a fallback file.")
  }
})

##### Load static data #####
ph <- readRDS("data/ph_0.001.rds")
levels(ph$mrca_name)[match("cellular organisms", levels(ph$mrca_name))] <- "Cellular Organisms"

emergence_df <- read.csv("data/Final_Taxa_Emergence_Timeline_with_Cellular_Organisms.csv")
# Reconcile minor naming inconsistencies between CSV and ph$mrca_name levels.
taxa_fixes <- c(
  "Cellular organisms" = "Cellular Organisms",
  "Haplorhini"         = "Haplorrhini",
  "Homoninae"          = "Homininae"
)
idx <- match(emergence_df$Taxa, names(taxa_fixes))
emergence_df$Taxa[!is.na(idx)] <- taxa_fixes[idx[!is.na(idx)]]

stopifnot("Taxa" %in% colnames(emergence_df))
mya_lookup <- emergence_df$Emergence..MYA.[match(levels(ph$mrca_name), emergence_df$Taxa)]
if (any(is.na(mya_lookup))) {
  warning("Some MRCA names did not match a row in the emergence timeline CSV: ",
          paste(levels(ph$mrca_name)[is.na(mya_lookup)], collapse = ", "))
}
levels(ph$mrca_name) <- paste0(levels(ph$mrca_name), ": ", mya_lookup, " MYA")

gene_uniprot_conv_df <- readRDS("data/uniprot_gene_conv_df_updated.rds")

ph <- merge(ph, gene_uniprot_conv_df, by.x = "qseqid", by.y = "Entry", all.x = TRUE)
ph$Gene <- ph$`Gene Names (primary)`
ph$gene <- ph$Gene

all_genes <- unique(ph$Gene[!is.na(ph$Gene)])

ph_filt <- ph |>
  dplyr::filter(!is.na(Gene)) |>
  dplyr::mutate(mrca_name = factor(mrca_name, levels = levels(ph$mrca_name), ordered = TRUE)) |>
  dplyr::arrange(Gene, mrca_name) |>
  dplyr::distinct(Gene, .keep_all = TRUE)

# Pre-build the BP annotation map once (used by GO enrichment).
annFUN.org.symbol <- topGO::annFUN.org(whichOnto = "BP",
                                       mapping   = "org.Hs.eg.db",
                                       ID        = "symbol")

my_theme <- ggplot2::theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    panel.border     = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(color = "black"),
    axis.text        = element_text(color = "black", size = 12),
    axis.title       = element_text(color = "black", size = 12),
    plot.title       = element_text(face = "bold", size = 14)
  )

##### UI #####
ui <- fluidPage(
  theme = shinytheme("yeti"),
  titlePanel(paste0("PhyloFastStrat v", VERSION)),
  sidebarLayout(
    sidebarPanel(
      shinyWidgets::switchInput(
        inputId  = "runmode",
        label    = "Run Mode",
        offLabel = "Preset Gene Lists",
        onLabel  = "Custom Gene Lists",
        value    = FALSE
      ),
      conditionalPanel(
        condition = "!input.runmode",
        selectInput(
          inputId   = "disease_of_interest",
          label     = "Disease-Gene List:",
          choices   = c(colnames(MasterGeneLists)),
          multiple  = FALSE,
          selected  = NULL,
          selectize = TRUE
        )
      ),
      conditionalPanel(
        condition = "input.runmode",
        textInput(
          inputId     = "other_disease_of_interest",
          label       = "Gene set name:",
          value       = "My Gene List",
          placeholder = "Gene List"
        ),
        selectizeInput(
          inputId  = "selectized_genes",
          label    = "Genes to query:",
          choices  = NULL,
          multiple = TRUE,
          selected = NULL,
          options  = list(
            splitOn   = I("(function() { return /[,; ]/; })()"),
            plugins   = list("remove_button"),
            delimiter = " ",
            persist   = TRUE
          )
        ),
        fileInput(
          inputId     = "custom_file",
          label       = "Custom Gene List:",
          accept      = ".csv",
          buttonLabel = "Upload"
        )
      ),
      headerPanel(""),
      htmlOutput("numGenesMappedUP"),
      headerPanel(""),
      shinyWidgets::switchInput(
        inputId  = "no_isoforms",
        label    = "Exclude alternative isoforms",
        onLabel  = "Yes",
        offLabel = "No",
        value    = FALSE
      ),
      htmlOutput("numIsoformsMappedUP"),
      headerPanel(""),
      actionButton("run", "Run analysis on selected genes")
    ),
    mainPanel(tabsetPanel(
      type = "tabs",
      tabPanel("Origin of genes",
               plotOutput("ph_plot", click = "ph_plot_click"),
               DT::dataTableOutput("selected_ph", height = "100%")
      ),
      tabPanel("Rate comparison",
               plotOutput("rate_plot"),
               "For each phylostratum, enrichment is assessed using a 2x2 contingency table comparing disease genes versus the overall human genome inside and outside that phylostratum.\n
    The plotted odds ratio is calculated from this contingency table. A value greater than 1 indicates enrichment of disease genes in that phylostratum, while a value less than 1 indicates depletion.\n
    Fisher's exact test is used for each phylostratum, and p-values are adjusted across phylostrata using the Benjamini-Hochberg false discovery rate procedure.\n
    The presence of an asterisk ('*') indicates FDR-adjusted Fisher p-value < 0.05."
      ),
      tabPanel("Results table",
               shinyWidgets::switchInput(
                 inputId  = "show_uniprot_ids",
                 label    = "Show individual isoform IDs?",
                 onLabel  = "Yes",
                 offLabel = "No",
                 value    = FALSE
               ),
               DT::dataTableOutput("phylo_table", height = "100%")
      ),
      tabPanel("Functional enrichment plots",
               numericInput(inputId = "labeltop", label = "Number of top GO-BPs to label:",
                            value = 5, min = 1, max = 10),
               numericInput(inputId = "cappval", label = "Max -log P-value to visualize",
                            value = 100, min = 1, max = 100),
               sliderInput(inputId = "pval_range", label = "P-value range to visualize",
                           min = 0, max = 1, value = c(0, 0.05), step = 0.001,
                           pre = "P-value range: "),
               plotOutput("GOplot_static", height = "700px", click = "go_plot_click") %>%
                 withSpinner(),
               plotOutput("selected_go", height = "700px") %>% withSpinner(),
               br(),
               downloadButton("download_go_table", "Download GO Table (CSV)")
      )
    ))
  )
)

##### SERVER #####
server <- function(input, output, session) {
  updateSelectizeInput(session, inputId = "selectized_genes",
                       choices = all_genes, server = TRUE, selected = NULL)

  ##### genes_of_interest #####
  genes_of_interest_reactive <- eventReactive(input$run, {
    if (input$runmode && length(input$selectized_genes) > 0) {
      genes_of_interest <- input$selectized_genes
    } else if (!input$runmode && !is.null(input$disease_of_interest)) {
      genes_of_interest <- MasterGeneLists[[input$disease_of_interest]]
    } else {
      req(input$custom_file)
      custom_gene_list <- readr::read_csv(file = input$custom_file$datapath)
      genes_of_interest <- dplyr::pull(custom_gene_list, 1)
    }

    genes_of_interest_input <- genes_of_interest[!is.na(genes_of_interest)]

    genes_of_interest <- genes_of_interest[!is.na(genes_of_interest)]
    genes_of_interest <- genes_of_interest[!grepl("RNU", genes_of_interest)]
    genes_of_interest <- unique(genes_of_interest)

    genes_of_interest_found <- unique(ph$Gene[ph$Gene %in% genes_of_interest])

    output$numGenesMappedUP <- renderUI({
      HTML(paste0(
        length(genes_of_interest_input), " input genes; ",
        length(genes_of_interest), " gene IDs were valid.", "<br/>",
        length(genes_of_interest_found), " genes correctly mapped to UniProt IDs."
      ))
    })

    output$numIsoformsMappedUP <- renderUI({
      if (input$no_isoforms) {
        HTML(paste0(
          "Alternative isoforms were excluded. ",
          nrow(ph_filt[ph_filt$gene %in% genes_of_interest, ]),
          " proteins are shown. Isoform with earliest emergence shown."
        ))
      } else {
        HTML(paste0(
          "Alternative isoforms were included. ",
          nrow(ph[ph$Gene %in% genes_of_interest, ]),
          " protein isoforms are shown. Any isoform emergence shown."
        ))
      }
    })

    genes_of_interest
  })

  ##### disease_name #####
  disease_name <- reactive({
    if (input$runmode) input$other_disease_of_interest else input$disease_of_interest
  })

  ##### phylo_table #####
  output$phylo_table <- DT::renderDataTable({
    ph_use <- if (input$no_isoforms) ph_filt else ph

    if (input$show_uniprot_ids) {
      ph_use <- ph_use %>%
        filter(gene %in% genes_of_interest_reactive()) %>%
        dplyr::select(Gene = gene, Phylostrata = mrca_name, `Isoform (UniProt)` = qseqid) %>%
        mutate(Phylostrata = factor(Phylostrata, levels = levels(ph$mrca_name)))
    } else {
      ph_use <- ph_use %>%
        filter(gene %in% genes_of_interest_reactive()) %>%
        group_by(Gene = gene, Phylostrata = mrca_name) %>%
        summarize(`Total Isoforms` = n(), .groups = "drop") %>%
        mutate(Phylostrata = factor(Phylostrata, levels = levels(ph$mrca_name)))
    }

    ph_use
  },
  server = TRUE,
  extensions = "Buttons",
  options = list(
    paging       = TRUE,
    searching    = TRUE,
    fixedColumns = TRUE,
    scrollX      = TRUE,
    scrollY      = 400,
    autoWidth    = TRUE,
    ordering     = TRUE,
    lengthMenu   = c(10, 25, 50, 100),
    dom          = "Bfrtip",
    buttons      = c("copy", "csv", "excel", "pageLength")
  ))

  ##### ph_plot: emergence of disease vs other human genes #####
  output$ph_plot <- renderPlot({
    ph_local <- if (input$no_isoforms) ph_filt else ph

    genes_of_interest <- genes_of_interest_reactive()
    ph_local$is_of_interest <- ph_local$Gene %in% genes_of_interest
    ph_local$annot_of_interest <- ifelse(
      ph_local$is_of_interest,
      disease_name(),
      "Other human genes"
    )

    ph_by_ps <- ph_local %>%
      dplyr::distinct() %>%
      dplyr::mutate(
        mrca_name         = as.factor(mrca_name),
        annot_of_interest = as.factor(annot_of_interest)
      ) %>%
      dplyr::count(mrca_name, annot_of_interest, .drop = FALSE) %>%
      dplyr::mutate(gene_relevance = factor(annot_of_interest,
                                            levels = c(disease_name(), "Other human genes")))

    ggplot(ph_by_ps,
           aes(x = mrca_name, y = log10(n), group = gene_relevance,
               fill = gene_relevance, label = round(n, digits = 2))) +
      geom_line(aes(color = gene_relevance), linewidth = 1.2) +
      geom_label(label.size = 0.05, alpha = 0.3, size = 2.5, vjust = -0.25,
                 show.legend = FALSE, nudge_y = 0.1) +
      labs(
        y = bquote(Number ~ of ~ novel ~ genes ~ log[10] ~ scaled),
        x = "Phylostrata",
        fill = "", color = "",
        title = paste0(
          "Emergence of <span style='color:#E41A1C;'><strong>",
          disease_name(),
          "</strong></span> vs <span style='color:#377EB8;'><strong>Other human genes</strong></span><br>"
        )
      ) +
      scale_color_brewer(palette = "Set1") +
      scale_fill_brewer(palette = "Set1") +
      ylim(c(0, 5)) +
      ggpubr::labs_pubr() +
      my_theme +
      theme(
        axis.text.x         = element_text(angle = 90, hjust = 1),
        panel.background    = element_blank(),
        plot.title.position = "plot",
        plot.title          = ggtext::element_markdown(size = 12, lineheight = 1, hjust = 0.5),
        axis.text           = element_text(size = 10, face = "bold"),
        legend.position     = "none"
      )
  })

  ##### ph plot click #####
  output$selected_ph <- DT::renderDataTable({
    click <- input$ph_plot_click
    if (is.null(click) || is.null(click$x) || is.null(click$y)) return(NULL)

    x_idx <- min(max(round(click$x), 1), length(levels(ph$mrca_name)))
    x <- levels(ph$mrca_name)[x_idx]

    src <- if (input$no_isoforms) ph_filt else ph
    src %>%
      filter(mrca_name == x & gene %in% genes_of_interest_reactive()) %>%
      dplyr::select(Gene = gene, `UniProt ID` = qseqid, Phylostrata = mrca_name)
  }, server = TRUE)

  ##### rate_plot: odds ratio per phylostratum #####
  output$rate_plot <- renderPlot({
    ph_local <- if (input$no_isoforms) ph_filt else ph

    genes_of_interest <- genes_of_interest_reactive()
    disease_label <- disease_name()

    ph_local$is_of_interest <- ph_local$Gene %in% genes_of_interest
    ph_local$annot_of_interest <- ifelse(
      ph_local$is_of_interest,
      disease_label,
      "All Human Genes"
    )

    ph_by_ps <- ph_local %>%
      dplyr::distinct() %>%
      dplyr::mutate(
        mrca_name         = factor(mrca_name, levels = levels(ph$mrca_name), ordered = TRUE),
        annot_of_interest = relevel(as.factor(annot_of_interest), ref = "All Human Genes")
      ) %>%
      dplyr::count(mrca_name, annot_of_interest, .drop = FALSE)

    ratio_plot_data <- ph_by_ps %>%
      tidyr::pivot_wider(
        id_cols     = mrca_name,
        names_from  = annot_of_interest,
        values_from = n,
        values_fill = 0
      )

    if (!disease_label %in% colnames(ratio_plot_data)) {
      return(
        ggplot() +
          annotate(geom = "text", x = 0, y = 0,
                   label = "No disease genes available for rate comparison.") +
          theme_void()
      )
    }

    ratio_plot_data <- ratio_plot_data %>%
      dplyr::mutate(
        Disease_In_Stratum        = .data[[disease_label]],
        Background_In_Stratum     = .data[["All Human Genes"]],
        Disease_Not_In_Stratum    = sum(Disease_In_Stratum) - Disease_In_Stratum,
        Background_Not_In_Stratum = sum(Background_In_Stratum) - Background_In_Stratum
      ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        fisher_p = fisher.test(matrix(c(
          Disease_In_Stratum,    Disease_Not_In_Stratum,
          Background_In_Stratum, Background_Not_In_Stratum
        ), nrow = 2, byrow = TRUE))$p.value,
        Odds_Ratio = dplyr::case_when(
          Disease_Not_In_Stratum == 0 | Background_In_Stratum == 0 ~ NA_real_,
          TRUE ~ (Disease_In_Stratum * Background_Not_In_Stratum) /
            (Disease_Not_In_Stratum * Background_In_Stratum)
        ),
        Disease_Rate          = Disease_In_Stratum / sum(Disease_In_Stratum),
        Background_Rate       = Background_In_Stratum / sum(Background_In_Stratum),
        Normalized_Rate_Ratio = dplyr::case_when(
          Background_Rate == 0 ~ NA_real_,
          TRUE ~ Disease_Rate / Background_Rate
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        fisher_fdr   = p.adjust(fisher_p, method = "BH"),
        Significance = dplyr::case_when(
          !is.na(fisher_fdr) & fisher_fdr < 0.05 ~ "*",
          TRUE                                   ~ ""
        )
      ) %>%
      dplyr::filter(!is.na(Odds_Ratio))

    ggplot(ratio_plot_data, aes(x = mrca_name, y = Odds_Ratio)) +
      geom_hline(yintercept = 1, color = "grey80") +
      geom_point() +
      geom_line(aes(group = 1), linewidth = 1.25) +
      geom_text(
        data = ratio_plot_data %>% dplyr::filter(Odds_Ratio >= 1),
        aes(label = Significance),
        vjust = 1, nudge_y = 0.1, size = 5, color = "red"
      ) +
      labs(
        x = "Phylostrata",
        y = "Odds Ratio of Enrichment\n(Disease Genes vs Overall Human Genome)",
        title = paste0(
          "Odds ratio of emergence of <br><span style='color:#E41A1C;'><strong>",
          disease_label,
          "</strong></span> relative to overall human genes<br>"
        )
      ) +
      ggpubr::labs_pubr() +
      my_theme +
      theme(
        axis.text.x         = element_text(angle = 90, hjust = 1),
        plot.title.position = "plot",
        plot.title          = ggtext::element_markdown(size = 12, lineheight = 1, hjust = 0.5)
      )
  })

  ##### GO enrichment via topGO Fisher (ORA) #####
  GO_df_reactive <- eventReactive(input$run, {
    ph_local <- if (isolate(input$no_isoforms)) ph_filt else ph

    ph_local$is_of_interest <- ph_local$Gene %in% genes_of_interest_reactive()

    ph_of_interest <- ph_local %>% filter(is_of_interest == TRUE)

    if (nrow(ph_of_interest) < 1) {
      return(data.frame())
    }

    gene_list_by_strata <- lapply(unique(ph_of_interest$ps), function(x) {
      ph_of_interest$gene[ph_of_interest$ps == x]
    })
    names(gene_list_by_strata) <- unique(ph_of_interest$ps)

    allGenes_universe <- unique(unlist(annFUN.org.symbol))

    GO_df_by_strata <- withProgress(
      message = "Running GO enrichment per phylostratum...",
      value = 0, {
        out <- lapply(seq_along(gene_list_by_strata), function(i) {
          incProgress(1 / length(gene_list_by_strata),
                      detail = names(gene_list_by_strata)[i])

          gene_set_to_enrich <- unique(gene_list_by_strata[[i]])
          gene_set_to_enrich <- gene_set_to_enrich[!is.na(gene_set_to_enrich)]
          if (length(gene_set_to_enrich) < 2) return(NULL)

          geneList <- factor(as.integer(allGenes_universe %in% gene_set_to_enrich))
          names(geneList) <- allGenes_universe

          if (sum(geneList == 1) < 1) return(NULL)

          GOdata <- suppressMessages(new(
            "topGOdata",
            ontology = "BP",
            allGenes = geneList,
            annot    = annFUN.org,
            mapping  = "org.Hs.eg.db",
            ID       = "symbol"
          ))
          resultFisher <- suppressMessages(
            runTest(GOdata, algorithm = "classic", statistic = "fisher")
          )
          allres <- GenTable(
            GOdata, classic = resultFisher,
            numChar  = 1000,
            topNodes = length(score(resultFisher))
          )

          allres$pval <- suppressWarnings(as.numeric(gsub("<", "", allres$classic)))
          allres$pval[is.na(allres$pval)] <- .Machine$double.eps

          allres$Significant <- as.numeric(allres$Significant)
          allres$Expected    <- as.numeric(allres$Expected)
          allres$Fold_Enrichment <- ifelse(
            allres$Expected > 0,
            allres$Significant / allres$Expected,
            NA_real_
          )
          allres$padj    <- p.adjust(allres$pval, method = "BH")
          allres$pathway <- allres$GO.ID
          allres
        })
        names(out) <- names(gene_list_by_strata)
        out
      })

    GO_df_by_strata <- GO_df_by_strata[!sapply(GO_df_by_strata, is.null)]
    if (length(GO_df_by_strata) == 0) return(data.frame())

    GO_res_df <- dplyr::bind_rows(
      lapply(GO_df_by_strata, as.data.frame), .id = "PS"
    ) %>%
      mutate(PS = factor(ph_local$mrca_name[match(PS, ph_local$ps)],
                         levels = levels(ph$mrca_name))) %>%
      filter(!is.na(pval))

    GO_res_df
  })

  GO_filtered_table <- reactive({
    req(GO_df_reactive())
    GO_df_reactive() %>%
      dplyr::filter(pval >= input$pval_range[1],
                    pval <= input$pval_range[2])
  })

  output$download_go_table <- downloadHandler(
    filename = function() paste0("GO_Enrichment_", Sys.Date(), ".csv"),
    content  = function(file) {
      go_df <- GO_filtered_table()
      if (is.data.frame(go_df) && nrow(go_df) > 0) {
        write.csv(go_df, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Note = "No GO enrichment results available"),
                  file, row.names = FALSE)
      }
    }
  )

  ##### Static GO plot #####
  output$GOplot_static <- renderPlot({
    GO_res_df <- GO_df_reactive()
    req(nrow(GO_res_df) > 0)

    cappval <- input$cappval
    if (cappval < 100) {
      GO_res_df$pval[-log10(GO_res_df$pval) > cappval] <- 10^(-1 * cappval)
    }

    min_lim <- input$pval_range[1]
    min_lim <- ifelse(min_lim < 10^(-1 * cappval), 10^(-1 * cappval), min_lim)
    max_lim <- input$pval_range[2]
    max_lim <- ifelse(max_lim == min_lim, 1, max_lim)

    GO_res_df <- GO_res_df %>%
      dplyr::filter(pval >= min_lim & pval <= max_lim)
    req(nrow(GO_res_df) > 0)

    GO_res_df_top <- GO_res_df %>%
      group_by(PS) %>%
      arrange(desc(pval)) %>%
      slice_min(n = input$labeltop, order_by = pval, with_ties = FALSE) %>%
      ungroup()

    palette <- scales::hue_pal()(length(unique(ph$mrca_name)))
    names(palette) <- levels(ph$mrca_name)

    ggplot(GO_res_df, aes(x = PS, y = as.numeric(pval), label = Term, color = PS)) +
      geom_point(show.legend = FALSE, position = position_jitter(0.1), size = 4) +
      ggrepel::geom_label_repel(
        data = GO_res_df_top,
        aes(label = Term, y = as.numeric(pval)),
        show.legend = FALSE, size = 3,
        hjust = 0.5, vjust = 1, max.overlaps = 30
      ) +
      scale_color_manual(values = palette) +
      scale_x_discrete(limits = levels(GO_res_df$PS)) +
      scale_y_continuous(transform = scales::new_transform(
        name      = "-log10",
        transform = function(x) -log10(x),
        inverse   = function(x) 10^(-1 * x)
      )) +
      my_theme +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(
        x = "Phylostrata (PS)",
        y = "P-value (Fisher exact, topGO)",
        title = paste0("Top enriched GO terms by phylostrata for ", disease_name())
      ) +
      ggpubr::labs_pubr()
  })

  ##### GO plot click #####
  output$selected_go <- renderPlot({
    click <- input$go_plot_click
    if (is.null(click) || is.null(click$x) || is.null(click$y)) return(NULL)

    x_idx <- min(max(round(click$x), 1), length(levels(ph$mrca_name)))
    x <- levels(ph$mrca_name)[x_idx]

    palette <- scales::hue_pal()(length(unique(ph$mrca_name)))
    names(palette) <- levels(ph$mrca_name)

    GO_res_df <- GO_df_reactive() %>%
      filter(PS == x) %>%
      dplyr::select(Term, `P-value` = pval,
                    `Fold Enrichment` = Fold_Enrichment,
                    `Adjusted P-value` = padj)
    req(nrow(GO_res_df) > 0)

    ggplot(GO_res_df, aes(x = `Fold Enrichment`, y = `Adjusted P-value`)) +
      geom_point(color = palette[x]) +
      ggrepel::geom_label_repel(
        data = GO_res_df %>% slice_min(`Adjusted P-value`, n = input$labeltop,
                                       with_ties = FALSE),
        aes(label = Term), max.overlaps = 30
      ) +
      scale_y_continuous(transform = scales::new_transform(
        name      = "-log10",
        transform = function(x) -log10(x),
        inverse   = function(x) 10^(-1 * x)
      )) +
      my_theme +
      labs(
        x = "Fold Enrichment (Significant / Expected)",
        y = "-log10 Adjusted P-value",
        title = paste0("Top enriched GO terms for ", disease_name(), " \nin PS ", x)
      ) +
      ggpubr::labs_pubr()
  })

  ##### STRING plot (currently not exposed in UI; kept for future use) #####
  output$string_plot <- renderPlot({
    ph_local <- if (input$no_isoforms) ph_filt else ph
    string_db <- STRINGdb::STRINGdb$new(species = 9606, version = "12",
                                        score_threshold = 200, input_directory = ".")
    ph_local$is_of_interest <- ph_local$Gene %in% genes_of_interest_reactive()
    ph_local$annot_of_interest <- ifelse(
      ph_local$is_of_interest, disease_name(), "Other human genes"
    )

    withProgress(message = "Pulling data from STRINGdb...", value = 0, {
      ph_mapped <- string_db$map(
        ph_local %>% filter(is_of_interest & mrca_name %in% input$ps) %>% as.data.frame(),
        "gene", removeUnmappedRows = TRUE
      )
      incProgress(amount = 0.5, message = "Generating plot...")
      string_db$plot_network(ph_mapped$STRING_id, required_score = 800)
      incProgress(amount = 0.5, message = "Done!")
    })
  })
}

##### RUN #####
shinyApp(ui = ui, server = server)
