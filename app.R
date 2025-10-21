# ==============================================================================
#
# Shiny App for Displaying Meta-analysis Data from glial RMA JCB Publication
#
# Author: Jeffrey Y Lee
# Date: Oct 2025
#
# Description:
# This Shiny application provides an interactive platform to explore the datasets
# associated with the glial localised RNA JCB publication. It includes an
# introduction to the study, a searchable table of the analysed dataset.
#
# --- DIRECTORY STRUCTURE (for shinylive) ---
#
# /app_root/
# |-- app.R         (this file)
# |-- data/
# |-- www/
#
# ==============================================================================


# 1. LOAD LIBRARIES
# ------------------------------------------------------------------------------
library(shiny)
library(DT)
library(dplyr)


# 2. LOAD DATA
# ------------------------------------------------------------------------------
# Load the main dataset
main_data <- tryCatch({
  readRDS("data/1740_main_table.RDS")
}, error = function(e) {
  # Create a placeholder dataframe if the file is missing
  data.frame(
    gene_name = c("pros", "GeneB", "GeneC"), # Ensure this column exists
    chromosome = c("chr1", "chr2", "chrX"),
    log2FC = c(1.5, -0.5, 2.1),
    p_value = c(0.001, 0.2, 0.0005),
    peak_sequence = c("AUCG", "GCUA", "UUUU"),
    notes = c("Note 1", "Note 2", "Note 3"),
    stringsAsFactors = FALSE
  )
})


do_data <- tryCatch({
  readRDS("data/disease-ontology.RDS")
}, error = function(e) {
  # Create a placeholder dataframe if the file is missing
  data.frame(
    gene_name = c("pros", "GeneB", "GeneC"), # Ensure this column exists
    chromosome = c("chr1", "chr2", "chrX"),
    log2FC = c(1.5, -0.5, 2.1),
    p_value = c(0.001, 0.2, 0.0005),
    peak_sequence = c("AUCG", "GCUA", "UUUU"),
    notes = c("Note 1", "Note 2", "Note 3"),
    stringsAsFactors = FALSE
  )
})


reactome_data <- tryCatch({
  readRDS("data/reactome.RDS")
}, error = function(e) {
  # Create a placeholder dataframe if the file is missing
  data.frame(
    gene_name = c("pros", "GeneB", "GeneC"), # Ensure this column exists
    chromosome = c("chr1", "chr2", "chrX"),
    log2FC = c(1.5, -0.5, 2.1),
    p_value = c(0.001, 0.2, 0.0005),
    peak_sequence = c("AUCG", "GCUA", "UUUU"),
    notes = c("Note 1", "Note 2", "Note 3"),
    stringsAsFactors = FALSE
  )
})


# Complementary non-localised disease ontology dataset
do_nonloc_data <- tryCatch({
  readRDS("data/disease-ontology-nonlocalised.RDS")
}, error = function(e) {
  data.frame(
    gene_name = c("pros", "GeneB", "GeneC"),
    chromosome = c("chr1", "chr2", "chrX"),
    log2FC = c(0.2, -0.1, 0.5),
    p_value = c(0.05, 0.3, 0.01),
    peak_sequence = c("AUGC", "GCAA", "UUAA"),
    notes = c("Nonloc Note 1", "Nonloc Note 2", "Nonloc Note 3"),
    stringsAsFactors = FALSE
  )
})


# Complementary non-localised reactome dataset
reactome_nonloc_data <- tryCatch({
  readRDS("data/reactome-nonlocalised.RDS")
}, error = function(e) {
  data.frame(
    gene_name = c("pros", "GeneB", "GeneC"),
    chromosome = c("chr1", "chr2", "chrX"),
    log2FC = c(0.3, -0.2, 0.6),
    p_value = c(0.02, 0.25, 0.005),
    peak_sequence = c("AUAA", "GCAA", "UUCC"),
    notes = c("Nonloc R Note 1", "Nonloc R Note 2", "Nonloc R Note 3"),
    stringsAsFactors = FALSE
  )
})


# 3. DEFINE UI (USER INTERFACE)
# ------------------------------------------------------------------------------
ui <- navbarPage(
  title = "Glial Localised RNA meta-analysis",
  collapsible = TRUE,

  # --- Main Page Tab ---
  tabPanel(
    "Data",

    # --- Custom CSS Styling ---
    tags$head(
      tags$style(
        HTML(
          "
        /* Target the full data table and its components */
        #main_data_table .table td,
        #main_data_table .table th {
          font-size: 90%; /* Adjust data cell and header font size */
        }

        #main_data_table .form-control {
          font-size: 90%; /* Adjust column filter input font size */
        }

        /* Target the wrapper for global search, pagination, etc. */
        .dataTables_wrapper {
          font-size: 90%; /* Adjust font size for table controls */
        }
      "
        )
      )
    ),

    sidebarLayout(
      # --- Sidebar Panel ---
      sidebarPanel(
        width = 2,

        # Publication Information Panel
        wellPanel(
          tags$p(
            tags$strong(
              "Murine glial protrusion transcripts predict localized Drosophila glial mRNAs involved in plasticity"
            )
          ),
          tags$p("Lee JY & Gala DS et al. (2024)"),
          tags$a(
            href = "https://rupress.org/jcb/article/223/10/e202306152/276869",
            target = "_blank",
            tags$img(src = "jcb_logo.png", height = "30px", style = "margin-right: 10px;"),
            tags$br(),
            "Link to the Article"
          ),
          tags$br(),
          tags$br(),
          tags$p(
            tags$a(
              href = "https://github.com/JYLeeSci/Gala2023_glia-localised-RNA",
              icon("github"),
              "GitHub source data",
              target = "_blank"
            )
          )
        ),

        # Contact Details Panel
        wellPanel(
          h4("Contact"),
          tags$p(
            "Jeffrey Y Lee",
            tags$a(href = "mailto:jeff.lee@glasgow.ac.uk", icon("envelope"), "Email", target = "_blank")
          ),
          tags$p(
            "Ilan Davis",
            tags$a(href = "mailto:ilan.davis@glasgow.ac.uk", icon("envelope"), "Email", target = "_blank")
          ),
          tags$br(),
          tags$p(tags$a(href = "https://ilandavis.com/", icon("globe"), "Lab website", target = "_blank")),
          tags$p(tags$a(href = "https://x.com/jefflee1103", icon("x-twitter"), "Twitter", target = "_blank"))
        )
      ),

      # --- Main Panel ---
      mainPanel(
        width = 9,
        tabsetPanel(
          id = "main_tabs",

          # Tab 1: Introduction
          tabPanel(
            "Introduction",
            div(
              style = "max-width: 85%;",
              h3("Glial Localised RNA Meta-analysis Data Explorer"),
              p(
                HTML(
                  "This web application provides a user-friendly interface to explore the data from our recent study on the glial protrusion localised RNAs. Here, you can browse the complete dataset, search for specific genes of interest, and identifiy their disease and reactome term association. The goal of this resource is to enhance the reproducibility and accessibility of our findings."
                )
              ),
              hr(),
              h4("Methodology"),
              p(
                HTML(
                  "We combined all available mammalian studies that profile RNAs in the far-reaching projections of glial cells and standardised the data so they could be compared. This revealed a shared set of more than 5,000 glial RNAs enriched in cell peripheries, many of which also appear in neuronal projections and in other protrusive cells, pointing to conserved transport and localisation mechanisms. To test conservation across species, we mapped high-confidence mammalian to Drosophila gene matches and cross-checked them against the Fly Cell Atlas for glial types that extend long processes in the peripheral nervous system. This predicted a curated list of 1,740 RNAs likely present in the projections of Drosophila peripheral glia."
                )
              ),
              tags$br(),
              fluidRow(
                column(
                  6,
                  tags$figure(
                    tags$img(
                      src = "JCB_figure1.jpg",
                      width = "100%",
                      style = "border: 1px solid #ddd; border-radius: 4px; padding: 5px;"
                    ),
                    tags$figcaption(
                      "Conserved transcripts localized to glial protrusions. (A) Graphical representation of the techniques used to separate the protrusion-localized transcripts in each of the studies analyzed here. (B) Bar graph showing the datasets included in this study. The high confidence expression cut-off was set at TPM >10, and transcripts detected in more than three independent datasets are shown in orange color. (C) Identification of high-confidence D. melanogaster orthologs of 5,028 mammalian genes that were detected in at least 7 datasets. DRSC Integrative Ortholog Prediction Tool (DIOPT) score of 8 was used as a cut-off. (D) Comparison of transcripts localized in the glial protrusion and neurites (von Kügelgen and Chekulaeva, 2020). (E) Upset plot showing a comparison of cellular protrusion localized transcripts between glia, neuron, and breast cancer cell line (MDA-MB231) (Mardakheh et al., 2015). (F) t-SNE plot of combined glial single-nuclei RNA-seq from the Fly Cell Atlas (Li et al., 2022). Perineurial, subperineurial, and ensheathing glial cell clusters are depicted in orange-red colors. (G) A schematic representing the distribution and location of glial cells in the Drosophila third instar larva. Six subtypes of glia exist in the larvae: perineurial and subperineurial glia, cortex glia, astrocyte-like and ensheathing glial cells and, finally, the peripheral nervous system (PNS) specific wrapping glia (Yildirim et al., 2019). Ensheathing glia of the CNS are continuous with the wrapping glia of the CNS. The perineurial and subperineurial glia have long and extensive projections which reach all the way to the neuromuscular junction (NMJ)."
                    )
                  )
                )
              ),
              tags$br(),
              hr(),
              tags$br(),
              tags$br(),
              tags$br()
            )
          ),

          # Tab 2: Explore Full Data Table
          tabPanel(
            "1,740 glial localised RNAs",
            h3("Meta-analysis output of 1,740 glial localised RNAs"),
            p("The table below contains the full set of 1,740 glial localised RNAs in Drosophila gene space (mouse homologs are also provided)."),
            hr(),
            div(
              style = "background-color: #f5f5f5; padding: 15px; border-radius: 5px; margin-bottom: 20px;",
              h4("Column Definitions"),
              p(HTML("<strong>dmel_gene_name:</strong> Drosophila gene name")),
              p(HTML("<strong>mmus_gene_symbols:</strong> Mouse gene symbols of homologs")),
              p(HTML("<strong>rna_in_protrusion:</strong> Number of datasets (out of 12) supporting RNA localisation to glial protrusions")),
              p(HTML("<strong>mean_tpm_across_dataset:</strong> Averge TPM across all glial RNAseq datasets (abundance)")),
              p(HTML("<strong>go_molecular_function:</strong> Gene ontology (Molecular Function) terms associated with gene")),
              p(HTML("<strong>go_biological_process:</strong> Gene ontology (Biological Process) terms associated with gene")),
              p(HTML("<strong>go_cellular_component:</strong> Gene ontology (Cellular Component) terms associated with gene"))
            ),
            hr(),
            DT::dataTableOutput("main_data_table"),
            downloadButton("download_main_data", "Download Full Table (.csv)", class = "btn-primary", style = "margin-top: 15px;"),
            tags$br(),
            tags$br(),
            tags$br()
          ),

          # Tab 3: Disease ontology Enrichment
          tabPanel(
            "Disease Ontology",
            h3("Human disease ontolgy enrichment analysis of 1,740 glial localised RNAs"),
            
            # Explanatory text for the plots
            div(style = "margin-bottom: 10px;",
                tags$p("Interactive volcano plots of enriched disease ontology terms. Left: predicted localised terms; Right: predicted NOT localised terms."),
                tags$p("Click on the term to see the associated genes.")
            ),
            hr(),
       # Two interactive plotly plots side-by-side
       fluidRow(
    column(6,
      plotly::plotlyOutput("do_plot", height = "590px")
    ),
    column(6,
      plotly::plotlyOutput("do_nonloc_plot", height = "590px")
    )
       ),
     # Single shared click-driven table below both plots
     fluidRow(
   column(12, DT::dataTableOutput("do_click_table"))
     ),
     tags$br(),
     tags$br(),
     # prominent separator between reactive table and following static content
     tags$hr(style = "border-top: 3px solid #444; margin-top: 10px; margin-bottom: 18px;"),
        # Move the descriptive header and table intro below the plots
        tags$br(),
            h3("Full table: Human disease ontolgy enrichment analysis of 1,740 glial localised RNAs"),
            p("The table below contains disease ontology enrichment analysis for each human disease term and the list of associated genes."),
            hr(),
            div(
              style = "background-color: #f5f5f5; padding: 15px; border-radius: 5px; margin-bottom: 20px;",
              h4("Column Definitions"),
              p(HTML("<strong>disease:</strong> Human disease term")),
              p(HTML("<strong>log2foldchange:</strong> log2 fold change of term enrichment in localised genes vs genome background")),
              p(HTML("<strong>adj_pvalue:</strong> Bonferroni adjusted p-value of term enrichment")),
              p(HTML("<strong>count_in_genome:</strong> Total number of genes associated with term in genome background")),
              p(HTML("<strong>count_in_localised_genes:</strong> Number of glial localised genes associated with term")),
              p(HTML("<strong>dmel_gene_list:</strong> List of Drosophila genes associated with the disease term")),
              p(HTML("<strong>mmus_gene_list:</strong> List of mouse genes - direct conversion of fly-to-mouse genes with high-homolog score"))
            ),
            hr(),
            DT::dataTableOutput("do_table"),
            downloadButton("download_do_data", "Download Full Table (.csv)", class = "btn-primary", style = "margin-top: 15px;"),
            tags$br(),
            tags$br(),
            tags$br()
          ),

          # Tab 4: Reactome Enrichment
          tabPanel(
            "Reactome",
            h3("Reactome enrichment analysis of 1,740 glial localised RNAs"),
            div(style = "margin-bottom: 10px;",
                tags$p("Interactive volcano plot. Click on the term to see the associated genes."),
            ),
            hr(),
            div(
              style = "margin-bottom: 10px;",
              plotly::plotlyOutput("reactome_plotly", height = "708px"),
              tags$br(),
              DT::dataTableOutput("reactome_click_table")
            ),
            tags$br(),
            tags$br(),
            tags$hr(style = "border-top: 3px solid #444; margin-top: 10px; margin-bottom: 18px;"),

            tags$br(),
            h3("Full table: Reactome enrichment analysis of 1,740 glial localised RNAs"),
            p("The table below contains Reactome enrichment analysis for each term and the list of associated genes."),
            hr(),
            div(
              style = "background-color: #f5f5f5; padding: 15px; border-radius: 5px; margin-bottom: 20px;",
              h4("Column Definitions"),
              p(HTML("<strong>name:</strong> Reactome term")),
              p(HTML("<strong>log2foldchange:</strong> log2 fold change of term enrichment in localised genes vs genome background")),
              p(HTML("<strong>adj_pvalue:</strong> Bonferroni adjusted p-value of term enrichment")),
              p(HTML("<strong>count_in_genome:</strong> Total number of genes associated with term in genome background")),
              p(HTML("<strong>count_in_localised_genes:</strong> Number of glial localised genes associated with term")),
              p(HTML("<strong>dmel_gene_list:</strong> List of Drosophila genes associated with the Reactome term")),
              p(HTML("<strong>mmus_gene_list:</strong> List of mouse genes - direct conversion of fly-to-mouse genes with high-homolog score"))
            ),
            hr(),
            DT::dataTableOutput("reactome_table"),
            downloadButton("download_reactome_data", "Download Full Table (.csv)", class = "btn-primary", style = "margin-top: 15px;"),
            tags$br(),
            tags$br(),
            tags$br()
          )
        )
      )
    )
  )
)


# 4. DEFINE SERVER LOGIC
# ------------------------------------------------------------------------------
server <- function(input, output, session) {

  # --- Server logic for Tab 2: Full Data Table ---
  output$main_data_table <- DT::renderDataTable({
    DT::datatable(
      main_data,
      options = list(pageLength = 15, autoWidth = TRUE, scrollX = TRUE),
      filter = 'top',
      rownames = FALSE,
      class = 'cell-border stripe'
    )
  })

  # --- Server logic for Tab 3: Disease ontology Table ---
  output$do_table <- DT::renderDataTable({
    DT::datatable(
      do_data,
      options = list(pageLength = 15, autoWidth = TRUE, scrollX = TRUE),
      filter = 'top',
      rownames = FALSE,
      class = 'cell-border stripe'
    )
  })

  # Plotly volcano-like plot for disease ontology (localised)
  output$do_plot <- plotly::renderPlotly({
    df <- do_data %>%
      filter(!is.infinite(log2foldchange)) %>%
      mutate(yval = -log10(adj_pvalue)) %>%
      dplyr::rename("count" = count_in_localised_genes)

    # compute shared ranges and color scale across both datasets
    df_nonloc_tmp <- do_nonloc_data %>% filter(!is.infinite(log2foldchange)) %>% mutate(yval = -log10(adj_pvalue))
    x_min <- min(c(df$log2foldchange, df_nonloc_tmp$log2foldchange), na.rm = TRUE)
    x_max <- max(c(df$log2foldchange, df_nonloc_tmp$log2foldchange), na.rm = TRUE)
    y_min <- min(c(df$yval, df_nonloc_tmp$yval), na.rm = TRUE)
    y_max <- max(c(df$yval, df_nonloc_tmp$yval), na.rm = TRUE)
    count_min <- min(c(df$count, df_nonloc_tmp$count), na.rm = TRUE)
    count_max <- max(c(df$count, df_nonloc_tmp$count), na.rm = TRUE)

    df <- df %>% mutate(put_label = if_else(adj_pvalue < 0.001 & log2foldchange >= 1.8, disease, NA_character_))

    # size mapping: scale counts into reasonable pixel sizes
    size_min <- 6
    size_max <- 18
    if (!all(is.na(df$count))) {
      sizes <- scales::rescale(df$count, to = c(size_min, size_max), from = c(count_min, count_max))
    } else {
      sizes <- rep((size_min + size_max)/2, nrow(df))
    }

    # base scatter with size mapped to count; hide automatic legend for trace
    p <- plotly::plot_ly(
      data = df,
      x = ~log2foldchange,
      y = ~yval,
      customdata = ~disease,
      source = "do_plot",
      type = 'scatter',
      mode = 'markers',
      marker = list(
        size = sizes, 
        opacity = 0.85, 
        line = list(width = 0), 
        cmin = count_min, 
        cmax = count_max, 
        colorbar = list(
          title = list(text = ""),
          orientation = 'h', 
          x = 0.5, 
          xanchor = 'center', 
          y = -0.25,
          yanchor = 'top',
          len = 0.7
        )
      ),
      color = ~count,
      colors = viridis::inferno(20),
      showlegend = FALSE
    )

    # add text labels separately (only where put_label not NA)
    labels_df <- df %>% filter(!is.na(put_label))
    if (nrow(labels_df) > 0) {
      p <- p %>% plotly::add_text(data = labels_df, x = ~log2foldchange, y = ~yval,
                                   text = ~put_label, textposition = 'right', showlegend = FALSE,
                                   textfont = list(color = 'gray40', size = 13))
    }

    # add dashed reference lines
    p <- p %>%
      plotly::layout(
        title = list(text = 'Predicted localised', x = 0.02, xanchor = 'left'),
        shapes = list(
          list(type = 'line', x0 = 0, x1 = 0, xref = 'x', y0 = 0, y1 = y_max, yref = 'y',
               line = list(dash = 'dash', color = 'gray80', width = 0.8)),
          list(type = 'line', x0 = 0.8, x1 = 3.5, xref = 'x',
               y0 = -log10(0.01), y1 = -log10(0.01), yref = 'y', line = list(dash = 'dash', color = 'gray80', width = 0.8))
        ),
        xaxis = list(title = 'log2FoldChange', range = c(0.8, 3.5)),
        yaxis = list(title = 'Adjusted p-value (-log10)', range = c(0, 49)),
        showlegend = FALSE,
        margin = list(l = 60, r = 20, b = 80, t = 40)
      )

    p
  })

  # Plotly volcano-like plot for disease ontology (non-localised)
  output$do_nonloc_plot <- plotly::renderPlotly({
    df <- do_nonloc_data %>%
      filter(!is.infinite(log2foldchange)) %>%
      mutate(yval = -log10(adj_pvalue)) %>%
      dplyr::rename("count" = count_in_localised_genes)

    # compute same shared ranges as used for the localised plot
    df_local_tmp <- do_data %>% filter(!is.infinite(log2foldchange)) %>% mutate(yval = -log10(adj_pvalue)) %>% dplyr::rename("count" = count_in_localised_genes)
    x_min <- min(c(df$log2foldchange, df_local_tmp$log2foldchange), na.rm = TRUE)
    x_max <- max(c(df$log2foldchange, df_local_tmp$log2foldchange), na.rm = TRUE)
    y_min <- min(c(df$yval, df_local_tmp$yval), na.rm = TRUE)
    y_max <- max(c(df$yval, df_local_tmp$yval), na.rm = TRUE)
    count_min <- min(c(df$count, df_local_tmp$count), na.rm = TRUE)
    count_max <- max(c(df$count, df_local_tmp$count), na.rm = TRUE)

    # label for non-localised: anything above y > 2
    df <- df %>% mutate(put_label = if_else(yval > 2, disease, NA_character_))

    # size mapping for nonlocalised
    size_min <- 6
    size_max <- 18
    if (!all(is.na(df$count))) {
      sizes <- scales::rescale(df$count, to = c(size_min, size_max), from = c(count_min, count_max))
    } else {
      sizes <- rep((size_min + size_max)/2, nrow(df))
    }

    p <- plotly::plot_ly(
      data = df,
      x = ~log2foldchange,
      y = ~yval,
      customdata = ~disease,
      source = "do_nonloc_plot",
      type = 'scatter',
      mode = 'markers',
      marker = list(
        size = sizes,
        opacity = 0.85,
        line = list(width = 0),
        cmin = count_min,
        cmax = count_max,
        colorbar = list(
          title = list(text = ""),
          orientation = 'h',
          x = 0.5,
          xanchor = 'center',
          y = -0.18
        )
      ),
      color = ~count,
      colors = viridis::inferno(20),
      showlegend = FALSE
    )

    labels_df <- df %>% filter(!is.na(put_label))
    if (nrow(labels_df) > 0) {
      p <- p %>% plotly::add_text(data = labels_df, x = ~log2foldchange, y = ~yval,
                                   text = ~put_label, textposition = 'right', showlegend = FALSE,
                                   textfont = list(color = 'gray40', size = 11))
    }

    p <- p %>%
      plotly::layout(
        title = list(text = 'Predicted NOT localised', x = 0.02, xanchor = 'left'),
        shapes = list(
          list(type = 'line', x0 = 0, x1 = 0, xref = 'x', y0 = 0, y1 = y_max, yref = 'y',
               line = list(dash = 'dash', color = 'gray80', width = 0.8)),
          list(type = 'line', x0 = 0.8, x1 = 3.5, xref = 'x',
               y0 = -log10(0.01), y1 = -log10(0.01), yref = 'y', line = list(dash = 'dash', color = 'gray80', width = 0.8))
        ),
        xaxis = list(title = 'log2FoldChange', range = c(0.8, 3.5)),
        yaxis = list(title = 'Adjusted p-value (-log10)', range = c(0, 49)),
        showlegend = FALSE,
        margin = list(l = 60, r = 20, b = 80, t = 40)
      )

    p
  })

  # Use reactiveValues and observers to capture clicks from either plot
  selected <- reactiveValues(data = NULL, source = NULL)

  observeEvent(plotly::event_data("plotly_click", source = "do_plot"), {
    ed <- plotly::event_data("plotly_click", source = "do_plot")
    if (!is.null(ed) && nrow(ed) > 0) {
      clicked_term <- ed$customdata[1]
      if (!is.null(clicked_term) && clicked_term != "") {
        selected$data <- do_data %>% filter(disease == clicked_term)
        selected$source <- "local"
      }
    }
  })

  observeEvent(plotly::event_data("plotly_click", source = "do_nonloc_plot"), {
    ed <- plotly::event_data("plotly_click", source = "do_nonloc_plot")
    if (!is.null(ed) && nrow(ed) > 0) {
      clicked_term <- ed$customdata[1]
      if (!is.null(clicked_term) && clicked_term != "") {
        selected$data <- do_nonloc_data %>% filter(disease == clicked_term)
        selected$source <- "nonloc"
      }
    }
  })

  output$do_click_table <- DT::renderDataTable({
    if (is.null(selected$data)) {
      DT::datatable(data.frame(message = "Click a point on either volcano to show the disease row"), options = list(dom = 't'), rownames = FALSE)
    } else {
      # Only show the requested columns if they exist
      wanted <- c("disease", "dmel_gene_list", "mmus_gene_list")
      present <- intersect(wanted, colnames(selected$data))
      if (length(present) == 0) {
        DT::datatable(data.frame(message = "Selected row has none of the requested columns (disease, dmel_gene_list, mmus_gene_list)"), options = list(dom = 't'), rownames = FALSE)
      } else {
        DT::datatable(selected$data %>% dplyr::select(all_of(present)), options = list(pageLength = 5, dom = 't'), rownames = FALSE, class = 'cell-border stripe')
      }
    }
  })

  # --- Server logic for Tab 4: Reactome Table ---
  output$reactome_table <- DT::renderDataTable({
    DT::datatable(
      reactome_data,
      options = list(pageLength = 15, autoWidth = TRUE, scrollX = TRUE),
      filter = 'top',
      rownames = FALSE,
      class = 'cell-border stripe'
    )
  })

  # Plotly volcano-like plot for Reactome enrichment
  output$reactome_plotly <- plotly::renderPlotly({
    # defensive: ensure expected columns exist (name, log2foldchange, adj_pvalue, count_in_gl)
    df <- reactome_data %>% dplyr::filter(!is.infinite(log2foldchange)) %>% dplyr::mutate(yval = -log10(adj_pvalue))

    if (!"count_in_gl" %in% colnames(df)) {
      # try common alternative names
      alt <- intersect(c("count_in_localised_genes", "count_in_localised_genes.x", "count"), colnames(df))
      if (length(alt) > 0) df$count_in_gl <- df[[alt[1]]] else df$count_in_gl <- 1
    }

    # label logic copied from ggplot pipeline
    df <- df %>%
      dplyr::mutate(put_label = if_else(yval > 15 | log2foldchange > 2.5, name, NA_character_)) %>%
      dplyr::mutate(put_label = if_else(yval < 5, NA_character_, put_label))

  # size mapping (increased ranges for larger symbols)
  # Previously 6-18; increase to 8-36 for greater visual impact
  size_min <- 5
  size_max <- 10
    count_min <- min(df$count_in_gl, na.rm = TRUE)
    count_max <- max(df$count_in_gl, na.rm = TRUE)
    if (!all(is.na(df$count_in_gl))) {
      sizes <- scales::rescale(df$count_in_gl, to = c(size_min, size_max), from = c(count_min, count_max))
    } else {
      sizes <- rep((size_min + size_max)/2, nrow(df))
    }

    y_max <- max(df$yval, na.rm = TRUE)

    # Use sizemode = 'diameter' and a sizeref to help plotly map data sizes to pixels
    # sizeref formula: 2.0 * max(size array) / (desired_max_marker_size_in_pixels^2)
    # We'll set desired_max_marker_size_in_pixels to size_max for a reasonable mapping
    desired_max_px <- size_max
    # avoid division by zero
    max_size_val <- ifelse(length(sizes) > 0, max(sizes, na.rm = TRUE), size_max)
    sizeref_val <- ifelse(max_size_val > 0, 2.0 * max_size_val / (desired_max_px^2), 1)

    p <- plotly::plot_ly(
      data = df,
      x = ~log2foldchange,
      y = ~yval,
      customdata = ~name,
      source = "reactome_plot",
      type = 'scatter',
      mode = 'markers',
      marker = list(size = sizes, opacity = 0.8, line = list(width = 0), sizemode = 'diameter', sizeref = sizeref_val),
      color = I('steelblue'),
      showlegend = FALSE
    )

    labels_df <- df %>% filter(!is.na(put_label))
    if (nrow(labels_df) > 0) {
      p <- p %>% plotly::add_text(data = labels_df, x = ~log2foldchange, y = ~yval,
                                   text = ~put_label, textposition = 'right', showlegend = FALSE,
                                   textfont = list(color = 'gray40', size = 11))
    }

    p <- p %>% plotly::layout(
      title = list(text = 'Reactome enrichment', x = 0.02, xanchor = 'left'),
      shapes = list(
        list(type = 'line', x0 = 1, x1 = 1, xref = 'x', y0 = 0, y1 = y_max, yref = 'y', line = list(dash = 'dash', color = 'gray80', width = 0.8)),
        list(type = 'line', x0 = 0.8, x1 = 3.5, xref = 'x', y0 = -log10(0.01), y1 = -log10(0.01), yref = 'y', line = list(dash = 'dash', color = 'gray80', width = 0.8))
      ),
      xaxis = list(title = list(text = 'log2FoldChange', font = list(size = 13)), range = c(1, 3.5)),
      yaxis = list(title = list(text = 'Adjusted p-value (-log10)', font = list(size = 13)), range = c(0, 70)),
      font = list(size = 13),
      margin = list(l = 60, r = 20, b = 80, t = 40)
    )

    # annotate p=0.01 label (moved to fit new x range)
  p <- p %>% plotly::add_annotations(x = 1.2, y = -log10(0.05) + 2, text = 'p=0.01', showarrow = FALSE, opacity = 0.5, font = list(size = 14))

    p
  })

  # Capture clicks on Reactome plot and show selected row
  reactome_selected <- reactiveValues(data = NULL)
  observeEvent(plotly::event_data("plotly_click", source = "reactome_plot"), {
    ed <- plotly::event_data("plotly_click", source = "reactome_plot")
    if (!is.null(ed) && nrow(ed) > 0) {
      clicked_name <- ed$customdata[1]
      if (!is.null(clicked_name) && clicked_name != "") {
        reactome_selected$data <- reactome_data %>% filter(name == clicked_name)
      }
    }
  })

  output$reactome_click_table <- DT::renderDataTable({
    if (is.null(reactome_selected$data)) {
      DT::datatable(data.frame(message = "Click a point on the Reactome plot to show the term row"), options = list(dom = 't'), rownames = FALSE)
    } else {
      wanted <- c("name", "dmel_gene_list", "mmus_gene_list")
      present <- intersect(wanted, colnames(reactome_selected$data))
      if (length(present) == 0) {
        DT::datatable(data.frame(message = "Selected row has none of the requested columns (name, dmel_gene_list, mmus_gene_list)"), options = list(dom = 't'), rownames = FALSE)
      } else {
        DT::datatable(reactome_selected$data %>% dplyr::select(all_of(present)), options = list(pageLength = 5, dom = 't'), rownames = FALSE, class = 'cell-border stripe')
      }
    }
  })

  # Download handler for the full dataset
  output$download_main_data <- downloadHandler(
    filename = function() {
      paste("Main_dataset_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(main_data, file, row.names = FALSE)
    }
  )

  # Download handler for the disease ontology dataset
  output$download_do_data <- downloadHandler(
    filename = function() {
      paste("Disease_ontology_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(do_data, file, row.names = FALSE)
    }
  )

  # Download handler for the reactome dataset
  output$download_reactome_data <- downloadHandler(
    filename = function() {
      paste("Reactome_enrichment_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(reactome_data, file, row.names = FALSE)
    }
  )

  # No additional server-side handlers required
}


# 5. RUN THE APPLICATION
# ------------------------------------------------------------------------------
shinyApp(ui = ui, server = server)
