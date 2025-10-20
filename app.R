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
        width = 3,

        # Publication Information Panel
        wellPanel(
          h4("Publication"),
          tags$p("This web resource accompanies the publication:"),
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
            h3("Human disease ontolgy enrichment analysis of of 1,740 glial localised RNAs"),
            p("The table below contains disease ontology enrichment analysis for each human disease term and the list of associated genes."),
            hr(),
            div(
              style = "margin-bottom: 10px;",
              tags$img(
                src = "do_plot.png",
                width = "40%",
                style = "object-fit: contain; border: none; display: block; margin-left: auto; margin-right: auto;"
              )
            ),
            div(
              style = "background-color: #f5f5f5; padding: 15px; border-radius: 5px; margin-bottom: 20px;",
              h4("Column Definitions"),
              p(HTML("<strong>disease:</strong> Human disease term")),
              p(HTML("<strong>log2foldchange:</strong> log2 fold change of term enrichment in localised genes vs genome background")),
              p(HTML("<strong>adj_pvalue:</strong> Bonferroni adjusted p-value of term enrichment")),
              p(HTML("<strong>count_in_genome:</strong> Total number of genes associated with term in genome background")),
              p(HTML("<strong>count_in_localised_genes:</strong> Number of glial localised genes associated with term")),
              p(HTML("<strong>dmel_gene_list:</strong> List of Drosophila genes associated with the disease term"))
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
            p("The table below contains Reactome enrichment analysis for each term and the list of associated genes."),
            hr(),
            div(
              style = "margin-bottom: 10px;",
              tags$img(
                src = "reactome_plot.png",
                width = "50%",
                style = "object-fit: contain; border: none; display: block; margin-left: auto; margin-right: auto;"
              )
            ),
            div(
              style = "background-color: #f5f5f5; padding: 15px; border-radius: 5px; margin-bottom: 20px;",
              h4("Column Definitions"),
              p(HTML("<strong>name:</strong> Reactome term")),
              p(HTML("<strong>log2foldchange:</strong> log2 fold change of term enrichment in localised genes vs genome background")),
              p(HTML("<strong>adj_pvalue:</strong> Bonferroni adjusted p-value of term enrichment")),
              p(HTML("<strong>count_in_genome:</strong> Total number of genes associated with term in genome background")),
              p(HTML("<strong>count_in_localised_genes:</strong> Number of glial localised genes associated with term")),
              p(HTML("<strong>dmel_gene_list:</strong> List of Drosophila genes associated with the Reactome term"))
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
