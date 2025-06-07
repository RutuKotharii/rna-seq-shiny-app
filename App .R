# Loading required libraries 
library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(pheatmap)
library(tidyr)

# Sourcing helper functions from an external file
source("functions.R")

# Increasing the maximum file upload size to handle large RNA-Seq datasets
options(shiny.maxRequestSize = 100 * 1024^2)
options(shiny.maxRequestSize = 100 * 1024^2)

# Defining the User Interface (UI) for the application
ui <- dashboardPage(
  # Creating the dashboard header with a title
  dashboardHeader(title = "RNA-Seq Analysis App"),
  
  # Adding a sidebar with navigation options for different tabs
  dashboardSidebar(
    sidebarMenu(
      menuItem("Sample Information", tabName = "sample_info", icon = icon("table")),
      menuItem("Counts Matrix", tabName = "counts_matrix", icon = icon("database")),
      menuItem("Differential Expression", tabName = "differential_expression", icon = icon("chart-line")),
      menuItem("Gene Set Enrichment", tabName = "gene_set", icon = icon("dna"))
    )
  ),
  
  # Adding content to the dashboard body with tab-specific layouts
  dashboardBody(
    tabItems(
      # Sample Information Tab
      tabItem(tabName = "sample_info",
              h2("Sample Information Exploration"),
              # Allowing the user to upload a sample information file
              fileInput("sample_info_file", "Upload Sample Info CSV"),
              tabsetPanel(
                # Displaying a summary of the sample data
                tabPanel("Summary", DTOutput("summary_table")),
                # Showing the sample data in a table format
                tabPanel("Data Table", DTOutput("data_table")),
                # Providing a custom visualization option
                tabPanel("Visualization",
                         selectInput("x_col", "Select a Continuous Variable:", choices = NULL),
                         selectInput("group_col", "Select a Grouping Variable:", choices = NULL),
                         plotOutput("custom_plot")
                )
              )
      ),
      
      # Counts Matrix Tab
      tabItem(tabName = "counts_matrix",
              h2("Counts Matrix Exploration"),
              # Allowing the user to upload a counts matrix file
              fileInput("counts_file", "Upload Counts Matrix CSV"),
              # Adding sliders for variance percentile and minimum non-zero samples
              sliderInput("var_threshold", "Variance Percentile", min = 0, max = 1, value = 0.8),
              sliderInput("min_nonzero", "Minimum Non-Zero Samples", min = 1, max = 69, value = 10),
              tabsetPanel(
                # Displaying a filtered summary of the counts matrix
                tabPanel("Filtered Summary", DTOutput("counts_summary_table")),
                # Adding diagnostic plots for exploratory analysis
                tabPanel("Diagnostic Plots",
                         selectInput("x_axis", "Select X-axis Metric:", choices = c("Median")),
                         selectInput("y_axis", "Select Y-axis Metric:", choices = c("Variance", "Zeros")),
                         checkboxInput("log_transform", "Log-transform Axes", value = TRUE),
                         plotOutput("diagnostic_plot")
                ),
                # Creating a heatmap visualization for top genes
                tabPanel("Heatmap",
                         checkboxInput("heatmap_log_transform", "Log-transform Counts", value = TRUE),
                         selectInput("color_palette", "Color Palette", 
                                     choices = c("Viridis", "Blue-Red", "Heat", "Coolwarm"), 
                                     selected = "Viridis"),
                         numericInput("top_genes", "Number of Top Genes (by Variance)", 
                                      value = 100, min = 10, max = 5000, step = 50),
                         plotOutput("heatmap")
                ),
                # Adding a PCA visualization tab
                tabPanel("PCA",
                         numericInput("top_genes_pca", "Number of Top Genes (by Variance)", value = 100, min = 10, step = 5),
                         selectInput("pc_x", "PC X", choices = NULL, selected = NULL),
                         selectInput("pc_y", "PC Y", choices = NULL, selected = NULL),
                         radioButtons("pca_plot_type", "Plot Type", 
                                      choices = c("Scatter Plot" = "scatter", "Beeswarm Plot" = "beeswarm"), 
                                      selected = "scatter"),
                         selectInput("num_components", "Number of PCA Components", choices = 5, selected = 5),
                         plotOutput("pca_plot")
                )
              )
      ),
      
      # Differential Expression Tab
      tabItem(tabName = "differential_expression",
              h2("Differential Expression Analysis"),
              # Allowing the user to upload a differential expression results file
              fileInput("de_file", "Upload Differential Expression Results (CSV):"),
              # Adding sliders for log2 fold change and adjusted p-value thresholds
              sliderInput("log2fc_threshold", "Log2 Fold Change Threshold:", min = 0, max = 5, value = 1.5, step = 0.05),
              sliderInput("sig_padj_threshold", "Adjusted P-Value Threshold:", min = 0, max = 0.1, value = 0.05, step = 0.01),
              tabsetPanel(
                # Displaying the results in a table
                tabPanel("Results Table", 
                         DTOutput("de_table")),
                # Showing histograms for p-values and log2 fold changes
                tabPanel("P-Value Histogram", 
                         plotOutput("pvalue_histogram")),
                tabPanel("Log2 Fold Change Histogram", 
                         plotOutput("logfc_histogram")),
                # Creating a volcano plot visualization
                tabPanel("Volcano Plot", 
                         plotOutput("volcano_plot")),
                # Displaying normalized counts for top genes
                tabPanel("Top Genes Normalized Counts",
                         plotOutput("top_genes_plot"))
              )
      ),
      
      # Gene Set Enrichment Tab
      tabItem(
        tabName = "gene_set",
        h2("Gene Set Enrichment Analysis"),
        # Allowing the user to upload FGSEA results
        fileInput("fgsea_upload", "Upload FGSEA Results (CSV/TSV)", accept = c(".csv", ".tsv")),
        # Adding sliders to filter pathways and thresholds
        sliderInput("top_n", "Number of Top Pathways to Display:", min = 5, max = 50, value = 10),
        sliderInput("padj_filter", "Adjusted p-value Threshold:", min = 0, max = 1, value = 0.05, step = 0.01),
        radioButtons("nes_filter", "Pathways by NES:", 
                     choices = c("All" = "all", "Positive" = "pos", "Negative" = "neg")),
        downloadButton("download_filtered", "Download Filtered Results"),
        
        tabsetPanel(
          # Displaying the top pathways as a barplot and table
          tabPanel("Top Pathways",
                   plotOutput("barplot_top", height = "500px"),
                   DTOutput("table_top")),
          # Showing filtered results in a table
          tabPanel("Filtered Table",
                   DTOutput("table_filtered")),
          # Adding a scatter plot visualization
          tabPanel("Scatter Plot",
                   plotOutput("scatter_plot", height = "500px"))
        )
      )
              
      )
    
  )
)

# Defining the server logic
server <- function(input, output, session) {
  # Loading and processing sample data upon file upload
  sample_data <- reactive({
    req(input$sample_info_file)
    load_sample_data(input$sample_info_file$datapath)
  })
  
  # Rendering a summary table for sample information
  output$summary_table <- renderDT({
    req(sample_data())
    summary_df <- generate_summary(sample_data())
    datatable(summary_df, options = list(scrollX = TRUE))
  })
  
  # Rendering a data table for displaying the uploaded sample information
  output$data_table <- renderDT({
    req(sample_data())
    datatable(sample_data(), options = list(pageLength = 25, scrollX = TRUE))
  })
  
  # Updating dropdown choices for columns
  observe({
    req(sample_data())
    data <- sample_data()
    
    # Restricting options for grouping and continuous variables
    grouping_choices <- c("Diagnosis", "Vonsattel_Grade")
    continuous_choices <- c("pmi", "Age_of_Death", "rin", "mrna_seq_reads", 
                            "Age_of_Onset", "Duration", "cag", 
                            "hv_Striatal_Score", "Cortical_Score")
    
    updateSelectInput(session, "x_col", choices = continuous_choices)
    updateSelectInput(session, "group_col", choices = grouping_choices)
  })
  
  # Rendering a custom plot for exploring sample information
  output$custom_plot <- renderPlot({
    req(sample_data(), input$x_col, input$group_col)
    create_density_plot(sample_data(), input$x_col, input$group_col)
  })
  
  #PART2
  # Loading and processing the counts matrix upon file upload
  counts_data <- reactive({
    req(input$counts_file)
    load_counts_matrix(input$counts_file$datapath)
  })
  
  # Rendering a filtered summary table of the counts matrix
  filtered_data <- reactive({
    req(counts_data())
    filter_counts(counts_data(), input$var_threshold, input$min_nonzero)
  })
  
  output$counts_summary_table <- renderDT({
    req(filtered_data())
    datatable(filtered_data()$stats, rownames = FALSE, options = list(dom = 't'))
  })
  
  # Rendering a diagnostic plot for exploratory data analysis of the counts matrix
  output$diagnostic_plot <- renderPlot({
    req(filtered_data())
    counts <- counts_data()
    filtered_indices <- rownames(filtered_data()$filtered_counts)
    
    num_zeros <- ncol(counts) - rowSums(counts > 0)
    
    create_diagnostic_plot(
      counts,
      input$x_axis,
      input$y_axis,
      rownames(counts) %in% filtered_indices,
      input$log_transform,
      num_zeros = num_zeros
    )
  })
  
  # Rendering a heatmap for the top variable genes in the counts matrix  
output$heatmap <- renderPlot({
  req(filtered_data())
  generate_heatmap(
    filtered_counts = filtered_data()$filtered_counts,
    top_genes = input$top_genes,
    log_transform = input$heatmap_log_transform,
    color_palette = input$color_palette
  )
})
  
# Reactive function for calculating PCA
pca_result <- reactive({
  req(filtered_data())
  
  # Performing PCA on filtered counts with the top selected genes
  perform_pca(filtered_data()$filtered_counts, top_genes = input$top_genes_pca)
})

# Updating PCA axis options dynamically based on the number of components
observe({
  req(pca_result())
  num_pcs <- length(pca_result()$pca$sdev)
  pc_choices <- paste0("PC", 1:num_pcs)
  
  updateSelectInput(session, "pc_x", choices = pc_choices, selected = "PC1")
  updateSelectInput(session, "pc_y", choices = pc_choices, selected = "PC2")
})

# Rendering PCA plot (Scatter or Beeswarm)
output$pca_plot <- renderPlot({
  req(pca_result(), input$pc_x)
  
  if (input$pca_plot_type == "scatter") {
    # Rendering scatter plot
    req(input$pc_y)  # Ensure both axes are selected
    create_pca_scatter(
      pca_result = pca_result()$pca,
      pc_x = input$pc_x,
      pc_y = input$pc_y,
      explained_variance = pca_result()$explained_variance
    )
  } else if (input$pca_plot_type == "beeswarm") {
    # Rendering beeswarm plot for multiple components
    components <- paste0("PC", 1:input$num_components)
    create_pca_beeswarm_multiple(
      pca_result = pca_result()$pca,
      components = components,
      explained_variance = pca_result()$explained_variance
    )
  }
})

#PART3
# Loading and processing differential expression results upon file upload
de_results <- reactive({
  req(input$de_file)
  read.csv(input$de_file$datapath) %>%
    process_deseq2_results()
})

# Rendering the results table for differential expression analysis
output$de_table <- renderDT({
  de_results() %>%
    datatable(
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
})

# Rendering a histogram of p-values for differential expression results
output$pvalue_histogram <- renderPlot({
  plot_pvalue_histogram(de_results())
})

# Rendering a histogram of log2 fold changes for differential expression results
output$logfc_histogram <- renderPlot({
  plot_logfc_histogram(de_results())
})

# Rendering a volcano plot for visualizing significant genes
output$volcano_plot <- renderPlot({
  plot_volcano(de_results())
})

# Rendering normalized counts for top genes identified from differential expression results
output$top_genes_plot <- renderPlot({
  req(de_results())
  topGenesData <- de_results()[order(de_results()$padj), ][1:10, ]
  topGenes <- topGenesData$Gene_ID
  plotData <- data.frame(
    Gene = rep(topGenes, each = 2),
    SampleType = rep(c("Control", "Huntington's Disease"), times = 10),
    NormCount = c(topGenesData$Control_mean, topGenesData$HD_mean)
  )
  
  ggplot(plotData, aes(x = Gene, y = log10(NormCount + 1), color = SampleType)) +
    geom_point(size = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Plot of Log10(normalized counts) for top ten DE genes", 
         x = "Gene", y = "log10(Normalized Counts)", color = "Sample Type")
})


#PART4
# Loading and processing FGSEA results upon file upload
fgsea_data <- reactive({
  req(input$fgsea_upload)
  ext <- tools::file_ext(input$fgsea_upload$name)
  if (ext == "csv") {
    read.csv(input$fgsea_upload$datapath)
  } else if (ext == "tsv") {
    read.delim(input$fgsea_upload$datapath)
  } else {
    stop("Please upload a .csv or .tsv file.")
  }
})

# Rendering a barplot of the top enriched pathways from FGSEA results
output$barplot_top <- renderPlot({
  req(fgsea_data())
  plot_top_pathways_by_nes(fgsea_data(), input$top_n)  # Using updated function
})

# Rendering a table of the top enriched pathways from FGSEA results
output$table_top <- renderDT({
  req(fgsea_data())
  fgsea_data() %>%
    arrange(desc(abs(NES))) %>%  # Sorting by absolute NES
    head(input$top_n) %>%
    datatable(options = list(pageLength = 10, scrollX = TRUE))
})

# Rendering a scatter plot for pathway enrichment analysis
filtered_results <- reactive({
  req(fgsea_data())
  filter_fgsea_results(fgsea_data(), input$padj_filter, input$nes_filter)
})

# Rendering a filtered table of FGSEA results
output$table_filtered <- renderDT({
  req(filtered_results())
  datatable(filtered_results(), options = list(pageLength = 10, scrollX = TRUE))
})

# Downloading filtered results
output$download_filtered <- downloadHandler(
  filename = function() {
    paste0("filtered_fgsea_results_", Sys.Date(), ".csv")
  },
  content = function(file) {
    save_filtered_results(filtered_results(), file)
  }
)

# Scatter plot
output$scatter_plot <- renderPlot({
  req(fgsea_data())
  plot_scatter(fgsea_data(), input$padj_filter)
})

}

# Running the application 
shinyApp(ui = ui, server = server)
