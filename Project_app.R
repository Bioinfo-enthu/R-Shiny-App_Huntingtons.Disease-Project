library(shiny)
library(DT)
library(tidyverse)
library(dplyr)
library(matrixStats)
library(gplots)
library(ggplot2)
library(gridExtra)
library(data.table)
library(beeswarm)
library(DESeq2)
library(EnhancedVolcano)
library(igraph)
library(colourpicker)
library(matrixStats)
library(viridis)
choices_list <- c("csv", "tsv", "text")

ui <- fluidPage(
  navbarPage("BF591: Final Project - Post-mortem Huntington’s Disease prefrontal cortex compared with neurologically healthy controls",
             tabPanel("Samples",
                      tabsetPanel(
                        tabPanel("Summary", DTOutput("sample_summary_table")),
                        tabPanel("Metadata", DTOutput("sample_info_table")),
                        tabPanel("Plots", 
                                 radioButtons("select_x", "Select X-axis:", choices = c("PMI", "Age of Death", "RIN", "mRNA-Seq reads")),
                                 plotOutput("sample_plots")
                        ),
                        fileInput("file", "Upload a file"),
                        conditionalPanel(
                          condition = "input.file != null"
                        )
                      )
             ),
             tabPanel("Counts",
                      fileInput("file1", "Upload file"),
                        tabsetPanel(
                          sliderInput("variance_slider", "Select the magnitude of genes with at least X percentile of variance:", min = 0, max = 100, value = 50),
                          sliderInput("non_zero_slider", "Select the magnitude of genes with at least X samples with non-zero values:", min = 0, max = 100, value = 50, step = 5),
                          tabPanel("Filter summary", dataTableOutput("counts_summary_table")),
                          tabPanel("Scatter Plot", plotOutput("scatter_plot_counts")),
                          tabPanel("HeatMap", plotOutput("heat_map_counts")),
                          tabPanel("PCA",
                                   sliderInput("selected_pcs", "Select top PC's to plot", min = 2, max = 50, value = 2),
                                   plotOutput("pca_scatter_counts"),
                                   DTOutput("pca_table")
                          
                        )
                      )
                      
             ),
             tabPanel("Differential Expression",
                      fileInput("csv_de", "Upload CSV file"),
                      
                        tabsetPanel(
                          tabPanel("DESeq2 Results", DTOutput("de_table")),
                          tabPanel("Volcano Plot", plotOutput("volcano_plot"),
                                   radioButtons("x_axis", "Choose the column for x axis", 
                                                choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
                                   ),
                                   radioButtons("y_axis", "Choose the column for y axis", 
                                                choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
                                   ),
                                   colourInput("color1", "Base point color:", value = "#DB349E"),
                                   colourInput("color2", "Highlight point color", value = "#3CDAE8"),
                                   
                          
                        )
                      )
             ),
             tabPanel("Correlation Network Analysis",
                      fileInput("csv_counts1", "Upload CSV file"),
                    
                        textInput("selected_genes", "Enter the gene ids of genes of your choice (comma-separated):"),
                        sliderInput("correlation_slider", "Select the minimum correlation for drawing an edge.", min = 0, max = 1, value = 0.5),
                        tabsetPanel(
                          tabPanel("Clustered HeatMap", plotOutput("clustered_heatmap")),
                          tabPanel("Correlation Network", plotOutput("correlation_network")),
                          tabPanel("Metrics", dataTableOutput("information_table"))
                        )
                      
             )
  )
)


server <- function(input, output) {
  
  #SAMPLE
  # Load datafile reactive function
  load_file <- reactive({
    req(input$file)
    dataf <- read.csv(input$file$datapath, header=TRUE)
    return(dataf)
  })
  
  # Summary table reactive function
  summarytable <- reactive({
    df <- load_file()
    summary_df <- data.frame(
      Column_name = names(df),
      Type = sapply(df, class),
      Mean = sapply(df, function(x) if (is.numeric(x)) mean(x, na.rm = TRUE) else NA)
    )
    return(summary_df)
  })
  
  # Render summary table
  output$sample_summary_table <- renderDT({
    summary_df <- summarytable()
    datatable(summary_df)
  })
  
  # Render metadata table
  output$sample_info_table <- renderDT({
    df <- load_file()
    datatable(df)
  })
  
  # Render plots
  output$sample_plots <- renderPlot({
    df <- load_file()
    
    # Check which variable is selected for the X-axis
    x_variable <- switch(input$select_x,
                         "PMI" = df$PMI,
                         "Age of Death" = df$Age_of_Death,
                         "RIN" = df$RIN,
                         "mRNA-Seq reads" = df$mRNA_Seq_reads)
    
    # Create the plot based on the selected variable
    ggplot(df, aes(x = x_variable)) +
      geom_histogram(bins = 50, alpha = 0.5, position = "identity", fill = "#995688", color = "#000000") +
      labs(x = input$select_x, y = "Count", title = paste("Histogram of", input$select_x))
  })
  
  #COUNTS  
  load_counts <- eventReactive(input$file1, {
    # Assuming your file is space-separated. If it's comma-separated, use sep = ","
    df <- read.table(input$file1$datapath, header = TRUE, sep = "\t")
   
    return(df)
  })
  
  
  load_counts_de <- eventReactive(input$csv_de, {
      df <- read.delim(input$csv_de$datapath)
      return(df)
  })
  
  output$de_table<- renderDT({
    load_counts_de()
  })
  
  output$volcano_plot <- renderPlot({
    create_volcano_plot()
  })
  
  filtered_genes <- reactive({
    counts_data <- load_counts()
    
    non_zero_threshold <- input$non_zero_slider
    variance_threshold <- input$variance_slider
    
    # Assuming your gene names are in the first column; adjust if needed
    gene_variances <- data.frame(Gene = counts_data$X, Variance = apply(counts_data[, -1], 1, var, na.rm = TRUE))
    
    # Filter genes based on variance
    genes_filtered_by_variance <- filter(gene_variances, Variance > quantile(Variance, variance_threshold / 100))
    
    # Filter genes based on non-zero sample counts
    non_zero_counts <- rowSums(counts_data[, -1] > 0, na.rm = TRUE)
    non_zero_genes <- filter(gene_variances, non_zero_counts >= non_zero_threshold)
    
    # Combine both filters
    final_filtered_genes <- intersect(genes_filtered_by_variance$Gene, non_zero_genes$Gene)
    
    # Return the filtered genes
    return(counts_data[counts_data$X %in% final_filtered_genes, ])
  })
  
  create_heatmap <- reactive({
  
    counts_data <- load_counts()
    filtered_data <- filtered_genes()
    genes_to_include <- filtered_data$X
    
    # Assuming your counts_data has a column for gene names (adjust if needed)
    counts_data_filtered <- counts_data[counts_data$X %in% genes_to_include, ]
    
    # Assuming "gene_column" is the name of the column in counts_data with gene names
    log_counts <- log2(counts_data_filtered[, -which(colnames(counts_data_filtered) %in% "X")] + 1)
    
    
    heatmap_output <- heatmap.2(as.matrix(log_counts),
                                col = viridis::viridis(75),
                                scale = "row",
                                cluster.cols = FALSE,
                                main = "Clustered Heatmap of Selected Genes",
                                labRow = counts_data_filtered$X)
    
    return(heatmap_output)
  })
  
  # Inside the create_pca_plot reactive function
  create_pca <- reactive({
    
    filtered_data <- filtered_genes()
    
    # Assuming "Gene" is the first column; adjust if needed
    filtered_data <- filtered_data
    
    # Remove the 'Gene' column
    filtered_data <- filtered_data[, -1]
    
    transposed_data <- transpose(data.frame(filtered_data))
    
    # Perform PCA directly using prcomp
    pca_result <- prcomp(transposed_data, scale. = TRUE)
    
    # Extract the top PCs based on the slider input
    selected_pcs <- seq_len(input$selected_pcs)
    
    # Initialize a list to store plots and tables for each PC
    pc_plots <- list()
    pc_tables <- list()
    
    for (pc in selected_pcs) {
      # Extract the specific PC
      top_pc <- pca_result$x[, pc]
      
      # Add % variance explained to the column names
      var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 2)
      label <- paste0("PC", pc, " (", var_explained[pc], "%)")
      
      # Create scatter plot
      beeswarm_plot <- beeswarm(
        as.matrix(top_pc),
        pch = 16,
        col = viridis::viridis(1, option = "D"),
        main = paste("Beeswarm Plot of PC", pc)
      )
      
      # Create data frame for the table
      table_data <- data.frame(
        "PC" = paste0("PC", pc),
        "Percent Variance Explained" = var_explained[pc]
      )
      
      # Add plots and tables to the lists
      pc_plots[[paste("pc_plot", pc, sep = "_")]] <- beeswarm_plot
      pc_tables <- rbind(pc_tables, data.frame(PC = label, "Percent Variance Explained" = var_explained[pc]))
    }
    
    return(list(pc_plots = pc_plots, pc_tables = pc_tables))
  })
  
  
  
  create_diagnostics <- reactive({
    counts_data <- load_counts()
    filtered_data <- filtered_genes()
    log_transform <- function(x) log2(x + 1)
    
    # Assuming your gene names are in the first column; adjust if needed
    counts_data$Gene <- counts_data$X
    
    # Assuming that numeric columns start from the second column; adjust if needed
    gene_column <- "Gene"
    
    # Assuming that numeric columns start from the second column; adjust if needed
    numeric_cols <- sapply(counts_data[, -which(colnames(counts_data) %in% gene_column)], is.numeric)
    
    if (sum(numeric_cols) > 0) {
      counts_data$median_count <- rowMedians(as.matrix(counts_data[, numeric_cols]), na.rm = TRUE)
      counts_data$log_variance <- log_transform(apply(as.matrix(counts_data[, numeric_cols]), 1, var, na.rm = TRUE))
    } else {
      # Handle the case where there are no numeric columns
      counts_data$median_count <- NA
      counts_data$log_variance <- NA
    }
    
    # Check if filtered_data has the "Gene" column
    counts_data$color <- ifelse(counts_data$Gene %in% filtered_data$X, "Passing Filter", "Filtered Out")
    
    # Calculate num_zeros
    counts_data$num_zeros <- rowSums(as.matrix(counts_data[, numeric_cols]) == 0, na.rm = TRUE)
    
    # Set a more visually appealing color palette
    my_palette <- c("Filtered Out" = "#CCCCCC", "Passing Filter" = "#3498db")
    
    p1 <- ggplot(counts_data, aes(x = median_count, y = log_variance, color = color)) +
      geom_point(size = 3, alpha = 0.7, shape = 16) +
      scale_color_manual(values = my_palette) +
      labs(title = "Median Count vs Variance", x = "Median Count", y = "Log Variance") +
      theme_minimal() +
      theme(legend.position = "right", plot.title = element_text(hjust = 0.5))
    
    p2 <- ggplot(counts_data, aes(x = median_count, y = num_zeros, color = color)) +
      geom_point(size = 3, alpha = 0.7, shape = 16) +
      scale_color_manual(values = my_palette) +
      labs(title = "Median Count vs Number of Zeros", x = "Median Count", y = "Number of Zeros") +
      theme_minimal() +
      theme(legend.position = "right", plot.title = element_text(hjust = 0.5))
    
    return(list(p1, p2))
  })
  
  output$counts_summary_table <- renderDataTable({
    counts_data <- load_counts()
    filtered_genes_data <- filtered_genes()
    
  
    
    total_genes <- nrow(counts_data)
    genes_passing_filter <- nrow(filtered_genes_data)
    genes_not_passing_filter <- total_genes - genes_passing_filter
    
    summary_df <- data.frame(
      "Variables" = c(
        "Number of Samples",
        "Total Number of Genes",
        "Number of Genes Passing Filter",
        "Percent of Genes Passing Filter",
        "Number of Genes Not Passing Filter",
        "Percent of Genes Not Passing Filter"
      ),
      "Values" = c(
        ncol(counts_data),
        total_genes,
        genes_passing_filter,
        round(100 * (genes_passing_filter / total_genes), 2),
        genes_not_passing_filter,
        round(100 * (genes_not_passing_filter / total_genes), 2)
      )
    )
    
    datatable(summary_df, caption = "Summary Table")
  })
  output$scatter_plot_counts <- renderPlot({
    scatter_plots <- create_diagnostics()
    grid.arrange(grobs = scatter_plots, ncol = 2)
  })
  
  output$heat_map_counts <- renderPlot({
    create_heatmap()
  })
  
  output$pca_scatter_counts <- renderPlot({
    create_pca()$combined_beeswarm_plot
  })
  
  output$pca_table <- renderDT({
    create_pca()$pc_tables
  })
  
  create_volcano_plot <- reactive({
    # Check if DE results are available
    de_results <- load_counts_de()
    
    # Get selected columns for x and y axis from radio buttons
    x_axis_col <- input$x_axis
    y_axis_col <- input$y_axis
    
    # Get color inputs
    base_color <- input$color1
    highlight_color <- input$color2
    
    # Create a volcano plot
    volcano_plot <- EnhancedVolcano(
      de_results,
      lab = de_results$symbol,
      x = x_axis_col,
      y = y_axis_col,
      xlim = c(-5, 5),
      title = 'Volcano Plot of Differential Expression',
      subtitle = 'Adjusted p-value vs. Log2 Fold Change',
      pCutoff = 0.05,
      FCcutoff = 2,
      pointSize = 3,
      labSize = 2.5,
      col = c("gray", base_color, highlight_color),
      legendLabels = c('Not Significant', 'Significant Up', 'Significant Down'),
      legendPosition = 'right'
    )
    
    return(volcano_plot)
  })
  
 
  load_counts_cna <- eventReactive(input$csv_counts1, {
   
      df <- read.delim(input$csv_counts1$datapath)
      print(df)
      return(df)
    }
  )
  create_correlation_network <- reactive({
    counts_data <- load_counts_cna()  # Assuming you want to load count data
    
    # Assuming your counts_data has a column for gene names (adjust if needed)
    selected_genes <- tools::toTitleCase(strsplit(input$selected_genes, ",")[[1]])
    print(selected_genes)
    not_found_genes <- setdiff(selected_genes, counts_data$X)
    if (length(not_found_genes) > 0) {
      stop(paste("Genes '", paste(not_found_genes, collapse = "', '"), "' not found in the counts data.", sep = ""))
    }
    
    counts_data_selected <- counts_data[counts_data$X %in% selected_genes, ]
    gene_ids_selected <- counts_data_selected$X
    threshold <- input$correlation_slider
    
    correlation_matrix <- cor(t(counts_data_selected[, -c(1)]))
    
    correlation_matrix[abs(correlation_matrix) < threshold] <- 0
    correlation_matrix[is.na(correlation_matrix)] <- 0
    
    graph <- graph.adjacency(correlation_matrix, weighted = TRUE, mode = "directed")
    
    V(graph)$name <- gene_ids_selected
    
    constant_color <- "pink"
    V(graph)$color <- constant_color
    
    p <- plot(
      graph,
      layout = layout_with_fr(graph),
      vertex.label = V(graph)$name,
      vertex.color = V(graph)$color,
      main = "Correlation Network"
    )
    return(list(p_obj = p, graph_obj = graph))
  })
  
  create_graph_information_table <- reactive({
    graph <- create_correlation_network()$graph_obj
    
    degree <- degree(graph)
    closeness <- closeness(graph)
    betweenness <- betweenness(graph)
    
    information_df <- data.frame(
      Gene = V(graph)$name,
      Degree = degree,
      ClosenessCentrality = closeness,
      BetweennessCentrality = betweenness
    )
    return(information_df)
  })
  
  create_clustered_heatmap <- reactive({
   
    counts_data <- load_counts_cna()
    
    selected_genes <- tools::toTitleCase(strsplit(input$selected_genes, ",")[[1]])
    
    not_found_genes <- setdiff(selected_genes, counts_data$X)
    if (length(not_found_genes) > 0) {
      stop(paste("Genes '", paste(not_found_genes, collapse = "', '"), "' not found in the counts data.", sep = ""))
    }
    
    counts_data_selected <- counts_data[counts_data$X %in% selected_genes, ]
    
    # Assuming "gene_column" is the name of the column in counts_data with gene names
    log_counts <- log2(counts_data_selected[, -c(1)] + 1)
    
    # Create clustered heatmap
    heatmap.2(as.matrix(log_counts),
              col = viridis::viridis(75),
              scale = "row",
              cluster.cols = FALSE,
              main = "Clustered Heatmap of Selected Genes",
              labRow = counts_data_selected$GeneID)
  })
  
  output$information_table <- renderDataTable({
    datatable(create_graph_information_table(), caption = "Metrics for Each Gene in the Graph Network")
  })
  
  output$correlation_network <- renderPlot({
    create_correlation_network()$p_obj
  })
  
  output$clustered_heatmap <- renderPlot({
    create_clustered_heatmap()
  })
  
  }



options(shiny.maxRequestSize = 50 * 1024^2)
shinyApp(ui,server)
