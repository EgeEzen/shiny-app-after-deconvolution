# app.R

library(shiny)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(ggplot2)
library(apeglm)
library(writexl)
library(ggrepel)
library(tidyr)
library(limma)
library(dplyr)
library(ggsignif)
library(tibble)
library(envalysis)

#> shiny::runApp("~/Desktop/PhD/2024 Second Term/shiny app after deconvolution")
#> rsconnect::deployApp("~/Desktop/PhD/2024 Second Term/shiny app after deconvolution",appFileManifest = "manifest.txt")

# Load data (adjust paths as necessary)
count_data <- readRDS("count_data_pauci_immune.rds")
metadata   <- readRDS("metadata_pauci_immune.rds")
nrm_counts <- readRDS("nrm_counts_pauci_immune.rds")

# Get a list of genes from the normalized counts matrix
gene_list <- rownames(nrm_counts)

ui <- fluidPage(
  titlePanel("Differential Expression Analysis by Selected Gene"),
  sidebarLayout(
    sidebarPanel(
      selectizeInput("selected_gene", "Select Gene Name:", choices = NULL, multiple = FALSE, 
                     options = list(placeholder = "Type to search gene...")),
      actionButton("run_analysis", "Run Analysis"),
      hr(),
      h4("Download Plot Settings:"),
      actionButton("toggleSettings", "Toggle Settings", icon = icon("cogs")),
      conditionalPanel(
        condition = "input.toggleSettings % 2 == 1",
        numericInput("plot_width", "Plot Width (inches)", value = 7, min = 3, max = 20, step = 0.5),
        numericInput("plot_height", "Plot Height (inches)", value = 7, min = 3, max = 20, step = 0.5)
      ),
      hr(),
      downloadButton("download_DE_results", "Download DE Results (Excel)"),
      downloadButton("download_nrm_counts","Download Normalized Counts")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("PCA Plot", 
                 plotOutput("pcaPlot"),
                 downloadButton("download_pca", "Download PCA Plot")),
        tabPanel("Volcano Plot", 
                 plotOutput("volcanoPlot"),
                 downloadButton("download_volcano", "Download Volcano Plot")),
        tabPanel("Heatmap", 
                 plotOutput("heatmapPlot"),
                 downloadButton("download_heatmap", "Download Heatmap")),
        tabPanel("Dot Plot", 
                 plotOutput("dotPlot"),
                 selectizeInput("dot_gene", "Select Gene for Dotplot:", choices = NULL, multiple = FALSE, 
                                options = list(placeholder = "Type to search gene...")),
                 actionButton("draw_dotplot", "Draw Dot Plot"),
                 downloadButton("download_dot", "Download Dot Plot"),
                 downloadButton("download_csv", "Download CSV"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  updateSelectizeInput(session, "selected_gene", choices = gene_list, server = TRUE)
  
  analysis_results <- eventReactive(input$run_analysis, {
    req(input$selected_gene)
    selected_gene <- input$selected_gene
    
    #quantile 
    selected_gene_expr <- nrm_counts[selected_gene, ]
    q25 <- quantile(selected_gene_expr, 0.25)
    q75 <- quantile(selected_gene_expr, 0.75)
    
    # groups
    metadata$selected_gene_expression <- ifelse(
      selected_gene_expr <= q25, "Low",
      ifelse(selected_gene_expr >= q75, "High", NA)
    )
    
    # subset metadata and countdata
    metadata_selected <- metadata[!is.na(metadata$selected_gene_expression), ]
    count_data_selected <- count_data[, colnames(count_data) %in% rownames(metadata_selected)]
    
    # create deseq2 object, sex as a covariant
    dds <- DESeqDataSetFromMatrix(
      countData = count_data_selected,
      colData   = metadata_selected,
      design    = ~ sex + selected_gene_expression
    )
    # low as a reference
    dds$selected_gene_expression <- relevel(dds$selected_gene_expression, ref = "Low")
    dds <- DESeq(dds)

    # shrink log fold change
    res_vs <- lfcShrink(dds, coef = "selected_gene_expression_High_vs_Low", type = "apeglm")
    
    # for visualization, batch effect remover
    vsd <- vst(dds, blind = FALSE)
    vsd_copy <- vsd
    mat_copy <- assay(vsd_copy)
    mm <- model.matrix(~ selected_gene_expression, colData(vsd_copy))
    mat_copy <- removeBatchEffect(mat_copy, batch = vsd_copy$sex, design = mm)
    assay(vsd_copy) <- mat_copy
    
    # pca data
    pcaData <- plotPCA(vsd_copy, intgroup = "selected_gene_expression", returnData = TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    
    # prepare heatmap -> only top 50 genes(according to padj) and foldchange of 1.5
    resSig <- res_vs[order(res_vs$padj), ]
    resSig <- subset(resSig, padj <= 0.05)
    resSig_fc <- subset(resSig, abs(log2FoldChange) >= 0.585)
    resSig_df <- as.data.frame(resSig_fc)
    resSig_merged <- merge(resSig_df, as.data.frame(mat_copy), by = "row.names")
    rownames(resSig_merged) <- resSig_merged$Row.names
    heatmap_data <- resSig_merged[order(resSig_merged$padj), ]
    heatmap_data <- heatmap_data[1:min(50, nrow(heatmap_data)), 7:ncol(heatmap_data), drop = FALSE]
    
    list(
      selected_gene     = selected_gene,
      dds               = dds,
      res_vs            = res_vs,
      pcaData           = pcaData,
      percentVar        = percentVar,
      heatmap_data      = heatmap_data,
      metadata_selected = metadata_selected
    )
    })
  
  updateSelectizeInput(session, "dot_gene", choices = gene_list, server = TRUE)
  
  #for dotplot drawing
  observeEvent(input$draw_dotplot, { 
    #prepare df for dotplot
    dot_gene <- input$dot_gene
    nrm_counts_filtered <- nrm_counts[dot_gene,rownames(analysis_results()$metadata_selected),drop=FALSE]
    nrm_counts_filtered <- as.data.frame(t(nrm_counts_filtered))
    nrm_counts_filtered <- rownames_to_column(nrm_counts_filtered, var = "Sample")
    nrm_counts_filtered$Condition <- analysis_results()$metadata_selected$selected_gene_expression
    colnames(nrm_counts_filtered) <- c("Sample","Counts","Condition") 
    nrm_counts_filtered$Condition <- factor(nrm_counts_filtered$Condition, labels = c("Low","High"))
    
    res_vs <- analysis_results()$res_vs
    checking_condition <- (res_vs[rownames(res_vs) == dot_gene,,drop=FALSE ])$padj
    if (nrow(nrm_counts_filtered) == 0) {
      output$dotPlot <- renderPlot({
        ggplot() + 
          geom_text(aes(x = 1, y = 1, label = "Gene not found!"), size = 6, color = "red") +
          theme_void()
      })
    } else {
      output$dotPlot <- renderPlot({
        if (!is.na(checking_condition) && checking_condition <= 0.05){ # with significance
          ggplot(nrm_counts_filtered, aes(x = Condition, y = Counts, fill = Condition)) +
            geom_violin(trim = TRUE, alpha = 0.8) + 
            geom_jitter(width = 0.2, size = 1, alpha = 0.6, color = "black") +
            labs(title = paste0("Gene Expression of\n", dot_gene ," in Pauci-immune RA with\n", analysis_results()$selected_gene," Expression Groups"), 
                 x = "Condition", y = "Normalized Counts") +
            theme_publish() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(size = 12),  
                  legend.position = "none") +
            ylim(NA,max(nrm_counts_filtered$Counts) * 1.2) +
            geom_signif(comparisons = list(c("Low", "High")),
                                 annotations = "*",  # Star for significance
                                 y_position = max(nrm_counts_filtered$Counts) * 1.05, 
                                 tip_length = 0.02, textsize = 5)
        }
        else {#without significance
              ggplot(nrm_counts_filtered, aes(x = Condition, y = Counts, fill = Condition)) +
                geom_violin(trim = TRUE, alpha = 0.8) + 
                geom_jitter(width = 0.2, size = 1, alpha = 0.6, color = "black") +
                labs(title = paste0("Gene Expression of\n", dot_gene," in Pauci-immune RA with\n", analysis_results()$selected_gene," Expression Groups"), 
                     x = "Condition", y = "Normalized Counts") +
                theme_publish() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(size = 12),  
                      legend.position = "none") +
                ylim(NA,max(nrm_counts_filtered$Counts) * 1.2)
            }
      
      })
    }
    
    output$download_csv <- downloadHandler(
      filename = function() {
        paste0("dotplot_",input$dot_gene,"_in_",analysis_results()$selected_gene, ".png")
      },
      content = function(file) {
        write.csv(nrm_counts_filtered, file, row.names = FALSE)
      }
    )
    
    
    })
    
  #### plot outputs
  
  output$pcaPlot <- renderPlot({
    req(analysis_results())
    pcaData <- analysis_results()$pcaData
    percentVar <- analysis_results()$percentVar
    ggplot(pcaData, aes(PC1, PC2, color = selected_gene_expression)) +
      geom_point(size = 3) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      scale_color_brewer(palette = "Set2") +
      coord_fixed() +
      labs(title = paste0("PCA Plot on\n", analysis_results()$selected_gene, " Expression")) +
      theme_publish()
  })
  
  output$volcanoPlot <- renderPlot({
    req(analysis_results())
    res_vs_df <- as.data.frame(analysis_results()$res_vs)
    res_vs_df$padj[is.na(res_vs_df$padj)] <- 1
    threshold_for_lfc <- 0.585
    res_vs_df$Classification <- with(res_vs_df, ifelse(
      padj <= 0.05 & log2FoldChange >= threshold_for_lfc, "Upregulated",
      ifelse(padj <= 0.05 & log2FoldChange <= -threshold_for_lfc, "Downregulated", "Not Significant")
    ))
    res_vs_df$log_adjP <- -log10(res_vs_df$padj)
    
    label_data <- subset(res_vs_df, (log2FoldChange >= threshold_for_lfc | log2FoldChange <= -threshold_for_lfc) & Classification != "Not Significant" )
    
    top_25_up <- label_data %>% filter(Classification == "Upregulated") %>% arrange(desc(log2FoldChange)) %>% head(25)
    top_25_down <- label_data %>% filter(Classification == "Downregulated") %>% arrange(log2FoldChange) %>% head(25)
    label_data <- bind_rows(top_25_up, top_25_down)
    
    ggplot(res_vs_df[rownames(res_vs_df) != analysis_results()$selected_gene, ], aes(x = log2FoldChange, y = log_adjP)) +
      geom_point(aes(color = Classification), size = 1,alpha=0.8) + #alpha 0.5-0.8 looks "cool"
      scale_color_manual(values = c("Upregulated" = "#EE4E4E", "Downregulated" = "#2A629A", "Not Significant" = "lightgray")) +
      geom_vline(xintercept = c(-threshold_for_lfc,threshold_for_lfc), linetype = "dashed", color = "black") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
      labs(
        title = paste0(analysis_results()$selected_gene ,"-Differentially Expressed Genes\nWithin Pauci-immune RA"), 
        x = "Log2 Fold Change",
        y = "-log10(Adjusted P-Value)",
        fill = "Classification"
      ) +   # Increase plot margins if necessary
      ylim(0,NA) +
      geom_label_repel(
        data = label_data,
        aes(label = rownames(label_data)),
        box.padding = unit(0.1, "lines"),
        segment.color = "#2B2A4C",
        max.overlaps= 50
      ) +
      guides(color = guide_legend(title = "Differential Expression", override.aes = list(size = 3)),alpha= "none") +
      theme_publish()
  })
  
  output$heatmapPlot <- renderPlot({
    req(analysis_results())
    custom_colors <- colorRampPalette(c("blue", "white", "red"))(50)
    pheatmap(
      analysis_results()$heatmap_data,
      cluster_rows   = TRUE,
      show_rownames  = TRUE,
      scale          = "row",
      color          = custom_colors,
      show_colnames  = FALSE,
      cluster_cols   = TRUE,
      annotation_col = analysis_results()$metadata_selected["selected_gene_expression"],
      main           = paste0("Within Pauci-immune RA\n", analysis_results()$selected_gene, " Differential Expression")
    )
  })
  
  #### download
  
  # ggsave for plots
  output$download_pca <- downloadHandler(
    filename = function() {
      paste0(analysis_results()$selected_gene, "_PCA_plot.png")
    },
    content = function(file) {
      ggsave(file, plot = last_plot(), device = "png", dpi = 300, 
             bg = "transparent", width = input$plot_width, height = input$plot_height)
    }
  )
  
  output$download_volcano <- downloadHandler(
    filename = function() {
      paste0(analysis_results()$selected_gene, "_volcano_plot.png")
    },
    content = function(file) {
      ggsave(file, plot = last_plot(), device = "png", dpi = 300, 
             bg = "transparent", width = input$plot_width, height = input$plot_height)
    }
  )
  
  output$download_dot <- downloadHandler(
    filename = function() {
      paste0("dotplot_",input$dot_gene,"_in_",analysis_results()$selected_gene, ".png")
    },
    content = function(file) {
      ggsave(file, plot = last_plot(), device = "png", dpi = 300, 
             bg = "transparent", width = input$plot_width, height = input$plot_height)
    }
  )
  
  # For the heatmap (grid object from pheatmap), we use png() and dev.off().
  output$download_heatmap <- downloadHandler(
    filename = function() {
      paste0(analysis_results()$selected_gene, "_heatmap.png")
    },
    content = function(file) {
      png(file, width = input$plot_width * 100, height = input$plot_height * 100, bg = "transparent")
      custom_colors <- colorRampPalette(c("blue", "white", "red"))(50)
      pheatmap(
        analysis_results()$heatmap_data,
        cluster_rows   = TRUE,
        show_rownames  = TRUE,
        scale          = "row",
        color          = custom_colors,
        show_colnames  = FALSE,
        cluster_cols   = TRUE,
        annotation_col = analysis_results()$metadata_selected["selected_gene_expression"],
        main           = paste0("Within Pauci-immune RA\n", analysis_results()$selected_gene, " Differential Expression")
      )
      dev.off()
    }
  )
  
  # Download DE results as an Excel file.
  output$download_DE_results <- downloadHandler(
    filename = function() {
      paste0(analysis_results()$selected_gene, "_DE_results.xlsx")
    },
    content = function(file) {
      results_to_save <- as.data.frame(analysis_results()$res_vs)
      results_to_save$gene <- rownames(results_to_save)
      writexl::write_xlsx(results_to_save, file)
    }
  )
  output$download_nrm_counts <- downloadHandler(
    filename = function() {
      "nrm_counts.csv"
    },
    content = function(file) {
      write.csv(nrm_counts,"normalized_counts.csv")
    }
  )
  
}

shinyApp(ui = ui, server = server)
