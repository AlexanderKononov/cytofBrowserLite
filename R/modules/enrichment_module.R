#' Enrichment Module for cytofBrowserLite
#'
#' This module encapsulates all server logic related to enrichment analysis,
#' including heatmap generation and deconvolution plots.
#'
#' @param input Shiny input object.
#' @param output Shiny output object.
#' @param session Shiny session object.
#' @param fcs_data Reactive values object containing main data.
#' @param plots Reactive values object for plots.
#' 
#' @return A list of reactive elements and observers for the enrichment module.

enrichment_module <- function(input, output, session, fcs_data, plots) {
  
  # Define reactive elements and observers specific to enrichment
  # This will include all the logic from the 'Enrichment' section of serverLite.R
  
  # Example of how to structure the module content:
  # (This is a simplified example, the actual content will be much larger)
  
  # UI to chose annotation for heatmap
  output$hm_rows_enrich_ui <- renderUI({
    if(is.null(fcs_data$cell_ann)){return(NULL)}
    col_names <- colnames(fcs_data$cell_ann)
    col_names <- col_names[!(col_names %in% c("all_cells", "tSNE1", "tSNE2", "UMAP1", "UMAP2"))]
    if(length(col_names) < 1){return(NULL)}
    sel <- col_names[1]
    if("clusters" %in% col_names){sel <- "clusters"}
    selectInput('hm_rows_enrich', label = h5("Heatmap rows"), choices = col_names, selected = sel)
  })

  # Drawing heatmap for enrichments
  output$hm_enrich <- renderPlot({
    if(is.null(input$hm_rows_enrich)){return(NULL)}
    hm_enrich_data <- get_hm_enrichments(exprs_data = fcs_data$exprs_data, use_markers = fcs_data$use_markers,
                                         groups = fcs_data$cell_ann[,input$hm_rows_enrich],
                                         summarize_mode = input$method_summarize_enrich)
    if(nrow(hm_enrich_data) < 2){
      showNotification("It was chosen less than 2 groups for heatmap", type = "warning")
      return(NULL)}
    hm_enrich_data <- get_scale_hm_enrichments(hm_enrich_data, method = input$row_scaling_enrich)
    color <- get_color_enrichments(pal = input$hm_enrich_palette, pal_reverse = input$hm_enrich_rev_pal)
    cluster_columns <- input$hm_enrich_col_dend
    if(any(is.na(hm_enrich_data))){
      cluster_columns <- FALSE
      showNotification("Some columns can't be scaled. Columns were not clustered.", type = "warning")
    }
    plots$hm_enrich <- ComplexHeatmap::Heatmap(hm_enrich_data, col = color,
                                               cluster_columns = cluster_columns,
                                               cluster_rows = input$hm_enrich_row_dend,
                                               show_heatmap_legend = input$hm_enrich_show_leg,
                                               heatmap_legend_param = list(title = "Enrichment"))
    return(plots$hm_enrich)
  })

  # Download clusters overlap plot in cluster section
  output$dwn_hm_enrich <- downloadHandler(
    filename = function() {
      ext <- input$dwn_hm_enrich_ext
      if(is.null(ext)){ext <- "pdf"}
      paste("Heatap_enrichments_plot", ext, sep = ".") },
    content = function(file) {
      ext <- input$dwn_hm_enrich_ext
      if(is.null(ext)){ext <- "pdf"}
      ggsave(file, plot = plots$hm_enrich, device = ext,
             width = input$dwn_hm_enrich_width,height = input$dwn_hm_enrich_height,
             units = input$dwn_hm_enrich_units,dpi =input$dwn_hm_enrich_dpi)
    }
  )

  # Create UI to choose marker for deconvolution plot
  output$mk_deconvol_enrich_ui <- renderUI({
    if(is.null(fcs_data$use_markers)){return(NULL)}
    selectInput("mk_deconvol_enrich", label = h5("Marker for deconvolution"),
                choices = names(fcs_data$use_markers),
                selected = 1)
  })
  # UI to chose annotation for deconvolution plot
  output$gp_deconvol_enrich_ui <- renderUI({
    if(is.null(fcs_data$cell_ann)){return(NULL)}
    col_names <- colnames(fcs_data$cell_ann)
    col_names <- col_names[!(col_names %in% c("all_cells", "tSNE1", "tSNE2", "UMAP1", "UMAP2", "samples"))]
    if(length(col_names) < 1){return(NULL)}
    sel <- col_names[1]
    if("clusters" %in% col_names){sel <- "clusters"}
    selectInput('gp_deconvol_enrich', label = h5("Annotation to plot"), choices = col_names, selected = sel)
  })


  # Drawing deconvolution of marker to sample-cluster plot
  output$deconvol_enrich <- renderPlot({
    if(is.null(input$gp_deconvol_enrich)){return(NULL)}
    if(is.null(input$mk_deconvol_enrich)){return(NULL)}
    if(is.null(fcs_data$exprs_data)){return(NULL)}
    deconvol_enrich_data <- get_deconvol_enrichments(exprs_data = fcs_data$exprs_data,
                                                     mk = fcs_data$use_markers[input$mk_deconvol_enrich],
                                                     samples = fcs_data$cell_ann$samples,
                                                     groups = fcs_data$cell_ann[, input$gp_deconvol_enrich],
                                                     summarize_mode = input$method_sum_deconvol_enrich)
    plots$deconvol_enrich <- ggplot(data = deconvol_enrich_data, aes(x=samples, y=groups, fill=enrich)) +
      geom_tile() +
      theme_bw()+
      ylab(input$gp_deconvol_enrich)+
      xlab("samples")+
      labs(fill = input$mk_deconvol_enrich) +
      geom_text(aes(label = round(cell_rate, 3)), size = 5, color = 'white') +
      theme(axis.text.x = element_text(angle = 90))
    return(plots$deconvol_enrich)
  })

  # Download clusters overlap plot in cluster section
  output$dwn_deconvol_enrich <- downloadHandler(
    filename = function() {
      ext <- input$dwn_deconvol_enrich_ext
      if(is.null(ext)){ext <- "pdf"}
      paste("Deconvolution_enrichments_plot", ext, sep = ".") },
    content = function(file) {
      ext <- input$dwn_deconvol_enrich_ext
      if(is.null(ext)){ext <- "pdf"}
      ggsave(file, plot = plots$deconvol_enrich, device = ext,
             width = input$dwn_deconvol_enrich_width,height = input$dwn_deconvol_enrich_height,
             units = input$dwn_deconvol_enrich_units,dpi =input$dwn_deconvol_enrich_dpi)
    }
  )
  
  # Return a list of reactive elements if needed
  # For now, this module primarily sets up observers and outputs
  return(list())
}