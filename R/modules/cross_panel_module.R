#' Cross-panel Module for cytofBrowserLite
#'
#' This module encapsulates all server logic related to cross-panel analysis,
#' including panel management and correlation analysis.
#'
#' @param input Shiny input object.
#' @param output Shiny output object.
#' @param session Shiny session object.
#' @param fcs_data Reactive values object containing main data.
#' @param plots Reactive values object for plots.
#' @param roots File system roots for shinyFiles.
#' 
#' @return A list of reactive elements and observers for the cross-panel module.

cross_panel_module <- function(input, output, session, fcs_data, plots, roots) {
  
  # Define reactive elements and observers specific to cross-panel analysis
  # This will include all the logic from the 'Cross-panel' section of serverLite.R
  
  # Example of how to structure the module content:
  # (This is a simplified example, the actual content will be much larger)
  
  # Creating reactive Values with panel data
  crosspanel <- reactiveValues(md_set = list(), fcs_raw_set = list(), panel_set = list(), use_marker_set = list(),
                               cell_clustering_list_set = list())

  shinyFileChoose(input, 'choose_fcs_cross_p', roots=roots)
  observeEvent(input$butt_upload_cross_p, {
    withProgress(message = "Clusters from fcs files", min =0, max = 4, value = 0,{
      panel_name <- input$panel_name_cross_p
      if(is.null(panel_name) | panel_name == ""){
        f_name <- basename(parseFilePaths(roots, input$choose_fcs_cross_p)$datapath[1])
        panel_name <- gsub(".fcs", "", f_name)
        if(grepl("CyBr_", f_name)){panel_name <- strsplit(f_name, split = "_")[[1]][2]}
      }
      crosspanel$md_set[[panel_name]] <- get_fcs_metadata(parseFilePaths(roots, input$choose_fcs_cross_p)$datapath)
      incProgress(1, detail = "Upload data" )
      crosspanel$fcs_raw_set[[panel_name]] <- get_fcs_raw(crosspanel$md_set[[panel_name]])
      incProgress(1, detail = "Extraction panel")
      crosspanel$panel_set[[panel_name]] <- get_fcs_panel(crosspanel$fcs_raw_set[[panel_name]])
      incProgress(1, detail = "Delete background markers" )
      crosspanel$use_marker_set[[panel_name]] <- get_use_marker(crosspanel$panel_set[[panel_name]])
      incProgress(1, detail = "Extract cluster info" )
      pattern <- "clust"
      if(!is.null(input$extr_clust_pattern_cross_p)){pattern <- input$extr_clust_pattern_cross_p}
      crosspanel$cell_clustering_list_set[[panel_name]] <- get_fcs_cluster_annotation(crosspanel$fcs_raw_set[[panel_name]], pattern = pattern)
    })
  })

  # Show the overlapped samples of the panels
  observe({
    if(length(crosspanel$md_set) == 0){return(NULL)}
    panel_content <- get_panel_content(crosspanel$fcs_raw_set)
    crosspanel$use_samples <- rownames(panel_content)[apply(panel_content, 1, function(x) {all(x == "presented")})]
  })

  # Show the overlapped markers of the panels
  observe({
    if(length(crosspanel$use_marker_set) == 0){return(NULL)}
    mk_panel_content <- get_common_markers(crosspanel$use_marker_set)
    crosspanel$use_common_mk <- rownames(mk_panel_content)[apply(mk_panel_content, 1, function(x) {all(x == "presented")})]
  })

  # Create ui for removing panels
  output$remove_panel_cross_p_ui <- renderUI({
    if(length(crosspanel$fcs_raw_set) < 1){return(NULL)}
    fluidRow(
      column(1),
      column(10,
             h4("Remove panels"),
             selectInput('panels_to_remove_cross_p', label = h5("Choose panels"),
                         choices = names(crosspanel$fcs_raw_set), multiple = TRUE),
             actionButton('remove_cross_p', label = "Remove")
      )
    )
  })

  # Removing panels
  observeEvent(input$remove_cross_p, {
    if(length(crosspanel$fcs_raw_set) < 1){return(NULL)}
    if(is.null(input$panels_to_remove_cross_p)){return(NULL)}
    stay_panel <- !(names(crosspanel$fcs_raw_set) %in% input$panels_to_remove_cross_p)
    crosspanel$md_set <- crosspanel$md_set[stay_panel]
    crosspanel$fcs_raw_set <- crosspanel$fcs_raw_set[stay_panel]
    crosspanel$panel_set <- crosspanel$panel_set[stay_panel]
    crosspanel$use_marker_set <- crosspanel$use_marker_set[stay_panel]
    crosspanel$cell_clustering_list_set <- crosspanel$cell_clustering_list_set[stay_panel]
  })

  # Show the sample content of the panels
  output$content_cross_p <- renderPlot({
    if(length(crosspanel$md_set) < 1){return(NULL)}
    panel_content <- get_panel_content(crosspanel$fcs_raw_set)
    ComplexHeatmap::Heatmap(panel_content, cluster_rows = F, cluster_columns = F,
                            row_names_side = "left", column_names_side = "top",
                            col = list("absent" = 'red3', "presented" = 'green3'), rect_gp = grid::gpar(col = "white", lwd = 2),
                            show_heatmap_legend = F)

  })

  # Show the marker content of the panels
  output$mk_content_cross_p <- renderPlot({
    if(length(crosspanel$use_marker_set) < 1){return(NULL)}
    mk_panel_content <- get_common_markers(crosspanel$use_marker_set)
    ComplexHeatmap::Heatmap(mk_panel_content, cluster_rows = F, cluster_columns = F,
                            row_names_side = "left", column_names_side = "top", row_names_gp = grid::gpar(fontsize = 8),
                            col = list("absent" = 'red3', "presented" = 'green3'), rect_gp = grid::gpar(col = "white", lwd = 2),
                            show_heatmap_legend = F)

  })

  # UI for choosing p-valueadjusting method for correlation cross-panel analysis
  output$abund_corr_padj_method_cross_p_ui <- renderUI({
    if(input$abund_corr_p_mode_cross_p == 'padj'){return((
      selectInput('abund_corr_padj_method_cross_p', "Select method for p-value adjusting",
                  choices = c("BH" = 'BH', "BY" = 'BY', "holm" ='holm', "hochberg" = 'hochberg', "hommel" = 'hommel',
                              "bonferroni" = 'bonferroni'), selected = 'BH')
    ))}
  })

  # Make abundance cross panel correlation analysis
  observeEvent(input$abund_corr_cross_p, {
    crosspanel$abund_cross_p <- get_abund_cross_panel(cell_clustering_list_set = crosspanel$cell_clustering_list_set,
                                                      use_samples = crosspanel$use_samples)
    crosspanel$abund_corr_info <- get_matrices_corr_info(crosspanel$abund_cross_p, method_cortest = input$abund_corr_method_cross_p,
                                                         method_padj = input$abund_corr_padj_method_cross_p)
  })

  # Drawing abundance correlation plot for cross panel analysis
  output$plot_abund_corr_cross_p <- renderPlot({
    if(is.null(crosspanel$abund_corr_info))return(NULL)
    signif_matrix <- crosspanel$abund_corr_info$padj
    if(input$abund_corr_p_mode_cross_p == 'pval'){signif_matrix <- crosspanel$abund_corr_info$p_value}
    corrplot::corrplot(crosspanel$abund_corr_info$corr_coef, method = "square", outline = T, order="hclust", type = 'lower',
                       p.mat = signif_matrix,
                       insig = 'label_sig', sig.level = c(.001, .01, .05), pch.cex = 0.7, pch.col = 'white',
                       addgrid.col = "darkgray", tl.col = "black", tl.cex = 1, cl.cex = 1,
                       col = colorRampPalette(c("midnightblue", "white", "darkred"))(100))
    gridGraphics::grid.echo()
    plots$abun_corr_cross_p <- grid::grid.grab()
    return(plots$abun_corr_cross_p)
  })

  # Download abundance correlation plot for cross panel analysis
  output$dwn_abund_corr_cross_p <- downloadHandler(
    filename = function() {
      ext <- input$dwn_abund_corr_cross_p_ext
      if(is.null(ext)){ext <- "pdf"}
      paste("Abunadnce_correlations_cross_panel", ext, sep = ".") },
    content = function(file) {
      ext <- input$dwn_abund_corr_cross_p_ext
      if(is.null(ext)){ext <- "pdf"}
      if(ext == "pdf"){
        pdf(file)
        grid::grid.draw(plots$abun_corr_cross_p)
        dev.off()
      }
      if(ext == "jpeg"){
        jpeg(file)
        grid::grid.draw(plots$abun_corr_cross_p)
        dev.off()
      }
      if(ext == "png"){
        png(file)
        grid::grid.draw(plots$abun_corr_cross_p)
        dev.off()
      }
      #ggsave(file, plot = plots$abun_corr_cross_p, device = ext)
    }
  )

  # Download abundance corelation tables  for cross-panel analusis
  output$dwn_table_abund_corr_cross_p <- downloadHandler(
    filename = function() {paste0("cross_panel_", input$dwn_table_abund_corr_cross_p_ext, ".csv")},
    content = function(file) {
      mode <- input$dwn_table_abund_corr_cross_p_ext
      if(is.null(mode)){mode <- 'corr_data'}
      if(mode == 'abund_data'){write.csv(crosspanel$abund_cross_p, file)}
      if(mode == 'corr_data'){write.csv(get_write_corr_info_cross_p(crosspanel$abund_corr_info), file)}
    }
  )

  # UI to choose marker for expressing cell fraction correlations
  output$mk_exp_cell_f_cross_p_ui <- renderUI({
    if(is.null(crosspanel$use_common_mk)){return(NULL)}
    selectInput('mk_exp_cell_f_cross_p', label = h5("Marker for cross-panel analysis"),
                choices = crosspanel$use_common_mk,
                selected = 1)
  })

  # UI to choose the threshold to expressing cell fraction allocation
  output$threshold_exp_cell_f_cross_p_ui <- renderUI({
    if(input$exp_cell_f_corr_method != 'threshold'){return(NULL)}
    if(is.null(crosspanel$fcs_raw_set)){return(NULL)}
    if(is.null(input$mk_exp_cell_f_cross_p)){return(NULL)}
    mk_exp_data_cross_p <- get_mk_exp_data_cross_p(fcs_raw_set = crosspanel$fcs_raw_set , use_samples = crosspanel$use_samples,
                                                   use_marker_set = crosspanel$use_marker_set, target_mk = input$mk_exp_cell_f_cross_p)
    mk_exp_data_cross_p <<- mk_exp_data_cross_p
    numericInput('threshold_exp_cell_f_cross_p', "threshold",min = min(mk_exp_data_cross_p), max = max(mk_exp_data_cross_p),
                 value = stats::median(mk_exp_data_cross_p), step = 0.1)
  })

  # UI for choosing p-value adjusting method for exp cell fraction cross-panel analysis
  output$exp_cell_f_corr_padj_method_cross_p_ui <- renderUI({
    if(input$exp_cell_f_corr_p_mode_cross_p == 'padj'){return((
      selectInput('exp_cell_f_corr_padj_method_cross_p', "Select method for p-value adjusting",
                  choices = c("BH" = 'BH', "BY" = 'BY', "holm" ='holm', "hochberg" = 'hochberg', "hommel" = 'hommel',
                              "bonferroni" = 'bonferroni'), selected = 'BH')
    ))}
  })

  # Make expressing cell fraction cross panel correlation analysis
  observeEvent(input$exp_cell_f_corr_cross_p, {
    exp_cell_f_corr_method <- input$exp_cell_f_corr_method
    if(is.null(exp_cell_f_corr_method)){exp_cell_f_corr_method <- 'clustering'}
    if(exp_cell_f_corr_method == 'clustering'){
      crosspanel$exp_cell_f_cross_p <- get_clustering_exp_cell_f_cross_panel(fcs_raw_set = crosspanel$fcs_raw_set,
                                                                             cell_clustering_list_set = crosspanel$cell_clustering_list_set,
                                                                             use_samples = crosspanel$use_samples,
                                                                             use_marker_set = crosspanel$use_marker_set,
                                                                             target_mk = input$mk_exp_cell_f_cross_p,
                                                                             min_cell_number = input$exp_cell_f_corr_min_cell_cross_p)
    }
    if(exp_cell_f_corr_method == 'threshold'){
      crosspanel$exp_cell_f_cross_p <- get_threshold_exp_cell_f_cross_panel(fcs_raw_set = crosspanel$fcs_raw_set,
                                                                            cell_clustering_list_set = crosspanel$cell_clustering_list_set,
                                                                            use_samples = crosspanel$use_samples,
                                                                            use_marker_set = crosspanel$use_marker_set,
                                                                            target_mk = input$mk_exp_cell_f_cross_p,
                                                                            threshold = input$threshold_exp_cell_f_cross_p,
                                                                            min_cell_number = input$exp_cell_f_corr_min_cell_cross_p)
    }
    crosspanel$exp_cell_f_corr_info <- get_matrices_corr_info(crosspanel$exp_cell_f_cross_p,
                                                              method_cortest = input$exp_cell_f_corr_method_cross_p,
                                                              method_padj = input$exp_cell_f_corr_padj_method_cross_p)
  })

  # Drawing abundance correlation plot for cross panel analysis
  output$plot_exp_cell_f_corr_cross_p <- renderPlot({
    if(is.null(crosspanel$exp_cell_f_corr_info))return(NULL)
    signif_matrix <- crosspanel$exp_cell_f_corr_info$padj
    if(input$exp_cell_f_corr_p_mode_cross_p == 'pval'){signif_matrix <- crosspanel$exp_cell_f_corr_info$p_value}
    corrplot::corrplot(crosspanel$exp_cell_f_corr_info$corr_coef, method = "square", outline = T, order="hclust", type = 'lower',
                       p.mat = signif_matrix,
                       insig = 'label_sig', sig.level = c(.001, .01, .05), pch.cex = 0.7, pch.col = 'white',
                       addgrid.col = "darkgray", tl.col = "black", tl.cex = 1, cl.cex = 1,
                       col = colorRampPalette(c("midnightblue", "white", "darkred"))(100))
    gridGraphics::grid.echo()
    plots$exp_cell_f_corr_cross_p <- grid::grid.grab()
    return(plots$exp_cell_f_corr_cross_p)
  })

  # Drawing expression density plot for a marker in cross-panel analysis
  output$mk_density_plot_cross_p <- renderPlot({
    if(is.null(crosspanel$fcs_raw_set)){return(NULL)}
    if(is.null(input$mk_exp_cell_f_cross_p)){return(NULL)}
    mk_exp_data_cross_p <- get_mk_exp_data_cross_p(fcs_raw_set = crosspanel$fcs_raw_set , use_samples = crosspanel$use_samples,
                                                   use_marker_set = crosspanel$use_marker_set, target_mk = input$mk_exp_cell_f_cross_p)
    ggplot(data.frame(expr = mk_exp_data_cross_p), aes(x = expr)) +
      geom_density(fill = 'black')
  })


  # Download abundance correlation plot for cross panel analysis
  output$dwn_exp_cell_f_corr_cross_p <- downloadHandler(
    filename = function() {
      ext <- input$dwn_exp_cell_f_corr_cross_p_ext
      if(is.null(ext)){ext <- "pdf"}
      paste0("Expr_cell_f_correlations_", input$mk_exp_cell_f_cross_p, "_cross_panel.", ext) },
    content = function(file) {
      ext <- input$dwn_exp_cell_f_corr_cross_p_ext
      if(is.null(ext)){ext <- "pdf"}
      if(ext == "pdf"){
        pdf(file)
        grid::grid.draw(plots$exp_cell_f_corr_cross_p)
        dev.off()
      }
      if(ext == "jpeg"){
        jpeg(file)
        grid::grid.draw(plots$exp_cell_f_corr_cross_p)
        dev.off()
      }
      if(ext == "png"){
        png(file)
        grid::grid.draw(plots$exp_cell_f_corr_cross_p)
        dev.off()
      }
      #ggsave(file, plot = plots$exp_cell_f_corr_cross_p, device = ext)
    }
  )

  # Download expressing cell fraction corelation tables for cross-panel analusis
  output$dwn_table_exp_cell_f_corr_cross_p <- downloadHandler(
    filename = function() {paste0("cross_panel_",input$mk_exp_cell_f_cross_p,"_", input$dwn_table_exp_cell_f_corr_cross_p_ext, ".csv")},
    content = function(file) {
      mode <- input$dwn_table_exp_cell_f_corr_cross_p_ext
      if(is.null(mode)){mode <- 'corr_data'}
      if(mode == 'exp_cell_f_data'){write.csv(crosspanel$exp_cell_f_cross_p, file)}
      if(mode == 'corr_data'){write.csv(get_write_corr_info_cross_p(crosspanel$exp_cell_f_corr_info), file)}
    }
  )
  
  # Return a list of reactive elements if needed
  # For now, this module primarily sets up observers and outputs
  # Also return the crosspanel reactive values object so it can be used in the main server function
  return(list(crosspanel = crosspanel))
}