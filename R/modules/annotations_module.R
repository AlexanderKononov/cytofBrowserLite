#' Annotations Module for cytofBrowserLite
#'
#' This module encapsulates all server logic related to annotation management
#' and visualization.
#'
#' @param input Shiny input object.
#' @param output Shiny output object.
#' @param session Shiny session object.
#' @param fcs_data Reactive values object containing main data.
#' @param clusters Reactive values object for clustering data.
#' @param data_prep_settings Reactive values object for data preparation settings.
#' @param plots Reactive values object for plots.
#' 
#' @return A list of reactive elements and observers for the annotations module.

annotations_module <- function(input, output, session, fcs_data, clusters, data_prep_settings, plots) {
  
  # Define reactive elements and observers specific to annotations
  # This will include all the logic from the 'Annotations' section of serverLite.R
  
  # Example of how to structure the module content:
  # (This is a simplified example, the actual content will be much larger)
  
  # UI to chose annotation for rename annotaions
  output$chose_managment_ann_ui <- renderUI({
    if(is.null(fcs_data$cell_ann)){return(NULL)}
    col_names <- colnames(fcs_data$cell_ann)
    col_names <- col_names[!(col_names %in% c("all_cells", "tSNE1", "tSNE2", "UMAP1", "UMAP2", "samples"))]
    if(length(col_names) < 1){return(NULL)}
    sel <- col_names[1]
    if("clusters" %in% col_names){sel <- "clusters"}
    selectInput('gp_managment_ann', label = h5("Annotation to manage"), choices = col_names, selected = sel)
  })

  # Dynamic UI for rename each group in chosen annotation
  output$managment_ann_ui <- renderUI({
    if(is.null(input$gp_managment_ann)){return(NULL)}
    lapply(unique(fcs_data$cell_ann[,input$gp_managment_ann]), function(i) {
      textInput(paste0('group', i), label = NULL, value = i)
    })
  })

  # change name in annotations
  observeEvent(input$rename_groups, {
    for (i in unique(fcs_data$cell_ann[,input$gp_managment_ann])){
      tmp <- fcs_data$cell_ann[,input$gp_managment_ann]
      tmp[tmp == i] <- input[[paste0('group', i)]]
      fcs_data$cell_ann[,input$gp_managment_ann] <- tmp
    }
  })



  # UI for advanced options annotatios page
  output$advanced_opt_ann_ui <- renderUI({
    if(is.null(input$method_plot_ann)){return(NULL)}
    if(input$method_plot_ann == 'tSNE'){
      ui <- fluidRow(
        column(1),
        column(10,
               numericInput("perplexity_ann", "tSNE Perplexity", min = 0, max = 200, value = 30, step = 5),
               numericInput("theta_ann", "tSNE Theta", min = 0, max = 1, value = 0.5, step = 0.1),
               numericInput("max_iter_ann", "tSNE Iterations", value = 1000, step = 500)
        )
      )
    }
    if(input$method_plot_ann == 'UMAP'){ui <- NULL}
    return(ui)
  })

  # Isolate number of plot from continiouse rection
  annotations <- reactiveValues(num_plot = 4)
  observeEvent(input$redraw_ann, {annotations$num_plot <- input$num_plot_ann})

  # UI to choose marker for scatter plots for annotating
  observe({
    lapply(1:annotations$num_plot, function(x){
      output[[paste0("mk_ann_ui", x)]] <- renderUI({
        selectInput(paste0("mk_target_ann", x), label = NULL,
                    choices = c(input$gp_managment_ann, names(fcs_data$use_markers)),selected = 1)
      })
    })
  })


  # Drawing scatter plots for annotating
  observe({
    if(!is.numeric(input$num_plot_ann)){return(NULL)}
    lapply(1:annotations$num_plot, function(x){
      output[[paste0("plots_ann", x)]] <- renderPlot({
        if(is.null(fcs_data$cell_ann$clusters)){return(NULL)}
        if(is.null(fcs_data$subset_coord)){return(NULL)}
        color_mk <- "clusters"
        if(!is.null(input[[paste0("mk_target_ann",x)]])){color_mk <- input[[paste0("mk_target_ann",x)]]}
        if(color_mk %in% names(fcs_data$use_markers)){color_mk <- fcs_data$use_markers[color_mk]}
        if(data_prep_settings$method == "tSNE"){
          plot_data <- data.frame(X = fcs_data$cell_ann[fcs_data$subset_coord, "tSNE1"],
                                  Y = fcs_data$cell_ann[fcs_data$subset_coord, "tSNE2"])}
        if(data_prep_settings$method == "UMAP"){
          plot_data <- data.frame(X = fcs_data$cell_ann[fcs_data$subset_coord, "UMAP1"],
                                  Y = fcs_data$cell_ann[fcs_data$subset_coord, "UMAP2"])}
        if(color_mk %in% colnames(fcs_data$cell_ann)){
          plot_data$mk <- as.factor(fcs_data$cell_ann[fcs_data$subset_coord, "clusters"])}
        if(!(color_mk %in% colnames(fcs_data$cell_ann))){
          plot_data$mk <- fcs_data$exprs_data[fcs_data$subset_coord, color_mk]}
        plt <- ggplot2::ggplot(plot_data,  aes(x = X, y = Y, color = mk)) +
          geom_point(size = input$point_size_clust)
        if(color_mk == "clusters"){plt <- plt + scale_color_manual(values = as.character(clusters$nodes$color))}
        if(color_mk != "clusters"){plt <- plt + scale_color_gradient2(midpoint=0.5, low='blue', mid='gray', high='red')}
        plt <- plt + theme_bw() + labs(color = color_mk)+
          guides(colour = guide_legend(override.aes = list(size=2)))
        return(plt)
      })
    })
  })


  output$plot_set_ann_ui <- renderUI({
    if(is.null(fcs_data$cell_ann)){return(NULL)}
    row_num <- (annotations$num_plot %/% 2)+(annotations$num_plot %% 2)
    ui <- fluidPage(
      column(6,
             lapply(1:row_num, function(x){
               fluidRow(
                 uiOutput(paste0("mk_ann_ui", x)),
                 plotOutput(paste0("plots_ann", x))
               )
             })
             ),
      column(6,
             lapply(((row_num + 1):annotations$num_plot), function(x){
               fluidRow(
                 uiOutput(paste0("mk_ann_ui", x)),
                 plotOutput(paste0("plots_ann", x))
               )
             })
             )
    )
    return(ui)
  })

  # Rewrite for scatter plot in Annotation
  observeEvent(input$redraw_ann, {
    withProgress(message = "Extraction data", min =0, max = 3, value = 0,{
      incProgress(1, detail = "Subsampling" )
      ## If subset was changed
      if(any(c(data_prep_settings$sampling_size, data_prep_settings$fuse) != c(input$sampling_size_ann, input$fuse_ann))){
        if(!is.null(input$sampling_size_clust)){data_prep_settings$sampling_size <- input$sampling_size_clust}
        if(!is.null(input$fuse_ann)){data_prep_settings$fuse <- input$fuse_ann}
        fcs_data$subset_coord <- get_subset_coord(cell_ann = fcs_data$cell_ann,
                                                  sampling_size = data_prep_settings$sampling_size,
                                                  fuse = data_prep_settings$fuse,
                                                  size_fuse = data_prep_settings$size_fuse)
      }
      ## Repeat dimension reducing
      incProgress(1, detail = "Dimention redicing" )
      dim_reduce_force <- FALSE
      if(!is.null(input$method_plot_ann)){data_prep_settings$method <- input$method_plot_ann; dim_reduce_force <- TRUE}
      if(!is.null(input$perplexity_ann)){data_prep_settings$perplexity <- input$perplexity_ann; dim_reduce_force <- TRUE}
      if(!is.null(input$theta_ann)){data_prep_settings$theta <- input$theta_ann; dim_reduce_force <- TRUE}
      if(!is.null(input$max_iter_ann)){data_prep_settings$max_iter <- input$max_iter_ann; dim_reduce_force <- TRUE}
      fcs_data$cell_ann <- get_dim_reduce(fcs_data$exprs_data, fcs_data$cell_ann,
                                          fcs_data$subset_coord, fcs_data$use_markers,force = dim_reduce_force,
                                          method = data_prep_settings$method, perplexity = data_prep_settings$perplexity,
                                          theta = data_prep_settings$theta, max_iter = data_prep_settings$max_iter,
                                          pca_param = FALSE, check_duplicates = FALSE, seed = 1234)
      incProgress(1, detail = "Plotting" )
    })
  })
  
  # Return a list of reactive elements if needed
  # For now, this module primarily sets up observers and outputs
  return(list())
}