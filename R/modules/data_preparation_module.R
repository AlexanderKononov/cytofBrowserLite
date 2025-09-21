#' Data Preparation Module for cytofBrowserLite
#'
#' This module encapsulates all server logic related to data preparation,
#' including file upload, preprocessing, transformations, and initial plotting.
#'
#' @param input Shiny input object.
#' @param output Shiny output object.
#' @param session Shiny session object.
#' @param fcs_data Reactive values object containing main data.
#' @param data_prep_settings Reactive values object for data preparation settings.
#' @param gates Reactive values object for gating data.
#' @param plots Reactive values object for plots.
#' @param roots File system roots for shinyFiles.
#' 
#' @return A list of reactive elements and observers for the data preparation module.

data_preparation_module <- function(input, output, session, fcs_data, data_prep_settings, gates, plots, roots) {
  
  # Define reactive elements and observers specific to data preparation
  # This will include all the logic from the 'Data preparation' section of serverLite.R
  
  # Example of how to structure the module content:
  # (This is a simplified example, the actual content will be much larger)
  
  # Showing selected fcs files
  output$selected_fcs_dp_ui <- renderUI({
    if((length(input$choose_fcs_dp) <= 1)){return(NULL)}
    fluidRow(
      column(1),
      column(10,
             h5("Selected files"),
             verbatimTextOutput('selected_fcs_dp')
      )
    )
  })
  output$selected_fcs_dp <- renderPrint({as.character(parseFilePaths(roots, input$choose_fcs_dp)$datapath)})

  # Upload data and automatic pre-processing steps
  observeEvent(input$butt_upload_dproc, {
    
    if((length(input$choose_fcs_dp) <= 1) & !input$test_data_upload_dproc){
      showNotification("Files were not chosen", type = "warning")
      return(NULL)}
    withProgress(message = "Extraction data", min =0, max = 12, value = 0,{
      ## Get row data fcs files
      if(input$test_data_upload_dproc){
        showNotification("Build-in data set were not uploaded", type = "warning")
        fcs_data$md <- get_test_fcs_metadata(input$test_data_dproc)}
      if(!input$test_data_upload_dproc){fcs_data$md <- get_fcs_metadata(parseFilePaths(roots, input$choose_fcs_dp)$datapath)}
      incProgress(1, detail = "Upload data" )
      fcs_data$fcs_raw <- get_fcs_raw(fcs_data$md)
      if(is.null(fcs_data$fcs_raw)){
        showNotification("Files have different channels", type = "error")
        return(NULL)}
      incProgress(1, detail = "Extraction ereachment" )
      fcs_data$exprs_data <- get_exprs_data(fcs_data$fcs_raw)
      incProgress(1, detail = "Extraction panel")
      fcs_data$panel <- get_fcs_panel(fcs_data$fcs_raw)
      incProgress(1, detail = "Markers processing" )
      fcs_data$use_markers <- get_use_marker(fcs_data$panel)
      fcs_data$entire_panel <- get_entire_panel(fcs_raw = fcs_data$fcs_raw)
      incProgress(1, detail = "Cell number calculation" )
      fcs_data$cell_number <- get_cell_number(fcs_data$fcs_raw)
      incProgress(1, detail = "Transformations" )
      fcs_data$trans <- get_extr_transformations(fcs_data$exprs_data, fcs_data$use_markers, mode = "dataset")
      fcs_data$trans_to_do <- check_transformation_list(fcs_data$trans_to_do, fcs_data$trans)
      fcs_data$exprs_data <- get_transformations(fcs_data$trans_to_do, fcs_data$exprs_data, fcs_data$use_markers)
      fcs_data$trans <- c(fcs_data$trans, fcs_data$trans_to_do)
      fcs_data$trans_to_do <- c()
      incProgress(1, detail = "Creating of cell annotation" )
      fcs_data$cell_ann <- data.frame(all_cells = rep("cell", sum(fcs_data$cell_number$cell_nmbr)),
                                samples = rep(fcs_data$cell_number$smpl, fcs_data$cell_number$cell_nmbr),
                                row.names = 1:sum(fcs_data$cell_number$cell_nmbr))
      if(input$extr_ann_dp){fcs_data$ann_panel <- get_extr_ann(exprs_data = fcs_data$exprs_data)}
      if(!is.null(fcs_data$ann_panel)){
        fcs_data$cell_ann <- cbind(fcs_data$cell_ann, fcs_data$exprs_data[, fcs_data$ann_panel])}

      incProgress(1, detail = "Subsampling" )
      fcs_data$subset_coord <- get_subset_coord(cell_ann = fcs_data$cell_ann,
                                                sampling_size = data_prep_settings$sampling_size,
                                                fuse = data_prep_settings$fuse,
                                                size_fuse = data_prep_settings$size_fuse)
      incProgress(1, detail = "Dimention redicing" )
      fcs_data$cell_ann <- get_dim_reduce(fcs_data$exprs_data, fcs_data$cell_ann,
                                          fcs_data$subset_coord, fcs_data$use_markers,
                                          method = data_prep_settings$method, perplexity = data_prep_settings$perplexity,
                                          theta = data_prep_settings$theta, max_iter = data_prep_settings$max_iter,
                                          pca_param = FALSE, check_duplicates = FALSE, seed = 1234)
      incProgress(1, detail = "Gate template" )
      ## Add primer column to cell type data frame and gates data frame
      gates$gates <- data.frame(all_cells = rep(TRUE, sum(fcs_data$cell_number$cell_nmbr)))
      gates$antology <- data.frame(name = "all_cells", parent = NA)
      rownames(gates$antology) <- gates$antology$name
      incProgress(1, detail = "Extraction fcs cluster info")
      incProgress(1)
    })
  })


  # Drawing the reactive scatter plot to data preparation
  output$scatter_plot_dp <- renderPlot({
    if(is.null(fcs_data$subset_coord)){return(NULL)}
    color_mk <- names(fcs_data$use_markers)[1]
    if(!is.null(input$mk_scatter_dp)){color_mk <- input$mk_scatter_dp}
    if(color_mk %in% names(fcs_data$use_markers)){color_mk <- fcs_data$use_markers[color_mk]}
    if(data_prep_settings$method == "tSNE"){
      plot_data <- data.frame(X = fcs_data$cell_ann[fcs_data$subset_coord, "tSNE1"],
                              Y = fcs_data$cell_ann[fcs_data$subset_coord, "tSNE2"])}
    if(data_prep_settings$method == "UMAP"){
      plot_data <- data.frame(X = fcs_data$cell_ann[fcs_data$subset_coord, "UMAP1"],
                              Y = fcs_data$cell_ann[fcs_data$subset_coord, "UMAP2"])}


    if(color_mk %in% colnames(fcs_data$cell_ann)){plot_data$mk <- as.factor(fcs_data$cell_ann[fcs_data$subset_coord, color_mk])}
    if(color_mk %in% colnames(fcs_data$exprs_data)){plot_data$mk <- fcs_data$exprs_data[fcs_data$subset_coord, color_mk]}

    plt <- ggplot2::ggplot(plot_data,  aes(x = X, y = Y, color = mk)) +
      geom_point(size = input$point_size)
    if(color_mk %in% colnames(fcs_data$cell_ann)){plt <- plt +
      #scale_color_manual(values = as.character(clusters$nodes$color))+
      guides(colour = guide_legend(override.aes = list(size=2)))}
    if(color_mk %in% colnames(fcs_data$exprs_data)){plt <- plt + scale_color_gradient2(midpoint=0.5, low='blue', mid='gray', high='red')}
    plt <- plt +
      labs(color = color_mk) +
      theme_bw()
    plots$scatter_dp <- plt
    return(plots$scatter_dp)
  })


  # Rewrite for scatter plot
  observeEvent(input$redraw_dp, {
    print(input$redraw_dp)
    print(str(input$redraw_dp))

    withProgress(message = "Extraction data", min =0, max = 3, value = 0,{
      incProgress(1, detail = "Subsampling" )
      print("--Repeat sabsampling--")
      ## If subset was changed
      #if(any(c(data_prep_settings$sampling_size, data_prep_settings$fuse) != c(input$sampling_size, input$fuse))){
      #}
      if(!is.null(input$sampling_size)){data_prep_settings$sampling_size <- input$sampling_size}
      if(!is.null(input$fuse)){data_prep_settings$fuse <- input$fuse}
      fcs_data$subset_coord <- get_subset_coord(cell_ann = fcs_data$cell_ann,
                                                sampling_size = data_prep_settings$sampling_size,
                                                fuse = data_prep_settings$fuse,
                                                size_fuse = data_prep_settings$size_fuse)
      ## Repeat dimension reducing
      print("--Repeat dim redicing--")
      incProgress(1, detail = "Dimention redicing" )
      dim_reduce_force <- FALSE
      if(!is.null(input$method_plot_dp)){data_prep_settings$method <- input$method_plot_dp; dim_reduce_force <- TRUE}
      if(!is.null(input$data_prep_perplexity)){data_prep_settings$perplexity <- input$data_prep_perplexity; dim_reduce_force <- TRUE}
      if(!is.null(input$data_prep_theta)){data_prep_settings$theta <- input$data_prep_theta; dim_reduce_force <- TRUE}
      if(!is.null(input$data_prep_max_iter)){data_prep_settings$max_iter <- input$data_prep_max_iter; dim_reduce_force <- TRUE}
      fcs_data$cell_ann <- get_dim_reduce(fcs_data$exprs_data, fcs_data$cell_ann,
                                          fcs_data$subset_coord, fcs_data$use_markers,force = dim_reduce_force,
                                          method = data_prep_settings$method, perplexity = data_prep_settings$perplexity,
                                          theta = data_prep_settings$theta, max_iter = data_prep_settings$max_iter,
                                          pca_param = FALSE, check_duplicates = FALSE, seed = 1234)
      incProgress(1, detail = "Plotting" )
    })
  })

  # Show the implemented transformations
  output$trans_dp <- renderPrint({fcs_data$trans})

  # Reaction to button "Transform" in "Data processing"
  observeEvent(input$butt_trans_dproc, {
    if(is.null(fcs_data$exprs_data)){return(NULL)}
    ## Transform row data to scaled data by set parameters
    withProgress(message = "Transformation", min =0, max = 3, value = 0,{
      incProgress(1)
      fcs_data$trans_to_do <- c(fcs_data$trans_to_do, input$transformation_list)
      fcs_data$trans_to_do <- check_transformation_list(fcs_data$trans_to_do, fcs_data$trans)
      incProgress(1)
      fcs_data$exprs_data <- get_transformations(fcs_data$trans_to_do, fcs_data$exprs_data, fcs_data$use_markers)
      fcs_data$trans <- c(fcs_data$trans, fcs_data$trans_to_do)
      incProgress(1)
    })
  })

  # Create UI to choose excluded markers
  output$mk_subset_dp_ui <- renderUI({
    if(is.null(fcs_data$use_markers)|is.null(fcs_data$entire_panel)){return(NULL)}
    fluidRow(
      column(1),
      column(10,
             selectInput("exclude_mk_dp", label = "Exclude markers",
                         choices = names(fcs_data$use_markers),
                         multiple = TRUE),
             selectInput("include_mk_dp", label = "Include markers",
                         choices = names(fcs_data$entire_panel[!(fcs_data$entire_panel %in% fcs_data$use_markers)]),
                         multiple = TRUE),
             actionButton("change_mk_button", label = "Change")
      )
    )
  })

  # Update reactive object use_markers after marker changing
  observeEvent(input$change_mk_button, {
    if(!is.null(input$exclude_mk_dp)){
      fcs_data$use_markers <- fcs_data$use_markers[!(names(fcs_data$use_markers) %in% input$exclude_mk_dp)]
    }
    if(!is.null(input$include_mk_dp)){
      fcs_data$use_markers <- fcs_data$entire_panel[c(names(fcs_data$use_markers), input$include_mk_dp)]
    }
  })

  # Create UI to choose excluded extracted annotations
  output$extr_ann_manage_dp_ui <- renderUI({
    if(is.null(fcs_data$ann_panel)){return(NULL)}
    fluidRow(
      column(1),
      column(10,
             selectInput("exclude_extr_ann_dp", label = "Exclude data from annotation",
                         choices = fcs_data$ann_panel,
                         multiple = TRUE),
             selectInput("include_extr_ann_dp", label = "Include data in annotation",
                         choices = fcs_data$entire_panel[!(fcs_data$entire_panel %in% fcs_data$ann_panel)],
                         multiple = TRUE),
             actionButton("change_extr_ann_button", label = "Change")
      )
    )
  })

  # Update annotation objects after extaract ann changing
  observeEvent(input$change_extr_ann_button, {
    if(!is.null(input$exclude_extr_ann_dp)){
      fcs_data$ann_panel <- fcs_data$ann_panel[!(fcs_data$ann_panel %in% input$exclude_extr_ann_dp)]
      fcs_data$cell_ann <- fcs_data$cell_ann[,!(colnames(fcs_data$cell_ann) %in% input$exclude_extr_ann_dp)]
    }
    if(!is.null(input$include_extr_ann_dp)){
      fcs_data$ann_panel <- c(fcs_data$ann_panel, input$include_extr_ann_dp)
      fcs_data$cell_ann <- cbind(fcs_data$cell_ann, fcs_data$exprs_data[, input$include_extr_ann_dp])
    }
  })


  # UI for saving cell annotation as fcs files
  output$save_cell_ann_dp_ui <- renderUI({
    if(is.null(fcs_data$cell_ann)){return(NULL)}
    selectInput('save_cell_ann_dp', label = h5("Adding cell annotation to FCS files"),
                choices = colnames(fcs_data$cell_ann), multiple = TRUE)
  })

  # First part of UI for subsetting data by annotation
  output$subset_dp_1ui <- renderUI({
    if(is.null(fcs_data$cell_ann)){return(NULL)}
    print(head(fcs_data$cell_ann))
    print(str(fcs_data$cell_ann))
    choices_set <- colnames(fcs_data$cell_ann)
    choices_set <- choices_set[!grepl("tSNE", choices_set)]
    choices_set <- choices_set[!grepl("UMAP", choices_set)]
    choices_set <- choices_set[!grepl("all_cells", choices_set)]
    if(length(choices_set) == 0){return(NULL)}
    fluidRow(
      column(1),
      column(10,
             radioGroupButtons(inputId = "subset_orient_dp", label = h5("Subsatting by annotaion"),
                               choices = c("exclude", "include"), selected = "exclude", justified = TRUE),
             selectInput("ann_to_subset", label = h5("Annotation to use "),
                         choices = choices_set,  multiple = FALSE)
      )
    )
  })

  # Second part of UI for subsetting data by annotation
  output$subset_dp_2ui <- renderUI({
    if(is.null(fcs_data$cell_ann)){return(NULL)}
    if(is.null(input$ann_to_subset)){return(NULL)}
    fluidRow(
      column(1),
      column(10,
             selectInput("ann_groups_to_subset", label = h5("Cell groups to use"),
                         choices = unique(fcs_data$cell_ann[,input$ann_to_subset]), multiple = TRUE),
             actionButton("butt_subset_dp", label = "Subset")
      )
    )
  })

  # Subsatting data by annotation
  observeEvent(input$butt_subset_dp, {
    if(is.null(fcs_data$cell_ann)){return(NULL)}
    if(is.null(input$ann_to_subset)){return(NULL)}
    if(is.null(input$ann_groups_to_subset)){return(NULL)}
    print("----Start subsetting-----")
    print(input$redraw_dp)
    print(str(input$redraw_dp))

    print("---1")
    stay_coord <- fcs_data$cell_ann[,input$ann_to_subset] %in% input$ann_groups_to_subset
    print("---2")
    if(input$subset_orient_dp == "exclude"){stay_coord <- !stay_coord}
    print("---3")
    fcs_data$exprs_data <- fcs_data$exprs_data[stay_coord,]
    print("---4")
    fcs_data$cell_ann <- fcs_data$cell_ann[stay_coord,]
    print("---5")
    fcs_data$subset_coord <- which(stay_coord[as.integer(fcs_data$subset_coord)])
    print("---6")
    if(!is.null(gates$gates)){
      new_gates <- as.data.frame(as.data.frame(gates$gates)[stay_coord,])
      colnames(new_gates) <- colnames(gates$gates)
      gates$gates <- as.data.frame(new_gates)
      }
    print("---7")
    if(!is.null(gates$subset_coord)){gates$subset_coord <- which(stay_coord[gates$subset_coord])}
    print("---7")
  })

  # Save data as fcs files
  shinyDirChoose(input, 'choose_dwn_folder_dp', roots = roots)
  observeEvent(input$dwn_fcs_dp, {
    if(is.null(fcs_data$fcs_raw) | is.null(input$choose_dwn_folder_dp)){return(NULL)}
    panel_folder_path <- parseDirPath(roots, input$choose_dwn_folder_dp)
    if(length(panel_folder_path) == 0){return(NULL)}
    if(!is.null(panel_folder_path) | length(panel_folder_path) !=0 ){
      filenames <- get_new_fcs_names(old_name = paste0(flowCore::sampleNames(fcs_data$fcs_raw), ".fcs"),
                                     folder_path = panel_folder_path, prefix = input$name_prefix_dp)
      new_fcs <- get_ann_adding_fcs_files(fcs_raw = fcs_data$fcs_raw, cell_ann = fcs_data$cell_ann,
                                         col_name = input$save_cell_ann_dp)
      flowCore::write.flowSet(new_fcs, outdir = as.character(panel_folder_path), filename = filenames)
    }
  })


  # UI for advanced options data preparation
  output$advanced_opt_dp_ui <- renderUI({
    if(is.null(input$method_plot_dp)){return(NULL)}
    if(is.null(input$method_plot_dp)){return(NULL)}
    if(input$method_plot_dp == 'tSNE'){
      ui <- fluidRow(
        column(1),
        column(10,
               numericInput("data_prep_perplexity", "tSNE Perplexity", min = 0, max = 200, value = 30, step = 5),
               numericInput("data_prep_theta", "tSNE Theta", min = 0, max = 1, value = 0.5, step = 0.1),
               numericInput("data_prep_max_iter", "tSNE Iterations", value = 1000, step = 500)
        )
      )
    }
    if(input$method_plot_dp == 'UMAP'){ui <- NULL}
    return(ui)
  })


  # UI to choose marker for scatter plot dp
  output$mk_scatter_dp_ui <- renderUI({
    if(is.null(fcs_data$use_markers)){return(NULL)}
    choices_set <- names(fcs_data$use_markers)
    if(!is.null(fcs_data$cell_ann)){
      choices_set <- c(colnames(fcs_data$cell_ann), choices_set)
      choices_set <- choices_set[!grepl("tSNE", choices_set)]
      choices_set <- choices_set[!grepl("UMAP", choices_set)]
      choices_set <- choices_set[!grepl("all_cells", choices_set)]
    }
    selected_1 <- names(fcs_data$use_markers)[1]
    selectInput('mk_scatter_dp', label = h4("Plotted marker"),
                choices = choices_set, selected = selected_1)
  })


  # Download scatter plot data preparation
  output$dwn_scatter_dp <- downloadHandler(
    filename = function() {
      ext <- input$dwn_scatter_dp_ext
      if(is.null(ext)){ext <- "pdf"}
      paste("Scatter_plot_data_preparation", ext, sep = ".") },
    content = function(file) {
      ext <- input$dwn_scatter_dp_ext
      if(is.null(ext)){ext <- "pdf"}
      ggplot2::ggsave(file, plot = plots$scatter_dp, device = ext, width = input$dwn_scatter_dp_width,
                      height = input$dwn_scatter_dp_height, units = input$dwn_scatter_dp_units,
                      dpi =input$dwn_scatter_dp_dpi)
    }
  )

  # Drawing the reactive plot of cell number
  output$smpl_hist_preparation <- renderPlot({
    if(is.null(fcs_data$cell_number)){return(NULL)}
    plots$smpl_hist <- ggplot(data = fcs_data$cell_number, aes(x = smpl, y = cell_nmbr))+
      geom_bar(stat="identity", fill = "black")+
      xlab("Samples")+
      ylab("Number of cells")+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    return(plots$smpl_hist)
  })

  # Download plot of cell number data preparation
  output$dwn_smpl_hist_dp <- downloadHandler(
    filename = function() {
      ext <- input$dwn_smpl_hist_dp_ext
      if(is.null(ext)){ext <- "pdf"}
      paste("Cell_number_plot_data_preparation", ext, sep = ".") },
    content = function(file) {
      ext <- input$dwn_smpl_hist_dp_ext
      if(is.null(ext)){ext <- "pdf"}
      ggsave(file, plot = plots$smpl_hist, device = ext, width = input$dwn_smpl_hist_dp_width,
             height = input$dwn_smpl_hist_dp_height, units = input$dwn_smpl_hist_dp_units,
             dpi =input$dwn_smpl_hist_dp_dpi)
    }
  )

  # Create UI to choose target marker for marker density plot
  output$mk_density_dp_ui <- renderUI({
    if(is.null(fcs_data$use_markers)){ return(NULL)}
    selectInput('mk_density_dp', label = h5("Plotted marker"),
                choices = names(fcs_data$use_markers), selected = 1)
  })

  # Drawing the reactive histogram plot of marker expression
  output$mk_density_plot_dp <- renderPlot({
    if(is.null(fcs_data$use_markers)){return(NULL)}
    color_mk <- names(fcs_data$use_markers)[1]
    if(!is.null(input$mk_density_dp)){color_mk <- input$mk_density_dp}
    plot_data <- data.frame(expr_data = fcs_data$exprs_data[,fcs_data$use_markers[color_mk]])
    plot_data[,"expr_data"] <- as.numeric(plot_data[,"expr_data"])
    if(input$mk_zero_del_density_plot_dp){
      plot_data <- data.frame(expr_data = plot_data$expr_data[plot_data$expr_data != 0])
    }
    plots$mk_density <- ggplot(plot_data, aes(x = expr_data))+
      geom_density(fill = 'black')+
      labs(x = color_mk)+
      theme_bw()
    return(plots$mk_density)
  })

  # Download plot of marker histogram data preparation
  output$dwn_mk_density_dp <- downloadHandler(
    filename = function() {
      ext <- input$dwn_mk_density_dp_ext
      if(is.null(ext)){ext <- "pdf"}
      print(paste("Marker_density_plot_data_preparation", ext, sep = "."))
      return(paste("Marker_density_plot_data_preparation", ext, sep = ".")) },
    content = function(file) {
      ext <- input$dwn_mk_density_dp_ext
      if(is.null(ext)){ext <- "pdf"}
      ggsave(file, plot = plots$mk_density, device = ext, width = input$dwn_mk_density_dp_width,
             height = input$dwn_mk_density_dp_height, units = input$dwn_mk_density_dp_units,
             dpi =input$dwn_mk_density_dp_dpi)
    }
  )

  # Reactive show current use_markers from "fcs_data" object
  output$mk_rested_dp <- renderPrint({fcs_data$use_markers})
  output$mk_excluded_dp <- renderPrint({
    fcs_data$entire_panel[!(fcs_data$entire_panel %in% fcs_data$use_markers)]
  })
  
  # Return a list of reactive elements if needed
  # For now, this module primarily sets up observers and outputs
  return(list())
}