#' Creation of server part of Shiny App
#' @description the function forms the body of backend part of Shiny App.
#' It also shapes the architecture of the pipeline. there are pipeline
#' blocks here which reflect the steps of CyTOF analysis. All computational
#' functions are placed in other files. The blocks of this function request
#' computational functions from that files by suitable order.
#' @param input shiny input object with from shiny UI to server
#' @param output shiny output object from server to UI and back to server
#'
#' @return
#'
#' @import shiny shinyFiles ggplot2
#' @importFrom flowCore sampleNames write.flowSet
#' @importFrom grDevices colorRampPalette
#'
cytofBrowser_server <- function(input, output){

  ########################
  ### Data preparation ###
  ########################
  roots <- c(Home = path.expand("~"), "R Installation" = R.home(), shinyFiles::getVolumes()())
  shinyFiles::shinyFileChoose(input, 'choose_fcs_dp', roots=roots, filetypes=c('', 'fcs'))

  ##### Create "fcs_data" as reactive object to store the CyTOF data
  fcs_data <-reactiveValues()
  gates <- reactiveValues()
  data_prep_settings <- reactiveValues(sampling_size = 0.5, fuse = TRUE, size_fuse = 3000, method = "tSNE",
                                       perplexity = 30, theta = 0.5, max_iter = 1000)
  plots <- reactiveValues()


  ##### Upload data and automatic pre-processing steps
  observeEvent(input$butt_upload_dproc, {
    if((length(input$choose_fcs_dp) <= 1)){return(NULL)}
    withProgress(message = "Extraction data", min =0, max = 11, value = 0,{
      ## Get row data fcs files
      fcs_data$md <- get_fcs_metadata(parseFilePaths(roots, input$choose_fcs_dp)$datapath)
      incProgress(1, detail = "Upload data" )
      fcs_data$fcs_raw <- get_fcs_raw(fcs_data$md)
      incProgress(1, detail = "Extraction ereachment" )
      fcs_data$exprs_data <- get_exprs_data(fcs_data$fcs_raw)
      incProgress(1, detail = "Extraction panel")
      fcs_data$panel <- get_fcs_panel(fcs_data$fcs_raw)
      incProgress(1, detail = "Markers processing" )
      fcs_data$use_markers <- get_use_marker(fcs_data$panel)
      fcs_data$entire_panel <- get_entire_panel(fcs_raw = fcs_data$fcs_raw)
      incProgress(1, detail = "Cell number calculation" )
      fcs_data$cell_number <- get_cell_number(fcs_data$fcs_raw)
      incProgress(1, detail = "Creating of cell annotation" )
      fcs_data$cell_ann <- data.frame(all_cells = rep("cell", sum(fcs_data$cell_number$cell_nmbr)),
                                samples = rep(fcs_data$cell_number$smpl, fcs_data$cell_number$cell_nmbr),
                                row.names = 1:sum(fcs_data$cell_number$cell_nmbr))
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


  ##### Drawing the reactive scatter plot to data preparation
  output$scatter_plot_dp <- renderPlot({
    if(is.null(fcs_data$subset_coord)){return(NULL)}
    color_mk <- names(fcs_data$use_markers)[1]
    if(!is.null(input$mk_scatter_dp)){color_mk <- input$mk_scatter_dp}
    if(data_prep_settings$method == "tSNE"){
      plot_data <- data.frame(X = fcs_data$cell_ann[fcs_data$subset_coord, "tSNE1"],
                              Y = fcs_data$cell_ann[fcs_data$subset_coord, "tSNE2"])}
    if(data_prep_settings$method == "UMAP"){
      plot_data <- data.frame(X = fcs_data$cell_ann[fcs_data$subset_coord, "UMAP1"],
                              Y = fcs_data$cell_ann[fcs_data$subset_coord, "UMAP2"])}
    plot_data$mk <- fcs_data$exprs_data[fcs_data$subset_coord, fcs_data$use_markers[color_mk]]
    plots$scatter_dp <- ggplot2::ggplot(plot_data,  aes(x = X, y = Y, color = mk)) +
      geom_point(size = input$point_size) +
      scale_color_gradient2(midpoint = 0.5, low = 'blue', mid = "gray",  high = 'red') +
      labs(color = color_mk) +
      theme_bw()
    return(plots$scatter_dp)
  })


  ##### Rewrite for scatter plot
  observeEvent(input$redraw_dp, {
    withProgress(message = "Extraction data", min =0, max = 3, value = 0,{
      incProgress(1, detail = "Subsampling" )
      ## If subset was changed
      if(any(c(data_prep_settings$sampling_size, data_prep_settings$fuse) != c(input$sampling_size, input$fuse))){
        if(!is.null(input$sampling_size)){data_prep_settings$sampling_size <- input$sampling_size}
        if(!is.null(input$fuse)){data_prep_settings$fuse <- input$fuse}
        fcs_data$subset_coord <- get_subset_coord(cell_ann = fcs_data$cell_ann,
                                                  sampling_size = data_prep_settings$sampling_size,
                                                  fuse = data_prep_settings$fuse,
                                                  size_fuse = data_prep_settings$size_fuse)
      }
      ## Repeat dimension reducing
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


  #### Reaction to button "Transform" in "Data processing"
  observeEvent(input$butt_trans_dproc, {
    if(is.null(fcs_data$exprs_data)){return(NULL)}
    ## Transform row data to scaled data by set parameters
    withProgress(message = "Transformation", min =0, max = 3, value = 0,{
      incProgress(1, detail = "Transformation" )
      if('asinh' %in% input$transformation_list){
        fcs_data$exprs_data <- exprs_asinh_transformation(exprs_data = fcs_data$exprs_data, use_markers = fcs_data$use_markers,
                                                       cofactor = input$cofactor)
      }
      incProgress(1, detail = "Outlier detection" )
      if('outlier_by_quantile' %in% isolate(input$transformation_list)){
        fcs_data$exprs_data <- exprs_outlier_squeezing(exprs_data = fcs_data$exprs_data, use_markers = fcs_data$use_markers,
                                                    quantile = input$quantile)
      }
      incProgress(1)
    })
  })


  ##### UI for saving cell annotation as fcs files
  output$save_cell_ann_ui <- renderUI({
    if(is.null(fcs_data$cell_ann)){return(NULL)}
    selectInput('save_cell_ann_dp', label = h5("Adding cell annotation to FCS files"),
                choices = colnames(fcs_data$cell_ann), multiple = TRUE)
  })


  ##### Save sample data with cluster info as panal files
  shinyDirChoose(input, 'choose_panel_clust', roots = roots)
  observeEvent(input$dwn_panel_clust, {
    if(is.null(fcs_data$fcs_raw)){return(NULL)}
    panel_folder_path <- parseDirPath(roots, input$choose_panel_clust)
    if(!is.null(panel_folder_path) | length(panel_folder_path) !=0 ){
      filenames <- get_new_fcs_names(old_name = paste0(flowCore::sampleNames(fcs_raw), ".fcs"),
                                     folder_path = panel_folder_path, prefix = input$name_prefix)
      new_fcs <- get_clustered_fcs_files(fcs_raw = fcs_data$fcs_raw, cell_ann = fcs_data$cell_ann,
                                         column_name = input$save_cell_ann_dp)
      flowCore::write.flowSet(new_fcs, outdir = as.character(panel_folder_path), filename = filenames)
    }
  })


  ##### UI for advanced options data preparation
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


  ##### UI to choose marker for scatter plot dp
  output$mk_scatter_dp_ui <- renderUI({
    if(is.null(fcs_data$use_markers)){return(NULL)}
    selectInput('mk_scatter_dp', label = h4("Plotted marker"),
                choices = names(fcs_data$use_markers), selected = 1)
  })


  ##### Download scatter plot data preparation
  output$dwn_scatter_dp <- downloadHandler(
    filename = function() {
      ext <- input$dwn_scatter_dp_ext
      if(is.null(ext)){ext <- "pdf"}
      paste("Scatter_plot_data_preparation", ext, sep = ".") },
    content = function(file) {
      ext <- input$dwn_scatter_dp_ext
      if(is.null(ext)){ext <- "pdf"}
      ggplot2::ggsave(file, plot = plots$scatter_dp, device = ext)
    }
  )


  ########################
  ###      Gating      ###
  ########################

  ##### Create UI to set the gating plot
  output$gating_api_ui <- renderUI({
    if(is.null(gates$gates)){return(NULL)}
    if(is.null(fcs_data$entire_panel)){return(NULL)}
    fluidRow(
      column(1),
      column(10,
             selectInput('gating_subset', label = h5("Data for gating"),
                         choices = colnames(gates$gates), selected = 1),
             selectInput('gating_mk1', label = h5("Marker for X"),
                         choices = names(fcs_data$entire_panel), selected = 1),
             selectInput('gating_mk2', label = h5("Marker for Y"),
                         choices = names(fcs_data$entire_panel), selected = 2)
      )
    )
  })


  ##### Create UI to convert gates to cell annotation
  output$mergeing_gates_ui <- renderUI({
    if(is.null(gates$gates)){return(NULL)}
    fluidRow(
      column(1),
      column(10,
             h4("Convert gates to cell annotation"),
             selectInput("gate_converting_method", label = h5("converting method"),
                         choices = list("Pure selection" = 'pure', "Gates squeezing" = 'squeeze'),selected = 1),
             selectInput("gates_to_convert_gating", label = h5("Choose gates"),
                         choices = gates$antology$name[-1], multiple = TRUE),
             actionButton("convert_gates", label = "Convert")
      )
    )
  })


  ##### Create UI to rename gates
  output$rename_gates_ui <- renderUI({
    if(is.null(input$gated_node_id)){return(NULL)}
    fluidRow(
      column(1),
      column(10,
             h4("Rename gates"),
             textInput('new_gate_name_gating', label = h5("Write new name"),
                       value = as.character(input$gated_node_id)),
             actionButton("rename_gates", label = "Rename")
      )
    )
  })


  ### Make subset of data to gatingplot
  observeEvent(input$butt_plot_for_gating, {
    if(is.null(gates$gates)){return(NULL)}
    if(is.null(input$gating_subset)){return(NULL)}
    if(is.null(input$gating_mk1)){return(NULL)}
    if(is.null(input$gating_mk2)){return(NULL)}
    withProgress(message = "Gate drawing", min =0, max = 3, value = 0,{
      incProgress(1)
      gates$gated_data_subset <- get_exprs_data_for_gating(exprs_data = fcs_data$exprs_data,
                                                     gating_subset = gates$gates[,input$gating_subset],
                                                     gating_mk1 = fcs_data$entire_panel[input$gating_mk1],
                                                     gating_mk2 = fcs_data$entire_panel[input$gating_mk2])
      incProgress(1)
      gates$gated_data_subset <- get_modif_sparse_data_gating(gates$gated_data_subset)
      incProgress(1)
    })
  })


  ### Drawing interactive dencity plot for gating
  output$scatter_plot_gating <- renderPlot({
    if(is.null(gates$gated_data_subset)){return(NULL)}
    colfunc <- grDevices::colorRampPalette(c("black", "red", "yellow"))
    min_mk1 <- min(gates$gated_data_subset$gating_mk1)
    max_mk1 <- max(gates$gated_data_subset$gating_mk1)
    min_mk2 <- min(gates$gated_data_subset$gating_mk2)
    max_mk2 <- max(gates$gated_data_subset$gating_mk2)
    plots$scatter_plot_gating <- ggplot(gates$gated_data_subset, aes(x = gating_mk1, y = gating_mk2)) +
      ylim((min_mk2 - 0.1*(max_mk2 - min_mk2)),
           (max_mk2 + 0.1*(max_mk2 - min_mk2))) +
      xlim((min_mk1 - 0.1*(max_mk1 - min_mk1)),
           (max_mk1 + 0.1*(max_mk1 - min_mk1))) +
      xlab(input$gating_mk1) +
      ylab(input$gating_mk2) +
      geom_point(size = 0.1) +
      stat_density_2d(aes(fill = ..level..), geom = "polygon") +
      scale_fill_gradientn(colours=colfunc(400))+
      geom_density2d(colour="black") +
      theme_bw() +
      theme(legend.position = "none")
    return(plots$scatter_plot_gating)
  })

  ### The object for the anthology graph of gates
  output$gate_antology_graph <- visNetwork::renderVisNetwork({
    if(is.null(gates$antology)){return(NULL)}
    nodes <- data.frame(id = gates$antology$name, label = gates$antology$name)
    edges <- data.frame(from = gates$antology$parent, to = gates$antology$name)
    edges <- edges[!is.na(edges$from),]
    visNetwork::visNetwork(nodes, edges) %>%
      visNetwork::visInteraction(hover = TRUE) %>%
      visNetwork::visEvents(select = "function(nodes) { Shiny.onInputChange('gated_node_id', nodes.nodes);}") %>%
      visNetwork::visHierarchicalLayout()
  })

  ### Allocation a new gate
  observeEvent(input$gete_chosen_cells, {
    if(is.null(input$brush_gating)){return(NULL)}
    if(is.null(gates$gated_data_subset)){return(NULL)}
    new_name <- paste0("Gate", as.character(ncol(gates$gates)+1), "_",
                       input$gating_mk1, "_", input$gating_mk2, "_from_", strsplit(input$gating_subset, "_")[[1]][1])
    if(!is.null(input$new_gate_name) & (input$new_gate_name != "marker1+/marker2+")){new_name <- input$new_gate_name}
    original_cell_coordinates <- brushedPoints(gates$gated_data_subset, input$brush_gating, xvar = "gating_mk1", yvar = "gating_mk2")
    original_cell_coordinates <- original_cell_coordinates$original_cell_coordinates
    gates$gates$new_gate <- FALSE
    gates$gates[original_cell_coordinates, "new_gate"] <- TRUE
    if(any(grepl(new_name, colnames(gates$gates)))){new_name <- paste0(new_name,"_",sum(grepl(new_name, colnames(gates$gates)))+1)}
    colnames(gates$gates)[which(colnames(gates$gates) ==  "new_gate")] <- new_name
    ## Add anthology note of gating for graph
    gates$antology <- rbind(gates$antology, data.frame(name = new_name, parent = input$gating_subset))
    rownames(gates$antology) <- gates$antology$name
  })

  ##### Converting gates to cell annotation data
  observeEvent(input$convert_gates, {
    fcs_data$cell_ann <- get_cell_type_from_gates(gates$gates, input$gates_to_convert_gating,
                                            fcs_data$cell_ann, method = input$gate_converting_method)
  })

  ##### Renew gate reactive object after gate rename
  observeEvent(input$rename_gates, {
    if(is.null(input$new_gate_name_gating)){return(NULL)}
    if(input$new_gate_name_gating == ""){return(NULL)}
    gates$antology$name[gates$antology$name == input$gated_node_id] <- input$new_gate_name_gating
    gates$antology$parent[gates$antology$parent == input$gated_node_id] <- input$new_gate_name_gating
    rownames(gates$antology) <- gates$antology$name
    colnames(gates$gates)[colnames(gates$gates) == input$gated_node_id] <- input$new_gate_name_gating
  })


}


