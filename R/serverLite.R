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
#' @import shiny shinyFiles ggplot2 visNetwork
#' @importFrom flowCore sampleNames write.flowSet
#' @importFrom grDevices colorRampPalette
#' @importFrom ComplexHeatmap Heatmap
#'
cytofBrowser_server <- function(input, output){

  ########################
  ### Data preparation ###
  ########################
  roots <- c(Home = path.expand("~"), "R Installation" = R.home(), shinyFiles::getVolumes()())
  shinyFiles::shinyFileChoose(input, 'choose_fcs_dp', roots=roots, filetypes=c('', 'fcs'))

  ##### Create reactive object to store the CyTOF data
  fcs_data <-reactiveValues(trans_to_do = c("asinh"))
  data_prep_settings <- reactiveValues(sampling_size = 0.5, fuse = TRUE, size_fuse = 3000, method = "tSNE",
                                       perplexity = 30, theta = 0.5, max_iter = 1000)
  gates <- reactiveValues()
  clusters <- reactiveValues()
  plots <- reactiveValues()

  ##### Showing selected fcs files
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

  ##### Upload data and automatic pre-processing steps
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


  ##### Drawing the reactive scatter plot to data preparation
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


  ##### Rewrite for scatter plot
  observeEvent(input$redraw_dp, {
    print(input$redraw_dp)
    print(str(input$redraw_dp))

    withProgress(message = "Extraction data", min =0, max = 3, value = 0,{
      incProgress(1, detail = "Subsampling" )
      print("--Repeat sabsampling--")
      ## If subset was changed
      #if(any(c(data_prep_settings$sampling_size, data_prep_settings$fuse) != c(input$sampling_size, input$fuse))){
      #}
      print("===1===")
      if(!is.null(input$sampling_size)){data_prep_settings$sampling_size <- input$sampling_size}
      print("===2===")
      if(!is.null(input$fuse)){data_prep_settings$fuse <- input$fuse}
      print("===3===")
      print(head(fcs_data$cell_ann))
      fcs_data$subset_coord <- get_subset_coord(cell_ann = fcs_data$cell_ann,
                                                sampling_size = data_prep_settings$sampling_size,
                                                fuse = data_prep_settings$fuse,
                                                size_fuse = data_prep_settings$size_fuse)
      print("===4===")
      ## Repeat dimension reducing
      print("--Repeat dim redicing--")
      incProgress(1, detail = "Dimention redicing" )
      dim_reduce_force <- FALSE
      if(!is.null(input$method_plot_dp)){data_prep_settings$method <- input$method_plot_dp; dim_reduce_force <- TRUE}
      print("===5===")
      if(!is.null(input$data_prep_perplexity)){data_prep_settings$perplexity <- input$data_prep_perplexity; dim_reduce_force <- TRUE}
      print("===6===")
      if(!is.null(input$data_prep_theta)){data_prep_settings$theta <- input$data_prep_theta; dim_reduce_force <- TRUE}
      print("===7===")
      if(!is.null(input$data_prep_max_iter)){data_prep_settings$max_iter <- input$data_prep_max_iter; dim_reduce_force <- TRUE}
      print(head(fcs_data$exprs_data))
      print(head(fcs_data$cell_ann))
      print(head(fcs_data$cell_ann))
      print(str(fcs_data$cell_ann))
      fcs_data$cell_ann <- get_dim_reduce(fcs_data$exprs_data, fcs_data$cell_ann,
                                          fcs_data$subset_coord, fcs_data$use_markers,force = dim_reduce_force,
                                          method = data_prep_settings$method, perplexity = data_prep_settings$perplexity,
                                          theta = data_prep_settings$theta, max_iter = data_prep_settings$max_iter,
                                          pca_param = FALSE, check_duplicates = FALSE, seed = 1234)
      print("===9===")
      print(head(fcs_data$cell_ann))
      incProgress(1, detail = "Plotting" )
    })
  })

  ##### Show the implemented transformations
  output$trans_dp <- renderPrint({fcs_data$trans})

  #### Reaction to button "Transform" in "Data processing"
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

  ##### Create UI to choose excluded markers
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

  ##### Update reactive object use_markers after marker changing
  observeEvent(input$change_mk_button, {
    if(!is.null(input$exclude_mk_dp)){
      fcs_data$use_markers <- fcs_data$use_markers[!(names(fcs_data$use_markers) %in% input$exclude_mk_dp)]
    }
    if(!is.null(input$include_mk_dp)){
      fcs_data$use_markers <- fcs_data$entire_panel[c(names(fcs_data$use_markers), input$include_mk_dp)]
    }
  })

  ##### Create UI to choose excluded extracted annotations
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

  ##### Update annotation objects after extaract ann changing
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


  ##### UI for saving cell annotation as fcs files
  output$save_cell_ann_dp_ui <- renderUI({
    if(is.null(fcs_data$cell_ann)){return(NULL)}
    selectInput('save_cell_ann_dp', label = h5("Adding cell annotation to FCS files"),
                choices = colnames(fcs_data$cell_ann), multiple = TRUE)
  })

  ##### First part of UI for subsetting data by annotation
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

  ##### Second part of UI for subsetting data by annotation
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

  ### Subsatting data by annotation
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

  ##### Save data as fcs files
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


  ##### Download scatter plot data preparation
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

  ##### Drawing the reactive plot of cell number
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

  ##### Download plot of cell number data preparation
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

  ##### Create UI to choose target marker for marker density plot
  output$mk_density_dp_ui <- renderUI({
    if(is.null(fcs_data$use_markers)){ return(NULL)}
    selectInput('mk_density_dp', label = h5("Plotted marker"),
                choices = names(fcs_data$use_markers), selected = 1)
  })

  ##### Drawing the reactive histogram plot of marker expression
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

  ##### Download plot of marker histogram data preparation
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

  ##### Reactive show current use_markers from "fcs_data" object
  output$mk_rested_dp <- renderPrint({fcs_data$use_markers})
  output$mk_excluded_dp <- renderPrint({
    fcs_data$entire_panel[!(fcs_data$entire_panel %in% fcs_data$use_markers)]
  })

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
             selectInput("gate_converting_method", label = h5("Converting method"),
                         choices = list("Pure selection" = 'pure', "Gates squeezing" = 'squeeze'),selected = 1),
             selectInput("gates_to_convert_gating", label = h5("Gates to convert"),
                         choices = gates$antology$name[-1], multiple = TRUE),
             textInput("gates_name_new_ann", label = h5("New annotation name"),value = NULL, placeholder = "annotation_from_gates"),
             actionButton("convert_gates", label = "Convert")
      )
    )
  })


  ##### Create UI to rename gates
  output$rename_gates_ui <- renderUI({
    if(is.null(input$gated_node_id)){return(NULL)}
    ui <- fluidRow(
      fluidRow(
        column(1),
        column(10,
               h4("Rename gate"),
               textInput('new_gate_name_gating', label = h5("Write new name"),
                         value = as.character(input$gated_node_id)),
               actionButton("rename_gates", label = "Rename")
        )
      ),
      fluidRow(
        column(1),
        column(6, h5(paste0("Delete gate ", input$gated_node_id))),
        column(4, actionButton("delete_gate", label = "Delete"))
      )
    )
    return(ui)
  })


  ### UI for make subset by gates
  output$subset_gates_ui <- renderUI({
    fluidRow(
      column(1),
      column(10,
             radioGroupButtons(inputId = "gates_subseting_orient", label = h5("Subsatting by gates"),
                               choices = c("exclude", "include"), selected = "exclude", justified = TRUE),
             selectInput("gates_to_subseting", label = h5("Gates to use "),
                         choices = gates$antology$name[-1], multiple = TRUE),
             actionButton("butt_subset_gating", label = "Subset")
      )
    )
  })

  ### Implement subsetting data by gates it all data objects
  observeEvent(input$butt_subset_gating,{
    if(is.null(gates$gates)){return(NULL)}
    if(is.null(input$gates_to_subseting)){return(NULL)}
    if(length(input$gates_to_subseting) == 1){stay_coord <- gates$gates[,input$gates_to_subseting]}
    if(length(input$gates_to_subseting) > 1){stay_coord <- apply(gates$gates[,input$gates_to_subseting], 1, any)}
    if(input$gates_subseting_orient == "exclude"){stay_coord <- !stay_coord}
    fcs_data$exprs_data <- fcs_data$exprs_data[stay_coord,]
    fcs_data$cell_ann <- fcs_data$cell_ann[stay_coord,]
    gates$gates <- gates$gates[stay_coord,]
    fcs_data$subset_coord <- which(stay_coord[as.integer(fcs_data$subset_coord)])
    gates$subset_coord <- which(stay_coord[gates$subset_coord])
  })


  ### Make subset of data to gatingplot
  observeEvent(input$butt_plot_for_gating, {
    if(is.null(gates$gates)){return(NULL)}
    if(is.null(input$gating_subset)){return(NULL)}
    if(is.null(input$gating_mk1)){return(NULL)}
    if(is.null(input$gating_mk2)){return(NULL)}
    withProgress(message = "Making of subset data", min =0, max = 3, value = 0,{
      incProgress(1)
      gates$gated_data_subset <- get_exprs_data_for_gating(exprs_data = fcs_data$exprs_data,
                                                     gating_subset = gates$gates[,input$gating_subset],
                                                     gating_mk1 = fcs_data$entire_panel[input$gating_mk1],
                                                     gating_mk2 = fcs_data$entire_panel[input$gating_mk2])
      incProgress(1)
      gates$gated_data_subset <- get_modif_sparse_data_gating(gates$gated_data_subset)
      gates$subset_coord <- c(1:nrow(gates$gated_data_subset))
      if(input$fuse_gate){gates$subset_coord <- get_size_subset_gating_data(gates$gated_data_subset)}
      incProgress(1)
    })
  })


  ### Drawing interactive density plot for gating
  output$scatter_plot_gating <- renderPlot({
    if(is.null(gates$gated_data_subset)){return(NULL)}
    withProgress(message = "Density drawing", min =0, max = 3, value = 0,{
      colfunc <- grDevices::colorRampPalette(c("black", "red", "yellow"))
      incProgress(1, detail = "Data preparation" )
      min_mk1 <- min(gates$gated_data_subset$gating_mk1)
      max_mk1 <- max(gates$gated_data_subset$gating_mk1)
      min_mk2 <- min(gates$gated_data_subset$gating_mk2)
      max_mk2 <- max(gates$gated_data_subset$gating_mk2)
      incProgress(1, detail = "Drawing")
      plots$scatter_plot_gating <- ggplot(gates$gated_data_subset[gates$subset_coord,], aes(x = gating_mk1, y = gating_mk2)) +
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
      incProgress(1)
    })
    return(plots$scatter_plot_gating)
  })

  ##### Download gates density plot
  output$dwn_density_gate <- downloadHandler(
    filename = function() {
      ext <- input$dwn_density_gate_ext
      if(is.null(ext)){ext <- "pdf"}
      paste("Gates_density_plot", ext, sep = ".") },
    content = function(file) {
      ext <- input$dwn_density_gate_ext
      if(is.null(ext)){ext <- "pdf"}
      ggsave(file, plot = plots$scatter_plot_gating, device = ext,
             width = input$dwn_density_gate_width,height = input$dwn_density_gate_height,
             units = input$dwn_density_gate_units,dpi =input$dwn_density_gate_dpi)
    }
  )

  ##### Drawing the gates overlap plot
  output$gates_overlap <- renderPlot({
    if(is.null(gates$gates)){return(NULL)}
    if(ncol(gates$gates) <= 1){return(NULL)}
    withProgress(message = "Gates overlap drawing", min =0, max = 3, value = 0,{
      incProgress(1, detail = "Overlap finding")
      plot_gates <- get_gates_overlap_data(gates = gates$gates)
      incProgress(1, detail = "Drawing")
      plots$gate_overlap <- ggplot(plot_gates, aes(Gates, Cells, fill= value))+
        geom_tile()+
        scale_fill_manual(values = c("white", "black"))+
        theme_bw()+
        theme(legend.position = "none")+
        theme(axis.text.x = element_text(angle = 90))
      incProgress(1)
    })
    return(plots$gate_overlap)
  })

  ##### Download gates overlap plot
  output$dwn_overlap_gate <- downloadHandler(
    filename = function() {
      ext <- input$dwn_overlap_gate_ext
      if(is.null(ext)){ext <- "pdf"}
      paste("Gates_overlap_plot", ext, sep = ".") },
    content = function(file) {
      ext <- input$dwn_overlap_gate_ext
      if(is.null(ext)){ext <- "pdf"}
      ggsave(file, plot = plots$gate_overlap, device = ext,
             width = input$dwn_overlap_gate_width,height = input$dwn_overlap_gate_height,
             units = input$dwn_overlap_gate_units,dpi =input$dwn_overlap_gate_dpi)
    }
  )

  ##### Drawing the annotation overlap plot in gate section
  output$ann_overlap_gate <- renderPlot({
    if(is.null(fcs_data$cell_ann)){return(NULL)}
    if(ncol(fcs_data$cell_ann) <= 1){return(NULL)}
    plot_ann <- get_ann_overlap_data(ann = fcs_data$cell_ann)
    pl <- ggplot(plot_ann, aes(fill=group, y=value, x = annotations)) +
      geom_bar(position="fill", stat="identity")+
      scale_fill_manual(values=plot_ann$color)+
      theme_bw() +
      ylab("")+
      theme(axis.text.x = element_text(angle = 90))
    if(input$legend_ann_overlap_plot_gate){pl <- pl + theme(legend.position="none")}
    plots$ann_overlap_gate <- pl
    return(plots$ann_overlap_gate)
  })

  ##### Download annotation overlap plot in gate section
  output$dwn_ann_overlap_gate <- downloadHandler(
    filename = function() {
      ext <- input$dwn_ann_overlap_gate_ext
      if(is.null(ext)){ext <- "pdf"}
      paste("Gates_overlap_plot", ext, sep = ".") },
    content = function(file) {
      ext <- input$dwn_ann_overlap_gate_ext
      if(is.null(ext)){ext <- "pdf"}
      ggsave(file, plot = plots$ann_overlap_gate, device = ext,
             width = input$dwn_ann_overlap_gate_width,height = input$dwn_ann_overlap_gate_height,
             units = input$dwn_ann_overlap_gate_units,dpi =input$dwn_ann_overlap_gate_dpi)
    }
  )

  ### The object for the anthology graph of gates
  output$gate_antology_graph <- visNetwork::renderVisNetwork({
    if(is.null(gates$antology)){return(NULL)}
    nodes <- data.frame(id = gates$antology$name, label = gates$antology$name)
    edges <- data.frame(from = gates$antology$parent, to = gates$antology$name)
    edges <- edges[!is.na(edges$from),]
    net <- visNetwork::visNetwork(nodes, edges)
    net <- visNetwork::visInteraction(net, hover = TRUE)
    net <- visNetwork::visEvents(net, select = "function(nodes) { Shiny.onInputChange('gated_node_id', nodes.nodes);}")
    visNetwork::visHierarchicalLayout(net)
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
    if(is.null(input$gate_converting_method)){return(NULL)}
    fcs_data$cell_ann <- get_cell_type_from_gates(gates$gates, input$gates_to_convert_gating,
                                            fcs_data$cell_ann, method = input$gate_converting_method,
                                            name_new_ann = input$gates_name_new_ann)
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

  ##### Delete gate from gate reactive object
  observeEvent(input$delete_gate, {
    if(is.null(input$gated_node_id)){return(NULL)}
    gates$antology$parent[gates$antology$parent == input$gated_node_id] <- gates$antology[input$gated_node_id, "parent"]
    gates$antology <- gates$antology[rownames(gates$antology) != input$gated_node_id,]
    gates$gates[,colnames(gates$gates) == input$gated_node_id] <- NULL
  })

  ########################
  ###    Clustering    ###
  ########################

  ##### Create UI to choose excluded markers from clusterisation
  output$mk_subset_clusters_ui <- renderUI({
    if(is.null(fcs_data$use_markers)){return(NULL)}
    selectInput('exclude_mk_clusters', label = "Exclude markers from clustering",
                choices = names(fcs_data$use_markers), multiple = TRUE)
  })

  ##### Action on the Clustering button
  observeEvent(input$start_clustering, {
    withProgress(message = "Clustering", min =0, max = 4, value = 0,{
      incProgress(1, detail = "clustering")
      ## Clustering with set parameters
      clusters$clust_markers <- fcs_data$use_markers[!(names(fcs_data$use_markers) %in% input$exclude_mk_clusters)]
      if(input$mode_k_choice == 1){k <- input$maxK}
      if(input$mode_k_choice == 2){k <- input$k}
      som <- get_som(fcs_data$exprs_data , fcs_data$fcs_raw, clusters$clust_markers, seed = 1234)
      incProgress(1, detail = "clustering")
      mc <- get_consensusClust(som, maxK = k, seed = 1234)
      incProgress(1, detail = "extration of cluster set")
      clusters$clusters <- get_all_consensusClust(som = som, mc = mc)
      if(input$mode_k_choice == 1){k <- get_optimal_clusters(mc, rate_var_expl = input$rate_var_explan)}
      fcs_data$cell_ann$clusters <- clusters$clusters[,as.character(k)]
      incProgress(1)
    })
  })

  ##### Update reactive objects euclid_dist
  observe({
    if(is.null(fcs_data$cell_ann$clusters)){return(NULL)}
    withProgress(message = "Cluster estimation", min =0, max = 2, value = 0,{
      incProgress(1, detail = "distance between clusters")
      clusters$clus_euclid_dist <- get_euclid_dist(exprs_data = isolate(fcs_data$exprs_data), use_markers = isolate(fcs_data$use_markers),
                                                   cell_clustering = fcs_data$cell_ann$clusters)
      incProgress(1, detail = "graph elements")
    })
  })

  ##### Update reactive objects with nodes and edges
  observe({
    if(is.null(fcs_data$cell_ann$clusters)){return(NULL)}
    withProgress(message = "Cluster estimation", min =0, max = 3, value = 0,{
      incProgress(1, detail = "graph elements")
      clusters$edges <- get_edges(clusters$clus_euclid_dist)
      incProgress(1, detail = "drawing scatter plot")
      clusters$nodes <- get_nodes(clusters$edges, fcs_data$cell_ann$clusters)
      incProgress(1)
    })
  })

  ##### Create UI to choose clusters to merge
  output$mergeing_clust_ui <- renderUI({
    if(is.null(fcs_data$cell_ann$clusters)){return(NULL)}
    fluidRow(
      column(1),
      column(10,
             selectInput('cluster_to_merge_clust', label = h5("Merged clusters"),
                         choices = unique(fcs_data$cell_ann$clusters), multiple = TRUE),
             actionButton("merge_clust", label = "Merge")
      )
    )
  })

  ##### Renew clustering reactive object after merging
  observeEvent(input$merge_clust, {
    fcs_data$cell_ann$clusters <- cluster_merging(fcs_data$cell_ann$clusters, input$cluster_to_merge_clust)
  })

  ##### Create UI to rename clusters
  output$rename_clust_ui <- renderUI({
    if(is.null(input$current_node_id)){return(NULL)}
    fluidRow(
      column(1),
      column(10,
             textInput("new_cluster_name_clust", label = h5("Rename chosen cluster"),
                       value = as.character(input$current_node_id)),
             actionButton("rename_clust", label = "Rename")
      )
    )
  })

  ##### Renew cluster reactive object after cluster rename
  observeEvent(input$rename_clust, {
    if(is.null(input$new_cluster_name_clust)){return(NULL)}
    if(input$new_cluster_name_clust == ""){return(NULL)}
    fcs_data$cell_ann$clusters[fcs_data$cell_ann$clusters==input$current_node_id] <- input$new_cluster_name_clust
  })

  ##### UI to choose marker for scatter plot with clusters
  output$mk_scatter_clust_ui <- renderUI({
    if(is.null(fcs_data$cell_ann$clusters)){return(NULL)}
    ann_groups <- colnames(fcs_data$cell_ann)
    ann_groups <- ann_groups[!grepl("tSNE", ann_groups)]
    ann_groups <- ann_groups[!grepl("UMAP", ann_groups)]
    ann_groups <- ann_groups[!grepl("all_cells", ann_groups)]
    selected_1 <- 1
    if("clusters" %in% ann_groups){selected_1 <- "clusters"}
    selectInput("mk_target_clusters", label = h5("Plotted marker"),
                choices = c(ann_groups, names(fcs_data$use_markers)),selected = selected_1)
  })

  ##### UI for advanced options for clustering page
  output$advanced_opt_clust_ui <- renderUI({
    if(is.null(input$method_plot_clust)){return(NULL)}
    if(is.null(input$method_plot_clust)){return(NULL)}
    if(input$method_plot_clust == 'tSNE'){
      ui <- fluidRow(
        column(1),
        column(10,
               numericInput("perplexity_clust", "tSNE Perplexity", min = 0, max = 200, value = 30, step = 5),
               numericInput("theta_clust", "tSNE Theta", min = 0, max = 1, value = 0.5, step = 0.1),
               numericInput("max_iter_clust", "tSNE Iterations", value = 1000, step = 500)
        )
      )
    }
    if(input$method_plot_clust == 'UMAP'){ui <- NULL}
    return(ui)
  })

  ##### Drawing the reactive and interactive UMAP plot
  output$scatter_plot_clust <- renderPlot({
    if(is.null(fcs_data$cell_ann$clusters)){return(NULL)}
    if(is.null(fcs_data$subset_coord)){return(NULL)}
    color_mk <- "clusters"
    if(!is.null(input$mk_target_clusters)){color_mk <- input$mk_target_clusters}
    if(color_mk %in% names(fcs_data$use_markers)){color_mk <- fcs_data$use_markers[color_mk]}
    if(data_prep_settings$method == "tSNE"){
      plot_data <- data.frame(X = fcs_data$cell_ann[fcs_data$subset_coord, "tSNE1"],
                              Y = fcs_data$cell_ann[fcs_data$subset_coord, "tSNE2"])}
    if(data_prep_settings$method == "UMAP"){
      plot_data <- data.frame(X = fcs_data$cell_ann[fcs_data$subset_coord, "UMAP1"],
                              Y = fcs_data$cell_ann[fcs_data$subset_coord, "UMAP2"])}
    if(color_mk %in% colnames(fcs_data$cell_ann)){plot_data$mk <- as.factor(fcs_data$cell_ann[fcs_data$subset_coord, color_mk])}
    if(color_mk %in% colnames(fcs_data$exprs_data)){plot_data$mk <- fcs_data$exprs_data[fcs_data$subset_coord, color_mk]}
    focus_node <- input$current_node_id
    plt <- ggplot2::ggplot(plot_data,  aes(x = X, y = Y, color = mk)) +
      geom_point(size = input$point_size_clust)
    if(color_mk %in% colnames(fcs_data$cell_ann)){plt <- plt + scale_color_manual(values = as.character(clusters$nodes$color))+
      guides(colour = guide_legend(override.aes = list(size=2)))}
    if(color_mk %in% colnames(fcs_data$exprs_data)){plt <- plt + scale_color_gradient2(midpoint=0.5, low='blue', mid='gray', high='red')}
    plt <- plt + geom_point(data = plot_data[fcs_data$cell_ann[fcs_data$subset_coord, "clusters"] == focus_node,],
                            colour = 'black', size = (input$point_size_clust*1.5))+
      labs(color = color_mk) +
      theme_bw()
    plots$scatter_clust <- plt
    return(plt)
  })

  ##### Rewrite for scatter plot in clusteringpage
  observeEvent(input$redraw_clust, {
    withProgress(message = "Extraction data", min =0, max = 3, value = 0,{
      incProgress(1, detail = "Subsampling" )
      ## If subset was changed
      #if(any(c(data_prep_settings$sampling_size, data_prep_settings$fuse) != c(input$sampling_size_clust, input$fuse_clust))){
      #}
      if(!is.null(input$sampling_size_clust)){data_prep_settings$sampling_size <- input$sampling_size_clust}
      if(!is.null(input$fuse_clust)){data_prep_settings$fuse <- input$fuse_clust}
      fcs_data$subset_coord <- get_subset_coord(cell_ann = fcs_data$cell_ann,
                                                sampling_size = data_prep_settings$sampling_size,
                                                fuse = data_prep_settings$fuse,
                                                size_fuse = data_prep_settings$size_fuse)
      ## Repeat dimension reducing
      incProgress(1, detail = "Dimention redicing" )
      dim_reduce_force <- FALSE
      if(!is.null(input$method_plot_clust)){data_prep_settings$method <- input$method_plot_clust; dim_reduce_force <- TRUE}
      if(!is.null(input$perplexity_clust)){data_prep_settings$perplexity <- input$perplexity_clust; dim_reduce_force <- TRUE}
      if(!is.null(input$theta_clust)){data_prep_settings$theta <- input$theta_clust; dim_reduce_force <- TRUE}
      if(!is.null(input$max_iter_clust)){data_prep_settings$max_iter <- input$max_iter_clust; dim_reduce_force <- TRUE}
      fcs_data$cell_ann <- get_dim_reduce(fcs_data$exprs_data, fcs_data$cell_ann,
                                          fcs_data$subset_coord, fcs_data$use_markers,force = dim_reduce_force,
                                          method = data_prep_settings$method, perplexity = data_prep_settings$perplexity,
                                          theta = data_prep_settings$theta, max_iter = data_prep_settings$max_iter,
                                          pca_param = FALSE, check_duplicates = FALSE, seed = 1234)
      incProgress(1, detail = "Plotting" )
    })
  })

  ##### Download scatter plot data preparation
  output$dwn_scatter_clust <- downloadHandler(
    filename = function() {
      ext <- input$dwn_scatter_clust_ext
      if(is.null(ext)){ext <- "pdf"}
      paste("Scatter_plot_clustering", ext, sep = ".") },
    content = function(file) {
      ext <- input$dwn_scatter_clust_ext
      if(is.null(ext)){ext <- "pdf"}
      ggplot2::ggsave(file, plot = plots$scatter_clust, device = ext,
                      width = input$dwn_scatter_clust_width,height = input$dwn_scatter_clust_height,
                      units = input$dwn_scatter_clust_units,dpi =input$dwn_scatter_clust_dpi)
    }
  )

  ##### Reactive show current use_markers from "fcs_data" object
  output$mk_clusted_clust <- renderPrint({clusters$clust_markers})
  output$mk_rested_clust <- renderPrint({fcs_data$use_markers})
  output$mk_excluded_clust <- renderPrint({
    fcs_data$entire_panel[!(fcs_data$entire_panel %in% fcs_data$use_markers)]
  })

  ##### Drawing the reactive abundance plot
  output$abundance_clust <- renderPlot({
    if(is.null(clusters$nodes$color)){return(NULL)}
    abundance_df <- get_abundance_dataframe(cell_ann = fcs_data$cell_ann,
                                            samples_col = "samples", clusters_col = "clusters")
    plots$abundance_clust <- ggplot(abundance_df, aes(x = samples, y = Freq, fill = clusters)) +
      geom_bar(stat = 'identity') +
      scale_fill_manual(values = as.character(clusters$nodes$color))+
      ylab("cell abundance")+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 90))
    return(plots$abundance_clust)
  })

  ##### Download abundance plot for clustering
  output$dwn_abundance_clust <- downloadHandler(
    filename = function() {
      ext <- input$dwn_abundance_clust_ext
      if(is.null(ext)){ext <- "pdf"}
      paste("Abundance_plot_clustering", ext, sep = ".") },
    content = function(file) {
      ext <- input$dwn_abundance_clust_ext
      if(is.null(ext)){ext <- "pdf"}
      ggsave(file, plot = plots$abundance_clust, device = ext,
             width = input$dwn_abundance_clust_width,height = input$dwn_abundance_clust_height,
             units = input$dwn_abundance_clust_units,dpi =input$dwn_abundance_clust_dpi)
    }
  )

  ##### Drawing the clusters overlap plot in cluster section
  output$clusters_overlap_clust <- renderPlot({
    if(is.null(clusters$clusters)){return(NULL)}
    if(ncol(clusters$clusters) <= 1){return(NULL)}
    plot_ann <- get_ann_overlap_data(ann = clusters$clusters)
    plots$clusters_overlap_clust <- ggplot(plot_ann, aes(fill=group, y=value, x = annotations)) +
      geom_bar(position="fill", stat="identity")+
      scale_fill_manual(values=plot_ann$color)+
      theme_bw() +
      ylab("")+
      xlab("clustering")+
      theme(axis.text.x = element_text(angle = 90))+
      theme(legend.position="none")
    return(plots$clusters_overlap_clust)
  })

  ##### Download clusters overlap plot in cluster section
  output$dwn_clusters_overlap_clust <- downloadHandler(
    filename = function() {
      ext <- input$dwn_clusters_overlap_clust_ext
      if(is.null(ext)){ext <- "pdf"}
      paste("Gates_overlap_plot", ext, sep = ".") },
    content = function(file) {
      ext <- input$dwn_clusters_overlap_clust_ext
      if(is.null(ext)){ext <- "pdf"}
      ggsave(file, plot = plots$clusters_overlap_clust, device = ext,
             width = input$dwn_clusters_overlap_clust_width,height = input$dwn_clusters_overlap_clust_height,
             units = input$dwn_clusters_overlap_clust_units,dpi =input$dwn_clusters_overlap_clust_dpi)
    }
  )

  ##### Drawing the reactive and interactive graph with clusters
  output$clusters_graph <- renderVisNetwork({
    if(is.null(clusters$nodes)){return(NULL)}
    edges_threshold <- input$edges_threshold_clusterisation
    if(is.null(input$edges_threshold_clusterisation)){edges_threshold <- 0.5}
    gravity <- input$gravity_clusterisation
    if(is.null(input$gravity_clusterisation)){gravity <- -40}
    edges <- filter_edges(clusters$edges, edges_threshold)

    net <- visNetwork::visNetwork(clusters$nodes, edges)
    net <- visNetwork::visInteraction(net, hover = TRUE)
    net <- visNetwork::visEvents(net, select = "function(nodes) { Shiny.onInputChange('current_node_id', nodes.nodes);}")
    net <- visNetwork::visPhysics(net, solver = "forceAtlas2Based", forceAtlas2Based = list(gravitationalConstant = gravity))
    return(net)
  })

  ########################
  ###   Enrichment     ###
  ########################

  ##### UI to chose annotation for heatmap
  output$hm_rows_enrich_ui <- renderUI({
    if(is.null(fcs_data$cell_ann)){return(NULL)}
    col_names <- colnames(fcs_data$cell_ann)
    col_names <- col_names[!(col_names %in% c("all_cells", "tSNE1", "tSNE2", "UMAP1", "UMAP2"))]
    if(length(col_names) < 1){return(NULL)}
    sel <- col_names[1]
    if("clusters" %in% col_names){sel <- "clusters"}
    selectInput('hm_rows_enrich', label = h5("Heatmap rows"), choices = col_names, selected = sel)
  })

  ##### Drawing heatmap for enrichments
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

  ##### Download clusters overlap plot in cluster section
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

  ##### Create UI to choose marker for deconvolution plot
  output$mk_deconvol_enrich_ui <- renderUI({
    if(is.null(fcs_data$use_markers)){return(NULL)}
    selectInput("mk_deconvol_enrich", label = h5("Marker for deconvolution"),
                choices = names(fcs_data$use_markers),
                selected = 1)
  })
  ##### UI to chose annotation for deconvolution plot
  output$gp_deconvol_enrich_ui <- renderUI({
    if(is.null(fcs_data$cell_ann)){return(NULL)}
    col_names <- colnames(fcs_data$cell_ann)
    col_names <- col_names[!(col_names %in% c("all_cells", "tSNE1", "tSNE2", "UMAP1", "UMAP2", "samples"))]
    if(length(col_names) < 1){return(NULL)}
    sel <- col_names[1]
    if("clusters" %in% col_names){sel <- "clusters"}
    selectInput('gp_deconvol_enrich', label = h5("Annotation to plot"), choices = col_names, selected = sel)
  })


  ##### Drawing deconvolution of marker to sample-cluster plot
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

  ##### Download clusters overlap plot in cluster section
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

  ########################
  ####  Annotations   ####
  ########################

  ##### UI to chose annotation for rename annotaions
  output$chose_managment_ann_ui <- renderUI({
    if(is.null(fcs_data$cell_ann)){return(NULL)}
    col_names <- colnames(fcs_data$cell_ann)
    col_names <- col_names[!(col_names %in% c("all_cells", "tSNE1", "tSNE2", "UMAP1", "UMAP2", "samples"))]
    if(length(col_names) < 1){return(NULL)}
    sel <- col_names[1]
    if("clusters" %in% col_names){sel <- "clusters"}
    selectInput('gp_managment_ann', label = h5("Annotation to manage"), choices = col_names, selected = sel)
  })

  ##### Dynamic UI for rename each group in chosen annotation
  output$managment_ann_ui <- renderUI({
    if(is.null(input$gp_managment_ann)){return(NULL)}
    lapply(unique(fcs_data$cell_ann[,input$gp_managment_ann]), function(i) {
      textInput(paste0('group', i), label = NULL, value = i)
    })
  })

  ##### change name in annotations
  observeEvent(input$rename_groups, {
    for (i in unique(fcs_data$cell_ann[,input$gp_managment_ann])){
      tmp <- fcs_data$cell_ann[,input$gp_managment_ann]
      tmp[tmp == i] <- input[[paste0('group', i)]]
      fcs_data$cell_ann[,input$gp_managment_ann] <- tmp
    }
  })



  ##### UI for advanced options annotatios page
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

  ##### Isolate number of plot from continiouse rection
  annotations <- reactiveValues(num_plot = 4)
  observeEvent(input$redraw_ann, {annotations$num_plot <- input$num_plot_ann})

  ##### UI to choose marker for scatter plots for annotating
  observe({
    lapply(1:annotations$num_plot, function(x){
      output[[paste0("mk_ann_ui", x)]] <- renderUI({
        selectInput(paste0("mk_target_ann", x), label = NULL,
                    choices = c(input$gp_managment_ann, names(fcs_data$use_markers)),selected = 1)
      })
    })
  })


  ##### Drawing scatter plots for annotating
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

  ##### Rewrite for scatter plot in Annotation
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

  ########################
  ####  Cross-panel   ####
  ########################

  ##### Creating reactive Values with panel data
  crosspanel <- reactiveValues(md_set = list(), fcs_raw_set = list(), panel_set = list(), use_marker_set = list(),
                               cell_clustering_list_set = list())

  shinyFileChoose(input, 'choose_fcs_cross_p', roots=roots, filetypes=c('', 'fcs'))
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

  ##### Show the overlapped samples of the panels
  observe({
    if(length(crosspanel$md_set) == 0){return(NULL)}
    panel_content <- get_panel_content(crosspanel$fcs_raw_set)
    crosspanel$use_samples <- rownames(panel_content)[apply(panel_content, 1, function(x) {all(x == "presented")})]
  })

  ##### Show the overlapped markers of the panels
  observe({
    if(length(crosspanel$use_marker_set) == 0){return(NULL)}
    mk_panel_content <- get_common_markers(crosspanel$use_marker_set)
    crosspanel$use_common_mk <- rownames(mk_panel_content)[apply(mk_panel_content, 1, function(x) {all(x == "presented")})]
  })

  ##### Create ui for removing panels
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

  ##### Removing panels
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

  ##### Show the sample content of the panels
  output$content_cross_p <- renderPlot({
    if(length(crosspanel$md_set) < 1){return(NULL)}
    panel_content <- get_panel_content(crosspanel$fcs_raw_set)
    ComplexHeatmap::Heatmap(panel_content, cluster_rows = F, cluster_columns = F,
                            row_names_side = "left", column_names_side = "top",
                            col = list("absent" = 'red3', "presented" = 'green3'), rect_gp = grid::gpar(col = "white", lwd = 2),
                            show_heatmap_legend = F)

  })

  ##### Show the marker content of the panels
  output$mk_content_cross_p <- renderPlot({
    if(length(crosspanel$use_marker_set) < 1){return(NULL)}
    mk_panel_content <- get_common_markers(crosspanel$use_marker_set)
    ComplexHeatmap::Heatmap(mk_panel_content, cluster_rows = F, cluster_columns = F,
                            row_names_side = "left", column_names_side = "top", row_names_gp = grid::gpar(fontsize = 8),
                            col = list("absent" = 'red3', "presented" = 'green3'), rect_gp = grid::gpar(col = "white", lwd = 2),
                            show_heatmap_legend = F)

  })

  ##### UI for choosing p-valueadjusting method for correlation cross-panel analysis
  output$abund_corr_padj_method_cross_p_ui <- renderUI({
    if(input$abund_corr_p_mode_cross_p == 'padj'){return((
      selectInput('abund_corr_padj_method_cross_p', "Select method for p-value adjusting",
                  choices = c("BH" = 'BH', "BY" = 'BY', "holm" ='holm', "hochberg" = 'hochberg', "hommel" = 'hommel',
                              "bonferroni" = 'bonferroni'), selected = 'BH')
    ))}
  })

  ##### Make abundance cross panel correlation analysis
  observeEvent(input$abund_corr_cross_p, {
    crosspanel$abund_cross_p <- get_abund_cross_panel(cell_clustering_list_set = crosspanel$cell_clustering_list_set,
                                                      use_samples = crosspanel$use_samples)
    crosspanel$abund_corr_info <- get_matrices_corr_info(crosspanel$abund_cross_p, method_cortest = input$abund_corr_method_cross_p,
                                                         method_padj = input$abund_corr_padj_method_cross_p)
  })

  #### Drawing abundance correlation plot for cross panel analysis
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

  ##### Download abundance correlation plot for cross panel analysis
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

  ##### Download abundance corelation tables  for cross-panel analusis
  output$dwn_table_abund_corr_cross_p <- downloadHandler(
    filename = function() {paste0("cross_panel_", input$dwn_table_abund_corr_cross_p_ext, ".csv")},
    content = function(file) {
      mode <- input$dwn_table_abund_corr_cross_p_ext
      if(is.null(mode)){mode <- 'corr_data'}
      if(mode == 'abund_data'){write.csv(crosspanel$abund_cross_p, file)}
      if(mode == 'corr_data'){write.csv(get_write_corr_info_cross_p(crosspanel$abund_corr_info), file)}
    }
  )

  ##### UI to choose marker for expressing cell fraction correlations
  output$mk_exp_cell_f_cross_p_ui <- renderUI({
    if(is.null(crosspanel$use_common_mk)){return(NULL)}
    selectInput('mk_exp_cell_f_cross_p', label = h5("Marker for cross-panel analysis"),
                choices = crosspanel$use_common_mk,
                selected = 1)
  })

  ##### UI to choose the threshold to expressing cell fraction allocation
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

  ##### UI for choosing p-value adjusting method for exp cell fraction cross-panel analysis
  output$exp_cell_f_corr_padj_method_cross_p_ui <- renderUI({
    if(input$exp_cell_f_corr_p_mode_cross_p == 'padj'){return((
      selectInput('exp_cell_f_corr_padj_method_cross_p', "Select method for p-value adjusting",
                  choices = c("BH" = 'BH', "BY" = 'BY', "holm" ='holm', "hochberg" = 'hochberg', "hommel" = 'hommel',
                              "bonferroni" = 'bonferroni'), selected = 'BH')
    ))}
  })

  ##### Make expressing cell fraction cross panel correlation analysis
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

  #### Drawing abundance correlation plot for cross panel analysis
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

  ##### Drawing expression density plot for a marker in cross-panel analysis
  output$mk_density_plot_cross_p <- renderPlot({
    if(is.null(crosspanel$fcs_raw_set)){return(NULL)}
    if(is.null(input$mk_exp_cell_f_cross_p)){return(NULL)}
    mk_exp_data_cross_p <- get_mk_exp_data_cross_p(fcs_raw_set = crosspanel$fcs_raw_set , use_samples = crosspanel$use_samples,
                                                   use_marker_set = crosspanel$use_marker_set, target_mk = input$mk_exp_cell_f_cross_p)
    ggplot(data.frame(expr = mk_exp_data_cross_p), aes(x = expr)) +
      geom_density(fill = 'black')
  })


  ##### Download abundance correlation plot for cross panel analysis
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

  ##### Download expressing cell fraction corelation tables for cross-panel analusis
  output$dwn_table_exp_cell_f_corr_cross_p <- downloadHandler(
    filename = function() {paste0("cross_panel_",input$mk_exp_cell_f_cross_p,"_", input$dwn_table_exp_cell_f_corr_cross_p_ext, ".csv")},
    content = function(file) {
      mode <- input$dwn_table_exp_cell_f_corr_cross_p_ext
      if(is.null(mode)){mode <- 'corr_data'}
      if(mode == 'exp_cell_f_data'){write.csv(crosspanel$exp_cell_f_cross_p, file)}
      if(mode == 'corr_data'){write.csv(get_write_corr_info_cross_p(crosspanel$exp_cell_f_corr_info), file)}
    }
  )


}


