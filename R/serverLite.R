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
cytofBrowser_server <- function(input, output){

  ########################
  ### Data preparation ###
  ########################
  roots <- c(Home = path.expand("~"), "R Installation" = R.home(), shinyFiles::getVolumes()())
  shinyFiles::shinyFileChoose(input, 'choose_fcs_dp', roots=roots, filetypes=c('', 'fcs'))

  ##### Create "fcs_data" as reactive object to store the CyTOF data
  fcs_data <-reactiveValues()
  data_prep_settings <- reactiveValues(sampling_size = 0.5, fuse = TRUE, size_fuse = 3000, method = "tSNE",
                                       perplexity = 30, theta = 0.5, max_iter = 1000)
  plots <- reactiveValues()


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


  ##### Upload data and automatic pre-processing steps
  observeEvent(input$butt_upload_dproc, {
    print("--0")
    if((length(input$choose_fcs_dp) <= 1)){return(NULL)}
    print("--0.1")
    withProgress(message = "Extraction data", min =0, max = 11, value = 0,{
      ## Get row data fcs files
      fcs_data$md <- get_fcs_metadata(parseFilePaths(roots, input$choose_fcs_dp)$datapath)
      print("--1")
      incProgress(1, detail = "Upload data" )
      fcs_data$fcs_raw <- get_fcs_raw(fcs_data$md)
      print("--2")
      incProgress(1, detail = "Extraction ereachment" )
      fcs_data$exprs_data <- get_exprs_data(fcs_data$fcs_raw)
      print("--3")
      incProgress(1, detail = "Extraction panel")
      fcs_data$panel <- get_fcs_panel(fcs_data$fcs_raw)
      print("--4")
      incProgress(1, detail = "Markers processing" )
      fcs_data$use_markers <- get_use_marker(fcs_data$panel)
      print("--5")
      fcs_data$entire_panel <- get_entire_panel(fcs_raw = fcs_data$fcs_raw)
      incProgress(1, detail = "Cell number calculation" )
      print("--6")
      fcs_data$cell_number <- get_cell_number(fcs_data$fcs_raw)
      incProgress(1, detail = "Creating of cell annotation" )
      print("--7")
      fcs_data$cell_ann <- data.frame(all_cells = rep("cell", sum(fcs_data$cell_number$cell_nmbr)),
                                samples = rep(fcs_data$cell_number$smpl, fcs_data$cell_number$cell_nmbr),
                                row.names = 1:sum(fcs_data$cell_number$cell_nmbr))
      incProgress(1, detail = "Subsampling" )
      print("--8")
      fcs_data$subset_coord <- get_subset_coord(cell_ann = fcs_data$cell_ann,
                                                sampling_size = data_prep_settings$sampling_size,
                                                fuse = data_prep_settings$fuse,
                                                size_fuse = data_prep_settings$size_fuse)
      print("--8.5")
      incProgress(1, detail = "Dimention redicing" )
      fcs_data$cell_ann <- get_dim_reduce(fcs_data$exprs_data, fcs_data$cell_ann,
                                          fcs_data$subset_coord, fcs_data$use_markers,
                                          method = data_prep_settings$method, perplexity = data_prep_settings$perplexity,
                                          theta = data_prep_settings$theta, max_iter = data_prep_settings$max_iter,
                                          pca_param = FALSE, check_duplicates = FALSE, seed = 1234)
      print("--9")
      print(dim(fcs_data$cell_ann[fcs_data$subset_coord,]))
      incProgress(1, detail = "Plotting" )

      incProgress(1, detail = "Extraction fcs cluster info")

      ## Add primer column to cell type data frame and gates data frame

      print("--10")

      incProgress(1)
    })
  })

  ##### UI to choose marker for scatter plot dp
  output$mk_scatter_dp_ui <- renderUI({
    if(is.null(fcs_data$use_markers)){return(NULL)}
    selectInput('mk_scatter_dp', label = h4("Plotted marker"),
                choices = names(fcs_data$use_markers), selected = 1)
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
    print(dim(plot_data))
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
      print("--8")
      ## If subset was changed
      print(c(data_prep_settings$sampling_size, data_prep_settings$fuse) != c(input$sampling_size, input$fuse))
      print(any(c(data_prep_settings$sampling_size, data_prep_settings$fuse) != c(input$sampling_size, input$fuse)))
      if(any(c(data_prep_settings$sampling_size, data_prep_settings$fuse) != c(input$sampling_size, input$fuse))){
        print(input$sampling_size)
        print(data_prep_settings$sampling_size)
        if(!is.null(input$sampling_size)){data_prep_settings$sampling_size <- input$sampling_size}
        if(!is.null(input$fuse)){data_prep_settings$fuse <- input$fuse}
        print(data_prep_settings$sampling_size)
        fcs_data$subset_coord <- get_subset_coord(cell_ann = fcs_data$cell_ann,
                                                  sampling_size = data_prep_settings$sampling_size,
                                                  fuse = data_prep_settings$fuse,
                                                  size_fuse = data_prep_settings$size_fuse)
        print(length(fcs_data$subset_coord))
      }
      ## Repeat dimention reducing
      print("--8.5")
      incProgress(1, detail = "Dimention redicing" )
      #if(any(c(data_prep_settings$method, data_prep_settings$perplexity,
      #     data_prep_settings$theta, data_prep_settings$max_iter) != c(input$method_plot_dp, input$data_prep_perplexity,
      #                                                                 input$data_prep_theta, input$data_prep_max_iter))){}
      print(input$data_prep_perplexity)
      print(data_prep_settings$perplexity)
      dim_reduce_force <- FALSE
      if(!is.null(input$method_plot_dp)){data_prep_settings$method <- input$method_plot_dp; dim_reduce_force <- TRUE}
      if(!is.null(input$data_prep_perplexity)){data_prep_settings$perplexity <- input$data_prep_perplexity; dim_reduce_force <- TRUE}
      if(!is.null(input$data_prep_theta)){data_prep_settings$theta <- input$data_prep_theta; dim_reduce_force <- TRUE}
      if(!is.null(input$data_prep_max_iter)){data_prep_settings$max_iter <- input$data_prep_max_iter; dim_reduce_force <- TRUE}
      print(data_prep_settings$perplexity)
      fcs_data$cell_ann <- get_dim_reduce(fcs_data$exprs_data, fcs_data$cell_ann,
                                          fcs_data$subset_coord, fcs_data$use_markers,force = dim_reduce_force,
                                          method = data_prep_settings$method, perplexity = data_prep_settings$perplexity,
                                          theta = data_prep_settings$theta, max_iter = data_prep_settings$max_iter,
                                          pca_param = FALSE, check_duplicates = FALSE, seed = 1234)
      print("--9")
      print(head(fcs_data$cell_ann))
      print(head(fcs_data$cell_ann[fcs_data$subset_coord,]))
      print(dim(fcs_data$cell_ann[fcs_data$subset_coord,]))
      incProgress(1, detail = "Plotting" )
    })
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
      ggsave(file, plot = plots$scatter_dp, device = ext)
    }
  )






}
