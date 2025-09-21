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

  # Call the data preparation module
  data_prep_module <- data_preparation_module(input, output, session, fcs_data, data_prep_settings, gates, plots, roots)
  
  ########################
  ###      Gating      ###
  ########################
  
  # Call the gating module
  gating_module <- gating_module(input, output, session, fcs_data, gates, plots, roots)
  
  ########################
  ###    Clustering    ###
  ########################
  
  # Call the clustering module
  clustering_module <- clustering_module(input, output, session, fcs_data, clusters, data_prep_settings, plots, roots)
  
  ########################
  ###   Enrichment     ###
  ########################
  
  # Call the enrichment module
  enrichment_module <- enrichment_module(input, output, session, fcs_data, plots)
  
  ########################
  ####  Annotations   ####
  ########################
  
  # Call the annotations module
  annotations_module <- annotations_module(input, output, session, fcs_data, clusters, data_prep_settings, plots)
  
  ########################
  ####  Cross-panel   ####
  ########################
  
  # Call the cross-panel module
  cross_panel_module_result <- cross_panel_module(input, output, session, fcs_data, plots, roots)
  # Assign the crosspanel reactive values from the module to the main scope
  crosspanel <- cross_panel_module_result$crosspanel
  
  # Any additional server logic that doesn't fit into the modules can go here
  
}


