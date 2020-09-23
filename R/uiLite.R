#' The main function to run cutofCore graphic interface
#' @description Run graphical user interface for cytofBrowser.
#' The function runs the Shiny App in browser. The current package
#' was created with assumption that it will be run with GUI mode
#' as Shiny App, That is why the function is the most appropriate
#' way to use the cytofBrowser.
#'
#' @return The function runs the Shiny App in browser.
#' @import shiny shinyFiles shinydashboard shinyWidgets
#' @export
#'
cytofBrowserGUI <-function(){
  cytofBrowser_ui <- shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(title = "cytofBrowser"),
    shinydashboard::dashboardSidebar(
      shinydashboard::sidebarMenu(
        shinydashboard::menuItem("Data", tabName = 'data_processing', icon = icon("bars")),
        shinydashboard::menuItem("Gating", tabName = 'gating', icon = icon("object-group")),
        shinydashboard::menuItem("Clustering", tabName = 'data_clustering', icon = icon("spinner")),
        shinydashboard::menuItem("Exploration", tabName = 'data_exploration', icon = icon("chart-area")),
        shinydashboard::menuItem("Correlation", tabName = 'data_correlation', icon = icon("braille")),
        shinydashboard::menuItem("Cross-panel", tabName = 'data_crosspanel', icon = icon("clone"))
      )
    ),
    shinydashboard::dashboardBody(
      shinydashboard::tabItems(
        ##################################
        ##### First tab content      #####
        ##################################
        shinydashboard::tabItem(tabName = 'data_processing',
                shiny::fluidRow(
                  shinydashboard::tabBox(
                    shiny::tabPanel("Uploading",
                                    uiOutput('selected_fcs_dp_ui'),
                                    fluidRow(
                                      column(1),
                                      column(5,shinyFiles::shinyFilesButton('choose_fcs_dp', label= 'Select FCS files',
                                                                            title='Please select FCS files', multiple=TRUE)),
                                      column(5, actionButton('butt_upload_dproc', label = "Upload files"))
                                    ),
                                    hr(),
                                    conditionalPanel(
                                      condition = "input.extr_clust_dproc == true",
                                      textInput("extr_clust_pattern_dproc",
                                                label = h5("full or part column name with clusters info (for cytofBrowser and cytofkit : <cluster>)"),
                                                value = "cluster")),

                    ),
                    shiny::tabPanel("Transforming",
                                    checkboxGroupInput("transformation_list", label = h4("Transformations"),
                                                       choices = list("asinh" = 'asinh', "outlier squeezing" = 'outlier_by_quantile',
                                                                      "extract cluster info" = 'extract_cluster_info'),
                                                       selected = c("asinh")),
                                    conditionalPanel(
                                      condition = "input.transformation_list.includes('asinh')",
                                      numericInput("cofactor", label = h4("Cofactor for dividing"), value = 5)
                                    ),
                                    conditionalPanel(
                                      condition = "input.transformation_list.includes('outlier_by_quantile')",
                                      numericInput("quantile", label = h4("Quantile for outlier detection"), value = 0.01)
                                    ),
                                    hr(),
                                    actionButton('butt_trans_dproc', label = "Transform")
                    ),
                    shiny::tabPanel("Markers",
                                    uiOutput('mk_subset_dp_ui')
                    ),
                    shiny::tabPanel("Save",
                                    h5("Saving the data as FCS files "),
                                    shinyDirButton('choose_dwn_folder_dp', label = "Folder choose", title = "Folder for saving FCS files"),
                                    hr(),
                                    h5("Saved set of files can be used as panel data in a cross-panel analysis"),
                                    textInput("name_prefix_dp", label = h5("Name prefix")),
                                    hr(),
                                    uiOutput('save_cell_ann_dp_ui'),
                                    actionButton('dwn_fcs_dp', label = "Save")
                    )
                  ),
                  shinydashboard::box(
                    fluidRow(
                      column(2, actionBttn(inputId = "redraw_dp", style = "material-circle", color = "default" ,icon = icon("redo"))),
                      column(2,
                             dropdownButton(
                               tags$h4("Options of plotting"),
                               numericInput('sampling_size', label = h5("Cell fraction to display"), value = 0.5, step = 0.1),
                               materialSwitch(inputId = 'fuse', label = h4("Size fuse"), value = TRUE),
                               selectInput("method_plot_dp", label = h5("Visualisation method"),
                                           choices = list("tSNE" = "tSNE", "UMAP" = "UMAP"),selected = "tSNE"),
                               numericInput('point_size', label = h5("Size of points"), value = 0.3, step = 0.1),
                               icon = icon("edit"), status = "primary", tooltip = tooltipOptions(title = "plot setting")
                             )
                      ),
                      column(2,
                             dropdownButton(
                               tags$h4("Advanced options"),
                               uiOutput('advanced_opt_dp_ui'),
                               icon = icon("gear"), status = "primary", tooltip = tooltipOptions(title = "plot setting")
                             )
                      ),
                      column(2,
                             dropdownButton(
                               selectInput('dwn_scatter_dp_ext', label = NULL,
                                           choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png",
                                                          'tiff' = "tiff", 'svg' = "svg", 'bmp' = "bmp")),
                               downloadButton('dwn_scatter_dp', ""),
                               hr(),
                               numericInput("dwn_scatter_dp_width", label = h5("width"), value = 10),
                               numericInput("dwn_scatter_dp_height", label = h5("height"), value = 8),
                               selectInput("dwn_scatter_dp_units", label = h5("units"),
                                           choices = list("in" = "in", "cm" = "cm", "mm" = "mm"), selected = "in"),
                               numericInput("dwn_scatter_dp_dpi", label = h5("dpi"), value = 300),
                               icon = icon("save"), status = "primary", tooltip = tooltipOptions(title = "save plot")
                             )
                      )
                    ),
                    plotOutput('scatter_plot_dp'),
                    uiOutput('mk_scatter_dp_ui')
                  ),
                  tabBox(
                    tabPanel("Cell number",
                             fluidRow(
                               column(10,plotOutput("smpl_hist_preparation")),
                               column(2,
                                      dropdownButton(
                                        selectInput('dwn_smpl_hist_dp_ext', label = NULL,
                                                    choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png",
                                                                   'tiff' = "tiff", 'svg' = "svg", 'bmp' = "bmp")),
                                        downloadButton('dwn_smpl_hist_dp', ""),
                                        hr(),
                                        numericInput("dwn_smpl_hist_dp_width", label = h5("width"), value = 10),
                                        numericInput("dwn_smpl_hist_dp_height", label = h5("height"), value = 8),
                                        selectInput("dwn_smpl_hist_dp_units", label = h5("units"),
                                                    choices = list("in" = "in", "cm" = "cm", "mm" = "mm"), selected = "in"),
                                        numericInput("dwn_smpl_hist_dp_dpi", label = h5("dpi"), value = 300),
                                        icon = icon("save"), status = "primary", tooltip = tooltipOptions(title = "save plot")
                                      )
                               )
                             )

                    ),
                    tabPanel("Expression",
                             fluidRow(
                               column(10,
                                      uiOutput('mk_density_dp_ui'),
                                      plotOutput('mk_density_plot_dp')
                                      ),
                               column(2,
                                      dropdownButton(
                                        tags$h4("Advanced options"),
                                        checkboxInput("mk_zero_del_density_plot_dp", label = "Remove zero expressions", value = FALSE),
                                        icon = icon("gear"), status = "primary", tooltip = tooltipOptions(title = "plot setting")
                                      ),
                                      dropdownButton(
                                        selectInput('dwn_mk_density_dp_ext', label = NULL,
                                                    choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png",
                                                                   'tiff' = "tiff", 'svg' = "svg", 'bmp' = "bmp")),
                                        downloadButton('dwn_mk_density_dp', ""),
                                        hr(),
                                        numericInput("dwn_mk_density_dp_width", label = h5("width"), value = 10),
                                        numericInput("dwn_mk_density_dp_height", label = h5("height"), value = 8),
                                        selectInput("dwn_mk_density_dp_units", label = h5("units"),
                                                    choices = list("in" = "in", "cm" = "cm", "mm" = "mm"), selected = "in"),
                                        numericInput("dwn_mk_density_dp_dpi", label = h5("dpi"), value = 300),
                                        icon = icon("save"), status = "primary", tooltip = tooltipOptions(title = "save plot")
                                      )
                               )
                             )
                    ),
                    tabPanel("Markers",
                             h4("Used markers"),
                             verbatimTextOutput('mk_rested_dp'),
                             h4("Unused markers"),
                             verbatimTextOutput('mk_excluded_dp')
                    )
                  ),
                )
        ),

        ###############################
        ##### Second tab content  #####
        ###############################
        shinydashboard::tabItem(tabName = 'gating',
                fluidRow(
                  tabBox(
                    tabPanel("Gating",
                             uiOutput('gating_api_ui'),
                             actionButton("butt_plot_for_gating", label = "Plot for gating")
                    ),
                    tabPanel("Gate management",
                             uiOutput("mergeing_gates_ui"),
                             hr(),
                             uiOutput("rename_gates_ui")
                    )
                  ),
                  shinydashboard::box(
                    plotOutput('scatter_plot_gating', brush = "brush_gating"),
                    fluidRow(
                      column(6, textInput('new_gate_name', label = h5("Name of new gated cells"), value = "marker1+/marker2+")),
                      column(6, actionButton("gete_chosen_cells", label = "Gate chosen cells"))
                    )

                  )
                ),
                fluidRow(
                  tabBox(
                    tabPanel("Gates overlap",
                             fluidRow(
                               column(10, plotOutput('gates_overlap')),
                               column(2,
                                      dropdownButton(
                                        selectInput('dwn_overlap_gate_ext', label = NULL,
                                                    choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png",
                                                                   'tiff' = "tiff", 'svg' = "svg", 'bmp' = "bmp")),
                                        downloadButton('dwn_overlap_gate', ""),
                                        hr(),
                                        numericInput("dwn_overlap_gate_width", label = h5("width"), value = 10),
                                        numericInput("dwn_overlap_gate_height", label = h5("height"), value = 8),
                                        selectInput("dwn_overlap_gate_units", label = h5("units"),
                                                    choices = list("in" = "in", "cm" = "cm", "mm" = "mm"), selected = "in"),
                                        numericInput("dwn_overlap_gate_dpi", label = h5("dpi"), value = 300),
                                        icon = icon("save"), status = "primary", tooltip = tooltipOptions(title = "save plot")
                                      )
                               )
                             )
                    ),
                    tabPanel("Annotations",
                             fluidRow(
                               column(10, plotOutput('ann_overlap_gate')),
                               column(2,
                                      dropdownButton(
                                        tags$h4("Advanced options"),
                                        checkboxInput("legend_ann_overlap_plot_gate", label = "Remove legend", value = TRUE),
                                        icon = icon("gear"), status = "primary", tooltip = tooltipOptions(title = "plot setting")
                                      ),
                                      dropdownButton(
                                        selectInput('dwn_ann_overlap_gate_ext', label = NULL,
                                                    choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png",
                                                                   'tiff' = "tiff", 'svg' = "svg", 'bmp' = "bmp")),
                                        downloadButton('dwn_ann_overlap_gate', ""),
                                        hr(),
                                        numericInput("dwn_ann_overlap_gate_width", label = h5("width"), value = 10),
                                        numericInput("dwn_ann_overlap_gate_height", label = h5("height"), value = 8),
                                        selectInput("dwn_ann_overlap_gate_units", label = h5("units"),
                                                    choices = list("in" = "in", "cm" = "cm", "mm" = "mm"), selected = "in"),
                                        numericInput("dwn_ann_overlap_gate_dpi", label = h5("dpi"), value = 300),
                                        icon = icon("save"), status = "primary", tooltip = tooltipOptions(title = "save plot")
                                      )
                               )
                             )
                             )

                  ),
                  shinydashboard::box(
                    visNetwork::visNetworkOutput('gate_antology_graph')
                  )
                )
        ),

        #############################
        ##### Third tab content #####
        #############################
        shinydashboard::tabItem(tabName = 'data_clustering',
                fluidRow(
                  tabBox(
                    tabPanel("clustering",
                             radioButtons("mode_k_choice", label = h4("Choose of number of clusters"),
                                          choices = list("Automatically detect optimum" = 1, "Manually choose" = 2),
                                          selected = 1),
                             conditionalPanel(
                               condition = "input.mode_k_choice == 1",
                               numericInput("rate_var_explan", label = h4("Rate of explained variance"), value = 0.9),
                               numericInput("maxK", label = h4("Max number of clusters"), value = 20)
                             ),
                             conditionalPanel(
                               condition = "input.mode_k_choice == 2",
                               numericInput("k", label = h4("Choose number of clusters"), value = 8)
                             ),
                             uiOutput("mk_subset_clusters_ui"),
                             actionButton('start_clustering', label = "Clustering")
                    ),
                    tabPanel("cluster management",
                    ),
                    tabPanel("Cells management",

                    ),
                    tabPanel("Save",
                    )
                  ),
                  shinydashboard::box(
                    fluidRow(
                      column(2, actionBttn(inputId = "redraw_clust", style = "material-circle", color = "default" ,icon = icon("redo"))),
                      column(2,
                             dropdownButton(
                               tags$h4("Options of plotting"),
                               numericInput('sampling_size_clust', label = h5("Cell fraction to display"), value = 0.5, step = 0.1),
                               materialSwitch(inputId = 'fuse_clust', label = h4("Size fuse"), value = TRUE),
                               selectInput("method_plot_clust", label = h5("Visualisation method"),
                                           choices = list("tSNE" = "tSNE", "UMAP" = "UMAP"),selected = "tSNE"),
                               numericInput('point_size_clust', label = h5("Size of points"), value = 0.3, step = 0.1),
                               icon = icon("edit"), status = "primary", tooltip = tooltipOptions(title = "plot setting")
                             )
                      ),
                      column(2,
                             dropdownButton(
                               tags$h4("Advanced options"),
                               uiOutput('advanced_opt_clust_ui'),
                               icon = icon("gear"), status = "primary", tooltip = tooltipOptions(title = "plot setting")
                             )

                      ),
                      column(2,
                             dropdownButton(
                               selectInput('dwn_scatter_clust_ext', label = NULL,
                                           choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png",
                                                          'tiff' = "tiff", 'svg' = "svg", 'bmp' = "bmp")),
                               downloadButton('dwn_scatter_clust', ""),
                               hr(),
                               numericInput("dwn_scatter_clust_width", label = h5("width"), value = 10),
                               numericInput("dwn_scatter_clust_height", label = h5("height"), value = 8),
                               selectInput("dwn_scatter_clust_units", label = h5("units"),
                                           choices = list("in" = "in", "cm" = "cm", "mm" = "mm"), selected = "in"),
                               numericInput("dwn_scatter_clust_dpi", label = h5("dpi"), value = 300),
                               icon = icon("save"), status = "primary", tooltip = tooltipOptions(title = "save plot")
                             )

                      )
                    ),
                    plotOutput('scatter_plot_clust'),
                    uiOutput('mk_scatter_clust_ui')
                  )
                ),
                fluidRow(
                  tabBox(
                    tabPanel("Markers",
                             h4("Used for clustering markers"),
                             verbatimTextOutput('mk_clusted_clust'),
                             h4("Used markers"),
                             verbatimTextOutput('mk_rested_clust'),
                             h4("Unused markers"),
                             verbatimTextOutput('mk_excluded_clust')
                    ),
                    tabPanel("Abundance",
                             fluidRow(
                               column(10, plotOutput('abundance_clust')),
                               column(2,
                                      dropdownButton(
                                        selectInput('dwn_abundance_clust_ext', label = NULL,
                                                    choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png",
                                                                   'tiff' = "tiff", 'svg' = "svg", 'bmp' = "bmp")),
                                        downloadButton('dwn_abundance_clust', ""),
                                        hr(),
                                        numericInput("dwn_abundance_clust_width", label = h5("width"), value = 10),
                                        numericInput("dwn_abundance_clust_height", label = h5("height"), value = 8),
                                        selectInput("dwn_abundance_clust_units", label = h5("units"),
                                                    choices = list("in" = "in", "cm" = "cm", "mm" = "mm"), selected = "in"),
                                        numericInput("dwn_abundance_clust_dpi", label = h5("dpi"), value = 300),
                                        icon = icon("save"), status = "primary", tooltip = tooltipOptions(title = "save plot")
                                      )
                               )
                             )
                    ),
                    tabPanel("Clusiering",
                             fluidRow(
                               column(10, plotOutput('clusters_overlap_clust')),



                               column(2,
                                      dropdownButton(
                                        tags$h4("Advanced options"),
                                        checkboxInput("legend_clusters_overlap_plot_clust", label = "Remove legend", value = TRUE),
                                        icon = icon("gear"), status = "primary", tooltip = tooltipOptions(title = "plot setting")
                                      ),
                                      dropdownButton(
                                        selectInput('dwn_clusters_overlap_clust_ext', label = NULL,
                                                    choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png",
                                                                   'tiff' = "tiff", 'svg' = "svg", 'bmp' = "bmp")),
                                        downloadButton('dwn_clusters_overlap_clust', ""),
                                        hr(),
                                        numericInput("dwn_clusters_overlap_clust_width", label = h5("width"), value = 10),
                                        numericInput("dwn_clusters_overlap_clust_height", label = h5("height"), value = 8),
                                        selectInput("dwn_clusters_overlap_clust_units", label = h5("units"),
                                                    choices = list("in" = "in", "cm" = "cm", "mm" = "mm"), selected = "in"),
                                        numericInput("dwn_clusters_overlap_clust_dpi", label = h5("dpi"), value = 300),
                                        icon = icon("save"), status = "primary", tooltip = tooltipOptions(title = "save plot")
                                      )
                               )
                             )
                             )


                  )
                )
        ),

        #############################
        #### Fourth tab content #####
        #############################
        tabItem(tabName = 'data_exploration',
                fluidRow(
                  shinydashboard::box(

                    width = 12
                  ),
                  shinydashboard::box(

                  )
                )
        ),

        # Fifth tab content
        tabItem(tabName = 'data_correlation',
                fluidRow(

                ),
                fluidRow(
                  shinydashboard::box(

                  ),
                  shinydashboard::box(
                    fluidRow(
                      column(2,

                      )
                    ),
                  ),
                  tabBox(
                    tabPanel("Abundance correlations",

                    ),
                    tabPanel("Marker correlations",
                    ),
                    width = 12
                  )
                )
        ),

        #############################
        ##### Sixth tab content #####
        #############################
        tabItem(tabName = 'data_crosspanel',
                fluidRow(
                  tabBox(
                    tabPanel("Upload panel",


                    ),
                    tabPanel("Remove panel",
                    )
                  ),
                  tabBox(
                    tabPanel("Samples content",
                    ),
                    tabPanel("Markers content",
                    )
                  ),
                  tabBox(
                    tabPanel("Abundance cross-panel correlation",
                             fluidRow(

                             ),
                    ),
                    tabPanel("Expressing cell fraction Correlation",
                             fluidRow(

                             ),
                             fluidRow(
                             )

                    ),
                    width = 12
                  )

                )
        )
      )
    )
  )

  shinyApp(ui = cytofBrowser_ui, server = cytofBrowser_server)
}
