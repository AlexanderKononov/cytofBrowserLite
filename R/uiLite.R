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
        shinydashboard::menuItem("Gating", tabName = 'gating', icon = icon("door-open")),
        shinydashboard::menuItem("Exploration", tabName = 'data_exploration', icon = icon("chart-area")),
        shinydashboard::menuItem("Correlation", tabName = 'data_correlation', icon = icon("braille")),
        shinydashboard::menuItem("Cross-panel", tabName = 'data_crosspanel', icon = icon("clone"))
      )
    ),
    shinydashboard::dashboardBody(
      shinydashboard::tabItems(
        # First tab content
        shinydashboard::tabItem(tabName = 'data_processing',
                shiny::fluidRow(
                  shinydashboard::tabBox(
                    shiny::tabPanel("Uploading",
                                    shinyFiles::shinyFilesButton('choose_fcs_dp', label= 'Select FCS files', title='Please select FCS files', multiple=TRUE),
                                    hr(),
                                    conditionalPanel(
                                      condition = "input.extr_clust_dproc == true",
                                      textInput("extr_clust_pattern_dproc",
                                                label = h5("full or part column name with clusters info (for cytofBrowser and cytofkit : <cluster>)"),
                                                value = "cluster")),
                                    hr(),
                                    actionButton('butt_upload_dproc', label = "Upload")
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

                    ),
                    shiny::tabPanel("Clustering",

                    ),
                    shiny::tabPanel("Cluster management",

                    ),
                    shiny::tabPanel("Save",
                                    h5("Saving the data as FCS files "),
                                    shinyDirButton('choose_panel_clust', label = "Folder choose", title = "Folder for saving FCS files"),
                                    hr(),
                                    h5("Saved set of files can be used as panel data in a cross-panel analysis"),
                                    textInput("name_prefix", label = h5("Name prefix")),
                                    hr(),
                                    uiOutput('save_cell_ann_ui'),
                                    actionButton('dwn_panel_clust', label = "Save")
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
                               icon = icon("save"), status = "primary", tooltip = tooltipOptions(title = "save plot")
                             )

                      )
                    ),
                    plotOutput('scatter_plot_dp'),
                    uiOutput('mk_scatter_dp_ui')
                  ),
                  tabBox(
                    tabPanel("Cell number",

                    ),
                    tabPanel("Expression",

                    ),
                    tabPanel("Markers",

                    ),
                    tabPanel("Abundance",
                    )
                  ),
                )
        ),

        # Second tab content
        tabItem(tabName = 'gating',
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
                    ),
                    tabPanel("Cells management",

                    ),
                    tabPanel("Save",
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
                  shinydashboard::box(

                  ),
                  shinydashboard::box(
                    visNetwork::visNetworkOutput('gate_antology_graph')
                  )
                )
        ),

        # Third tab content
        tabItem(tabName = 'data_exploration',
                fluidRow(
                  shinydashboard::box(

                    width = 12
                  ),
                  shinydashboard::box(

                  )
                )
        ),

        # Fourth tab content
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

        # Fifth tab content
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
