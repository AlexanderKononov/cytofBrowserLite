#' The main function to run cutofCore graphic interface
#' @description Run graphical user interface for cytofBrowser.
#' The function runs the Shiny App in browser. The current package
#' was created with assumption that it will be run with GUI mode
#' as Shiny App, That is why the function is the most appropriate
#' way to use the cytofBrowser.
#'
#' @return The function runs the Shiny App in browser.
#' @import shiny shinyFiles shinydashboard shinyWidgets
#' @importFrom visNetwork visNetworkOutput
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
        shinydashboard::menuItem("Annotations", tabName = 'data_annotations', icon = icon("buromobelexperte")),
        shinydashboard::menuItem("Enrichments", tabName = 'data_enrichments', icon = icon("chart-area")),
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
                    shiny::tabPanel("Upload",
                                    uiOutput('selected_fcs_dp_ui'),
                                    fluidRow(
                                      column(1),
                                      column(5,shinyFiles::shinyFilesButton('choose_fcs_dp', label= 'Select FCS files',
                                                                            title='Please select FCS files', multiple=TRUE)),
                                      column(5, actionButton('butt_upload_dproc', label = "Upload files"))
                                    ),
                                    hr(),
                                    materialSwitch(inputId = 'test_data_upload_dproc', label = h4("build-in dataset"), value = FALSE),
                                    conditionalPanel(
                                      condition = "input.test_data_upload_dproc == true",
                                      selectInput('test_data_dproc', label = NULL, choices = c("Test data" = 'test_data'), selected = "Test data")
                                    ),
                                    hr(),
                                    materialSwitch(inputId = 'extr_ann_dp', label = h5("Extract annotations"), value = TRUE)

                    ),
                    shiny::tabPanel("Transform",
                                    h4("Implemented transformations"),
                                    verbatimTextOutput('trans_dp'),
                                    hr(),
                                    checkboxGroupInput("transformation_list", label = h4("Transformations"),
                                                       choices = list("asinh" = 'asinh',
                                                                      "log-transfomation" = 'log',
                                                                      "z-score" = 'z-score',
                                                                      "cytofAsinh" = 'cytofAsinh',
                                                                      "cytofLog10" = 'cytofLog10',
                                                                      "outlier squeezing" = 'outlier_squeezing'),
                                                       selected = c("asinh")),
                                    conditionalPanel(
                                      condition = "input.transformation_list.includes('asinh')",
                                      numericInput("cofactor", label = h4("Cofactor for dividing"), value = 5)
                                    ),
                                    conditionalPanel(
                                      condition = "input.transformation_list.includes('outlier_squeezing')",
                                      numericInput("quantile", label = h4("Quantile for outlier detection"), value = 0.01)
                                    ),
                                    hr(),
                                    actionButton('butt_trans_dproc', label = "Transform")
                    ),
                    shiny::tabPanel("Markers",
                                    uiOutput('mk_subset_dp_ui'),
                                    uiOutput('extr_ann_manage_dp_ui')
                    ),
                    shiny::tabPanel("Subset",
                                    uiOutput("subset_dp_1ui"),
                                    uiOutput("subset_dp_2ui")
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
                    ),
                    tabPanel("Gate subseting",
                             uiOutput("subset_gates_ui")
                    )
                  ),
                  shinydashboard::box(
                    fluidRow(
                      column(7),
                      column(2,
                             dropdownButton(
                               tags$h4("Options of plotting"),
                               materialSwitch(inputId = 'fuse_gate', label = h4("Size fuse"), value = TRUE),
                               icon = icon("edit"), status = "primary", tooltip = tooltipOptions(title = "plot setting")
                             )
                      ),
                      column(2,
                             dropdownButton(
                               selectInput('dwn_density_gate_ext', label = NULL,
                                           choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png",
                                                          'tiff' = "tiff", 'svg' = "svg", 'bmp' = "bmp")),
                               downloadButton('dwn_density_gate', ""),
                               hr(),
                               numericInput("dwn_density_gate_width", label = h5("width"), value = 10),
                               numericInput("dwn_density_gate_height", label = h5("height"), value = 8),
                               selectInput("dwn_density_gate_units", label = h5("units"),
                                           choices = list("in" = "in", "cm" = "cm", "mm" = "mm"), selected = "in"),
                               numericInput("dwn_density_gate_dpi", label = h5("dpi"), value = 300),
                               icon = icon("save"), status = "primary", tooltip = tooltipOptions(title = "save plot")
                             )
                      )
                    ),
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
                    tabPanel("Clustering",
                             radioButtons("mode_k_choice", label = h5("Choose of number of clusters"),
                                          choices = list("Automatically detect optimum" = 1, "Manually choose" = 2),
                                          selected = 1),
                             conditionalPanel(
                               condition = "input.mode_k_choice == 1",
                               numericInput("rate_var_explan", label = h5("Rate of explained variance"), value = 0.9),
                               numericInput("maxK", label = h5("Max number of clusters"), value = 20)
                             ),
                             conditionalPanel(
                               condition = "input.mode_k_choice == 2",
                               numericInput("k", label = h5("Choose number of clusters"), value = 8)
                             ),
                             uiOutput("mk_subset_clusters_ui"),
                             actionButton('start_clustering', label = "Clustering")
                    ),
                    tabPanel("Cluster management",
                             uiOutput("mergeing_clust_ui"),
                             hr(),
                             uiOutput("rename_clust_ui")
                    )
                  ),
                  shinydashboard::box(
                    fluidRow(
                      column(2, actionBttn(inputId = "redraw_clust", style = "material-circle", color = "default" ,icon = icon("redo"))),
                      column(2,
                             dropdownButton(
                               tags$h4("Options of plotting"),
                               numericInput('sampling_size_clust', label = h5("Cell fraction to display"), value = 0.5, step = 0.1),
                               materialSwitch(inputId = 'fuse_clust', label = h5("Size fuse"), value = TRUE),
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
                  ),
                  box(
                    visNetwork::visNetworkOutput("clusters_graph"),
                    sliderInput('edges_threshold_clusterisation', "Edge weight threshold for graph",
                                min =0, max = 1, value = 0.5, step = 0.01),
                    sliderInput('gravity_clusterisation', "Gravity for graph",
                                min = -100, max = 0, value = -40, step = 1)
                  )
                )
        ),

        ##################################
        ##### Fourth tab content     #####
        ##################################
        tabItem(tabName = 'data_annotations',
                shinydashboard::box(
                  fluidRow(
                    column(1),
                    column(11,
                           uiOutput('chose_managment_ann_ui'),
                           uiOutput('managment_ann_ui'),
                           actionButton('rename_groups', label = "Rename")
                           ),                  ),
                  width = 3
                ),
                shinydashboard::box(
                  fluidRow(
                    column(3, numericInput('num_plot_ann', label = h5("Number of plots"), value = 4, step = 1)),
                    column(2, actionBttn(inputId = "redraw_ann", style = "material-circle", color = "default" ,icon = icon("redo"))),
                    column(2,
                           dropdownButton(
                             tags$h5("Options of plotting"),
                             numericInput('sampling_size_ann', label = h5("Cell fraction to display"), value = 0.5, step = 0.1),
                             materialSwitch(inputId = 'fuse_ann', label = h5("Size fuse"), value = TRUE),
                             selectInput("method_plot_ann", label = h5("Visualisation method"),
                                         choices = list("tSNE" = "tSNE", "UMAP" = "UMAP"),selected = "tSNE"),
                             numericInput('point_size_ann', label = h5("Size of points"), value = 0.3, step = 0.1),
                             icon = icon("edit"), status = "primary", tooltip = tooltipOptions(title = "plot setting")
                           )
                    ),
                    column(2,
                           dropdownButton(
                             tags$h4("Advanced options"),
                             uiOutput('advanced_opt_ann_ui'),
                             icon = icon("gear"), status = "primary", tooltip = tooltipOptions(title = "plot setting")
                           )
                    )
                  ),
                  fluidRow(
                    uiOutput('plot_set_ann_ui')
                  ),
                  width = 9
                ),
                ),

        #############################
        #### Fifth tab content  #####
        #############################
        tabItem(tabName = 'data_enrichments',
                fluidRow(
                  shinydashboard::box(
                    fluidRow(
                      fluidRow(
                        column(1),
                        column(1, actionBttn(inputId = 'redraw_expression', style = "material-circle", color = "default" ,icon = icon("redo"))),
                        column(3, uiOutput('hm_rows_enrich_ui')),
                        column(1,
                               dropdownButton(
                                 selectInput('method_summarize_enrich', label = h5("Summarise method"),
                                             choices = list('median' = "median", 'mean' = "mean"),
                                             selected = 'median'),
                                 selectInput('row_scaling_enrich', label = h5("Rows scaling method"),
                                             choices = list("Z-score Standardization" = 'zscore',
                                                            "Min-max Normalization" = 'norm', "No scaling" = 'none'),
                                             selected = 'zscore'),
                                 icon = icon("edit"), status = "primary", tooltip = tooltipOptions(title = "plot setting")
                               )
                        ),
                        column(1,
                               dropdownButton(
                                 checkboxInput("hm_enrich_row_dend", label = "Rows clustering", value = TRUE),
                                 checkboxInput("hm_enrich_col_dend", label = "Columns clustering", value = TRUE),
                                 selectInput("hm_enrich_palette", label = h5("Palette"),
                                             choices = list('PiYG'="PiYG", 'PRGn'="PRGn", 'PuOr'="PuOr", 'RdBu'="RdBu", 'RdGy'="RdGy",
                                                            'RdYlBu'="RdYlBu", 'RdYlGn'="RdYlGn", 'Spectral'="Spectral",
                                                            'Blues'="Blues", 'BuGn'="BuGn", 'BuPu'="BuPu", 'GnBu'="GnBu", 'Greens'="Greens",
                                                            'Greys'="Greys", 'Oranges'="Oranges",'OrRd'="OrRd", 'PuBu'="PuBu",
                                                            'PuBuGn'="PuBuGn", 'PuRd'="PuRd", 'Purples'="Purples", 'RdPu'="RdPu",
                                                            'Reds'="Reds", 'YlGn'="YlGn", 'YlGnBu'="YlGnBu", 'YlOrBr'="YlOrBr", 'YlOrRd'="YlOrRd"),
                                             selected = "RdYlBu"),
                                 checkboxInput("hm_enrich_rev_pal", label = "Reverse palette", value = TRUE),
                                 checkboxInput("hm_enrich_show_leg", label = "Show legend", value = TRUE),

                                 icon = icon("gear"), status = "primary", tooltip = tooltipOptions(title = "plot setting")
                               )
                               ),
                        column(1,
                               dropdownButton(
                                 selectInput('dwn_hm_enrich_ext', label = NULL,
                                             choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png",
                                                            'tiff' = "tiff", 'svg' = "svg", 'bmp' = "bmp")),
                                 downloadButton('dwn_hm_enrich', ""),
                                 hr(),
                                 numericInput("dwn_hm_enrich_width", label = h5("width"), value = 10),
                                 numericInput("dwn_hm_enrich_height", label = h5("height"), value = 8),
                                 selectInput("dwn_hm_enrich_units", label = h5("units"),
                                             choices = list("in" = "in", "cm" = "cm", "mm" = "mm"), selected = "in"),
                                 numericInput("dwn_hm_enrich_dpi", label = h5("dpi"), value = 300),
                                 icon = icon("save"), status = "primary", tooltip = tooltipOptions(title = "save plot")
                               )
                        )
                      )
                    ),
                    fluidRow(column(1),
                             column(11, plotOutput('hm_enrich'))
                             ),
                    width = 12
                  ),
                  shinydashboard::box(
                    fluidRow(
                      column(1),
                      column(4, uiOutput('mk_deconvol_enrich_ui')),
                      column(4, uiOutput('gp_deconvol_enrich_ui')),
                      column(1,
                             dropdownButton(
                               selectInput('method_sum_deconvol_enrich', label = h5("Summarise method"),
                                           choices = list('median' = "median", 'mean' = "mean"),
                                           selected = 'median'),
                               icon = icon("edit"), status = "primary", tooltip = tooltipOptions(title = "plot setting")
                             )
                             ),
                      column(1,
                             dropdownButton(
                               selectInput('dwn_deconvol_enrich_ext', label = NULL,
                                           choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png",
                                                          'tiff' = "tiff", 'svg' = "svg", 'bmp' = "bmp")),
                               downloadButton('dwn_deconvol_enrich', ""),
                               hr(),
                               numericInput("dwn_deconvol_enrich_width", label = h5("width"), value = 10),
                               numericInput("dwn_deconvol_enrich_height", label = h5("height"), value = 8),
                               selectInput("dwn_deconvol_enrich_units", label = h5("units"),
                                           choices = list("in" = "in", "cm" = "cm", "mm" = "mm"), selected = "in"),
                               numericInput("dwn_deconvol_enrich_dpi", label = h5("dpi"), value = 300),
                               icon = icon("save"), status = "primary", tooltip = tooltipOptions(title = "save plot")
                             )
                      )
                    ),
                    plotOutput('deconvol_enrich'),
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
                             h4("FCS files"),
                             h5("Choose FSC files to upload them as one panel"),
                             shinyFilesButton('choose_fcs_cross_p', label= h5("Select FCS files"),
                                              title='Please select clustered FCS files to upload as panel', multiple=TRUE),
                             hr(),
                             h4("Name the panel"),
                             textInput("panel_name_cross_p", label = h5("You can name the new panel (try to use a short name or even one or few letters)")),
                             hr(),
                             textInput("extr_clust_pattern_cross_p",
                                       label = h5("full or part column name with clusters info (for cytofBrowser and cytofkit : <cluster>)"), value = "cluster"),
                             hr(),
                             actionButton('butt_upload_cross_p', label = "Upload")

                    ),
                    tabPanel("Remove panel",
                             uiOutput('remove_panel_cross_p_ui')
                    )
                  ),
                  tabBox(
                    tabPanel("Samples content",
                             plotOutput('content_cross_p')
                    ),
                    tabPanel("Markers content",
                             plotOutput('mk_content_cross_p')
                    )
                  ),
                  tabBox(
                    tabPanel("Abundance cross-panel correlation",
                             fluidRow(
                               column(2, actionBttn(inputId = "abund_corr_cross_p", style = "material-circle", color = "default" ,icon = icon("play"))),
                               column(2,
                                      dropdownButton(
                                        tags$h4("Advanced options"),
                                        selectInput('abund_corr_method_cross_p', "Select correlation method:", c("spearman" = 'spearman', "pearson" = 'pearson')),
                                        selectInput('abund_corr_p_mode_cross_p', "p-value", c("p-value" = "pval", "adjusted p-value" = "padj"), selected = 'padj'),
                                        uiOutput('abund_corr_padj_method_cross_p_ui'),
                                        icon = icon("gear"), status = "primary", tooltip = tooltipOptions(title = "Correlation test options")
                                      )
                               ),
                               column(2,
                                      dropdownButton(
                                        selectInput('dwn_abund_corr_cross_p_ext', label = NULL,
                                                    choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png")),
                                        downloadButton('dwn_abund_corr_cross_p', ""),
                                        hr(),
                                        selectInput('dwn_table_abund_corr_cross_p_ext', label = NULL,
                                                    choices = list("abundance data" = 'abund_data', "correlation data" = 'corr_data'),
                                                    selected = 'corr_data'),
                                        downloadButton('dwn_table_abund_corr_cross_p', ""),
                                        icon = icon("save"), status = "primary", tooltip = tooltipOptions(title = "save plot")
                                      )
                               )
                             ),
                             plotOutput('plot_abund_corr_cross_p')
                    ),
                    tabPanel("Expressing cell fraction Correlation",
                             fluidRow(
                               column(2, actionBttn(inputId = "exp_cell_f_corr_cross_p", style = "material-circle", color = "default" ,icon = icon("play"))),
                               column(4, uiOutput('mk_exp_cell_f_cross_p_ui')),

                               column(2,
                                      dropdownButton(
                                        tags$h4("Settings"),
                                        selectInput("exp_cell_f_corr_method", label =  "Choose the method to get expression cells fraction",
                                                    choices =c("clustering" = 'clustering', "threshold" = 'threshold'), selected = 'clustering'),
                                        uiOutput('threshold_exp_cell_f_cross_p_ui'),
                                        icon = icon("edit"), status = "primary", tooltip = tooltipOptions(title = "Cell fraction allocation methods")
                                      )
                               ),
                               column(2,
                                      dropdownButton(
                                        tags$h4("Advanced options"),
                                        selectInput('exp_cell_f_corr_method_cross_p', "Select correlation method:", c("spearman" = 'spearman', "pearson" = 'pearson')),
                                        selectInput('exp_cell_f_corr_p_mode_cross_p', "p-value", c("p-value" = "pval", "adjusted p-value" = "padj"), selected = 'padj'),
                                        uiOutput('exp_cell_f_corr_padj_method_cross_p_ui'),
                                        numericInput('exp_cell_f_corr_min_cell_cross_p', "Min expressing cell number ", min = 0, value = 5, step = 1),
                                        icon = icon("gear"), status = "primary", tooltip = tooltipOptions(title = "Correlation test options")
                                      )

                               ),
                               column(2,
                                      dropdownButton(
                                        selectInput('dwn_exp_cell_f_corr_cross_p_ext', label = NULL,
                                                    choices = list('pdf' = "pdf", 'jpeg' = "jpeg", 'png' = "png")),
                                        downloadButton('dwn_exp_cell_f_corr_cross_p', ""),
                                        hr(),
                                        selectInput('dwn_table_exp_cell_f_corr_cross_p_ext', label = NULL,
                                                    choices = list("expressing cell fractions" = 'exp_cell_f_data', "correlation data" = 'corr_data'),
                                                    selected = 'corr_data'),
                                        downloadButton('dwn_table_exp_cell_f_corr_cross_p', ""),
                                        icon = icon("save"), status = "primary", tooltip = tooltipOptions(title = "save plot")
                                      )
                               )
                             ),
                             fluidRow(
                               column(7, plotOutput('plot_exp_cell_f_corr_cross_p')),
                               column(5, plotOutput('mk_density_plot_cross_p'))
                             )

                    ),
                    width = 12
                  )

                )
        )
        ################
      )
    )
  )

  shinyApp(ui = cytofBrowser_ui, server = cytofBrowser_server)
}
