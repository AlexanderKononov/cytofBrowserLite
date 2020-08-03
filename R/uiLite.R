#' The main function to run cutofCore graphic interface
#' @description Run graphical user interface for cytofBrowser.
#' The function runs the Shiny App in browser. The current package
#' was created with assumption that it will be run with GUI mode
#' as Shiny App, That is why the function is the most appropriate
#' way to use the cytofBrowser.
#'
#' @return The function runs the Shiny App in browser.
#' @import shiny shinydashboard
#' @export
#'
#' @examples
#' \dontrun{
#' cytofBrowserGUI()
#' }
cytofBrowserGUI <-function(){
  cytofBrowser_ui <- dashboardPage(
    dashboardHeader(title = "cytofBrowser"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Data", tabName = 'data_processing', icon = icon("bars")),
        menuItem("Gating", tabName = 'gating', icon = icon("door-open")),
        menuItem("Exploration", tabName = 'data_exploration', icon = icon("chart-area")),
        menuItem("Correlation", tabName = 'data_correlation', icon = icon("braille")),
        menuItem("Cross-panel", tabName = 'data_crosspanel', icon = icon("clone"))
      )
    ),
    dashboardBody(
      tabItems(
        # First tab content
        tabItem(tabName = 'data_processing',
                fluidRow(
                  tabBox(
                    tabPanel("Uploading",

                    ),
                    tabPanel("Transforming",

                    ),
                    tabPanel("Markers",

                    ),
                    tabPanel("Clustering",

                    ),
                    tabPanel("Cluster management",

                    ),
                    tabPanel("Save",

                    )
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
                    ),
                    tabPanel("Gate management",

                    ),
                    tabPanel("Cells management",

                    ),
                    tabPanel("Save",
                    )
                  ),
                  shinydashboard::box(

                  )
                ),
                fluidRow(
                  shinydashboard::box(

                  ),
                  shinydashboard::box(

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
