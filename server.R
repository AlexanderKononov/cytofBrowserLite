library(shiny)

# Source all necessary files
source("R/serverLite.R")
source("R/Data_preparation.R")
source("R/Gating.R")
source("R/Clustering.R")
source("R/Enrichments.R")
source("R/Cross_panel_analysis.R")

# Also source the modules
module_files <- list.files("R/modules", pattern = "\\.R$", full.names = TRUE)
for (file in module_files) {
  source(file)
}

# Load required libraries
library(shinydashboard)
library(shinyWidgets)
library(flowCore)
library(Rtsne)
library(umap)
library(ggplot2)
library(shinyFiles)
library(grDevices)
library(visNetwork)
library(FlowSOM)
library(ConsensusClusterPlus)
library(RColorBrewer)
library(reshape2)
library(ComplexHeatmap)

function(input, output) {
  cytofBrowser_server(input, output)
}