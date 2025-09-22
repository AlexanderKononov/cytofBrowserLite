# Test script to verify the application can be loaded
library(shiny)

# Test UI
source("ui.R")
cat("UI loaded successfully\n")

# Test server
source("server.R")
cat("Server loaded successfully\n")

cat("Application dependencies loaded successfully!\n")