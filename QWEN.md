# cytofBrowserLite Project Context

## Project Overview

This directory contains the `cytofBrowserLite` R package. It is a lite version of `cytofBrowser`, designed for the analysis and visualization of Cytometry by Time-Of-Flight (CyTOF) proteomic data. The package provides a streamlined pipeline from raw FCS files to downstream analysis including cell population identification, clustering, and correlation analysis. It is primarily intended to be used via its R Shiny graphical user interface (GUI).

Key features:
- Starts analysis directly from FCS files.
- Includes data preprocessing and transformation.
- Facilitates cell population identification using markers and clustering.
- Provides tools for visualizing and interacting with data at each step.
- Features a user-friendly Shiny GUI as the main interaction point.

## Technologies and Architecture

- **Language:** R
- **Framework:** R Shiny (for the GUI)
- **Key Dependencies (from DESCRIPTION):**
  - `shiny`, `shinydashboard`, `shinyWidgets`, `shinyFiles`: For the graphical interface.
  - `flowCore`: For handling FCS files.
  - `Rtsne`, `umap`: For dimensionality reduction.
  - `ggplot2`, `visNetwork`, `ComplexHeatmap`: For data visualization.
  - `FlowSOM`, `ConsensusClusterPlus`: For clustering algorithms.
- **Architecture:**
  - The Shiny app is defined in `R/uiLite.R` (UI) and `R/serverLite.R` (Server logic).
  - Core analysis functions are modularized into separate files in the `R/` directory (e.g., `Data_preparation.R`, `Gating.R`, `Clustering.R`, `Enrichments.R`).
  - The main entry point is the `cytofBrowserGUI()` function.

## Development and Running

### Installation

To install the development version:

1. Ensure you have `BiocManager` and `devtools` installed:
   ```R
   if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   if (!requireNamespace("devtools", quietly = TRUE))
       install.packages("devtools")
   ```
2. Install `cytofBrowserLite` from GitHub:
   ```R
   devtools::install_github("AlexanderKononov/cytofBrowserLite")
   ```

### Running the Application

The application is designed to run as an R Shiny app. Execute the following command in your R console:

```R
cytofBrowser::cytofBrowserGUI()
```

Alternatively, for local development, you can run the script `launch_file_for_tests.R` which sources the necessary files and starts the app.

### Dependencies

The list of required R packages is specified in the `DESCRIPTION` file under the `Imports` section. The `launch_file_for_tests.R` script also lists these libraries explicitly.

### Troubleshooting

Common installation issues can often be resolved by:
- Checking and manually installing dependencies listed in `DESCRIPTION`.
- Updating R and RStudio if compatibility issues arise.
- Manually updating or downgrading the `shinyFiles` package if GUI elements become unresponsive.
