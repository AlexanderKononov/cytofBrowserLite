#' Clustering Module for cytofBrowserLite
#'
#' This module encapsulates all server logic related to clustering,
#' including cluster creation, management, and visualization.
#'
#' @param input Shiny input object.
#' @param output Shiny output object.
#' @param session Shiny session object.
#' @param fcs_data Reactive values object containing main data.
#' @param clusters Reactive values object for clustering data.
#' @param data_prep_settings Reactive values object for data preparation settings.
#' @param plots Reactive values object for plots.
#' @param roots File system roots for shinyFiles.
#' 
#' @return A list of reactive elements and observers for the clustering module.

clustering_module <- function(input, output, session, fcs_data, clusters, data_prep_settings, plots, roots) {
  
  # Define reactive elements and observers specific to clustering
  # This will include all the logic from the 'Clustering' section of serverLite.R
  
  # Example of how to structure the module content:
  # (This is a simplified example, the actual content will be much larger)
  
  # Create UI to choose excluded markers from clusterisation
  output$mk_subset_clusters_ui <- renderUI({
    if(is.null(fcs_data$use_markers)){return(NULL)}
    selectInput('exclude_mk_clusters', label = "Exclude markers from clustering",
                choices = names(fcs_data$use_markers), multiple = TRUE)
  })

  # Action on the Clustering button
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

  # Update reactive objects euclid_dist
  observe({
    if(is.null(fcs_data$cell_ann$clusters)){return(NULL)}
    withProgress(message = "Cluster estimation", min =0, max = 2, value = 0,{
      incProgress(1, detail = "distance between clusters")
      clusters$clus_euclid_dist <- get_euclid_dist(exprs_data = isolate(fcs_data$exprs_data), use_markers = isolate(fcs_data$use_markers),
                                                   cell_clustering = fcs_data$cell_ann$clusters)
      incProgress(1, detail = "graph elements")
    })
  })

  # Update reactive objects with nodes and edges
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

  # Create UI to choose clusters to merge
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

  # Renew clustering reactive object after merging
  observeEvent(input$merge_clust, {
    fcs_data$cell_ann$clusters <- cluster_merging(fcs_data$cell_ann$clusters, input$cluster_to_merge_clust)
  })

  # Create UI to rename clusters
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

  # Renew cluster reactive object after cluster rename
  observeEvent(input$rename_clust, {
    if(is.null(input$new_cluster_name_clust)){return(NULL)}
    if(input$new_cluster_name_clust == ""){return(NULL)}
    fcs_data$cell_ann$clusters[fcs_data$cell_ann$clusters==input$current_node_id] <- input$new_cluster_name_clust
  })

  # UI to choose marker for scatter plot with clusters
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

  # UI for advanced options for clustering page
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

  # Drawing the reactive and interactive UMAP plot
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

  # Rewrite for scatter plot in clusteringpage
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

  # Download scatter plot data preparation
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

  # Reactive show current use_markers from "fcs_data" object
  output$mk_clusted_clust <- renderPrint({clusters$clust_markers})
  output$mk_rested_clust <- renderPrint({fcs_data$use_markers})
  output$mk_excluded_clust <- renderPrint({
    fcs_data$entire_panel[!(fcs_data$entire_panel %in% fcs_data$use_markers)]
  })

  # Drawing the reactive abundance plot
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

  # Download abundance plot for clustering
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

  # Drawing the clusters overlap plot in cluster section
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

  # Download clusters overlap plot in cluster section
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

  # Drawing the reactive and interactive graph with clusters
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
  
  # Return a list of reactive elements if needed
  # For now, this module primarily sets up observers and outputs
  return(list())
}