#' Gating Module for cytofBrowserLite
#'
#' This module encapsulates all server logic related to gating,
#' including gate creation, management, and visualization.
#'
#' @param input Shiny input object.
#' @param output Shiny output object.
#' @param session Shiny session object.
#' @param fcs_data Reactive values object containing main data.
#' @param gates Reactive values object for gating data.
#' @param plots Reactive values object for plots.
#' @param roots File system roots for shinyFiles.
#' 
#' @return A list of reactive elements and observers for the gating module.

gating_module <- function(input, output, session, fcs_data, gates, plots, roots) {
  
  # Define reactive elements and observers specific to gating
  # This will include all the logic from the 'Gating' section of serverLite.R
  
  # Example of how to structure the module content:
  # (This is a simplified example, the actual content will be much larger)
  
  # Create UI to set the gating plot
  output$gating_api_ui <- renderUI({
    if(is.null(gates$gates)){return(NULL)}
    if(is.null(fcs_data$entire_panel)){return(NULL)}
    fluidRow(
      column(1),
      column(10,
             selectInput('gating_subset', label = h5("Data for gating"),
                         choices = colnames(gates$gates), selected = 1),
             selectInput('gating_mk1', label = h5("Marker for X"),
                         choices = names(fcs_data$entire_panel), selected = 1),
             selectInput('gating_mk2', label = h5("Marker for Y"),
                         choices = names(fcs_data$entire_panel), selected = 2)
      )
    )
  })


  # Create UI to convert gates to cell annotation
  output$mergeing_gates_ui <- renderUI({
    if(is.null(gates$gates)){return(NULL)}
    fluidRow(
      column(1),
      column(10,
             h4("Convert gates to cell annotation"),
             selectInput("gate_converting_method", label = h5("Converting method"),
                         choices = list("Pure selection" = 'pure', "Gates squeezing" = 'squeeze'),selected = 1),
             selectInput("gates_to_convert_gating", label = h5("Gates to convert"),
                         choices = gates$antology$name[-1], multiple = TRUE),
             textInput("gates_name_new_ann", label = h5("New annotation name"),value = NULL, placeholder = "annotation_from_gates"),
             actionButton("convert_gates", label = "Convert")
      )
    )
  })


  # Create UI to rename gates
  output$rename_gates_ui <- renderUI({
    if(is.null(input$gated_node_id)){return(NULL)}
    ui <- fluidRow(
      fluidRow(
        column(1),
        column(10,
               h4("Rename gate"),
               textInput('new_gate_name_gating', label = h5("Write new name"),
                         value = as.character(input$gated_node_id)),
               actionButton("rename_gates", label = "Rename")
        )
      ),
      fluidRow(
        column(1),
        column(6, h5(paste0("Delete gate ", input$gated_node_id))),
        column(4, actionButton("delete_gate", label = "Delete"))
      )
    )
    return(ui)
  })


  # UI for make subset by gates
  output$subset_gates_ui <- renderUI({
    fluidRow(
      column(1),
      column(10,
             radioGroupButtons(inputId = "gates_subseting_orient", label = h5("Subsatting by gates"),
                               choices = c("exclude", "include"), selected = "exclude", justified = TRUE),
             selectInput("gates_to_subseting", label = h5("Gates to use "),
                         choices = gates$antology$name[-1], multiple = TRUE),
             actionButton("butt_subset_gating", label = "Subset")
      )
    )
  })

  # Implement subsetting data by gates it all data objects
  observeEvent(input$butt_subset_gating,{
    if(is.null(gates$gates)){return(NULL)}
    if(is.null(input$gates_to_subseting)){return(NULL)}
    if(length(input$gates_to_subseting) == 1){stay_coord <- gates$gates[,input$gates_to_subseting]}
    if(length(input$gates_to_subseting) > 1){stay_coord <- apply(gates$gates[,input$gates_to_subseting], 1, any)}
    if(input$gates_subseting_orient == "exclude"){stay_coord <- !stay_coord}
    fcs_data$exprs_data <- fcs_data$exprs_data[stay_coord,]
    fcs_data$cell_ann <- fcs_data$cell_ann[stay_coord,]
    gates$gates <- gates$gates[stay_coord,]
    fcs_data$subset_coord <- which(stay_coord[as.integer(fcs_data$subset_coord)])
    gates$subset_coord <- which(stay_coord[gates$subset_coord])
  })


  # Make subset of data to gatingplot
  observeEvent(input$butt_plot_for_gating, {
    if(is.null(gates$gates)){return(NULL)}
    if(is.null(input$gating_subset)){return(NULL)}
    if(is.null(input$gating_mk1)){return(NULL)}
    if(is.null(input$gating_mk2)){return(NULL)}
    withProgress(message = "Making of subset data", min =0, max = 3, value = 0,{
      incProgress(1)
      gates$gated_data_subset <- get_exprs_data_for_gating(exprs_data = fcs_data$exprs_data,
                                                     gating_subset = gates$gates[,input$gating_subset],
                                                     gating_mk1 = fcs_data$entire_panel[input$gating_mk1],
                                                     gating_mk2 = fcs_data$entire_panel[input$gating_mk2])
      incProgress(1)
      gates$gated_data_subset <- get_modif_sparse_data_gating(gates$gated_data_subset)
      gates$subset_coord <- c(1:nrow(gates$gated_data_subset))
      if(input$fuse_gate){gates$subset_coord <- get_size_subset_gating_data(gates$gated_data_subset)}
      incProgress(1)
    })
  })


  # Drawing interactive density plot for gating
  output$scatter_plot_gating <- renderPlot({
    if(is.null(gates$gated_data_subset)){return(NULL)}
    withProgress(message = "Density drawing", min =0, max = 3, value = 0,{
      colfunc <- grDevices::colorRampPalette(c("black", "red", "yellow"))
      incProgress(1, detail = "Data preparation" )
      min_mk1 <- min(gates$gated_data_subset$gating_mk1)
      max_mk1 <- max(gates$gated_data_subset$gating_mk1)
      min_mk2 <- min(gates$gated_data_subset$gating_mk2)
      max_mk2 <- max(gates$gated_data_subset$gating_mk2)
      incProgress(1, detail = "Drawing")
      plots$scatter_plot_gating <- ggplot(gates$gated_data_subset[gates$subset_coord,], aes(x = gating_mk1, y = gating_mk2)) +
        ylim((min_mk2 - 0.1*(max_mk2 - min_mk2)),
             (max_mk2 + 0.1*(max_mk2 - min_mk2))) +
        xlim((min_mk1 - 0.1*(max_mk1 - min_mk1)),
             (max_mk1 + 0.1*(max_mk1 - min_mk1))) +
        xlab(input$gating_mk1) +
        ylab(input$gating_mk2) +
        geom_point(size = 0.1) +
        stat_density_2d(aes(fill = ..level..), geom = "polygon") +
        scale_fill_gradientn(colours=colfunc(400))+
        geom_density2d(colour="black") +
        theme_bw() +
        theme(legend.position = "none")
      incProgress(1)
    })
    return(plots$scatter_plot_gating)
  })

  # Download gates density plot
  output$dwn_density_gate <- downloadHandler(
    filename = function() {
      ext <- input$dwn_density_gate_ext
      if(is.null(ext)){ext <- "pdf"}
      paste("Gates_density_plot", ext, sep = ".") },
    content = function(file) {
      ext <- input$dwn_density_gate_ext
      if(is.null(ext)){ext <- "pdf"}
      ggsave(file, plot = plots$scatter_plot_gating, device = ext,
             width = input$dwn_density_gate_width,height = input$dwn_density_gate_height,
             units = input$dwn_density_gate_units,dpi =input$dwn_density_gate_dpi)
    }
  )

  # Drawing the gates overlap plot
  output$gates_overlap <- renderPlot({
    if(is.null(gates$gates)){return(NULL)}
    if(ncol(gates$gates) <= 1){return(NULL)}
    withProgress(message = "Gates overlap drawing", min =0, max = 3, value = 0,{
      incProgress(1, detail = "Overlap finding")
      plot_gates <- get_gates_overlap_data(gates = gates$gates)
      incProgress(1, detail = "Drawing")
      plots$gate_overlap <- ggplot(plot_gates, aes(Gates, Cells, fill= value))+
        geom_tile()+
        scale_fill_manual(values = c("white", "black"))+
        theme_bw()+
        theme(legend.position = "none")+
        theme(axis.text.x = element_text(angle = 90))
      incProgress(1)
    })
    return(plots$gate_overlap)
  })

  # Download gates overlap plot
  output$dwn_overlap_gate <- downloadHandler(
    filename = function() {
      ext <- input$dwn_overlap_gate_ext
      if(is.null(ext)){ext <- "pdf"}
      paste("Gates_overlap_plot", ext, sep = ".") },
    content = function(file) {
      ext <- input$dwn_overlap_gate_ext
      if(is.null(ext)){ext <- "pdf"}
      ggsave(file, plot = plots$gate_overlap, device = ext,
             width = input$dwn_overlap_gate_width,height = input$dwn_overlap_gate_height,
             units = input$dwn_overlap_gate_units,dpi =input$dwn_overlap_gate_dpi)
    }
  )

  # Drawing the annotation overlap plot in gate section
  output$ann_overlap_gate <- renderPlot({
    if(is.null(fcs_data$cell_ann)){return(NULL)}
    if(ncol(fcs_data$cell_ann) <= 1){return(NULL)}
    plot_ann <- get_ann_overlap_data(ann = fcs_data$cell_ann)
    pl <- ggplot(plot_ann, aes(fill=group, y=value, x = annotations)) +
      geom_bar(position="fill", stat="identity")+
      scale_fill_manual(values=plot_ann$color)+
      theme_bw() +
      ylab("")+
      theme(axis.text.x = element_text(angle = 90))
    if(input$legend_ann_overlap_plot_gate){pl <- pl + theme(legend.position="none")}
    plots$ann_overlap_gate <- pl
    return(plots$ann_overlap_gate)
  })

  # Download annotation overlap plot in gate section
  output$dwn_ann_overlap_gate <- downloadHandler(
    filename = function() {
      ext <- input$dwn_ann_overlap_gate_ext
      if(is.null(ext)){ext <- "pdf"}
      paste("Gates_overlap_plot", ext, sep = ".") },
    content = function(file) {
      ext <- input$dwn_ann_overlap_gate_ext
      if(is.null(ext)){ext <- "pdf"}
      ggsave(file, plot = plots$ann_overlap_gate, device = ext,
             width = input$dwn_ann_overlap_gate_width,height = input$dwn_ann_overlap_gate_height,
             units = input$dwn_ann_overlap_gate_units,dpi =input$dwn_ann_overlap_gate_dpi)
    }
  )

  # The object for the anthology graph of gates
  output$gate_antology_graph <- visNetwork::renderVisNetwork({
    if(is.null(gates$antology)){return(NULL)}
    nodes <- data.frame(id = gates$antology$name, label = gates$antology$name)
    edges <- data.frame(from = gates$antology$parent, to = gates$antology$name)
    edges <- edges[!is.na(edges$from),]
    net <- visNetwork::visNetwork(nodes, edges)
    net <- visNetwork::visInteraction(net, hover = TRUE)
    net <- visNetwork::visEvents(net, select = "function(nodes) { Shiny.onInputChange('gated_node_id', nodes.nodes);}")
    visNetwork::visHierarchicalLayout(net)
  })

  # Allocation a new gate
  observeEvent(input$gete_chosen_cells, {
    if(is.null(input$brush_gating)){return(NULL)}
    if(is.null(gates$gated_data_subset)){return(NULL)}
    new_name <- paste0("Gate", as.character(ncol(gates$gates)+1), "_",
                       input$gating_mk1, "_", input$gating_mk2, "_from_", strsplit(input$gating_subset, "_")[[1]][1])
    if(!is.null(input$new_gate_name) & (input$new_gate_name != "marker1+/marker2+")){new_name <- input$new_gate_name}
    original_cell_coordinates <- brushedPoints(gates$gated_data_subset, input$brush_gating, xvar = "gating_mk1", yvar = "gating_mk2")
    original_cell_coordinates <- original_cell_coordinates$original_cell_coordinates
    gates$gates$new_gate <- FALSE
    gates$gates[original_cell_coordinates, "new_gate"] <- TRUE
    if(any(grepl(new_name, colnames(gates$gates)))){new_name <- paste0(new_name,"_",sum(grepl(new_name, colnames(gates$gates)))+1)}
    colnames(gates$gates)[which(colnames(gates$gates) ==  "new_gate")] <- new_name
    ## Add anthology note of gating for graph
    gates$antology <- rbind(gates$antology, data.frame(name = new_name, parent = input$gating_subset))
    rownames(gates$antology) <- gates$antology$name
  })

  # Converting gates to cell annotation data
  observeEvent(input$convert_gates, {
    if(is.null(input$gate_converting_method)){return(NULL)}
    fcs_data$cell_ann <- get_cell_type_from_gates(gates$gates, input$gates_to_convert_gating,
                                            fcs_data$cell_ann, method = input$gate_converting_method,
                                            name_new_ann = input$gates_name_new_ann)
  })

  # Renew gate reactive object after gate rename
  observeEvent(input$rename_gates, {
    if(is.null(input$new_gate_name_gating)){return(NULL)}
    if(input$new_gate_name_gating == ""){return(NULL)}
    gates$antology$name[gates$antology$name == input$gated_node_id] <- input$new_gate_name_gating
    gates$antology$parent[gates$antology$parent == input$gated_node_id] <- input$new_gate_name_gating
    rownames(gates$antology) <- gates$antology$name
    colnames(gates$gates)[colnames(gates$gates) == input$gated_node_id] <- input$new_gate_name_gating
  })

  # Delete gate from gate reactive object
  observeEvent(input$delete_gate, {
    if(is.null(input$gated_node_id)){return(NULL)}
    gates$antology$parent[gates$antology$parent == input$gated_node_id] <- gates$antology[input$gated_node_id, "parent"]
    gates$antology <- gates$antology[rownames(gates$antology) != input$gated_node_id,]
    gates$gates[,colnames(gates$gates) == input$gated_node_id] <- NULL
  })
  
  # Return a list of reactive elements if needed
  # For now, this module primarily sets up observers and outputs
  return(list())
}