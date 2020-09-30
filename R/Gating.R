
#' Take data to gating scatterplot
#'
#' @param fcs_raw flowSet object with data
#' @param gating_subset vector of logic values for expression matrix extracted from flowSet object
#' @param gating_mk1 marker for plotting X axis of density plot for gating
#' @param gating_mk2 marker for plotting Y axis of density plot for gating
#'
#' @return
#' @importClassesFrom flowCore flowSet
#' @importFrom flowCore fsApply exprs "exprs<-"
#'
get_data_for_gating <- function(fcs_raw, gating_subset, gating_mk1, gating_mk2){
  raw_df <- flowCore::fsApply(fcs_raw, function(x) flowCore::exprs(x))
  raw_df <- as.data.frame(raw_df)
  raw_df$original_cell_coordinates <- 1:nrow(raw_df)
  raw_df <- raw_df[gating_subset, c(gating_mk1, gating_mk2, "original_cell_coordinates")]
  colnames(raw_df) <- c("gating_mk1", "gating_mk2", "original_cell_coordinates")
  return(raw_df)
}

### test
#raw_df <- get_data_for_gating(fcs_raw, rep(TRUE, 13402), "I127Di", "Tm169Di")


#' Extract data table for plot to gate
#'
#' @param exprs_data matrix with expression data
#' @param gating_subset subset to plot to gating
#' @param gating_mk1 marker for X axis for gating plot
#' @param gating_mk2 marker for Y axis for gating plot
#'
#' @return
get_exprs_data_for_gating <- function(exprs_data, gating_subset, gating_mk1, gating_mk2){
  raw_df <- as.data.frame(exprs_data[,c(gating_mk1, gating_mk2)])
  raw_df$original_cell_coordinates <- 1:nrow(exprs_data)
  raw_df <- raw_df[gating_subset,]
  colnames(raw_df) <- c("gating_mk1", "gating_mk2", "original_cell_coordinates")
  return(raw_df)
}

#' Check and adjusting sparse data for density plot for gating
#'
#' @param gated_data_subset matrix of data for adjusting
#' @param mk column in matrix for adjusting and plotting
#'
#' @return
#' @importFrom stats quantile
get_modif_sparse_data_gating <- function(gated_data_subset, mk = c("gating_mk1", "gating_mk2")){
  new_gated_data_subset <- gated_data_subset
  for (i in mk){
    data_vector <- new_gated_data_subset[,i]
    qu <- stats::quantile(data_vector)
    print(qu)
    if((qu["50%"] == qu["75%"]) & (qu["25%"] == qu["50%"])){
      n_modif <- (as.integer(length(data_vector)*0.25)+1)-sum(data_vector > qu["50%"])
      adj_modif <- (max(data_vector) - min(data_vector)) *0.001
      mid <- which(data_vector == qu["50%"])
      data_vector[sample(mid,n_modif)] <- qu["50%"] + adj_modif
    }
    new_gated_data_subset[,i] <- data_vector
  }
  return(new_gated_data_subset)
}



#' Function adds set of gates to cell annotation object
#'
#' @param gates data frame with allocated gates with logical values in it
#' @param gate_list name of the gates, column names from gates data frame
#' @param cell_annotation cell annotation data frame to save new annotation from gates
#' @param method methd of converting gates to cell annotation "pure" or "squeeze"
#'
#' @return
get_cell_type_from_gates <- function(gates, gate_list, cell_annotation, method = 'pure'){
  target_gates <- as.matrix(gates[,as.character(unlist(gate_list))])
  populations <- apply(target_gates, 1, function(cell){paste(gate_list[cell], sep = "", collapse = "_&_")})
  populations[populations == ""] <- "untyped"
  new_cell_annotation <- cell_annotation
  if(method == 'pure'){
    new_name <- "pure_gated_cell_type"
  }
  if(method == 'squeeze'){
    purity <- rowSums(target_gates)
    populations[!(purity == 1)] <- "untyped"
    new_name <- "squeeze_gated_cell_type"
  }
  new_cell_annotation$new_gate <- populations
  if(any(grepl(new_name, colnames(new_cell_annotation)))){
    new_name <- paste0(new_name,"_",sum(grepl(new_name, colnames(new_cell_annotation)))+1)}
  colnames(new_cell_annotation)[colnames(new_cell_annotation) == "new_gate"] <- new_name
  return(new_cell_annotation)
}


#' Extract annotation fron fcs_raw object to cell annotation object
#'
#' @param entire_panel wide range of column names from floeSet without any filtration
#' @param fcs_raw flowSet object with data
#' @param cell_annotation cell annotation data frame to save new annotation from gates
#'
#' @return
#' @importClassesFrom flowCore flowSet
#' @importFrom flowCore fsApply exprs "exprs<-"

get_add_cell_annotation_from_data <- function(entire_panel, fcs_raw, cell_annotation){
  raw_df <- flowCore::fsApply(fcs_raw, function(x) flowCore::exprs(x)[,entire_panel])
  raw_df <- as.data.frame(raw_df)
  cell_annotation <- cbind(cell_annotation, raw_df)
  return(cell_annotation)
}


#' Adding of cluster-info to flowSet object
#'
#' @param fcs_raw flowSet object with data
#' @param cell_clustering_list vector of cluster annotation for cells in order of expression data
#'
#' @return
#' @importFrom flowCore sampleNames fr_append_cols
#' @importClassesFrom flowCore flowSet
#'
get_cell_annotation_fcs_files <- function(fcs_raw, cell_annotation, column_names = NULL){
  if(is.null(column_names)){column_names <- colnames(cell_annotation)[-1]}
  if(!all(flowCore::sampleNames(fcs_raw) %in% unique(cell_annotation$samples))){
    print("Cluster and data samples does not match")
    return(NULL)}
  clustered_fcs <- lapply(flowCore::sampleNames(fcs_raw), function(s) {
    addition_col <- as.matrix(cell_annotation[cell_annotation$samples == s, column_names])
    addition_col <- apply(addition_col, 2, function(x) as.integer(as.factor(x)))
    colnames(addition_col) <- column_names
    fs_d <- flowCore::fr_append_cols(fcs_raw[[s]], addition_col)
    return(fs_d)
  })
  names(clustered_fcs) <- flowCore::sampleNames(fcs_raw)
  clustered_fcs <- as(clustered_fcs, 'flowSet')
  return(clustered_fcs)
}

#' Create tidy data frame to order overlaped gates to plotting
#'
#' @param gates dataframe with gates data
#'
#' @return
#' @importFrom reshape2 melt
get_gates_overlap_data <- function(gates){
  plot_gates <- gates
  order_tags <- apply(plot_gates, 1,function(x){paste(x,collapse = "_")})
  plot_gates <- plot_gates[order(order_tags),]
  plot_gates$order_tags <- 1:nrow(plot_gates)
  plot_gates <- reshape2::melt(plot_gates, id = "order_tags")
  colnames(plot_gates) <- c("Cells", "Gates", "value")
  return(plot_gates)
}

#' Create tidy data frame to plot portions of groups from annotations
#'
#' @param ann annotation data frame
#' @param col_names vector with column names to plot
#'
#' @return
#' @importFrom RColorBrewer brewer.pal
#'
get_ann_overlap_data <- function(ann,col_names = NULL){
  if(is.null(col_names)){
    exeptions <- c("all_cells", "tSNE1", "tSNE2", "UMAP1", "UMAP2")
    col_names <- colnames(ann)
    col_names <- col_names[!(col_names %in% exeptions)]
  }
  plot_ann <- do.call(rbind, lapply(col_names, function(x){
    count_data <- table(ann[,x])
    data.frame(annotations = x, group = names(count_data), value = as.numeric(count_data) )
  }))
  pal_order <- c("Dark2", "Set1", "Set2", "Paired", "Accent", "Set3", "Pastel1", "Pastel2")
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  ann_colour <- as.character(unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals[pal_order, 'maxcolors'],
                                               pal_order))[1:nrow(plot_ann)])
  plot_ann$color <- ann_colour
  return(plot_ann)
}
