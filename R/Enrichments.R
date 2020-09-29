

##### Get dataframe for plotting heatmap
#' Get dataframe for plotting heatmap
#'
#' @param exprs_data matrix with enrichment data
#' @param use_markers list of used markers
#' @param groups vector with grouped info for each cell (clusters, samples or other)
#' @param summarize_mode method to calculate enrichment of cells groups(median or mean)
#'
#' @return
#' @importFrom stats median
#'
get_hm_enrichments <- function(exprs_data, use_markers, groups, summarize_mode = 'median'){

  exprs_group <- lapply(unique(groups), function(x){exprs_data[groups == x,use_markers]})
  names(exprs_group) <- unique(groups)
  if(summarize_mode == 'median'){
    hm_enrich_data <- lapply(exprs_group, function(x){apply(x, 2, stats::median)})
  }
  if(summarize_mode == 'mean'){
    hm_enrich_data <- lapply(exprs_group, colMeans)
  }
  hm_enrich_data <- as.data.frame(hm_enrich_data)
  names(hm_enrich_data) <- names(exprs_group)
  hm_enrich_data <- t(hm_enrich_data)
  hm_enrich_data <- hm_enrich_data[,use_markers]
  colnames(hm_enrich_data) <- names(use_markers)
  return(hm_enrich_data)
}

#' Scalin data table for enrichments heatmap
#'
#' @param hm_enrich_data table with enrichments data
#' @param method method to scale data ("zscore", "norm" or "none")
#'
#' @return
get_scale_hm_enrichments <- function(hm_enrich_data, method = 'zscore'){
  if(method == 'none'){return(hm_enrich_data)}
  if(method == 'zscore'){scaled_hm_enrich <- apply(hm_enrich_data, 2, function(x){(x-mean(x))/sd(x)})}
  if(method == 'norm'){scaled_hm_enrich <- apply(hm_enrich_data, 2, function(x){(x-mean(x))/(max(x)-min(x))})}
  scaled_hm_enrich <- apply(scaled_hm_enrich, 2, function(x){
    new_x <- x
    new_x[is.nan(new_x)] <- NA
    new_x[is.infinite(new_x)] <- NA
    return(new_x)
  })
  return(scaled_hm_enrich)
}


#' Get vector of colors for heatmap
#'
#' @param pal name of palette for brewer.pal function (RdYlBu, RdBu, Spectral, PiYG etc.)
#' @param pal_reverse reverse or not the order of colors in palette
#' @param range number of demanded colours
#'
#' @return
#'
get_color_enrichments <- function(pal = "RdYlBu", pal_reverse = TRUE, range = 100){
  color <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(
    RColorBrewer::brewer.pal.info[pal, "maxcolors"], pal))(range)
  if(pal_reverse){color <- rev(color)}
  return(color)
}


#' get tidy data frame for deconvolution heatmap plot
#'
#' @param exprs_data table with enrichment data
#' @param mk marker for deconvolution
#' @param samples vector with sample information for each cell of enrichment data
#' @param groups vector with group information for each cell of enrichment data
#' @param summarize_mode method to calculate enrichment of cells groups("median" or "mean")
#'
#' @return
get_deconvol_enrichments <- function(exprs_data, mk, samples, groups, summarize_mode = 'median'){
  deconv_expr <- data.frame(samples = samples, groups = groups, enrich = exprs_data[,mk])
  col_names <- unique(deconv_expr$groups)
  row_names <- unique(deconv_expr$samples)
  deconv_expr <- as.data.frame(lapply(col_names, function(g){
    sapply(row_names , function(s){
      enrich <- deconv_expr$enrich[(deconv_expr$samples == s) & (deconv_expr$groups == g)]
      if(summarize_mode == 'median'){return(stats::median(enrich))}
      if(summarize_mode == 'mean'){return(mean(enrich))}
    })
  }))
  colnames(deconv_expr) <- col_names
  deconv_expr$samples <- row_names
  deconv_expr <-reshape2::melt(deconv_expr, id.vars = "samples", variable.name = "groups", value.name = "enrich")
  deconv_expr$cell_rate <- apply(deconv_expr, 1, function(x){
    n1 <- sum((samples == x[1]) & (groups == x[2]))
    n2 <- sum((samples == x[1]))
    n1/n2
    })
  return(deconv_expr)
}
