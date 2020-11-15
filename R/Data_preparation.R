
#' Create basic metadata for flowSet object
#' @description The function creates two-column dataframe with paths
#' and names of files by the list of paths to fcs files with CyTOF data.
#' this data frame use as base metadata or phenoData in flowSet object
#' from flowCore package.
#'
#' @param fcs_files list with paths
#' @return The dataframe with paths to fcs files within one column and names of these files in another column
#'
get_fcs_metadata <- function(fcs_files){
  md <- data.frame(path = fcs_files)
  md$path <- as.character(md$path)
  md$file_name <- basename(md$path)
  samples <- gsub(".fcs", "", md$file_name)
  md$short_name <- unlist(lapply(1:length(samples), function(x){
    tmp <- unlist(strsplit(samples[x], ""))
    if(length(tmp) <= 12){return(samples[x])}
    return(paste0(paste0(tmp[1:7], collapse = ""), "_smpl_", as.character(x)))
  }))

  return(md)
}

### test
#md <- get_fcs_metadata(c("./../test_data/c13_20190704_hnp_perf_11_0_Alex2.fcs","./../test_data/c14_20190704_hnp_perf_11_0_Alex2.fcs",
#                       "./../test_data/c13_20190704_hnp_perf_11_0_Alez1.fcs","./../test_data/c14_20190704_hnp_perf_11_0_Alez1.fcs"))
#md <- get_fcs_metadata("./../../Toni_data/EC_200117_Freshly labelled PBMCs _1_0/Activation_Activation full panel unstim TILs_033.fcs")
#md <- get_fcs_metadata("./../../cytofBrowser_project/Tape_data/Figure-1_S2_S3_raw/Figure-1_S2_S3_processed.fcs")


#' extract metadata for build-in test dataset
#'
#' @param test_data_dproc flag to identifier build-in data
#'
#' @return
#'
get_test_fcs_metadata <- function(test_data_dproc){
  if(test_data_dproc == 'test_data'){
    test_files <- c("KPC1_stroma.fcs", "KPC2_stroma.fcs", "KPC3_stroma.fcs", "KPC4_stroma.fcs",
                    "KPC5_stroma.fcs", "KPC6_stroma.fcs", "KPC7_stroma.fcs")
    test_files <- system.file("extdata",test_files,package = "cytofBrowserLite")
    md <- get_fcs_metadata(test_files)
  }
  return(md)
}


#' Creat FlowSet object
#' @description The function takes dataframe with column contained  paths to
#' fcs files and creates flowSet object from flowCore package.
#'
#' @param md The input metadata should be dataframe format with
#' a column named path. The function assumed to process the result
#' of the function \dontrun{
#' get_fcs_metadata
#' }
#'
#' @return flowSet object from flowCore packag
#' @importFrom flowCore read.FCS read.flowSet sampleNames "sampleNames<-"
#'
get_fcs_raw <- function(md){
  pathes <- as.vector(md$path)
  fcs_list <- lapply(pathes, flowCore::read.FCS)
  fcs_colnames <- lapply(fcs_list, colnames)
  colnames_test <- unlist(lapply(fcs_colnames, function(x){all(unique(unlist(fcs_colnames)) %in% x)}))
  if(any(!colnames_test)){return(NULL)}
  fcs_raw <- as(fcs_list, "flowSet")
  #fcs_raw <- flowCore::read.flowSet(pathes, transformation = FALSE, truncate_max_range = FALSE)
  #sampleNames(fcs_raw) <- unlist(gsub(".fcs", "", flowCore::sampleNames(fcs_raw)))
  sampleNames(fcs_raw) <- basename(pathes)
  return(fcs_raw)
}

### test
#fcs_raw <- get_fcs_raw(md)


#' Get enrichment data from flowFrame object to matrix format
#'
#' @param fcs_raw flowSet object
#' @importFrom flowCore fsApply exprs
#'
#' @return
#'
get_exprs_data <- function(fcs_raw){
  exprs_data <- flowCore::fsApply(fcs_raw, flowCore::exprs)
  rownames(exprs_data) <- 1:nrow(exprs_data)
  return(exprs_data)
}

### test
#exprs_data <- get_exprs_data(fcs_raw)


#' Create panel data to fowSet object
#' @description extract panel data about markers from flowSet
#' object from flowCore package. The function tries to assign
#' a description to markers and its type (technical or not)
#'
#' @param fcs_raw flowSet object from flowCore package
#'
#' @return data frame with information about the markers
#' @importFrom flowCore pData "pData<-" parameters "parameters<-"
#' @importClassesFrom flowCore flowSet
#'
get_fcs_panel <- function(fcs_raw){
  tech_patterns <-list(computational_tech = c("Time", "Event", "length","Center", "Offset", "Width", "Residual", "tSNE", "<NA>"),
                       marker_tech = c("_BC", "BCKG", "DNA"))
  panel <- as.data.frame(flowCore::pData(flowCore::parameters(fcs_raw[[1]]))[,c("name", "desc")])
  rownames(panel) <- NULL
  panel <- panel[sapply(panel$name, function(x) !any(sapply(tech_patterns$computational_tech, function(y) grepl(y,x)))),]
  panel <- panel[sapply(panel$desc, function(x) !any(sapply(tech_patterns$computational_tech, function(y) grepl(y,x)))),]
  panel <- panel[!is.na(panel$desc),]
  panel$antigen <- sapply(strsplit(panel$desc, "_"), function(x){
    if(length(x) <= 1){return(x)}
    return(paste(x[-c(1)], sep = "_", collapse = "_"))
  })
  panel$marker_class <- "type"
  panel$marker_class[sapply(panel$desc, function(x) any(sapply(tech_patterns$marker_tech, function(y) grepl(y,x))))] <- "state"
  return(panel)
}

### test
#panel <- get_fcs_panel(fcs_raw)


#' Filtering out technical markers
#' @description the function formed list of markers which have class
#' mentioned as "type" in the panel. Other markers with class "state"
#' are considered as technical signals
#'
#' @param panel dataframe with marker panel data
#' @return the list value of which are marker names and names of the
#' list element are names of antibody, relevant to the markers
#'
get_use_marker <- function(panel){
  use_markers <- as.character(panel$name[panel$marker_class == "type"])
  names(use_markers) <- panel$antigen[panel$marker_class == "type"]
  return(use_markers)
}

### test
#use_markers <- get_use_marker(panel)


#' Create entire named vector with obsermations: markers,clusters, qc etc.
#'
#' @param fcs_raw flowSet object
#'
#' @return
#' @importFrom flowCore pData parameters
#'
get_entire_panel <- function(fcs_raw){
  tmp <- as.data.frame(flowCore::pData(flowCore::parameters(fcs_raw[[1]]))[,c("name", "desc")])
  rownames(tmp) <- NULL
  tmp$antigen <- sapply(strsplit(tmp$desc, "_"), function(x){
    if(length(x) <= 1){return(x)}
    return(paste(x[-c(1)], sep = "_", collapse = "_"))
  })
  entire_panel <- as.character(tmp$name)
  names(entire_panel) <- tmp$antigen
  return(entire_panel)
}

### test
#entire_panel <- get_entire_panel(fcs_raw)


#' Detect columns with annotations
#'
#' @param fcs_raw flowDet object
#'
#' @return
get_extr_ann <-function(exprs_data){
  col_names <- colnames(exprs_data)
  expectations <- c("tSNE", "UMAP", "clust", "tsne", "umap", "Clust", "Sampl")
  extr_ann <- unlist(lapply(expectations, function(x){grep(x,col_names)}))
  if(length(extr_ann) == 0){return(NULL)}
  extr_ann <- col_names[extr_ann]
  return(extr_ann)
}

add_extr_ann_to_cell_ann <- function(exprs_data, extr_ann, cell_ann){

}

#' Detecting of implemented transformation by data
#'
#' @param exprs_data data frame with enrichment data
#' @param use_markers list of markers to consider
#' @param mode mode of running the function ("dataset": to consider all data  together, "markers": to consider each marker independently)
#'
#' @return
get_extr_transformations <- function(exprs_data, use_markers, mode = "dataset"){
  trans <- c()
  max_values <- apply(exprs_data[,use_markers], 2, max)
  min_values <- apply(exprs_data[,use_markers], 2, min)
  mean_values <- colMeans(exprs_data[,use_markers])
  if(mode == "dataset"){
    max_values <- max(max_values)
    min_values <- min(min_values)
    mean_values <- mean(mean_values)
  }
  if((max_values > 100) & (min_values >= 0)){return(trans)}
  if((min_values >= 0) & (mean_values > 0) & (mean_values < 2)){trans <- c(trans, "log-like")}
  if((max_values <= 1) & (min_values >= 0)){trans <- c(trans, "normalised")}
  if((min_values <= 0) & (max_values <= 20)){trans <- c(trans, "z-score")}
  return(trans)
}



#' Transformation exprs matrix by asinh fransformation
#'
#' @param exprs_data matrix with expression data
#' @param cofactor digit with used as denominator to transform data
#' @param use_marker list of marker for transformation
#' @return flowSet object with transform data by asinh function and divided by cofactor
#'
exprs_asinh_transformation <- function(exprs_data, use_markers, cofactor = 5){
  exprs_asinh <- exprs_data
  exprs_asinh[,use_markers] <- asinh(exprs_data[,use_markers] / cofactor)
  return(exprs_asinh)
}


#' Detect outliers from exprs data and squeeze from 0 to 1
#'
#' @param exprs_data matrix with expression data
#' @param quantile number from 0 to 1 to set the quantile which is of cutoff threshold
#' @param use_marker list of markers for transformation
#'
#' @return
#' @importFrom stats quantile
#'
exprs_outlier_squeezing <- function(exprs_data, use_markers, quantile = 0.01){
  exprs_outlier_squeeze <- exprs_data
  rng <- stats::quantile(exprs_data[,use_markers], probs = c(quantile, 1-quantile))
  exprs_outlier_squeeze[,use_markers] <- t((t(exprs_data[,use_markers]) - rng[1]) / (rng[2] - rng[1]))
  exprs_outlier_squeeze[,use_markers] <- apply(exprs_outlier_squeeze[,use_markers], 2, function(x){
    new_x <- x
    new_x[new_x < 0] <- 0
    new_x[new_x > 1] <- 1
    return(new_x)
  })
  return(exprs_outlier_squeeze)
}

### test
#exprs_data <- exprs_outlier_squeezing(exprs_data, use_markers)

#' check does the transformation which will be implemented already implemented
#'
#' @param trans_to_do list of transformation to implementation
#' @param trans list of transformations which already implemented
#'
#' @return

check_transformation_list <- function(trans_to_do, trans){
  if(any(c("log","asinh","log-like") %in% trans)){trans <- c(trans, "log","asinh","log-like")}
  trans_to_do <- trans_to_do[!(trans_to_do %in% trans)]
  return(trans_to_do)
}


#' Implement the transformation from trans_to_do if they are not in trans list
#'
#' @param trans_to_do list of transformations for implementation
#' @param exprs_data data frame with enrichment data
#' @param use_markers list of considered markers
#' @param cofactor co-factior for asinh transformation
#' @param quantile quantile for outlier detection
#'
#' @return

get_transformations <- function(trans_to_do, exprs_data, use_markers,
                                cofactor = 5, quantile = 0.01){
  if("log" %in% trans_to_do){
    exprs_data[,use_markers] <- log(exprs_data[,use_markers])}
  if("asinh" %in% trans_to_do){
    exprs_data <- exprs_asinh_transformation(exprs_data, use_markers, cofactor = cofactor)}
  if("outlier_squeezing" %in% trans_to_do){
    exprs_data <- exprs_outlier_squeezing(exprs_data, use_markers, quantile = quantile)}
  if("z-score" %in% trans_to_do){
    exprs_data[,use_markers] <- apply(exprs_data[,use_markers], 2, scale)}
  return(exprs_data)
}

#' Extract cell number
#'
#' @param fcs_raw flowSet object
#'
#' @return
#' @importFrom flowCore sampleNames "sampleNames<-" fsApply
#'
get_cell_number <- function(fcs_raw){
  cell_number <- data.frame(smpl = flowCore::sampleNames(fcs_raw), cell_nmbr = flowCore::fsApply(fcs_raw, nrow))
  colnames(cell_number) <- c('smpl', 'cell_nmbr')
  return(cell_number)
}

### test
#cell_number <- get_cell_number(fcs_raw)

###################
### Subsampling ###
###################

#' Subsempling data from exprs data format with sample and fuze considering
#'
#' @param sampling_size fraction of cell amount to plotting
#' @param cell_ann cell annotation data
#' @param fuse logic value to use fuse for data size
#' @param size_fuse number of cells to as max number as fuze for data size
#'
#' @return
#'
get_subset_coord <- function(cell_ann, sampling_size = 0.5, fuse = TRUE, size_fuse = 3000){
  if(!("samples" %in% colnames(cell_ann))){print("samples info are not in cell_ann")}
  ## How many cells to sampling per-sample
  smpl_cell_amount <- as.integer((table(cell_ann$samples) + 0.1) * sampling_size)
  if(fuse & (sum(smpl_cell_amount) > size_fuse)){
    smpl_cell_amount <- as.integer((smpl_cell_amount/sum(smpl_cell_amount))*size_fuse)}
  names(smpl_cell_amount) <- names(table(cell_ann$samples))
  ## Get subsample indices
  set.seed(1234)
  subset_coord <- unlist(lapply(names(smpl_cell_amount), function(i){
    sample(rownames(cell_ann[cell_ann$samples == i,]), smpl_cell_amount[i], replace = FALSE)}))
  subset_coord <- subset_coord[order(as.integer(subset_coord))]
  return(subset_coord)
}

### test
#cell_ann <- data.frame(all_cells = rep("cell", sum(cell_number$cell_nmbr)),
#                       samples = rep(cell_number$smpl, cell_number$cell_nmbr),
#                       row.names = 1:sum(cell_number$cell_nmbr))

#subset_coord <- get_subset_coord(cell_ann, sampling_size = 0.5, fuse = TRUE, size_fuse = 3000)


############################
###  Dimensional reduce  ###
############################

#' Preparing data to dimensional reduce plot
#'
#' @param exprs_data matrix with enrichment data
#' @param cell_ann cell annotation data
#' @param subset_coord coordinates of subset
#' @param method method to dimensional reducing (tSNE or UMAP)
#' @param perplexity perplexity for tSNE
#' @param theta theta for tSNE
#' @param max_iter max iterations for tSNE
#' @param use_markers list of markers for dimensional reducing
#' @param force flag to make dimensional reducing
#' @param pca_param pca parameters for Rtsne function
#' @param check_duplicates check_duplicates parameters for Rtsne function
#' @param seed seed for randomisation of Rtsne function
#'
#' @return
#' @importFrom flowCore fsApply sampleNames "sampleNames<-" exprs "exprs<-"
#' @importFrom Rtsne Rtsne
#' @importFrom umap umap
#'
get_dim_reduce <- function(exprs_data, cell_ann, subset_coord, use_markers, force = FALSE,
                           method = "tSNE", perplexity = 30, theta = 0.5, max_iter = 1000,
                           pca_param = FALSE, check_duplicates = FALSE, seed = 1234){
  col_position <- grepl(method, colnames(cell_ann))
  if((any(col_position) & !force) & all(!(is.na(unlist(cell_ann[subset_coord,which(col_position)]))))){
    return(cell_ann)
  }
  new_cell_ann <- cell_ann

  if(method == "tSNE"){
    new_cell_ann$tSNE1 <- NA
    new_cell_ann$tSNE2 <- NA
    ##### Run t-SNE
    set.seed(seed)
    tsne_result <- Rtsne::Rtsne(exprs_data[subset_coord, use_markers], check_duplicates = check_duplicates, pca = pca_param,
                                perplexity = perplexity, theta = theta, max_iter = max_iter)
    new_cell_ann[subset_coord, "tSNE1"] <- tsne_result$Y[, 1]
    new_cell_ann[subset_coord, "tSNE2"] <- tsne_result$Y[, 2]
  }

  if(method == "UMAP"){
    new_cell_ann$UMAP1 <- NA
    new_cell_ann$UMAP2 <- NA
    ##### Run UMAP
    umap_out <- umap::umap(exprs_data[subset_coord, use_markers])
    new_cell_ann[subset_coord, "UMAP1"] <- umap_out$layout[, 1]
    new_cell_ann[subset_coord, "UMAP2"] <- umap_out$layout[, 2]
  }
  return(new_cell_ann)
}


#' Get new names with additional prefix for samples to save
#'
#' @param old_name vector of assumed names to save
#' @param folder_path directory to save file
#' @param prefix additional prefix to add to names
#'
#' @return
#'
get_new_fcs_names <- function(old_name, folder_path, prefix = ""){
  if(prefix == "" | is.null(prefix)){prefix <- "CyBr"}
  new_names <- paste(prefix, old_name, sep = "_")
  ovrlp <- new_names %in% list.files(folder_path)
  quantity <- 1
  while (any(ovrlp)){
    quantity <- quantity +1
    new_names <- paste(paste0(prefix, quantity), old_name, sep = "_")
    ovrlp <- new_names %in% list.files(folder_path)
  }
  return(new_names)
}

### test
#get_new_fcs_names(old_name, folder_path)


#' Adding of annotation info to flowSet object
#'
#' @param fcs_raw flowSet object
#' @param cell_ann cell annotation table
#' @param col_name vetor of column names from cell annotation data which we would like to add to fcs
#'
#' @return
#' @importFrom flowCore sampleNames fr_append_cols
#' @importClassesFrom flowCore flowSet
#'
get_ann_adding_fcs_files <- function(fcs_raw, cell_ann, col_name = NULL){
  if(is.null(col_name) | any(!(col_name %in% colnames(cell_ann)))){
    print(paste0("Annotations ", col_name[!(col_name %in% colnames(cell_ann))], " didn't add to files."))
    return(fcs_raw)
  }
  format_cell_ann <- lapply(col_name, function(x){
    as.integer(as.factor(cell_ann[,x]))
  })
  format_cell_ann <- as.data.frame(format_cell_ann)
  colnames(format_cell_ann) <- col_name
  new_fcs <- lapply(flowCore::sampleNames(fcs_raw), function(s) {
    addition_col <- as.matrix(format_cell_ann[cell_ann$samples == s, col_name])
    colnames(addition_col) <- col_name
    flowCore::fr_append_cols(fcs_raw[[s]], addition_col)
  })
  names(new_fcs) <- flowCore::sampleNames(fcs_raw)
  new_fcs <- as(new_fcs, 'flowSet')
  return(new_fcs)
}
