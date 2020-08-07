
##### Create metadata
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
#md <- get_fcs_metadata(c("./test_data/c13_20190704_hnp_perf_11_0_Alex2.fcs","./test_data/c14_20190704_hnp_perf_11_0_Alex2.fcs",
#                        "./test_data/c13_20190704_hnp_perf_11_0_Alez1.fcs","./test_data/c14_20190704_hnp_perf_11_0_Alez1.fcs"))
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
    test_files <- system.file("extdata",test_files,package = "cytofBrowser")
    md <- get_fcs_metadata(test_files)
  }
  return(md)
}


##### creat FlowSet object
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
#' @importFrom flowCore read.flowSet sampleNames "sampleNames<-"
#'

get_fcs_raw <- function(md){
  pathes <- as.vector(md$path)
  fcs_raw <- flowCore::read.flowSet(pathes, transformation = FALSE, truncate_max_range = FALSE)
  sampleNames(fcs_raw) <- unlist(gsub(".fcs", "", flowCore::sampleNames(fcs_raw)))
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
  #exprs_data <- cbind(exprs_data, as.matrix(data.frame(Sample = rep(1:length(fcs_raw), fsApply(fcs_raw, nrow)))))
  rownames(exprs_data) <- 1:nrow(exprs_data)
  return(exprs_data)
}

##### Create panel data
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
  #panel$antigen <- sapply(strsplit(panel$desc, "_"), function(x) x[length(x)])
  #panel$antigen <- gsub(" \\(v)", "", panel$antigen)
  panel$marker_class <- "type"
  panel$marker_class[sapply(panel$desc, function(x) any(sapply(tech_patterns$marker_tech, function(y) grepl(y,x))))] <- "state"
  #panel <- as.data.frame(apply(panel, c(1,2), function(x) gsub(" ", "_", x)))
  return(panel)
}

### test
#panel <- get_fcs_panel(fcs_raw)


##### Create use_marker
#' filtering out technical markers
#' @description the function formed list of markers which have class
#' mentioned as "type" in the panel. Other markers with class "state"
#' are considered as technical signals
#'
#' @param panel dataframe with marker panel data
#' @return the list value of which are marker names and names of the
#' list element are names of antibody, relevant to the markers
#'
#'

get_use_marker <- function(panel){
  use_markers <- as.character(panel$name[panel$marker_class == "type"])
  names(use_markers) <- panel$antigen[panel$marker_class == "type"]
  #names(use_markers) <- gsub("_(v)", "", names(use_markers), fixed = T)          ### It should be fixed (problem with "_(v)")
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


##### Upload data from fcs files
upload_fcs_data <- function(fcs_files){
  md <- get_fcs_metadata(fcs_files)
  fcs_raw <- get_fcs_raw(md)
  panel <- get_fcs_panel(fcs_raw)
  use_markers <- get_use_marker(panel)
  return(list(fcs_raw,md, panel, use_markers))
}


#### Transformation from "count" to "asinh" data
#' Transformation from "count" to "asinh" data
#'
#' @param fcs_raw flowSet object
#' @param cofactor digit with used as denominator to transform data
#' @param use_marker list of marker for transformation
#' @return flowSet object with transform data by asinh function and divided
#' by cofactor
#' @importClassesFrom flowCore flowSet
#' @importFrom flowCore fsApply exprs
#'

asinh_transformation <- function(fcs_raw, cofactor = 5, use_markers = NULL){
  markers <- use_markers
  if(is.null(use_markers)){
    computational_tech <- c("Time", "Event", "length","Center", "Offset", "Width", "Residual", "tSNE", "clust")
    markers <- colnames(fcs_raw[[1]])[sapply(colnames(fcs_raw[[1]]), function(x) !any(sapply(computational_tech, function(y) grepl(y,x))))]
  }
  fcs_asinh <- flowCore::fsApply(fcs_raw, function(x, cf = cofactor, mk = markers){
    exprs(x)[,mk] <- asinh(exprs(x)[,mk] / cf)
    return(x)})
  return(fcs_asinh)
}

### test
#fcs_raw <- asinh_transformation(fcs_raw, 5)


##### Transformation to a from 0 to 1 variable and removing outliers
#' Transformation to a from 0 to 1 variable and removing outliers
#'
#' @param fcs_raw flowSet object
#' @param quantile number from 0 to 1 to set the quantile which is of cutoff threshold
#' @param use_marker list of markers for transformation
#'
#' @return
#' @importClassesFrom flowCore flowSet
#' @importFrom flowCore fsApply exprs "exprs<-"
#' @importFrom stats quantile
#'

outlier_by_quantile_transformation <- function(fcs_raw, quantile = 0.01, use_markers = NULL){
  markers <- use_markers
  if(is.null(use_markers)){
    computational_tech <- c("Time", "Event", "length","Center", "Offset", "Width", "Residual", "tSNE")
    markers <- colnames(fcs_raw[[1]])[sapply(colnames(fcs_raw[[1]]), function(x) !any(sapply(computational_tech, function(y) grepl(y,x))))]
  }
  fcs_outlier_by_quantile <- flowCore::fsApply(fcs_raw, function(x, ql = quantile, mk = markers){
    rng <- stats::quantile(exprs(x)[,mk], probs = c(ql, 1-ql))
    expr_data <- t((t(exprs(x)[,mk]) - rng[1]) / (rng[2] - rng[1]))
    expr_data[expr_data < 0] <- 0
    expr_data[expr_data > 1] <- 1
    exprs(x)[,mk] <- expr_data
    return(x)
  })
  return(fcs_outlier_by_quantile)
}

### test
#fcs_raw <- outlier_by_quantile_transformation(fcs_raw, 0.01)

##### Extract cell number
#' Extract cell number
#'
#' @param fcs_raw flowSet object
#'
#' @return
#' @importFrom flowCore sampleNames "sampleNames<-" fsApply


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
  print("-------0")
  if(!("samples" %in% colnames(cell_ann))){print("samples info are not in cell_ann")}
  print("-------1")
  ## How many cells to sampling per-sample
  smpl_cell_amount <- as.integer((table(cell_ann$samples) + 0.1) * sampling_size)
  print("-------2")
  print(smpl_cell_amount)
  if(fuse & (sum(smpl_cell_amount) > size_fuse)){
    smpl_cell_amount <- as.integer((smpl_cell_amount/sum(smpl_cell_amount))*size_fuse)}
  print("-------3")
  print(str(smpl_cell_amount))
  print(table(cell_ann$samples))
  names(smpl_cell_amount) <- names(table(cell_ann$samples))

  print("-------4")
  ## Get subsample indices
  set.seed(1234)
  subset_coord <- unlist(lapply(names(smpl_cell_amount), function(i){
    sample(rownames(cell_ann[cell_ann$samples == i,]), smpl_cell_amount[i], replace = FALSE)}))
  print("-------5")
  subset_coord <- subset_coord[order(as.integer(subset_coord))]
  print(head(subset_coord))
  print(length(subset_coord))
  print("-------6")
  return(subset_coord)
}

### test
#cell_ann <- data.frame(all_cells = rep("cell", sum(cell_number$cell_nmbr)),
#                       samples = rep(cell_number$smpl, cell_number$cell_nmbr),
#                       row.names = 1:sum(cell_number$cell_nmbr))

#subset_coord <- get_subset_coord(cell_ann, sampling_size = 0.5, fuse = TRUE, size_fuse = 3000)

############################
###  Dimentional reduce  ###
############################

##### Preparing data to dimentional reduce plot
#' Preparing data to dimentional reduce plot
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
  print("--------0")
  col_position <- grepl(method, colnames(cell_ann))
  print("--------1")
  print(col_position)
  print(head(cell_ann[subset_coord,which(col_position)]))
  print(head(unlist(cell_ann[subset_coord,which(col_position)])))
  print(table(is.na(unlist(cell_ann[subset_coord,which(col_position)]))))
  print(table(!(is.na(unlist(cell_ann[subset_coord,which(col_position)])))))
  print(all(!(is.na(unlist(cell_ann[subset_coord,which(col_position)])))))
  if((any(col_position) & !force) & all(!(is.na(unlist(cell_ann[subset_coord,which(col_position)]))))){
    return(cell_ann)
  }
  print("--------2")
  print(head(cell_ann[,col_position]))

  new_cell_ann <- cell_ann

  if(method == "tSNE"){
    new_cell_ann$tSNE1 <- NA
    new_cell_ann$tSNE2 <- NA
    print("--------3")
    print(head(new_cell_ann))
    print(perplexity)
    ##### Run t-SNE
    set.seed(seed)
    tsne_result <- Rtsne::Rtsne(exprs_data[subset_coord, use_markers], check_duplicates = check_duplicates, pca = pca_param,
                                perplexity = perplexity, theta = theta, max_iter = max_iter)
    print("--------4")
    print(head(new_cell_ann))
    print(dim(exprs_data[subset_coord, use_markers]))
    print(head(exprs_data[subset_coord, use_markers]))
    print(str(tsne_result))
    print(dim(tsne_result$Y))
    print(head(tsne_result$Y))
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
  print("--------5")
  return(new_cell_ann)
}


#################
### tSNE plot ###
#################

##### Preparing data to tSNE
#' Preparing data to tSNE
#'
#' @param fcs_raw flowSet object
#' @param use_markers list of markers for tSNE
#' @param sampling_size fraction of cell amount for analysis
#'
#' @return
#' @importFrom flowCore fsApply sampleNames "sampleNames<-" exprs "exprs<-"
#' @importFrom Rtsne Rtsne
#' @importFrom umap umap
#'

scatter_plot_data_prep <- function(fcs_raw, use_markers, sampling_size = 0.5, method = "tSNE",
                         perplexity = 30, theta = 0.5, max_iter = 1000, size_fuse = 5000){
  #sampling_size <- as.integer(sampling_size/length(fcs_raw))
  expr <- flowCore::fsApply(fcs_raw[,use_markers], flowCore::exprs)
  sample_ids <- rep(flowCore::sampleNames(fcs_raw), flowCore::fsApply(fcs_raw, nrow))

  ## Find and skip duplicates
  dups <- which(!duplicated(expr[, use_markers]))

  ## Data subsampling: create indices by sample
  inds <- split(1:length(sample_ids), sample_ids)

  ## How many cells to downsample per-sample
  #tsne_ncells <- pmin(table(sample_ids), sampling_size)             ################ Number of cells to ploting
  tsne_ncells <- as.integer((table(sample_ids) + 0.1) * sampling_size)
  if((!is.null(size_fuse) & !is.na(size_fuse)) & (sum(tsne_ncells) > size_fuse)){
    tsne_ncells <- as.integer((tsne_ncells/sum(tsne_ncells))*size_fuse)}
  names(tsne_ncells) <- names(table(sample_ids))

  ## Get subsampled indices
  set.seed(1234)
  tsne_inds <- lapply(names(inds), function(i){
    s <- sample(inds[[i]], tsne_ncells[i], replace = FALSE)
    intersect(s, dups)
  })

  tsne_inds <- unlist(tsne_inds)
  tsne_expr <- expr[tsne_inds, use_markers]

  if(method == "tSNE"){
    ##### Run t-SNE
    set.seed(1234)
    tsne_result <- Rtsne::Rtsne(tsne_expr, check_duplicates = FALSE, pca = FALSE,
                         perplexity = perplexity, theta = theta, max_iter = max_iter)
    #tsne_out <- data.frame(tSNE1 = tsne_result$Y[, 1], tSNE2 = tsne_result$Y[, 2])
    tsne_out <- data.frame(tSNE1 = tsne_result$Y[, 1], tSNE2 = tsne_result$Y[, 2], expr[tsne_inds, use_markers])
    colnames(tsne_out) <- c("tSNE1", "tSNE2", use_markers)
  }

  if(method == "UMAP"){
    ##### Run UMAP
    umap_out <- umap::umap(tsne_expr)
    tsne_out <- data.frame(tSNE1 = umap_out$layout[, 1], tSNE2 = umap_out$layout[, 2], expr[tsne_inds, use_markers])
    colnames(tsne_out) <- c("tSNE1", "tSNE2", use_markers)
  }

  colnames(tsne_out)[match(use_markers, colnames(tsne_out))] <- names(use_markers)
  return(tsne_out)
}

#tSNE <- scatter_plot_data_prep(fcs_raw, use_markers, sampling_size = 0.1, method = "tSNE")
#ggplot(tSNE,  aes(x = tSNE1, y = tSNE2, color = tSNE[,names(use_markers)[10]])) +
#  geom_point(size = 0.2) +
#  labs(color = names(use_markers)[10])
