#' softThreshold.Seurat
#'
#' @rdname softThreshold
#' @method softThreshold Seurat
#'
#' @export
#'
softThreshold.Seurat <- function(
  object = NULL,
  sealing = 10
) {
  if(!is.null(object)){
    # determine the soft threshold.
    if(exists(object$percent.mt)){
      mito <- object$percent.mt
      mito <- sort(mito)
      mito.cut <- mito[floor(length(mito)*0.95)]
      if(mito.cut > sealing){
        mito.cut = sealing
      }
      return(mito.cut)
    }else{
      stop("There is no percent.mt information in your Seurat object! Please run PercentageFeatureSet function in Seurat, or create your object use createObject function first!")
    }
  } else {
    stop("Please provide Seurat object!")
  }
}


#' softThreshold.default
#'
#' @rdname softThreshold
#' @method softThreshold default
#'
#' @export
#'
softThreshold.default <- function(
  object = NULL,
  sealing = 10
) {
  if(!is.null(object)){
    # determine the soft threshold.
    mito <- object
    mito <- sort(mito)
    mito.cut <- mito[floor(length(mito)*0.95)]
    if(mito.cut > sealing){
      mito.cut = sealing
    }
    return(mito.cut)
  } else {
    stop("Please provide Seurat object!")
  }
}


#' removeGene
#'
#' remove unwanted genes by gene name pattern (e.g. Ig genes)
#' For some cells, the expression of some specific genes may not related to the real cell popultion but will affect the cell clustering results. User can remove those genes from raw data before perform further analysis
#' An example is Ig genes for B cells.
#'
#' @param object Seurat object
#' @param pattern the pattern (regular expression) of unwanted genes (e.g. '^IG[HKL]' for IgH, IgK, IgL genes). For more information about regular expression, please visit https://en.wikipedia.org/wiki/Regular_expression
#'
#' @export
#'
removeGene <- function(
  object = NULL,
  pattern = NULL
) {
  if(!is.null(object) & !is.null(pattern)){
    # remove some unwanted genes from RNA data
    cat("remove genes from your RNA assay data... \n")
    assay.data <- GetAssayData(object = object, slot = "counts",assay = "RNA")
    remove.index <- grep(pattern = pattern, x = rownames(assay.data), value = FALSE)
    object@assays[["RNA"]]@counts <- object@assays[["RNA"]]@counts[-remove.index,]
    object@assays[["RNA"]]@data <- object@assays[["RNA"]]@data[-remove.index,]
    object@assays[["RNA"]]@meta.features <- object@assays[["RNA"]]@meta.features[-remove.index,]
    return(object)
  }else if (is.null(object) ) {
    stop("Please provide an Seurat object!")
  } else {
    warning("You didn't provide any pattern, do nothing to your object")
    return(object)
  }
}


#' removeProtein
#'
#' remove unwanted surface protein signals by protein names
#'
#' @param object Seurat object
#' @param names the name of proteins you want to remove from current analysis
#'
#' @export
removeProtein <- function(
  object = NULL,
  names = NULL
) {
  if(!is.null(object) & !is.null(names)){
    # remove some unwanted genes from ADT data
    cat("remove proteins from your ADT assay data... \n")
    assay.data <- GetAssayData(object = object, slot = "counts",assay = "ADT")
    remove.index <- which(rownames(assay.data) %in% names)
    object@assays[["ADT"]]@counts <- object@assays[["ADT"]]@counts[-remove.index,]
    object@assays[["ADT"]]@data <- object@assays[["ADT"]]@data[-remove.index,]
    object@assays[["ADT"]]@meta.features <- object@assays[["ADT"]]@meta.features[-remove.index,]
    return(object)
  }else if (is.null(object) ) {
    stop("Please provide an Seurat object!")
  } else {
    warning("You didn't provide any protein names, do nothing to your object")
    return(object)
  }
}


#' dataNormalization
#'
#' data normalization for both RNA and ADT data
#' a key step for pre-processing. Call seurat function
#'
#' @param object Seurat object
#' @param scale.factor Sets the scale factor for cell-level normalization. Default is 10000
#' @param margin If performing CLR normalization, normalize across features (1) or cells (2)
#' @param verbose display progress bar for normalization procedure
#'
#' @export
dataNormalization <- function(
  object = NULL,
  scale.factor = 10000,
  margin = 1,
  verbose = TRUE
) {
  if(!is.null(object)) {
    assays <- names(object@assays)
    for (assay in assays) {
      if(assay == "RNA") {
        object <- NormalizeData(object, assay = assay, normalization.method = "LogNormalize", scale.factor = scale.factor, margin = margin, verbose = verbose)
      } else if (assay == "ADT") {
        object <- NormalizeData(object, assay = assay, normalization.method = "CLR", margin = margin, verbose = verbose)
      } else if (assay == "HTO") {
        object <- NormalizeData(object, assay = assay, normalization.method = "CLR", margin = margin, verbose = verbose)
      } else {
        cat("We don't pre-define normalization method for your assay named ",assay,", please run Seurat function NormalizeData for this assay individually!")
      }
    }
    return(object)
  } else {
    stop("Please provide an Seurat object!")
  }
}


#' dataScaling
#'
#' data scaling for both RNA and ADT (for later Heatmap viosualization)
#'
#' @param object = NULL, A Seurat object
#' @param features = NULL, features on data scaling
#' @param assay = NULL,  assay name, RNA, ADT, HTO
#' @param vars.to.regress = NULL, vars.to.regress, same as Seurat,  Variables to regress out (previously latent.vars in RegressOut). For example, nUMI, or percent.mito.
#' @param model.use = "linear", same as Seurat, Use a linear model or generalized linear model (poisson, negative binomial) for the regression. Options are 'linear' (default), 'poisson', and 'negbinom'
#' @param use.umi = FALSE, same as Seurat, Regress on UMI count data. Default is FALSE for linear modeling, but automatically set to TRUE if model.use is 'negbinom' or 'poisson'
#' @param do.scale = TRUE, same as Seurat, Whether to scale the data.
#' @param do.center = TRUE, same as Seurat, Whether to center the data.
#' @param scale.max = 10,  same as Seurat,  Max value to return for scaled data. The default is 10. Setting this can help reduce the effects of feautres that are only expressed in a very small number of cells.
#' @param block.size = 1000, same as Seurat, Default size for number of feautres to scale at in a single computation. Increasing block.size may speed up calculations but at an additional memory cost.
#' @param min.cells.to.block = 3000, same as Seurat, If object contains fewer than this number of cells, don't block for scaling calculations.
#' @param verbose = TRUE, same as Seurat, Displays a progress bar for scaling procedure
#'
#' @export
dataScaling <- function(
  object = NULL,
  features = NULL,
  assay = NULL,
  vars.to.regress = NULL,
  model.use = "linear",
  use.umi = FALSE,
  do.scale = TRUE,
  do.center = TRUE,
  scale.max = 10,
  block.size = 1000,
  min.cells.to.block = 3000,
  verbose = TRUE
) {
  if(!is.null(object)) {
    assays <- names(object@assays)
    for (assay in assays) {
      cat("Scaling for assay",assay, "...\n")
      object <- ScaleData(
        object = object,
        assay = assay,
        features = features,
        vars.to.regress = vars.to.regress,
        model.use = model.use,
        use.umi = use.umi,
        do.scale = do.scale,
        do.center = do.center,
        scale.max = scale.max,
        block.size = block.size,
        min.cells.to.block = min.cells.to.block,
        verbose = verbose
      )
    }
    return(object)
  } else {
    stop("Please provide an Seurat object!")
  }
}
