#' jointDistance
#'
#' If input seurat object, will calculate cell-cell pairwise distances for RNA and ADT separately, then calculate the joint cell-cell pairwise distances. It alao can calculate joint distances for any two dist objects
#'
#' @param object (For Seurat) Seurat object
#' @param dims (For Seurat) number of PCs used for RNA data. Default is 20
#' @param alpha (For Seurat) use alpha to balence contributions from RNA and ADT in the joint distances. We suggest use 0.5 for initial analysis, and then adjust alpha for better results. User can set alpha = null, then we will automatically determine the value
#' @param dist1 (For default) dist object of first assay
#' @param dist2 (For default) dist object of second assay
#' @param alpha (For default) use alpha to balence contributions from two assays in joint distances. We suggest use 0.5 for initial analysis, and then adjust alpha for better results. User can set alpha = null, then we will automatically determine the value
#'
#' @export
#'
jointDistance <- function(object, ...) {
  UseMethod(generic = 'jointDistance', object = object)
}


#' umapFromDistane
#'
#' Run UMAP to project single cells into 2D space using cell-cell pairwise distances
#'
#' @param object (For Seurat) Seurat object
#' @param assay (For Seurat) run UMAP for which assay, choose from RNA, ADT, Joint or All
#' @param dist (For default) dist object or dist matrix
#' @param seed see number. default is 42
#' @param method could be "naive" or "umap-learn"
#' @param n.neighbors integer; number of nearest neighbors
#' @param n.components  integer; dimension of target (output) space
#' @param metric character or function; determines how distances between data points are computed. When using a string, available metrics are: euclidean, manhattan. Other available generalized metrics are: cosine, pearson, pearson2. Note the triangle inequality may not be satisfied by some generalized metrics, hence knn search may not be optimal. When using metric.function as a function, the signature must be function(matrix, origin, target) and should compute a distance between the origin column and the target columns
#' @param verbose logical or integer; determines whether to show progress messages
#' @param n.epochs  integer; number of iterations performed during layout optimization
#' @param min.dist  numeric; determines how close points appear in the final layout
#' @param spread numeric; used during automatic estimation of a/b parameters.
#' @param set.op.mix.ratio numeric in range [0,1]; determines who the knn-graph is used to create a fuzzy simplicial graph
#' @param local.connectivity  numeric; used during construction of fuzzy simplicial set
#' @param negative.sample.rate  integer; determines how many non-neighbor points are used per point and per iteration during layout optimization
#'
#' @export
#'
umapFromDistane <- function(object, ...) {
  UseMethod(generic = 'umapFromDistane', object = object)
}


#' tsneFromDistane
#'
#' Run t-SNE to project single cells into 2D map using cell-cell pairwise distances
#'
#' @param object (For Seurat) Seurat object
#' @param assay (For Seurat) run t-SNE for which assay, choose from RNA, ADT, Joint or All
#' @param dist (For default) dist object or dist matrix
#' @param perplexity numeric; Perplexity parameter (should not be bigger than 3 * perplexity < nrow(X) - 1, see details for interpretation)
#' @param dim integer; Output dimensionality (default: 2)
#' @param seed integer; seed for reproducible results.
#' @param theta numeric; Speed/accuracy trade-off (increase for less accuracy), set to 0.0 for exact TSNE (default: 0.5)
#'
#' @export
#'
tsneFromDistane <- function(object, ...) {
  UseMethod(generic = 'tsneFromDistane', object = object)
}

#' buildMST
#'
#' If input Seurat object, build Minimum spanning tree (MST) based on user's choice (RNA, ADT, Joint or All). Also can build MST for any input dist object
#'
#' @param object Seurat object, or dist object, or N*N distance martix
#' @param assay If object = Seurat object, user choose the assay data they want to build a MST based on, can be RNA, ADT, Joint or All (All means calculate all three)
#'
#' @export
#'
buildMST <- function(object, ...) {
  UseMethod(generic = 'buildMST', object = object)
}


#' buildPG
#'
#' build PG based on user's choice (RNA, ADT, Joint or All) MST
#'
#' @param object (For Seurat) Seurat object
#' @param assay (For Seurat) user choose the assay data they want to build a MST based on, can be RNA, ADT, Joint or All (All means calculate all three)
#' @param reduction.prefix (For Seurat) the prefix of reduction. e.g. "tsne_" or "umap_"
#' @param object (For default) cell embedding
#' @param mst (For default) MST for building PG
#' @param pg.nodes parameter for computing Ecipal Graph. number of nodes you want to build PG. if NULL, we will calculated a number based on your data size
#' @param pg.min.nodes parameter for computing Elastic Principal Graph. min value of pg.nodes
#' @param pg.Lambda parameter for computing Elastic Principal Graph, real, the lambda parameter used the compute the elastic energy
#' @param pg.Mu parameter for computing Elastic Principal Graph, real, the lambda parameter used the compute the elastic energy
#' @param trimming.radius parameter for computing Elastic Principal Graph, real, maximal distance of point from a node to affect its embedment
#' @param final.energy parameter for computing Elastic Principal Graph, string indicating the final elastic emergy associated with the configuration. Currently it can be "Base" or "Penalized"
#' @param initial.MST use MST as initial for PG or not. default is FALSE because use MST as initial is time consuming for large dataset
#'
#' @export
#'
buildPG <- function(object, ...) {
  UseMethod(generic = 'buildPG', object = object)
}


#' buildKNN
#'
#' build knn network based on user's choice (RNA, ADT, Joint or All)
#'
#' @param object Seurat object, or dist object
#' @param assay (For Seurat) user choose the assay data they want to build a MST based on, can be RNA, ADT, Joint or All (All means calculate all three)
#' @param k k value for KNN graph. Default = 10
#'
#' @export buildKNN
#'
buildKNN <- function(object, ...) {
  UseMethod(generic = 'buildKNN', object = object)
}


#' softThreshold
#'
#' determine the soft threshold of percentage of mitochondrial gene. The criteria is: for each sample, we only keep the 95 percent cells with lowest percent of mito genes. If the cell quality is too bad, we use 10 percent as a sealing threshold
#'
#' @title softThreshold
#' @param object The object
#' @param ... Arguments passed to other methods (ignored for now)
#'
#' @param object (For Seurat) Seurat object
#' @param sealing (For Seurat) the sealing point of mitochondrial gene threshold, for 10X data we suggest use 10 percent
#' @param data (For default) Array of mito percent
#' @param sealing (For default) the sealing point of mitochondrial gene threshold, for 10X data we suggest use 10 percent
#'
#' @export
#'
softThreshold <- function(object, ...) {
  UseMethod(generic = 'softThreshold', object = object)
}


