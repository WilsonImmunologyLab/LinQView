#' @importFrom Rcpp evalCpp
#'
NULL

#' clusteringFromDistance
#'
#' Run clustering method (implemented by Seurat package) to identify cell populations using cell-cell pairwise distances
#'
#' @param object Seurat object
#' @param assay run UMAP for which assay, choose from RNA, ADT, Joint or All
#' @param resolution resolution for 1) Joint, 2) RNA and 3) ADT clustering. if assay = All, user should provide all three solutions. otherwise user only need to provide one resolution for selected assay.
#' @param graph.k.param parameter for Seurat function FindNeighbors. Defines k for the k-nearest neighbor algorithm
#' @param graph.compute.SNN parameter for Seurat function FindNeighbors. also compute the shared nearest neighbor graph
#' @param graph.prune.SNN parameter for Seurat function FindNeighbors. Sets the cutoff for acceptable Jaccard index when computing the neighborhood overlap for the SNN construction. Any edges with values less than or equal to this will be set to 0 and removed from the SNN graph. Essentially sets the strigency of pruning (0 — no pruning, 1 — prune everything).
#' @param graph.nn.eps parameter for Seurat function FindNeighbors. Error bound when performing nearest neighbor seach using RANN; default of 0.0 implies exact nearest neighbor search
#' @param graph.force.recalc parameter for Seurat function FindNeighbors. Force recalculation of SNN.
#' @param cluster.modularity.fxn parameter for Seurat function FindClusters.Modularity function (1 = standard; 2 = alternative).
#' @param cluster.initial.membership parameter for Seurat function FindClusters.Parameters to pass to the Python leidenalg function.
#' @param cluster.node.sizes parameter for Seurat function FindClusters.Parameters to pass to the Python leidenalg function.
#' @param cluster.algorithm parameter for Seurat function FindClusters.Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python.
#' @param cluster.n.start parameter for Seurat function FindClusters.Number of random starts.
#' @param cluster.n.iter parameter for Seurat function FindClusters.Maximal number of iterations per random start.
#' @param cluster.random.seed parameter for Seurat function FindClusters.Seed of the random number generator.
#' @param cluster.group.singletons parameter for Seurat function FindClusters.Group singletons into nearest cluster. If FALSE, assign all singletons to a "singleton" group
#' @param cluster.temp.file.location parameter for Seurat function FindClusters.Directory where intermediate files will be written. Specify the ABSOLUTE path.
#' @param cluster.edge.file.name parameter for Seurat function FindClusters.Edge file to use as input for modularity optimizer jar.
#'
#' @export
clusteringFromDistance <- function(
  object,
  assay = "All",
  resolution = c(0.6, 0.6, 0.6),
  graph.k.param = 20,
  graph.compute.SNN = TRUE,
  graph.prune.SNN = 1/15,
  graph.nn.eps = 0,
  graph.force.recalc = FALSE,
  cluster.modularity.fxn = 1,
  cluster.initial.membership = NULL,
  cluster.node.sizes = NULL,
  cluster.algorithm = 1,
  cluster.n.start = 10,
  cluster.n.iter = 10,
  cluster.random.seed = 0,
  cluster.group.singletons = TRUE,
  cluster.temp.file.location = NULL,
  cluster.edge.file.name = NULL
) {
  if(assay == "All") {
    if(length(resolution) < 3) {
      cat("You set assay = All, but didn't provide 3 resolutions, will use default value 0.6 for all three\n")
      resolution = c(0.6,0.6,0.6)
    } else {
      cat("You set resolution for RNA = ",resolution[2],", For ADT = ",resolution[3],", For joint = ",resolution[1], "...\n")
    }
  }
  if(!is.null(object)){
    if(assay == "All") {
      # for RNA
      cat("Start run clustering for RNA using cell-cell distances...\n")

      # identify cell clusters by calling Seurat function
      object[["rna_snn"]] <- FindNeighbors(object = object@misc[['RNA']][['dist']], k.param = graph.k.param, compute.SNN = graph.compute.SNN, prune.SNN = graph.prune.SNN, nn.eps = graph.nn.eps, verbose = FALSE, force.recalc = graph.force.recalc)$snn
      object <- FindClusters(object = object, resolution = resolution[2], graph.name = "rna_snn", modularity.fxn = cluster.modularity.fxn, initial.membership = cluster.initial.membership, node.sizes = cluster.node.sizes, algorithm = cluster.algorithm, n.start = cluster.n.start, n.iter = cluster.n.iter, random.seed = cluster.random.seed, group.singletons = cluster.group.singletons, temp.file.location = cluster.temp.file.location, edge.file.name = cluster.edge.file.name, verbose = FALSE)
      object[["rnaClusterID"]] <- Idents(object = object)
      object@misc[['RNA']][['cluster']] <- "rnaClusterID"

      # for ADT
      cat("Start run clustering for ADT using cell-cell distances...\n")

      # identify cell clusters by calling Seurat function
      object[["adt_snn"]] <- FindNeighbors(object = object@misc[['ADT']][['dist']], k.param = graph.k.param, compute.SNN = graph.compute.SNN, prune.SNN = graph.prune.SNN, nn.eps = graph.nn.eps, verbose = FALSE, force.recalc = graph.force.recalc)$snn
      object <- FindClusters(object = object, resolution = resolution[3], graph.name = "adt_snn", modularity.fxn = cluster.modularity.fxn, initial.membership = cluster.initial.membership, node.sizes = cluster.node.sizes, algorithm = cluster.algorithm, n.start = cluster.n.start, n.iter = cluster.n.iter, random.seed = cluster.random.seed, group.singletons = cluster.group.singletons, temp.file.location = cluster.temp.file.location, edge.file.name = cluster.edge.file.name, verbose = FALSE)
      object[["adtClusterID"]] <- Idents(object = object)
      object@misc[['ADT']][['cluster']] <- "adtClusterID"

      # for Joint
      cat("Start run clustering for Joint data using cell-cell distances...\n")

      # identify cell clusters by calling Seurat function
      object[["joint_snn"]] <- FindNeighbors(object = object@misc[['Joint']][['dist']], k.param = graph.k.param, compute.SNN = graph.compute.SNN, prune.SNN = graph.prune.SNN, nn.eps = graph.nn.eps, verbose = FALSE, force.recalc = graph.force.recalc)$snn
      object <- FindClusters(object = object, resolution = resolution[1], graph.name = "joint_snn", modularity.fxn = cluster.modularity.fxn, initial.membership = cluster.initial.membership, node.sizes = cluster.node.sizes, algorithm = cluster.algorithm, n.start = cluster.n.start, n.iter = cluster.n.iter, random.seed = cluster.random.seed, group.singletons = cluster.group.singletons, temp.file.location = cluster.temp.file.location, edge.file.name = cluster.edge.file.name, verbose = FALSE)
      object[["jointClusterID"]] <- Idents(object = object)
      object@misc[['Joint']][['cluster']] <- "jointClusterID"

    } else if (assay == "Joint") {
      # for Joint
      cat("Start run clustering for Joint data using cell-cell distances...\n")

      # identify cell clusters by calling Seurat function
      object[["joint_snn"]] <- FindNeighbors(object = object@misc[['Joint']][['dist']], k.param = graph.k.param, compute.SNN = graph.compute.SNN, prune.SNN = graph.prune.SNN, nn.eps = graph.nn.eps, verbose = FALSE, force.recalc = graph.force.recalc)$snn
      object <- FindClusters(object = object, resolution = resolution[1], graph.name = "joint_snn", modularity.fxn = cluster.modularity.fxn, initial.membership = cluster.initial.membership, node.sizes = cluster.node.sizes, algorithm = cluster.algorithm, n.start = cluster.n.start, n.iter = cluster.n.iter, random.seed = cluster.random.seed, group.singletons = cluster.group.singletons, temp.file.location = cluster.temp.file.location, edge.file.name = cluster.edge.file.name, verbose = FALSE)
      object[["jointClusterID"]] <- Idents(object = object)
      object@misc[['Joint']][['cluster']] <- "jointClusterID"

    } else if (assay == "ADT") {
      # for ADT
      cat("Start run clustering for ADT using cell-cell distances...\n")

      # identify cell clusters by calling Seurat function
      object[["adt_snn"]] <- FindNeighbors(object = object@misc[['ADT']][['dist']], k.param = graph.k.param, compute.SNN = graph.compute.SNN, prune.SNN = graph.prune.SNN, nn.eps = graph.nn.eps, verbose = FALSE, force.recalc = graph.force.recalc)$snn
      object <- FindClusters(object = object, resolution = resolution[1], graph.name = "adt_snn", modularity.fxn = cluster.modularity.fxn, initial.membership = cluster.initial.membership, node.sizes = cluster.node.sizes, algorithm = cluster.algorithm, n.start = cluster.n.start, n.iter = cluster.n.iter, random.seed = cluster.random.seed, group.singletons = cluster.group.singletons, temp.file.location = cluster.temp.file.location, edge.file.name = cluster.edge.file.name, verbose = FALSE)
      object[["adtClusterID"]] <- Idents(object = object)
      object@misc[['ADT']][['cluster']] <- "adtClusterID"

    } else if (assay == "RNA") {
      # for RNA
      cat("Start run clustering for RNA using cell-cell distances...\n")

      # identify cell clusters by calling Seurat function
      object[["rna_snn"]] <- FindNeighbors(object = object@misc[['RNA']][['dist']], k.param = graph.k.param, compute.SNN = graph.compute.SNN, prune.SNN = graph.prune.SNN, nn.eps = graph.nn.eps, verbose = FALSE, force.recalc = graph.force.recalc)$snn
      object <- FindClusters(object = object, resolution = resolution[1], graph.name = "rna_snn", modularity.fxn = cluster.modularity.fxn, initial.membership = cluster.initial.membership, node.sizes = cluster.node.sizes, algorithm = cluster.algorithm, n.start = cluster.n.start, n.iter = cluster.n.iter, random.seed = cluster.random.seed, group.singletons = cluster.group.singletons, temp.file.location = cluster.temp.file.location, edge.file.name = cluster.edge.file.name, verbose = FALSE)
      object[["rnaClusterID"]] <- Idents(object = object)
      object@misc[['RNA']][['cluster']] <- "rnaClusterID"

    } else {
      stop("Please provide correct assay name! choose from RNA, ADT, Joint, All")
    }
    return(object)
  } else {
    stop("Please provide a Seurat Object!")
  }
}

#' umapFromDistane.Seurat
#'
#' @importFrom umap umap umap.defaults
#' @rdname umapFromDistane
#' @export
#'
umapFromDistane.Seurat <- function(
  object,
  assay = "All",
  seed = 42,
  method = "umap-learn",
  n.neighbors = 15,
  n.components = 2,
  metric = "euclidean",
  verbose = FALSE,
  n.epochs = 200,
  min.dist = 0.1,
  spread = 1,
  set.op.mix.ratio = 1,
  local.connectivity = 1L,
  negative.sample.rate = 5L
) {
  if(!is.null(object)){
    set.seed(seed = seed)

    my.umap.conf <- umap.defaults
    my.umap.conf$input <- "dist"
    my.umap.conf$n_neighbors <- n.neighbors
    my.umap.conf$n_components <- n.components
    my.umap.conf$metric <- metric
    my.umap.conf$verbose <- verbose
    my.umap.conf$n_epochs <- n.epochs
    my.umap.conf$min_dist <- min.dist
    my.umap.conf$spread <- spread
    my.umap.conf$set_op_mix_ratio <- set.op.mix.ratio
    my.umap.conf$local_connectivity <- local.connectivity
    my.umap.conf$negative_sample_rate <- negative.sample.rate

    if(assay == "All") {
      # for RNA
      cat("Start run UMAP for RNA using cell-cell distances...\n")

      rna.dist.matrix <- as.matrix(object@misc[['RNA']][['dist']])
      # run umap with pairwise distances

      my.umap <- umap(rna.dist.matrix,my.umap.conf,method = method)
      colnames(my.umap$layout) <- c('rnaUMAP_1','rnaUMAP_2')
      umap.reduction <- CreateDimReducObject(embeddings = my.umap$layout, key = "rnaUMAP_", assay = "RNA")
      object[["umap_rna"]] <- umap.reduction
      rownames(object@reductions[["umap_rna"]]@cell.embeddings) <- rownames(object@reductions[["pca"]]@cell.embeddings)
      object@misc[['RNA']][['reduction']] <- "umap_rna"

      # for ADT
      cat("Start run UMAP for ADT using cell-cell distances...\n")

      adt.dist.matrix <- as.matrix(object@misc[['ADT']][['dist']])
      # run umap with pairwise distances
      my.umap <- umap(adt.dist.matrix, my.umap.conf, method = method)
      colnames(my.umap$layout) <- c('adtUMAP_1','adtUMAP_2')
      umap.reduction <- CreateDimReducObject(embeddings = my.umap$layout, key = "adtUMAP_", assay = "ADT")
      object[["umap_adt"]] <- umap.reduction
      rownames(object@reductions[["umap_adt"]]@cell.embeddings) <- rownames(object@reductions[["pca"]]@cell.embeddings)
      object@misc[['ADT']][['reduction']] <- "umap_adt"

      # for Joint
      cat("Start run UMAP for Joint data using cell-cell distances...\n")

      joint.dist.matrix <- as.matrix(object@misc[['Joint']][['dist']])
      my.umap <- umap(joint.dist.matrix,my.umap.conf,method = method)
      colnames(my.umap$layout) <- c('jointUMAP_1','jointUMAP_2')
      umap.reduction <- CreateDimReducObject(embeddings = my.umap$layout,key = "jointUMAP_",assay = "ADT")
      object[["umap_joint"]] <- umap.reduction
      rownames(object@reductions[["umap_joint"]]@cell.embeddings) <- rownames(object@reductions[["pca"]]@cell.embeddings)
      object@misc[['Joint']][['reduction']] <- "umap_joint"
    } else if (assay == "Joint") {
      # for Joint
      cat("Start run UMAP for Joint data using cell-cell distances...\n")

      joint.dist.matrix <- as.matrix(object@misc[['Joint']][['dist']])
      my.umap <- umap(joint.dist.matrix,my.umap.conf,method = method)
      colnames(my.umap$layout) <- c('jointUMAP_1','jointUMAP_2')
      umap.reduction <- CreateDimReducObject(embeddings = my.umap$layout,key = "jointUMAP_",assay = "ADT")
      object[["umap_joint"]] <- umap.reduction
      rownames(object@reductions[["umap_joint"]]@cell.embeddings) <- rownames(object@reductions[["pca"]]@cell.embeddings)
      object@misc[['Joint']][['reduction']] <- "umap_joint"
    } else if (assay == "ADT") {
      # for ADT
      cat("Start run UMAP for ADT using cell-cell distances...\n")

      adt.dist.matrix <- as.matrix(object@misc[['ADT']][['dist']])
      # run umap with pairwise distances
      my.umap <- umap(adt.dist.matrix, my.umap.conf, method = method)
      colnames(my.umap$layout) <- c('adtUMAP_1','adtUMAP_2')
      umap.reduction <- CreateDimReducObject(embeddings = my.umap$layout, key = "adtUMAP_", assay = "ADT")
      object[["umap_adt"]] <- umap.reduction
      rownames(object@reductions[["umap_adt"]]@cell.embeddings) <- rownames(object@reductions[["pca"]]@cell.embeddings)
      object@misc[['ADT']][['reduction']] <- "umap_adt"
    } else if (assay == "RNA") {
      # for RNA
      cat("Start run UMAP for RNA using cell-cell distances...\n")

      rna.dist.matrix <- as.matrix(object@misc[['RNA']][['dist']])
      # run umap with pairwise distances

      my.umap <- umap(rna.dist.matrix,my.umap.conf,method = method)
      colnames(my.umap$layout) <- c('rnaUMAP_1','rnaUMAP_2')
      umap.reduction <- CreateDimReducObject(embeddings = my.umap$layout, key = "rnaUMAP_", assay = "RNA")
      object[["umap_rna"]] <- umap.reduction
      rownames(object@reductions[["umap_rna"]]@cell.embeddings) <- rownames(object@reductions[["pca"]]@cell.embeddings)
      object@misc[['RNA']][['reduction']] <- "umap_rna"
    } else {
      stop("Please provide correct assay name! choose from RNA, ADT, Joint, All")
    }
    return(object)
  } else {
    stop("Please provide a Seurat Object!")
  }
}

#' tsneFromDistane.Seurat
#'
#' @importFrom Rtsne Rtsne
#' @rdname tsneFromDistane
#' @export
#'

tsneFromDistane.Seurat <- function(
  object,
  assay = "All",
  perplexity = 30,
  dim = 2,
  seed = 42,
  theta = 0.5
) {
  if(!is.null(object)){
    set.seed(seed = seed)
    if(assay == "All") {
      # for RNA
      cat("Start run t-SNE for RNA using cell-cell distances...\n")

      rna.dist <- object@misc[['RNA']][['dist']]
      # run tsne with pairwise distances
      my.tsne <- Rtsne(rna.dist, perplexity = perplexity, is_distance = TRUE, dims = dim, theta = theta)

      colnames(my.tsne$Y) <- c("rnatsne_1", "rnatsne_2")
      tsne.reduction <- CreateDimReducObject(embeddings = my.tsne$Y, key = "rnatsne_", assay = "RNA")
      object[["tsne_rna"]] <- tsne.reduction
      rownames(object@reductions[["tsne_rna"]]@cell.embeddings) <- rownames(object@reductions[["pca"]]@cell.embeddings)
      object@misc[['RNA']][['reduction']] <- "tsne_rna"

      # for ADT
      cat("Start run t-SNE for ADT using cell-cell distances...\n")

      adt.dist <- object@misc[['ADT']][['dist']]
      # run tsne with pairwise distances
      my.tsne <- Rtsne(adt.dist, perplexity = perplexity, is_distance = TRUE, dims = dim, theta = theta)

      colnames(my.tsne$Y) <- c("adttsne_1", "adttsne_2")
      tsne.reduction <- CreateDimReducObject(embeddings = my.tsne$Y, key = "adttsne_", assay = "ADT")
      object[["tsne_adt"]] <- tsne.reduction
      rownames(object@reductions[["tsne_adt"]]@cell.embeddings) <- rownames(object@reductions[["pca"]]@cell.embeddings)
      object@misc[['ADT']][['reduction']] <- "tsne_adt"

      # for Joint
      cat("Start run t-SNE for Joint data using cell-cell distances...\n")

      joint.dist <- object@misc[['Joint']][['dist']]
      # run tsne with pairwise distances
      my.tsne <- Rtsne(joint.dist, perplexity = perplexity, is_distance = TRUE, dims = dim, theta = theta)

      colnames(my.tsne$Y) <- c("jointtsne_1", "jointtsne_2")
      tsne.reduction <- CreateDimReducObject(embeddings = my.tsne$Y, key = "jointtsne_", assay = "ADT")
      object[["tsne_joint"]] <- tsne.reduction
      rownames(object@reductions[["tsne_joint"]]@cell.embeddings) <- rownames(object@reductions[["pca"]]@cell.embeddings)
      object@misc[['Joint']][['reduction']] <- "tsne_joint"

    } else if (assay == "Joint") {
      # for Joint
      cat("Start run t-SNE for Joint data using cell-cell distances...\n")

      joint.dist <- object@misc[['Joint']][['dist']]
      # run tsne with pairwise distances
      my.tsne <- Rtsne(joint.dist, perplexity = perplexity, is_distance = TRUE, dims = dim, theta = theta)

      colnames(my.tsne$Y) <- c("jointtsne_1", "jointtsne_2")
      tsne.reduction <- CreateDimReducObject(embeddings = my.tsne$Y, key = "jointtsne_", assay = "ADT")
      object[["tsne_joint"]] <- tsne.reduction
      rownames(object@reductions[["tsne_joint"]]@cell.embeddings) <- rownames(object@reductions[["pca"]]@cell.embeddings)
      object@misc[['Joint']][['reduction']] <- "tsne_joint"
    } else if (assay == "ADT") {
      # for ADT
      cat("Start run t-SNE for ADT using cell-cell distances...\n")

      adt.dist <- object@misc[['ADT']][['dist']]
      # run tsne with pairwise distances
      my.tsne <- Rtsne(adt.dist, perplexity = perplexity, is_distance = TRUE, dims = dim, theta = theta)

      colnames(my.tsne$Y) <- c("adttsne_1", "adttsne_2")
      tsne.reduction <- CreateDimReducObject(embeddings = my.tsne$Y, key = "adttsne_", assay = "ADT")
      object[["tsne_adt"]] <- tsne.reduction
      rownames(object@reductions[["tsne_adt"]]@cell.embeddings) <- rownames(object@reductions[["pca"]]@cell.embeddings)
      object@misc[['ADT']][['reduction']] <- "tsne_adt"
      object@misc[['ADT']][['reduction']] <- "tsne_adt"
    } else if (assay == "RNA") {
      # for RNA
      cat("Start run t-SNE for RNA using cell-cell distances...\n")

      rna.dist <- object@misc[['RNA']][['dist']]
      # run tsne with pairwise distances
      my.tsne <- Rtsne(rna.dist, perplexity = perplexity, is_distance = TRUE, dims = dim, theta = theta)

      colnames(my.tsne$Y) <- c("rnatsne_1", "rnatsne_2")
      tsne.reduction <- CreateDimReducObject(embeddings = my.tsne$Y, key = "rnatsne_", assay = "RNA")
      object[["tsne_rna"]] <- tsne.reduction
      rownames(object@reductions[["tsne_rna"]]@cell.embeddings) <- rownames(object@reductions[["pca"]]@cell.embeddings)
      object@misc[['RNA']][['reduction']] <- "tsne_rna"
    } else {
      stop("Please provide correct assay name! choose from RNA, ADT, Joint, All")
    }
    return(object)
  } else {
    stop("Please provide a Seurat Object!")
  }
}


#' tsneFromDistane.default
#'
#' @importFrom Rtsne Rtsne
#' @rdname tsneFromDistane
#' @export
#'

tsneFromDistane.default <- function(
  object,
  perplexity = 30,
  dim = 2,
  seed = 42,
  theta = 0.5
) {
  if(!is.null(object)){
    set.seed(seed = seed)

    cat("Start run t-SNE for gaven cell-cell distances...\n")

    if(!is.matrix(object)) {
      dist <- as.matrix(object)
    }

    # run tsne with pairwise distances
    my.tsne <- Rtsne(dist, perplexity = perplexity, is_distance = TRUE, dims = dim, theta = theta)
    return(my.tsne)
  } else {
    stop("Please provide a dist object or dist matrix!")
  }
}


#' umapFromDistane.default
#'
#' @importFrom umap umap umap.defaults
#' @rdname umapFromDistane
#' @export
#'
umapFromDistane.default <- function(
  object,
  seed = 42,
  method = "umap-learn",
  n.neighbors = 15,
  n.components = 2,
  metric = "euclidean",
  verbose = TRUE,
  n.epochs = 200,
  min.dist = 0.1,
  spread = 1,
  set.op.mix.ratio = 1,
  local.connectivity = 1L,
  negative.sample.rate = 5L
) {
  if(!is.null(object)){
    set.seed(seed = seed)

    my.umap.conf <- umap.defaults
    my.umap.conf$input <- "dist"
    my.umap.conf$n_neighbors <- n.neighbors
    my.umap.conf$n_components <- n.components
    my.umap.conf$metric <- metric
    my.umap.conf$verbose <- verbose
    my.umap.conf$n_epochs <- n.epochs
    my.umap.conf$min_dist <- min.dist
    my.umap.conf$spread <- spread
    my.umap.conf$set_op_mix_ratio <- set.op.mix.ratio
    my.umap.conf$local_connectivity <- local.connectivity
    my.umap.conf$negative_sample_rate <- negative.sample.rate

    cat("Start run UMAP for gaven using cell-cell distances...\n")

    if(!is.matrix(object)) {
      dist <- as.matrix(object)
    }
    my.umap <- umap(dist,my.umap.conf,method = method)
    return(my.umap)
  } else {
    stop("Please provide a dist object or dist matrix!")
  }
}


#' jointDistance.Seurat
#'
#' @rdname jointDistance
#'
#' @export
jointDistance.Seurat <- function(
  object,
  dims = 20,
  beta = 0.5,
  model = "LP",
  keep.rna = TRUE,
  keep.adt = TRUE,
  sigmoid.n = 10,
  sigmoid.k = 0.5,
  precision.mode = TRUE,
  c = NULL
) {
  cat("Start working...\n")

  if(isTRUE(precision.mode)){
    adt.num <- dim(object@assays[["ADT"]]@data)[1]
    if(is.null(c)){
      c <- rep(1,adt.num)
    } else {
      if(length(c) != adt.num){
        stop("Number of C should be equal to the number of ADT features! Please check your input!")
      }
    }
  }

  # for RNA
  cat("Start calculate cell-cell pairwise distances for RNA...\n")
  DefaultAssay(object = object) <- "RNA"
  # calculate cell-cell pairwise distances for RNA using top n PCs
  rna.pca <- object@reductions[["pca"]]
  if(is.null(rna.pca)) {
    stop("Please do PCA analysis for RNA data first!")
  }
  rna.pca <- object@reductions[["pca"]]@cell.embeddings[,1:dims]
  rna.dist <- dist(x = rna.pca)

  # for ADT
  cat("Start calculate cell-cell pairwise distances for ADT...\n")
  DefaultAssay(object = object) <- "ADT"
  # calculate cell-cell pairwise distances for ADT directly using normalized data
  adt.data <- as.matrix(GetAssayData(object, slot = "data"))

  if(isTRUE(precision.mode)){
    adt.dist <- scaleDistUpdateCpp(adt.data, sigmoid.n, sigmoid.k, c)
  } else {
    adt.dist <- scaleDistCpp(adt.data, sigmoid.n, sigmoid.k)
  }
  rownames(adt.dist) <- colnames(adt.data)
  colnames(adt.dist) <- colnames(adt.data)
  adt.dist <- as.dist(adt.dist)

  # Joint cell-cell distance
  cat("Start calculate joint cell-cell pairwise distances... \n")
  if(!is.null(model)) {
    cat("Scale ADT and RNA distances into same level... \n")
    learning.rate <- 1
    low.threshold <- 1e-5
    max.iter <- 10000
    X <- rna.dist + adt.dist
    Y <- adt.dist
    #optimal <- gradientDescent(X, Y, alpha = 0, learning.rate, low.threshold)
    optimal <- gradientDescentCpp(X, Y, alpha = 0, learning.rate, low.threshold, max.iter)
    rna.dist.scale <- rna.dist*optimal
    adt.dist.scale <- adt.dist*(1-optimal)

    if(isTRUE(keep.rna)) {
      if(!is.null(object@misc[['RNA']])) {
        object@misc[['RNA']][['dist']] <- rna.dist.scale
      } else {
        object@misc[['RNA']] <- list()
        object@misc[['RNA']][['dist']] <- rna.dist.scale
      }
    }

    if(isTRUE(keep.adt)) {
      if(!is.null(object@misc[['ADT']])) {
        object@misc[['ADT']][['dist']] <- adt.dist.scale
      } else {
        object@misc[['ADT']] <- list()
        object@misc[['ADT']][['dist']] <- adt.dist.scale
      }
    }

    if(model == "LP") {
      cat("Calculate joint distance using L-infinite model ... \n")

      joint.dist <- pmax(rna.dist.scale,adt.dist.scale)

      c.rna <- joint.dist - rna.dist.scale
      c.adt <- joint.dist - adt.dist.scale
      n.rna <- length(which(c.rna == 0))
      n.adt <- length(which(c.adt == 0))

      total <- n.rna + n.adt
      p.rna <- n.rna/total*100
      p.adt <- n.adt/total*100

      cat("Contribution of ADT is",p.adt,"% \n")
      cat("Contribution of RNA is",p.rna,"% \n")

    } else if (model == "L1") {
      cat("Calculate joint distance using L-1 model ... \n")
      if(!is.numeric(beta)) {
        stop("beta should be a number!\n")
      } else {
        joint.dist <- rna.dist.scale*beta + adt.dist.scale*(1-beta)
        p.rna <- 1
        p.adt <- 1
      }
    } else {
      stop("Can not recognize your model! Please set model, 'LP' or 'L1'! \n")
    }

    if(!is.null(object@misc[['Joint']])) {
      object@misc[['Joint']][['dist']] <- joint.dist
      object@misc[['Joint']][['alpha']] <- optimal
      object@misc[['Joint']][['model']] <- model
      object@misc[['Joint']][['contribution']][['rna']] <- p.rna
      object@misc[['Joint']][['contribution']][['adt']] <- p.adt
    } else {
      object@misc[['Joint']] <- list()
      object@misc[['Joint']][['dist']] <- joint.dist
      object@misc[['Joint']][['alpha']] <- optimal
      object@misc[['Joint']][['model']] <- model
      object@misc[['Joint']][['contribution']][['rna']] <- p.rna
      object@misc[['Joint']][['contribution']][['adt']] <- p.adt
    }
  } else {
    stop("Please set model, 'LP' or 'L1'! \n")
  }
  return(object)
}


#' jointDistance.default
#'
#' @rdname jointDistance
#' @export
#'
jointDistance.default <- function(
  object,
  dist1 = NULL,
  dist2 = NULL,
  beta = 1,
  model = "LP"
) {
  cat("Start calculate joint cell-cell pairwise distances... \n")
  if(!is.null(model)) {
    cat("Scale two distances into same level... \n")
    learning.rate <- 1
    low.threshold <- 1e-5
    max.iter <- 10000
    X <- dist1 + dist2
    Y <- dist2
    #optimal <- gradientDescent(X, Y, alpha = 0, learning.rate, low.threshold)
    optimal <- gradientDescentCpp(X, Y, alpha = 0, learning.rate, low.threshold, max.iter)
    dist1.scale <- dist1*optimal
    dist2.scale <- dist2*(1-optimal)
    if(model == "LP") {
      cat("Calculate joint distance using L-infinite model ... \n")

      joint.dist <- pmax(dist1.scale,dist2.scale)

      c1 <- joint.dist - dist1.scale
      c2 <- joint.dist - dist2.scale
      n1 <- length(which(c1 == 0))
      n2 <- length(which(c2 == 0))

      total <- n1 + n2
      p1 <- n1/total*100
      p2 <- n2/total*100

      cat("Contribution of dist1 is",p2,"% \n")
      cat("Contribution of dist2 is",p1,"% \n")

    } else if (model == "L1") {
      cat("Calculate joint distance using L-1 model ... \n")
      if(!is.numeric(beta)) {
        stop("beta should be a number!\n")
      } else {
        joint.dist <- dist1.scale*beta + dist2.scale*(1-beta)
      }
    } else {
      stop("Can not recognize your model! Please set model, 'LP' or 'L1'! \n")
    }
  } else {
    stop("Please set model, 'LP' or 'L1'! \n")
  }
  return(joint.dist)
}

#' errorFunction
#'
#' objective function value
#'
#' @param alpha current alpha
#' @param X dist of X
#' @param Y dist of Y
#'

errorFunction <- function(alpha, X, Y) {
  diff <- sum((alpha*X - Y)^2)/(2*length(Y))
  return(diff)
}

#' gradientFunction
#'
#' gradient value
#'
#' @param alpha current alpha
#' @param X dist of X
#' @param Y dist of Y
#'

gradientFunction <- function(alpha, X, Y) {
  diff <- sum(alpha*X - Y)/length(Y)
  return(diff)
}

#' gradientDescent
#'
#' gradient descent algorithm
#'
#' @param alpha starting alpha
#' @param X dist of X
#' @param Y dist of Y
#' @param learning.rate step size of each iteration
#' @param low.threshold low threshold of gradient
#'
#' @export
gradientDescent <- function(X, Y, alpha, learning.rate, low.threshold) {
  gradient <- gradientFunction(alpha, X, Y)
  iter <- 1
  alpha.last <- alpha
  while (abs(gradient) > low.threshold) {
    alpha <- alpha - learning.rate*gradient
    gradient <- gradientFunction(alpha, X, Y)
    #cat("iter = ",iter, ",current alpha = ", alpha, ",current gradiant = ", gradient, ",current learning rate = ", learning.rate,"\n")
    iter <- iter + 1
    if((alpha < 0)||(alpha > 1)) {
      # the learning rate is too high, reduce the learning rate
      #cat("the learning rate is too high, reducing the learning rate\n")
      if(abs(alpha.last) < abs(alpha)) {
        learning.rate <- learning.rate/2
        alpha.last <- alpha
        #reset alpha
        alpha <- 0
      }
    }
  }
  return(alpha)
}


#' scoreRNA
#'
#' RNA score function. use ROGUE score for RNA data
#'
#' @importFrom ROGUE rogue
#'
#' @param object single cell data object
#' @param group.by name of clustering
#' @param platform platform of rna sequencing
#' @param span parameter span
#' @param samples name of sample information, can be blank
#'
#' @export
#'
scoreRNA <- function(
  object = NULL,
  group.by = NULL,
  platform = 'UMI',
  span = 0.6,
  samples = NULL
){
  if(is.null(object)) {
    stop('Please provide data object!')
  }
  if(is.null(group.by)) {
    stop('Please provide clustering name! e.g. jointClusterID, rnaClusterID...')
  }
  if(is.null(samples)) {
    samples <- rep('EntireDataset', length(object@meta.data[[group.by]]))
  } else {
    samples <- object@meta.data[[samples]]
  }

  rna_expr <- object@assays[["RNA"]]@counts

  rogue.res.rna <- rogue(rna_expr, labels = object@meta.data[[group.by]], samples = samples, platform = "UMI", span = 0.6)
  rogue.res.rna <- rogue.res.rna[,sort(colnames(rogue.res.rna))]
  return(rogue.res.rna)
}


#' scoreADT
#'
#' ADT score function
#'
#'
#' @param object single cell data object
#' @param group.by name of clustering
#' @param samples name of sample information, can be blank
#' @param k parameter k option. set "auto" will use the sum of SDs of entire dataset
#'
#' @export

scoreADT <- function(
  object = NULL,
  group.by = NULL,
  k = 'auto',
  samples = NULL,
  rank = 2
){
  if(is.null(object)) {
    stop('Please provide data object!')
  }
  if(is.null(group.by)) {
    stop('Please provide clustering name! e.g. jointClusterID, rnaClusterID...')
  }
  adt_expr_norm <- object@assays[["ADT"]]@data
  if(!is.numeric(k)){
    k = 0
    for (i in c(1:dim(adt_expr_norm)[1])) {
      value <- sd(adt_expr_norm[i,])
      k <- k + value^rank
    }
  }
  if(is.null(samples)) {
    adt.score <- scoreADTfun(adt_expr_norm, labels = object@meta.data[[group.by]], k = k, rank = rank)
  } else {
    adt.score <- NULL
    samples <- as.character(object@meta.data[[samples]])
    samples.name <- unique(samples)
    labels <- as.character(object@meta.data[[group.by]])
    labels.name <- unique(labels)
    for (sample.name in samples.name) {
      cur_index <- which(samples == sample.name)
      cur_adt_expr_norm <- adt_expr_norm[,cur_index]
      cur_labels <- object@meta.data[[group.by]]
      cur_labels <- cur_labels[cur_index]
      cur_adt.score <- scoreADTfun(cur_adt_expr_norm, labels = cur_labels,labels.name = labels.name, k = k, rank = rank)
      rownames(cur_adt.score) <- c(sample.name)
      if(is.null(adt.score)) {
        adt.score <- cur_adt.score
      } else {
        adt.score <- rbind(adt.score, cur_adt.score)
      }
    }
  }
  adt.score <- adt.score[,sort(colnames(adt.score))]
  return(adt.score)
}

#' scoreADTfun
#'
#' ADT score function
#'
#' @param object single cell data object
#' @param labels name of clustering
#' @param labels.name a list of all labels
#' @param k parameter k
#'
#'

scoreADTfun <- function(
  object = NULL,
  labels = NULL,
  labels.name = NULL,
  k = 10,
  rank = 2
) {
  # determine k value based on current dataset
  if(!is.numeric(k)){
    k = 0
    adt_expr <- object
    for (i in c(1:dim(adt_expr)[1])) {
      value <- sd(adt_expr[i,])
      k <- k + value^rank
    }
  }

  labels <- as.character(labels)
  if(is.null(labels.name)) {
    labels.name <- unique(labels)
  }
  scores <- c()
  # for each cluster
  for(label.name in labels.name){
    index <- which(labels == label.name)
    if(length(index) == 0){
      cat('Current sample do not have members in this cluster, set to 0...')
      plus.sd <- 0
      scores <- c(scores, plus.sd)
    } else {
      adt_expr <- object[,index]
      plus.sd <- 0
      for (i in c(1:dim(adt_expr)[1])) {
        value <- sd(adt_expr[i,])
        plus.sd <- plus.sd + value^rank
      }
      plus.sd <- k/(plus.sd + k)
      scores <- c(scores, plus.sd)
    }
  }
  names(scores) <- labels.name
  scores <- t(as.data.frame(scores))
  scores <- as.data.frame(scores)
  return(scores)
}

#' CombineScore
#'
#' Combine RNA and ADT scores
#'
#' @param rna.score RNA score, output from scoreRNA()
#' @param adt.score ADT score, output from scoreADT()
#'
#' @export

CombineScore <- function(
  rna.score = NULL,
  adt.score = NULL
) {
  rna.score <- rna.score[,sort(colnames(rna.score))]
  adt.score <- adt.score[,sort(colnames(adt.score))]
  combine.score <- c()
  if (length(rna.score) == length(adt.score)) {
    for (i in c(1:length(rna.score))) {
      current.score <- sqrt(rna.score[i]*adt.score[i])
      combine.score <- c(combine.score, current.score)
    }
    combine.score <- as.data.frame(combine.score)
    colnames(combine.score) <- colnames(rna.score)
    rownames(combine.score) <- c('Purity score')
    return(combine.score)
  } else {
    stop('Your ADT score and RNA score do not have same length! Check your input!')
  }
}
