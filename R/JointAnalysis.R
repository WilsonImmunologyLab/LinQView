#' jointClusteing
#' 1) Calculate cell-cell pairwise distances for RNA and ADT separately, then calculate the joint cell-cell pairwise distances
#' 2) Run UMAP to project single cells into 2D map using cell-cell pairwise distances
#' 3) Run clustering method (call Seurat function) to identify cell populations using cell-cell pairwise distances
#' Since trajectory analysis require distances, this function is not recommended
#'
#' @param object Seurat object
#' @param run_seprate also calculate cell-cell pairwise distances, umap and clustering for RNA and ADT separately
#' @param reduction_method reduction method. Now only support UMAP. Will support t-SNE later...
#' @param dims number of PCs used for RNA data. Default is top 20
#' @param resolution resolution for 1) Joint, 2) RNA and 3) ADT clustering. if run_seprate = TRUE, user should provide all three solutions, otherwise will use default value 0.2 for all three datasets
#' @param alpha use alpha to balence contributions from RNA and ADT in the joint distances. We suggest use 0.5 for initial analysis, and then adjust alpha for better results
#'
#' @export
jointClusteing <- function(
  object,
  run_seprate = FALSE,
  reduction_method = "UMAP",
  dims = 20,
  resolution = c(0.2,0.2,0.2),
  alpha = 0.5
) {
  cat("Start working...\n")
  if(isTRUE(run_seprate)) {
    if(length(resolution < 3)) {
      cat("You set run_seprate = TRUE, but didn't provide 3 resolutions, will use default value 0.2 for all three \n")
      resolution = c(0.2,0.2,0.2)
    }
  }

  # for RNA
  cat("Start calculate cell-cell pairwise distances for RNA...\n")
  DefaultAssay(object = object) <- "RNA"
  # calculate cell-cell pairwise distances for RNA using top n PCs
  rna.pca <- object@reductions[["pca"]]@cell.embeddings[,1:dims]
  rna.dist <- dist(x = rna.pca)
  object@misc[['rnaDist']] <- rna.dist

  if(isTRUE(run_seprate)) {
    rna.dist.matrix <- as.matrix(rna.dist)
    # run umap with pairwise distances
    my.umap.conf <- umap.defaults
    my.umap.conf$input <- "dist"
    my.umap <- umap(rna.dist.matrix,my.umap.conf,method = "umap-learn")
    umap.reduction <- CreateDimReducObject(embeddings = my.umap$layout, key = "rnaUMAP_", assay = "RNA")
    object[["umap_rna"]] <- umap.reduction
    rownames(object@reductions[["umap_rna"]]@cell.embeddings) <- rownames(object@reductions[["pca"]]@cell.embeddings)
    #colnames(object@reductions[["umap_rna"]]@cell.embeddings) <- c("UMAP1","UMAP2")

    # identify cell clusters by calling Seurat function
    object[["rna_snn"]] <- FindNeighbors(object = rna.dist)$snn
    object <- FindClusters(object = object, resolution = resolution[2], graph.name = "rna_snn")
    object[["rnaClusterID"]] <- Idents(object = object)
  }

  # for ADT
  cat("Start calculate cell-cell pairwise distances for ADT...\n")
  DefaultAssay(object = object) <- "ADT"
  # calculate cell-cell pairwise distances for ADT directly using normalized data
  adt.data <- t(as.matrix(GetAssayData(object, slot = "data")))
  adt.dist <- dist(x = adt.data)
  object@misc[['adtDist']] <- adt.dist

  if(isTRUE(run_seprate)) {
    adt.dist.matrix <- as.matrix(adt.dist)
    # run umap with pairwise distances
    my.umap.conf <- umap.defaults
    my.umap.conf$input <- "dist"
    my.umap <- umap(adt.dist.matrix, my.umap.conf, method = "umap-learn")
    umap.reduction <- CreateDimReducObject(embeddings = my.umap$layout, key = "adtUMAP_", assay = "ADT")
    object[["umap_adt"]] <- umap.reduction
    rownames(object@reductions[["umap_adt"]]@cell.embeddings) <- rownames(object@reductions[["pca"]]@cell.embeddings)
    #colnames(object@reductions[["umap_adt"]]@cell.embeddings) <- c("UMAP1","UMAP2")

    # identify cell clusters by calling Seurat function
    object[["adt_snn"]] <- FindNeighbors(object = adt.dist)$snn
    object <- FindClusters(object = object, resolution = resolution[3], graph.name = "adt_snn")
    object[["adtClusterID"]] <- Idents(object = object)
  }

  # Joint cell-cell distance
  cat("Start calculate joint cell-cell pairwise distances... alpha = ",alpha,"\n")
  # calculate scale factor to scale two distances into the same level
  scale.factor <- sum(adt.dist)/sum(rna.dist)
  joint.dist <- adt.dist*alpha + rna.dist*(1-alpha)*scale.factor
  object@misc[['jointDist']] <- joint.dist

  joint.dist.matrix <- as.matrix(joint.dist)
  my.umap.conf <- umap.defaults
  my.umap.conf$input <- "dist"
  my.umap <- umap(joint.dist.matrix,my.umap.conf,method = "umap-learn")
  umap.reduction <- CreateDimReducObject(embeddings = my.umap$layout,key = "jointUMAP_",assay = "ADT")
  object[["umap_joint"]] <- umap.reduction
  rownames(object@reductions[["umap_joint"]]@cell.embeddings) <- rownames(object@reductions[["pca"]]@cell.embeddings)
  #colnames(object@reductions[["umap_joint"]]@cell.embeddings) <- c("UMAP1","UMAP2")

  # identify cell clusters by calling Seurat function
  object[["joint_snn"]] <- FindNeighbors(object = joint.dist)$snn
  object <- FindClusters(object = object, resolution = resolution[1], graph.name = "joint_snn")
  object[["jointClusterID"]] <- Idents(object = object)

  return(object)
}


#' clusteringFromDistance
#'
#' Run clustering method (implemented by Seurat package) to identify cell populations using cell-cell pairwise distances
#'
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
#' @param cluster.weights parameter for Seurat function FindClusters.Parameters to pass to the Python leidenalg function.
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
  cluster.weights = NULL,
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
      cat("You set resolution for RNA = ",resolution[2],", For ADT = ",resolution[3],", For joint = ",resolution[1], "...\n ")
    }
  }
  if(!is.null(object)){
    if(assay == "All") {
      # for RNA
      cat("Start run clustering for RNA using cell-cell distances...\n")

      # identify cell clusters by calling Seurat function
      object[["rna_snn"]] <- FindNeighbors(object = object@misc[['rnaDist']], k.param = graph.k.param, compute.SNN = graph.compute.SNN, prune.SNN = graph.prune.SNN, nn.eps = graph.nn.eps, verbose = FALSE, force.recalc = graph.force.recalc)$snn
      object <- FindClusters(object = object, resolution = resolution[2], graph.name = "rna_snn", modularity.fxn = cluster.modularity.fxn, initial.membership = cluster.initial.membership, weights = cluster.weights, node.sizes = cluster.node.sizes, algorithm = cluster.algorithm, n.start = cluster.n.start, n.iter = cluster.n.iter, random.seed = cluster.random.seed, group.singletons = cluster.group.singletons, temp.file.location = cluster.temp.file.location, edge.file.name = cluster.edge.file.name, verbose = FALSE)
      object[["rnaClusterID"]] <- Idents(object = object)

      # for ADT
      cat("Start run clustering for ADT using cell-cell distances...\n")

      # identify cell clusters by calling Seurat function
      object[["adt_snn"]] <- FindNeighbors(object = object@misc[['adtDist']], k.param = graph.k.param, compute.SNN = graph.compute.SNN, prune.SNN = graph.prune.SNN, nn.eps = graph.nn.eps, verbose = FALSE, force.recalc = graph.force.recalc)$snn
      object <- FindClusters(object = object, resolution = resolution[3], graph.name = "adt_snn", modularity.fxn = cluster.modularity.fxn, initial.membership = cluster.initial.membership, weights = cluster.weights, node.sizes = cluster.node.sizes, algorithm = cluster.algorithm, n.start = cluster.n.start, n.iter = cluster.n.iter, random.seed = cluster.random.seed, group.singletons = cluster.group.singletons, temp.file.location = cluster.temp.file.location, edge.file.name = cluster.edge.file.name, verbose = FALSE)
      object[["adtClusterID"]] <- Idents(object = object)

      # for Joint
      cat("Start run clustering for Joint data using cell-cell distances...\n")

      # identify cell clusters by calling Seurat function
      object[["joint_snn"]] <- FindNeighbors(object = object@misc[['jointDist']], k.param = graph.k.param, compute.SNN = graph.compute.SNN, prune.SNN = graph.prune.SNN, nn.eps = graph.nn.eps, verbose = FALSE, force.recalc = graph.force.recalc)$snn
      object <- FindClusters(object = object, resolution = resolution[1], graph.name = "joint_snn", modularity.fxn = cluster.modularity.fxn, initial.membership = cluster.initial.membership, weights = cluster.weights, node.sizes = cluster.node.sizes, algorithm = cluster.algorithm, n.start = cluster.n.start, n.iter = cluster.n.iter, random.seed = cluster.random.seed, group.singletons = cluster.group.singletons, temp.file.location = cluster.temp.file.location, edge.file.name = cluster.edge.file.name, verbose = FALSE)
      object[["jointClusterID"]] <- Idents(object = object)

    } else if (assay == "Joint") {
      # for Joint
      cat("Start run clustering for Joint data using cell-cell distances...\n")

      # identify cell clusters by calling Seurat function
      object[["joint_snn"]] <- FindNeighbors(object = object@misc[['jointDist']], k.param = graph.k.param, compute.SNN = graph.compute.SNN, prune.SNN = graph.prune.SNN, nn.eps = graph.nn.eps, verbose = FALSE, force.recalc = graph.force.recalc)$snn
      object <- FindClusters(object = object, resolution = resolution[1], graph.name = "joint_snn", modularity.fxn = cluster.modularity.fxn, initial.membership = cluster.initial.membership, weights = cluster.weights, node.sizes = cluster.node.sizes, algorithm = cluster.algorithm, n.start = cluster.n.start, n.iter = cluster.n.iter, random.seed = cluster.random.seed, group.singletons = cluster.group.singletons, temp.file.location = cluster.temp.file.location, edge.file.name = cluster.edge.file.name, verbose = FALSE)
      object[["jointClusterID"]] <- Idents(object = object)

    } else if (assay == "ADT") {
      # for ADT
      cat("Start run clustering for ADT using cell-cell distances...\n")

      # identify cell clusters by calling Seurat function
      object[["adt_snn"]] <- FindNeighbors(object = object@misc[['adtDist']], k.param = graph.k.param, compute.SNN = graph.compute.SNN, prune.SNN = graph.prune.SNN, nn.eps = graph.nn.eps, verbose = FALSE, force.recalc = graph.force.recalc)$snn
      object <- FindClusters(object = object, resolution = resolution[1], graph.name = "adt_snn", modularity.fxn = cluster.modularity.fxn, initial.membership = cluster.initial.membership, weights = cluster.weights, node.sizes = cluster.node.sizes, algorithm = cluster.algorithm, n.start = cluster.n.start, n.iter = cluster.n.iter, random.seed = cluster.random.seed, group.singletons = cluster.group.singletons, temp.file.location = cluster.temp.file.location, edge.file.name = cluster.edge.file.name, verbose = FALSE)
      object[["adtClusterID"]] <- Idents(object = object)

    } else if (assay == "RNA") {
      # for RNA
      cat("Start run clustering for RNA using cell-cell distances...\n")

      # identify cell clusters by calling Seurat function
      object[["rna_snn"]] <- FindNeighbors(object = object@misc[['rnaDist']], k.param = graph.k.param, compute.SNN = graph.compute.SNN, prune.SNN = graph.prune.SNN, nn.eps = graph.nn.eps, verbose = FALSE, force.recalc = graph.force.recalc)$snn
      object <- FindClusters(object = object, resolution = resolution[1], graph.name = "rna_snn", modularity.fxn = cluster.modularity.fxn, initial.membership = cluster.initial.membership, weights = cluster.weights, node.sizes = cluster.node.sizes, algorithm = cluster.algorithm, n.start = cluster.n.start, n.iter = cluster.n.iter, random.seed = cluster.random.seed, group.singletons = cluster.group.singletons, temp.file.location = cluster.temp.file.location, edge.file.name = cluster.edge.file.name, verbose = FALSE)
      object[["rnaClusterID"]] <- Idents(object = object)

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

      rna.dist.matrix <- as.matrix(object@misc[['rnaDist']])
      # run umap with pairwise distances

      my.umap <- umap(rna.dist.matrix,my.umap.conf,method = method)
      umap.reduction <- CreateDimReducObject(embeddings = my.umap$layout, key = "rnaUMAP_", assay = "RNA")
      object[["umap_rna"]] <- umap.reduction
      rownames(object@reductions[["umap_rna"]]@cell.embeddings) <- rownames(object@reductions[["pca"]]@cell.embeddings)
      #colnames(object@reductions[["umap_rna"]]@cell.embeddings) <- c("UMAP1","UMAP2")

      # for ADT
      cat("Start run UMAP for ADT using cell-cell distances...\n")

      adt.dist.matrix <- as.matrix(object@misc[['adtDist']])
      # run umap with pairwise distances
      my.umap <- umap(adt.dist.matrix, my.umap.conf, method = method)
      umap.reduction <- CreateDimReducObject(embeddings = my.umap$layout, key = "adtUMAP_", assay = "ADT")
      object[["umap_adt"]] <- umap.reduction
      rownames(object@reductions[["umap_adt"]]@cell.embeddings) <- rownames(object@reductions[["pca"]]@cell.embeddings)
      #colnames(object@reductions[["umap_adt"]]@cell.embeddings) <- c("UMAP1","UMAP2")

      # for Joint
      cat("Start run UMAP for Joint data using cell-cell distances...\n")

      joint.dist.matrix <- as.matrix(object@misc[['jointDist']])
      my.umap <- umap(joint.dist.matrix,my.umap.conf,method = method)
      umap.reduction <- CreateDimReducObject(embeddings = my.umap$layout,key = "jointUMAP_",assay = "ADT")
      object[["umap_joint"]] <- umap.reduction
      rownames(object@reductions[["umap_joint"]]@cell.embeddings) <- rownames(object@reductions[["pca"]]@cell.embeddings)
      #colnames(object@reductions[["umap_joint"]]@cell.embeddings) <- c("UMAP1","UMAP2")
    } else if (assay == "Joint") {
      # for Joint
      cat("Start run UMAP for Joint data using cell-cell distances...\n")

      joint.dist.matrix <- as.matrix(object@misc[['jointDist']])
      my.umap <- umap(joint.dist.matrix,my.umap.conf,method = method)
      umap.reduction <- CreateDimReducObject(embeddings = my.umap$layout,key = "jointUMAP_",assay = "ADT")
      object[["umap_joint"]] <- umap.reduction
      rownames(object@reductions[["umap_joint"]]@cell.embeddings) <- rownames(object@reductions[["pca"]]@cell.embeddings)
      #colnames(object@reductions[["umap_joint"]]@cell.embeddings) <- c("UMAP1","UMAP2")
    } else if (assay == "ADT") {
      # for ADT
      cat("Start run UMAP for ADT using cell-cell distances...\n")

      adt.dist.matrix <- as.matrix(object@misc[['adtDist']])
      # run umap with pairwise distances
      my.umap <- umap(adt.dist.matrix, my.umap.conf, method = method)
      umap.reduction <- CreateDimReducObject(embeddings = my.umap$layout, key = "adtUMAP_", assay = "ADT")
      object[["umap_adt"]] <- umap.reduction
      rownames(object@reductions[["umap_adt"]]@cell.embeddings) <- rownames(object@reductions[["pca"]]@cell.embeddings)
      #colnames(object@reductions[["umap_adt"]]@cell.embeddings) <- c("UMAP1","UMAP2")
    } else if (assay == "RNA") {
      # for RNA
      cat("Start run UMAP for RNA using cell-cell distances...\n")

      rna.dist.matrix <- as.matrix(object@misc[['rnaDist']])
      # run umap with pairwise distances
      my.umap <- umap(rna.dist.matrix,my.umap.conf,method = method)
      umap.reduction <- CreateDimReducObject(embeddings = my.umap$layout, key = "rnaUMAP_", assay = "RNA")
      object[["umap_rna"]] <- umap.reduction
      rownames(object@reductions[["umap_rna"]]@cell.embeddings) <- rownames(object@reductions[["pca"]]@cell.embeddings)
      #colnames(object@reductions[["umap_rna"]]@cell.embeddings) <- c("UMAP1","UMAP2")
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

      rna.dist <- object@misc[['rnaDist']]
      # run tsne with pairwise distances
      my.tsne <- Rtsne(rna.dist, perplexity = perplexity, is_distance = TRUE, dims = dim, theta = theta)

      tsne.reduction <- CreateDimReducObject(embeddings = my.tsne$Y, key = "rna_tsne_", assay = "RNA")
      object[["tsne_rna"]] <- tsne.reduction
      rownames(object@reductions[["tsne_rna"]]@cell.embeddings) <- rownames(object@reductions[["pca"]]@cell.embeddings)

      # for ADT
      cat("Start run t-SNE for ADT using cell-cell distances...\n")

      adt.dist <- object@misc[['adtDist']]
      # run tsne with pairwise distances
      my.tsne <- Rtsne(adt.dist, perplexity = perplexity, is_distance = TRUE, dims = dim, theta = theta)

      tsne.reduction <- CreateDimReducObject(embeddings = my.tsne$Y, key = "adt_tsne_", assay = "ADT")
      object[["tsne_adt"]] <- tsne.reduction
      rownames(object@reductions[["tsne_adt"]]@cell.embeddings) <- rownames(object@reductions[["pca"]]@cell.embeddings)

      # for Joint
      cat("Start run t-SNE for Joint data using cell-cell distances...\n")

      joint.dist <- object@misc[['jointDist']]
      # run tsne with pairwise distances
      my.tsne <- Rtsne(joint.dist, perplexity = perplexity, is_distance = TRUE, dims = dim, theta = theta)

      tsne.reduction <- CreateDimReducObject(embeddings = my.tsne$Y, key = "joint_tsne_", assay = "ADT")
      object[["tsne_joint"]] <- tsne.reduction
      rownames(object@reductions[["tsne_joint"]]@cell.embeddings) <- rownames(object@reductions[["pca"]]@cell.embeddings)

    } else if (assay == "Joint") {
      # for Joint
      cat("Start run t-SNE for Joint data using cell-cell distances...\n")

      joint.dist <- object@misc[['jointDist']]
      # run tsne with pairwise distances
      my.tsne <- Rtsne(joint.dist, perplexity = perplexity, is_distance = TRUE, dims = dim, theta = theta)

      tsne.reduction <- CreateDimReducObject(embeddings = my.tsne$Y, key = "joint_tsne_", assay = "ADT")
      object[["tsne_joint"]] <- tsne.reduction
      rownames(object@reductions[["tsne_joint"]]@cell.embeddings) <- rownames(object@reductions[["pca"]]@cell.embeddings)
    } else if (assay == "ADT") {
      # for ADT
      cat("Start run t-SNE for ADT using cell-cell distances...\n")

      adt.dist <- object@misc[['adtDist']]
      # run tsne with pairwise distances
      my.tsne <- Rtsne(adt.dist, perplexity = perplexity, is_distance = TRUE, dims = dim, theta = theta)

      tsne.reduction <- CreateDimReducObject(embeddings = my.tsne$Y, key = "adt_tsne_", assay = "ADT")
      object[["tsne_adt"]] <- tsne.reduction
      rownames(object@reductions[["tsne_adt"]]@cell.embeddings) <- rownames(object@reductions[["pca"]]@cell.embeddings)
    } else if (assay == "RNA") {
      # for RNA
      cat("Start run t-SNE for RNA using cell-cell distances...\n")

      rna.dist <- object@misc[['rnaDist']]
      # run tsne with pairwise distances
      my.tsne <- Rtsne(rna.dist, perplexity = perplexity, is_distance = TRUE, dims = dim, theta = theta)

      tsne.reduction <- CreateDimReducObject(embeddings = my.tsne$Y, key = "rna_tsne_", assay = "RNA")
      object[["tsne_rna"]] <- tsne.reduction
      rownames(object@reductions[["tsne_rna"]]@cell.embeddings) <- rownames(object@reductions[["pca"]]@cell.embeddings)
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
#' @export
jointDistance.Seurat <- function(
  object,
  dims = 20,
  alpha = NULL
) {
  cat("Start working...\n")

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
  object@misc[['rnaDist']] <- rna.dist

  # for ADT
  cat("Start calculate cell-cell pairwise distances for ADT...\n")
  DefaultAssay(object = object) <- "ADT"
  # calculate cell-cell pairwise distances for ADT directly using normalized data
  adt.data <- t(as.matrix(GetAssayData(object, slot = "data")))
  adt.dist <- dist(x = adt.data)
  object@misc[['adtDist']] <- adt.dist

  # Joint cell-cell distance
  cat("Start calculate joint cell-cell pairwise distances... \n")
  if(is.null(alpha)) {
    cat("Will automatically determine an alpha to calculate the joint distance \n")
    learning.rate <- 1
    low.threshold <- 1e-5
    alpha <- 0
    X <- rna.dist + adt.dist
    Y <- adt.dist
    optimal <- gradientDescent(X, Y, alpha, learning.rate, low.threshold)
    cat("alpha is set to",optimal," \n")
    joint.dist <- rna.dist*optimal + adt.dist*(1-optimal)
    object@misc[['jointDist']] <- joint.dist
  } else {
    cat("will use user's alpha value",alpha," to calculate the joint distance \n")
    # calculate scale factor to scale two distances into the same level
    scale.factor <- sum(adt.dist)/sum(rna.dist)
    joint.dist <- adt.dist*alpha + rna.dist*(1-alpha)*scale.factor
    object@misc[['jointDist']] <- joint.dist
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
  alpha = NULL
) {
  if(is.null(alpha)) {
    cat("will automatically determine an alpha to calculate the joint distance \n")
    learning.rate <- 1
    low.threshold <- 1e-5
    alpha <- 0
    X <- dist1 + dist2
    Y <- dist2
    optimal <- gradientDescent(X, Y, alpha, learning.rate, low.threshold)

    joint.dist <- dist1*optimal + dist2*(1-optimal)
  } else if (is.numeric(alpha)) {
    cat("will use user's alpha value",alpha," to calculate the joint distance \n")
    # calculate scale factor to scale two distances into the same level
    scale.factor <- sum(dist1)/sum(dist2)
    joint.dist <- dist1*alpha + dist2*(1-alpha)*scale.factor
  } else {
    stop("Please provide an alpha between 0 and 1, or you can leave alpha = NULL, we will determine for you \n")
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


