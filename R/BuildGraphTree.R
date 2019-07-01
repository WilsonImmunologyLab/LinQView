#' buildMST.Seurat
#'
#' @rdname buildMST
#' @export
#'
buildMST.Seurat <- function(
  object,
  assay = "Joint"
) {
  if(!is.null(object)) {
    if(assay == "All") {
      rna.dist <- object@misc[['rnaDist']]
      adt.dist <- object@misc[['adtDist']]
      joint.dist <- object@misc[['jointDist']]

      # build MST based on cell-cell distances
      joint.tree <- mlpack_mst(as.matrix(joint.dist))
      rna.tree <- mlpack_mst(as.matrix(rna.dist))
      adt.tree <- mlpack_mst(as.matrix(adt.dist))

      object@misc[["rnaTree"]] <- rna.tree
      object@misc[["adtTree"]] <- adt.tree
      object@misc[["jointTree"]] <- joint.tree
    } else if (assay == "Joint") {
      joint.dist <- object@misc[['jointDist']]

      # build MST based on cell-cell distances
      joint.tree <- mlpack_mst(as.matrix(joint.dist))
      object@misc[["jointTree"]] <- joint.tree
    } else if (assay == "RNA") {
      rna.dist <- object@misc[['rnaDist']]

      # build MST based on cell-cell distances
      rna.tree <- mlpack_mst(as.matrix(rna.dist))
      object@misc[["rnaTree"]] <- rna.tree
    } else if (assay == "ADT") {
      adt.dist <- object@misc[['adtDist']]

      # build MST based on cell-cell distances
      adt.tree <- mlpack_mst(as.matrix(adt.dist))
      object@misc[["adtTree"]] <- adt.tree
    } else {
      stop("Invalid assay name! Assay name can only be RNA, ADT, Joint or All")
    }
    return(object)
  } else {
    stop("Please provide an Seurat object!")
  }
}

#' buildMST.dist
#'
#' @rdname buildMST
#' @export
#'
buildMST.dist <- function(
  object
) {
  if(!is.null(object)) {
    # build MST based on cell-cell distances
    mst <- mlpack_mst(as.matrix(object))
    return(mst)
  } else {
    stop("Please provide an Seurat object!")
  }
}


#' buildMST.matrix
#'
#' @rdname buildMST
#' @export
#'
buildMST.matrix <- function(
  object
) {
  if(!is.null(object)) {
    # build MST based on cell-cell distances
    mst <- mlpack_mst(object)
    return(mst)
  } else {
    stop("Please provide an Seurat object!")
  }
}

#' buildMST.default
#'
#' @rdname buildMST
#' @export
#'
buildMST.default <- function(
  object
) {
  if(!is.null(object)) {
    # build MST based on cell-cell distances
    mst <- mlpack_mst(as.matrix(object))
    return(mst)
  } else {
    stop("Please provide an Seurat object!")
  }
}



#' buildPG.Seurat
#'
#' @rdname buildPG
#' @export
#'
buildPG.Seurat <- function(
  object,
  assay = "Joint",
  reduction.prefix = "tsne_",
  pg.nodes = NULL,
  pg.min.nodes = 30,
  pg.Lambda = 0.03,
  pg.Mu = 0.01,
  trimming.radius = Inf,
  final.energy = "Penalized",
  initial.MST = FALSE
) {
  if(!is.null(object)) {
    if(assay == "All") {
      # build PG for joint data
      cat("build PG for joint data ... \n")
      joint.mst <- object@misc[['jointTree']]
      joint.reduction.name <- paste0(reduction.prefix,"joint")
      joint.embedding <- object@reductions[[joint.reduction.name]]@cell.embeddings

      object@misc[['jointPG']] <- buildPG(object = joint.embedding, mst = joint.mst, pg.nodes = pg.nodes, pg.min.nodes = pg.min.nodes, pg.Lambda = pg.Lambda, pg.Mu = pg.Mu, initial.MST = initial.MST, trimming.radius = trimming.radius, final.energy = final.energy)

      # build PG for RNA data
      cat("build PG for RNA data ... \n")
      rna.mst <- object@misc[['rnaTree']]
      rna.reduction.name <- paste0(reduction.prefix,"rna")
      rna.embedding <- object@reductions[[rna.reduction.name]]@cell.embeddings

      object@misc[['rnaPG']] <- buildPG(object = rna.embedding, mst = rna.mst, pg.nodes = pg.nodes, pg.min.nodes = pg.min.nodes, pg.Lambda = pg.Lambda, pg.Mu = pg.Mu, initial.MST = initial.MST, trimming.radius = trimming.radius, final.energy = final.energy)

      # build PG for ADT data
      cat("build PG for ADT data ... \n")
      adt.mst <- object@misc[['adtTree']]
      adt.reduction.name <- paste0(reduction.prefix,"adt")
      adt.embedding <- object@reductions[[adt.reduction.name]]@cell.embeddings

      object@misc[['adtPG']] <- buildPG(object = adt.embedding, mst = adt.mst, pg.nodes = pg.nodes, pg.min.nodes = pg.min.nodes, pg.Lambda = pg.Lambda, pg.Mu = pg.Mu, initial.MST = initial.MST, trimming.radius = trimming.radius, final.energy = final.energy)
    } else if (assay == "Joint") {
      # build PG for joint data
      cat("build PG for joint data ... \n")
      joint.mst <- object@misc[['jointTree']]
      joint.reduction.name <- paste0(reduction.prefix,"joint")
      joint.embedding <- object@reductions[[joint.reduction.name]]@cell.embeddings

      object@misc[['jointPG']] <- buildPG(object = joint.embedding, mst = joint.mst, pg.nodes = pg.nodes, pg.min.nodes = pg.min.nodes, pg.Lambda = pg.Lambda, pg.Mu = pg.Mu, initial.MST = initial.MST, trimming.radius = trimming.radius, final.energy = final.energy)
    } else if (assay == "RNA") {
      # build PG for RNA data
      cat("build PG for RNA data ... \n")
      rna.mst <- object@misc[['rnaTree']]
      rna.reduction.name <- paste0(reduction.prefix,"rna")
      rna.embedding <- object@reductions[[rna.reduction.name]]@cell.embeddings

      object@misc[['rnaPG']] <- buildPG(object = rna.embedding, mst = rna.mst, pg.nodes = pg.nodes, pg.min.nodes = pg.min.nodes, pg.Lambda = pg.Lambda, pg.Mu = pg.Mu, initial.MST = initial.MST, trimming.radius = trimming.radius, final.energy = final.energy)
    } else if (assay == "ADT") {
      # build PG for ADT data
      cat("build PG for ADT data ... \n")
      adt.mst <- object@misc[['adtTree']]
      adt.reduction.name <- paste0(reduction.prefix,"adt")
      adt.embedding <- object@reductions[[adt.reduction.name]]@cell.embeddings

      object@misc[['adtPG']] <- buildPG(object = adt.embedding, mst = adt.mst, pg.nodes = pg.nodes, pg.min.nodes = pg.min.nodes, pg.Lambda = pg.Lambda, pg.Mu = pg.Mu, initial.MST = initial.MST, trimming.radius = trimming.radius, final.energy = final.energy)
    } else {
      stop("Invalid assay name! Assay name can only be RNA, ADT, Joint or All")
    }
    return(object)
  } else {
    stop("Please provide an Seurat object!")
  }
}




#' buildPG.default
#'
#' @rdname buildPG
#' @export
#'
buildPG.default <- function(
  object,
  mst = NULL,
  pg.nodes = NULL,
  pg.min.nodes = 30,
  pg.Lambda = 0.03,
  pg.Mu = 0.01,
  trimming.radius = Inf,
  final.energy = "Penalized",
  initial.MST = FALSE
) {
  if(!is.null(object)) {
    embedding <- object
    if(!is.null(pg.nodes)) {
      nodes <- pg.nodes
    } else {
      n.point <- dim(mst)[2]
      nodes <- floor(n.point/100)
      nodes <- ifelse(test = nodes > pg.min.nodes, yes = nodes, no = pg.min.nodes)
    }

    mst <- t(mst)[,1:2]
    if(isTRUE(initial.MST)) {
      if(is.null(mst)) {
        stop("You set 'initial.MST = TRUE', please provide an MST using 'mst =' paramater \n")
      }
      # infer PG from cell embedding and use trimmed MST as initial
      trimed <- trimMST(mst = mst, embedding = embedding)
      PrincipalTree <- computeElasticPrincipalTree(X = embedding,
                                                   NumNodes = nodes,
                                                   Lambda = pg.Lambda,
                                                   Mu = pg.Mu,
                                                   InitNodePositions = trimed[[2]],
                                                   InitEdges = trimed[[1]],
                                                   FinalEnergy = final.energy,
                                                   TrimmingRadius= trimming.radius,
                                                   Do_PCA = FALSE,
                                                   verbose = FALSE,
                                                   CenterData=FALSE,
                                                   drawAccuracyComplexity = FALSE,
                                                   drawEnergy = FALSE,
                                                   drawPCAView = FALSE)
    } else {
      # infer PG from cell embedding without initial
      PrincipalTree <- computeElasticPrincipalTree(X = embedding,
                                                   NumNodes = nodes,
                                                   Lambda = pg.Lambda,
                                                   Mu = pg.Mu,
                                                   InitNodePositions = NULL,
                                                   InitEdges = NULL,
                                                   FinalEnergy = final.energy,
                                                   TrimmingRadius= trimming.radius,
                                                   Do_PCA = FALSE,
                                                   verbose = FALSE,
                                                   CenterData=FALSE,
                                                   drawAccuracyComplexity = FALSE,
                                                   drawEnergy = FALSE,
                                                   drawPCAView = FALSE)
    }
    return(PrincipalTree)
  } else {
    stop("Please provide an Embedding!")
  }
}




#' buildKNN.Seurat
#'
#' @rdname buildKNN
#' @export
#'
buildKNN.Seurat <- function(
  object,
  assay = "Joint",
  k = 10
) {
  if(!is.null(object)) {
    if(assay == "All") {
      rna.dist <- object@misc[['rnaDist']]
      adt.dist <- object@misc[['adtDist']]
      joint.dist <- object@misc[['jointDist']]

      # build KNN graph based on cell-cell distances
      cat("building kNN graph from joint distance... k = ",k,"\n")
      joint.graph <- fastKNN(joint.dist, k = k)
      cat("building kNN graph from RNA distance... k = ",k,"\n")
      rna.graph <- fastKNN(rna.dist, k = k)
      cat("building kNN graph from ADT distance... k = ",k,"\n")
      adt.graph <- fastKNN(adt.dist, k = k)

      object@misc[["rnaGraph"]] <- rna.graph
      object@misc[["adtGraph"]] <- adt.graph
      object@misc[["jointGraph"]] <- joint.graph
    } else if (assay == "Joint") {
      joint.dist <- object@misc[['jointDist']]

      # build KNN graph based on cell-cell distances
      cat("building kNN graph from joint distance... k = ",k,"\n")
      joint.graph <- fastKNN(joint.dist, k = k)
      object@misc[["jointGraph"]] <- joint.graph
    } else if (assay == "RNA") {
      rna.dist <- object@misc[['rnaDist']]

      # build KNN graph based on cell-cell distances
      cat("building kNN graph from RNA distance... k = ",k,"\n")
      rna.graph <- fastKNN(rna.dist, k = k)
      object@misc[["rnaGraph"]] <- rna.graph
    } else if (assay == "ADT") {
      adt.dist <- object@misc[['adtDist']]

      # build KNN graph based on cell-cell distances
      cat("building kNN graph from ADT distance... k = ",k,"\n")
      adt.graph <- fastKNN(adt.dist, k = k)
      object@misc[["adtGraph"]] <- adt.graph
    } else {
      stop("Invalid assay name! Assay name can only be RNA, ADT, Joint or All")
    }
    return(object)
  } else {
    stop("Please provide an Seurat object!")
  }
}


#' buildKNN.default
#'
#' @rdname buildKNN
#' @export
#'
buildKNN.default <- function(
  object,
  assay,
  k = 10
) {
  if(!is.null(object)) {
    # build KNN graph based on cell-cell distances
    cat("building kNN graph from distance... k = ",k,"\n")
    graph <- fastKNN(object, k = k)

    return(graph)
  } else {
    stop("Please provide an Seurat object!")
  }
}



#' fastKNN
#'
#' run fastKNN method
#'
#' @param dist dist object or matrix
#' @param k k value for KNN graph
#'
#' @export
fastKNN <- function(
  dist = NULL,
  k = NULL
) {
  dist <- as.matrix(dist)
  n = dim(dist)[1]

  knn = matrix(0,n,k) # n x k
  for (i in 1:n) {
    knn[i,] = k.nearest.neighbors(i, dist, k = k)
  }
  return(knn)
}

#' batchPseudoTime
#'
#' calculate pseudotime using a given start point from MST or KNN. This function only can be used if you know which cell is the root. Not suggested.
#'
#' @param object Seurat object
#' @param method could be MST or KNN
#' @param mst user choose the assay data MST they want to use for pseudotime, can be RNA, ADT, Joint or All (All means calculate all three), use also can determine any mst (for example, mst = "rnaTree", the rnaTree must be located in object[at]misc[['rnaTree']]
#' @param knn user choose the assay data KNN they want to use for pseudotime, can be RNA, ADT, Joint or All (All means calculate all three), use also can determine any knn by name (for example, knn = "rnaGraph", the rnaTree must be located in object[at]misc[['rnaGraph']]
#' @param dist if use choose their own KNN or MST, user should determine the dist object or distance matrix for the knn graph or MST. the dist can be dist object, or the name of the dist object that is located in object[at]misc[[dist]]. for example, user can set dist = object[at]misc[["userDist"]] or knn.dist = "userDist"
#' @param root user choose the root cell
#' @export
batchPseudoTime <- function(
  object = NULL,
  method = c("MST", "KNN"),
  mst = "Joint",
  knn = "Joint",
  dist = NULL,
  root = NULL
) {
  if(!is.null(object)) {
    if(method == "MST") {
      if(mst == "All") {
        # MST
        rna.tree <- object@misc[['rnaTree']]
        adt.tree <- object@misc[['adtTree']]
        joint.tree <- object@misc[['jointTree']]

        # order cells
        cat("calculating pseudotime for RNA data...\n")
        rna.time <- orderCellsMST(data = rna.tree,root = root)
        cat("calculating pseudotime for ADT data...\n")
        adt.time <- orderCellsMST(data = adt.tree,root = root)
        cat("calculating pseudotime for joint data...\n")
        joint.time <- orderCellsMST(data = joint.tree,root = root)

        # save
        object$rnaMSTTime <- rna.time
        object$adtMSTTime <- adt.time
        object$jointMSTTime <- joint.time

      } else if (mst == "Joint") {
        # MST
        joint.tree <- object@misc[['jointTree']]

        # order cells
        cat("calculating pseudotime for joint data...\n")
        joint.time <- orderCellsMST(data = joint.tree,root = root)

        # save
        object$jointMSTTime <- joint.time
      } else if (mst == "RNA") {
        # MST
        rna.tree <- object@misc[['rnaTree']]

        # order cells
        cat("calculating pseudotime for RNA data...\n")
        rna.time <- orderCellsMST(data = rna.tree,root = root)

        # save
        object$rnaMSTTime <- rna.time
      } else if (mst == "ADT") {
        # MST
        adt.tree <- object@misc[['adtTree']]

        # order cells
        cat("calculating pseudotime for ADT data...\n")
        adt.time <- orderCellsMST(data = adt.tree,root = root)

        # save
        object$adtMSTTime <- adt.time
      } else {
        # MST
        cur.tree <- object@misc[[mst]]

        # order cells
        cat("calculating pseudotime for your data...\n")
        cur.time <- orderCellsMST(data = adt.tree,root = root)
        timename <- paste0(mst,"Time")

        # save
        object@meta.data[[timename]] <- cur.time
      }
      return(object)
    } else if (method == "KNN") {
      if(knn == "All") {
        # KNN
        rna.graph <- object@misc[['rnaGraph']]
        adt.graph <- object@misc[['adtGraph']]
        joint.graph <- object@misc[['jointGraph']]

        rna.dist <- object@misc[['rnaDist']]
        adt.dist <- object@misc[['adtDist']]
        joint.dist <- object@misc[['jointDist']]

        # order cells
        cat("calculating pseudotime for RNA data...\n")
        rna.time <- orderCellsKNN(data = rna.graph,dist = rna.dist,root = root)
        cat("calculating pseudotime for ADT data...\n")
        adt.time <- orderCellsKNN(data = adt.graph,dist = adt.dist,root = root)
        cat("calculating pseudotime for Joint data...\n")
        joint.time <- orderCellsKNN(data = joint.graph,dist = joint.dist,root = root)

        # save
        object$rnaKNNTime <- rna.time
        object$adtKNNTime <- adt.time
        object$jointKNNTime <- joint.time

      } else if (knn == "Joint") {
        # KNN
        joint.graph <- object@misc[['jointGraph']]
        joint.dist <- object@misc[['jointDist']]

        # order cells
        cat("calculating pseudotime for Joint data...\n")
        joint.time <- orderCellsKNN(data = joint.graph,dist = joint.dist,root = root)

        # save
        object$jointKNNTime <- joint.time
      } else if (knn == "RNA") {
        # KNN
        rna.graph <- object@misc[['rnaGraph']]
        rna.dist <- object@misc[['rnaDist']]

        # order cells
        cat("calculating pseudotime for RNA data...\n")
        rna.time <- orderCellsKNN(data = rna.graph,dist = rna.dist,root = root)

        # save
        object$rnaTime <- rna.time
      } else if (knn == "ADT") {
        # KNN
        adt.graph <- object@misc[['adtGraph']]
        adt.dist <- object@misc[['adtDist']]

        # order cells
        cat("calculating pseudotime for ADT data...\n")
        adt.time <- orderCellsKNN(data = adt.graph,dist = adt.dist,root = root)

        # save
        object$adtKNNTime <- adt.time
      } else {
        # KNN
        cur.graph <- object@misc[[knn]]
        if(is.character(dist)) {
          cur.dist <- object@misc[[dist]]
        } else {
          cur.dist <- dist
        }

        # order cells
        cat("calculating pseudotime for your data...\n")
        cur.time <- orderCellsKNN(data = cur.graph,dist = cur.dist,root = root)
        timename <- paste0(graph,"Time")

        # save
        object@meta.data[[timename]] <- cur.time
      }
      return(object)
    } else {
      stop("Please provide a valid method, MST or KNN!")
    }
  } else {
    stop("Please provide an Seurat object!")
  }
}


#' pseudoTime
#'
#' calculate pseudotime using a given start point from MST or KNN.
#'
#' @param object Seurat object
#' @param method could be MST or KNN
#' @param group.by cell clustering results.
#' @param graph kNN or MST graph. use can determine any knn by name.for example, knn = "rnaGraph", the rnaTree must be located in object[at]misc[['rnaGraph']]
#' @param dist if use choose their own KNN or MST, user should determine the dist object or distance matrix for the knn graph or MST. the dist can be dist object, or the name of the dist object that is located in object[at]misc[[dist]]. for example, user can set dist = object[at]misc[["userDist"]] or knn.dist = "userDist"
#' @param root.cluster user choose the root cell cluster
#'
#' @export
#'
pseudoTime <- function(
  object,
  method = c("MST", "KNN"),
  group.by = NULL,
  graph = NULL,
  dist = NULL,
  root.cluster = NULL
) {
  if(!is.null(object)) {
    if(is.character(dist)) {
      cur.dist <- object@misc[[dist]]
    } else {
      cur.dist <- dist
    }

    if(root.cluster %in% levels(object@meta.data[[group.by]])) {
      index1 <- which(object@meta.data[[group.by]] == root.cluster)
      index2 <- which(object@meta.data[[group.by]] != root.cluster)
      root.id <- findRootNodeID(dist = cur.dist, root.cluster.index = index1, other.cluster.index = index2)
    } else {
      stop("Your root.cluster name ",root.cluster," can not be found in the group.by ", group.by, "\n")
    }

    if(method == "MST") {
      # MST
      tree <- object@misc[[graph]]

      # order cells
      cat("calculating pseudotime for given data...",root.id,"\n")
      time <- orderCellsMST(data = tree, root = root.id)

      time.name <- paste0(graph,"MSTTime")
      # save
      object@meta.data[[time.name]] <- time

    } else if (method == "KNN") {
      # KNN
      knn <- object@misc[[graph]]

      # order cells
      cat("calculating pseudotime for ADT data...\n")
      time <- orderCellsKNN(data = knn,dist = cur.dist,root = root.id)

      time.name <- paste0(graph,"KNNTime")
      # save
      object@meta.data[[time.name]] <- time
    } else {
      stop("Please provide a valid method, MST or KNN!")
    }
    return(object)
  } else {
    stop("Please provide an Seurat object!")
  }
}

#' findRootNodeID
#'
#' calculate pseduotime using a given start point from MST
#'
#' @param dist a MST
#' @param root.cluster.index user choose the a root node
#'
#' @export
#'
findRootNodeID<- function(
  dist = NULL,
  root.cluster.index = NULL,
  other.cluster.index = NULL,
  option = 1
) {
  # we calculate distance from each node of candicate cluster to all nodes in other clusters. The node with largest overall distance will be the root node
  dist <- as.matrix(dist)
  if(option == 1) {
    dist.within.cluster <- dist[root.cluster.index, root.cluster.index]
    dist.within.cluster <- rowSums(dist.within.cluster)
    root.id <- which(dist.within.cluster == min(dist.within.cluster))
    root.id <- root.cluster.index[root.id]
  } else {
    dist.to.other <- dist[root.cluster.index, other.cluster.index]
    dist.to.other <- rowSums(dist.to.other)
    root.id <- which(dist.to.other == max(dist.to.other))
    root.id <- root.cluster.index[root.id]
  }
  return(root.id)
}


#' orderCellsMST
#'
#' calculate pseduotime using a given start point from MST
#'
#' @param data a MST
#' @param root user choose the a root node
#'
#' @export

orderCellsMST <- function(
  data = NULL,
  root = NULL
) {
  if(!is.null(data)) {
    origin1 <- data
    origin2 <- origin1
    origin2[1,] <- origin1[2,]
    origin2[2,] <- origin1[1,]
    origin <- t(cbind(origin1,origin2))
    origin <- as.data.frame(origin)
    colnames(origin) <- c("from","to","distance")

    # temp pool waiting for process
    pool <- data.frame(node = c(NA), distance = c(NA))
    pool <- pool[-1,]

    # start
    StartNode <- root
    pastDist <- 0

    # result data frame
    result <- data.frame(node = c(StartNode), distance = c(pastDist))

    # ugly code start ... try to find the distance from root node to all other nodes on the MST
    tmp <- origin[origin$from == StartNode,]
    for (i in 1:dim(tmp)[1]) {
      dis <- tmp$distance[i] + pastDist
      row <- data.frame(node=tmp$to[i],distance=dis)

      pool <- rbind(pool, row)
      result <- rbind(result, row)
    }

    a <- (origin$from == StartNode)
    b <- (origin$to == StartNode)
    row_a <- as.numeric(rownames(origin[a,]))
    row_b <- as.numeric(rownames(origin[b,]))
    row <- c(row_a,row_b)
    if(length(row) > 0){
      origin[row,] <- c(0,0,0)
    }

    while (dim(pool)[1] > 0) {
      StartNode <- pool$node[1]
      pastDist <- pool$distance[1]
      pool <- pool[-1,]
      if(dim(pool)[1] > 0){
        rownames(pool) <- 1:dim(pool)[1]
      }

      tmp <- origin[origin$from == StartNode,]
      if(dim(tmp)[1] > 0){
        for (i in 1:dim(tmp)[1]) {
          dis <- tmp$distance[i] + pastDist
          row <- data.frame(node=tmp$to[i],distance=dis)

          pool <- rbind(pool, row)
          result <- rbind(result, row)
        }
      }
      a <- (origin$from == StartNode)
      b <- (origin$to == StartNode)
      row_a <- as.numeric(rownames(origin[a,]))
      row_b <- as.numeric(rownames(origin[b,]))
      row <- c(row_a,row_b)
      if(length(row) > 0){
        origin[row,] <- c(0,0,0)
      }
    }
    # ugly code end
    result <- result[order(result$node),]
    pseduotime <- result$distance/max(result$distance)

    return(pseduotime)
  } else {
    stop("Please provide MST!")
  }
}


#' orderCellsKNN
#'
#' calculate pseduotime using a given start point from KNN
#'
#' @param data a KNN
#' @param dist dist matrix the KNN calculated from
#' @param root user choose the a root node
#'
#' @export

orderCellsKNN <- function(
  data = NULL,
  dist = NULL,
  root = NULL
) {
  if(!is.null(data)) {
    # generate a matrix from KNN
    n <- dim(data)[1]
    dist <- as.matrix(dist)
    x <- matrix(0, n, n)
    diag(x) <- 0
    for (i in 1:n) {
      x[i,data[i,]] <- dist[i,data[i,]]
    }

    # Create graphs from adjacency matrices using igraph
    my.graph <- graph_from_adjacency_matrix(x, mode = "undirected",weighted = TRUE, diag = FALSE)

    # Calculate shortest path
    #short.path <- shortest_paths(my.graph, from = root, to = V(my.graph))

    # Calculate shortest distances
    short.dist <- shortest.paths(my.graph, v=root, to=V(my.graph))

    time <- as.numeric(t(short.dist))

    # scale to [0,1]
    pseduotime <- time/max(time)

    return(pseduotime)
  } else {
    stop("Please provide KNN!")
  }
}



#' trimMST
#'
#' trim MST, remove leaf. This step create initial tree for principle graph method.
#'
#' @param mst a MST
#' @param embedding cell embedding
#'
#' @export

trimMST <- function(
  mst = NULL,
  embedding = NULL
) {
  mst.all.nodes <- c(mst[,1],mst[,2])
  mst.count <- table(mst.all.nodes)
  mst.leaf.node <- which(mst.count == 1)

  index.1 <- which(mst[,1] %in% mst.leaf.node)
  index.2 <- which(mst[,2] %in% mst.leaf.node)
  index <- unique(c(index.1,index.2))
  mst.trimmed <- mst[-index,]
  embedding.trimmed <- embedding[-mst.leaf.node,]

  # make a projection between old node numbering and new node numbering
  ori.names <- rownames(embedding)
  new.names <- rownames(embedding.trimmed)

  location <- ori.names %in% new.names
  projection <- rep(NA,length(location))
  new.index <- 0
  for (i in 1:length(location)) {
    if(isTRUE(location[i])) {
      new.index <- new.index + 1
      projection[i] <- new.index
    }
  }

  mst.trimmed[,1] <- projection[mst.trimmed[,1]]
  mst.trimmed[,2] <- projection[mst.trimmed[,2]]

  result <- list()
  result[[1]] <- mst.trimmed
  result[[2]] <- embedding.trimmed

  return(result)
}
