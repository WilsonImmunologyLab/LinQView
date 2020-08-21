#' buildMST.Seurat
#' @importFrom emstreeR mlpack_mst
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
      rna.dist <- object@misc[['RNA']][['dist']]
      adt.dist <- object@misc[['ADT']][['dist']]
      joint.dist <- object@misc[['Joint']][['dist']]

      # build MST based on cell-cell distances
      joint.tree <- mlpack_mst(as.matrix(joint.dist))
      rna.tree <- mlpack_mst(as.matrix(rna.dist))
      adt.tree <- mlpack_mst(as.matrix(adt.dist))

      object@misc[["RNA"]][['mst']] <- rna.tree
      object@misc[["ADT"]][['mst']] <- adt.tree
      object@misc[["Joint"]][['mst']] <- joint.tree
    } else if (assay == "Joint") {
      joint.dist <- object@misc[['Joint']][['dist']]

      # build MST based on cell-cell distances
      joint.tree <- mlpack_mst(as.matrix(joint.dist))
      object@misc[["Joint"]][['mst']] <- joint.tree
    } else if (assay == "RNA") {
      rna.dist <- object@misc[['RNA']][['dist']]

      # build MST based on cell-cell distances
      rna.tree <- mlpack_mst(as.matrix(rna.dist))
      object@misc[["RNA"]][['mst']] <- rna.tree
    } else if (assay == "ADT") {
      adt.dist <- object@misc[['ADT']][['dist']]

      # build MST based on cell-cell distances
      adt.tree <- mlpack_mst(as.matrix(adt.dist))
      object@misc[["ADT"]][['mst']] <- adt.tree
    } else {
      stop("Invalid assay name! Assay name can only be RNA, ADT, Joint or All")
    }
    return(object)
  } else {
    stop("Please provide an Seurat object!")
  }
}

#' buildMST.dist
#' @importFrom emstreeR mlpack_mst
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
#' @importFrom emstreeR mlpack_mst
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
#' @importFrom emstreeR mlpack_mst
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
#' @importFrom ElPiGraph.R computeElasticPrincipalTree
#'
#' @rdname buildPG
#' @export
#'
buildPG.Seurat <- function(
  object,
  assay = "Joint",
  reduction = NULL,
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
      if(!is.null(reduction)) {
        stop("reduction name can only be used for single Assay! Please set reduction to NULL is you want run ALL Assay. \n")
      }
      # build PG for joint data
      cat("build PG for joint data ... \n")
      joint.mst <- object@misc[['Joint']][['mst']]
      joint.reduction.name <- object@misc[['Joint']][['reduction']]
      joint.embedding <- object@reductions[[joint.reduction.name]]@cell.embeddings

      object@misc[['Joint']][['pg']] <- buildPG(object = joint.embedding, mst = joint.mst, pg.nodes = pg.nodes, pg.min.nodes = pg.min.nodes, pg.Lambda = pg.Lambda, pg.Mu = pg.Mu, initial.MST = initial.MST, trimming.radius = trimming.radius, final.energy = final.energy)

      # build PG for RNA data
      cat("build PG for RNA data ... \n")
      rna.mst <- object@misc[['RNA']][['mst']]
      rna.reduction.name <- object@misc[['RNA']][['reduction']]
      rna.embedding <- object@reductions[[rna.reduction.name]]@cell.embeddings

      object@misc[['RNA']][['pg']] <- buildPG(object = rna.embedding, mst = rna.mst, pg.nodes = pg.nodes, pg.min.nodes = pg.min.nodes, pg.Lambda = pg.Lambda, pg.Mu = pg.Mu, initial.MST = initial.MST, trimming.radius = trimming.radius, final.energy = final.energy)

      # build PG for ADT data
      cat("build PG for ADT data ... \n")
      adt.mst <- object@misc[['ADT']][['mst']]
      adt.reduction.name <- object@misc[['ADT']][['reduction']]
      adt.embedding <- object@reductions[[adt.reduction.name]]@cell.embeddings

      object@misc[['ADT']][['pg']] <- buildPG(object = adt.embedding, mst = adt.mst, pg.nodes = pg.nodes, pg.min.nodes = pg.min.nodes, pg.Lambda = pg.Lambda, pg.Mu = pg.Mu, initial.MST = initial.MST, trimming.radius = trimming.radius, final.energy = final.energy)
    } else if (assay == "Joint") {
      # build PG for joint data
      cat("build PG for joint data ... \n")
      joint.mst <- object@misc[['Joint']][['mst']]
      if(!is.null(reduction)) {
        joint.reduction.name <- reduction
      } else {
        joint.reduction.name <- object@misc[['Joint']][['reduction']]
      }
      joint.embedding <- object@reductions[[joint.reduction.name]]@cell.embeddings

      object@misc[['Joint']][['pg']] <- buildPG(object = joint.embedding, mst = joint.mst, pg.nodes = pg.nodes, pg.min.nodes = pg.min.nodes, pg.Lambda = pg.Lambda, pg.Mu = pg.Mu, initial.MST = initial.MST, trimming.radius = trimming.radius, final.energy = final.energy)
    } else if (assay == "RNA") {
      # build PG for RNA data
      cat("build PG for RNA data ... \n")
      rna.mst <- object@misc[['RNA']][['mst']]
      if(!is.null(reduction)) {
        rna.reduction.name <- reduction
      } else {
        rna.reduction.name <- object@misc[['RNA']][['reduction']]
      }
      rna.embedding <- object@reductions[[rna.reduction.name]]@cell.embeddings

      object@misc[['RNA']][['pg']] <- buildPG(object = rna.embedding, mst = rna.mst, pg.nodes = pg.nodes, pg.min.nodes = pg.min.nodes, pg.Lambda = pg.Lambda, pg.Mu = pg.Mu, initial.MST = initial.MST, trimming.radius = trimming.radius, final.energy = final.energy)
    } else if (assay == "ADT") {
      # build PG for ADT data
      cat("build PG for ADT data ... \n")
      adt.mst <- object@misc[['ADT']][['mst']]
      if(!is.null(reduction)) {
        adt.reduction.name <- reduction
      } else {
        adt.reduction.name <- object@misc[['ADT']][['reduction']]
      }
      adt.embedding <- object@reductions[[adt.reduction.name]]@cell.embeddings

      object@misc[['ADT']][['pg']] <- buildPG(object = adt.embedding, mst = adt.mst, pg.nodes = pg.nodes, pg.min.nodes = pg.min.nodes, pg.Lambda = pg.Lambda, pg.Mu = pg.Mu, initial.MST = initial.MST, trimming.radius = trimming.radius, final.energy = final.energy)
    } else {
      stop("Invalid assay name! Assay name can only be RNA, ADT, Joint or All")
    }
    return(object)
  } else {
    stop("Please provide an Seurat object!")
  }
}




#' buildPG.default
#' @importFrom ElPiGraph.R computeElasticPrincipalTree
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
      rna.dist <- object@misc[['RNA']][['dist']]
      adt.dist <- object@misc[['ADT']][['dist']]
      joint.dist <- object@misc[['Joint']][['dist']]

      # build KNN graph based on cell-cell distances
      cat("building kNN graph from joint distance... k = ",k,"\n")
      joint.graph <- fastKNN(joint.dist, k = k)
      cat("building kNN graph from RNA distance... k = ",k,"\n")
      rna.graph <- fastKNN(rna.dist, k = k)
      cat("building kNN graph from ADT distance... k = ",k,"\n")
      adt.graph <- fastKNN(adt.dist, k = k)

      object@misc[["RNA"]][['knn']] <- rna.graph
      object@misc[["ADT"]][['knn']] <- adt.graph
      object@misc[["Joint"]][['knn']] <- joint.graph
    } else if (assay == "Joint") {
      joint.dist <- object@misc[['Joint']][['dist']]

      # build KNN graph based on cell-cell distances
      cat("building kNN graph from joint distance... k = ",k,"\n")
      joint.graph <- fastKNN(joint.dist, k = k)
      object@misc[["Joint"]][['knn']] <- joint.graph
    } else if (assay == "RNA") {
      rna.dist <- object@misc[['RNA']][['dist']]

      # build KNN graph based on cell-cell distances
      cat("building kNN graph from RNA distance... k = ",k,"\n")
      rna.graph <- fastKNN(rna.dist, k = k)
      object@misc[["RNA"]][['knn']] <- rna.graph
    } else if (assay == "ADT") {
      adt.dist <- object@misc[['ADT']][['dist']]

      # build KNN graph based on cell-cell distances
      cat("building kNN graph from ADT distance... k = ",k,"\n")
      adt.graph <- fastKNN(adt.dist, k = k)
      object@misc[["ADT"]][['knn']] <- adt.graph
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
#' @importFrom FastKNN k.nearest.neighbors
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


#' pseudoTime
#'
#' calculate pseudotime using a given start point from MST or KNN.
#'
#' @param object Seurat object
#' @param method could be MST
#' @param assay Assay name. Can be RNA, ADT or Joint
#' @param root.cluster user choose the root cell cluster
#' @param percentile.cutoff setup a percentile cutoff for scaling distance to pseudoTime, can avoid some outliers
#'
#' @export
#'
pseudoTime <- function(
  object,
  method = "MST",
  assay = "Joint",
  root.cluster = NULL,
  percentile.cutoff = 1
) {
  if(!is.null(object)) {
    group.by <- object@misc[[assay]][['cluster']]
    cur.dist <- object@misc[[assay]][['dist']]
    # identify root cell
    if(root.cluster %in% levels(object@meta.data[[group.by]])) {
      index1 <- which(object@meta.data[[group.by]] == root.cluster)
      index2 <- which(object@meta.data[[group.by]] != root.cluster)
      root.id <- findRootNodeID(dist = cur.dist, root.cluster.index = index1, other.cluster.index = index2)
    } else {
      stop("Your root.cluster name ",root.cluster," can not be found in the group.by ", group.by, "\n")
    }

    if(method == "MST") {
      # MST
      tree <- object@misc[[assay]][['mst']]

      # order cells
      cat("calculating pseudotime for ",assay," data...\n")
      time <- orderCellsMST(data = tree, root = root.id, cutoff = percentile.cutoff)

      time.name <- paste0(tolower(assay),"MSTtime")
      # save
      object@meta.data[[time.name]] <- time
      object@misc[[assay]][['time']][['mst']] <- time.name
    } else if (method == "KNN") {
      stop("your method can not be regconized! Only can be MST!")
      # KNN
      knn <- object@misc[[assay]][['knn']]

      # order cells
      cat("calculating pseudotime for ",assay," data...\n")
      time <- orderCellsKNN(data = knn,dist = cur.dist,root = root.id, cutoff = percentile.cutoff)

      time.name <- paste0(tolower(assay),"KNNtime")
      # save
      object@meta.data[[time.name]] <- time
      object@misc[[assay]][['time']][['knn']] <- time.name
    } else {
      stop("your method can not be regconized! Only can be MST!")
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
#' @param cutoff setup a percentile cutoff for scaling distance to pseudoTime, can avoid some outliers
#'
#' @export

orderCellsMST <- function(
  data = NULL,
  root = NULL,
  cutoff = 1
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

    sort.dist <- sort(result$distance)
    p95.dist <- sort.dist[floor(length(sort.dist)*cutoff)]

    pseduotime <- result$distance/p95.dist
    pseduotime[which(pseduotime > 1)] <- 1

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
#' @param cutoff setup a percentile cutoff for scaling distance to pseudoTime, can avoid some outliers
#'
#' @importFrom igraph graph_from_adjacency_matrix shortest.paths V
#'
#' @export

orderCellsKNN <- function(
  data = NULL,
  dist = NULL,
  root = NULL,
  cutoff = 1
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
    sort.dist <- sort(time)
    p95.dist <- sort.dist[floor(length(sort.dist)*cutoff)]

    pseduotime <- time/p95.dist
    pseduotime[which(pseduotime > 1)] = 1

    return(pseduotime)
  } else {
    stop("Please provide KNN!")
  }
}


#' trimMST
#'
#' trim MST, remove all 1) leaf nodes (nodes only connected to one node), 2) connecting nodes (nodes connected to two nodes after remove all leaf nodes).  This step create initial tree for principle graph method.
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

  # remove all leaf nodes
  mst.leaf.node <- mst.count[which(mst.count == 1)]
  mst.leaf.node <- as.integer(names(mst.leaf.node))
  index.1 <- which(mst[,1] %in% mst.leaf.node)
  index.2 <- which(mst[,2] %in% mst.leaf.node)
  index <- unique(c(index.1,index.2))
  mst.trimmed <- mst[-index,]

  # remove all twoway nodes
  mst.all.nodes <- c(mst.trimmed[,1],mst.trimmed[,2])
  mst.count <- table(mst.all.nodes)
  mst.twoway.node <- mst.count[which(mst.count == 2)]
  mst.twoway.node <- as.integer(names(mst.twoway.node))

  while (length(which(mst.count == 2)) > 0) {
    cur.twoway.node <- mst.count[which(mst.count == 2)]
    cur.twoway.node <- as.integer(names(cur.twoway.node))

    cur.twoway.node <- cur.twoway.node[1]

    index.1 <- which(mst.trimmed[,1] %in% cur.twoway.node)
    index.2 <- which(mst.trimmed[,2] %in% cur.twoway.node)
    if(length(index.1) == 2) {
      node1 <- mst.trimmed[index.1[1],2]
      node2 <- mst.trimmed[index.1[2],2]
    } else if (length(index.1) == 1) {
      node1 <- mst.trimmed[index.1,2]
      node2 <- mst.trimmed[index.2,1]
    } else if (length(index.2) == 2){
      node1 <- mst.trimmed[index.2[1],1]
      node2 <- mst.trimmed[index.2[2],1]
    } else {
      stop("Error!")
    }

    index <- unique(c(index.1,index.2))
    mst.trimmed <- mst.trimmed[-index,]

    mst.trimmed <- rbind(mst.trimmed,c(node2, node1))

    mst.all.nodes <- c(mst.trimmed[,1],mst.trimmed[,2])
    mst.count <- table(mst.all.nodes)
  }

  # remove all leaf nodes and two-way nodes from embedding
  embedding.trimmed <- embedding[-c(mst.leaf.node,mst.twoway.node),]

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
