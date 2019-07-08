
#' colorGradient
#'
#' get gradient colors
#'
#' @param x value vector
#' @param colors color plate
#' @param colsteps color steps
#'
colorGradient <- function(
  x,
  colors=c("red","yellow","green","blue"),
  colsteps=100
) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}


#' trajectoryPlotKNN
#'
#' plot cells with KNN method based trajectory
#'
#' @param object seurat object
#' @param assay Assay name
#' @param line.size width of cell trajectory
#' @param line.color color of cell trajectory
#'
#' @export
trajectoryPlotKNN <- function(
  object = NULL,
  assay = 'Joint',
  line.size = 0.5,
  line.color = "gray"
) {
  if(!is.null(object)) {
    group.namae <- as.character(object@misc[[assay]][['cluster']])
    group <- as.character(object@meta.data[[group.namae]])
    graph <- object@misc[[assay]][['knn']]
    reduction <- object@misc[[assay]][['reduction']]
    dim1 = object@reductions[[reduction]]@cell.embeddings[, 1]
    dim2 = object@reductions[[reduction]]@cell.embeddings[, 2]

    # trajectory
    data <- data.frame(
      x=as.numeric(dim1),
      y=as.numeric(dim2),
      cluster=group
    )

    nrow <- dim(graph)[1]
    ncol <- dim(graph)[2]
    row <- rep(1:nrow,ncol)
    graph.frame <- data.frame(row,value =as.numeric(graph))
    graph.frame <- graph.frame[order(row),]
    graph.frame <- as.matrix(graph.frame)

    # initial path.x and path.y
    path.x <- rep(0, dim(graph.frame)[1]*3)
    path.y <- rep(0, dim(graph.frame)[1]*3)

    for (x in 1:dim(graph.frame)[1]) {
      i <- graph.frame[x,1]
      j <- graph.frame[x,2]

      xxx <- 3*(x-1) + 1
      path.x[xxx] <- dim1[i]
      path.x[xxx+1] <- dim1[j]
      path.x[xxx+2] <- 0

      path.y[xxx] <- dim2[i]
      path.y[xxx+1] <- dim2[j]
      path.y[xxx+2] <- NA

    }

    path <- data.frame(x=path.x,y=path.y)
    p <- ggplot(data, aes(x,y)) +
      geom_path(data = path,mapping = aes(x,y), colour = line.color, size = line.size) +
      geom_point(aes(colour = cluster)) + xlab("dim1") + ylab("dim2")

    return(p)
  } else {
    stop("Please provide a Seurat object!")
  }
}



#' pseudoTimePlotKNN
#'
#' plot cells with trajectory and pseudoTime
#'
#' @param object seurat object
#' @param assay Assay name
#' @param colors gradient colors for pseudotime. color order is low to high.
#' @param line.size line width for kNN network
#' @param line.color line color for kNN network
#'
#' @export
#'
pseudoTimePlotKNN <- function(
  object = NULL,
  assay = 'Joint',
  colors = c("blue", "green","yellow", "red"),
  line.size = 0.5,
  line.color = "gray"
) {
  if(!is.null(object)) {

    graph <- object@misc[[assay]][['knn']]
    time.name <- object@misc[[assay]][['time']][['knn']]
    pseudotime <- object@meta.data[[time.name]]
    reduction <- object@misc[[assay]][['reduction']]
    dim1 = object@reductions[[reduction]]@cell.embeddings[, 1]
    dim2 = object@reductions[[reduction]]@cell.embeddings[, 2]

    # joint trajectory
    data <- data.frame(
      x=as.numeric(dim1),
      y=as.numeric(dim2),
      time=as.numeric(pseudotime)
    )

    nrow <- dim(graph)[1]
    ncol <- dim(graph)[2]
    row <- rep(1:nrow,ncol)
    graph.frame <- data.frame(row,value =as.numeric(graph))
    graph.frame <- graph.frame[order(row),]
    graph.frame <- as.matrix(graph.frame)

    # initial path.x and path.y
    path.x <- rep(0, dim(graph.frame)[1]*3)
    path.y <- rep(0, dim(graph.frame)[1]*3)

    for (x in 1:dim(graph.frame)[1]) {
      i <- graph.frame[x,1]
      j <- graph.frame[x,2]

      xxx <- 3*(x-1) + 1
      path.x[xxx] <- dim1[i]
      path.x[xxx+1] <- dim1[j]
      path.x[xxx+2] <- 0

      path.y[xxx] <- dim2[i]
      path.y[xxx+1] <- dim2[j]
      path.y[xxx+2] <- NA

    }

    path <- data.frame(x=path.x,y=path.y)

    p <- ggplot(data, aes(x,y)) +
      geom_path(data = path, mapping = aes(x,y), colour = line.color, size = line.size) +
      geom_point(aes(colour = time)) +
      scale_color_gradientn(colours = colors) + xlab("dim1") + ylab("dim2")
    return(p)
  } else {
    stop("Please provide a Seurat object!")
  }
}

#' trajectoryPlotMST
#'
#' plot cells with trajectory (MST based)
#'
#' @param object seurat object
#' @param assay Assay name
#' @param line.size width of cell trajectory
#' @param line.color color of cell trajectory
#'
#' @export
trajectoryPlotMST <- function(
  object = NULL,
  assay = 'Joint',
  line.size = 0.5,
  line.color = "gray"
) {
  if(!is.null(object)) {
    group.namae <- as.character(object@misc[[assay]][['cluster']])
    group <- as.character(object@meta.data[[group.namae]])
    pg <- object@misc[[assay]][['pg']]
    reduction <- object@misc[[assay]][['reduction']]
    dim1 = object@reductions[[reduction]]@cell.embeddings[, 1]
    dim2 = object@reductions[[reduction]]@cell.embeddings[, 2]

    # trajectory
    data <- data.frame(
      x=as.numeric(dim1),
      y=as.numeric(dim2),
      cluster=group
    )

    node.pg <- pg[[1]]$NodePositions
    tree.pg <- pg[[1]][["Edges"]][["Edges"]]


    path.x <- rep(0, dim(tree.pg)[1]*3)
    path.y <- rep(0, dim(tree.pg)[1]*3)

    for (x in 1:dim(tree.pg)[1]) {
      i <- tree.pg[x,1]
      j <- tree.pg[x,2]

      xxx <- 3*(x-1) + 1
      path.x[xxx] <- node.pg[i,1]
      path.x[xxx+1] <- node.pg[j,1]
      path.x[xxx+2] <- 0

      path.y[xxx] <- node.pg[i,2]
      path.y[xxx+1] <- node.pg[j,2]
      path.y[xxx+2] <- NA
    }

    path <- data.frame(x=path.x,y=path.y)

    p <- ggplot(data, aes(x,y)) +
      geom_point(aes(colour = cluster)) +
      geom_path(data = path, mapping = aes(x,y), colour = line.color, size = line.size) +
      xlab("dim1") + ylab("dim2")

    return(p)
  } else {
    stop("Please provide a Seurat object!")
  }
}



#' pseudoTimePlotMST
#'
#' plot cells with trajectory and pseudoTime (MST based)
#'
#' @param object seurat object
#' @param assay reduction name
#' @param colors gradient colors for pseudotime. color order is low to high.
#' @param line.size line width for trajectory
#' @param line.color line color for trajectory
#'
#' @export
pseudoTimePlotMST <- function(
  object = NULL,
  assay = NULL,
  colors = c("blue", "green","yellow", "red"),
  line.size = 0.5,
  line.color = "gray"
) {
  if(!is.null(object)) {
    pg <- object@misc[[assay]][['pg']]
    reduction <- object@misc[[assay]][['reduction']]
    dim1 = object@reductions[[reduction]]@cell.embeddings[, 1]
    dim2 = object@reductions[[reduction]]@cell.embeddings[, 2]
    time.name <- object@misc[[assay]][['time']][['mst']]
    pseudotime <- object@meta.data[[time.name]]

    # joint trajectory
    data <- data.frame(
      x=as.numeric(dim1),
      y=as.numeric(dim2),
      time=as.numeric(pseudotime)
    )

    node.pg <- pg[[1]]$NodePositions
    tree.pg <- pg[[1]][["Edges"]][["Edges"]]


    path.x <- rep(0, dim(tree.pg)[1]*3)
    path.y <- rep(0, dim(tree.pg)[1]*3)

    for (x in 1:dim(tree.pg)[1]) {
      i <- tree.pg[x,1]
      j <- tree.pg[x,2]

      xxx <- 3*(x-1) + 1
      path.x[xxx] <- node.pg[i,1]
      path.x[xxx+1] <- node.pg[j,1]
      path.x[xxx+2] <- 0

      path.y[xxx] <- node.pg[i,2]
      path.y[xxx+1] <- node.pg[j,2]
      path.y[xxx+2] <- NA
    }

    path <- data.frame(x=path.x,y=path.y)

    p <- ggplot(data, aes(x,y)) +
      geom_point(aes(colour = time)) + scale_color_gradientn(colours = colors) +
      geom_path(data = path, mapping = aes(x,y), colour = line.color, size = line.size) +
      xlab("dim1") + ylab("dim2")

    return(p)
  } else {
    stop("Please provide a Seurat object!")
  }
}


#' clusterConnectionPlot
#'
#' plot cell clusters with connections
#'
#' @param object seurat object
#' @param assay Assay name
#' @param cutoff only connections between two cluster have connection factor larger than this cutoff will be shown. Default value = 0, user can adjust the cutoff to hide weak connections
#' @param line.color line color of the connections between cell clusters
#' @param point.size control dot size of cell clusters
#' @param path.size control path size of cell connections
#'
#' @export
clusterConnectionPlot <- function(
  object = NULL,
  assay = 'Joint',
  cutoff = 0,
  line.color = "black",
  point.size = 2,
  path.size = 1.5
) {
  if(!is.null(object)) {
    group.namae <- as.character(object@misc[[assay]][['cluster']])
    group <- as.character(object@meta.data[[group.namae]])
    graph <- object@misc[[assay]][['knn']]
    reduction <- object@misc[[assay]][['reduction']]
    dim1 = object@reductions[[reduction]]@cell.embeddings[, 1]
    dim2 = object@reductions[[reduction]]@cell.embeddings[, 2]

    data <- data.frame(
      x=as.numeric(dim1),
      y=as.numeric(dim2),
      cluster=group
    )

    clusters <- sort(unique(group))
    nCluster <- length(unique(group))
    cluster.center <- data.frame(dim1 = 0, dim2 = 0, cells = 0)

    for (i in 1:nCluster) {
      cur.cluster <- clusters[i]
      sub.data <- data[which(data$cluster == cur.cluster),]
      cluster.center[i,] <- c(mean(sub.data$x), mean(sub.data$y), dim(sub.data)[1])
    }
    cluster.center$cluster <- clusters
    cluster.center$dim1 <- as.numeric(cluster.center$dim1)
    cluster.center$dim2 <- as.numeric(cluster.center$dim2)
    cluster.center$cells <- as.numeric(cluster.center$cells)

    # evaluate connections among cell clusters
    knn.matrix <- matrix(data = 0, nrow = dim(graph)[1], ncol = dim(graph)[1])
    for (i in 1:dim(graph)[1]) {
      j.index <- graph[i,]
      knn.matrix[i,j.index] <- 1
    }

    cluster.path <- data.frame(dim1 = 0, dim2 = 0, connection = 0)
    cluster.path.index <- 1

    for (i in 1:(nCluster-1)) {
      c1 <- clusters[i]
      index.c1 <- which(group == c1)
      for (j in (i+1):nCluster) {
        c2 <- clusters[j]
        index.c2 <- which(group == c2)
        within.c1 <- length(which(knn.matrix[index.c1,index.c1] > 0))
        within.c2 <- length(which(knn.matrix[index.c2,index.c2] > 0))
        across <- length(which(knn.matrix[index.c1,index.c2] > 0))

        connection.factor <- across/sqrt(within.c1*within.c2)

        if(connection.factor > cutoff) {
          where.c1 <- which(cluster.center$cluster == c1)
          where.c2 <- which(cluster.center$cluster == c2)

          cluster.path[cluster.path.index,] <- c(cluster.center$dim1[where.c1], cluster.center$dim2[where.c1], connection.factor)
          cluster.path.index <- cluster.path.index + 1
          cluster.path[cluster.path.index,] <- c(cluster.center$dim1[where.c2], cluster.center$dim2[where.c2], connection.factor)
          cluster.path.index <- cluster.path.index + 1
          cluster.path[cluster.path.index,] <- c(0, NA, 0)
          cluster.path.index <- cluster.path.index + 1
        }
      }
    }

    p <- ggplot(cluster.center,aes(x=dim1,y=dim2)) +
      geom_path(data = cluster.path, size= path.size*40*cluster.path$connection, colour = line.color) +
      geom_point(aes(colour = cluster), size = point.size*log2(cluster.center$cells)) +
      geom_text(aes(label = cluster))

    return(p)
  } else {
    stop("Please provide a Seurat object!")
  }
}



#' heatMapPlot
#'
#' plot heatmap of ADT and selected RNA features of cell clusters
#'
#' @param object seurat object
#' @param adt.feature ADT features used for plot
#' @param rna.feature RNA features used for plot
#' @param group.by cell clusters
#' @param adt.label annotate ADT features on HEATmap (TRUE or FALSE)
#' @param rna.label annotate RNA features on HEATmap (TRUE or FALSE)
#' @param label.hjust hjust for ADT labels, number between 0 and 1
#' @param label.size font size of ADT labels
#' @param height.rel reletive height of RNA compare to ADT. Default 1 means equal height. set to null , the height will be set according to number of features
#' @param legend show legend or not. TRUE or FALSE
#'
#' @export

heatMapPlot <- function(
  object = NULL,
  adt.feature = NULL,
  rna.feature = NULL,
  group.by = NULL,
  adt.label = TRUE,
  rna.label = FALSE,
  label.hjust = 0.1,
  label.size = 5,
  height.rel = 1,
  legend = FALSE
) {
  if(!is.null(object)) {
    if(!is.null(group.by)){
      if(is.null(adt.feature)) {
        cat("User didn't provide ADT feature list, will use all!\n")
        adt.feature = rownames(object@assays[["ADT"]]@scale.data)
      }

      if(is.null(rna.feature)) {
        cat("User didn't provide RNA feature list, will choose top 10 for each cluster!\n")

        Idents(object = object) <- object[[group.by]]
        rna.markers <- FindAllMarkers(object, max.cells.per.ident = 100, min.diff.pct = 0.3, only.pos = TRUE,assay = 'RNA')
        rna.top10 <- rna.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
        rna.feature <- rna.top10
        rna.feature <- unique(rna.feature$gene)
      } else {
        rna.feature <- unique(rna.feature$gene)
      }

      ph.adt <- DoHeatmap(object, features = adt.feature, assay = "ADT", angle = 90, group.by = group.by) + theme(axis.text.y = element_blank())
      if(!isTRUE(legend)) {
        ph.adt <- ph.adt + NoLegend()
      }
      if(isTRUE(adt.label)) {
        if(!is.null(label.hjust)) {
          if((label.hjust > 1) || (label.hjust < 0)) {
            cat("label.hjust only can be number between 0 and 1, will use default value 0.1!\n")
            label.hjust = 0.1
          }
          label.hjust <- dim(object)[2]*label.hjust
          ph.adt <- ph.adt + annotate("text",label = adt.feature, x = label.hjust, y = length(adt.feature):1,colour="white", fontface =2, size = label.size)
        } else {
          label.hjust <- dim(object)[2]*0.1
          ph.adt <- ph.adt + annotate("text",label = adt.feature, x = label.hjust, y = length(adt.feature):1,colour="white", fontface =2, size = label.size)
        }
      }

      ph.rna <- DoHeatmap(object, features = rna.feature, assay = "RNA",group.by = group.by, group.bar = FALSE, label = FALSE) + theme(axis.text.y = element_blank())
      if(!isTRUE(legend)) {
        ph.rna <- ph.rna + NoLegend()
      }
      rna.feature <- levels(ph.rna[["data"]][["Feature"]])
      if(isTRUE(rna.label)) {
        if(!is.null(label.hjust)) {
          if((label.hjust > 1) || (label.hjust < 0)) {
            cat("label.hjust only can be number between 0 and 1, will use default value 0.1!\n")
            label.hjust = 0.1
          }
          label.hjust <- dim(object)[2]*label.hjust
          ph.rna <- ph.rna + annotate("text",label = rna.feature, x = label.hjust, y = length(rna.feature):1,colour="white", fontface =2, size = label.size)
        } else {
          label.hjust <- dim(object)[2]*0.1
          ph.rna <- ph.rna + annotate("text",label = rna.feature, x = label.hjust, y = length(rna.feature):1,colour="white", fontface =2, size = label.size)
        }
      }
      if(is.null(height.rel)) {
        height.rel <- length(rna.feature)/length(adt.feature)
      }
      rel_size_height <- c(1,height.rel)
      p <- plot_grid(ph.adt,ph.rna,ncol = 1,rel_heights = rel_size_height)

      return(p)
    } else {
      stop("Please determine cell cluster!\n")
    }
  } else {
    stop("Please provide a seurat object!\n")
  }
}

#' gridDimPlot
#'
#' grid dim plot for RNA, ADT and joint analysis
#'
#' @param object seurat object
#' @param assays assays user want to plot and compare
#' @param darkTheme switch for darkTheme (TRUE or FALSE)
#' @param legend switch for showing legend in each sub plot (TRUE or FALSE)
#' @param figure.label.size font size for figure labels
#' @param cluster.label swicth for showing cluster label on each sub plot (TRUE or FALSE)
#' @param cluster.lable.size font size for cluster labels
#' @param wide.rel relative width of figure label to the left. size of plots is 3.
#' @param height.rel relative height of figure label on the top. size of plots is 3.
#'
#' @export

gridDimPlot <- function(
  object = NULL,
  assays = c("RNA","ADT","Joint"),
  darkTheme = TRUE,
  legend = TRUE,
  figure.label.size = 7,
  cluster.label = TRUE,
  cluster.lable.size = 5,
  wide.rel = 1.5,
  height.rel = 1
) {
  if(!is.null(object)) {
    empty <- ggplot() + annotate("text",label = "", x = 0, y = 0, size = figure.label.size) + theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line = element_blank()
    )
    p.list <- list(empty)
    # for label on the top
    for (assay in assays) {
      cur.label <- paste0(assay, " based clustering")
      top.label.figure <- ggplot() + annotate("text",label = cur.label, x = 0, y = 0, size = figure.label.size, colour = "black") + theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_blank()
      )
      cur.index <- length(p.list) + 1
      p.list[[cur.index]] <- top.label.figure
    }

    for (assay.x in assays) {
      # for label to the left
      cur.label <- paste0(assay.x, " based map")
      top.label.figure <- ggplot() + annotate("text",label = cur.label, x = 0, y = 0, size = figure.label.size, colour = "black") + theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_blank()
      )
      cur.index <- length(p.list) + 1
      p.list[[cur.index]] <- top.label.figure

      # dim plots for assay.x
      for (assay.y in assays) {
        cur.cluster <- paste0(tolower(assay.y),cluster.suffix)
        cur.map <- paste0(reduction.prefix,tolower(assay.x))

        cur.p <- DimPlot(object = object, group.by = cur.cluster,reduction = cur.map) + xlab(paste0(reduction.prefix,"1")) + ylab(paste0(reduction.prefix,"2"))
        if(!isTRUE(legend)) {
          cur.p <- cur.p + NoLegend()
        }
        if(isTRUE(darkTheme)) {
          cur.p <- cur.p + DarkTheme()
        }
        if(isTRUE(cluster.label)) {
          if(isTRUE(darkTheme)) {
            cur.p <- LabelClusters(plot = cur.p, id = cur.cluster, size = cluster.lable.size, color = "white")
          } else {
            cur.p <- LabelClusters(plot = cur.p, id = cur.cluster, size = cluster.lable.size, color = "black")
          }
        }
        cur.index <- length(p.list) + 1
        p.list[[cur.index]] <- cur.p
      }
    }

    # determine the relevent size of figures (set label:plot = 1:3)
    rel_size <- c(3,3,3,3,3,3)
    rel_size <- rel_size[1:length(assays)]
    rel_size_wide <- c(wide.rel,rel_size)
    rel_size_height <- c(height.rel,rel_size)

    plot <- plot_grid(plotlist = p.list, ncol=length(assays) + 1, rel_widths = rel_size_wide, rel_heights = rel_size_height )
    return(plot)
  } else {
    stop("Please provide a seurat object!\n")
  }
}



#' generateGridDimPlot
#'
#' generate a grid dim plot object list for RNA, ADT and joint analysis
#'
#' @param object seurat object
#' @param assays assays user want to plot and compare
#' @param darkTheme switch for darkTheme (TRUE or FALSE)
#' @param legend switch for showing legend in each sub plot (TRUE or FALSE)
#' @param cluster.label swicth for showing cluster label on each sub plot (TRUE or FALSE)
#' @param cluster.lable.size font size for cluster labels
#'
#' @export

generateGridDimPlot <- function(
  object = NULL,
  assays = c("RNA","ADT","Joint"),
  darkTheme = TRUE,
  legend = TRUE,
  cluster.label = TRUE,
  cluster.lable.size = 5
) {
  if(!is.null(object)) {
    p.list <- list()

    for (assay.x in assays) {
      # dim plots for assay.x
      for (assay.y in assays) {
        cur.cluster <- object@misc[[assay.y]][['cluster']]
        cur.map <- object@misc[[assay.x]][['reduction']]
        reductionType <- gsub(pattern = "_.+", replacement = "", x = cur.map)
        cur.p <- DimPlot(object = object, group.by = cur.cluster,reduction = cur.map) + xlab(paste0(reductionType," 1")) + ylab(paste0(reductionType," 2"))
        if(!isTRUE(legend)) {
          cur.p <- cur.p + NoLegend()
        }
        if(isTRUE(darkTheme)) {
          cur.p <- cur.p + DarkTheme()
        }
        if(isTRUE(cluster.label)) {
          if(isTRUE(darkTheme)) {
            cur.p <- LabelClusters(plot = cur.p, id = cur.cluster, size = cluster.lable.size, color = "white")
          } else {
            cur.p <- LabelClusters(plot = cur.p, id = cur.cluster, size = cluster.lable.size, color = "black")
          }
        }

        info <- list("ident" = assay.y, "map" = assay.x)
        cur.p$info <- info
        cur.index <- length(p.list) + 1
        p.list[[cur.index]] <- cur.p
      }
    }
    return(p.list)
  } else {
    stop("Please provide a seurat object!\n")
  }
}

#' listPlot
#'
#' plot grid dim plots. the object should be generated by generateGridDimPlot() function
#'
#' @param object list of dim plots
#' @param align parameter of plot_grid() from cowplot package.(optional) Specifies whether graphs in the grid should be horizontally ("h") or vertically ("v") aligned. Options are "none" (default), "hv" (align in both directions), "h", and "v".
#' @param axis parameter of plot_grid() from cowplot package.(optional) Specifies whether graphs should be aligned by the left ("l"), right ("r"), top ("t"), or bottom ("b") margins. Options are "none" (default), or a string of any combination of l, r, t, and b in any order (e.g. "tblr" or "rlbt" for aligning all margins). Must be specified if any of the graphs are complex (e.g. faceted) and alignment is specified and desired. See align_plots() for details.
#' @param ncol parameter of plot_grid() from cowplot package.(optional) Number of columns in the plot grid.
#' @param fig.id Numerical ids for selected plots. could be number or vector
#' @param fig.ident ident property of selected plots. For example, set to "RNA" will select and plot all plots that have RNA clusters
#' @param fig.map map property of selected plots. For example, set to "RNA" will select and plot all plots that have RNA maps
#' @param rel_widths parameter of plot_grid() from cowplot package. (optional) Numerical vector of relative columns widths. For example, in a two-column grid, rel_widths = c(2, 1) would make the first column twice as wide as the second column.
#' @param rel_heights parameter of plot_grid() from cowplot package.(optional) Numerical vector of relative columns heights. Works just as rel_widths does, but for rows rather than columns.
#' @param labels parameter of plot_grid() from cowplot package.(optional) List of labels to be added to the plots. You can also set labels="AUTO" to auto-generate upper-case labels or labels="auto" to auto-generate lower-case labels.
#' @param label_size parameter of plot_grid() from cowplot package.(optional) Numerical value indicating the label size. Default is 14.
#' @param label_fontfamily parameter of plot_grid() from cowplot package.(optional) Font family of the plot labels. If not provided, is taken from the current theme.
#' @param label_fontface parameter of plot_grid() from cowplot package.(optional) Font face of the plot labels. Default is "bold".
#' @param label_colour parameter of plot_grid() from cowplot package.(optional) Color of the plot labels. If not provided, is taken from the current theme.
#' @param label_x parameter of plot_grid() from cowplot package.(optional) Single value or vector of x positions for plot labels, relative to each subplot. Defaults to 0 for all labels. (Each label is placed all the way to the left of each plot.)
#' @param label_y parameter of plot_grid() from cowplot package.(optional) Single value or vector of y positions for plot labels, relative to each subplot. Defaults to 1 for all labels. (Each label is placed all the way to the top of each plot.)
#' @param hjust parameter of plot_grid() from cowplot package.Adjusts the horizontal position of each label. More negative values move the label further to the right on the plot canvas. Can be a single value (applied to all labels) or a vector of values (one for each label). Default is -0.5.
#' @param vjust parameter of plot_grid() from cowplot package.Adjusts the vertical position of each label. More positive values move the label further down on the plot canvas. Can be a single value (applied to all labels) or a vector of values (one for each label). Default is 1.5.
#' @param scale parameter of plot_grid() from cowplot package.Individual number or vector of numbers greater than 0. Enables you to scale the size of all or select plots. Usually it's preferable to set margins instead of using scale, but scale can sometimes be more powerful.
#'
#' @export

listPlot <- function(
  object = NULL,
  align = "none",
  axis = "none",
  fig.id = NULL,
  fig.ident = NULL,
  fig.map = NULL,
  ncol = NULL,
  rel_widths = 1,
  rel_heights = 1,
  labels = NULL,
  label_size = 10,
  label_fontfamily = NULL,
  label_fontface = "bold",
  label_colour = NULL,
  label_x = 0.1,
  label_y = 1,
  hjust = -0.5,
  vjust = 1,
  scale = 1
) {
  if(!is.null(object)) {
    if(!is.list(object)) {
      stop("Please provide a list of plot object!\n")
    } else {
      if(!is.null(fig.id)) {
        object <- object[fig.id]
      }

      if(!is.null(fig.ident)) {
        temp <- list()
        for (plot in object) {
          if(plot$info[["ident"]] == fig.ident) {
            index <- length(temp) + 1
            temp[[index]] <- plot
          }
        }
        object <- temp
      }

      if(!is.null(fig.map)) {
        temp <- list()
        for (plot in object) {
          if(plot$info[["map"]] == fig.map) {
            index <- length(temp) + 1
            temp[[index]] <- plot
          }
        }
        object <- temp
      }

      if(is.null(ncol)) {
        n <- length(object)
        n <- floor(sqrt(n))
        ncol <- n
      }

      if(is.null(labels)) {
        labels <- c()
        for (p in object) {
          info <- p$info
          title <- paste0(info[["ident"]]," cluster + ",info[["map"]]," map")
          labels <- c(labels, title)
        }
      }

      plots <- plot_grid(
        plotlist = object,
        ncol = ncol,
        align = align,
        axis = axis,
        rel_widths = rel_widths,
        rel_heights = rel_heights,
        labels = labels,
        label_size = label_size,
        label_fontfamily = label_fontfamily,
        label_fontface = label_fontface,
        label_colour = label_colour,
        label_x = label_x,
        label_y = label_y,
        hjust = hjust,
        vjust = vjust,
        scale = scale
        )
      return(plots)
    }
  } else {
    stop("Please provide a list of plot object!\n")
  }
}


#' plotInfo
#'
#' print figure information for all plots
#'
#' @param object list of ggplot2 plots
#'
#' @export


plotInfo <- function (
  object
) {
  if(!is.null(object)) {
    info <- data.frame(
      "index" = NA,
      "ident" = NA,
      "map" = NA
        )
    i = 1
    for (p in object) {
      temp <- c(i ,p$info[["ident"]] ,p$info[["map"]])
      info[i,] <- temp
      i <- i + 1
    }
    return(info)
  } else {
    stop("Please provide a list of plot object!\n")
  }
}


#' highCorrelatedGenePlot
#'
#' order all cells by the gene expression level of a specific gene, and explorer the gene expression profiles for highly corrolated genes or a group of user choosed genes
#'
#' @param object Seurat object
#' @param target.gene target gene user want focus on
#' @param coef coef methods: spearman, pearson or kendall
#' @param genes user's choice gene list
#' @param topx number of genes that highly correlated with target gene
#'
#' @export

highCorrelatedGenePlot <- function(
  object = NULL,
  target.gene = NULL,
  coef = "spearman",
  genes = NULL,
  topx = 20,
  group.by = "ident"
) {
  if(!is.null(object)) {
    expression <- GetAssayData(object = object, slot = "scale.data", assay = "RNA")

    # Order cells according to target gene
    target.gene.express <- expression[target.gene,]
    cell.order <- names(sort(target.gene.express))
    target.gene.express <- target.gene.express[cell.order]
    target.gene.express <- as.data.frame(target.gene.express)
    colnames(target.gene.express) <- c("exp")
    target.gene.express$order <- 1:length(cell.order)

    if(!is.null(genes)) {
      # make corrosponding RNA expression data
      gene.list <- genes[which(genes %in% rownames(expression))]
      data <- expression[gene.list,]
      list <- setdiff(genes,rownames(expression))
      zero.matrix <- as.data.frame(matrix(0,length(list),length(cell.order)))
      rownames(zero.matrix) <- list
      colnames(zero.matrix) <- colnames(data)
      data <- rbind(data,zero.matrix)
      data <- data[,cell.order]
      data <- t(data)
    } else {
      all.coef <- c()
      for (i in 1:dim(expression)[1]) {
        all.coef[i] <- cor(as.numeric(expression[target.gene,]),as.numeric(expression[i,]),method = coef)
      }

      coef.data <- data.frame(gene=rownames(expression),cc=all.coef)
      coef.data <- coef.data[order(coef.data$cc, decreasing = TRUE),]
      coef.data <- coef.data[1:topx,]

      genes <- as.character(coef.data$gene)

      gene.list <- genes[which(genes %in% rownames(expression))]
      data <- expression[gene.list,]
      list <- setdiff(genes,rownames(expression))
      zero.matrix <- as.data.frame(matrix(0,length(list),length(cell.order)))
      rownames(zero.matrix) <- list
      colnames(zero.matrix) <- colnames(data)
      data <- rbind(data,zero.matrix)
      data <- data[,cell.order]
      data <- t(data)
    }

    if(group.by == "ident") {
      cluster <- as.data.frame(object@active.ident)
      colnames(cluster) <- c("cluster")
    } else {
      cluster <- as.data.frame(object@meta.data[[group.by]])
      colnames(cluster) <- c("cluster")
    }

    # plot
    p_line <- ggplot(data=target.gene.express, aes(x=order,y=exp)) + geom_line(arrow = arrow(),color="red") + labs(x="", y = "Gene Expression") + theme(axis.text.x = element_blank())
    colfunc <- colorRampPalette(c("gray", "red"))
    # expression
    plot <- Seurat:::SingleRasterMap(data = data, raster = TRUE, colors = colfunc(20),
                            disp.min = -2.5, disp.max = 2.5, feature.order = genes,
                            cell.order = cell.order) + NoLegend()

    # plot group bars
    pbuild <- ggplot_build(plot = plot)
    cols <- hue_pal()(length(x = levels(x = cluster$cluster)))
    names(x = cols) <- levels(x = cluster$cluster)
    y.pos <- max(pbuild$layout$panel_params[[1]]$y.range)*1.01
    y.max <- max(pbuild$layout$panel_params[[1]]$y.range)*1.03
    plot <- plot + annotation_raster(raster = t(x = cols[cluster$cluster]), xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) + coord_cartesian(ylim = c(0, y.max), clip = "off")
    figure <- plot_grid(p_line, plot,ncol = 1)
    return(figure)
  } else {
    stop("Please provide Seurat Object!")
  }
}



# Set a default value if an object is null
#
# @param lhs An object to set if it's null
# @param rhs The value to provide if x is null
#
# @return rhs if lhs is null, else lhs
#
# @author Hadley Wickham
# @references https://adv-r.hadley.nz/functions.html#missing-arguments
#
`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

# Set a default value if an object is NOT null
#
# @param lhs An object to set if it's NOT null
# @param rhs The value to provide if x is NOT null
#
# @return lhs if lhs is null, else rhs
#
# @author Hadley Wickham
# @references https://adv-r.hadley.nz/functions.html#missing-arguments
#
`%iff%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(rhs)
  } else {
    return(lhs)
  }
}

# Melt a data frame
#
# @param x A data frame
#
# @return A molten data frame
#
Melt <- function(x) {
  if (!is.data.frame(x = x)) {
    x <- as.data.frame(x = x)
  }
  return(data.frame(
    rows = rep.int(x = rownames(x = x), times = ncol(x = x)),
    cols = unlist(x = lapply(X = colnames(x = x), FUN = rep.int, times = nrow(x = x))),
    vals = unlist(x = x, use.names = FALSE)
  ))
}
