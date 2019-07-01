#' readDataFrom10X
#'
#' read data from 10X (cell ranger V3 output) by calling Seurat Read10X function
#'
#' @param dir path of 10X cell ranger output folder
#' @param gene.column Specify which column of genes.tsv or features.tsv to use for gene names; default is 2
#' @param unique.features Make feature names unique (default TRUE)
#'
#' @export
readDataFrom10X <- function(
  dir = NULL,
  gene.column = 2,
  unique.features = TRUE
  ) {
  if(!is.null(dir))
  {
    if(dir.exists(dir)){
      raw.data <- Read10X(data.dir = dir, gene.column = gene.column, unique.features = unique.features)
      rownames(x = raw.data[["Antibody Capture"]]) <- gsub(pattern = "_.+", replacement = "", x = rownames(x = raw.data[["Antibody Capture"]]))
      names(raw.data) <- c("GeneExpression","AntibodyCapture")
      return(raw.data)
    } else {
      stop("The file path you provided is not exist! Please check your input!")
    }
  } else {
    stop("Please provide the path of cell ranger output!")
  }
}


#' createObject
#'
#' Create Seurat V3 object from raw data by calling Seurat CreateSeuratObject function and CreateAssayObject function
#'
#' @param data data of 10X cell ranger output
#' @param project project name for current dataset
#' @param assay assay name, default is RNA
#' @param min.cells Same as CreateSeuratObject function in Seurat. Include features detected in at least this many cells. Will subset the counts matrix as well. To reintroduce excluded features, create a new object with a lower cutoff.
#' @param min.features Same as CreateSeuratObject function in Seurat. Include cells where at least this many features are detected.
#' @param names.field Same as CreateSeuratObject function in Seurat. For the initial identity class for each cell, choose this field from the cell's name. E.g. If your cells are named as BARCODE_CLUSTER_CELLTYPE in the input matrix, set names.field to 3 to set the initial identities to CELLTYPE.
#' @param names.delim Same as CreateSeuratObject function in Seurat. For the initial identity class for each cell, choose this delimiter from the cell's column name. E.g. If your cells are named as BARCODE-CLUSTER-CELLTYPE, set this to "-" to separate the cell name into its component parts for picking the relevant field.
#' @param meta.data Same as CreateSeuratObject function in Seurat. Additional cell-level metadata to add to the Seurat object. Should be a data frame where the rows are cell names and the columns are additional metadata fields.
#' @param min.cells.adt Same as CreateAssayObject function in Seurat. Include ADT features detected in at least this many cells. Will subset the counts matrix as well. To reintroduce excluded features, create a new object with a lower cutoff.
#' @param min.features.adt Same as CreateAssayObject function in Seurat. Include cells where at least this many ADT features are detected.
#'
#' @export
createObject <- function(
  data = NULL,
  project = "SeuratProject",
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_",
  meta.data = NULL,
  min.cells.adt = 0,
  min.features.adt = 0,
  rna = NULL,
  adt = NULL
  ) {
  if(!is.null(data)) {
    if(is.null(data[["AntibodyCapture"]])){
      # Create Seurat object from RNA data matrix (for single modal data)
      object <- CreateSeuratObject(counts = data, project = project, assay = "RNA", min.cells = min.cells, min.features = min.features, names.field = names.field, names.delim = names.delim, meta.data = meta.data)
    } else {
      # Create Seurat object from RNA and ADT data (for multi modal data, cell ranger V3 output)
      object <- CreateSeuratObject(counts = data[["GeneExpression"]], project = project, assay = "RNA", min.cells = min.cells, min.features = min.features, names.field = names.field, names.delim = names.delim, meta.data = meta.data)
      object[["ADT"]] <- CreateAssayObject(data[["AntibodyCapture"]][, colnames(x = object)],min.cells = min.cells.adt, min.features = min.features.adt)
    }
    # calculate the percentage of mito genes for later cellQC analysis
    object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
    return(object)
  } else {
    if( !is.null(rna) && !is.null(adt) ) {
      # sss
      index <- intersect(colnames(rna), colnames(adt))
      rna <- as.matrix(rna[, index])
      adt <- as.matrix(adt[, index])

      object <- CreateSeuratObject(counts = rna, project = project, assay = "RNA", min.cells = min.cells, min.features = min.features, names.field = names.field, names.delim = names.delim, meta.data = meta.data)
      object[["ADT"]] <- CreateAssayObject(adt,min.cells = min.cells.adt, min.features = min.features.adt)
      object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
    } else if (!is.null(rna)) {
      object <- CreateSeuratObject(counts = rna, project = project, assay = "RNA", min.cells = min.cells, min.features = min.features, names.field = names.field, names.delim = names.delim, meta.data = meta.data)
      object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
    } else if (!is.null(adt)) {
      object <- CreateSeuratObject(counts = rna, project = project, assay = "ADT", min.cells = min.cells, min.features = min.features, names.field = names.field, names.delim = names.delim, meta.data = meta.data)
    } else {
      stop("Please provide input data through either data = or rna = , adt = ,hto = \n")
    }
    return(object)
  }
}



#' loadHashTag
#'
#' load cell hash tag information to a Seurat object from a CSV or TXT file
#'
#' @param object Seurat object
#' @param file hash tag csv file
#'
#' @export

loadHashTag <- function(
  object = NULL,
  file = NULL
) {
  if(is.null(object)) {
    stop("Please provide Seurat objects!\n")
  }
  if(!is.null(file)) {
    table <- read.csv(file = file, header = TRUE, row.names = 1)
    remove.index <- which(rownames(table) %in% c("bad_struct","no_match", "total_reads"))
    table <- table[-remove.index,]

    share.names <- intersect(colnames(object), colnames(table))
    object <- object[,share.names]
    table <- table[,share.names]

    object[["HTO"]] <- CreateAssayObject(table)
    return(object)
  } else {
    stop("Please provide HTO file path!\n")
  }
}


#' loadSgRNA
#'
#' load SgRNA information to a Seurat object from a CSV or TXT file
#'
#' @param object Seurat object
#' @param file hash tag csv file
#'
#' @export

loadSgRNA <- function(
  object = NULL,
  file = NULL
) {
  if(is.null(object)) {
    stop("Please provide Seurat objects!")
  }
  if(!is.null(file)) {
    table <- read.csv(file = file, header = TRUE, row.names = 1)
    remove.index <- which(rownames(table) %in% c("bad_struct","no_match", "total_reads"))
    table <- table[-remove.index,]

    share.names <- intersect(colnames(object), colnames(table))
    object <- object[,share.names]
    table <- table[,share.names]

    object[["GDO"]] <- CreateAssayObject(table)
    return(object)
  } else {
    stop("Please provide GDO file path!")
  }
}


#' loadADT
#'
#' load ADT information to a Seurat object from a CSV or TXT file
#'
#' @param object Seurat object
#' @param file hash tag csv file
#'
#' @export

loadADT <- function(
  object = NULL,
  file = NULL
) {
  if(is.null(object)) {
    stop("Please provide Seurat objects!")
  }
  if(!is.null(file)) {
    table <- read.csv(file = file, header = TRUE, row.names = 1)
    remove.index <- which(rownames(table) %in% c("bad_struct","no_match", "total_reads"))
    table <- table[-remove.index,]

    share.names <- intersect(colnames(object), colnames(table))
    object <- object[,share.names]
    table <- table[,share.names]

    object[["ADT"]] <- CreateAssayObject(table)
    return(object)
  } else {
    stop("Please provide ADT file path!")
  }
}


#' readTCR
#'
#' read TCR/BCR information from a CSV file to a table
#'
#' @param file hash tag csv file
#'
#' @export

readTCR <- function(
  file = NULL
) {
  if(!is.null(file)) {
    table <- read.csv(file = file, header = TRUE)
    barcodes <- table$barcode
    table$barcode <- sub("-1","",barcodes)
    return(table)
  } else {
    stop("Please provide TCR/BCR csv file path!")
  }
}


#' loadTCR (not finished)
#'
#' load TCR/BCR information from a table to a seurat object. Only information in user selected fields will be loaded into seurat object.
#'
#' @param object Seurat object
#' @param file TCR/BCR information file path
#' @param table TCR/BCR information table
#' @param field hash tag csv file
#'
#' @export

loadTCR <- function(
  object = NULL,
  file = NULL,
  table = NULL,
  field = NULL
) {
  if(is.null(object)) {
    stop("Please provide Seurat objects!")
  }
  if(!is.null(field)) {
    if(!is.null(file) && !is.null(table)) {
      cat("You provided both path and table, will use information from table!")
      table <- table
    } else if (!is.null(file)) {
      table <- read.csv(file = file, header = TRUE)
      barcodes <- table$barcode
      table$barcode <- sub("-1","",barcodes)
    } else if (!is.null(table)) {

    }else {
      stop("Please provide TCR/BCR csv file path or information table!")
    }
    table
  } else {
    stop("Please provide at least one field name!")
  }
}


#' attachModal (not finished)
#'
#' read modalities and attach to cells
#'
#' @param object Seurat object
#' @param data data table contains modality information with barcode
#' @param target.col.name ssss
#' @param type could be "bool" or "value"
#' @param barcode.index index of barcode column in data table
#'
#' @export

attachModal <- function(
  object = NULL,
  data = NULL,
  target.col.name = NULL,
  type = c("bool","value"),
  barcode.index = 1
) {
  if(!is.null(object)) {
    if(!is.null(data)) {
      barcodes.all <- rownames(object@meta.data)
      if(type == "bool") {
        new.meta <- rep("Negative",length(barcodes.all))
      } else if (type == "value") {
        new.meta <- rep(NA,length(barcodes.all))
      } else {
        stop("type can only be bool or value")
      }
    } else {
      stop("Please provide data table!")
    }
  } else {
    stop("Please provide Seurat object!")
  }
}

