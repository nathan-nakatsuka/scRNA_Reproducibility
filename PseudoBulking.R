# This is a function to give back labels to an object after you pseudobulked it.
AverageMetaData <- function(orig, avg, f = NULL) {
  f <- f %||% Idents(object = orig)
  f <- as.character(x = f)
  if (!all(colnames(x = avg) %in% unique(x = f))) {
    stop("Not all averaged groups present in original object", call. = FALSE)
  }
  groups <- colnames(x = avg)
  md <- orig[[]]
  for (col in names(x = md)) {
    m <- split(x = md[[col]], f = f)[groups]
    m <- sapply(X = m, FUN = unique, simplify = FALSE, USE.NAMES = TRUE)
    if (!all(vapply(X = m, FUN = length, FUN.VALUE = integer(length = 1L)) == 1L)) {
      warning(
        "Meta data ",
        sQuote(x = col),
        " is not unique across groups",
        call. = FALSE,
        immediate. = TRUE
      )
      next
    }
    message("Adding meta data ", sQuote(x = col), " to averaged object")
    m <- data.frame(unlist(x = m), row.names = names(x = m))
    names(x = m) <- col
    avg[[col]] <- m
  }
  return(avg)
}

# This is a function to pseudobulk objects and give them back their original labels.
# CellTypeLevel indicates the cell type level resolution (e.g. "predicted.celltype.l1")
# Note, this requires the individual ID to be labeled with the column name "patient".
# Note, this is an aggregate pseudobulk for DESeq2. If you want to get the mean use the PseudobulkSeuratObject_Mean function
PseudobulkSeuratObject_Aggregate <- function(SeuratObject, CellTypeLevel){
   DefaultAssay(SeuratObject) <- "RNA"
   SeuratObject <- NormalizeData(SeuratObject)
   avg_exp_SeuratObject <- Seurat:::PseudobulkExpression(SeuratObject, slot="counts",pb.method='aggregate', group.by = c(CellTypeLevel, "patient"), return.seurat = T)
   #Give back labels to the object
   Idents(SeuratObject) <- CellTypeLevel
   avg_exp_SeuratObject <- AverageMetaData(SeuratObject, avg_exp_SeuratObject, f=paste(as.character(x=Idents(object=SeuratObject)),SeuratObject$patient,sep='_'))
   return(avg_exp_SeuratObject)
}

PseudobulkSeuratObject_Mean <- function(SeuratObject, CellTypeLevel){
   DefaultAssay(SeuratObject) <- "RNA"
   SeuratObject <- NormalizeData(SeuratObject)
   avg_exp_SeuratObject <- AverageExpression(SeuratObject, group.by = c(CellTypeLevel, "patient"), return.seurat = T)
   #Give back labels to the object
   Idents(SeuratObject) <- CellTypeLevel
   avg_exp_SeuratObject <- AverageMetaData(SeuratObject, avg_exp_SeuratObject, f=paste(as.character(x=Idents(object=SeuratObject)),SeuratObject$patient,sep='_'))
   return(avg_exp_SeuratObject)
}
