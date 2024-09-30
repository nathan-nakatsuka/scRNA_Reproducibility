
library(Seurat)
library(dplyr)
library(Azimuth)
library(SeuratData)
library(presto)
library(car)


##### Step 1: Obtain cell types by mapping to Azimuth reference
# Download Azimuth human motor cortex reference from here: https://zenodo.org/records/4546932

# Load Azimuth human motor cortex
CortexRef = readRDS("ref.Rds")

### Table of all Individuals and Metadata
FileNames = read.table("IDs.txt",header=F)
Individuals = unique(FileNames[,1])
IndividualTable = data.frame(1:length(Individuals))
IndividualTable[,1]=Individuals
Metadata = read.table("MetaData.txt",header=T)

## Goal: create a merged Seurat object with all individuals.
# Read in Seurat files for each individual and name them IndividualID_Seurat.
### Merge all data and then map them to Azimuth reference (integration beforehand is not necessary if mapping to the same reference)
objs <- list()
for(i in 1:nrow(IndividualTable)){
  objs[i] <- get(paste(IndividualTable[i,1],"_Seurat",sep=""))
}
Total_Seurat <- merge(objs[[1]],objs[2:nrow(IndividualTable)])

####QC (example; note: metrics are chosen based on the individual datasets)
Total_Seurat[["percent.mt"]] <- PercentageFeatureSet(Total_Seurat, pattern = "^MT-",assay="RNA")
Total_Seurat <- subset(Total_Seurat, subset = percent.mt <20 & nCount_RNA<65000 & nCount_RNA>200 & nFeature_RNA<9000 & nFeature_RNA>200)

Total_Seurat <- SCTransform(Total_Seurat, vst.flavor = "v2")
Total_Seurat <- RunPCA(Total_Seurat, npcs = 30)
Total_Seurat <- RunUMAP(Total_Seurat, reduction = "pca", dims = 1:30) 

## Find Anchors
Total_Seurat.anchors <- FindTransferAnchors(reference = CortexRef,
                                             query = Total_Seurat, 
                                             dims = 1:30, 
                                             reference.reduction = "refDR", 
                                             features = rownames(CortexRef[['refDR']]@feature.loadings),
                                             k.filter = NA)
#Map to Azimuth cortex reference
Total_Seurat <- MapQuery(
  anchorset = Total_Seurat.anchors,
  query = Total_Seurat,
  reference = CortexRef,
  refdata = list(
    celltype.l2 = "subclass",
    celltype.l1 = "class"
  ), 
  reduction.model = "refUMAP"
)

Total_Seurat$predicted.celltype.l3 <- Total_Seurat$predicted.celltype.l2
Total_Seurat$predicted.celltype.l3 <- car::recode(Total_Seurat$predicted.celltype.l3, "c('Astro') = c('Astrocytes')")
Total_Seurat$predicted.celltype.l3 <- car::recode(Total_Seurat$predicted.celltype.l3, "c('Endo') = c('Endothelial_Cell')")
Total_Seurat$predicted.celltype.l3 <- car::recode(Total_Seurat$predicted.celltype.l3, "c('Oligo') = c('Oligodendrocytes')")
Total_Seurat$predicted.celltype.l3 <- car::recode(Total_Seurat$predicted.celltype.l3, "c('L2/3 IT','L5 IT','L5 ET','L5/6 NP','L6 CT','L6 IT','L6 IT Car3','L6b') = c('Glut_ExcitatoryNeuron')")
Total_Seurat$predicted.celltype.l3 <- car::recode(Total_Seurat$predicted.celltype.l3, "c('Vip','Sst','Sst Chodl','Sncg','Pvalb','Lamp5') = c('GABA_InhibitoryNeuron')")
Total_Seurat$predicted.celltype.l3 <- car::recode(Total_Seurat$predicted.celltype.l3, "c('Micro-PVM') = c('Microglia')")
Total_Seurat$predicted.celltype.l3 <- car::recode(Total_Seurat$predicted.celltype.l3, "c('OPC') = c('Oligodendrocyte_Precursor_Cell')")

Metadata = read.table("MetaData.txt",header=T)
row.names(Metadata)=Metadata$Individual
Metadata2 <- Metadata
Metadata2$Individual <- NULL
md <- Metadata[as.character(x = Total_Seurat$Individual), ]
row.names(x = md) <- colnames(Total_Seurat)

Total_Seurat <- AddMetaData(
  object = Total_Seurat,
  metadata = md
)
Total_Seurat$SEX = Total_Seurat$Sex

saveRDS(Total_Seurat,file="Total_Seurat_filteredmapped_annotated.rds")
