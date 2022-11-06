# Pipeline for single cell sequencing workflow 

### Before getting started 
Before getting started, it might be helpful to watch a few videos on single cell sequencing and the goal behind it. I'd recommend these two videos, although basic, 
they provide a good overview of the ultimate goal:
  - [Single-cell sequencing explained in 2 minutes](https://www.youtube.com/watch?v=6UVOdCc1Q7I&ab_channel=Sanbomics)
  - [Introduction to single-cell RNA-Seq and Seurat | Bioinformatics for beginners](https://www.youtube.com/watch?v=xbX49h7BiUU&ab_channel=Bioinformagician)

Overall, the goal of single cell RNA sequencing is to understand what types of cells are present, in what frequency, and how subpopulations of cells differ througout different conditions. In general, I would recommend the two above youtube accounts for help with bioinformatics analysis. [This video](https://www.youtube.com/watch?v=5HBzgsz8qyk&ab_channel=Bioinformagician) is a workflow for single cell sequencing that I'd suggest to watch as a beginner. 

Finally, the package that our lab mainly uses is [Seurat](https://satijalab.org/seurat/). Seurat is a package designed to study single cell RNA sequencing on R. The creators of Seurat have documented many vignettes on their website; these are tutorials that you can follow during your analysis. You should know that there are hundreds of packages designed to analyze single cell single cell sequencing, including some on python. Any of these packages can be used for analysis, based on your comfortability level. 

## Step 1: Getting FASTQ files 
After sequencing is done, the output of sequencing is provided to us in a FASTQ format. Typically, there are two FASTQ reads per sample (Read 1 and Read 2) - you need both! At our lab, these are typically provided on a platform called [Globus Connect Personal](https://www.globus.org/globus-connect-personal). Globus allows for transfer of large data files from computer to computer. [Here](https://docs.globus.org/how-to/) are tutorials for using Globus. You should note that since the data files we are transferring are very large, downloading files from Globus can take hours to days to complete. 

## Step 2: FASTQ to barcoding 
At this point, we have our FASTQ files. FASTQ files contain sequence information, but we are unable to analyze data from this file format. To utlize FASTQ files, they have to be allgned. This can be done manually, but in our lab we use [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger). Cell Ranger is a pipeline that processes single cell data (FASTQ format) and alligns the reads for you. Additionally, it creates a set of feature-barcode matrices. This is used for our downstream analysis. 

How to use Cell Ranger? Cell Ranger can be utilized on linux (terminal) by connecting to a remote server (USC uses [CARC](https://www.carc.usc.edu/)). However, for beginners, the easiest way to use Cell Ranger is through [Cell Ranger Cloud](https://www.10xgenomics.com/products/cloud-analysis). Cell Ranger Cloud is a user friendly UI where you won't have to learn linux syntax. To use Cell Ranger Cloud, you just need to upload the two reads of your sample (the software automatically detects its from the same sample, if the sample is labeled correctly). Then, you can select which type of data yours is (for example, Single Cell 3' Gene Expression) and it will automatically run for you. 

Things to note with Cell Ranger Cloud:
1. Analyses take time -- depending on file size, it can take multiple hours for this pipeline to run 
2. Cell Ranger CLoud is a free service, up to a point. They allow for free data download for all your outputs for one single time. I suggest downloading all your outputs at moment of creation. After a certain number of days, Cell Ranger Cloud will delete all your data. They are not a data storage service. 

What is our desired output of Cell Ranger Cloud? If studying single cell RNA sequencing, the major output that we are looking for are the 'filtered.feature.bc.matrix.h5' files. Download all files, even if you don't think you will be using them for now. 

#### Aside: What are filtered.feature.bc.matrix.h5 files?
'Bc' stands for barcode. These files tell us the features per barcode. Barcodes are cells. Overall, this type of file gives you a table of cells (denoted by barcode) and gene counts (in umi). An .h5 file allows for storage for lots of data. 

## Step 3: Analysis with Seurat 
Now, you can start your analysis! 

For analysis, we use R. If you don't have R downloaded, you can download it [here](https://posit.co/). R is programming language for statistical computing. 
Next, the package we will start off using is [Seurat](https://satijalab.org/seurat/). You will need to install Seurat to your R computing environment -- you can do this through CRAN: 
1. Go to 'Tools'
2. 'Install packages'
3. Type 'Seurat'. It should automatically start downloading. 

We need a few other packages that can be similarily installed:
1. dplyr
2. patchwork

Once you have these packages installed, you need to load them in prior to utilizing it. Note that you need to do this each time you restart R. Just because you have the package downloaded does not mean the environment has it loaded in. The term "library" loads in the specific package. To do this, type:   
`library(Seurat)`   
`library(dplyr)`   
`library(patchwork)`

Once these three packages are loaded in, we can begin. 

### Creating a seurat object
Seurat's package uses 'seurat objects'. These hold the information found in the filtered.feature.bc.matrix.h5 files, and additionally store any analysis data that you perform on the seurat object. In lines of code, "so" stands for seurat object. 

First, we will load in the filtered.feature.bc.matrix.h5 file:   
`file <- Read10X_h5("/path/to/file.filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)`

Then, we can turn it into a seurat object:   
```
so <- CreateSeuratObject(counts = file$`Gene Expression`, project = "projectname")
```

To check if a seurat object was created, we can:   
`so` 
This will show you what type of object it is, the number of features, samples, and assays. 

Additionally, to check the metadata of the seurat object, you can:   
`so[[]]`
Metadata shows you all the information that seurat object holds.  

#### Aside on merging samples
Most likely, you will be analyzing more than one sample. You can introduce all the samples in your dataset similarly to above (1. load file and 2. create seurat object). From here, you can merge all the samples together. To do this:    
`merged_dataset <- merge(so1, y = c(so2, so3, so4), project = "projectname")`

You can merge as many or as less samples as you want. One thing to note, if you are merging samples, when you are creating the individual seurat object at the start, under project you can label the sample. For example, if you are comparing diseased to healthy, some samples would have a `project = "Healthy"`. This identity will be reintroduced later on in the analysis. 

Additionally, please note in all downstream code I am writing "so". The following code can be applied to the "merged_dataset" or seurat object. 

### Preliminary QC
There are three quality control metrics we typically examine when analyzing single cell RNA sequencing data:
1. nCount_RNA: total number of molecules detected in a cell
2. nFeature_RNA: total number of genes detected in a cell
3. Mitochondrial percentage of a cell: indicative of poor sample quality 

Your seurat object will automatically have the first two metrics, however, mitochondrial percentage is something you will need to add. To add this metric to your dataset, type:   
`so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")`

Now, we can visualize all three of our QC metrics in a violin plot:    
`VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)`

From this violin plot, you should be able to visualize these QC metrics. When examining the mitochondrial percentage, if it is too high (typically, "too high" means larger than 20 or 25% -- you can read about this in the literature if you want), this is indicative of low sample quality. High mitochondrial percentage will overshadow any future analysis, so to deal with this, we need to filter out the cells with a high mitochondiral percentage:   
`so <- subset(so, subset = percent.mt <= 25)`
You can adjust for the amount of cut off based on your sample set. Note that this cut off is happening on the same seurat object we were using (i.e. we did not create a new seurat object). You are welcome to create a new seurat object if you prefer. Here, once again, you can visualize on a violin plot to double check that the cut off did indeed happen. 

### Analysis Pipeline
Here is a basic analysis pipeline:   

First, normalizing and scaling the data. This normalized expression measurements and scales the expression of each gene.      
`so <- NormalizeData(so)
all.genes <- rownames(so)
so <- ScaleData(so, features = all.genes)`

Find variable features. From here, you can also examine top variable genes if you'd like, but this step is necessary for downstream analysis.    
`so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)`

Run a PCA. PCA = principal component analysis. PCA is a way to reduce the dimensionality of the dataset.
`so <- RunPCA(so, features = VariableFeatures(object = so))`
You can visualize the PCA:
`DimPlot(so, reduction = "pca")`

Determine dimensionality of dataset. There are two ways to do this: jack straw plot and elbow plot. The jack straw plot is more computationally expensive compared to the elbow plot. If you have a large dataset, the jack straw plot can take upwards of 30 minutes to run. Examining these metrics is useful for determining the optimal number of clusters in your dataset.

Jack straw:   
`so <- JackStraw(so, num.replicate = 100)
so <- ScoreJackStraw(so, dims = 1:20)
JackStrawPlot(so, dims = 1:15)`

Elbow:    
`ElbowPlot(so)`

Find neighbors and clusters. Here, you can optimize the "dims" and "resolution" to get the optimal number of clusters as defined by the previous step. This "optimal" number is not a given, it will likely change once you explore into the contents of the clusters.   
`so <- FindNeighbors(so, dims = 1:10)
so <- FindClusters(so, resolution = 0.5)`

Run UMAP and visualize:
`so <- RunUMAP(so, dims = 1:10)
DimPlot(so)`

### Cluster identification
At this stage, you will have a UMAP with clusters. Now the next question is, what are these clusters? We can do cluster identification in two ways: manually or algorithmically (but if you decide to choose algorithmically, you should always double check manually).

1. Manually: Find the top differentially expressed genes for each cluster, and manually identify the cell type those genes belong to. To do this:
`markers1 <- FindMarkers(so, ident.1 = 1, min.pct = 0.25)
head(markers1, n = 10)`

Here, I am examining the top ten markers (n = 10) for the cluster #1 (ident.1 = 1). Similarily, you can get the top markers for every cluster, and you can adjust the number of markers you want to visualize as well. Note that this step can take time to run, especially if you have a large sample set. 

Now, you have the top differentially expressed genes. I would suggest using the [Human Protein Atlas](https://www.proteinatlas.org/) to identify cell types. You are welcome to use whatever site you prefer for this. 

2. Algorithmically. There are two algorithms that I think work best, but you are welcome to explore more. 
   1. [ScType](https://github.com/IanevskiAleksandr/sc-type)
   I prefer this for scRNA data. For ScType, you will run this after the last "Run UMAP and visualize step". I've added the code for this at the bottom, as it is quite long(*). 
   3. [Azimuth](https://azimuth.hubmapconsortium.org/)
   I prefer this for scATAC data. For Azimuth, you'll need to install and load it in: `library(Azimuth)` and after creating the seurat object, `so <- RunAzimuth(so, reference = "lungref")`. 


### Re-introduction of original identity
Once you know what cluster is what cell type, you can re-introduce the original identiy of the cell. This can mean, bringing back whether that cell is a "diseased" or "healthy" cell. To do this,    
`DimPlot(so, reduction = "umap", group.by = "orig.ident")`
From here, you can make your analyses! 


### Single cell ATAC sequencing
scATAC sequencing is a newer technology, and there aren't a lot of options out there to study it. We use [Signac](https://stuartlab.org/signac/), which is created by the same people who created Seurat. 

You'll need to install and load "Signac". Here is a sample workflow for scATAC analysis (under the impression that these are multimodal samples): 

```
inputdata.object1 <- Read10X_h5("/path/to/file.filtered_feature_bc_matrix.h5")
object1 <- CreateSeuratObject(counts = inputdata.object1$`Gene Expression`, project = "projectName")
object1 <- RunAzimuth(object1, reference = "lungref")
object1[["ATAC"]] <- CreateAssayObject(counts = inputdata.object1$Peaks)
object1[["percent.mt"]] <- PercentageFeatureSet(object1, pattern = "MT")
object1 <- subset(x = object1, percent.mt < 20 )
```

Here, we loaded in the file, created an object, brought in the "peaks" information (ATAC) and filtered. We also ran Azimuth for downstream analysis. 

You can do this for how many ever samples you have, and merge: `mergedobject <- merge(object1, y = c(object2, object3), project = "projectName")`

Then, 
```
DefaultAssay(mergedobject) <- "RNA"
mergedobject <- SCTransform(mergedobject, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

DefaultAssay(mergedobject) <- "ATAC"
mergedobject <- RunTFIDF(mergedobject)
mergedobject <- FindTopFeatures(mergedobject, min.cutoff = 'q0')
mergedobject <- RunSVD(mergedobject)

mergedobject <- RunUMAP(mergedobject, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

mergedobject <- FindMultiModalNeighbors(mergedobject, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
mergedobject <- RunUMAP(mergedobject, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
mergedobject <- FindClusters(mergedobject, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(mergedobject, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(mergedobject, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(mergedobject, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

p4 <- DimPlot(mergedobject, reduction = "umap.rna", group.by = "predicted.ann_level_4", label = TRUE, label.size = 3)
p5 <- DimPlot(mergedobject, reduction = "umap.atac", group.by = "predicted.ann_level_4", label = TRUE, label.size = 3)
p6 <- DimPlot(mergedobject, reduction = "wnn.umap", group.by = "predicted.ann_level_4", label = TRUE, label.size = 3)

p4 + p5 + p6
```
Note that this can take upwards of an hour to run, based on the number of samples. 


### Random tips
1. If you get stuck on anything, google it! You can try to find vignettes, youtube tutorials, or look on stack overflow for help with code.
2. Remember to save your code! Especially with a lot of samples, the program can crash at times. 
3. Save your images in a powerpoint. 
4. If you run out of memory, you can remove old objects you don't need by typing `rm(object)`. This opens up memory for your computer to run. 

(*) Code for ScType:    
```
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Lung" 

gs_list = gene_sets_prepare(db_, tissue)

es.max = sctype_score(scRNAseqData = **so**[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

cL_resutls = do.call("rbind", lapply(unique(**so**@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(**so**@meta.data[**so**@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(**so**@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores) 

sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

**so**@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  **so**@meta.data$customclassif[**so**@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(**so**, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')  
```
