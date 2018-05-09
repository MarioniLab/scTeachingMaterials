---
title: Preparing to correct batch effects in single-cell RNA-seq data
author: Aaron T. L. Lun
output: 
  BiocStyle::html_document:
    fit_caption: false
---



# Setting up the Nestorowa data

## Data input

We download the counts from NCBI GEO corresponding to https://doi.org/10.1182/blood-2016-05-716480.


```r
fnameN <- "GSE81682_HTSeq_counts.txt.gz"
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81682&format=file&file=GSE81682%5FHTSeq%5Fcounts%2Etxt%2Egz", 
    fnameN)
```

We load in the counts:


```r
dataN <- read.table(fnameN, header=TRUE, row.names=1, check.names=FALSE)
dataN <- as.matrix(dataN)
```

We convert this to a sparse matrix to save space:


```r
library(Matrix)
dataN <- as(dataN, "dgCMatrix")
dim(dataN)
```

```
## [1] 46175  1920
```

... and we get rid of non-genes at the end:



```r
dataN <- dataN[!grepl("^_", rownames(dataN)),]
dim(dataN)
```

```
## [1] 46170  1920
```

We construct a `SingleCellExperiment` object:


```r
library(SingleCellExperiment)
sceN <- SingleCellExperiment(list(counts=dataN))
```

... and throw in some spike-in information:


```r
isSpike(sceN, "ERCC") <- grepl("^ERCC", rownames(sceN))
```

We also add information about the cell type:


```r
sceN$CellType <- sub("_.*", "", colnames(sceN))
```

## Normalizing for cell-specific biases

The authors have already done the quality control for us, so we skip straight to the normalization.
First, we break up the dataset with `quickCluster`:


```r
library(scran)
clustersN <- quickCluster(sceN, method="igraph")
table(clustersN)
```

```
## clustersN
##   1   2   3 
## 353 792 775
```

Then we apply the deconvolution method:


```r
sceN <- computeSumFactors(sceN, clusters=clustersN)
summary(sizeFactors(sceN))
```

```
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##  0.008336  0.364923  0.738963  1.000000  1.305771 15.033123
```

We plot the size factors against the library sizes as a sanity check:


```r
plot(Matrix::colSums(counts(sceN)), sizeFactors(sceN), log="xy")
```

<img src="prepare_files/figure-html/unnamed-chunk-10-1.png" width="100%" />

We also generate size factors for the spike-ins separately.
Remember, it's fine that the spike-in size factors are not well-correlated to the deconvolution size factors.


```r
sceN <- computeSpikeFactors(sceN, general.use=FALSE, type="ERCC")
plot(sizeFactors(sceN), sizeFactors(sceN, "ERCC"), log="xy")
```

<img src="prepare_files/figure-html/unnamed-chunk-11-1.png" width="100%" />

Finally, we compute log-normalized expression values:


```r
sceN <- normalize(sceN)
```

## Detecting highly variable genes

Let's pick out some highly variable genes.


```r
fitN <- trendVar(sceN, parametric=TRUE, loess.args=list(span=0.2))
decN <- decomposeVar(sceN, fitN)
head(decN)
```

```
## DataFrame with 6 rows and 6 columns
##                                   mean              total
##                              <numeric>          <numeric>
## ENSMUSG00000000001    6.34334751550413   7.93088831876345
## ENSMUSG00000000003 0.00979007329151634 0.0253290312805928
## ENSMUSG00000000028     3.9902456058454   9.35478320944155
## ENSMUSG00000000031   0.111949614544398  0.460678594211501
## ENSMUSG00000000037    1.08606588143884   4.13900720759496
## ENSMUSG00000000049  0.0214270738805416 0.0572629876677554
##                                    bio               tech
##                              <numeric>          <numeric>
## ENSMUSG00000000001      1.211919636143   6.71896868262045
## ENSMUSG00000000003 -0.0268192853388259 0.0521483166194187
## ENSMUSG00000000028  -0.577730972579486   9.93251418202104
## ENSMUSG00000000031  -0.132253854687009  0.592932448898509
## ENSMUSG00000000037  -0.815115814425082   4.95412302202004
## ENSMUSG00000000049  -0.056871585117566  0.114134572785321
##                                 p.value                  FDR
##                               <numeric>            <numeric>
## ENSMUSG00000000001 6.00580100025002e-08 4.75490203590242e-06
## ENSMUSG00000000003                    1                    1
## ENSMUSG00000000028    0.966156176729986                    1
## ENSMUSG00000000031    0.999999999999966                    1
## ENSMUSG00000000037    0.999999965828181                    1
## ENSMUSG00000000049                    1                    1
```

Having a look at them on the plot:


```r
plot(decN$mean, decN$total, xlab="Mean log-expression",
    ylab="Variance of log-expression")
points(fitN$mean, fitN$var, col="red", pch=16)
curve(fitN$trend(x), col="red", add=TRUE)
```

<img src="prepare_files/figure-html/unnamed-chunk-14-1.png" width="100%" />


```
##            used  (Mb) gc trigger   (Mb)  max used   (Mb)
## Ncells  5119599 273.5    9630424  514.4   9630424  514.4
## Vcells 77674854 592.7  211447095 1613.3 214450664 1636.2
```

# Setting up the Paul data

## Data input

We download the counts from NCBI GEO corresponding to https://doi.org/10.1016/j.cell.2015.11.013.


```r
fnameP <- "umitab_Amit.txt.gz"
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE72857&format=file&file=GSE72857%5Fumitab%2Etxt%2Egz", 
    fnameP)
```

We load in the counts:


```r
dataP <- read.table(fnameP, header=TRUE, row.names=1, check.names=FALSE)
dataP <- as.matrix(dataP)
```

We convert this to a sparse matrix, and take a random set of 2000 cells for the sake of speed:


```r
set.seed(0)
chosen <- sample(ncol(dataP), 2000)
dataP <- as(dataP[,chosen], "dgCMatrix")
dim(dataP)
```

```
## [1] 27297  2000
```

We convert the gene names back to Ensembl:


```r
library(org.Mm.eg.db)
symbols <- sub(".*;", "", rownames(dataP))
ens <- mapIds(org.Mm.eg.db, keys=symbols, keytype="SYMBOL", column="ENSEMBL")
```

We remove duplicates or genes without any Ensembl:


```r
keep <- !is.na(ens) & !duplicated(ens)
summary(keep)
```

```
##    Mode   FALSE    TRUE 
## logical    9205   18092
```

```r
dataP <- dataP[keep,]
rownames(dataP) <- unname(ens[keep])
```

We now construct a `SingleCellExperiment` object:


```r
sceP <- SingleCellExperiment(list(counts=dataP))
```

## Removing low-quality cells

We will remove some low quality cells here:


```r
library(scater)
sceP <- calculateQCMetrics(sceP)
multiplot(cols=2,
    plotColData(sceP, "log10_total_counts"),
    plotColData(sceP, "log10_total_features_by_counts"))
```

<img src="prepare_files/figure-html/unnamed-chunk-22-1.png" width="100%" />

We'll apply our outlier strategy:


```r
low.total <- isOutlier(sceP$log10_total_counts, type="lower", nmads=3)
low.nexpr <- isOutlier(sceP$log10_total_features_by_counts, type="lower", nmads=3)
discard <- low.total | low.nexpr
data.frame(LowLib=sum(low.total), LowNexprs=sum(low.nexpr), Lost=sum(discard))
```

```
##   LowLib LowNexprs Lost
## 1     98       138  138
```

... and discard the offending cells:


```r
sceP <- sceP[,!discard]
```

## Normalizing for cell-specific biases

Again, we break up the dataset with `quickCluster` - note `min.mean` is set lower than the default, as we are dealing with UMI counts.


```r
library(scran)
clustersP <- quickCluster(sceP, method="igraph", min.mean=0.1)
table(clustersP)
```

```
## clustersP
##   1   2 
## 885 977
```

Then we apply the deconvolution method:


```r
sceP <- computeSumFactors(sceP, clusters=clustersP, min.mean=0.1)
summary(sizeFactors(sceP))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.02028 0.36357 0.78597 1.00000 1.37441 8.10157
```

We plot the size factors against the library sizes as a sanity check:


```r
plot(Matrix::colSums(counts(sceP)), sizeFactors(sceP), log="xy")
```

<img src="prepare_files/figure-html/unnamed-chunk-27-1.png" width="100%" />

Finally, we compute log-normalized expression values:


```r
sceP <- normalize(sceP)
```

## Detecting highly variable genes

Let's pick out some highly variable genes. 
There aren't any spike-ins, so we'll have to use `use.spikes=FALSE`.


```r
fitP <- trendVar(sceP, use.spikes=FALSE, loess.args=list(span=0.01))
decP <- decomposeVar(sceP, fitP)
head(decP)
```

```
## DataFrame with 6 rows and 6 columns
##                                  mean              total
##                             <numeric>          <numeric>
## ENSMUSG00000007777 0.0371305230862186 0.0400754051582754
## ENSMUSG00000024442  0.150766095383802  0.174716159359473
## ENSMUSG00000078880                  0                  0
## ENSMUSG00000093989  0.956864963377254  0.918615249367187
## ENSMUSG00000107002   0.48964365559762  0.509101611842776
## ENSMUSG00000058706 0.0644151149273843 0.0768957720586398
##                                     bio               tech           p.value
##                               <numeric>          <numeric>         <numeric>
## ENSMUSG00000007777 -0.00983619395322893 0.0499115991115043 0.999999999943339
## ENSMUSG00000024442  -0.0121316752377102  0.186847834597183 0.978003138518491
## ENSMUSG00000078880                    0                  0                 1
## ENSMUSG00000093989  -0.0115400651780048  0.930155314545192 0.643969515694305
## ENSMUSG00000107002  -0.0420411886549673  0.551142800497743   0.9912918542818
## ENSMUSG00000058706 -0.00969232651222003 0.0865880985708598 0.999799862578407
##                          FDR
##                    <numeric>
## ENSMUSG00000007777         1
## ENSMUSG00000024442         1
## ENSMUSG00000078880         1
## ENSMUSG00000093989         1
## ENSMUSG00000107002         1
## ENSMUSG00000058706         1
```

Having a look at them on the plot:


```r
plot(decP$mean, decP$total, xlab="Mean log-expression",
    ylab="Variance of log-expression")
curve(fitP$trend(x), col="blue", add=TRUE)
```

<img src="prepare_files/figure-html/unnamed-chunk-30-1.png" width="100%" />

# Cleaning up

We clean up after ourselves by deleting the two source files.


```r
unlink(fnameN)
unlink(fnameP)
```

We also save the objects for further use.


```r
save(file="objects.Rda", sceN, sceP, decN, decP, fitN, fitP)
```

Finally, as this procedure was somewhat involved, we will report the session information as well.


```r
sessionInfo()
```

```
## R version 3.5.0 Patched (2018-04-30 r74681)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 16.04.4 LTS
## 
## Matrix products: default
## BLAS: /home/cri.camres.org/lun01/Software/R/R-3-5-branch-release/lib/libRblas.so
## LAPACK: /home/cri.camres.org/lun01/Software/R/R-3-5-branch-release/lib/libRlapack.so
## 
## locale:
##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] scater_1.8.0                ggplot2_2.2.1              
##  [3] org.Mm.eg.db_3.6.0          AnnotationDbi_1.42.1       
##  [5] scran_1.8.1                 SingleCellExperiment_1.2.0 
##  [7] SummarizedExperiment_1.10.0 DelayedArray_0.6.0         
##  [9] BiocParallel_1.14.1         matrixStats_0.53.1         
## [11] Biobase_2.40.0              GenomicRanges_1.32.2       
## [13] GenomeInfoDb_1.16.0         IRanges_2.14.5             
## [15] S4Vectors_0.18.1            BiocGenerics_0.26.0        
## [17] Matrix_1.2-14               knitr_1.20                 
## [19] BiocStyle_2.8.0            
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-6             bit64_0.9-7             
##  [3] rprojroot_1.3-2          dynamicTreeCut_1.63-1   
##  [5] tools_3.5.0              backports_1.1.2         
##  [7] R6_2.2.2                 irlba_2.3.2             
##  [9] DT_0.4                   vipor_0.4.5             
## [11] DBI_1.0.0                lazyeval_0.2.1          
## [13] colorspace_1.3-2         gridExtra_2.3           
## [15] bit_1.1-12               compiler_3.5.0          
## [17] labeling_0.3             bookdown_0.7            
## [19] scales_0.5.0             stringr_1.3.0           
## [21] digest_0.6.15            rmarkdown_1.9           
## [23] XVector_0.20.0           pkgconfig_2.0.1         
## [25] htmltools_0.3.6          limma_3.36.1            
## [27] htmlwidgets_1.2          rlang_0.2.0             
## [29] RSQLite_2.1.1            FNN_1.1                 
## [31] shiny_1.0.5              DelayedMatrixStats_1.2.0
## [33] bindr_0.1.1              dplyr_0.7.4             
## [35] RCurl_1.95-4.10          magrittr_1.5            
## [37] GenomeInfoDbData_1.1.0   Rcpp_0.12.16            
## [39] ggbeeswarm_0.6.0         munsell_0.4.3           
## [41] Rhdf5lib_1.2.0           viridis_0.5.1           
## [43] stringi_1.2.2            yaml_2.1.19             
## [45] edgeR_3.22.1             zlibbioc_1.26.0         
## [47] rhdf5_2.24.0             plyr_1.8.4              
## [49] blob_1.1.1               grid_3.5.0              
## [51] promises_1.0.1           shinydashboard_0.7.0    
## [53] lattice_0.20-35          cowplot_0.9.2           
## [55] locfit_1.5-9.1           pillar_1.2.2            
## [57] igraph_1.2.1             rjson_0.2.18            
## [59] reshape2_1.4.3           glue_1.2.0              
## [61] evaluate_0.10.1          data.table_1.11.2       
## [63] httpuv_1.4.2             gtable_0.2.0            
## [65] assertthat_0.2.0         xfun_0.1                
## [67] mime_0.5                 xtable_1.8-2            
## [69] later_0.7.2              viridisLite_0.3.0       
## [71] tibble_1.4.2             memoise_1.1.0           
## [73] beeswarm_0.2.3           tximport_1.8.0          
## [75] bindrcpp_0.2.2           statmod_1.4.30
```
