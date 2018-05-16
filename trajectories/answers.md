---
title: Creating trajectories from scRNA-seq data
author: Aaron Lun
date: "2018-05-16"
output:
  BiocStyle::html_document:
    fig_caption: false
---



# Introduction

For continuous phenomena, it is less appropriate to summarize the data in terms of discrete clusters.
Rather, we want to characterize the cell "trajectories", i.;e., curves in the high-dimensional expression space on which the cells lie.
These trajectories (hopefully) correspond to continuous biological processes such as differentiation.

To demonstrate, we will be using the Nestorowa dataset from the previous workshop. 
This involves capturing cells as they progress throughout haematopoietic differentiation, and is a good candidate for playing with trajectories:


```r
load("../unbatch/objects.Rda")
sceN
```

```
## class: SingleCellExperiment 
## dim: 46170 1920 
## metadata(1): log.exprs.offset
## assays(2): counts logcounts
## rownames(46170): ENSMUSG00000000001 ENSMUSG00000000003 ...
##   ERCC-00170 ERCC-00171
## rowData names(0):
## colnames(1920): HSPC_007 HSPC_013 ... Prog_852 Prog_810
## colData names(1): CellType
## reducedDimNames(0):
## spikeNames(1): ERCC
```

We'll perform a PCA, and do our trajectory calculations on the PCs instead of the original expression values.


```r
library(scran)
sceN <- denoisePCA(sceN, decN, approximate=TRUE)
rd <- reducedDim(sceN, "PCA")
ncol(rd)
```

```
## [1] 5
```

# Using *[destiny](http://bioconductor.org/packages/destiny)*

## Creating visualizations

First we will use the *[destiny](http://bioconductor.org/packages/destiny)* package, which models trajectories through diffusion maps.

1. We calculate distances between all pairs of cells (or in practice, just the nearest neighbours).
2. The probability of moving from one cell to another is dependent on the distance (larger distance = lower probability).
3. We build a transition probability matrix, normalized by the density of cells.
4. We define diffusion components from the eigenvectors of the transition matrix.


```r
library(destiny)
dm <- DiffusionMap(rd) # cells are rows here!
dm
```

```
## DiffusionMap (20 Diffusion components and 1920 observations)
## eigenvalues:    num [1:20] 0.996 0.983 0.978 0.975 0.968 ...
## eigenvectors:   num [1:1920, 1:20] -0.00258 -0.00908 -0.00778 -0.01031 -0.01189 ...
##   ..colnames:   chr [1:20] "DC1" "DC2" "DC3" "DC4" ...
## optimal_sigma:  num [1:1920] 7.05 6.53 7.23 5.69 7.6 ...
## distance:       chr "euclidean"
```

We extract the first two "diffusion components":


```r
plot(dm@eigenvectors[,1], dm@eigenvectors[,2], xlab="DC1", ylab="DC2")
```

<img src="answers_files/figure-html/unnamed-chunk-4-1.png" width="100%" />

We can store this in the object for more general plotting with *[scater](http://bioconductor.org/packages/scater)*:


```r
library(scater)
reducedDim(sceN, "DiffusionMap") <- dm@eigenvectors
plotDiffusionMap(sceN, colour_by="ENSMUSG00000031162") # Gata1
```

<img src="answers_files/figure-html/unnamed-chunk-5-1.png" width="100%" />

More plotting:


```r
plotDiffusionMap(sceN, colour_by="ENSMUSG00000015355") # Cd48
```

<img src="answers_files/figure-html/unnamed-chunk-6-1.png" width="100%" />

More plotting:


```r
plotDiffusionMap(sceN, colour_by="ENSMUSG00000006389") # Mpl
```

<img src="answers_files/figure-html/unnamed-chunk-7-1.png" width="100%" />

<div class="alert alert-warning">
**Exercise:** 

Direct slot access is **Bad**:

- Internal fields are generally subject to changes without prior notice.
- Data representation and scientific meaning may not be equivalent.

Avoid if possible!
</div>

## Getting pseudotime orderings

A visualization is nice, but to perform analyses we need the location of cells along the trajectory.
This is referred to as "pseudotime", and can be obtained using the `DPT` function.


```r
dp.out <- DPT(dm)
```

We can see the branch assignments from pulling out the `branch` slot.


```r
sceN$branch <- factor(dp.out@branch[,1])
plotDiffusionMap(sceN, colour_by="branch")
```

<img src="answers_files/figure-html/unnamed-chunk-9-1.png" width="100%" />

We can also get pseudotime values:


```r
sceN$pseudotime <- dp.out[,1]
plotDiffusionMap(sceN, colour_by="pseudotime")
```

<img src="answers_files/figure-html/unnamed-chunk-10-1.png" width="100%" />

# Testing for DE along pseudotime

We can use the pseudotime orderings to test for DE along pseudotime, using methods from the *[limma](http://bioconductor.org/packages/limma)* package:


```r
library(limma)
design <- model.matrix(~sceN$pseudotime)
fit <- lmFit(logcounts(sceN), design)
fit <- eBayes(fit, robust=TRUE, trend=TRUE)
summary(decideTests(fit))    
```

```
##        (Intercept) sceN$pseudotime
## Down            87            9263
## NotSig       10323           25946
## Up           35760           10961
```

```r
topTable(fit)   
```

```
##                         logFC   AveExpr         t       P.Value
## ENSMUSG00000049775 -0.4244947 9.3474010 -54.49938  0.000000e+00
## ENSMUSG00000009687 -0.4943123 8.3126575 -53.96799  0.000000e+00
## ENSMUSG00000030707 -0.4887822 8.6538472 -47.38497  0.000000e+00
## ENSMUSG00000030220 -0.3484593 9.7785415 -47.34377  0.000000e+00
## ENSMUSG00000086438  0.3023960 0.8189963  45.70263 3.513250e-309
## ENSMUSG00000041237  0.4319794 1.6029404  45.40354 2.453450e-306
## ENSMUSG00000028581 -0.3784900 8.8925078 -45.15360 5.846834e-304
## ENSMUSG00000008843  0.5042432 2.7213878  41.14017 8.780261e-266
## ENSMUSG00000006567  0.4359989 2.1404924  40.79271 1.749830e-262
## ENSMUSG00000001249  0.4574981 2.1317109  40.66732 2.711703e-261
##                        adj.P.Val        B
## ENSMUSG00000049775  0.000000e+00 887.9938
## ENSMUSG00000009687  0.000000e+00 876.6112
## ENSMUSG00000030707  0.000000e+00 733.9716
## ENSMUSG00000030220  0.000000e+00 733.0716
## ENSMUSG00000086438 3.244135e-305 697.1794
## ENSMUSG00000041237 1.887929e-302 690.6309
## ENSMUSG00000028581 3.856404e-300 685.1575
## ENSMUSG00000008843 5.067308e-262 597.2603
## ENSMUSG00000006567 8.976628e-259 589.6640
## ENSMUSG00000001249 1.251993e-257 586.9238
```

What about changes along one branch?


```r
on.branch.3 <- which(sceN$branch==3) # NA protection
sceN3 <- sceN[,on.branch.3]
design3 <- model.matrix(~sceN3$pseudotime)
fit3 <- lmFit(logcounts(sceN3), design3)
fit3 <- eBayes(fit3, robust=TRUE, trend=TRUE)
summary(decideTests(fit3))    
```

```
##        (Intercept) sceN3$pseudotime
## Down            61              323
## NotSig       34996            45036
## Up           11113              811
```

```r
topTable(fit3)   
```

```
##                         logFC   AveExpr         t      P.Value    adj.P.Val
## ENSMUSG00000022584  4.8078510  6.545605  14.13867 3.836339e-33 1.771238e-28
## ENSMUSG00000057729  1.7553070 11.120422  13.62944 1.856351e-31 4.285385e-27
## ENSMUSG00000031722  3.8414705  5.085402  12.71399 1.904671e-28 2.931288e-24
## ENSMUSG00000009350  1.7700128 13.366380  12.27534 5.158000e-27 5.953622e-23
## ENSMUSG00000044258 -2.7150372  4.871581 -12.08714 2.098604e-26 1.937851e-22
## ENSMUSG00000086567 -0.7718498  6.616103 -11.86062 1.134018e-25 8.726267e-22
## ENSMUSG00000020125  4.2823108  9.135763  11.84033 1.325954e-25 8.745615e-22
## ENSMUSG00000079018  2.5190199  2.689835  11.74306 2.714192e-25 1.566428e-21
## ENSMUSG00000025130  0.8628608 10.115648  11.46723 2.085667e-24 1.069947e-20
## ENSMUSG00000029322  0.8034596 11.749411  11.33138 5.667987e-24 2.616909e-20
##                           B
## ENSMUSG00000022584 63.48868
## ENSMUSG00000057729 59.61379
## ENSMUSG00000031722 52.69192
## ENSMUSG00000009350 49.40069
## ENSMUSG00000044258 48.00090
## ENSMUSG00000086567 46.31855
## ENSMUSG00000020125 46.16281
## ENSMUSG00000079018 45.44843
## ENSMUSG00000025130 43.41583
## ENSMUSG00000029322 42.41957
```

<div class="alert alert-warning">
**Exercise:** 

What about non-linear changes over time?


```r
spl <- splines::ns(sceN$pseudotime, df=4)
design.spl <- model.matrix(~spl)
fit <- lmFit(logcounts(sceN), design.spl)
fit <- eBayes(fit, robust=TRUE, trend=TRUE)
topTable(fit)   
```

```
##                         spl1       spl2       spl3      spl4  AveExpr
## ENSMUSG00000001249 -4.599053 -1.0162186 -2.8248708  6.558791 2.131711
## ENSMUSG00000006567 -3.334459  1.2718687  0.3696295  5.872700 2.140492
## ENSMUSG00000008843 -6.402860 -1.2870068 -4.6580211  6.448468 2.721388
## ENSMUSG00000009687  2.432390 -4.4706492 -3.8065762 -6.124916 8.312658
## ENSMUSG00000023926 -6.133549 -0.3172904 -3.4821534  6.130932 2.804653
## ENSMUSG00000028132 -3.669380 -2.2724342 -2.6382883  6.829254 1.832446
## ENSMUSG00000028581  2.585791 -0.4656129 -1.1513364 -5.562100 8.892508
## ENSMUSG00000028825 -6.670475  0.6288379 -3.9719149  5.542931 3.260629
## ENSMUSG00000030000 -4.471329 -0.2366134 -2.1926073  5.822867 1.991842
## ENSMUSG00000030220  2.054157 -0.1522402 -0.9984370 -5.466239 9.778542
##                            F P.Value adj.P.Val
## ENSMUSG00000001249  631.0278       0         0
## ENSMUSG00000006567  571.8284       0         0
## ENSMUSG00000008843  787.4870       0         0
## ENSMUSG00000009687 1005.5271       0         0
## ENSMUSG00000023926  730.6586       0         0
## ENSMUSG00000028132  581.6446       0         0
## ENSMUSG00000028581  782.4958       0         0
## ENSMUSG00000028825  707.8036       0         0
## ENSMUSG00000030000  624.9676       0         0
## ENSMUSG00000030220  798.1835       0         0
```
</div>

# Using *[monocle](http://bioconductor.org/packages/monocle)*

Next, we will use the *[destiny](http://bioconductor.org/packages/destiny)* package.
This _used_ to use a minimum spanning tree on the cell expression profiles, now it uses reverse graph embedding.
The idea is to cluster cells in low-dimensional space and build a MST on the cluster centroids, while preserving similarities in the high-dimensional space.


```r
library(monocle)
cds <- convertTo(sceN, type="monocle")
cds <- setOrderingFilter(cds, rownames(decN)[which(decN$FDR<=0.05)])
```

We do dimensionality reduction for visualization:


```r
cds <- reduceDimension(cds, max_components=2, method='DDRTree')
```

And we define the cell orderings.
This gives us the `Pseudotime` and the `State` variables.


```r
cds <- orderCells(cds)    
head(phenoData(cds)$Pseudotime)
```

```
## [1] 19.01258 14.70975 18.72591 15.39094 26.47259 20.59203
```

```r
head(phenoData(cds)$State)
```

```
## [1] 4 3 4 3 7 5
## Levels: 1 2 3 4 5 6 7
```

We can plot this using the `plot` method:


```r
plot_cell_trajectory(cds)
```

<img src="answers_files/figure-html/unnamed-chunk-17-1.png" width="100%" />

... or the pseudotimes.


```r
plot_cell_trajectory(cds, color_by="Pseudotime")
```

<img src="answers_files/figure-html/unnamed-chunk-18-1.png" width="100%" />

Again, this information can be used in DE testing with *[limma](http://bioconductor.org/packages/limma)* or other packages.

# Session information


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
##  [1] splines   parallel  stats4    stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] bindrcpp_0.2.2              monocle_2.8.0              
##  [3] DDRTree_0.1.5               irlba_2.3.2                
##  [5] VGAM_1.0-5                  Matrix_1.2-14              
##  [7] limma_3.36.1                scater_1.8.0               
##  [9] ggplot2_2.2.1               destiny_2.10.0             
## [11] scran_1.8.1                 SingleCellExperiment_1.2.0 
## [13] SummarizedExperiment_1.10.1 DelayedArray_0.6.0         
## [15] matrixStats_0.53.1          Biobase_2.40.0             
## [17] GenomicRanges_1.32.2        GenomeInfoDb_1.16.0        
## [19] IRanges_2.14.9              S4Vectors_0.18.1           
## [21] BiocGenerics_0.26.0         BiocParallel_1.14.1        
## [23] knitr_1.20                  BiocStyle_2.8.0            
## 
## loaded via a namespace (and not attached):
##   [1] readxl_1.1.0             backports_1.1.2         
##   [3] RcppEigen_0.3.3.4.0      plyr_1.8.4              
##   [5] igraph_1.2.1             lazyeval_0.2.1          
##   [7] sp_1.2-7                 shinydashboard_0.7.0    
##   [9] densityClust_0.3         fastICA_1.2-1           
##  [11] digest_0.6.15            htmltools_0.3.6         
##  [13] viridis_0.5.1            magrittr_1.5            
##  [15] cluster_2.0.7-1          openxlsx_4.0.17         
##  [17] docopt_0.4.5             xts_0.10-2              
##  [19] colorspace_1.3-2         ggrepel_0.8.0           
##  [21] haven_1.1.1              xfun_0.1                
##  [23] dplyr_0.7.4              sparsesvd_0.1-4         
##  [25] RCurl_1.95-4.10          tximport_1.8.0          
##  [27] bindr_0.1.1              zoo_1.8-1               
##  [29] glue_1.2.0               gtable_0.2.0            
##  [31] zlibbioc_1.26.0          XVector_0.20.0          
##  [33] car_3.0-0                Rhdf5lib_1.2.0          
##  [35] DEoptimR_1.0-8           abind_1.4-5             
##  [37] VIM_4.7.0                scales_0.5.0            
##  [39] pheatmap_1.0.8           edgeR_3.22.1            
##  [41] ggthemes_3.5.0           Rcpp_0.12.16            
##  [43] viridisLite_0.3.0        xtable_1.8-2            
##  [45] laeken_0.4.6             foreign_0.8-70          
##  [47] proxy_0.4-22             DT_0.4                  
##  [49] vcd_1.4-4                htmlwidgets_1.2         
##  [51] FNN_1.1                  RColorBrewer_1.1-2      
##  [53] pkgconfig_2.0.1          nnet_7.3-12             
##  [55] locfit_1.5-9.1           dynamicTreeCut_1.63-1   
##  [57] labeling_0.3             rlang_0.2.0             
##  [59] reshape2_1.4.3           later_0.7.2             
##  [61] munsell_0.4.3            cellranger_1.1.0        
##  [63] tools_3.5.0              evaluate_0.10.1         
##  [65] stringr_1.3.1            yaml_2.1.19             
##  [67] robustbase_0.93-0        RANN_2.5.1              
##  [69] mime_0.5                 slam_0.1-43             
##  [71] compiler_3.5.0           beeswarm_0.2.3          
##  [73] curl_3.2                 e1071_1.6-8             
##  [75] smoother_1.1             tibble_1.4.2            
##  [77] statmod_1.4.30           stringi_1.2.2           
##  [79] forcats_0.3.0            lattice_0.20-35         
##  [81] HSMMSingleCell_0.114.0   pillar_1.2.2            
##  [83] combinat_0.0-8           lmtest_0.9-36           
##  [85] data.table_1.11.2        cowplot_0.9.2           
##  [87] bitops_1.0-6             httpuv_1.4.3            
##  [89] R6_2.2.2                 bookdown_0.7            
##  [91] promises_1.0.1           gridExtra_2.3           
##  [93] rio_0.5.10               vipor_0.4.5             
##  [95] boot_1.3-20              MASS_7.3-50             
##  [97] assertthat_0.2.0         rhdf5_2.24.0            
##  [99] rprojroot_1.3-2          rjson_0.2.18            
## [101] qlcMatrix_0.9.7          GenomeInfoDbData_1.1.0  
## [103] grid_3.5.0               class_7.3-14            
## [105] rmarkdown_1.9            DelayedMatrixStats_1.2.0
## [107] carData_3.0-1            Rtsne_0.13              
## [109] TTR_0.23-3               scatterplot3d_0.3-41    
## [111] shiny_1.0.5              ggbeeswarm_0.6.0
```
