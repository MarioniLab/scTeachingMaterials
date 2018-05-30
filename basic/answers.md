---
title: A simple analysis of single-cell RNA seq data with Bioconductor packages
author: Aaron Lun 
date: "2018-05-30"
output: 
  BiocStyle::html_document:
    fig_caption: false
---



# Introduction

This practical will perform some of the initial steps in a basic scRNA-seq data analysis, including: 

- Quality control on the cells
- Cell cycle phase assignment
- Normalization for cell-specific biases
- Modelling technical noise
- Dimensionality reduction
- Some visualization and clustering
- Differential gene detection

This will use R version 3.5.0 or higher, with a number of Bioconductor packages.
If you haven't downloaded and installed them already, you can do so by running the code below.
**This only needs to be done once** - the packages will be on your computer once installed, and can be loaded with `library`.


```r
source("https://bioconductor.org/biocLite.R")
biocLite(c("knitr", "BiocStyle", "org.Mm.eg.db", "scater", "Rtsne", 
    "TxDb.Mmusculus.UCSC.mm10.ensGene", "scran", "pheatmap"))
```

# Setting up the data 

## Reading in the counts

To demonstrate, we'll use a subset of the data from a study of mouse embryonic stem cells (mESCs).
These mESCs were cultured under various conditions - 2i (ground state), lif (serum) and a2i (alternative ground state).
To keep things simple, I've only taken only cells in the single batch that contains spike-in transcripts.
A more detailed description of the study is available at http://www.ebi.ac.uk/teichmann-srv/espresso/.

Our first task is to read in the counts and the associated metadata.
Here, both files are stored in the CSV format so we can just use the `read.csv` command.
(The `.gz` suffix just indicates that it is compressed, which is automatically handled by `read.csv`.)


```r
count.data <- read.csv("es_data.csv.gz", header=TRUE, row.names=1)
head(count.data[,1:10])
```

```
##                    ola_mES_2i_3_10.counts ola_mES_2i_3_11.counts
## ENSMUSG00000000001                    603                    140
## ENSMUSG00000000003                      0                      0
## ENSMUSG00000000028                    346                    407
## ENSMUSG00000000031                      0                      0
## ENSMUSG00000000037                      0                      0
## ENSMUSG00000000049                      0                      0
##                    ola_mES_2i_3_12.counts ola_mES_2i_3_13.counts
## ENSMUSG00000000001                    421                    173
## ENSMUSG00000000003                      0                      0
## ENSMUSG00000000028                    708                    276
## ENSMUSG00000000031                      0                      0
## ENSMUSG00000000037                      0                      0
## ENSMUSG00000000049                      0                      0
##                    ola_mES_2i_3_15.counts ola_mES_2i_3_16.counts
## ENSMUSG00000000001                    406                    296
## ENSMUSG00000000003                      0                      0
## ENSMUSG00000000028                    777                     17
## ENSMUSG00000000031                      0                      0
## ENSMUSG00000000037                      0                      0
## ENSMUSG00000000049                      0                      0
##                    ola_mES_2i_3_17.counts ola_mES_2i_3_18.counts
## ENSMUSG00000000001                     92                   1407
## ENSMUSG00000000003                      0                      0
## ENSMUSG00000000028                      3                   1933
## ENSMUSG00000000031                      0                      0
## ENSMUSG00000000037                      0                      0
## ENSMUSG00000000049                      0                      0
##                    ola_mES_2i_3_19.counts ola_mES_2i_3_2.counts
## ENSMUSG00000000001                    199                    48
## ENSMUSG00000000003                      0                     0
## ENSMUSG00000000028                    613                   866
## ENSMUSG00000000031                      0                     0
## ENSMUSG00000000037                      0                     0
## ENSMUSG00000000049                      0                     0
```

The original study used _HTSeq_ to assign reads to genes to obtain gene counts in each cell.
It's worth pointing out that _HTSeq_ puts some gunk at the end of the count matrix, e.g., number of unassigned reads.


```r
tail(count.data[,1:10])
```

```
##                        ola_mES_2i_3_10.counts ola_mES_2i_3_11.counts
## ERCC-00171                              68675                  60871
## __no_feature                          1402733                1352141
## __ambiguous                            104932                 141055
## __too_low_aQual                             0                      0
## __not_aligned                          645271                 599720
## __alignment_not_unique                      0                      0
##                        ola_mES_2i_3_12.counts ola_mES_2i_3_13.counts
## ERCC-00171                              32067                  35563
## __no_feature                          1162808                1341302
## __ambiguous                            106137                 111237
## __too_low_aQual                             0                      0
## __not_aligned                          459533                 528995
## __alignment_not_unique                      0                      0
##                        ola_mES_2i_3_15.counts ola_mES_2i_3_16.counts
## ERCC-00171                              27481                  59013
## __no_feature                          1051578                1495291
## __ambiguous                            109613                 138660
## __too_low_aQual                             0                      0
## __not_aligned                          465062                 695486
## __alignment_not_unique                      0                      0
##                        ola_mES_2i_3_17.counts ola_mES_2i_3_18.counts
## ERCC-00171                              88843                 135546
## __no_feature                          1210200                1835213
## __ambiguous                            155099                 221749
## __too_low_aQual                             0                      0
## __not_aligned                          738637                 821338
## __alignment_not_unique                      0                      0
##                        ola_mES_2i_3_19.counts ola_mES_2i_3_2.counts
## ERCC-00171                              32120                 87124
## __no_feature                          1218745               1818393
## __ambiguous                            120198                123300
## __too_low_aQual                             0                     0
## __not_aligned                          575089                592157
## __alignment_not_unique                      0                     0
```

We throw these out because they're not counts for actual genes.


```r
count.data <- count.data[!grepl("^_", rownames(count.data)),]
tail(count.data[,1:10])
```

```
##            ola_mES_2i_3_10.counts ola_mES_2i_3_11.counts
## ERCC-00163                    417                    388
## ERCC-00164                      0                      0
## ERCC-00165                   3323                   2862
## ERCC-00168                      0                      0
## ERCC-00170                   1181                    487
## ERCC-00171                  68675                  60871
##            ola_mES_2i_3_12.counts ola_mES_2i_3_13.counts
## ERCC-00163                    329                    167
## ERCC-00164                      2                      0
## ERCC-00165                   1089                   1186
## ERCC-00168                      0                     97
## ERCC-00170                    379                    201
## ERCC-00171                  32067                  35563
##            ola_mES_2i_3_15.counts ola_mES_2i_3_16.counts
## ERCC-00163                     64                    483
## ERCC-00164                      0                      0
## ERCC-00165                    636                   1726
## ERCC-00168                      0                     14
## ERCC-00170                    199                    333
## ERCC-00171                  27481                  59013
##            ola_mES_2i_3_17.counts ola_mES_2i_3_18.counts
## ERCC-00163                    431                    954
## ERCC-00164                      0                      0
## ERCC-00165                   2331                   4625
## ERCC-00168                      0                      0
## ERCC-00170                    379                   1714
## ERCC-00171                  88843                 135546
##            ola_mES_2i_3_19.counts ola_mES_2i_3_2.counts
## ERCC-00163                    127                   915
## ERCC-00164                      0                     0
## ERCC-00165                   1289                  2169
## ERCC-00168                      0                     0
## ERCC-00170                    130                   471
## ERCC-00171                  32120                 87124
```

<div class="alert alert-warning">
**Exercise:**

Regular expressions: what does the caret do?


```r
grepl("^_", "__UNMAPPED")
```

```
## [1] TRUE
```


```r
grepl("_", "I_AM_A_GENE")
```

```
## [1] TRUE
```


```r
grepl("^_", "I_AM_A_GENE")
```

```
## [1] FALSE
```
</div>


## Organizing the cell-based metadata

We read in the metadata:


```r
metadata <- read.csv("metadata.csv", header=TRUE, row.names=1)
head(metadata)
```

```
##                        Culture Batch
## ola_mES_2i_3_10.counts      2i     3
## ola_mES_2i_3_11.counts      2i     3
## ola_mES_2i_3_12.counts      2i     3
## ola_mES_2i_3_13.counts      2i     3
## ola_mES_2i_3_15.counts      2i     3
## ola_mES_2i_3_16.counts      2i     3
```

It's a good idea to check that the metadata actually matches up with the count data.
The `match()` command below ensures that the ordering of metadata rows are the same as the ordering of count columns.


```r
m <- match(colnames(count.data), rownames(metadata))
metadata <- metadata[m,]
head(metadata)
```

```
##                        Culture Batch
## ola_mES_2i_3_10.counts      2i     3
## ola_mES_2i_3_11.counts      2i     3
## ola_mES_2i_3_12.counts      2i     3
## ola_mES_2i_3_13.counts      2i     3
## ola_mES_2i_3_15.counts      2i     3
## ola_mES_2i_3_16.counts      2i     3
```

This information can be stored alongside the counts in a `SingleCellExperiment` object.
By storing everything together, we avoid book-keeping errors, e.g., when one matrix is subsetted and the other is not.


```r
library(SingleCellExperiment)
sce <- SingleCellExperiment(list(counts=as.matrix(count.data)), 
    colData=metadata)
sce    
```

```
## class: SingleCellExperiment 
## dim: 38653 204 
## metadata(0):
## assays(1): counts
## rownames(38653): ENSMUSG00000000001 ENSMUSG00000000003 ...
##   ERCC-00170 ERCC-00171
## rowData names(0):
## colnames(204): ola_mES_2i_3_10.counts ola_mES_2i_3_11.counts ...
##   ola_mES_lif_3_95.counts ola_mES_lif_3_96.counts
## colData names(2): Culture Batch
## reducedDimNames(0):
## spikeNames(0):
```

<div class="alert alert-warning">
**Exercise:** 

What column metadata are available for this dataset?


```r
colData(sce)
```

```
## DataFrame with 204 rows and 2 columns
##                          Culture     Batch
##                         <factor> <integer>
## ola_mES_2i_3_10.counts        2i         3
## ola_mES_2i_3_11.counts        2i         3
## ola_mES_2i_3_12.counts        2i         3
## ola_mES_2i_3_13.counts        2i         3
## ola_mES_2i_3_15.counts        2i         3
## ...                          ...       ...
## ola_mES_lif_3_90.counts      lif         3
## ola_mES_lif_3_92.counts      lif         3
## ola_mES_lif_3_94.counts      lif         3
## ola_mES_lif_3_95.counts      lif         3
## ola_mES_lif_3_96.counts      lif         3
```

How can you get column metadata values for each cell?


```r
sce$Culture
```

```
##   [1] 2i  2i  2i  2i  2i  2i  2i  2i  2i  2i  2i  2i  2i  2i  2i  2i  2i  2i 
##  [19] 2i  2i  2i  2i  2i  2i  2i  2i  2i  2i  2i  2i  2i  2i  2i  2i  2i  2i 
##  [37] 2i  2i  2i  2i  2i  2i  2i  2i  2i  2i  2i  2i  2i  2i  2i  2i  2i  2i 
##  [55] 2i  2i  2i  2i  2i  a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i
##  [73] a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i
##  [91] a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i
## [109] a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i a2i lif
## [127] lif lif lif lif lif lif lif lif lif lif lif lif lif lif lif lif lif lif
## [145] lif lif lif lif lif lif lif lif lif lif lif lif lif lif lif lif lif lif
## [163] lif lif lif lif lif lif lif lif lif lif lif lif lif lif lif lif lif lif
## [181] lif lif lif lif lif lif lif lif lif lif lif lif lif lif lif lif lif lif
## [199] lif lif lif lif lif lif
## Levels: 2i a2i lif
```

How do you get the counts?


```r
# str() => structure of the object
str(counts(sce))
```

```
##  int [1:38653, 1:204] 603 0 346 0 0 0 0 0 342 0 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:38653] "ENSMUSG00000000001" "ENSMUSG00000000003" "ENSMUSG00000000028" "ENSMUSG00000000031" ...
##   ..$ : chr [1:204] "ola_mES_2i_3_10.counts" "ola_mES_2i_3_11.counts" "ola_mES_2i_3_12.counts" "ola_mES_2i_3_13.counts" ...
```
</div>

## Adding gene-based annotation

We pull out annotation from *[org.Mm.eg.db](http://bioconductor.org/packages/org.Mm.eg.db)* to relate the ENSEMBL identifiers to the gene symbols.
The `mapIds` call just ensures that only one gene symbol is used if two symbols map to the same Ensembl ID.


```r
library(org.Mm.eg.db)
symbols <- mapIds(org.Mm.eg.db, keys=rownames(sce), keytype="ENSEMBL", column="SYMBOL")
anno <- data.frame(ENSEMBL=rownames(sce), SYMBOL=symbols, stringsAsFactors=FALSE)
head(anno)
```

```
##                               ENSEMBL SYMBOL
## ENSMUSG00000000001 ENSMUSG00000000001  Gnai3
## ENSMUSG00000000003 ENSMUSG00000000003   Pbsn
## ENSMUSG00000000028 ENSMUSG00000000028  Cdc45
## ENSMUSG00000000031 ENSMUSG00000000031    H19
## ENSMUSG00000000037 ENSMUSG00000000037  Scml2
## ENSMUSG00000000049 ENSMUSG00000000049   Apoh
```

To identify which rows correspond to mitochondrial genes, we need to use extra annotation describing the genomic location of each gene.
For Ensembl, this involves using the *[TxDb.Mmusculus.UCSC.mm10.ensGene](http://bioconductor.org/packages/TxDb.Mmusculus.UCSC.mm10.ensGene)* package.


```r
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
location <- mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene, keys=rownames(sce),
    column="CDSCHROM", keytype="GENEID")
anno$CHR <- location
table(anno$CHR)
```

```
## 
##                 chr1 chr1_GL456210_random chr1_GL456211_random 
##                 1237                    5                    6 
## chr1_GL456212_random chr1_GL456221_random                chr10 
##                    2                    7                 1039 
##                chr11                chr12                chr13 
##                 1629                  809                  867 
##                chr14                chr15                chr16 
##                 1028                  794                  667 
##                chr17                chr18                chr19 
##                 1090                  492                  719 
##                 chr2                 chr3                 chr4 
##                 1834                 1050                 1349 
## chr4_GL456216_random chr4_GL456350_random chr4_JH584292_random 
##                    1                    7                    1 
## chr4_JH584293_random chr4_JH584294_random chr4_JH584295_random 
##                   12                    8                    1 
##                 chr5 chr5_GL456354_random chr5_JH584296_random 
##                 1303                    5                    3 
## chr5_JH584297_random chr5_JH584298_random chr5_JH584299_random 
##                    3                    3                   12 
##                 chr6                 chr7 chr7_GL456219_random 
##                 1284                 2044                    2 
##                 chr8                 chr9                 chrM 
##                 1072                 1263                   13 
##       chrUn_JH584304                 chrX chrX_GL456233_random 
##                    1                  945                    4 
##                 chrY chrY_JH584303_random 
##                  377                    1
```

We add all of this information to our `SingleCellExperiment` object.


```r
rowData(sce) <- anno
rowData(sce)
```

```
## DataFrame with 38653 rows and 3 columns
##                  ENSEMBL      SYMBOL         CHR
##              <character> <character> <character>
## 1     ENSMUSG00000000001       Gnai3        chr3
## 2     ENSMUSG00000000003        Pbsn        chrX
## 3     ENSMUSG00000000028       Cdc45       chr16
## 4     ENSMUSG00000000031         H19          NA
## 5     ENSMUSG00000000037       Scml2        chrX
## ...                  ...         ...         ...
## 38649         ERCC-00164          NA          NA
## 38650         ERCC-00165          NA          NA
## 38651         ERCC-00168          NA          NA
## 38652         ERCC-00170          NA          NA
## 38653         ERCC-00171          NA          NA
```

<div class="alert alert-warning">
**Exercise:**

How do you recover the row metadata for each gene?


```r
head(rowData(sce)$SYMBOL)
```

```
## [1] "Gnai3" "Pbsn"  "Cdc45" "H19"   "Scml2" "Apoh"
```

How do you add an extra row metadata variable?


```r
BLAH <- runif(nrow(sce))
rowData(sce)$BLAH <- BLAH
rowData(sce)
```

```
## DataFrame with 38653 rows and 4 columns
##                  ENSEMBL      SYMBOL         CHR               BLAH
##              <character> <character> <character>          <numeric>
## 1     ENSMUSG00000000001       Gnai3        chr3  0.185841772938147
## 2     ENSMUSG00000000003        Pbsn        chrX 0.0107907785568386
## 3     ENSMUSG00000000028       Cdc45       chr16   0.67441970971413
## 4     ENSMUSG00000000031         H19          NA   0.47338613960892
## 5     ENSMUSG00000000037       Scml2        chrX  0.531499963486567
## ...                  ...         ...         ...                ...
## 38649         ERCC-00164          NA          NA  0.155990440864116
## 38650         ERCC-00165          NA          NA  0.615467643830925
## 38651         ERCC-00168          NA          NA 0.0438267334830016
## 38652         ERCC-00170          NA          NA  0.282002961263061
## 38653         ERCC-00171          NA          NA  0.540163872996345
```
</div>

Identification of rows that correspond to spike-in transcripts is much easier, given that the ERCC spike-ins were used.
(If you're doing this on the gene symbols, beware of the human gene family that also starts with "ERCC".)


```r
is.spike <- grepl("^ERCC", rownames(sce))
sum(is.spike)
```

```
## [1] 92
```

Note that we need to explicitly indicate that the ERCC set is, in fact, a spike-in set.
This is necessary as spike-ins require special treatment in some downstream steps such as variance estimation and normalization.


```r
isSpike(sce, "ERCC") <- is.spike
sce
```

```
## class: SingleCellExperiment 
## dim: 38653 204 
## metadata(0):
## assays(1): counts
## rownames(38653): ENSMUSG00000000001 ENSMUSG00000000003 ...
##   ERCC-00170 ERCC-00171
## rowData names(4): ENSEMBL SYMBOL CHR BLAH
## colnames(204): ola_mES_2i_3_10.counts ola_mES_2i_3_11.counts ...
##   ola_mES_lif_3_95.counts ola_mES_lif_3_96.counts
## colData names(2): Culture Batch
## reducedDimNames(0):
## spikeNames(1): ERCC
```

To make things easier to interpret, we'll use the gene symbols as row names.
This requires some fiddling to avoid non-unique gene symbols (in which case we paste the Ensembl ID after it) or missing gene symbols (in which case we use the Ensembl ID).
*[scater](http://bioconductor.org/packages/scater)* provides the `uniquifyFeatureNames()` function to do this conveniently:


```r
library(scater)
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ENSEMBL, rowData(sce)$SYMBOL)
head(rownames(sce), 50)
```

```
##  [1] "Gnai3"    "Pbsn"     "Cdc45"    "H19"      "Scml2"    "Apoh"    
##  [7] "Narf"     "Cav2"     "Klf6"     "Scmh1"    "Cox5a"    "Tbx2"    
## [13] "Tbx4"     "Zfy2"     "Ngfr"     "Wnt3"     "Wnt9a"    "Fer"     
## [19] "Xpo6"     "Tfe3"     "Axin2"    "Brat1"    "Gna12"    "Slc22a18"
## [25] "Itgb2l"   "Igsf5"    "Pih1d2"   "Dlat"     "Sdhd"     "Fgf23"   
## [31] "Fgf6"     "Ccnd2"    "Gpr107"   "Nalcn"    "Btbd17"   "Slfn4"   
## [37] "Th"       "Ins2"     "Scnn1g"   "Drp2"     "Tspan32"  "Lhx2"    
## [43] "Clec2g"   "Gmpr"     "Glra1"    "Mid2"     "Trim25"   "Dgke"    
## [49] "Scpep1"   "Mnt"
```

**<rant>**

If you're releasing your own data with a count table, be sure to put a stable ID as the row name, e.g., Ensembl or Entrez.
People can always easily convert from these IDs to gene symbols; it is much harder to go from symbols to umambiguous IDs.
(And put _just_ the ID; don't paste the ID with something else, as this makes life difficult for people who want to use your data.)

**</rant>**

## Quality control on the cells 

We need to remove low-quality cells to ensure that technical effects do not distort downstream analyses.
We have a number of simple metrics that we use to assess quality:

- Library size 
- Number of expressed features
- Proportion of spike-in reads
- Proportion of mitochondrial reads (if spike-ins are not available)

Together, these catch failures in cDNA capture or sequencing, especially for the endogenous transcripts.
Spike-in information is automatically extracted from the `SingleCellExperiment` object, but we need to explicitly specify genes on the mitochondrial genome:


```r
is.mito <- which(rowData(sce)$CHR == "chrM")
is.mito
```

```
##  [1] 18158 18162 18168 18171 18173 18174 18175 18177 18179 18183 18184 18186
## [13] 19142
```

<div class="alert alert-warning">
**Exercise:**

Why which()? See this neat trick:


```r
Y <- LETTERS[1:5]
X <- c(TRUE, NA, FALSE, NA, TRUE)

Y[X]
```

```
## [1] "A" NA  NA  "E"
```

```r
Y[which(X)]
```

```
## [1] "A" "E"
```

Why use it with 'is.mito'?


```r
sum(is.na(rowData(sce)$CHR))
```

```
## [1] 15664
```
</div>

For each cell, we calculate quality control metrics such as the total number of counts or the proportion of counts in mitochondrial genes or spike-in transcripts.
These are stored in the `colData` of the `SingleCellExperiment` for future reference.


```r
library(scater)
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=is.mito)) 
head(colnames(colData(sce)))
```

```
## [1] "Culture"                        "Batch"                         
## [3] "is_cell_control"                "total_features_by_counts"      
## [5] "log10_total_features_by_counts" "total_counts"
```

We create plots of these statistics below for each culture condition.


```r
multiplot(cols=2,
    plotColData(sce, x="Culture", y="log10_total_counts"),
    plotColData(sce, x="Culture", y="total_features"),
    plotColData(sce, x="Culture", y="pct_counts_ERCC"),
    plotColData(sce, x="Culture", y="pct_counts_Mt")
)
```

<img src="answers_files/figure-html/qualplot-1.png" width="960" />

We assume that, in each culture condition, most cells are of high quality.
Cells with outlier values for each of the QC metrics are identified based on some number of MADs from the median value.


```r
low.lib <- isOutlier(sce$total_counts, log=TRUE, nmads=3, type="lower", batch=sce$Culture) # using log-values, here.
low.nfeatures <- isOutlier(sce$total_features, log=TRUE, nmads=3, type="lower", batch=sce$Culture)
high.ercc <- isOutlier(sce$pct_counts_ERCC, nmads=3, type="higher", batch=sce$Culture)
data.frame(LowLib=sum(low.lib), LowNgenes=sum(low.nfeatures), HighSpike=sum(high.ercc))
```

```
##   LowLib LowNgenes HighSpike
## 1      2         5         9
```

<div class="alert alert-warning">
**Exercise:**

What does `batch=` do in `isOutlier()`?


```r
x <- c(1,2,0,1,3,2,1,3,0,1,3,2,4,2,0,1,1,5,2,50,52,53,51)
b <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2)

isOutlier(x, batch=b) # last four are outliers
```

```
##  [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [13] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
```

```r
isOutlier(x, batch=b) # last four are not outliers
```

```
##  [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [13] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
```
</div>

We only retain cells that pass all of the specified criteria.
Of course, this involves some assumptions about independence from biology. 
(For example, don't use the mitochondrial proportions if the number/activity of mitochondria changes between cell types.)


```r
discard <- low.lib | low.nfeatures | high.ercc
data.frame(TotalLost=sum(discard), TotalLeft=sum(!discard))
```

```
##   TotalLost TotalLeft
## 1        12       192
```

We toss out the cells that we consider to be low-quality, and keep the rest.
Here, most cells are retained (which makes sense, as some QC was already applied to the published data).


```r
sce <- sce[,!discard]
ncol(sce)
```

```
## [1] 192
```

Of course, more sophisticated QC procedures can be used, e.g., *[cellity](http://bioconductor.org/packages/cellity)*.
We stick to the simple stuff because it's easier to interpret and troubleshoot.

# Classification of cell cycle phase 

We use the `cyclone` method to classify cells into different cell cycle phases (http://dx.doi.org/10.1016/j.ymeth.2015.06.021).
This was previously trained on some mouse data with known cell cycle classifications - we can get the trained classifier with `mm.pairs`.


```r
library(scran)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
str(mm.pairs)
```

```
## List of 3
##  $ G1 :'data.frame':	12052 obs. of  2 variables:
##   ..$ first : chr [1:12052] "ENSMUSG00000000001" "ENSMUSG00000000001" "ENSMUSG00000000001" "ENSMUSG00000000001" ...
##   ..$ second: chr [1:12052] "ENSMUSG00000001785" "ENSMUSG00000005470" "ENSMUSG00000012443" "ENSMUSG00000015120" ...
##  $ S  :'data.frame':	6459 obs. of  2 variables:
##   ..$ first : chr [1:6459] "ENSMUSG00000000001" "ENSMUSG00000000001" "ENSMUSG00000000001" "ENSMUSG00000000001" ...
##   ..$ second: chr [1:6459] "ENSMUSG00000002014" "ENSMUSG00000004771" "ENSMUSG00000007656" "ENSMUSG00000015749" ...
##  $ G2M:'data.frame':	9981 obs. of  2 variables:
##   ..$ first : chr [1:9981] "ENSMUSG00000000001" "ENSMUSG00000000001" "ENSMUSG00000000001" "ENSMUSG00000000001" ...
##   ..$ second: chr [1:9981] "ENSMUSG00000014402" "ENSMUSG00000017499" "ENSMUSG00000022432" "ENSMUSG00000027285" ...
```

The method uses a randomization step, so we use `set.seed()` to obtain consistent results from different runs.


```r
set.seed(100)
assignments <- cyclone(sce, mm.pairs, gene.names=rowData(sce)$ENSEMBL)
plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16)
```

<img src="answers_files/figure-html/unnamed-chunk-31-1.png" width="100%" />

Cells are classified as being in G1 phase if the G1 score is above 0.5 and greater than the G2/M score; 
    in G2/M phase if the G2/M score is above 0.5 and greater than the G1 score; 
    and in S phase if neither score is above 0.5.
It seems that most cells here are in the S phase.

<div class="alert alert-warning">
**Exercise:**


```r
# How do we figure out the specific calls?
table(assignments$phases)
```

```
## 
##  G1 G2M   S 
##  30  31 131
```
</div>

We store the cell cycle phases for future use, e.g., to check whether any downstream results are driven by cell cycle effects.
Phase assignment is difficult so I generally don't use the assignments for anything more than diagnostics.


```r
sce$cycle_phase <- assignments$phases
tail(colnames(colData(sce)))
```

```
## [1] "pct_counts_ERCC"                    "pct_counts_in_top_50_features_ERCC"
## [3] "total_features_ERCC"                "log10_total_features_ERCC"         
## [5] "pct_counts_top_50_features_ERCC"    "cycle_phase"
```

# Examining the genes

We inspect the distribution of log-mean counts across all genes.
The peak represents the bulk of moderately expressed genes while the rectangular component corresponds to lowly expressed genes.


```r
ave.counts <- calcAverage(sce)
hist(log10(ave.counts), breaks=100, main="", col="grey80",
    xlab=expression(Log[10]~"average count"))
```

<img src="answers_files/figure-html/unnamed-chunk-34-1.png" width="100%" />

We also look at the identities of the most highly expressed genes.
This should generally be dominated by constitutively expressed transcripts, such as those for ribosomal or mitochondrial proteins.
The presence of other classes of features may be cause for concern if they are not consistent with expected biology.


```r
plotHighestExprs(sce, n=50)
```

<img src="answers_files/figure-html/unnamed-chunk-35-1.png" width="100%" />

We can also have a look at the number of cells expressing each gene. 
This is usually well-correlated to the average expression of each gene.


```r
numcells <- nexprs(sce, byrow=TRUE)
smoothScatter(log2(ave.counts), numcells, 
    xlab=expression(Log[2]~"average count"), 
    ylab="Number of expressing cells")
```

<img src="answers_files/figure-html/unnamed-chunk-36-1.png" width="100%" />

We discard genes that are not expressed in any cell, as these are obviously uninformative.


```r
sce <- sce[numcells > 0,]
summary(numcells > 0)
```

```
##    Mode   FALSE    TRUE 
## logical   10359   28294
```

# Normalization of cell-specific biases

## For the endogenous genes

We apply the deconvolution method (https://dx.doi.org/10.1186/s13059-016-0947-7) to eliminate biases in the counts for the endogenous transcripts.
This computes size factors for each cell representing the scaling bias in the counts.
Some filtering is required to remove low-abundance genes prior to normalization, see `min.mean=`.


```r
sce <- computeSumFactors(sce, min.mean=1) 
summary(sizeFactors(sce))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.2418  0.6464  0.9461  1.0000  1.2742  2.3574
```

<div class="alert alert-warning">
**Exercise:**

It's a good idea to plot the size factors against the library sizes as a sanity check.
There should be some positive correlation.


```r
plot(sizeFactors(sce), sce$total_counts/1e6, 
    log="xy", ylab="Library size (millions)", xlab="Size factor")
```

<img src="answers_files/figure-html/unnamed-chunk-39-1.png" width="100%" />
</div>

For highly heterogeneous data sets with multiple cell types, we should cluster first with `quickCluster()`, and then normalize within each cluster with `cluster=`.
This avoids trying to pool very different cells together during the deconvolution process.
We don't bother doing this here, as all the cells are mESCs.

## Computing separate size factors for spike-in transcripts

We need to normalize spike-in transcripts separately as they are not subject to all biases affecting endogenous transcripts - in particular, total RNA content.
Applying the endogenous size factors would over-normalize, so we define separate size factors for the spike-ins.
This is simply defined for each cell as the total count across all transcripts in the spike-in set.


```r
sce <- computeSpikeFactors(sce, type="ERCC", general.use=FALSE)
summary(sizeFactors(sce, "ERCC"))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.1790  0.5174  0.8289  1.0000  1.3044  4.2152
```

These size factors are stored in a separate field of the `SingleCellExperiment` object by setting `general.use=FALSE` in `computeSpikeFactors`.
This ensures that they will only be used with the spike-in transcripts but not the endogenous genes.

<div class="alert alert-warning">
**Exercise:**

The two sets of size factors tend to agree less due to the effects of heterogeneity in total RNA content between cells - this is expected.


```r
plot(sizeFactors(sce, 'ERCC'), sizeFactors(sce),
    log="xy", xlab="Size factor (ERCC)", ylab="Size factor (genes)")
```

<img src="answers_files/figure-html/unnamed-chunk-41-1.png" width="100%" />

Note that if you _do_ want to use the spike-in size factors to normalize all genes, set `general.use=TRUE` instead.
</div>

## Applying the size factors to normalize gene expression

Counts are transformed into normalized log-expression values for use in downstream analyses.
Each value is defined as the log-ratio of each count to the size factor for the corresponding cell, after adding a pseudo-count of 1 to avoid undefined values at zero counts.
Division of the counts for each gene by its appropriate size factor ensures that any cell-specific biases are removed.


```r
sce <- normalize(sce)
sce
```

```
## class: SingleCellExperiment 
## dim: 28294 192 
## metadata(1): log.exprs.offset
## assays(2): counts logcounts
## rownames(28294): Gnai3 Pbsn ... ERCC-00170 ERCC-00171
## rowData names(15): ENSEMBL SYMBOL ... n_cells_counts
##   pct_dropout_counts
## colnames(192): ola_mES_2i_3_10.counts ola_mES_2i_3_11.counts ...
##   ola_mES_lif_3_95.counts ola_mES_lif_3_96.counts
## colData names(60): Culture Batch ... pct_counts_top_50_features_ERCC
##   cycle_phase
## reducedDimNames(0):
## spikeNames(1): ERCC
```

The log-transformation provides some measure of variance stabilization, so that high-abundance genes with large variances do not dominate downstream analyses.
The computed values are stored as an `logcounts` matrix in addition to the other assay elements.

<div class="alert alert-warning">
**Exercise:**

How do we find the available assays?


```r
assayNames(sce)
```

```
## [1] "counts"    "logcounts"
```

How to we get a particular assay?


```r
str(assay(sce, 'logcounts'))
```

```
##  num [1:28294, 1:192] 9.74 0 8.94 0 0 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:28294] "Gnai3" "Pbsn" "Cdc45" "H19" ...
##   ..$ : chr [1:192] "ola_mES_2i_3_10.counts" "ola_mES_2i_3_11.counts" "ola_mES_2i_3_12.counts" "ola_mES_2i_3_13.counts" ...
```

How do we get the logcounts?


```r
str(logcounts(sce))
```

```
##  num [1:28294, 1:192] 9.74 0 8.94 0 0 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:28294] "Gnai3" "Pbsn" "Cdc45" "H19" ...
##   ..$ : chr [1:192] "ola_mES_2i_3_10.counts" "ola_mES_2i_3_11.counts" "ola_mES_2i_3_12.counts" "ola_mES_2i_3_13.counts" ...
```
</div>

# Modelling the technical component of variation

We identify HVGs to focus on the genes that are driving heterogeneity across the population of cells.
This requires estimation of the variance in expression for each gene, followed by decomposition of the variance into biological and technical components.


```r
var.fit <- trendVar(sce, method="loess", loess.args=list(span=0.2))
var.out <- decomposeVar(sce, var.fit)
head(var.out)    
```

```
## DataFrame with 6 rows and 6 columns
##                      mean               total                 bio
##                 <numeric>           <numeric>           <numeric>
## Gnai3    8.41578552934507    4.96161788483393    3.90814886233558
## Pbsn  0.00427098707057064 0.00350233546694046 -0.0271161719198812
## Cdc45    7.47181265135016    7.83090605228553    5.15755613484333
## H19      1.93713736788656    11.9647715438154    4.47693950090266
## Scml2   0.674288582888631    3.32629495103183 -0.0637453227711968
## Narf     1.99404374962859    8.57829279442094   0.965781826445298
##                     tech              p.value                  FDR
##                <numeric>            <numeric>            <numeric>
## Gnai3   1.05346902249835 2.80080880349774e-92 3.49518631349765e-90
## Pbsn  0.0306185073868216                    1                    1
## Cdc45    2.6733499174422 7.60355389782798e-38 2.96192031188457e-36
## H19     7.48783204291271 2.85999179397056e-07 3.45588468574771e-06
## Scml2   3.39004027380303    0.559913621002494                    1
## Narf    7.61251096797564    0.110341829204011    0.531053004955754
```

We can have a look at the fitted trend to the spike-in variances.
Some tinkering may be required to get a good fit, usually by modifying `span=`.


```r
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
    ylab="Variance of log-expression")
curve(var.fit$trend(x), add=TRUE, col="dodgerblue", lwd=2)
cur.spike <- isSpike(sce)
points(var.out$mean[cur.spike], var.out$total[cur.spike], col="red", pch=16)
```

<img src="answers_files/figure-html/hvgplothsc-1.png" width="100%" />

If you don't have spike-ins, you can fit the trend to the variances of the genes with `use.spikes=FALSE` (but this probably overestimates the technical component).

<div class="alert alert-warning">
**Exercise:**

It's wise to check the distribution of expression values for the top HVGs to ensure that the variance estimate is not being dominated by one or two outlier cells.


```r
top <- rownames(var.out)[order(var.out$bio, decreasing=TRUE)[1:10]]
plotExpression(sce, features = top)
```

<img src="answers_files/figure-html/unnamed-chunk-47-1.png" width="100%" />
</div>

An alternative approach is to use `technicalCV2`, which implements the method described by Brennecke _et al._ (http://dx.doi.org/10.1038/nmeth.2645).
This has some pros and cons compared to the log-variance method described above.

# Dimensionality reduction based on the technical noise 

We use all genes with a positive biological component in `denoisePCA`.
This performs a principal components analysis on the expression profiles, choosing the number of PCs to retain based on the total technical noise in the data set.
The idea is to discard later PCs that contain random technical noise, thus enriching for early biological signal (and also reducing work in downstream steps).


```r
sce <- denoisePCA(sce, technical=var.fit$trend)
sce
```

```
## class: SingleCellExperiment 
## dim: 28294 192 
## metadata(1): log.exprs.offset
## assays(2): counts logcounts
## rownames(28294): Gnai3 Pbsn ... ERCC-00170 ERCC-00171
## rowData names(15): ENSEMBL SYMBOL ... n_cells_counts
##   pct_dropout_counts
## colnames(192): ola_mES_2i_3_10.counts ola_mES_2i_3_11.counts ...
##   ola_mES_lif_3_95.counts ola_mES_lif_3_96.counts
## colData names(60): Culture Batch ... pct_counts_top_50_features_ERCC
##   cycle_phase
## reducedDimNames(1): PCA
## spikeNames(1): ERCC
```

We can have a look at the PCs directly:


```r
pcs <- reducedDim(sce) # stored in the object
dim(pcs) # Cells are rows, PCs are columns
```

```
## [1] 192  20
```

Creating pairwise plots between the first four PCs:


```r
plotPCA(sce, ncomponents=4, colour_by="Culture")
```

<img src="answers_files/figure-html/pcaplothsc-1.png" width="960" />

<div class="alert alert-warning">
**Exercise:**

What if we want to look at the proportion of variance explained by each PC?


```r
prop.var <- attr(pcs, "percentVar")
plot(prop.var, xlab="PC", ylab="Proportion of variance explained")
```

<img src="answers_files/figure-html/unnamed-chunk-50-1.png" width="100%" />
</div>

We also use _t_-SNE, which is very good at displaying distinct clusters of cells and resolving complex structure.
Note that we use the PCA results as "denoised expression values" for input into downstream functions like _t_-SNE.
This is valid as it is the distance between cells that is important.


```r
sce <- runTSNE(sce, use_dimred="PCA", perplexity=30, rand_seed=100)
plotTSNE(sce, colour_by="Culture")
```

<img src="answers_files/figure-html/unnamed-chunk-51-1.png" width="100%" />

However, _t_-SNE is stochastic and more complicated than PCA.
Testing different settings of the "perplexity" parameter is recommended, as well as running multiple times to check that the conclusions are the same.
(Check out http://distill.pub/2016/misread-tsne/ for examples of odd _t_-SNE behaviour.)

# Clustering into putative subpopulations

Clearly there's some structure here.
Here, we know that they're associated with the different culture conditions, but if we didn't we'd have to cluster the cells.
A quick and dirty approach with hierarchical clustering on Euclidean distances:


```r
my.dist <- dist(pcs) # Making a distance matrix
my.tree <- hclust(my.dist, method="ward.D2") # Building a tree
plot(my.tree)
```

<img src="answers_files/figure-html/unnamed-chunk-52-1.png" width="100%" />

Cutting the tree into three clusters, and seeing how these correspond to the known culture conditions:


```r
my.clusters <- cutree(my.tree, k=3) 
table(my.clusters, sce$Culture)
```

```
##            
## my.clusters 2i a2i lif
##           1 51   0   0
##           2  6  63   0
##           3  0   0  72
```

<div class="alert alert-warning">
**Exercise:**

How do I make a _t_-SNE plot coloured by cluster?


```r
sce$Cluster <- factor(my.clusters)
plotTSNE(sce, colour_by='Cluster')
```

<img src="answers_files/figure-html/unnamed-chunk-54-1.png" width="100%" />
</div>

We use a heatmap to visualize the expression profiles of the top 100 genes with the largest biological components.
We see "blocks" in expression that correspond nicely to the known culture conditions.
We possibly could have subclustered further, in which case we would subset and repeat the above process.


```r
plotHeatmap(sce, 
    features=rownames(var.out)[order(var.out$bio, decreasing=TRUE)[1:100]],
    cluster_cols=my.tree, colour_columns_by="Culture",
    center=TRUE, symmetric=TRUE)
```

<img src="answers_files/figure-html/unnamed-chunk-55-1.png" width="100%" />

When clustering, it is often useful to look at silhouette plots to assess cluster separatedness.
Each bar corresponds to a cell, and is proportional to the difference in the average distances to all other cells in the same cluster versus cells in the nearest neighbouring cluster.
A good gauge for the number of clusters is that which maximizes the average silhouette width.


```r
library(cluster)
par(mfrow=c(2,2))
for (k in 2:5) { 
    example.clusters <- cutree(my.tree, k=k)
    sil <- silhouette(example.clusters, dist=my.dist)
    plot(sil, col=rainbow(k)[sort(sil[,1])])
}
```

<img src="answers_files/figure-html/silhouette-1.png" width="960" />

Other options for clustering are:

- Use the `cutreeDynamic()` function in the *[dynamicTreeCut](https://CRAN.R-project.org/package=dynamicTreeCut)* package, for toplogy-aware cutting of the tree.
- Use graph-based methods such as `buildSNNGraph()` or `buildKNNGraph()`, followed by clustering methods from *[igraph](https://CRAN.R-project.org/package=igraph)*.
- Use methods with pre-specified number of clusters, e.g., k-means with `kmeans()` and *[SC3](http://bioconductor.org/packages/SC3)*, self-organizing maps in *[flowSOM](http://bioconductor.org/packages/flowSOM)*. 

# Identifying marker genes between subpopulations

We use the `findMarkers` function to detect differences between clusters.
This will perform pairwise DE analyses between clusters, and consolidate the results into a single table of marker genes per cluster.
The tricky part is how to summarize results from many pairwise comparisons into a single ranking of genes.


```r
out <- findMarkers(sce, clusters=my.clusters)
```

Consider the results of cluster 2:


```r
marker.set <- out[[2]] # Marker set for cluster 2
head(marker.set)
```

```
## DataFrame with 6 rows and 4 columns
##                          Top                  FDR           logFC.1
##                    <integer>            <numeric>         <numeric>
## Lamb1                      1 1.39090399942238e-39 -1.22364166542808
## ENSMUSG00000064193         1  9.6141469969557e-22   1.6251199688414
## Rps18                      2 3.09529239052511e-35 0.245633601063361
## Trh                        2  3.5543418027953e-20  6.97830927490577
## Peg10                      3 1.00759539494572e-34 0.668232809427393
## Lin28a                     3 4.24938699170752e-20  6.85247776991814
##                              logFC.3
##                            <numeric>
## Lamb1               -6.8640626903952
## ENSMUSG00000064193 -0.72710300905236
## Rps18              0.857216361608522
## Trh                 5.38279175581989
## Peg10              -5.98502687483942
## Lin28a             -1.41528446825544
```

<div class="alert alert-warning">
**Exercise:**

What does the `Top` parameter represent?


```r
marker.set[marker.set$Top <= 1,]
```

```
## DataFrame with 2 rows and 4 columns
##                          Top                  FDR           logFC.1
##                    <integer>            <numeric>         <numeric>
## Lamb1                      1 1.39090399942238e-39 -1.22364166542808
## ENSMUSG00000064193         1  9.6141469969557e-22   1.6251199688414
##                              logFC.3
##                            <numeric>
## Lamb1               -6.8640626903952
## ENSMUSG00000064193 -0.72710300905236
```

Try it yourself:


```r
marker.set[marker.set$Top <= 10,]
```

```
## DataFrame with 20 rows and 4 columns
##                          Top                  FDR             logFC.1
##                    <integer>            <numeric>           <numeric>
## Lamb1                      1 1.39090399942238e-39   -1.22364166542808
## ENSMUSG00000064193         1  9.6141469969557e-22     1.6251199688414
## Rps18                      2 3.09529239052511e-35   0.245633601063361
## Trh                        2  3.5543418027953e-20    6.97830927490577
## Peg10                      3 1.00759539494572e-34   0.668232809427393
## ...                      ...                  ...                 ...
## Adam23                     8 2.59691517742769e-19    5.26922077985474
## ENSMUSG00000089944         9 3.48726577391293e-32 -0.0406834298381531
## Aqp3                       9 1.40092859566935e-18   -5.05053153350963
## Rpl18a                    10 9.91614462407516e-32  -0.347637768265264
## Ddx3y                     10 5.89941962770506e-18   -6.38529774474118
##                               logFC.3
##                             <numeric>
## Lamb1                -6.8640626903952
## ENSMUSG00000064193  -0.72710300905236
## Rps18               0.857216361608522
## Trh                  5.38279175581989
## Peg10               -5.98502687483942
## ...                               ...
## Adam23             -0.387397863543287
## ENSMUSG00000089944   1.99159019432366
## Aqp3                  2.7132371801288
## Rpl18a              0.955866793863532
## Ddx3y                0.58157788618842
```
</div>

We can visualize this more clearly with a heatmap of the top 20 genes.


```r
plotHeatmap(sce, 
    features=rownames(marker.set)[marker.set$Top <= 10],
    cluster_cols=my.tree, colour_columns_by="Cluster",
    center=TRUE, symmetric=TRUE)
```

<img src="answers_files/figure-html/unnamed-chunk-60-1.png" width="100%" />

A valid alternative strategy is to detect marker genes that are uniquely up-regulated or down-regulated in each cluster (set `pval.type="all"`).
However, be aware that no such genes may exist.
For example, in a mixed population of T cells, you could have CD4^+^ T cells, CD8^+^ T cells, double negative and double positive cells.
If each of these formed a cluster, and we only looked for unique genes, neither CD4 or CD8 would be detected!

# Additional comments

It's a good idea to save the `SingleCellExperiment` object to file with the `saveRDS` function.
The object can then be easily restored into new R sessions using the `readRDS` function.


```r
saveRDS(file="data.rds", sce)
copy <- readRDS("data.rds")
copy
```

```
## class: SingleCellExperiment 
## dim: 28294 192 
## metadata(1): log.exprs.offset
## assays(2): counts logcounts
## rownames(28294): Gnai3 Pbsn ... ERCC-00170 ERCC-00171
## rowData names(15): ENSEMBL SYMBOL ... n_cells_counts
##   pct_dropout_counts
## colnames(192): ola_mES_2i_3_10.counts ola_mES_2i_3_11.counts ...
##   ola_mES_lif_3_95.counts ola_mES_lif_3_96.counts
## colData names(61): Culture Batch ... cycle_phase Cluster
## reducedDimNames(2): PCA TSNE
## spikeNames(1): ERCC
```

Data within it can be extracted and used for more complex analyses.

- Droplet-based or UMI data analysis. 
- Batch correction, see `?mnnCorrect`.
- Trajectory reconstruction, see *[destiny](http://bioconductor.org/packages/destiny)* and *[monocle](http://bioconductor.org/packages/monocle)*.

See https://www.bioconductor.org/packages/devel/workflows/vignettes/simpleSingleCell/inst/doc/work-0-intro.html for more details.

Meanwhile, show the session information for record-keeping:


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
##  [1] cluster_2.0.7-1                       
##  [2] scran_1.8.2                           
##  [3] scater_1.8.0                          
##  [4] ggplot2_2.2.1                         
##  [5] TxDb.Mmusculus.UCSC.mm10.ensGene_3.4.0
##  [6] GenomicFeatures_1.32.0                
##  [7] org.Mm.eg.db_3.6.0                    
##  [8] AnnotationDbi_1.42.1                  
##  [9] SingleCellExperiment_1.2.0            
## [10] SummarizedExperiment_1.10.1           
## [11] DelayedArray_0.6.0                    
## [12] BiocParallel_1.14.1                   
## [13] matrixStats_0.53.1                    
## [14] Biobase_2.40.0                        
## [15] GenomicRanges_1.32.3                  
## [16] GenomeInfoDb_1.16.0                   
## [17] IRanges_2.14.10                       
## [18] S4Vectors_0.18.2                      
## [19] BiocGenerics_0.26.0                   
## [20] knitr_1.20                            
## [21] BiocStyle_2.8.1                       
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-6             bit64_0.9-7             
##  [3] RColorBrewer_1.1-2       progress_1.1.2          
##  [5] httr_1.3.1               rprojroot_1.3-2         
##  [7] dynamicTreeCut_1.63-1    tools_3.5.0             
##  [9] backports_1.1.2          DT_0.4                  
## [11] R6_2.2.2                 KernSmooth_2.23-15      
## [13] vipor_0.4.5              DBI_1.0.0               
## [15] lazyeval_0.2.1           colorspace_1.3-2        
## [17] tidyselect_0.2.4         gridExtra_2.3           
## [19] prettyunits_1.0.2        bit_1.1-13              
## [21] compiler_3.5.0           labeling_0.3            
## [23] rtracklayer_1.40.2       bookdown_0.7            
## [25] scales_0.5.0             stringr_1.3.1           
## [27] digest_0.6.15            Rsamtools_1.32.0        
## [29] rmarkdown_1.9            XVector_0.20.0          
## [31] pkgconfig_2.0.1          htmltools_0.3.6         
## [33] limma_3.36.1             htmlwidgets_1.2         
## [35] rlang_0.2.0              RSQLite_2.1.1           
## [37] FNN_1.1                  shiny_1.1.0             
## [39] DelayedMatrixStats_1.2.0 bindr_0.1.1             
## [41] dplyr_0.7.5              RCurl_1.95-4.10         
## [43] magrittr_1.5             GenomeInfoDbData_1.1.0  
## [45] Matrix_1.2-14            Rcpp_0.12.17            
## [47] ggbeeswarm_0.6.0         munsell_0.4.3           
## [49] Rhdf5lib_1.2.1           viridis_0.5.1           
## [51] stringi_1.2.2            yaml_2.1.19             
## [53] edgeR_3.22.2             zlibbioc_1.26.0         
## [55] Rtsne_0.13               rhdf5_2.24.0            
## [57] plyr_1.8.4               grid_3.5.0              
## [59] blob_1.1.1               promises_1.0.1          
## [61] shinydashboard_0.7.0     lattice_0.20-35         
## [63] cowplot_0.9.2            Biostrings_2.48.0       
## [65] locfit_1.5-9.1           pillar_1.2.3            
## [67] igraph_1.2.1             rjson_0.2.19            
## [69] reshape2_1.4.3           biomaRt_2.36.1          
## [71] XML_3.98-1.11            glue_1.2.0              
## [73] evaluate_0.10.1          data.table_1.11.2       
## [75] httpuv_1.4.3             gtable_0.2.0            
## [77] purrr_0.2.4              assertthat_0.2.0        
## [79] xfun_0.1                 mime_0.5                
## [81] xtable_1.8-2             later_0.7.2             
## [83] viridisLite_0.3.0        pheatmap_1.0.10         
## [85] tibble_1.4.2             GenomicAlignments_1.16.0
## [87] beeswarm_0.2.3           memoise_1.1.0           
## [89] tximport_1.8.0           bindrcpp_0.2.2          
## [91] statmod_1.4.30
```
