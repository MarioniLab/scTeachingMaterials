# Teaching materials for single-cell RNA-seq data analysis

This repository contains some standard teaching materials for scRNA-seq data analysis practicals.
It currently relies on packages from Bioconductor version 3.6, and will be updated as required.

- `basic/` uses read count data with ERCC spike-ins from the Fluidigm system.
This covers the basic concepts of scRNA-seq data analysis and takes approximately 2 hours to complete.
- `droplet/` uses the PBMC dataset from 10X Genomics.
This describes the analysis of droplet-based scRNA-seq data and takes approximately 1 hour to complete.

Obviously, timings are dependent on the number of questions from the audience, and should be treated appropriately.

# Instructions for uploaders

Do **not** upload slides here unless you are uploading the raw TeX files. 

Do **not** upload images here unless you are uploading SVG files.

As a general rule, files should only be uploaded if they can be compiled anywhere.

# Instructions for instructors

If you use these materials, please cite their ancestral source:

> Lun ATL, McCarthy DJ and Marioni JC (2016). “A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor.” F1000Res., 5, pp. 2122.

Additional information can be found in the Bioconductor workflow, from which the workflow is most recently derived:

https://bioconductor.org/packages/3.7/workflows/html/simpleSingleCell.html
