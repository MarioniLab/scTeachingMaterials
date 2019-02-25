# Teaching materials for single-cell RNA-seq data analysis

This repository contains some standard teaching materials for scRNA-seq data analysis practicals.
It currently relies on packages from Bioconductor version 3.7, and will be updated as required.

- `basic/` uses read count data with ERCC spike-ins from the Fluidigm system.
This covers the basic concepts of scRNA-seq data analysis and takes approximately 2.5 hours to complete.
- `droplet/` uses the publicly available Pan T cell dataset from 10X Genomics.
This describes the analysis of droplet-based scRNA-seq data and takes approximately 1 hour to complete.
- `unbatch/` uses read and UMI count data from studies of the haemtopoietic lineage.
This describes how to perform batch correction and takes approximately 15 minutes to complete.

Obviously, timings are dependent on the number of questions from the audience, and should be treated appropriately.

# Instructions for uploaders

Do **not** upload slides here unless you are uploading the raw TeX files. 

Do **not** upload images here unless you are uploading SVG files.

As a general rule, files should only be uploaded if they can be compiled anywhere.
Required binary content should be hosted at https://jmlab-gitlab.cruk.cam.ac.uk/teaching/scTeachingFiles instead.

# Instructions for instructors

To compile the lectures, run `download.sh` to obtain the image files.
It is then straightforward to run `pdflatex` on each set of slides.
The source SVGs for some images are available in the `pics` subdirectory and can be modified as necessary.

To obtain the data files for each practical, run `prepare.Rmd` in the corresponding directory.
The practical itself can be executed with `test.sh`.
This will update the `answers.md` file for a quick comparison to the reference results.

To create a tarball for distribution:

1. Clone this repository into a fresh directory.
2. Execute the top-level `build.sh`.
This will download all necessary data files.
It will also create two versions of each practical, one with the answers and one without.
3. Delete any directories that will not be used in the workshop.

It is **strongly** recommended that you try compiling each `answers.Rmd` file to ensure that the code is up to date!

# Citation

If you use these materials, please cite their ancestral source:

> Lun ATL, McCarthy DJ and Marioni JC (2016). “A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor.” F1000Res., 5, pp. 2122.

Additional information can be found in the Bioconductor workflow, from which these workshops are most recently derived:

https://bioconductor.org/packages/release/workflows/html/simpleSingleCell.html
