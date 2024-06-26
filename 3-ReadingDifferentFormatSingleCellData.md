R Notebook
================

# Reading Different Formats of Single Cell Data

you can read different input formats into Seurat Object

| Name                                    | Extension |
|-----------------------------------------|-----------|
| 10x hdf5                                | .hdf5     |
| R data format                           | .rds      |
| AnnData Format                          | .h5ad     |
| Loom                                    | .loom     |
| text based market exchange format (MEX) | .mtx      |

## Load the libraries

``` r
library(Seurat)
```

    ## Warning: package 'Seurat' was built under R version 4.3.3

    ## Loading required package: SeuratObject

    ## Warning: package 'SeuratObject' was built under R version 4.3.3

    ## Loading required package: sp

    ## Warning: package 'sp' was built under R version 4.3.3

    ## 
    ## Attaching package: 'SeuratObject'

    ## The following object is masked from 'package:base':
    ## 
    ##     intersect

``` r
library(SeuratDisk)
```

    ## Registered S3 method overwritten by 'SeuratDisk':
    ##   method            from  
    ##   as.sparse.H5Group Seurat

## Read the datasets

``` r
# .RDS format
rds_obj <- readRDS("ependymal_cells.rds")

# hdf5 format

hdf5_obj <-Read10X_h5(filename = "20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5",
                      use.names = TRUE,
                      unique.features = TRUE)
seurat_hdf5 <- CreateSeuratObject(counts = hdf5_obj)

# .mtx files
mtx_obj <- ReadMtx(mtx = "raw_feature_bc_matrix/matrix.mtx.gz",
                   features = "raw_feature_bc_matrix/features.tsv.gz",
                   cells = "raw_feature_bc_matrix/barcodes.tsv.gz")
seurat_mtx <- CreateSeuratObject(counts = mtx_obj)

# .loom format
loom_obj <- Connect(filename = "CryoPancreatic-human-pancreas-SS2.loom", mode = 'r')
seurat_loom <- as.Seurat(loom_obj)

# AnnData format

Convert("adata_SS2_for_download.h5ad", dest = "h5seurat", overwrite = TRUE)
seurat_anndata <- LoadH5Seurat("adata_SS2_for_download.h5seurat")
```

# Download sequencing data from SRA NCBI

fastq files

Still use the GSE183947 dataset :
[link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183947)

click on SRA Run selector on the bottom of the link.

You can save all the srr items by pressing the accession list below
“Download”
