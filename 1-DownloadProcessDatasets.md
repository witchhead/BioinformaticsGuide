R Notebook
================

# Basics of Single cell RNA-seq

## Bulk RNA-Seq

- estimate average expression level for each gene across a population of
  cells

- use to characterize expression signatures between two groups or
  conditions ex) healthy / diseased

- use to find and annotate new genes, gene isoforms, and other
  transcripts

## Single Cell RNA-Seq

- allows comparison of trascriptomes in individual cells

- assess transcriptional similarity and difference between cells

- can find rare cell populations

- find developmental relationships between heterogeneous cell states

# Key terms

## UMI (Unique Molecular Identifiers)

molecular tags that can be applied to detect and quantify the unique
transcripts

## Features = genes

## Barcodes

Single cell sequencing experiments use short “DNA barcode” tags to
identify reads that originate from the **same cell**

## Count Matrix / Feature-Barcode Matrix / Gene-Barcode Matrix

A Matrix of counts representing the number of unique observations of
each feature within each cell barcode

## Doublets

**Two cells** are encapsulated into one reaction volume -\> should take
them out

# Packages for single-cell RNA Seq Analysis

## R packages

- SingleCellExperiment

- Seurat

- Monocle3

- Scater

## Python

- Scanpy

# Download Gene expression data from NCBI GEO

# Manipulate Gene experession data with dplyr

## Load Libraries

``` r
library(dplyr)
library(tidyverse)
library(GEOquery)
```

## Read the data

``` r
data <- read.csv(".gitignore/GSE183947_fpkm.csv")
```

## Get the metadata using GEO Query package

Geoquery practice =
<https://www.bioconductor.org/packages/devel/bioc/vignettes/GEOquery/inst/doc/GEOquery.html>

``` r
gse <- getGEO(GEO = "GSE183947", GSEMatrix = TRUE)
```

    ## Found 1 file(s)

    ## GSE183947_series_matrix.txt.gz

``` r
gse
```

    ## $GSE183947_series_matrix.txt.gz
    ## ExpressionSet (storageMode: lockedEnvironment)
    ## assayData: 0 features, 60 samples 
    ##   element names: exprs 
    ## protocolData: none
    ## phenoData
    ##   sampleNames: GSM5574685 GSM5574686 ... GSM5574744 (60 total)
    ##   varLabels: title geo_accession ... tissue:ch1 (41 total)
    ##   varMetadata: labelDescription
    ## featureData: none
    ## experimentData: use 'experimentData(object)'
    ##   pubMedIds: 35046993 
    ## Annotation: GPL11154

Get Phenotype data

``` r
metadata <- pData(phenoData(gse[[1]]))
head(metadata)
```

    ##                 title geo_accession                status submission_date
    ## GSM5574685 tumor rep1    GSM5574685 Public on Sep 15 2021     Sep 11 2021
    ## GSM5574686 tumor rep2    GSM5574686 Public on Sep 15 2021     Sep 11 2021
    ## GSM5574687 tumor rep3    GSM5574687 Public on Sep 15 2021     Sep 11 2021
    ## GSM5574688 tumor rep4    GSM5574688 Public on Sep 15 2021     Sep 11 2021
    ## GSM5574689 tumor rep5    GSM5574689 Public on Sep 15 2021     Sep 11 2021
    ## GSM5574690 tumor rep6    GSM5574690 Public on Sep 15 2021     Sep 11 2021
    ##            last_update_date type channel_count source_name_ch1 organism_ch1
    ## GSM5574685      Sep 15 2021  SRA             1          breast Homo sapiens
    ## GSM5574686      Sep 15 2021  SRA             1          breast Homo sapiens
    ## GSM5574687      Sep 15 2021  SRA             1          breast Homo sapiens
    ## GSM5574688      Sep 15 2021  SRA             1          breast Homo sapiens
    ## GSM5574689      Sep 15 2021  SRA             1          breast Homo sapiens
    ## GSM5574690      Sep 15 2021  SRA             1          breast Homo sapiens
    ##             characteristics_ch1 characteristics_ch1.1 characteristics_ch1.2
    ## GSM5574685 tissue: breast tumor       metastasis: yes         donor: 102548
    ## GSM5574686 tissue: breast tumor       metastasis: yes         donor: 104338
    ## GSM5574687 tissue: breast tumor       metastasis: yes         donor: 105094
    ## GSM5574688 tissue: breast tumor        metastasis: no         donor: 109745
    ## GSM5574689 tissue: breast tumor        metastasis: no        donor: 1906415
    ## GSM5574690 tissue: breast tumor       metastasis: yes        donor: 1912627
    ##            molecule_ch1
    ## GSM5574685    total RNA
    ## GSM5574686    total RNA
    ## GSM5574687    total RNA
    ## GSM5574688    total RNA
    ## GSM5574689    total RNA
    ## GSM5574690    total RNA
    ##                                                                                                            extract_protocol_ch1
    ## GSM5574685 Total RNA was isolated and purified using TRIzol (Life, cat.265709, CA, USA) following the manufacturer's procedure.
    ## GSM5574686 Total RNA was isolated and purified using TRIzol (Life, cat.265709, CA, USA) following the manufacturer's procedure.
    ## GSM5574687 Total RNA was isolated and purified using TRIzol (Life, cat.265709, CA, USA) following the manufacturer's procedure.
    ## GSM5574688 Total RNA was isolated and purified using TRIzol (Life, cat.265709, CA, USA) following the manufacturer's procedure.
    ## GSM5574689 Total RNA was isolated and purified using TRIzol (Life, cat.265709, CA, USA) following the manufacturer's procedure.
    ## GSM5574690 Total RNA was isolated and purified using TRIzol (Life, cat.265709, CA, USA) following the manufacturer's procedure.
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           extract_protocol_ch1.1
    ## GSM5574685 After the quality inspection of Agilent 2100 Bioanalyzer (Agilent, cat.G2939AA, CA, USA) and NanoPhotometer (Implen, cat.N60, Munich, Germany), mRNA with poly(A) is purified from 1μg total RNA using VAHTS mRNA Capture Beads with Oligo (dT) (Vazyme, cat.N401-01, Nanjing, China) through two rounds of purification. Subsequently, mRNA fragment was interrupted using VAHTS Universal V6 RNA-seq Library Prep Kit (Vazyme, cat.NR604, Nanjing, China) under 94℃ 8min and reversed transcription into cDNA which would use to synthesise U-labeled second-stranded DNAs. An A-base was added to the blunt ends of each strand to ligase the indexed adapters which contains a T-base at the tail end. After UDG enzyme treatment of the U-labeled double-strand DNA, size selection was performed with VAHTS DNA Clean Beads (Vazyme, cat.N411, Nanjing, China). Then the ligated products are amplified with PCR by the following conditions: initial denaturation at 98℃ for 5 min; 12-17 cycles of denaturation at 98℃ for 10 sec, annealing at 60℃ for 30 sec, and extension at 72℃ for 30 sec; final extension at 72℃ for 5 min. The average insert size of cDNA library was 280±80 bp. After purification by VAHTS DNA Clean Beads (Vazyme, cat.N411-02, Nanjing, China), quality control of concentration and fragment size is performed by Agilent 2100 Bioanalyzer (Agilent, cat.G2939AA, CA, USA) and Qubit assay tubes (Life, cat. 1604220, CA, USA).
    ## GSM5574686 After the quality inspection of Agilent 2100 Bioanalyzer (Agilent, cat.G2939AA, CA, USA) and NanoPhotometer (Implen, cat.N60, Munich, Germany), mRNA with poly(A) is purified from 1μg total RNA using VAHTS mRNA Capture Beads with Oligo (dT) (Vazyme, cat.N401-01, Nanjing, China) through two rounds of purification. Subsequently, mRNA fragment was interrupted using VAHTS Universal V6 RNA-seq Library Prep Kit (Vazyme, cat.NR604, Nanjing, China) under 94℃ 8min and reversed transcription into cDNA which would use to synthesise U-labeled second-stranded DNAs. An A-base was added to the blunt ends of each strand to ligase the indexed adapters which contains a T-base at the tail end. After UDG enzyme treatment of the U-labeled double-strand DNA, size selection was performed with VAHTS DNA Clean Beads (Vazyme, cat.N411, Nanjing, China). Then the ligated products are amplified with PCR by the following conditions: initial denaturation at 98℃ for 5 min; 12-17 cycles of denaturation at 98℃ for 10 sec, annealing at 60℃ for 30 sec, and extension at 72℃ for 30 sec; final extension at 72℃ for 5 min. The average insert size of cDNA library was 280±80 bp. After purification by VAHTS DNA Clean Beads (Vazyme, cat.N411-02, Nanjing, China), quality control of concentration and fragment size is performed by Agilent 2100 Bioanalyzer (Agilent, cat.G2939AA, CA, USA) and Qubit assay tubes (Life, cat. 1604220, CA, USA).
    ## GSM5574687 After the quality inspection of Agilent 2100 Bioanalyzer (Agilent, cat.G2939AA, CA, USA) and NanoPhotometer (Implen, cat.N60, Munich, Germany), mRNA with poly(A) is purified from 1μg total RNA using VAHTS mRNA Capture Beads with Oligo (dT) (Vazyme, cat.N401-01, Nanjing, China) through two rounds of purification. Subsequently, mRNA fragment was interrupted using VAHTS Universal V6 RNA-seq Library Prep Kit (Vazyme, cat.NR604, Nanjing, China) under 94℃ 8min and reversed transcription into cDNA which would use to synthesise U-labeled second-stranded DNAs. An A-base was added to the blunt ends of each strand to ligase the indexed adapters which contains a T-base at the tail end. After UDG enzyme treatment of the U-labeled double-strand DNA, size selection was performed with VAHTS DNA Clean Beads (Vazyme, cat.N411, Nanjing, China). Then the ligated products are amplified with PCR by the following conditions: initial denaturation at 98℃ for 5 min; 12-17 cycles of denaturation at 98℃ for 10 sec, annealing at 60℃ for 30 sec, and extension at 72℃ for 30 sec; final extension at 72℃ for 5 min. The average insert size of cDNA library was 280±80 bp. After purification by VAHTS DNA Clean Beads (Vazyme, cat.N411-02, Nanjing, China), quality control of concentration and fragment size is performed by Agilent 2100 Bioanalyzer (Agilent, cat.G2939AA, CA, USA) and Qubit assay tubes (Life, cat. 1604220, CA, USA).
    ## GSM5574688 After the quality inspection of Agilent 2100 Bioanalyzer (Agilent, cat.G2939AA, CA, USA) and NanoPhotometer (Implen, cat.N60, Munich, Germany), mRNA with poly(A) is purified from 1μg total RNA using VAHTS mRNA Capture Beads with Oligo (dT) (Vazyme, cat.N401-01, Nanjing, China) through two rounds of purification. Subsequently, mRNA fragment was interrupted using VAHTS Universal V6 RNA-seq Library Prep Kit (Vazyme, cat.NR604, Nanjing, China) under 94℃ 8min and reversed transcription into cDNA which would use to synthesise U-labeled second-stranded DNAs. An A-base was added to the blunt ends of each strand to ligase the indexed adapters which contains a T-base at the tail end. After UDG enzyme treatment of the U-labeled double-strand DNA, size selection was performed with VAHTS DNA Clean Beads (Vazyme, cat.N411, Nanjing, China). Then the ligated products are amplified with PCR by the following conditions: initial denaturation at 98℃ for 5 min; 12-17 cycles of denaturation at 98℃ for 10 sec, annealing at 60℃ for 30 sec, and extension at 72℃ for 30 sec; final extension at 72℃ for 5 min. The average insert size of cDNA library was 280±80 bp. After purification by VAHTS DNA Clean Beads (Vazyme, cat.N411-02, Nanjing, China), quality control of concentration and fragment size is performed by Agilent 2100 Bioanalyzer (Agilent, cat.G2939AA, CA, USA) and Qubit assay tubes (Life, cat. 1604220, CA, USA).
    ## GSM5574689 After the quality inspection of Agilent 2100 Bioanalyzer (Agilent, cat.G2939AA, CA, USA) and NanoPhotometer (Implen, cat.N60, Munich, Germany), mRNA with poly(A) is purified from 1μg total RNA using VAHTS mRNA Capture Beads with Oligo (dT) (Vazyme, cat.N401-01, Nanjing, China) through two rounds of purification. Subsequently, mRNA fragment was interrupted using VAHTS Universal V6 RNA-seq Library Prep Kit (Vazyme, cat.NR604, Nanjing, China) under 94℃ 8min and reversed transcription into cDNA which would use to synthesise U-labeled second-stranded DNAs. An A-base was added to the blunt ends of each strand to ligase the indexed adapters which contains a T-base at the tail end. After UDG enzyme treatment of the U-labeled double-strand DNA, size selection was performed with VAHTS DNA Clean Beads (Vazyme, cat.N411, Nanjing, China). Then the ligated products are amplified with PCR by the following conditions: initial denaturation at 98℃ for 5 min; 12-17 cycles of denaturation at 98℃ for 10 sec, annealing at 60℃ for 30 sec, and extension at 72℃ for 30 sec; final extension at 72℃ for 5 min. The average insert size of cDNA library was 280±80 bp. After purification by VAHTS DNA Clean Beads (Vazyme, cat.N411-02, Nanjing, China), quality control of concentration and fragment size is performed by Agilent 2100 Bioanalyzer (Agilent, cat.G2939AA, CA, USA) and Qubit assay tubes (Life, cat. 1604220, CA, USA).
    ## GSM5574690 After the quality inspection of Agilent 2100 Bioanalyzer (Agilent, cat.G2939AA, CA, USA) and NanoPhotometer (Implen, cat.N60, Munich, Germany), mRNA with poly(A) is purified from 1μg total RNA using VAHTS mRNA Capture Beads with Oligo (dT) (Vazyme, cat.N401-01, Nanjing, China) through two rounds of purification. Subsequently, mRNA fragment was interrupted using VAHTS Universal V6 RNA-seq Library Prep Kit (Vazyme, cat.NR604, Nanjing, China) under 94℃ 8min and reversed transcription into cDNA which would use to synthesise U-labeled second-stranded DNAs. An A-base was added to the blunt ends of each strand to ligase the indexed adapters which contains a T-base at the tail end. After UDG enzyme treatment of the U-labeled double-strand DNA, size selection was performed with VAHTS DNA Clean Beads (Vazyme, cat.N411, Nanjing, China). Then the ligated products are amplified with PCR by the following conditions: initial denaturation at 98℃ for 5 min; 12-17 cycles of denaturation at 98℃ for 10 sec, annealing at 60℃ for 30 sec, and extension at 72℃ for 30 sec; final extension at 72℃ for 5 min. The average insert size of cDNA library was 280±80 bp. After purification by VAHTS DNA Clean Beads (Vazyme, cat.N411-02, Nanjing, China), quality control of concentration and fragment size is performed by Agilent 2100 Bioanalyzer (Agilent, cat.G2939AA, CA, USA) and Qubit assay tubes (Life, cat. 1604220, CA, USA).
    ##            taxid_ch1 description
    ## GSM5574685      9606   CA.102548
    ## GSM5574686      9606   CA.104338
    ## GSM5574687      9606   CA.105094
    ## GSM5574688      9606   CA.109745
    ## GSM5574689      9606  CA.1906415
    ## GSM5574690      9606  CA.1912627
    ##                                                                                                                                                                                                                                                               data_processing
    ## GSM5574685 In order to remove technical sequences, including adapters, polymerase chain reaction (PCR) primers, or fragments thereof, and quality of bases lower than 20, pass filter data of fastq format were processed by Cutadapt (V1.9.1) to be high quality clean data.
    ## GSM5574686 In order to remove technical sequences, including adapters, polymerase chain reaction (PCR) primers, or fragments thereof, and quality of bases lower than 20, pass filter data of fastq format were processed by Cutadapt (V1.9.1) to be high quality clean data.
    ## GSM5574687 In order to remove technical sequences, including adapters, polymerase chain reaction (PCR) primers, or fragments thereof, and quality of bases lower than 20, pass filter data of fastq format were processed by Cutadapt (V1.9.1) to be high quality clean data.
    ## GSM5574688 In order to remove technical sequences, including adapters, polymerase chain reaction (PCR) primers, or fragments thereof, and quality of bases lower than 20, pass filter data of fastq format were processed by Cutadapt (V1.9.1) to be high quality clean data.
    ## GSM5574689 In order to remove technical sequences, including adapters, polymerase chain reaction (PCR) primers, or fragments thereof, and quality of bases lower than 20, pass filter data of fastq format were processed by Cutadapt (V1.9.1) to be high quality clean data.
    ## GSM5574690 In order to remove technical sequences, including adapters, polymerase chain reaction (PCR) primers, or fragments thereof, and quality of bases lower than 20, pass filter data of fastq format were processed by Cutadapt (V1.9.1) to be high quality clean data.
    ##                                                                                                                                                                                                                                                                                  data_processing.1
    ## GSM5574685 Firstly, reference genome sequences and gene model annotation files of relative species were downloaded from genome ENSEMBL website, Tophat2 (v2.1.0) was used to index reference genome sequence. Finally, clean data were aligned to reference genome via software Tophat2 (v2.1.02).
    ## GSM5574686 Firstly, reference genome sequences and gene model annotation files of relative species were downloaded from genome ENSEMBL website, Tophat2 (v2.1.0) was used to index reference genome sequence. Finally, clean data were aligned to reference genome via software Tophat2 (v2.1.02).
    ## GSM5574687 Firstly, reference genome sequences and gene model annotation files of relative species were downloaded from genome ENSEMBL website, Tophat2 (v2.1.0) was used to index reference genome sequence. Finally, clean data were aligned to reference genome via software Tophat2 (v2.1.02).
    ## GSM5574688 Firstly, reference genome sequences and gene model annotation files of relative species were downloaded from genome ENSEMBL website, Tophat2 (v2.1.0) was used to index reference genome sequence. Finally, clean data were aligned to reference genome via software Tophat2 (v2.1.02).
    ## GSM5574689 Firstly, reference genome sequences and gene model annotation files of relative species were downloaded from genome ENSEMBL website, Tophat2 (v2.1.0) was used to index reference genome sequence. Finally, clean data were aligned to reference genome via software Tophat2 (v2.1.02).
    ## GSM5574690 Firstly, reference genome sequences and gene model annotation files of relative species were downloaded from genome ENSEMBL website, Tophat2 (v2.1.0) was used to index reference genome sequence. Finally, clean data were aligned to reference genome via software Tophat2 (v2.1.02).
    ##                                                                                                                                                                                                                                                 data_processing.2
    ## GSM5574685 In the beginning transcripts in fasta format are converted from known gff annotation file and indexed properly. Then, with the file as a reference gene file, RSEM (v1.3.1) estimated gene and isoform expression levels from the pair-end clean data.
    ## GSM5574686 In the beginning transcripts in fasta format are converted from known gff annotation file and indexed properly. Then, with the file as a reference gene file, RSEM (v1.3.1) estimated gene and isoform expression levels from the pair-end clean data.
    ## GSM5574687 In the beginning transcripts in fasta format are converted from known gff annotation file and indexed properly. Then, with the file as a reference gene file, RSEM (v1.3.1) estimated gene and isoform expression levels from the pair-end clean data.
    ## GSM5574688 In the beginning transcripts in fasta format are converted from known gff annotation file and indexed properly. Then, with the file as a reference gene file, RSEM (v1.3.1) estimated gene and isoform expression levels from the pair-end clean data.
    ## GSM5574689 In the beginning transcripts in fasta format are converted from known gff annotation file and indexed properly. Then, with the file as a reference gene file, RSEM (v1.3.1) estimated gene and isoform expression levels from the pair-end clean data.
    ## GSM5574690 In the beginning transcripts in fasta format are converted from known gff annotation file and indexed properly. Then, with the file as a reference gene file, RSEM (v1.3.1) estimated gene and isoform expression levels from the pair-end clean data.
    ##                              data_processing.3
    ## GSM5574685 Genome_build: Homo sapiens (GRCh37)
    ## GSM5574686 Genome_build: Homo sapiens (GRCh37)
    ## GSM5574687 Genome_build: Homo sapiens (GRCh37)
    ## GSM5574688 Genome_build: Homo sapiens (GRCh37)
    ## GSM5574689 Genome_build: Homo sapiens (GRCh37)
    ## GSM5574690 Genome_build: Homo sapiens (GRCh37)
    ##                                                                                                            data_processing.4
    ## GSM5574685 Supplementary_files_format_and_content: fpkm.csv: Comma-separated text file includes FPKM values for each Sample.
    ## GSM5574686 Supplementary_files_format_and_content: fpkm.csv: Comma-separated text file includes FPKM values for each Sample.
    ## GSM5574687 Supplementary_files_format_and_content: fpkm.csv: Comma-separated text file includes FPKM values for each Sample.
    ## GSM5574688 Supplementary_files_format_and_content: fpkm.csv: Comma-separated text file includes FPKM values for each Sample.
    ## GSM5574689 Supplementary_files_format_and_content: fpkm.csv: Comma-separated text file includes FPKM values for each Sample.
    ## GSM5574690 Supplementary_files_format_and_content: fpkm.csv: Comma-separated text file includes FPKM values for each Sample.
    ##            platform_id contact_name
    ## GSM5574685    GPL11154  Jianhua,,Xu
    ## GSM5574686    GPL11154  Jianhua,,Xu
    ## GSM5574687    GPL11154  Jianhua,,Xu
    ## GSM5574688    GPL11154  Jianhua,,Xu
    ## GSM5574689    GPL11154  Jianhua,,Xu
    ## GSM5574690    GPL11154  Jianhua,,Xu
    ##                                                                                        contact_institute
    ## GSM5574685 Department of Laboratory Medicine Shunde Hospital of Guangzhou University of Chinese Medicine
    ## GSM5574686 Department of Laboratory Medicine Shunde Hospital of Guangzhou University of Chinese Medicine
    ## GSM5574687 Department of Laboratory Medicine Shunde Hospital of Guangzhou University of Chinese Medicine
    ## GSM5574688 Department of Laboratory Medicine Shunde Hospital of Guangzhou University of Chinese Medicine
    ## GSM5574689 Department of Laboratory Medicine Shunde Hospital of Guangzhou University of Chinese Medicine
    ## GSM5574690 Department of Laboratory Medicine Shunde Hospital of Guangzhou University of Chinese Medicine
    ##                                                                  contact_address
    ## GSM5574685 No. 12, Jinsha Avenue, Shunde District, Foshan, Guangdong, P.R. China
    ## GSM5574686 No. 12, Jinsha Avenue, Shunde District, Foshan, Guangdong, P.R. China
    ## GSM5574687 No. 12, Jinsha Avenue, Shunde District, Foshan, Guangdong, P.R. China
    ## GSM5574688 No. 12, Jinsha Avenue, Shunde District, Foshan, Guangdong, P.R. China
    ## GSM5574689 No. 12, Jinsha Avenue, Shunde District, Foshan, Guangdong, P.R. China
    ## GSM5574690 No. 12, Jinsha Avenue, Shunde District, Foshan, Guangdong, P.R. China
    ##            contact_city contact_state contact_zip/postal_code contact_country
    ## GSM5574685       Fushan     Guangdong                  528300           China
    ## GSM5574686       Fushan     Guangdong                  528300           China
    ## GSM5574687       Fushan     Guangdong                  528300           China
    ## GSM5574688       Fushan     Guangdong                  528300           China
    ## GSM5574689       Fushan     Guangdong                  528300           China
    ## GSM5574690       Fushan     Guangdong                  528300           China
    ##            data_row_count    instrument_model library_selection library_source
    ## GSM5574685              0 Illumina HiSeq 2000              cDNA transcriptomic
    ## GSM5574686              0 Illumina HiSeq 2000              cDNA transcriptomic
    ## GSM5574687              0 Illumina HiSeq 2000              cDNA transcriptomic
    ## GSM5574688              0 Illumina HiSeq 2000              cDNA transcriptomic
    ## GSM5574689              0 Illumina HiSeq 2000              cDNA transcriptomic
    ## GSM5574690              0 Illumina HiSeq 2000              cDNA transcriptomic
    ##            library_strategy
    ## GSM5574685          RNA-Seq
    ## GSM5574686          RNA-Seq
    ## GSM5574687          RNA-Seq
    ## GSM5574688          RNA-Seq
    ## GSM5574689          RNA-Seq
    ## GSM5574690          RNA-Seq
    ##                                                                  relation
    ## GSM5574685 BioSample: https://www.ncbi.nlm.nih.gov/biosample/SAMN21395376
    ## GSM5574686 BioSample: https://www.ncbi.nlm.nih.gov/biosample/SAMN21395377
    ## GSM5574687 BioSample: https://www.ncbi.nlm.nih.gov/biosample/SAMN21395378
    ## GSM5574688 BioSample: https://www.ncbi.nlm.nih.gov/biosample/SAMN21395379
    ## GSM5574689 BioSample: https://www.ncbi.nlm.nih.gov/biosample/SAMN21394912
    ## GSM5574690 BioSample: https://www.ncbi.nlm.nih.gov/biosample/SAMN21394913
    ##                                                        relation.1
    ## GSM5574685 SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRX12143676
    ## GSM5574686 SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRX12143617
    ## GSM5574687 SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRX12143618
    ## GSM5574688 SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRX12143619
    ## GSM5574689 SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRX12143620
    ## GSM5574690 SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRX12143621
    ##            supplementary_file_1 donor:ch1 metastasis:ch1   tissue:ch1
    ## GSM5574685                 NONE    102548            yes breast tumor
    ## GSM5574686                 NONE    104338            yes breast tumor
    ## GSM5574687                 NONE    105094            yes breast tumor
    ## GSM5574688                 NONE    109745             no breast tumor
    ## GSM5574689                 NONE   1906415             no breast tumor
    ## GSM5574690                 NONE   1912627            yes breast tumor

## Use dplyr to process data

``` r
metadata_subset <- metadata %>% 
  select(c(1, 10, 11, 17)) %>%
  rename(tissue = characteristics_ch1,
         metastasis = characteristics_ch1.1) %>%
  mutate(tissue = gsub("tissue: ", "", tissue),
         metastasis = gsub("metastasis: ", "", metastasis))
```

## join both datasets

``` r
data_long <- data %>% rename(gene = X) %>%
  gather(key = 'samples', value = 'fpkm', -gene)
data_long <- data_long %>% left_join(metadata_subset, by = c("samples" = "description"))
head(data_long)
```

    ##       gene   samples fpkm      title       tissue metastasis
    ## 1   TSPAN6 CA.102548 0.93 tumor rep1 breast tumor        yes
    ## 2     TNMD CA.102548 0.00 tumor rep1 breast tumor        yes
    ## 3     DPM1 CA.102548 0.00 tumor rep1 breast tumor        yes
    ## 4    SCYL3 CA.102548 5.78 tumor rep1 breast tumor        yes
    ## 5 C1orf112 CA.102548 2.83 tumor rep1 breast tumor        yes
    ## 6      FGR CA.102548 4.80 tumor rep1 breast tumor        yes

## Summarize the data

``` r
data_long %>% filter(gene == "BRCA1" | gene == "BRCA2") %>%
  group_by(gene, tissue) %>%
  summarise(avg_fpkm = mean(fpkm),
            med_fpkm = median(fpkm)) %>%
  arrange(-avg_fpkm)
```

    ## `summarise()` has grouped output by 'gene'. You can override using the
    ## `.groups` argument.

    ## # A tibble: 4 × 4
    ## # Groups:   gene [2]
    ##   gene  tissue               avg_fpkm med_fpkm
    ##   <chr> <chr>                   <dbl>    <dbl>
    ## 1 BRCA1 breast tumor            10.0      6.96
    ## 2 BRCA1 normal breast tissue     7.70     6.45
    ## 3 BRCA2 normal breast tissue     3.05     1.25
    ## 4 BRCA2 breast tumor             2.04     1.6
