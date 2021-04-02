ERR_finder.R is the script that searches for enhancer-promoter connections using GTEx v7 eQTL data.

The script requires R packages 'tidyverse' and bioconductor package 'GenomicRanges' for running.

Annotation files required for running the script are included in the folder 'ERR-datafiles' can be found at https://www.dropbox.com/sh/zmlpnru9nk70mlr/AAC05e0PBHDFBy2PxIozqL7za?dl=0:

Mandatory input files:
-gencode_info.txt :- File containing GENCODE gene annotations
-combinedEnhancerAnnotations.txt :- Enhancer annotations from HEDD database (http://zdzlab.einsteinmed.org/1/hedd.php)
-v7_eqtl_data :- Tissue specific eQTL data from GTEx v7 release
