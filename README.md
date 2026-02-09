==INTEGRATION PIPELINE===

1_run_sample_processing.sh - SLRUM script that calls 1_run_preprocessing.py. First step in the QC process. Removes ambient RNA present in the filtered feature matrix, and then uses median absolute deviation to remove outlier cells. Lastly, calls the python scripts scDblFinder.py and run_doubletdetector.py to remove doublets for the same dataset. 
2_run_downstream_parallel.py - Reads in the adata objects with classified doublets, and takes the intersection of doublet identified by both tools for removal. The python wrapper for SingleR is used to annotate coarse immune cell types T and B cells.  These T and B cells are filtered based on their abundance in each cluster so as to ensure effective subsetting.
3_anndata_to_rds.sh - Converts adata object to RDS object for analysis in R
4_run_rira_annotation.py - Uses the RIRA R package to annotate T/B cell subtypes from a rhesus macaque atlas for validation
5_run_vdj_overall.py - Incorporates existing VDJ-pairing metadata into the atlas
7_run_scvi.sh - SLURM script that calls the python script "run_scvi.py". scVI is used to remove batch effects in the T cells and B cells separately. 
8_anndata_to_rds_integration.sh - Converts integrated adata object into RDS object
