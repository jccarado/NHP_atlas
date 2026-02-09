## Integration Pipeline

1. **1_run_sample_processing.sh**  
   SLURM script that calls `1_run_preprocessing.py`.  
   First step in the QC process:
   - Removes ambient RNA from the filtered feature matrix
   - Uses median absolute deviation to remove outlier cells
   - Calls `scDblFinder.py` and `run_doubletdetector.py` to remove doublets

2. **2_run_downstream_parallel.py**  
   - Reads in `adata` objects with classified doublets  
   - Removes cells identified as doublets by *both* tools  
   - Uses the Python wrapper for **SingleR** to annotate coarse immune cell types (T and B cells)  
   - Filters T and B cells based on cluster abundance to ensure effective subsetting

3. **3_anndata_to_rds.sh**  
   Converts `adata` objects to RDS format for downstream analysis in R.

4. **4_run_rira_annotation.py**  
   Uses the **RIRA** R package to annotate T/B cell subtypes using a rhesus macaque atlas for validation.

5. **5_run_vdj_overall.py**  
   Incorporates existing VDJ-pairing metadata into the atlas.

6. **7_run_scvi.sh**  
   SLURM script that calls `run_scvi.py`.  
   **scVI** is used to remove batch effects in T cells and B cells separately.

7. **8_anndata_to_rds_integration.sh**  
   Converts the integrated `adata` object into an RDS object.
