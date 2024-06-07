# Misc Notes

## Steps needed for updating cells

### Updating SC tab data

- s1a_gen_sstab_data_v2.py
  - change output zarr file name (generates scZarr_YYMMDD.zarr)
  - Ensure only needed cells are present in the output
- s2_gen_aggr_sstab_h5ad.R
  - change input file names
  - prepare input for s2b_aggr_sstab_zarr.py
- s3_gen_sstab_zarr.py
  - change input and output file names (generates aggedSCdata_YYMMDD.zarr)

### Updating CS tab data

- s1c_beadxcell_zarr.py (generates cellscores_cshl_YYMMDD/)
  - change input (one) and output (two) paths that include a date
- s1d_region_agg.py
  - change input (one) and output (one) path that include a date
  - prepares input for s2c_sorted_counts.py
- s2c_sorted_counts.py (generates s2c_cell_jsons_YYMMDD/)
  - change input (one) and output (one) paths that include a date

### Updating SC tab data dendrogram component

- s2_gen_regtocell.py (generates s2_regtocell_YYMMDD/)
  - create output directory in advance
  - requires cell_scores_cshl_YYMMDD/ generated in s1c_beadxcell_zarr.py
