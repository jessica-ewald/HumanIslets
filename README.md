# HumanIslets

This repository contains all of the Rscripts necessary for:

1. processing and normalizing the bulk omics data
2. processing and normalizing the single-cell omics data
3. formatting the omics data into SQLite and HDF5 files for use by the HumanIslets web-tool
4. estmating cell type proportions via deconvolution analysis of the bulk proteomics data
5. precomputing omics:metadata associations for the Feature Search page on the HumanIslets web-tool
6. performing summary statistics and figures for the HumanIslets v1 web-tool manuscript

R scripts relevant to each of the tasks described above are contained in a folder with the same number. The "set_paths.R" script in the root directory contains local paths for the HumanIslets "home" directory (other.tables.path), SQLite directory (sqlite.path), HDF5 directory (hdf5.path), and XiaLab REST API repository (restapi.path). If this is your first time running scripts from this repository, add your local filepaths to set_paths.R. If you believe that you need access to the XiaLab REST API or the HumanIslets home directory, please contact either Dr. Jeff Xia (McGill University) or Dr. Patrick MacDonald (University of Alberta). 
