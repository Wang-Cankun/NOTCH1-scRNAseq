# NOTCH1-scRNAseq

This is the repository for the manuscript: TODO. 

The original paper was published on TODO

If you have any questions or feedback regarding this notebook, please contact Cankun Wang <cankun.wang@osumc.edu>.

## Directory structure

- 10x_code: HPC scripts to run CellRanger workflow
- velocyto_code: HPC scripts to run velocyto to generate velocity matrix in loom format
- 10x_qc_report: web_summary.html output from 10x CellRanger workflow
- velocity_loom: the authors use this folder to put velocity loom files and did not upload to GitHub due to size limit
- rmarkdown: analysis code from raw read counts
  1. Load 10x read counts
  2. Quality control
  3. Dimension reduction
  4. Cell type annotation
- jupyter: analysis that run using python packages
  1. RNA velocity analysis using scvelo

## Author

- [Cankun Wang](https://github.com/Wang-Cankun)

