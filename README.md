

Repo for imputing sex, cell line, and other metadata labels from refine-bio data and examining the results for the paper [Large-scale Labeling and Assessment of Sex Bias in Publicly Available Expression data](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04070-2) by Flynn, Chang, and Altman.

Data files are not included for space limits, but final inferred labels are located [here](https://figshare.com/s/985621c1705043421962) and will be available in refine.bio. This analysis focuses on human and mouse data; in principle it should be applicable to other organisms previously but this has yet to be examined.

The `code/` directory is setup as follows:
- `01_metadata/` - code for extracting metadata labels from refine-bio, cleaning, and mapping 
- `02_expression/` - code for extracting microarray and RNA-seq data from compendia 
- `03_models/` - code for building models 
- `04_apply_models/` - code for applying models 
- `05_analysis/` - preliminary analyses

Within directory, scripts are numbered according to order to run. 
Most scripts are setup as: \
  `Rscript <script_name> <organism> `
  
Many also take as a final argument an idx for chunking/array scripts.

R package dependencies: `tidyverse`, `data.table`, `glmnet`. \
For microarray extraction, data is converted to `gctx` format using the `cmappy` Python package. 

Please contact erflynn -AT- stanford -DOT- edu with questions.

