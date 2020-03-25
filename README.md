

Initial repo for imputing sex, cell line, and other metadata labels from refine-bio data and examining the results.

For now, data files are not included for space limits. This analysis focuses on human, mouse, and rat data; in principle it should be applicable to other organisms previously but this has to be examined.

The `code/` directory is setup as follows:
- `01_metadata/` - code for extracting metadata labels from refine-bio [c] , cleaning [p/c], and mapping them [p]
- `02_expression/` - code for extracting microarray and RNA-seq data from compendia [c]
- `03_models/` - code for building models [p/c]
- `04_apply_models/` - code for applying models [p/c]
- `05_analysis/` - preliminary analyses [p]

Within directory, scripts are numbered according to order to run. This is still in progress. For each part, [p] indicates in progress/preliminary and [c] indicates complete portion of the pipeline.

Most scripts are setup as: \
  `Rscript <script_name> <organism> `
  
Many also take as a final argument an idx for chunking/array scripts.

R package dependencies: `tidyverse`, `data.table`, `glmnet`. \
For microarray extraction, data is converted to `gctx` format using the `cmappy` Python package. 
