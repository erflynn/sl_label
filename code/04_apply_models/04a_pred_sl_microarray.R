require('cmapR')
require('tidyverse')
require('glmnet')
args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]
idx <- args[2]
ds <- "sex"
load(file=sprintf("data/07_model_dat/fit_%s_microarray_%s.RData", prefix, ds))
ds_path <- sprintf("data/microarray/%s/03_gctx/%s_%s.gctx", prefix, prefix, idx)
my_ds <- parse_gctx(ds_path)

# load the XY genes
if (ds=="sex"){
  xy_genes <- read_csv(sprintf("data/microarray/%s/04_sl_input/xy_genes.csv",prefix))
  xy_genes_l <- xy_genes$overlap.genes
}

gene_names <- my_ds@rid
overlap.genes <- intersect(gene_names,xy_genes_l)

my_mat <- my_ds@mat
my_dat <- my_mat[overlap.genes,]


preds <- predict(fit, newx=t(my_dat), s=my.lambda, type="response")
preds2 <- data.frame("id"=rownames(preds), "pred"=preds[,1])
preds2 %>% write_csv(sprintf("data/08_model_out/%s_microarray_%s_%s_sl.csv", prefix, ds, idx))
