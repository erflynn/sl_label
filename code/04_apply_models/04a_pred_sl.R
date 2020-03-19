require('cmapR')
require('tidyverse')
require('glmnet')
args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]
idx <- args[2]
load(file=sprintf("data/%s/04_sl_input/cvfit.RData", prefix))
ds_path <- sprintf("data/%s/03_gctx/%s_%s.gctx", prefix, prefix, idx)
my_ds <- parse_gctx(ds_path)


load(sprintf("../sex_labeling/geo_pipeline/gpl_ref/%s_gene_map.RData", prefix))
xy_genes <- ref_dat %>% 
  filter(chromosome_name %in% c("X", "Y")) %>% 
  select(chromosome_name, ensembl_gene_id) %>% 
  unique()

gene_names <- my_ds@rid
overlap.genes <- intersect(gene_names,xy_genes$ensembl_gene_id)

my_mat <- my_ds@mat
my_dat <- my_mat[overlap.genes,]

my_dat2 <- my_dat 
new_dat <- as.matrix(t(my_dat2))
new_dat <- apply(new_dat, c(1,2),as.numeric)
preds <- predict(cvfit, newx=t(my_dat), s="lambda.min", type="response")
preds2 <- data.frame("id"=rownames(preds), "pred"=preds[,1])
preds2 %>% write_csv(sprintf("data/%s/05_sl_output/%s_%s_sl.csv", prefix, prefix, idx))
