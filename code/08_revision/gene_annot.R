# ------- CONSTRUCT GENE DATASETS --------- #
library('biomaRt')
library('readxl')
library('tidyverse')

# Location of PAR region from https://www.ncbi.nlm.nih.gov/assembly

# listDatasets(mart) %>% filter(str_detect(dataset, "hsapiens")) # for version
hmart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl") # GRCh38.p13
#PAR#1 	X 	10,001 	2,781,479
#PAR#2 	X 	155,701,383 	156,030,895
#PAR#1 	Y 	10,001 	2,781,479
#PAR#2 	Y 	56,887,903 	57,217,415

mmart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl") # GRCm38.p6
#PAR 	X 	169,969,759 	170,931,299
#PAR 	Y 	90,745,845 	91,644,698


add_par_human <- function(df) {
  df %>% dplyr::mutate(region=case_when(
    chromosome_name=="X" & 
      start_position >=10001 &
      start_position <=	2781479 ~ "PAR1",
    chromosome_name=="X"  &
      start_position >= 155701383 &
      start_position <= 156030895 ~ "PAR2",
    chromosome_name=="Y" &
      start_position >= 10001  &
      start_position <= 2781479 ~ "PAR1",
    chromosome_name=="Y" &
      start_position >= 56887903  &
      start_position <= 57217415 ~ "PAR2",
    TRUE ~ "non-PAR"))
}





add_par_mouse <- function(df){
  df %>% dplyr::mutate(region=case_when(
    chromosome_name=="X" &
      start_position >= 169969759 &
      start_position <= 170931299  ~ "PAR",
    chromosome_name=="Y" &
      start_position >= 90745845 & 
      start_position <= 91644698 ~ "PAR",
    TRUE ~ "non-PAR"))
}

### load a list of mouse + human genes
human_microarray <- getBM(
  attributes=c("ensembl_gene_id", "hgnc_symbol",
               "start_position", "end_position", "chromosome_name"),
  filters="chromosome_name",
  values=c("X", "Y"),
  mart=hmart) %>%
  as_tibble() %>%
  add_par_human()

human_rnaseq <- getBM(
  attributes=c("ensembl_gene_id", "ensembl_transcript_id","hgnc_symbol",
               "start_position", "end_position", "chromosome_name"),
  filters="chromosome_name",
  values=c("X", "Y"),
  mart=hmart) %>%
  as_tibble() %>%
  add_par_human()


mouse_microarray <- getBM(
  attributes=c("ensembl_gene_id","external_gene_name",
               "hsapiens_homolog_associated_gene_name",
               "start_position", "end_position", "chromosome_name"),
  filters="chromosome_name",
  values=c("X", "Y"),
  mart=mart) %>%
  as_tibble() %>%
  add_par_mouse()

mouse_rnaseq <- getBM(
  attributes=c("ensembl_gene_id", "ensembl_transcript_id","external_gene_name",
               "hsapiens_homolog_associated_gene_name",
               "start_position", "end_position", "chromosome_name"),
  filters="chromosome_name",
  values=c("X", "Y"),
  mart=mart) %>%
  as_tibble() %>%
  add_par_mouse()

save(human_microarray, human_rnaseq, mouse_microarray, mouse_rnaseq, 
     file="data/00_reference/genes/xy_genes_annot_par_biomart.RData")

# ---- add XCI information ----- #

# mouse XCI from ...
# 13 genes escape:
#  Yang F, Babak T, Shendure J, Disteche CM. Global survey of escape from X inactivation by RNA-sequencing in 
# mouse. Genome Res. 2010 May;20(5):614-22. doi: 10.1101/gr.103200.109. Epub 2010 Apr 2.
# PMID: 20363980; PMCID: PMC2860163.
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2860163/
yang_escape_mouse <-c(
  "Xist",
  "Mid1",
  "Kdm6a",
  "Ddx3x",
  "Eif2s3x",
  "Shroom4",
  "Kdm5c",
  "Car5b",
  "Bgn",
  "6720401G13Rik",
  "2610029G23Rik",
  "1810030O07Rik",
  "BC022960")

# Berletch JB, Ma W, Yang F, et al. Escape from X inactivation varies in mouse tissues. 
# PLoS Genet. 2015;11(3):e1005079. Published 2015 Mar 18. doi:10.1371/journal.pgen.1005079
berletch <- mouse_genes <- read_xlsx("data/00_reference/genes/pgen.1005079.s009.xlsx", skip=1)
# pull: `Gene Symbol`, `Escapes XCI in`, `Sex bias`


berletch2 <- berletch %>%
  group_by(`Gene Symbol`) %>%
  dplyr::rename(gene_name=`Gene Symbol`, sex_bias=`Sex bias`) %>%
  mutate(num_tissues_escape=length(str_split(`Escapes XCI in`, ",")[[1]]),
         sex_bias=tolower(sex_bias)) %>%
  ungroup() %>%
  mutate(xci_status="escape", reference="Berletch et al. 2015") %>%
  dplyr::select(gene_name, xci_status, num_tissues_escape, sex_bias, reference)

yang <- tibble("gene_name"=yang_escape_mouse)
yang2 <- yang %>% mutate(xci_status="escape", reference="Yang et al. 2010", num_tissues_escape=NA, sex_bias=NA) 
mouse_escape <- berletch2 %>% 
  bind_rows(yang2 %>% dplyr::select(colnames(berletch2))) %>%
  group_by(gene_name) %>%
  summarize(reference=paste(reference, collapse="; "),
            xci_status=unique(xci_status),
            num_tissues_escape=unique(num_tissues_escape),
            sex_bias=unique(sex_bias)) %>%
  mutate(sex_bias=ifelse(sex_bias=="na", NA, sex_bias))

mouse_escape %>% 
  write_csv("data/00_reference/genes/aggreg_mouse_escape.csv")




# human XCI from Tukiainen et al 2017 https://www.nature.com/articles/nature24265#Sec26
tukiainen <- read_xlsx("data/00_reference/genes/Suppl.Table.13.xlsx")

tukiainen2 <- tukiainen %>% 
  dplyr::rename(gene_name=`Gene name`, sex_bias=`Sex-bias in GTEx`, xci_status=`Reported XCI status`,
                single_cell_xci=`XCI in single cells`, tissue_xci=`XCI across tissues`) %>% 
  dplyr::select(gene_name, xci_status, sex_bias, tissue_xci, single_cell_xci) %>%
  mutate(reference="Tukiainen et al 2017", xci_status=tolower(xci_status), 
         tissue_xci=tolower(tissue_xci), single_cell_xci=tolower(single_cell_xci), 
         sex_bias=str_squish(tolower(str_replace_all(sex_bias, "bias|-", "")) )) %>%
  mutate(tissue_xci=str_replace_all(tissue_xci, ",", " -")) %>%
  mutate(across(everything(), ~ifelse(.=="na", NA, .)))

tukiainen2 %>%
  write_csv("data/00_reference/genes/aggreg_human_escape.csv")

# To write out
# gene id, gene, chr, location, region, xci, sex-bias, resource
# <transcript>, <human homolog>

# ---- filter for present genes ---- #


# --- write this out --- #
