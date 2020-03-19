# Grab a microarray expression data frame from the .gctx chunks
# based on a list of IDs in a file.
# The file must contain the fields:
#  'acc' (the sample accession) 
#  'f_idx' (the file index) 
#  'idx' (the index of the column within the file)

from cmapPy.pandasGEXpress.parse import parse
import pandas as pd
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--infile', help='the file with column names')
parser.add_argument('--prefix', help='the prefix for the files (e.g. human, mouse)')
parser.add_argument('--outfile', help='the prefix of file to write out')

args = parser.parse_args()

prefix=args.prefix
infile=args.infile
out_prefix=args.outfile

sample_list = pd.read_csv(infile)
sample_sm = sample_list[['acc', 'f_idx', 'idx']]
list_idx = sample_sm['f_idx'].unique()
list_idx.sort()

dfs = []

for i in list_idx:
  idx_df = sample_sm[sample_sm['f_idx']==i]
  df_cols = idx_df['acc'].tolist()
  res=parse("data/microarray/%s/03_gctx/%s_%s.gctx" %(prefix, prefix, i), cid=df_cols)
  dfs.append(res.data_df)

# concat and write it out
my_df = pd.concat(dfs, axis=1)
my_df.to_csv("data/microarray/%s/04_sl_input/%s_expr.csv" %(prefix, out_prefix))



