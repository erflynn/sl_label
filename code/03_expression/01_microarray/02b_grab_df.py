# Grab a microarray expression data frame from the .gctx chunks
# based on a list of IDs in a file.
# The file must contain the fields:
#  'acc' (the sample accession) 
#  'f_idx' (the file index) 
#  'idx' (the index of the column within the file)

# ml python/3.6.1
# ml py-pandas/1.0.3_py36 

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

samples = pd.read_csv(infile)
# //TODO FIX THIS
#sample_list = samples[samples['idx']%20000!=2]
sample_sm = samples[['sample_acc', 'f_idx', 'idx']]
list_idx = sample_sm['f_idx'].unique()
list_idx.sort()

dfs = []

for i in list_idx:
  idx_df = sample_sm[sample_sm['f_idx']==i]
  if i==0:
    idx_df2=idx_df[idx_df['rem_idx']!=1]
  else:
    idx_df2=idx_df[idx_df['rem_idx']!=0]
  df_cols = idx_df2['sample_acc'].tolist()
  res=parse("data/03_expression/microarray/%s/03_gctx/compendia_%s.gctx" %(prefix, i), cid=df_cols)
  dfs.append(res.data_df)

# concat and write it out
my_df = pd.concat(dfs, axis=1)
my_df.to_csv("data/04_df/%s_%s_expr.csv" %(prefix, out_prefix))



