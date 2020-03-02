# grab from expression data frame using a file
from cmapPy.pandasGEXpress.parse import parse
import pandas as pd
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--infile', help='the file with column names')
parser.add_argument('--prefix', help='the prefix for the files')
parser.add_argument('--outfile', help='the prefix of file to write out')

args = parser.parse_args()

prefix=args.prefix
infile=args.infile
out_prefix=args.outfile

sample_list = pd.read_csv(infile)
sample_sm = sample_list[['acc', 'sex', 'f_idx']]
list_idx = sample_sm['f_idx'].unique()
list_idx.sort()

# iterate through

dfs = []
sex_dfs = []
for idx in list_idx:
  idx_df = sample_sm[sample_sm['f_idx']==idx]
  idx_df = idx_df.sort_values(by=['acc'])
  df_cols = idx_df['acc'].tolist()
  df_sex = idx_df[['acc', 'sex']]
  res=parse("data/%s/03_gctx/%s_%s.gctx" %(prefix, prefix, idx), cid=df_cols)
  dfs.append(res.data_df)
  sex_dfs.append(df_sex)

# concat and write it out
my_df = pd.concat(dfs, axis=1)
my_df.to_csv("data/%s/04_sl_input/%s_expr.csv" %(prefix, out_prefix))

# write out the sex labels
sex_df = pd.concat(sex_dfs, axis=0)
sex_df.to_csv("data/%s/04_sl_input/%s_sex_lab.csv" %(prefix,out_prefix))


