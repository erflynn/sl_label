# Code for reformatting sex labels
# "PATO:0000047"
#  "PATO:0000384" - male
#  "PATO:0000383" - female

import json
import pandas as pd
import argparse


# parse out human argument
parser = argparse.ArgumentParser()
parser.add_argument('--prefix', help='the prefix for the files')
parser.add_argument('--datatype', help='the prefix for the files')

args = parser.parse_args()
prefix=args.prefix
data_type=args.datatype

metadata = pd.read_csv("data/02_labeled_data/%s_%s_sl.csv" %(prefix, data_type ))

print(metadata.shape)
metadata_sm = metadata[["id", "pred"]]

row_dicts = []
for index,row in metadata_sm.iterrows():
  sample_accession = row['id']
  #sex = row['sex_lab']
  prob = row['pred']
  if prob < 0.5:
    sex = "PATO:0000383"
    prob = 1-prob
  else:
    sex = "PATO:0000384"
  row_dict = {}
  row_dict['sample_accession'] = sample_accession
  row_dict['attributes']=[{"PATO:0000047" : {"value" : sex, "probability" : prob}}]
  row_dicts.append(row_dict)

# json dump the whole thing
with open("data/02_labeled_data/%s_%s_sex_lab.json" %(prefix, data_type), 'w') as f:
	row_str = json.dumps(row_dicts)
	f.write(row_str)
