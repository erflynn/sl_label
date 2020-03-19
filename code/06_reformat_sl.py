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
args = parser.parse_args()
prefix=args.prefix

metadata = pd.read_csv("data/%s_metadata2.csv" %prefix )
print(metadata.shape)
metadata_sm = metadata[["acc", "sex_lab", "pred"]]

row_dicts = []
for index,row in metadata_sm.iterrows():
  sample_accession = row['acc']
  sex = row['sex_lab']
  prob = row['pred']
  if sex=="female":
    sex = "PATO:0000383"
    prob = 1-prob
  else:
    sex = "PATO:0000384"
  row_dict = {}
  row_dict['sample_accession'] = sample_accession
  row_dict['attributes']=[{"PATO:0000047" : {"value" : sex, "probability" : prob}}]
  row_dicts.append(row_dict)

# json dump the whole thing
with open("data/%s_sex_lab.json" %prefix, 'w') as f:
	row_str = json.dumps(row_dicts)
	f.write(row_str)
