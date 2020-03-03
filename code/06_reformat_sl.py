
# "PATO:0000047"
#  "PATO:0000384" - male
#  "PATO:0000383" - female

import json
import pandas as pd
metadata = pd.read_csv("data/human_metadata2.csv")
metadata_sm = metadata[["acc", "sex_lab", "pred"]]

row_dicts = []
for index,row in metadata_sm.head(500).iterrows():
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
with open("data/human_500_sex_lab.json", 'w') as f:
	row_str = json.dumps(row_dicts)
	f.write(row_str)
