import xml.etree.ElementTree as ET
import json
from os import listdir


sra_dict={} # dictionary for dates
list_acc = listdir()

missing_acc=[]

for my_acc in list_acc:
  #print(my_acc)
  tree=ET.parse(my_acc)
  root= tree.getroot()
  study = root.find("STUDY")
  if (study is None):
    print(my_acc)
    missing_acc.append(my_acc)
    continue
  accession = study.attrib["accession"]
  study_attr = study.find("STUDY_ATTRIBUTES").findall("STUDY_ATTRIBUTE")
  first_public=""
  last_update=""
  for my_attr in study_attr:
    tag_text = my_attr.find("TAG").text
    if tag_text=="ENA-FIRST-PUBLIC":
      first_public = my_attr.find("VALUE").text
    if tag_text=="ENA-LAST-UPDATE":
      last_update = my_attr.find("VALUE").text
  sra_dict[accession]= {"first_public":first_public, "last_update":last_update}
  
with open("../sra_date.json", 'w') as f:
	exp_str = json.dumps(sra_dict)
	f.write(exp_str)
