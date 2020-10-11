import xml.etree.ElementTree as ET
import json
from os import listdir

exp_dict={} # dictionary for dates
list_acc = listdir()
for my_acc in list_acc:
  print(my_acc)
  tree=ET.parse(my_acc)
  root=tree.getroot()
  experiment_list=root.findall("experiment")
  assert(len(experiment_list)==1)
  release_date = experiment_list[0].find("releasedate").text
  update_date = experiment_list[0].find("lastupdatedate").text
  exp_dict[my_acc]= {"release_date":release_date, "update_date":update_date}
  
with open("../ae_date.json", 'w') as f:
	exp_str = json.dumps(exp_dict)
	f.write(exp_str)

