# Do this with ENA data


# for a given SRR (run),
#   - grab all SRS
#   - grab all SRX
#
# fields to extract:
# - attributes
import urllib.request
import xml.etree.ElementTree as ET
import json
import sys
import time


missing_acc = set()
list_samples = set()
list_experiments = set()

def parse_obj(my_acc, acc_type):
	my_f = "https://www.ebi.ac.uk/ena/data/view/%s&display=xml" %(my_acc)
	with urllib.request.urlopen(my_f) as f:
		tree=ET.parse(f)
		root= tree.getroot()
		obj = root.find(acc_type)
		obj_dict = {}
		if (obj is None):
			print(my_acc)
			missing_acc.append(my_acc)
			return []
		accession = obj.attrib["accession"]
		if accession != my_acc:
			print("Error mismatched accessions %s %s" % (my_acc, accession))
			return []
		title = obj.find("TITLE").text
		obj_dict["title"]=title
		if acc_type == "RUN":
			exp_ref = obj.find("EXPERIMENT_REF")
			if exp_ref is not None:
				experiment = exp_ref.attrib["accession"]
				list_experiments.add(experiment)
				experiments = []
			else: experiments = []
		else:
			experiments = []
		# find friends!
		samples = []
		experiments = []
		runs = []
		obj_links = [r.find("XREF_LINK") for r in obj.find("%s_LINKS" %acc_type).findall("%s_LINK" %acc_type)]
		for obj_link in obj_links:
			link_type = obj_link.find("DB").text # TODO, double check multiples
			link_id = obj_link.find("ID").text
			if link_type=="ENA-SAMPLE":
				samples.append(link_id)
				list_samples.add(link_id)
			elif link_type=="ENA-RUN":
				runs.append(link_id)
			elif link_type=="ENA-EXPERIMENT":
				experiments.append(link_id)
				list_experiments.add(link_id)
			else:
				t = "other link"   
		obj_dict["experiments"]=experiments
		obj_dict["runs"]=runs
		obj_dict["samples"]=samples
		# find attributes
		attr_dict = {}
		my_attr = obj.find("%s_ATTRIBUTES" %acc_type).findall("%s_ATTRIBUTE" %acc_type)
		for attr in my_attr:
			tag_text = attr.find("TAG").text
			value_text = attr.find("VALUE").text
			attr_dict[tag_text]= value_text
		obj_dict["attr"]=attr_dict
		return obj_dict
	return []


#setup
idx = sys.argv[1]
with open("data/rnaseq_runs_%s.csv" %idx, 'r') as f:
	list_runs = [run.rstrip("\n") for run in f]

start_time = time.time()
# parse the runs 
run_dict = {}
for my_acc in list_runs:
	run_dict[my_acc] = parse_obj(my_acc, "RUN")

with open("data/sra_out/run_chunk_%s.json" %(idx), 'w') as f:
	run_str = json.dumps(run_dict)
	f.write(run_str)
# parse the samples
sample_dict = {}
for my_acc in list_samples:
	sample_dict[my_acc] = parse_obj(my_acc, "SAMPLE")

with open("data/sra_out/sample_chunk_%s.json" %(idx), 'w') as f:
	sample_str = json.dumps(sample_dict)
	f.write(sample_str)
# parse the experiments
experiment_dict = {}
for my_acc in list_experiments:
	experiment_dict[my_acc] = parse_obj(my_acc, "EXPERIMENT")

with open("data/sra_out/experiment_chunk_%s.json" %(idx), 'w') as f:
	experiment_str = json.dumps(experiment_dict)
	f.write(experiment_str)

print("--- %s seconds ---" % (time.time() - start_time))




