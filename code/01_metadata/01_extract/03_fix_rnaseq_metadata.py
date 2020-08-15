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

def parse_obj(my_acc, acc_type, obj_f, attr_f):
	my_f = "https://www.ebi.ac.uk/ena/data/view/%s&display=xml" %(my_acc)
	with urllib.request.urlopen(my_f) as f:
		tree=ET.parse(f)
		root= tree.getroot()
		obj = root.find(acc_type)
		if (obj is None):
			return -1
		accession = obj.attrib["accession"]
		if accession != my_acc:
			print("Error mismatched accessions %s %s" % (my_acc, accession))
			return -1
		title = obj.find("TITLE").text
		experiments = []
		if acc_type == "RUN":
			exp_ref = obj.find("EXPERIMENT_REF")
			if exp_ref is not None:
				experiment = exp_ref.attrib["accession"]
				list_experiments.add(experiment)
				experiments.append(experiment)
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
		obj_f.write("%s\t%s\t%s\t%s\t%s\n" %(my_acc, title, ";".join(runs), ";".join(samples), ";".join(experiments)))
		# find attributes
		#attr_dict = {}
		my_attr = obj.find("%s_ATTRIBUTES" %acc_type).findall("%s_ATTRIBUTE" %acc_type)
		for attr in my_attr:
			tag_text = attr.find("TAG").text
			value_text = attr.find("VALUE").text
			attr_f.write("%s\t%s\t%s\n" %(my_acc, tag_text, value_text))
			#attr_dict[tag_text]= value_text
		#obj_dict["attr"]=attr_dict
		#return obj_dict
		return 1
	return -1


#setup
idx = sys.argv[1]
with open("data/rnaseq_runs_%s.csv" %idx, 'r') as f:
	list_runs = [run.rstrip("\n") for run in f]


def parse_write(acc_type, acc_list):
	obj_f = open("data/sra_out/%s_info_%s.tsv" %(acc_type.lower(), idx), 'w')
	attr_f = open("data/sra_out/%s_attr_%s.tsv" %(acc_type.lower(), idx), 'w')
	# parse the runs 
	for my_acc in acc_list:
		ret_val = parse_obj(my_acc, acc_type, obj_f, attr_f)
		if ret_val==-1:
			missing_acc.add(my_acc)
	attr_f.close()
	obj_f.close()

start_time = time.time()

parse_write("RUN", list_runs)
parse_write("SAMPLE", list_samples)
parse_write("EXPERIMENT", list_experiments)

# write out missing accessions
with open("data/sra_out/missing_acc_%s.tsv" %(idx), 'w') as f:
	f.write("\n".join(missing_acc))

print("--- %s seconds ---" % (time.time() - start_time))


