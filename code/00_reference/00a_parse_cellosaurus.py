# parse_cellosaurus.py
# E Flynn
# 3/7/2019
#
# Parse the cellosaurus XML into json file.
#  Extract: age, sex, accesion, synonyms
#  also: disease of origin, MESH/GEO accessions

import xml.etree.ElementTree as ET
import json

cellFile = "/Users/eflynn/Downloads/cellosaurus.xml" # 3/28/2020
tree = ET.parse(cellFile)
root = tree.getroot() 

cl_list = root.find("cell-line-list")
cls = cl_list.findall("cell-line") # 118785
cl_dict = {} # dictionary for holding output

# to look at these we need examples, to get the examples, we need the accession

# acc_list = []
# for cl in cls:
#   accs = cl.find("accession-list").findall("accession")
#   primary_acc=""
#   for acc in accs:
#     acc_num = acc.text
#     acc_type = acc.get("type") # primary or secondary
#     if acc_type=="primary":
#       primary_acc=acc_num
#   acc_list.append(primary_acc)
# 
# cl_parent = cls[acc_list.index("CVCL_SH91")]
# cl_marker = cls[acc_list.index("CVCL_7175")]

for cl in cls:
	# initialize variables
	cl_sex, cl_category, cl_age, identifier, mesh_acc = "","","","",""
	primary_acc, secondary_acc, synonyms, dzs, gex_list, species, alleles, derivs, origins = [],[],[],[],[],[],[],[],[]
  
  # get age/sex/accessions
	cl_sex = cl.get("sex")
	cl_category = cl.get("category")
	cl_age = cl.get("age")
	accs = cl.find("accession-list").findall("accession")
	for acc in accs:
		acc_num = acc.text
		acc_type = acc.get("type") # primary or secondary
		if acc_type=="primary":
			primary_acc.append(acc_num)
		else:
			secondary_acc.append(acc_num)
	acc_info = {"primary" : primary_acc, "secondary" : secondary_acc}
	
	# grab identifiers and synonyms
	names = cl.find("name-list").findall("name")
	for name in names:
		name_text = name.text
		name_type = name.get("type") # identifier or synonym
		if name_type == "identifier":
			if identifier != "":
				print("Error multiple identifiers: %s, %s" %(identifier, name_text)) 
			else:
				identifier = name_text
		else:
			synonyms.append(name_text)
	if identifier == "":
		print("Error no identifier for cell line")
		continue
	
	# grab parents/sibs
	deriv_list = cl.find("derived-from")
	if deriv_list is not None:
		derivs=[deriv.get("accession") for deriv in deriv_list.findall("cv-term")]
	origin_list = cl.find("same-origin-as")
	if origin_list is not None:
		origins =[origin.get("accession") for origin in origin_list.findall("cv-term")]
	
	# grab cross-reference
	atcc_acc = []
	xref_list = cl.find("xref-list")
	if xref_list is not None:
		xrefs = xref_list.findall("xref")
		for xref in xrefs:
			database = xref.get("database")
			if database == "ATCC":
				atcc_acc=xref.get("accession")
	
	# grab markers and their sources
	str_list = cl.find("str-list")
	if str_list is not None:
		marker_list = str_list.find("marker-list")
		if marker_list is not None:
			markers = marker_list.findall("marker")
			for marker in markers:
				marker_id = marker.get("id")
				if marker_id=="Amelogenin":
					marker_data_el = marker.find("marker-data-list").findall("marker-data")
					alleles = { }
					for m in marker_data_el:
						my_alleles = m.find("alleles").text
						source_list = m.find("source-list")
						if source_list is not None:
							my_sources = [my_source.text for my_source in source_list.findall("source")]
							ref_list = source_list.findall("reference-list")
							if len(ref_list) > 0 and ref_list is not None:
								my_refs = []
								for refs in ref_list:
									my_refs.append([ref.text for ref in refs.findall("reference")])
							else:
								my_refs=""
						else:
							my_sources=""
							my_refs=""
						alleles[my_alleles]= {"src": my_sources, "ref":my_refs}

	# extract the disease of origin - this may be helpeful
	dz_list = cl.find("disease-list")
	if dz_list is not None:
		dzs =[dz.text for dz in dz_list.findall("cv-term")]
	
	species_list = cl.find("species-list")
	if species_list is not None:
		species =[spec.text for spec in species_list.findall("cv-term")]
	
	# extract database mappings - Gene Expression databases
	xref_l = cl.find("xref-list")
	if xref_l is not None:
		xrefs = xref_l.findall("xref")
		for xref in xrefs:
			cat = xref.get("category")
			if cat == "Gene expression databases":
				gex_acc = xref.get("accession")
				if gex_acc is not None:
					gex_list.append(gex_acc)
			if cat == "Ontologies":
				if xref.get("database")=="MeSH":
					mesh_acc = xref.get("accession")
	
	# format the output
	cl_dict[identifier] = {"dz":dzs, "species":species, "sex": cl_sex, "category": cl_category, "age" : cl_age, \
	 "synonyms": synonyms, "accession":acc_info, "alleles": alleles,  "derivs": derivs, "origins": origins, \
	 "atcc_acc":atcc_acc, "mesh":mesh_acc, "gex": gex_list}


# write out the json
with open("data/00_db_data/cellosaurus_v3.json", 'w') as f:
	cl_str = json.dumps(cl_dict)
	f.write(cl_str)
