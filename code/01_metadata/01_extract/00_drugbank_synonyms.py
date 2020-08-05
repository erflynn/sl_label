# drugbank_synonyms.py
# E Flynn
# 2/27/2019
#
# Get drugbank ID, name, synonym, and ATCCC code from the drugbank XML.
#
# TODO:
# - may want to parse other OTC, FDA approved fields


import xml.etree.ElementTree as ET
import json

drugBankFile = "data/00_db_data/DrugBank.xml" # TODO add to utils
tree = ET.parse(drugBankFile) # slowest step
root = tree.getroot()
treeDrugs = root.findall("{http://www.drugbank.ca}drug")

drug_data = [] # empty list

for drug in treeDrugs:
  # initialize variables
	dbID, drugName, drugCas, drugUnii, chebi_id, cat_label, meshid, atc_id = "","","","","","","",""
	syn_list, category_list, group_list = [],[], []
	
	# get the primary ID and name
	drugBankIDs = drug.findall("{http://www.drugbank.ca}drugbank-id")
	for drugBankID in drugBankIDs:
		if drugBankID.get("primary")=="true":
			dbID = drugBankID.text
	nameField = drug.find("{http://www.drugbank.ca}name")
	if nameField is not None:
		drugName = nameField.text

  # grab mapping to CAS, UNII, and CHEBI if present
	cas = drug.find("{http://www.drugbank.ca}cas-number")
	if cas is not None:
		drugCas = cas.text
	unii = drug.find("{http://www.drugbank.ca}unii")
	if unii is not None:
		drugUnii = unii.text

	groups = drug.findall("{http://www.drugbank.ca}groups")[0].findall("{http://www.drugbank.ca}group")
	for group in groups:
		if groups is not None:
			group_list.append(group.text)

	# chebi mapping is in external identifiers
	external_ids = drug.findall("{http://www.drugbank.ca}external-identifiers")[0].findall("{http://www.drugbank.ca}external-identifier")
	for external_id in external_ids:
		resource = external_id.find("{http://www.drugbank.ca}resource")
		if resource is not None:
			if resource.text=="ChEBI":
				chebi_id = external_id.find("{http://www.drugbank.ca}identifier").text
	
  # grab synonyms
	synonyms = drug.findall("{http://www.drugbank.ca}synonyms")[0].findall("{http://www.drugbank.ca}synonym")
	for synonym in synonyms:
		if synonym is not None:
			syn_list.append(synonym.text)

  # grab the category, meshid 
	categories = drug.findall("{http://www.drugbank.ca}categories")[0].findall("{http://www.drugbank.ca}category")
	for category in categories:
		cat=category.find("{http://www.drugbank.ca}category")
		if cat is not None:
			cat_label = cat.text
		mesh = category.find("{http://www.drugbank.ca}mesh-id")
		if mesh is not None:
			meshid = mesh.text
		category_list.append((cat_label, meshid))
		
	# grab ATC code
	atc_code = drug.findall("{http://www.drugbank.ca}atc-codes")[0].find("{http://www.drugbank.ca}atc-code")
	if atc_code is not None:
		atc_id = atc_code.get("code")
  
  # drug info json
	drug_info = {"dbID":dbID, "name":drugName, "group": group_list, "chebi": chebi_id, "cas":drugCas, "unii":drugUnii, "synonyms" :syn_list, "categories":category_list, "ATC":atc_id}
	drug_data.append(drug_info)

# json dump the whole thing
with open("data/00_db_data/drugbank_info.json", 'w') as f:
	drug_str = json.dumps(drug_data)
	f.write(drug_str)
