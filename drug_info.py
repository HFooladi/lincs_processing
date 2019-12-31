'''
Provided are annotations for 6,125 drug and tool compounds 
(2,369 FDA-approved drugs, 1,619 drugs that reached phases 1-3 of clinical development, 
96 compounds that were previously approved but withdrawn from use, and 2,041 preclinical
 or tool compounds). Annotations include compound name, chemical structure, clinical 
 trial status, mechanism of action, protein targets, disease areas, approved indications 
 (where applicable), purity of the purchased sample, and vendor ID.
'''

import pandas as pd
from collections import Counter


__author__ = "Hosein Fooladi"
__email__ = "fooladi.hosein@gmail.com"

def print_drug_statistics(drug_info_dir):
	"""
	This function takes the directory of drug_info.txt
	and print some information about perturbations.
	
	Input:
		Mandatory:
		-drug_info_dir (str): The directory of drug_info file. E.g.,
		'./Data/repurposing_drugs_20180907.txt'
		

	"""
	
	assert isinstance(drug_info_dir, str), "The dataset_dir must be a string object"
	
	drug_info = pd.read_csv(drug_info_dir, sep = '\t', skiprows=9, encoding='latin-1')
	
	
	print("=================================================================")
	print("Data Statistics\n")
	print("Number of available drugs in the datasets: {}".format(drug_info.shape[0]))
	print("Number of Unique mechanism of actions: {}".format(drug_info.moa.unique().shape[0]))
	
	print("Information about Number of drugs in different clinical phases: {}".format(Counter(drug_info.clinical_phase)))
	
	
	
def drug_pert_retrieval(drug_info_dir, pert_info_dir, pert_type='trt_cp'):
	"""
	This function takes the directory of drug_info.txt, pert_info.txt, 
	and pert_type; and return drug information correspondant to
	perturbation information.
	
	
	Input:
		Mandatory:
		-drug_info_dir (str): The directory of drug_info file. E.g.,
		'./Data/repurposing_drugs_20180907.txt'
		-pert_info_dir (str): The directory of pert_info file. E.g.,
		'./Data/pert_info.txt'
		
		Optional:
		-pert_type (str): perturbation type that you want to extract
		information about it. Default='trt_cp'
		
	Output:
		-pert_supp_info(pandas.DataFrame): Pandas dataframe which contains
		supplementary information about touchstone compounds such as 
		mode of actions, tarfets and ...
		pert_list (List): List of strings which indicates the compounds that
		we have additional information for them.
		Example: ['abiraterone', 'ABT-737', 'ABT-751', 'AC-55649',...]
		

	"""
	
	assert isinstance(drug_info_dir, str), "The dataset_dir must be a string object"
	assert isinstance(pert_info_dir, str), "The dataset_dir must be a string object"
	assert isinstance(pert_type, str), "The pert_type must be a string object"
	
	drug_info = pd.read_csv(drug_info_dir, sep = '\t', skiprows=9, encoding='latin-1')
	pert_info = pd.read_csv(pert_info_dir, sep = '\t')
	
	assert pert_type in pert_info.pert_type.unique(), "pert_type should be in the list of available perturbations"
	
	x = pert_info[(pert_info.pert_type=='trt_cp') & (pert_info.is_touchstone)]
	
	print("=================================================================")
	print("Number of Touchstone of {}: {}".format(pert_type, x.shape[0]))
	
	pert_supp_info = drug_info[drug_info.pert_iname.isin(x.pert_iname)]
	pert_list = list(y.pert_iname)
	
	print("Number of {} that we have additional information about them: {}".format(pert_type, len(pert_list)))

	return pert_supp_info, pert_list
