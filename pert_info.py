import pandas as pd
from collections import Counter


__author__ = "Hosein Fooladi"
__email__ = "fooladi.hosein@gmail.com"

def print_pert_statistics(pert_info_dir, pert_type='trt_cp'):
	"""
	This function takes the directory of pert_info.txt
	file and perturbation type, and print some information about perturbations.
	
	Input:
		Mandatory:
		-pert_info_dir (str): The directory of pert_info file. E.g.,
		'./Data/pert_info.txt'
		
		Optional:
		-pert_type (str): perturbation type that you want to extract
		information about it. Default='trt_cp'
	"""
	
	assert isinstance(pert_info_dir, str), "The dataset_dir must be a string object"
	assert isinstance(pert_type, str), "The pert_type must be a string object"
	
	pert_info = pd.read_csv(pert_info_dir, sep = '\t')
	
	assert pert_type in pert_info.pert_type.unique(), "pert_type should be in the list of available perturbations"
	
	print("Data Statistics\n")
	print("Number of available perturbations: {}".format(pert_info.shape[0]))
	print("Number of all Touchstone perturbations: {}".format(pert_info.is_touchstone.sum()))
	
	x = pert_info[pert_info.pert_type==pert_type]
	print("Number of available {}: {}".format(pert_type, x.shape[0]))
	print("Number of Touchstone of {} perturbations: {}".format(pert_type, x.is_touchstone.sum()))
	
	
def pert_touchstone(pert_info_dir, pert_type='trt_cp'):
	"""
	This function takes the directory of pert_info.txt
	file and perturbation type, and return a dictionary that 
	maps perturbations to 0 and 1 (whether they are touchstone or not)
	For example:
	{'BRD-A00100033':1, BRD-A00150179:0, ...}
	
	
	Input:
		Mandatory:
		-pert_info_dir (str): The directory of pert_info file. E.g.,
		'./Data/pert_info.txt'
		
		Optional:
		-pert_type (str): perturbation type that you want to extract
		information about it. Default='trt_cp'
		
	Output:
		-pert_dict_id (Dictionary): A dictionary that maps perturbation ids
		to 0 and 1 (whether they are touchstone or not)
		Keys are pert_id (str) and values are 0 or 1.
		For example:
		{'BRD-A00100033':1, BRD-A00150179:0, ...}
		
		-pert_dict_iname (Dictionary): A dictionary that maps perturbation names
		to 0 and 1 (whether they are touchstone or not)
		Keys are pert_iname (str) and values are 0 or 1.
		For example:
		{'BRD-A00100033':1, BRD-A00150179:0, ...}
		
	"""

	assert isinstance(pert_info_dir, str), "The dataset_dir must be a string object"
	assert isinstance(pert_type, str), "The pert_type must be a string object"
	
	pert_info = pd.read_csv(pert_info_dir, sep = '\t')
	
	assert pert_type in pert_info.pert_type.unique(), "pert_type should be in the list of available perturbations"
	
	x = pert_info[pert_info.pert_type==pert_type]
	
	pert_dict_id = dict(zip(x.pert_id, x.is_touchstone))
	pert_dict_iname = dict(zip(x.pert_iname, x.is_touchstone))
	
	return pert_dict_id, pert_dict_iname
	
	


def duplicate_pert_name(pert_info_dir):
	"""
	There are some perturbation name that maps to more than one
	perturbation ID in Lincs dataset. Here, I am going to find those 
	perturbations that have more than one IDs.
	You can find more information about this by following this link:
	https://clue.io/connectopedia/some_perts_have_over_one_brdid
	
	
	Input:
		Mandatory:
		-pert_info_dir (str): The directory of pert_info file. E.g.,
		'./Data/pert_info.txt'
		

		
	Output:
		-duplicate_list (List): A list of pert_inames that have 
		multiple pert_ids.
		
	"""

	assert isinstance(pert_info_dir, str), "The dataset_dir must be a string object"
	
	pert_info = pd.read_csv(pert_info_dir, sep = '\t')
		
	duplicate_list = [name for name, count in Counter(pert_info.pert_iname).items() if count>1]
	print("Number of pert_iname that have multiple pert_ids: {}".format(len(duplicate_list)))
	
	return duplicate_list
	