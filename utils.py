from __future__ import unicode_literals, print_function, division

import pickle
import random
import numpy as np
from tqdm import tqdm
from collections import Counter


__author__ = "Hosein Fooladi"
__email__ = "fooladi.hosein@gmail.com"



def print_statistics(dataset_dir, data=None):
	"""
	This function takes the directory of dataset and
	returns some useful statistics about the data.
	
	Input:
		Mandatory:
		-:param dataset_dir (str): It must be string file that shows the directory of the dataset.
		dataset should be a pickle file. e.g., valid argument is something like this:
		'./Data/level3_trt_cp_landmark.pkl'
		
		Optional:
		-:param data (list): It must be a list with the following format:
		line[0]:(cell_line, drug, drug_type, does, does_type, time, time_type)
		line[1]: 978 or 12328-dimensional Vector(Gene_expression_profile)
		
	Note:
	If you provide the data argument, the function igonres the dataset_dir argument 
	and returns output based on the provided data. Otherwise, it returns output 
	based on dataset_dir.
	"""
	
	assert isinstance(dataset_dir, str), "The dataset_dir must be a string object"
	
	
	print("=================================================================")
	print("Data Loading..")
	
	if data is None: 
		with open(dataset_dir, "rb") as f:
			train = pickle.load(f)
	else:
		assert isinstance(data, list), "The data must be a list object"
		train = data
		
		
	print("Data Statistics\n")
	print("Number of Train Data: {}".format(len(train)))
	
	
	print("Please wait while we are retriving information ...")
	cell_lines = []
	compounds = []
	doses = []
	times = []
	for i in range(len(train)):
		cell_lines.append(train[i][0][0])
		compounds.append(train[i][0][1])
		doses.append(train[i][0][3])
		times.append(train[i][0][5])
		
	print("Number of unique Cell Lines: {}".format(len(set(cell_lines))))
	print("Number of unique Compounds: {}".format(len(set(compounds))))
	print("Number of unique doses: {}".format(len(set(doses))))
	print("Number of unique times: {}".format(len(set(times))))
		
		

def print_most_frequent(dataset_dir, n=3, data=None):
	"""
	This function takes the directory of dataset and integer n
	and returns The n most frequent cell lines, compounds and
	doses in the dataset.
	
	Input:
		Mandatory:
		-:param dataset_dir (str): It must be string file that shows the directory of the dataset.
		dataset should be a pickle file. e.g., valid argument is something like this:
		'./Data/level3_trt_cp_landmark.pkl'
		
		Optional:
		-:param n (int): An integern which determine number of frequent statistics we want 
		to retrieve. Default=3.
		-:param data (list): It must be a list with the following format:
		line[0]:(cell_line, drug, drug_type, does, does_type, time, time_type)
		line[1]: 978 or 12328-dimensional Vector(Gene_expression_profile)
		
	Note:
	If you provide the data argument, the function igonres the dataset_dir argument 
	and returns output based on the provided data. Otherwise, it returns output 
	based on dataset_dir.
	"""
	
	assert isinstance(dataset_dir, str), "The dataset_dir must be a string object"
	assert isinstance(n, int), "The parameter n must be an integer"
	
	
	
	print("=================================================================")
	print("Data Loading..")
	
	if data is None: 
		with open(dataset_dir, "rb") as f:
			train = pickle.load(f)
	else:
		assert isinstance(data, list), "The data must be a list object"
		train = data
		
	
	print("Please wait while we are retriving information ...")
	cell_lines = []
	compounds = []
	doses = []
	for i in tqdm(range(len(train))):
		cell_lines.append(train[i][0][0])
		compounds.append(train[i][0][1])
		doses.append(train[i][0][3])
		
	print("loop finished !!!")
				
	print("Most frequent Cell Lines: {}".format(Counter(cell_lines).most_common(n)))
	print("Most frequent Compounds: {}".format(Counter(compounds).most_common(n)))
	print("Most frequent Doses: {}".format(Counter(doses).most_common(n)))
	
	
	
def cell_line_frequent(dataset_dir, n=3, data=None):
	"""
	This function takes the directory of dataset and integer n,
	and parse the data to keep only the data that belongs to n 
	most frequent cell lines.
	
	Input:
		Mandatory:
		-:param dataset_dir (str): It must be string file that shows the directory of the dataset.
		dataset should be a pickle file. e.g., valid argument is something like this:
		'./Data/level3_trt_cp_landmark.pkl'
		
		Optional:
		-:param n (int): An integern which determine number of frequent statistics we want 
		to retrieve. Default=3.
		-:param data (list): It must be a list with the following format:
		line[0]:(cell_line, drug, drug_type, does, does_type, time, time_type)
		line[1]: 978 or 12328-dimensional Vector(Gene_expression_profile)
		
	Output:
		-:param parse_data (list): A list containing data that belongs to n most frequent
		cell lines.
	
	Note:
	If you provide the data argument, the function igonres the dataset_dir argument 
	and returns output based on the provided data. Otherwise, it returns output 
	based on dataset_dir.
		
	"""
	
	assert isinstance(dataset_dir, str), "The dataset_dir must be a string object"
	assert isinstance(n, int), "The parameter n must be an integer"
	
		
	print("=================================================================")
	print("Data Loading..")
	if data is None: 
		with open(dataset_dir, "rb") as f:
			train = pickle.load(f)
	else:
		assert isinstance(data, list), "The data must be a list object"
		train = data
		
	
	print("Please wait while we are retriving information ...")
	cell_lines = []

	for i in tqdm(range(len(train))):
		cell_lines.append(train[i][0][0])
		
	print("Number of unique Cell Lines: {}".format(len(set(cell_lines))))	
	print("Most frequent Cell Lines: {}".format(Counter(cell_lines).most_common(n)))
	
	assert n < len(set(cell_lines)), "n is out of valid range!"
	
	## List of n most frequent cell lines
	x = list(map(lambda x : x[0], Counter(cell_lines).most_common(n)))
	
	parse_data = [line for line in train if line[0][0] in x] 
	
	return parse_data
	


def cell_line_list(dataset_dir, cells = ['MCF7'], data=None):
	"""
	This function takes the directory of dataset and a list cells,
	and parse the data to keep only the data that belongs to cells list.
	
	Input:
		Mandatory:
		-:param dataset_dir (str): It must be string file that shows the directory of the dataset.
		dataset should be a pickle file. e.g., valid argument is something like this:
		'./Data/level3_trt_cp_landmark.pkl'
		
		Optional:
		-:param cells (list of strings): list of cell lines that we want to keep their data 
		to retrieve. Default=['MCF7']
		-:param data (list): It must be a list with the following format:
		line[0]:(cell_line, drug, drug_type, does, does_type, time, time_type)
		line[1]: 978 or 12328-dimensional Vector(Gene_expression_profile)
		
	Output:
		-:param parse_data (list): A list containing data that belongs to desired list.	

		
	Note:
	If you provide the data argument, the function igonres the dataset_dir argument 
	and returns output based on the provided data. Otherwise, it returns output 
	based on dataset_dir.
	"""
	
	assert isinstance(dataset_dir, str), "The dataset_dir must be a string object"
	assert isinstance(cells, list), "The parameter cells must be a list"
	
		
	print("=================================================================")
	print("Data Loading..")
	if data is None: 
		with open(dataset_dir, "rb") as f:
			train = pickle.load(f)
	else:
		assert isinstance(data, list), "The data must be a list object"
		train = data
		

	print("Number of Train Data: {}".format(len(train)))
	
	parse_data = [line for line in train if line[0][0] in cells] 
	
	print("Number of Data after parsing: {}".format(len(parse_data)))
	return parse_data
	
	
	

def parse_list(dataset_dir, indicator=0, query=['MCF7'], data=None):
	"""
	This function takes the directory of dataset, indicator that indicates
	whether you want to subset the data based on cell line, compound, dose, or time
	and a list which shows what part of the data you want to keep.
	The output will be a list of desired parsed dataset. 
	
	
	Input:
		Mandatory:
		-:param dataset_dir (str): It must be string file that shows the directory of the dataset.
		dataset should be a pickle file. e.g., valid argument is something like this:
		'./Data/level3_trt_cp_landmark.pkl'
		
		Optional:
		-:params indicator (int): it must be an integer from 0 1 2 and 3 that shows whether
		we want to retrieve the data based on cells, compound or dose.
		0: cell_lines   
		1:compounds
		2:doses	
		3:time
		Default=0 (cell_lines)
		-:params query (list): list of cells or compounds or doses that we want to retrieve.
		The list depends on the indicator. If the indicator is 0, you should enter the
		list of desired cell lines and so on. Default=['MCF7']
		-:param data (list): It must be a list with the following format:
		line[0]:(cell_line, drug, drug_type, does, does_type, time, time_type)
		line[1]: 978 or 12328-dimensional Vector(Gene_expression_profile)
	
	Output:
		-:params parse_data (list): A list containing data that belongs to desired list.	
		
	Note:
	If you provide the data argument, the function igonres the dataset_dir argument 
	and returns output based on the provided data. Otherwise, it returns output 
	based on dataset_dir.
	"""
	
	assert isinstance(dataset_dir, str), "The dataset_dir must be a string object"
	assert isinstance(indicator, int), "The indicator must be an int object"
	assert indicator in [0, 1, 2, 3], "You should choose indicator from 0, 1, 2 range"
	assert isinstance(query, list), "The parameter query must be a list"
		
	print("=================================================================")
	print("Data Loading..")
	if data is None: 
		with open(dataset_dir, "rb") as f:
			train = pickle.load(f)
	else:
		assert isinstance(data, list), "The data must be a list object"
		train = data
		
	mapping = {0:0, 1:1, 2:3, 3:5}
	k = mapping[indicator]
	mapping_name = {0:'cell_lines', 1:'compounds', 2:'doses', 3:'time'}
	
	print("Number of Train Data: {}".format(len(train)))
	print("You are parsing the data base on {}". format(mapping_name[indicator]))
	
	parse_data = [line for line in train if line[0][k] in query] 
	
	print("Number of Data after parsing: {}".format(len(parse_data)))
	return parse_data
	

def parse_most_frequent(dataset_dir, indicator=0, n=3, data=None):
	"""
	This function takes the directory of dataset, indicator that indicates
	whether you want to subset the data based on cell line, compound, dose, or time
	and a n which how much frequent items you want to keep.
	The output will be a list of desired parsed dataset. 
	
	Input:
		Mandatory:
		-:param dataset_dir (str): It must be string file that shows the directory of the dataset.
		dataset should be a pickle file. e.g., valid argument is something like this:
		'./Data/level3_trt_cp_landmark.pkl'
		
		Optional:
		-:params indicator (int): it must be an integer from 0 1 2 and 3 that shows whether
		we want to retrieve the data based on cells, compound or dose.
		0: cell_lines   
		1:compounds
		2:doses
		3:time
		Default=0	
		-:params n (int): number of most frequent cells or compounds or doses that we want to retrieve.
		The list depends on the indicator. If the indicator is 0, you should enter the
		number of desired cell lines and so on. Default=3 
		-:param data (list): It must be a list with the following format:
		line[0]:(cell_line, drug, drug_type, does, does_type, time, time_type)
		line[1]: 978 or 12328-dimensional Vector(Gene_expression_profile)
	
	Output:
		-:params parse_data: A list containing data that belongs to desired list.

	Note:
	If you provide the data argument, the function igonres the dataset_dir argument 
	and returns output based on the provided data. Otherwise, it returns output 
	based on dataset_dir.		
	"""
	
	assert isinstance(dataset_dir, str), "The dataset_dir must be a string object"
	assert isinstance(indicator, int), "The indicator must be an int object"
	assert indicator in [0, 1, 2], "You should choose indicator from 0, 1, 2 range"
	assert isinstance(n, int), "The parameter n must be an integer"
		
	print("=================================================================")
	print("Data Loading..")
	if data is None: 
		with open(dataset_dir, "rb") as f:
			train = pickle.load(f)
	else:
		assert isinstance(data, list), "The data must be a list object"
		train = data
		
	mapping = {0:0, 1:1, 2:3, 3:5}
	k = mapping[indicator]
	mapping_name = {0:'cell_lines', 1:'compounds', 2:'doses', 3:'time'}
	
	mylist = []

	for i in tqdm(range(len(train))):
		mylist.append(train[i][0][k])
		
	print("Number of unique {}: {}".format(mapping_name[indicator], len(set(mylist))))	
	print("Most frequent {}: {}".format(mapping_name[indicator], Counter(mylist).most_common(n)))
	
	assert n < len(set(mylist)), "n is out of valid range!"
	
	## List of n most frequent cell lines
	y = list(map(lambda x : x[0], Counter(mylist).most_common(n)))
	
	parse_data = [line for line in train if line[0][k] in y] 
	
	return parse_data




def parse_chunk_frequent(dataset_dir, indicator=0, start=0, end=3, data=None):
	"""
	This function takes the directory of dataset, indicator that indicates
	whether you want to subset the data based on cell line, compound, dose, or time
	and a start and end which shows what chunk of data is desirable.
	E.g., if start=0 and end=3, you are subsetting 3 most frequent data.
	The output will be a list of desired parsed dataset. 
	
	Input:
		Mandatory:
		-:param dataset_dir (str): It must be string file that shows the directory of the dataset.
		dataset should be a pickle file. e.g., valid argument is something like this:
		'./Data/level3_trt_cp_landmark.pkl'
		
		Optional:
		-:params indicator (int): it must be an integer from 0 1 2 and 3 that shows whether
		we want to retrieve the data based on cells, compound or dose.
		0: cell_lines   
		1:compounds
		2:doses
		3:time
		Default=0	
		-:params start: indicates the start of the list you want to subset. Default=0
		-:params end: indicates the end of the list you want to subset. Default=3
		
		-:param data (list): It must be a list with the following format:
		line[0]:(cell_line, drug, drug_type, does, does_type, time, time_type)
		line[1]: 978 or 12328-dimensional Vector(Gene_expression_profile)
	
	Output:
		-:params parse_data: A list containing data that belongs to desired list.	
		
	Note:
	If you provide the data argument, the function igonres the dataset_dir argument 
	and returns output based on the provided data. Otherwise, it returns output 
	based on dataset_dir.
	"""
	
	assert isinstance(dataset_dir, str), "The dataset_dir must be a string object"
	assert isinstance(indicator, int), "The indicator must be an int object"
	assert indicator in [0, 1, 2], "You should choose indicator from 0, 1, 2 range"
	assert isinstance(start, int), "The parameter start must be an integer"
	assert isinstance(end, int), "The parameter end must be an integer"
	assert start <= end, "The start should be less than the end!!"


		
	print("=================================================================")
	print("Data Loading..")
	if data is None: 
		with open(dataset_dir, "rb") as f:
			train = pickle.load(f)
	else:
		assert isinstance(data, list), "The data must be a list object"
		train = data
		
	mapping = {0:0, 1:1, 2:3, 3:5}
	k = mapping[indicator]
	mapping_name = {0:'cell_lines', 1:'compounds', 2:'doses', 3:'time'}
	
	mylist = []

	for i in range(len(train)):
		mylist.append(train[i][0][k])
		
	print("Number of unique {}: {}".format(mapping_name[indicator], len(set(mylist))))	
	
	assert end < len(set(mylist)), "end is out of valid range!"
	
	## List of n most frequent cell lines
	y = list(map(lambda x : x[0], Counter(mylist).most_common()))[start:end]
	
	print("Desired {}: {}".format(mapping_name[indicator], y))

	parse_data = [line for line in train if line[0][k] in y] 
	
	return parse_data
	
	
	
def parse_dose_range(dataset_dir, dose_min=0, dose_max=5, data=None):
	"""
	This function takes the directory of dataset minimum and maximum dose
	and return a list of data that are within the desired range. 
	
	Input:
		Mandatory:
		-:param dataset_dir: It must be string file that shows the directory of the dataset.
		dataset should be a pickle file. e.g., valid argument is something like this:
		'./Data/level3_trt_cp_landmark.pkl'
		
		Optional:
		-:params dose_min (int): minimum dose. Default=0 
		-:params dose_max (int): maximum_dose. Default=5
		
		-:param data (list): It must be a list with the following format:
		line[0]:(cell_line, drug, drug_type, does, does_type, time, time_type)
		line[1]: 978 or 12328-dimensional Vector(Gene_expression_profile)
		
	Output:
		-:params parse_data: A list containing data that belongs to desired list (
		Desired range of doses).

	Note:
	If you provide the data argument, the function igonres the dataset_dir argument 
	and returns output based on the provided data. Otherwise, it returns output 
	based on dataset_dir.		
	"""
	
	assert isinstance(dataset_dir, str), "The dataset_dir must be a string object"
	assert isinstance(dose_min, int), "The parameter dose_min must be an integer"
	assert isinstance(dose_max, int), "The parameter dose_max must be an integer"
	assert dose_min < dose_max , "The minimum dose must be less than the maximum dose !!"

	
		
	print("=================================================================")
	print("Data Loading..")
	if data is None: 
		with open(dataset_dir, "rb") as f:
			train = pickle.load(f)
	else:
		assert isinstance(data, list), "The data must be a list object"
		train = data
		

	print("Number of Train Data: {}".format(len(train)))
	
	parse_data = [line for line in train if line[0][3]>dose_min and line[0][3]<dose_max] 
	
	print("Number of Data after parsing: {}".format(len(parse_data)))
	return parse_data
	
	
	
	
def parse_list_v2(dataset_dir, indicator=0, query=['MCF7'], data=None):
	"""
	This function takes the directory of dataset, indicator that indicates
	whether you want to subset the data based on cell line, compound, dose, time, touchstone,
	clinical phase, MOA or target. Moreover, it takes a list which shows what part of the data you want to keep.
	The output will be a list of desired parsed dataset. 
	
	
	Input:
		Mandatory:
		-:param dataset_dir (str): It must be string file that shows the directory of the dataset.
		dataset should be a pickle file. e.g., valid argument is something like this:
		'./Data/level3_trt_cp_landmark_allinfo.pkl'
		
		The pickle file should be as the following:
		list Format:
		line[0]:(cell_line,
					drug, 
					drug_type, 
					does, 
					does_type, 
					time, 
					time_type,
					touchstone,
					clinical phase,
					moa,
					target)
		line[1]: 978 or 12328-dimensional Vector(Gene_expression_profile)
		
		Optional:
		-:params indicator (int): it must be an integer from 0 1 2 3 4 5 6 7 that shows whether
		we want to retrieve the data based on cells, compound, dose, touchstone, clinical phase, moa or target.
		0: cell_lines   
		1: compounds
		2: doses
		3: time		
		4: touchstone
		5: clinical phase
		6: moa
		7: target
		Default=0 (cell_lines)
		-:params query (list): list of cells or compounds or doses or time or touchstone or clinical phase or MOA or target that we want to retrieve.
		The list depends on the indicator. If the indicator is 0, you should enter the
		list of desired cell lines and so on. Default=['MCF7']
		
		-:param data (list): It must be a list with the following format:
		line[0]:(cell_line, drug, drug_type, does, does_type, time, time_type, 
		touchstone, clinical phase, moa, target)
		line[1]: 978 or 12328-dimensional Vector(Gene_expression_profile)
	
	Output:
		-:params parse_data (list): A list containing data that belongs to desired list.
		
	Note:
	If you provide the data argument, the function igonres the dataset_dir argument 
	and returns output based on the provided data. Otherwise, it returns output 
	based on dataset_dir.
	"""
	
	assert isinstance(dataset_dir, str), "The dataset_dir must be a string object"
	assert isinstance(indicator, int), "The indicator must be an int object"
	assert indicator in [0, 1, 2, 3, 4, 5, 6, 7], "You should choose indicator from 0, 1, 2, 3, 4, 5, 6, 7 range"
	assert isinstance(query, list), "The parameter query must be a list"
		
	print("=================================================================")
	print("Data Loading..")
	if data is None: 
		with open(dataset_dir, "rb") as f:
			train = pickle.load(f)
	else:
		assert isinstance(data, list), "The data must be a list object"
		train = data
		
	mapping = {0:0, 1:1, 2:3, 3:5, 4:7, 5:8, 6:9, 7:10}
	k = mapping[indicator]
	mapping_name = {0:'cell_lines', 1:'compounds', 2:'doses', 3:'time', 4:'tochstone', 5:'clinical_phase', 6:'moa', 7:'target'}
	
	print("Number of Train Data: {}".format(len(train)))
	print("You are parsing the data base on {}". format(mapping_name[indicator]))
	
	parse_data =[]
	if indicator in [0, 1, 2, 3, 4]:
		parse_data = [line for line in train if line[0][k] in query] 
		
	elif indicator in [5, 6, 7]:
		for line in train:
			tmp = line[0][k][0].split('|')
			for a in tmp:
				if a in query:
					parse_data.append(line)
					break
	
	print("Number of Data after parsing: {}".format(len(parse_data)))
	return parse_data
	