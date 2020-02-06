import argparse
import numpy as np
import pickle
import random
from tqdm import tqdm
from collections import Counter


__author__ = "Hosein Fooladi"
__email__ = "fooladi.hosein@gmail.com"


parser = argparse.ArgumentParser(description='Parsing LINCS')
parser.add_argument('--dataset_dir', type=str, default='Data/level3_trt_cp_landmark.pkl')
parser.add_argument('--data', type=list, default=None)
parser.add_argument('--cells', type=str, nargs='+', default=None)
parser.add_argument('--compounds', type=str, nargs='+', default=None)
parser.add_argument('--doses', type=float, nargs='+', default=None)
parser.add_argument('--times', type=int, nargs='+', default=None)
parser.add_argument('--output_dir', type=str, default='Data/after_parsing.pkl')


flags = parser.parse_args()

assert isinstance(flags.dataset_dir, str), "The dataset_dir must be a string object"

print("=================================================================")
print("Data Loading..")

## If you enter the data (which is a list), it overrides the dataset_dir and ignore it.	
if flags.data is None: 
	with open(flags.dataset_dir, "rb") as f:
		train = pickle.load(f)
else:
	assert isinstance(flags.data, list), "The data must be a list object"
	train = flags.data
	
print("Number of Train Data: {}".format(len(train)))

if flags.cells is None:
	pass
else:
	train = [line for line in train if line[0][0] in flags.cells]
	print("Number of training data after parsing based on cell lines: {}".format(len(train)))

	
if flags.compounds is None:
	pass
else:
	train = [line for line in train if line[0][1] in flags.compounds]
	print("Number of training data after parsing based on compounds: {}".format(len(train)))
	
if flags.doses is None:
	pass
else:
	train = [line for line in train if line[0][3] in flags.doses]
	print("Number of training data after parsing based on doses: {}".format(len(train)))
	
if flags.times is None:
	pass
else:
	train = [line for line in train if line[0][5] in flags.times]
	print("Number of training data after parsing based on times: {}".format(len(train)))
	
print("Number of final training data after parsing: {}".format(len(train)))
	
with open(flags.output_dir, 'wb') as f:
	pickle.dump(train, f)