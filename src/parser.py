from __future__ import unicode_literals, print_function, division

import numpy as np
import pandas as pd
from collections import Counter
from cmapPy.pandasGEXpress.parse import parse
import cmapPy.pandasGEXpress.write_gctx as wg
from typing import List, Tuple, Optional

__author__ = "Hosein Fooladi"
__email__ = "fooladi.hosein@gmail.com"


def parsing_level3_cp(dataset_dir: str,
                      inst_info_dir: str,
                      gene_info_dir: str,
                      pert_type: str = "trt_cp",
                      landmarks: bool = True) -> List[List]:
  """Parsing the data to keep desired sig_ids
  
  This function takes the directory of dataset, perturbation type, and
  whether we want to only keep landmark genes or not. It returns a list
  based on the inputs.


  Parameters
  ----------
  dataset_dir: str
    It must be string file that shows the directory of the dataset.
    dataset should be a gctx file. e.g., valid argument is something like this:
    './Data/Level3_INF_mlr12k_n1319138x12328.gctx'
  param inst_info_dir: str 
    directory of inst_info. It contains the information about the
    experiment, perturbation type and cell line. For example:
    './Data/inst_info.txt'
  gene_info_dir: str 
    directory of gene_info. It contains the information about the genes.
    For example: './Data/gene_info.txt'
  pert_type: str (default= "trt_cp")
    String object that determine which perturbation type you want to parse.
    Default='trt_cp'
  landmarks: bool
    boolean which determines whether you want to just keep landmark genes
    after parsing or you want to keep all the genes. Default=True

  Returns
  ------
  parse_list: List
    Output list (Train, Validation, Test) Format:
    line[0]:(cell_line,
    drug,
    drug_type,
    does,
    does_type,
    time,
    time_type)
    line[1]: 978 or 12328-dimensional Vector(Gene_expression_profile)

  """

  assert isinstance(pert_type, str), "pert_type must be a string object"
  assert isinstance(landmarks, bool), "landmarks must be a boolean object"

  gene_info = pd.read_csv(gene_info_dir, sep="\t", dtype=str)
  print("Number of measured genes in the dataset: {}".format(
      gene_info.shape[0]))

  landmark_gene_row_ids = gene_info["gene_id"][gene_info["is_lm"] == "1"]
  print("Number of landmark genes in the dataset: {}".format(
      landmark_gene_row_ids.shape[0]))

  inst_info = pd.read_csv(inst_info_dir, sep="\t")
  print("Number of availbale gene expression profiles: {}".format(
      inst_info.shape[0]))
  print("Unique perturbation types: {}".format(inst_info.pert_type.unique()))

  assert pert_type in inst_info.pert_type.unique(), "pert_type is not valid!!"

  query_trt = inst_info[inst_info["pert_type"] == pert_type]
  print("Number of availbale gene expression profiles of {}: {}".format(
      pert_type, query_trt.shape[0]))

  print((query_trt == '-666').sum())
  print("Number of different pert_dose_units: {}".format(
      Counter(query_trt.pert_dose_unit)))

  ## It needs to be better. I am supposed to modify this part!
  if pert_type == 'trt_cp':
    query_trt = query_trt[query_trt.pert_dose_unit != '-666']
    query_trt = query_trt[query_trt.pert_dose_unit == 'um']
  else:
    pass

  query_ids = query_trt.inst_id
  print("Number of samples at the end: {}".format(query_ids.shape[0]))

  print("=================================================================")
  print("Please wait while we are parsing the data ...")

  if landmarks:
    query_gctoo = parse(dataset_dir, rid=landmark_gene_row_ids, cid=query_ids)
  else:
    query_gctoo = parse(dataset_dir, cid=query_ids)

  print("Parse Completed")
  print("Size of the data after parsing: {}".format(query_gctoo.data_df.shape))

  query_trt = query_trt.set_index(query_trt.inst_id)
  query_trt = query_trt.reindex(query_gctoo.data_df.columns)

  parse_list = []
  for i in range(query_trt.shape[0]):
    parse_list.append([
        (query_trt.cell_id[i], query_trt.pert_id[i], query_trt.pert_type[i],
         query_trt.pert_dose[i], query_trt.pert_dose_unit[i],
         query_trt.pert_time[i], query_trt.pert_time_unit[i]),
        np.array(query_gctoo.data_df.iloc[:, i])
    ])

  return parse_list


def parsing_level5_cp(dataset_dir: str,
                      sig_info_dir: str,
                      gene_info_dir: str,
                      pert_type: str = 'trt_cp',
                      landmarks: bool = True,
                      cell_line: Optional[str] = None) -> List[List]:
  """Parsing the data to keep desired sig_ids
  
  This function takes the directory of dataset, perturbation type, and
  whether we want to only keep landmark genes or not. It returns a list
  based on the inputs.


  Parameters
  ----------
  dataset_dir: str
    It must be string file that shows the directory of the dataset.
    dataset should be a gctx file. e.g., valid argument is something like this:
    './Data/Level3_INF_mlr12k_n1319138x12328.gctx'
  sig_info_dir: str 
    directory of sig_info. It contains the information about the
    experiment, perturbation type and cell line. For example:
    './Data/sig_info.txt'
  gene_info_dir: str
    directory of gene_info. It contains the information about the genes.
    For example: './Data/gene_info.txt'
  pert_type: str (default="trt_cp") 
    String object that determine which perturbation type you want to parse.
    Default='trt_cp'
  landmarks: bool (default=True)
    boolean which determines whether you want to just keep landmark genes
    after parsing or you want to keep all the genes. Default=True
  cell_line: str (default=None)
    Whether you want to select a particular cell_line and parse data just
    for that cell line or not. Default=None Which means parse information of all the cell lines.

  Returns
  -------
  parse_list: List
    Output list (Train, Validation, Test) Format:
    line[0]:(cell_line,
    drug,
    drug_type,
    does,
    does_type,
    time,
    time_type)
    line[1]: 978 or 12328-dimensional Vector(Gene_expression_profile)

  """

  assert isinstance(pert_type, str), "pert_type must be a string object"
  assert isinstance(landmarks, bool), "landmarks must be a boolean object"

  gene_info = pd.read_csv(gene_info_dir, sep="\t", dtype=str)
  print("Number of measured genes in the dataset: {}".format(
      gene_info.shape[0]))

  landmark_gene_row_ids = gene_info["gene_id"][gene_info["is_lm"] == "1"]
  print("Number of landmark genes in the dataset: {}".format(
      landmark_gene_row_ids.shape[0]))

  sig_info = pd.read_csv(sig_info_dir, sep="\t")
  print("Number of availbale gene expression profiles: {}".format(
      sig_info.shape[0]))
  print("Unique perturbation types: {}".format(sig_info.pert_type.unique()))

  assert pert_type in sig_info.pert_type.unique(), "pert_type is not valid!!"

  query_trt = sig_info[sig_info["pert_type"] == pert_type]
  print("Number of availbale gene expression profiles of {}: {}".format(
      pert_type, query_trt.shape[0]))

  print((query_trt == '-666').sum())
  print("Number of different pert_dose_units: {}".format(
      Counter(query_trt.pert_dose_unit)))

  ## It needs to be better. I am supposed to modify this part!
  if pert_type == 'trt_cp':
    query_trt = query_trt[query_trt.pert_dose_unit != '-666']
    #query_trt = query_trt[query_trt.pert_dose_unit == 'um']
  else:
    pass

  if cell_line is not None:
    query_trt = query_trt[sig_info["cell_id"] == cell_line]
  else:
    pass

  query_ids = query_trt.sig_id
  print("Number of samples at the end: {}".format(query_ids.shape[0]))

  print("=================================================================")
  print("Please wait while we are parsing the data ...")

  if landmarks:
    query_gctoo = parse(dataset_dir, rid=landmark_gene_row_ids, cid=query_ids)
  else:
    query_gctoo = parse(dataset_dir, cid=query_ids)

  print("Parse Completed")
  print("Size of the data after parsing: {}".format(query_gctoo.data_df.shape))

  query_trt = query_trt.set_index(query_trt.sig_id)
  query_trt = query_trt.reindex(query_gctoo.data_df.columns)

  parse_list = []
  for i in range(query_trt.shape[0]):
    parse_list.append([
        (query_trt.cell_id[i], query_trt.pert_id[i], query_trt.pert_type[i],
         float(query_trt.pert_dose[i]), query_trt.pert_dose_unit[i],
         query_trt.pert_time[i], query_trt.pert_time_unit[i]),
        np.array(query_gctoo.data_df.iloc[:, i])
    ])

  return parse_list
