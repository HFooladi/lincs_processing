import pandas as pd
from collections import Counter

__author__ = "Hosein Fooladi"
__email__ = "fooladi.hosein@gmail.com"


def sig_info_augment(sig_info_dir: str):
  """ Unification between GSE70138 and GSE92742

  This function has been written for working with GSE70138.
  Unfortunately, format of sig_info (number of columns) differs between
  GSE92742 and GSE70138. So, I have written this function for unification between 
  these two dataset. Particularly, I am going to add 4 columns (pert_dose, pert_dose_unit,
  pert_time, pert_time_unit) to sig_info of GSE70138.
	
  Parameters
  ----------
  sig_info_dir: str  
    The directory of sig_info file. E.g., './Data/sig_info.txt'
			
  Returns
  -------
  sig_info_v1: pd.DataFrame 
    A dataframe. It is like the input file, except it has four more columns.
  """

  assert isinstance(sig_info_dir,
                    str), "The dataset_dir must be a string object"

  sig_info = pd.read_csv(sig_info_dir, sep='\t')

  print("Data Statistics\n")
  print("Number of available gene expression signature: {}".format(
      sig_info.shape[0]))
  print("Number of available columns: {}".format(sig_info.shape[1]))

  dose = [
      x.split() if x != '-666' else ['-666', '-666']
      for x in sig_info.pert_idose
  ]
  time = [
      x.split() if x != '-666' else ['-666', '-666']
      for x in sig_info.pert_itime
  ]

  pert_dose = list(map(lambda x: float(x[0]), dose))
  pert_dose_unit = list(map(lambda x: x[1], dose))
  pert_time = list(map(lambda x: float(x[0]), time))
  pert_time_unit = list(map(lambda x: x[1], time))

  adding = pd.DataFrame({
      'pert_dose': pert_dose,
      'pert_dose_unit': pert_dose_unit,
      'pert_time': pert_time,
      'pert_time_unit': pert_time_unit
  })

  sig_info_v1 = pd.concat([sig_info, adding], axis=1)
  print("Number of available columns after augmentation: {}".format(
      sig_info_v1.shape[1]))

  return sig_info_v1
