"""
Test that utils file works.
"""
import unittest

from ..utils import load_pickle
from ..utils import cell_line_list

data = load_pickle("data/sample.pkl")


class TestUtils(unittest.TestCase):
  """
  Tests that utils functions work as expected.
  """

  def test_cell_line_list(self):
    """Check the number of cell lines in output"""
    cells = ['HL60']
    parse_data = cell_line_list("dummy", cells, data=data)
    output_cells = [line[0][0] for line in parse_data]
    self.assertEqual(len(cells), len(list(set(output_cells))))
    
    

if __name__ == '__main__':
    unittest.main()
    

