"""
use
python -m unittest discover
or
python -m unittest test.test_table
from the main folder
"""

from geecq.table import Table
from geecq.fastqc import Fastqc
from geecq.sam import Sam
from geecq.meta import Meta
import unittest
import os

ROOTDIR = os.path.dirname(__file__)

class TestTable(unittest.TestCase):
    def test_make_csv(self):
        fastqc1 = Fastqc()
        fastqc1.load_from_file(ROOTDIR + '/fastqc_v10.txt')
        fastqc2 = Fastqc()
        fastqc2.load_from_file(ROOTDIR + '/fastqc_v11.txt')
        sam = Sam()
        sam.load_from_file(ROOTDIR + '/samstat_v2.html')
        meta = Meta('dummy')
        short_table, long_table = Table([[fastqc1, fastqc2, sam, meta]]).make_tables()
        self.assertTrue(all(len(row) == len(short_table[0]) for row in short_table))
        self.assertTrue(all(len(row) == len(long_table[0]) for row in long_table))
