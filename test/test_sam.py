"""
use
python -m unittest discover
or
python -m unittest test.test_sam
from the main folder
"""

from geecq.sam import Sam
import unittest
import os

ROOTDIR = os.path.dirname(__file__)

class TestSam(unittest.TestCase):

    def is_loaded(self, fastqc):
        self.assertTrue(fastqc.name)
        self.assertTrue(fastqc.mapq)

    def test_load_from_string(self):
        sam = Sam()
        sam.load_from_string(open(ROOTDIR + '/samstat_v1.html', 'r').read())
        self.is_loaded(sam)

    def test_load_from_file(self):
        sam = Sam()
        sam.load_from_file(ROOTDIR + '/samstat_v1.html')
        self.is_loaded(sam)

    def test_load_v2(self):
        sam = Sam()
        sam.load_from_file(ROOTDIR + '/samstat_v2.html')
        self.is_loaded(sam)
