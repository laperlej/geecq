"""
use
python -m unittest discover
or
python -m unittest test.test_fastqc
from the main folder
"""

from geecq.fastqc import Fastqc
import unittest
import os

ROOTDIR = os.path.dirname(__file__)

class TestFastqc(unittest.TestCase):

    def is_loaded(self, fastqc):
        self.assertTrue(fastqc.name)
        self.assertTrue(fastqc.version)
        self.assertTrue(fastqc.nb_sequences)
        self.assertTrue(fastqc.gc_content)
        self.assertTrue(fastqc.pos_quality)
        self.assertTrue(fastqc.seq_length)
        self.assertTrue(fastqc.qual)
        self.assertTrue(fastqc.dup)

    def test_load_from_string(self):
        fastqc = Fastqc()
        fastqc.load_from_string(open(ROOTDIR + '/fastqc_v10.txt', 'r').read())
        self.is_loaded(fastqc)

    def test_load_from_file(self):
        fastqc = Fastqc()
        fastqc.load_from_file(ROOTDIR + '/fastqc_v10.txt')
        self.is_loaded(fastqc)

    def test_load_from_zip(self):
        fastqc = Fastqc()
        fastqc.load_from_zip(ROOTDIR + '/fastqc_v10.zip')
        self.is_loaded(fastqc)
