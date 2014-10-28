"""
Created on January 25, 2014

@author: Jonathan Laperle(jonathan.laperle@usherbrooke.ca)
"""

import re

class Sam(object):
    """Used to extract and contain data from a .sam.html file

    Attributes:
        MAPQ: list of lists with each element of the form [%, nb sequences]
            follows the following order >=30, >=20, >=10, >=3, < 3, Unmapped
        name: name of the fastq file
    """
    def __init__(self):
        self.mapq = []
        self.name = ''

    def convert_name(self, raw_name):
        """Converts the name to be coherent with the other files

        Args:
            rawName: unconverted name(string)
        """
        self.name = raw_name.replace('_1', '')
        self.name = raw_name.replace('_2', '.1')

    def load_from_string(self, file_txt):
        """Extracts data from a .sam.html file

        Args:
            fileTxt: raw text from a .sam.html file(string)
        """
        for line in file_txt.split('\n'):
            if '<title>Library' in line:
                self.convert_name(line.split()[1].replace('.sam</title>', ''))
                self.name = self.name.replace('.sam</title>', '')
                self.convert_name(self.name)

            if 'MAPQ' in line or 'Unmapped' in  line:
                line = re.search('\(.*\((.*?)% , (.*?)\)', line)
                if line:
                    self.mapq.append([line.group(1), line.group(2)])

    def load_from_file(self, file_name):
        """Extracts data from a .sam.html file

        Args:
            fileName: path to a .sam.html file
        """
        try:
            self.load_from_string(open(file_name, 'r').read())
        except IOError:
            print 'Could not open ' + file_name
