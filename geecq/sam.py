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

    def convert_name_(self, raw_name):
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
        version = 2
        lines = file_txt.split('\n')
        for line in file_txt.split('\n'):
            if '<title>' in line:
                raw_name = re.search('<title>(.*)</title>', line).group(1)
                if 'Library ' in raw_name:#if old samstat file
                    raw_name = raw_name.replace('Library ', '')
                    version = 1
                self.convert_name_(raw_name)
                self.name = self.name.replace('.sam', '')
                self.convert_name_(self.name)
                break

        if not self.name:
            print "ERROR: Invalid sam file, missing a <title> tag"
            exit(1)

        if version == 1:
            for line in file_txt.split('\n'):
                if 'MAPQ' in line or 'Unmapped' in  line:
                    line = re.search('\(.*\((.*?)% , (.*?)\)', line)
                    if line:
                        self.mapq.append([line.group(1), line.group(2)])
        else:
            for line in file_txt.split('tr>'):
                if '<td style="background-color: rgba(' in line and "Total" not in line:
                    line = re.search('<td>(.*?)\.0</td>\n<td>(.*?)</td>', line)
                    self.mapq.append([line.group(2), line.group(1)])


    def load_from_file(self, file_name):
        """Extracts data from a .sam.html file

        Args:
            fileName: path to a .sam.html file
        """
        try:
            self.load_from_string(open(file_name, 'r').read())
        except IOError:
            print 'Could not open ' + file_name
