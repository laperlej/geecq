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
        self.MAPQ = []
        self.name = ''
    
    def convertName(self, rawName):
        """Converts the name to be coherent with the other files
        
        Args:
            rawName: unconverted name(string)
        """
        self.name = self.name.replace('_1', '')
        self.name = self.name.replace('_2', '.1')
    
    def loadFromString(self, fileTxt):
        """Extracts data from a .sam.html file
        
        Args: 
            fileTxt: raw text from a .sam.html file(string)
        """
        for line in fileTxt.split('\n'):
            if '<title>Library' in line:
                self.convertName(line.split()[1].replace('.sam</title>', ''))
                self.name = self.name.replace('.sam</title>', '')
                self.convertName(self.name)
                
            if 'MAPQ' in line or 'Unmapped' in  line:
                line = re.search('\(.*\((.*?)% , (.*?)\)', line)
                if line:
                    self.MAPQ.append([line.group(1), line.group(2)])
        
    def loadFromFile(self, fileName):
        """Extracts data from a .sam.html file
        
        Args: 
            fileName: path to a .sam.html file
        """
        try:
            self.loadFromString(open(fileName, 'r').read())
        except IOError:
            print('Could not open ' + fileName)