"""
Created on January 25, 2014

@author: Jonathan Laperle(jonathan.laperle@usherbrooke.ca)
"""

import re

class Meta(object):
    """Used to extract and contain data from a meta data file
    
    Attributes:
        meta: dictionnary where the key is the row name and the value is the value of the row
            example: meta['Instrument'] = 'Illumina Genome Analyzer II'
        name:
            name of the fastq file
    """
    def __init__(self, name):
        self.meta = {}
        self.name = name

    def loadFromString(self, fileTxt):
        """Extracts data from a meta data file's content
        
        Args: 
            fileTxt: raw text from a meta data file(string)
        """
        for line in fileTxt.split('\n'):
            if '<td>GEO</td>' in line:
                string = filter(None, re.sub(r'<[^>]+>', '|', line).split('|'))
                self.meta[string[2]]=string[3] #Instrument
                self.meta[string[6]]=string[7] #Alias
                self.meta[string[12]]=string[13] #Organism
                self.meta[string[16]]=string[17] #Library name
                self.meta[string[18]]=string[19] #Library strategy
                self.meta[string[24]]=string[25] #Library layout
                self.meta[string[32]]=string[33] #Submissions
    
    def loadFromFile(self, fileName):
        """Extracts data from a meta data file
        
        Args: 
            fileName: path to a meta data file
        """
        try:
            self.loadFromString(open(fileName, 'r').read())
        except IOError:
            print('Could not open ' + fileName)