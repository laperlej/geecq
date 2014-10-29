"""
Created on January 25, 2014

@author: Jonathan Laperle(jonathan.laperle@usherbrooke.ca)
"""

import re

class MetaOld(object):
    """Used to extract and contain data from a GEO meta data file

    Attributes:
        meta: dictionnary where the key is the row name and the value
              is the value of the row
            example: meta['Instrument'] = 'Illumina Genome Analyzer II'
        name:
            name of the fastq file
    """
    def __init__(self, name):
        self.meta = {}
        self.name = name

    def load_from_string(self, file_txt):
        """Extracts data from a meta data file's content

        Args:
            fileTxt: raw text from a meta data file(string)
        """
        for line in file_txt.split('\n'):
            if '<td>GEO</td>' in line:
                string = filter(None, re.sub(r'<[^>]+>', '|', line).split('|'))
                self.meta[string[2]] = string[3] #Instrument
                self.meta[string[6]] = string[7] #Alias
                self.meta[string[12]] = string[13] #Organism
                self.meta[string[16]] = string[17] #Library name
                self.meta[string[18]] = string[19] #Library strategy
                self.meta[string[24]] = string[25] #Library layout
                self.meta[string[32]] = string[33] #Submissions

    def load_from_file(self, file_name):
        """Extracts data from a meta data file

        Args:
            fileName: path to a meta data file
        """
        try:
            self.load_from_string(open(file_name, 'r').read())
        except IOError:
            print 'Could not open ' + file_name

class Meta(object):
    """Used to extract and contain data from a meta data file

    Attributes:
        meta: dictionnary where the key is the row name and the value
              is the value of the row
            example: meta['Instrument'] = 'Illumina Genome Analyzer II'
        name:
            name of the fastq file
    """
    def __init__(self, name):
        self.meta = {}
        self.name = name

    def load_from_string(self, file_txt):
        """Extracts data from a meta data file's content

        Args:
            fileTxt: raw text from a meta data file(string)
        """
        for line in file_txt.split('\n'):
            prog = re.compile("!Sample_(.*?) = (.*?)\r\n?|\n")
            if line.startswith("!Sample_"):
                info = prog.match(line)
                if info.groups()[0] in self.meta:
                    self.meta[info.groups()[0]] += info.groups()[1]
                else:
                    self.meta[info.groups()[0]] = info.groups()[1]
            elif line.startswith("^SAMPLE"):
                matched_name = re.match(r"\^SAMPLE = (.*?)\r\n?|\n", line)
                self.name = matched_name.groups()[0]

    def load_from_file(self, file_name):
        """Extracts data from a meta data file

        Args:
            fileName: path to a meta data file
        """
        try:
            self.load_from_string(open(file_name, 'r').read())
        except IOError:
            print 'Could not open ' + file_name
