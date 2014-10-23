"""
Created on January 25, 2014

@author: Jonathan Laperle(jonathan.laperle@usherbrooke.ca)
"""

import zipfile

class Fastqc():
    """Used to extract and contain data from a fastqc_data.txt file
    
    Attributes:
        fileName: name of the fastq file
        nbSequences: total number of sequences
        gcContent: content in gc (%)
        posQuality: list of numbers with index = the position
            and the value = average quality for that position
        seqLength: List of lists of format [Length, nb sequences]
        qual: list of numbers with index = quality score and value = nb sequences
        dup: list of numbers with index = duplication level and value = nb sequences
            [0] is Total
    """
    def __init__(self):
        self.fileName = ''
        self.nbSequences = 0
        self.gcContent = 0
        self.posQuality = []
        self.seqLength = []
        self.qual = []
        self.dup = []

    def loadBasic(self, module):
        """Extracts data from the basic statistics section of a fastqc_data.txt file
        
        Args:
            module: the raw text of the section(delimited by '>>'
                which should not be included in the string)
        """
        for line in module.split('\n'):
            if line.startswith('Filename'):
                self.fileName = line.split()[-1]
            if line.startswith('Total Sequences'):
                self.nbSequences = int(line.split()[-1])
            if line.startswith('%GC'):
                self.gcContent = int(line.split()[-1])

    def loadPosQuality(self, module):
        """Extracts the data from the Per base sequence quality section of a fastqc_data.txt file
        
        Args:
            module: the raw text of the section(delimited by '>>'
                which should not be included in the string)
        """
        for line in module.split('\n'):
            if line and line[0].isdigit():
                line = line.split()
                self.posQuality.append(float(line[2]))

    def loadQual(self, module):
        """Extracts data from the Per sequence quality scores section of a fastqc_data.txt file
        
        Args:
            module: the raw text of the section(delimited by '>>'
                which should not be included in the string)
        """
        for line in module.split('\n'):
            if line and line[0].isdigit():
                line = line.split()
                while int(line[0]) > len(self.qual):
                    self.qual.append(0.0)
                self.qual.append(float(line[1]))

    def loadLength(self, module):
        """Extracts data from the Sequence Length Distribution section of a fastqc_data.txt file
        
        Args:
            module: the raw text of the section(delimited by '>>'
                which should not be included in the string)
        """
        for line in module.split('\n'):
            if line and line[0].isdigit():
                line = line.split()
                if '-' in line[0]:
                    line2 = line[0].split('-')
                    for x in line2:
                        self.seqLength.append([int(x), float(line[1])/2])
                else:
                    x = int(line[0])
                    y = float(line[1])
                    self.seqLength.append([x, y])

    def loadDup(self, module):
        """Extracts data from the Sequence Duplication Levels of a fastqc_data.txt file
        
        Args:
            module: the raw text of the section(delimited by '>>'
                which should not be included in the string)
        """
        for line in module.split('\n'):
            if line.startswith('#Total Duplicate Percentage'):
                y = float(line.split()[-1])
                self.dup.append(y)
            if line and line[0].isdigit():
                line = line.replace('+', '')
                line = line.split()
                y = float(line[1])
                self.dup.append(y)

    def loadFromString(self, fileTxt):
        """Extracts data from a string
        
        Args:
            fileTxt: raw text from a fastqc_data.txt file(string)
        """
        for module in fileTxt.split('>>'):
            if 'Basic Statistics' in module:
                self.loadBasic(module)
            if 'Per base sequence quality' in module:
                self.loadPosQuality(module)
            if 'Per sequence quality scores' in module:
                self.loadQual(module)
            if 'Sequence Length Distribution' in module:
                self.loadLength(module)
            if 'Sequence Duplication Levels' in module:
                self.loadDup(module)
        
    def loadFromFile(self, fileName):
        """Extracts data from a fastqc_data.txt file
        
        Args: 
            fileName: path to a fastqc_data.txt file
        """
        try:
            self.loadFromString(open(fileName, 'r').read())
        except IOError:
            print('Could not open ' + fileName)
    
    def loadFromZip(self, zipFile):
        """Extracts data from a zip archive
        
        Args: 
            fileName: path to a zip archive containing a fastqc_data.txt file
        """
        try:
            z=zipfile.ZipFile(zipFile)
            for filename in z.namelist():
                if filename.endswith('fastqc_data.txt'):
                    f = z.read(filename)
                    self.loadFromString(f)
            z.close()
        except zipfile.BadZipfile:
            print('Could not open ' + zipFile)