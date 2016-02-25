"""
Created on January 25, 2014

@author: Jonathan Laperle(jonathan.laperle@usherbrooke.ca)
"""

import zipfile

class Fastqc(object):
    """Used to extract and contain data from a fastqc_data.txt file

    Attributes:
        fileName: name of the fastq file
        nbSequences: total number of sequences
        gcContent: content in gc (%)
        posQuality: list of numbers with index = the position
            and the value = average quality for that position
        seqLength: List of lists of format [Length, nb sequences]
        qual: list of numbers with index = quality score and
              value = nb sequences
        dup: list of numbers with index = duplication level and
             value = nb sequences, [0] is Total
    """
    def __init__(self):
        self.name = ''
        self.nb_sequences = 0
        self.gc_content = 0
        self.pos_quality = []
        self.seq_length = []
        self.qual = []
        self.dup = []

    def load_basic(self, module):
        """Extracts data from the basic statistics section of a
        fastqc_data.txt file

        Args:
            module: the raw text of the section(delimited by '>>'
                which should not be included in the string)
        """
        for line in module.split('\n'):
            if line.startswith('Filename'):
                self.name = line.split()[-1]
            if line.startswith('Total Sequences'):
                self.nb_sequences = int(line.split()[-1])
            if line.startswith('%GC'):
                self.gc_content = int(line.split()[-1])

    def load_pos_quality(self, module):
        """Extracts the data from the Per base sequence quality section
        of a fastqc_data.txt file

        Args:
            module: the raw text of the section(delimited by '>>'
                which should not be included in the string)
        """
        for line in module.split('\n'):
            if line and line[0].isdigit():
                line = line.split()
                self.pos_quality.append(float(line[2]))

    def load_qual(self, module):
        """Extracts data from the Per sequence quality scores section
        of a fastqc_data.txt file

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

    def load_length(self, module):
        """Extracts data from the Sequence Length Distribution
        section of a fastqc_data.txt file

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
                        self.seq_length.append([int(x), float(line[1])/2])
                else:
                    x = int(line[0])
                    y = float(line[1])
                    self.seq_length.append([x, y])

    def load_dup(self, module):
        """Extracts data from the Sequence Duplication Levels of a
        fastqc_data.txt file

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

    def load_from_string(self, file_txt):
        """Extracts data from a string

        Args:
            fileTxt: raw text from a fastqc_data.txt file(string)
        """
        for module in file_txt.split('>>'):
            if 'Basic Statistics' in module:
                self.load_basic(module)
            if 'Per base sequence quality' in module:
                self.load_pos_quality(module)
            if 'Per sequence quality scores' in module:
                self.load_qual(module)
            if 'Sequence Length Distribution' in module:
                self.load_length(module)
            if 'Sequence Duplication Levels' in module:
                self.load_dup(module)

    def load_from_file(self, file_name):
        """Extracts data from a fastqc_data.txt file

        Args:
            fileName: path to a fastqc_data.txt file
        """
        try:
            self.load_from_string(open(file_name, 'r').read())
        except IOError:
            print 'Could not open ' + file_name

    def load_from_zip(self, zip_file):
        """Extracts data from a zip archive

        Args:
            fileName: path to a zip archive containing a fastqc_data.txt file
        """
        try:
            z = zipfile.ZipFile(zip_file)
            for filename in z.namelist():
                if filename.endswith('fastqc_data.txt'):
                    f = z.read(filename)
                    self.load_from_string(f)
            z.close()
        except (zipfile.BadZipfile, IOError):
            print 'Could not open ' + zip_file
