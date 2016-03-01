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
        self.version = 0
        self.nb_sequences = 0
        self.gc_content = 0
        self.pos_quality = []
        self.seq_length = []
        self.qual = []
        self.dup = []

    def load_version_(self, module):
        """find the version number of the fastqc file, extracts the 2nd number
        note: 1.10 is very different from 1.11 and requires new code
        """
        for line in module.split('\n'):
            if line.startswith('##FastQC'):
                self.version = int(line.split()[1].split('.')[1])

    def load_basic_(self, module):
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

    def load_pos_quality_(self, module):
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

    def load_qual_(self, module):
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

    def load_length_(self, module):
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

    def load_dup_(self, module):
        """Extracts data from the Sequence Duplication Levels of a
        fastqc_data.txt file

        Args:
            module: the raw text of the section(delimited by '>>'
                which should not be included in the string)
        """
        for line in module.split('\n'):
            if line.startswith('#Total Duplicate Percentage') or line.startswith('#Total Deduplicated Percentage'):
                y = float(line.split()[-1])
                self.dup.append(y)
            line = line.replace('>','')
            line = line.replace('+', '')
            line = line.replace('k', '000')
            if line and line[0].isdigit():
                line = line.split()
                if int(line[0]) > 10:
                    self.dup[-1] = self.dup[-1] + float(line[1])
                else:
                    y = float(line[1])
                    self.dup.append(y)

    def load_from_string(self, file_txt):
        """Extracts data from a string

        Args:
            fileTxt: raw text from a fastqc_data.txt file(string)
        """
        for module in file_txt.split('>>'):
            if module.startswith('##FastQC'):
                self.load_version_(module)
            if module.startswith('Basic Statistics'):
                self.load_basic_(module)
            if module.startswith('Per base sequence quality'):
                self.load_pos_quality_(module)
            if module.startswith('Per sequence quality scores'):
                self.load_qual_(module)
            if module.startswith('Sequence Length Distribution'):
                self.load_length_(module)
            if module.startswith('Sequence Duplication Levels'):
                self.load_dup_(module)

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
            zip_file = zipfile.ZipFile(zip_file)
            for filename in zip_file.namelist():
                if filename.endswith('fastqc_data.txt'):
                    txt_file = zip_file.read(filename)
                    self.load_from_string(txt_file)
            zip_file.close()
        except (zipfile.BadZipfile, IOError):
            print 'Could not open ' + zip_file
