"""
Created on January 25, 2014

@author: Jonathan Laperle(jonathan.laperle@usherbrooke.ca)
"""
import others

class Table(object):
    """Used to generate .TAB tables from fastqc/sam/meta files

    Attributes:
        matchedList:
        path: output path
    """
    def __init__(self, matchedLists):
        self.matched_lists = matchedLists
        self.short_csv = []
        self.long_csv = []

    def max_qual_(self):
        """Maxmimal quality value for all fastqc objects
        """
        return max([max([len(qc.qual) for qc in \
            [row[0] for row in self.matched_lists]]), \
            max([len(qc.qual) for qc in \
            [row[1] for row in self.matched_lists]])])

    def min_length_(self):
        """Minimal length for all fastqc objects
        """
        temp = []
        for i in range(2):
            for fastqc in [row[i] for row in self.matched_lists]:
                if fastqc.name:
                    temp.append(fastqc.seq_length[0][0])
        return min(temp)

    def max_length_(self):
        """Maximal length for all fastqc objects in fastqcList
        """
        temp = []
        for i in range(2):
            for fastqc in [row[i] for row in self.matched_lists]:
                if fastqc.name:
                    temp.append(fastqc.seq_length[-1][0])
        return max(temp)

    def header_short_(self):
        """Writes the header for the short table file
        """
        short_header = ['File Name',
                        'Library Name',
                        'nb Sequences',
                        'nb Sequences (trimmed)',
                        'GC content(%)',
                        'GC content(%) (trimmed)',
                        'Mean quality score of first 10 bases',
                        'Mean quality score of last 10 bases',
                        'Mean quality score of first 10 bases (trimmed)',
                        'Mean quality score of last 10 bases (trimmed)',
                        'mode for average quality',
                        '% of sequences with the mode',
                        'mode for average quality (trimmed)',
                        '% of sequences with the mode (trimmed)',
                        'mode for sequence length',
                        '% of sequences with the mode',
                        'mode for sequence length (trimmed)',
                        '% of sequences with the mode (trimmed)',
                        '% Duplication',
                        '% Duplication (trimmed)',
                        'MAPQ >= 20 (%)',
                        'MAPQ >= 20',
                        '20 > MAPQ >= 3 (%)',
                        '20 > MAPQ >= 3',
                        'MAPQ < 3 or unmapped (%)',
                        'MAPQ < 3 or unmapped']

        self.short_csv.append(short_header)

    def header_long_(self):
        """writes the header for the long table file
        """
        self.long_csv.append(self.first_header_level_())
        self.long_csv.append(self.second_header_level_())
        self.long_csv.append(self.third_header_level_())

    def first_header_level_(self):
        """Generates the first level for the long table's header

        Returns:
            the header as a string
        """
        first_header = []
        first_header.append('General information')
        first_header += empty_slot(9)
        for _ in range(2):
            first_header.append('Per base sequence quality')
            first_header += empty_slot(19)
        for _ in range(2):
            first_header.append('Number of sequences per quality score')
            first_header += empty_slot(self.max_qual_()-1)
        for _ in range(2):
            first_header.append('Sequence Length Distribution')
            first_header += empty_slot(self.max_length_() - self.min_length_())
        for _ in range(2):
            first_header.append('Sequences duplication levels')
            first_header += empty_slot(10)
        first_header.append('MAPQ')
        first_header += empty_slot(11)
        return first_header

    def second_header_level_(self):
        """Generates the second level for the long table's header

        Returns:
            the header as a string
        """
        #General information
        second_header = ['File Name',
                         'nb Sequences',
                         'nb Sequences (trimmed)',
                         'GC content(%)',
                         'GC content(%) (trimmed)',
                         'Instrument',
                         'Organism',
                         'Library Name',
                         'Library Strategy',
                         'Characteristics']
        #in order:
        #Per base sequence quality
        #Per sequence quality scores
        #Sequence Length Distribution
        #Sequences duplication levels
        for nb_empty in [19,
                         self.max_qual_() - 1,
                         self.max_length_() - self.min_length_(),
                         10]:
            second_header += ['Untrimmed'] + empty_slot(nb_empty)
            second_header += ['Trimmed'] + empty_slot(nb_empty)

        #MAPQ
        second_header += empty_slot(12)

        return second_header

    def third_header_level_(self):
        """Generates the third level for the long table's header

        Returns:
            the header as a string
        """
        third_header = []

        #General information
        third_header += empty_slot(10)

        #Per base sequence quality
        for _ in range(2):
            third_header += [str(i+1) for i in range(10)]
            third_header += [str(i-10) for i in range(10)]

        #Per sequence quality scores
        for _ in range(2):
            third_header += [str(i) for i in range(self.max_qual_())]

        #Sequence Length Distribution
        for _ in range(2):
            third_header += [str(i) for i in range(self.min_length_(),
                                                   self.max_length_()+1)]

        #Sequences duplication levels
        for i in range(2):
            third_header.append('Total')
            third_header += [str(i+1) for i in range(9)]
            third_header.append('10+')

        #MAPQ
        third_header += ["MAPQ >= 30 (%)",
                         "MAPQ >= 30",
                         "MAPQ >= 20 (%)",
                         "MAPQ >= 20",
                         "MAPQ >= 10 (%)",
                         "MAPQ >= 10",
                         "MAPQ >= 3 (%)",
                         "MAPQ >= 3",
                         "MAPQ < 3 (%)",
                         "MAPQ < 3",
                         "Unmapped (%)",
                         "Unmapped"]

        return third_header

    def table_short_(self, matched_list):
        """Appends the content of the short table for matching
        before/after/sam/meta files

        Args:
            qcbefore:Fastqc object before trimming
            qcafter: Fastqc object after trimming
            sam: Sam object
            meta: Meta object
        """
        fastqcs = [matched_list[0], matched_list[1]]
        sam = matched_list[2]
        meta = matched_list[3]
        empty = ['-']
        #File name
        output = [search_name(matched_list)]

        #Library name
        if meta.meta:
            output.append(meta.meta['title'])
        else:
            output += empty

        #nb sequences
        for fastqc in fastqcs:
            if fastqc.nb_sequences:
                output.append(str(fastqc.nb_sequences))
            else:
                output += empty

        #GC content
        for fastqc in fastqcs:
            if fastqc.gc_content:
                output.append(str(fastqc.gc_content))
            else:
                output += empty

        #Average base sequence quality for first and last 10 bases
        for fastqc in fastqcs:
            if fastqc.pos_quality:
                output.append(str(others.mean(fastqc.pos_quality[:10])))
                output.append(str(others.mean(fastqc.pos_quality[-10:])))
            else:
                output += empty * 2

        #Mode for sequence quality scores and how many sequences ahve the mode
        #in absolute and %
        for fastqc in fastqcs:
            if fastqc.qual:
                i = others.argmax(fastqc.qual)
                output += [str(i)]
                output += [str(fastqc.qual[i]/sum(fastqc.qual))]
            else:
                output += empty * 2

        #Mode for sequence lengthand how many sequences ahve the mode
        #in absolute and %
        for fastqc in fastqcs:
            if fastqc.seq_length:
                i = others.argmax([row[1] for row in fastqc.seq_length])
                output += [str(fastqc.seq_length[i][0])]
                output += [str(fastqc.seq_length[i][1]/\
                           sum([row[1] for row in fastqc.seq_length]))]
            else:
                output += empty * 2

        #Total duplication levels (%)
        for fastqc in fastqcs:
            if fastqc.dup:
                output += [str(fastqc.dup[0])]
            else:
                output += empty * 2

        #MAPQ, joins each category 2 at a time
        if sam.mapq:
            for i in range(3):
                output += [str(float(sam.mapq[2*i][0]) + \
                           float(sam.mapq[2*i+1][0]))]
                output += [str(int(sam.mapq[2*i][1]) + \
                           int(sam.mapq[2*i+1][1]))]
        else:
            output += empty * 6

        self.short_csv.append(output)

    def table_long_(self, matched_list):
        """Appends the content of the short table for matching
        before/after/sam/meta files

        Args:
            qcbefore:Fastqc object before trimming
            qcafter: Fastqc object after trimming
            sam: Sam object
            meta: Meta object
        """
        fastqcs = [matched_list[0], matched_list[1]]
        sam = matched_list[2]
        meta = matched_list[3]
        empty = ['-']

        #File name
        output = [search_name(matched_list)]

        #nb sequences
        for fastqc in fastqcs:
            if fastqc.nb_sequences:
                output += [str(fastqc.nb_sequences)]
            else:
                output += empty

        #GC content
        for fastqc in fastqcs:
            if fastqc.gc_content:
                output += [str(fastqc.gc_content)]
            else:
                output += empty

        #Meta data
        """
        #code for MetaOld
        if meta.meta:
            output += [meta.meta['Instrument']]
            output += [meta.meta['Alias']]
            output += [meta.meta['Organism']]
            output += [meta.meta['Library name']]
            output += [meta.meta['Library strategy']]
            output += [meta.meta['Library layout']]
            output += [meta.meta['Submissions']]
        else:
            output += empty * 7
        """

        if meta.meta:
            output += [meta.meta['instrument_model']]
            output += [meta.meta['organism_ch1']]
            output += [meta.meta['title']]
            output += [meta.meta['library_strategy']]
            output += [meta.meta['characteristics_ch1']]
        else:
            output += empty * 5

        #Per base sequence quality
        for fastqc in fastqcs:
            if fastqc.pos_quality:
                for i in fastqc.pos_quality[:10]:
                    output += [str(i)]
                for i in fastqc.pos_quality[-10:]:
                    output += [str(i)]
            else:
                output += empty * 20

        #Per sequence quality scores
        for fastqc in fastqcs:
            for i in range(self.max_qual_()):
                if fastqc.qual:
                    if i < len(fastqc.qual):
                        output += [str(fastqc.qual[i])]
                    else:
                        output += ['0']
                else:
                    output += empty

        #Sequence Length Distribution
        for fastqc in fastqcs:
            for i in range(self.min_length_(), self.max_length_()+1):
                if fastqc.seq_length:
                    if i in [row[0] for row in fastqc.seq_length]:
                        temp = [row[0] for row in fastqc.seq_length]
                        temp2 = fastqc.seq_length[temp.index(i)][1]
                        output += [str(temp2)]
                    else:
                        output += ['0']
                else:
                    output += empty

        #Sequences duplication levels
        for fastqc in fastqcs:
            if fastqc.dup:
                for i in fastqc.dup:
                    output += [str(i)]
            else:
                output += empty * 17

        #MAPQ >=30, >=20, >=10, >=3, < 3, Unmapped  (nb)
        if sam.mapq:
            for i in range(6):
                output += [float(sam.mapq[i][0]), int(sam.mapq[i][1])]
        else:
            for i in range(6):
                output += empty*2

        self.long_csv.append(output)

    def make_tables(self):
        """Matches the different lists with one another, based on filename
            and generates the tables
        """
        self.header_short_()
        self.header_long_()
        for matched_list in self.matched_lists:
            self.table_short_(matched_list)
            self.table_long_(matched_list)
        return self.short_csv, self.long_csv


def search_name(matched_list):
    """Finds a filename to be used by tables
    the order is the same as the order in MatchedLists sublists

    Args:
        matched_list:
            list of the following items:
                qcbefore:Fastqc object before trimming
                qcafter: Fastqc object after trimming
                sam: Sam object
                meta: Meta object

    Returns: a filename
    """
    name = ''
    for i in range(4):
        if not name and matched_list[i].name:
            name = matched_list[i].name
    if not name:
        name = '-'
    return name


def empty_slot(nb_empty):
    """Return a list of n empty strings
    """
    return ['']*nb_empty
