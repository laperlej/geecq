"""
Created on January 25, 2014

@author: Jonathan Laperle(jonathan.laperle@usherbrooke.ca)
"""
import csv
import geffq.others as others

class Table(object):
    """Used to generate .TAB tables from fastqc/sam/meta files

    Attributes:
        matchedList:
        path: output path
    """
    def __init__(self, matchedLists, out_path):
        self.matched_lists = matchedLists
        self.out_path = out_path
        self.short_csv = self.open_csv('tableShort.tab')
        self.long_csv = self.open_csv('tableLong.tab')

    def open_csv(self, file_name):
        """opens a file and returns a csv writer
        """
        return csv.writer(open(self.out_path + file_name, 'w'),
                          dialect='excel-tab')

    def max_qual(self):
        """Maxmimal quality value for all fastqc objects
        """
        return max([max([len(qc.qual) for qc in \
            [row[0] for row in self.matched_lists]]), \
            max([len(qc.qual) for qc in \
            [row[1] for row in self.matched_lists]])])

    def min_length(self):
        """Minimal length for all fastqc objects
        """
        temp = []
        for i in range(2):
            for fastqc in [row[i] for row in self.matched_lists]:
                if fastqc.name:
                    temp.append(fastqc.seq_length[0][0])
        return min(temp)

    def max_length(self):
        """Maximal length for all fastqc objects in fastqcList
        """
        temp = []
        for i in range(2):
            for fastqc in [row[i] for row in self.matched_lists]:
                if fastqc.name:
                    temp.append(fastqc.seq_length[-1][0])
        return max(temp)

    def write_header_short(self):
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

        self.short_csv.writerow(short_header)

    def write_header_long(self):
        """writes the header for the long table file
        """
        self.long_csv.writerow(self.first_header_level())
        self.long_csv.writerow(self.second_header_level())
        self.long_csv.writerow(self.third_header_level())

    def first_header_level(self):
        """Generates the first level for the long table's header

        Returns:
            the header as a string
        """
        first_header = []
        first_header.append('General information')
        first_header += empty_slot(10)
        for _ in range(2):
            first_header.append('Per base sequence quality')
            first_header += empty_slot(19)
        for _ in range(2):
            first_header.append('Number of sequences per quality score')
            first_header += empty_slot(self.max_qual()-1)
        for _ in range(2):
            first_header.append('Sequence Length Distribution')
            first_header += empty_slot(self.max_length() - self.min_length())
        for _ in range(2):
            first_header.append('Sequences duplication levels')
            first_header += empty_slot(10)
        first_header.append('MAPQ')
        first_header += empty_slot(11)
        return first_header

    def second_header_level(self):
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
                         'Characteristics',
                         'Submissions']
        #in order:
        #Per base sequence quality
        #Per sequence quality scores
        #Sequence Length Distribution
        #Sequences duplication levels
        for nb_empty in [19,
                         self.max_qual() - 1,
                         self.max_length() - self.min_length(),
                         10]:
            second_header += ['Untrimmed'] + empty_slot(nb_empty)
            second_header += ['Trimmed'] + empty_slot(nb_empty)

        #MAPQ
        second_header += empty_slot(12)

        return second_header

    def third_header_level(self):
        """Generates the third level for the long table's header

        Returns:
            the header as a string
        """
        third_header = []

        #General information
        third_header += empty_slot(11)

        #Per base sequence quality
        for _ in range(2):
            third_header += [str(i+1) for i in range(10)]
            third_header += [str(i-10) for i in range(10)]

        #Per sequence quality scores
        for _ in range(2):
            third_header += [str(i) for i in range(self.max_qual())]

        #Sequence Length Distribution
        for _ in range(2):
            third_header += [str(i) for i in range(self.min_length(),
                                                   self.max_length()+1)]

        #Sequences duplication levels
        for i in range(2):
            third_header.append('Total')
            third_header += [str(i+1) for i in range(9)]
            third_header.append('10++')

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

    def write_table_short(self, matched_list):
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
            output.append(meta.meta['Library name'])
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

        self.short_csv.writerow(output)

    def write_table_long(self, matched_list):
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
            output += [meta.meta['Submissions']]
        else:
            output += empty * 6

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
            for i in range(self.max_qual()):
                if fastqc.qual:
                    if i < len(fastqc.qual):
                        output += [str(fastqc.qual[i])]
                    else:
                        output += ['0']
                else:
                    output += empty

        #Sequence Length Distribution
        for fastqc in fastqcs:
            for i in range(self.min_length(), self.max_length()+1):
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
                output += empty * 11

        #MAPQ >=30, >=20, >=10, >=3, < 3, Unmapped  (nb)
        if sam.mapq:
            for i in range(6):
                output += [float(sam.mapq[i][0]), int(sam.mapq[i][1])]
        else:
            for i in range(6):
                output += empty*2

        self.long_csv.writerow(output)

    def main(self):
        """Matches the different lists with one another, based on filename
            and generates the tables
        """
        self.write_header_short()
        self.write_header_long()

        for matched_list in self.matched_lists:
            self.write_table_short(matched_list)
            self.write_table_long(matched_list)


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
