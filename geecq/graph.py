"""
Created on January 25, 2014

@author: Jonathan Laperle(jonathan.laperle@usherbrooke.ca)
"""
import geecq.others as others
from itertools import product
import string
import subprocess
import os
import os.path

class GraphMaker(object):
    """Used to write R scripts in order to generate graphs

    the graphs are in .png format

    Attributes:
        labels: list of 2 letter strings in alphabetical order,
            used as labels in R data-frames
        fastqcList: list of all fastqc objects to be included in the graphs
        path: output path
    """
    def __init__(self, fastqc_list, path):
        self.labels = [''.join(word) for word in
                       product(string.ascii_lowercase, string.ascii_lowercase)]
        self.fastqc_list = []
        for fastqc in fastqc_list:
            if fastqc.name:
                self.fastqc_list.append(fastqc)
        self.path = path

    def max_qual(self):
        """Maxmimal quality value for all fastqc objects in fastqcList

        Returns:
            Maxmimal quality value for all fastqc objects in fastqcList

        """
        return max([len(qc.qual) for qc in self.fastqc_list])

    def min_length(self):
        """Minimal length for all fastqc objects in fastqcList

        Returns:
            Minimal length for all fastqc objects in fastqcList
        """
        return min([qc.seq_length[0][0] for qc in self.fastqc_list])

    def max_length(self):
        """Maximal length for all fastqc objects in fastqcList

        Returns:
            Maximal length for all fastqc objects in fastqcList
        """
        return max([qc.seq_length[-1][0] for qc in self.fastqc_list])

    def generate_graph(self, file_name, rscript):
        """Generates the graph with an rscript

        Args:
            fileName: name of the file to create(no extention)
            rscript: text of the rscript used to generate a graph
        """
        others.write_file(self.path + file_name + '.R', rscript)
        subprocess.call(['Rscript', self.path + file_name + '.R'])
        os.remove(self.path + file_name + '.R')

    def make_per_base_qual_graph(self):
        """Generates the per base quality graph
        """
        script = 'data<-data.frame(\n\t'
        for x_value in range(10):
            script += "'%s' = c(" % (self.labels[x_value])
            for y_value in [fastqc.pos_quality for fastqc in self.fastqc_list]:
                script += '%g, ' % (y_value[x_value])
            script = script[:-2] + '),\n\t'
        for x_value in range(10):
            script += "'%s' = c(" % (self.labels[x_value+10])
            for y_value in [fastqc.pos_quality for fastqc in self.fastqc_list]:
                script += '%g, ' % (y_value[x_value-10])
            script = script[:-2] + '),\n\t'
        script = script[:-3] + '\n)\n'
        script += 'data<-data[,order(names(data))]\n'
        script += "png('" + self.path + "per_base_quality.png')\n"
        script += "boxplot(data, ylab='Quality score', xlab='Position', " + \
                  "names=c('1','2','3','4','5','6','7','8','9','10'," + \
                  "'-10','-9','-8','-7','-6','-5','-4','-3','-2','-1'))\n"
        script += "dev.off()"
        self.generate_graph('per_base_quality', script)

    def make_per_sequence_qual_graph(self):
        """Generates the Per sequence quality scores graph
        """
        script = 'data<-data.frame('
        for i in range(self.max_qual()+1):
            script += "'%s' = c(" % (self.labels[i])
            for fastqc in self.fastqc_list:
                if len(fastqc.qual) > i:
                    script += '%g, ' % (fastqc.qual[i] / fastqc.nb_sequences)
                else:
                    script += '0, '
            script = script[:-2] + '),\n\t'
        script = script[:-3] + '\n)\n'

        script += 'data<-data[,order(names(data))]\n'
        script += "png('"+ self.path +"per_sequence_quality.png')\n"
        script += "boxplot(data, ylab='nb sequences / total sequences', " + \
                  "xlab='Quality score', names=c("
        for i in range(self.max_qual()+1):
            script += "'%s', " % (i)
        script = script[:-2] + "))\n"
        script += "dev.off()"
        self.generate_graph('per_sequence_quality', script)

    def make_seq_len_graph(self):
        """Generates the Sequence Length Distribution graph
        """
        script = 'data<-data.frame('
        for i in range(self.min_length(), self.max_length()+1):
            script += "'%s' = c(" % (self.labels[i])
            for fastqc in self.fastqc_list:
                if i in [row[0] for row in fastqc.seq_length]:
                    temp = fastqc.seq_length[i-fastqc.seq_length[0][0]][1]
                    script += '%g, ' % (temp / fastqc.nb_sequences)
                else:
                    script += '0, '
            script = script[:-2] + '),\n\t'
        script = script[:-3] + '\n)\n'

        script += 'data<-data[,order(names(data))]\n'
        script += "png('"+self.path+"sequence_length.png')\n" + \
                  "boxplot(data, ylab='nb sequences / total sequences', " + \
                  "xlab='Length(nucl)', names=c("
        for i in range(self.min_length(), self.max_length()+1):
            script += "'%s', " % (i)
        script = script[:-2] + "))\n"
        script += "dev.off()"
        self.generate_graph('sequence_length', script)

    def make_duplication_graph(self):
        """Generates the Sequence Duplication Levels graph
        """
        fastq_version = self.fastqc_list[0].version
        if fastq_version >= 11 and not True:
            xlabels = "c('0','1','2','3','4','5','6','7','8','9','>10','>50','>100','>500', '>1k', '>5k', '>10k+')"
            length = 17
        else:
            xlabels = "c('0','1','2','3','4','5','6','7','8','9','10+')"
            length = 11
        script = 'data<-data.frame('
        for x_value in range(length):
            script += "'%s' = c(" % (self.labels[x_value])
            for y_value in [fastqc.dup for fastqc in self.fastqc_list]:
                script += '%g, ' % (y_value[x_value])
            script = script[:-2] + '),\n\t'
        script = script[:-3] + '\n)\n'
        script += 'data<-data[,order(names(data))]\n'
        script += "png('"+self.path+"sequence_duplication.png')\n"
        script += "boxplot(data, ylab='nb sequences(%)', " + \
                  "xlab='Sequence duplication level', " + \
                  "names=" + xlabels + ")\n"
        script += "dev.off()"
        self.generate_graph('sequence_duplication', script)

    def generate_all(self):
        """Generates all graphs in path specified by self.path
        """
        self.make_per_base_qual_graph()
        self.make_per_sequence_qual_graph()
        self.make_seq_len_graph()
        self.make_duplication_graph()
