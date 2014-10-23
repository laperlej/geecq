"""
Created on January 25, 2014

@author: Jonathan Laperle(jonathan.laperle@usherbrooke.ca)
"""
import geffq.others as others
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
    def __init__(self, fastqcList, path):
        self.labels = [''.join(word) for word in product(string.ascii_lowercase, string.ascii_lowercase)]
        self.fastqcList = []
        for fastqc in fastqcList:
            if fastqc.fileName:
                self.fastqcList.append(fastqc)
        self.path = path
    
    def maxQual(self):
        """Maxmimal quality value for all fastqc objects in fastqcList
        
        Returns:
            Maxmimal quality value for all fastqc objects in fastqcList
            
        """
        return max([len(qc.qual) for qc in self.fastqcList])
    
    def minLength(self):
        """Minimal length for all fastqc objects in fastqcList
        
        Returns:
            Minimal length for all fastqc objects in fastqcList
        """
        return min([qc.seqLength[0][0] for qc in self.fastqcList])
    
    def maxLength(self):
        """Maximal length for all fastqc objects in fastqcList
        
        Returns:
            Maximal length for all fastqc objects in fastqcList
        """
        return max([qc.seqLength[-1][0] for qc in self.fastqcList])
    
    def generateGraph(self, fileName, rscript):
        """Generates the graph with an rscript
        
        Args:
            fileName: name of the file to create(no extention)
            rscript: text of the rscript used to generate a graph
        """
        others.writeFile(self.path + fileName + '.R', rscript)
        subprocess.call(['Rscript', self.path + fileName + '.R'])
        os.remove(self.path + fileName + '.R')
        
    def makePerBaseQualGraph(self):
        """Generates the per base quality graph
        """
        s = 'data<-data.frame(\n\t'
        for x in range(10):
            s += "'%s' = c(" % (self.labels[x])
            for y in [fastqc.posQuality for fastqc in self.fastqcList]:
                s += '%g, ' % (y[x])
            s = s[:-2] + '),\n\t'
        for x in range(10):
            s += "'%s' = c(" % (self.labels[x+10])
            for y in [fastqc.posQuality for fastqc in self.fastqcList]:
                s += '%g, ' % (y[x-10])
            s = s[:-2] + '),\n\t'
        s = s[:-3] + '\n)\n'
        s += 'data<-data[,order(names(data))]\n'
        s += "png('"+ self.path +"per_base_quality.png')\nboxplot(data, ylab='Quality score', xlab='Position', names=c('1','2','3','4','5','6','7','8','9','10','-10','-9','-8','-7','-6','-5','-4','-3','-2','-1'))\ndev.off()"
        self.generateGraph('per_base_quality', s)

    def makePerSequenceQualGraph(self):
        """Generates the Per sequence quality scores graph
        """
        s = 'data<-data.frame('
        for i in range(self.maxQual()+1):
            s += "'%s' = c(" % (self.labels[i])
            for fastqc in self.fastqcList:
                if len(fastqc.qual) > i:
                    s += '%g, ' % (fastqc.qual[i] / fastqc.nbSequences)
                else:
                    s += '0, '
            s = s[:-2] + '),\n\t'
        s = s[:-3] + '\n)\n'
        
        s += 'data<-data[,order(names(data))]\n'
        s += "png('"+ self.path +"per_sequence_quality.png')\nboxplot(data, ylab='nb sequences / total sequences', xlab='Quality score', names=c("
        for i in range(self.maxQual()+1):
            s += "'%s', " % (i)
        s = s[:-2] + '))\ndev.off()'
        self.generateGraph('per_sequence_quality', s)

    def makeSeqLenGraph(self):
        """Generates the Sequence Length Distribution graph
        """
        s = 'data<-data.frame('
        for i in range(self.minLength(), self.maxLength()+1):
            s += "'%s' = c(" % (self.labels[i])
            for fastqc in self.fastqcList:
                if i in [row[0] for row in fastqc.seqLength] :
                    s += '%g, ' % (fastqc.seqLength[i-fastqc.seqLength[0][0]][1] / fastqc.nbSequences)
                else:
                    s += '0, '
            s = s[:-2] + '),\n\t'
        s = s[:-3] + '\n)\n'
        
        s += 'data<-data[,order(names(data))]\n'
        s += "png('"+self.path+"sequence_length.png')\nboxplot(data, ylab='nb sequences / total sequences', xlab='Length(nucl)', names=c("
        for i in range(self.minLength(), self.maxLength()+1):
            s += "'%s', " % (i)
        s = s[:-2] + '))\ndev.off()'
        self.generateGraph('sequence_length', s)

    def makeDuplicationGraph(self):
        """Generates the Sequence Duplication Levels graph
        """
        s = 'data<-data.frame('
        for x in range(11):
            s += "'%s' = c(" % (self.labels[x])
            for y in [fastqc.dup for fastqc in self.fastqcList]:
                s += '%g, ' % (y[x])
            s = s[:-2] + '),\n\t'
        s = s[:-3] + '\n)\n'
        s += 'data<-data[,order(names(data))]\n'
        s += "png('"+self.path+"sequence_duplication.png')\nboxplot(data, ylab='nb sequences(%)', xlab='Sequence duplication level', names=c('0','1','2','3','4','5','6','7','8','9','10+'))\ndev.off()"
        self.generateGraph('sequence_duplication', s)
        
    def generateAll(self):
        """Generates all graphs in path specified by self.path
        """
        self.makePerBaseQualGraph()
        self.makePerSequenceQualGraph()
        self.makeSeqLenGraph()
        self.makeDuplicationGraph()