"""
Created on January 25, 2014

@author: Jonathan Laperle(jonathan.laperle@usherbrooke.ca)
"""
import geffq.others as others

class Table(object):
    """Used to generate .TAB tables from fastqc/sam/meta files
    
    Attributes:
        matchedList:
        path: output path
    """
    def __init__(self, matchedLists, path):
        self.matchedLists = matchedLists
        self.path = path

    def main(self):
        """Matches the different lists with one another, based on filename
            and generates the tables
        """
        self.writeHeaderShort()
        self.writeHeaderLong()
        
        for matchedList in self.matchedLists:
            self.writeTableShort(matchedList[0], matchedList[1], matchedList[2], matchedList[3])
            self.writeTableLong(matchedList[0], matchedList[1], matchedList[2], matchedList[3])

    def maxQual(self):
        """Maxmimal quality value for all fastqc objects
        """
        return max([max([len(qc.qual) for qc in [row[0] for row in self.matchedLists]]),
                max([len(qc.qual) for qc in [row[1] for row in self.matchedLists]])])

    def minLength(self):
        """Minimal length for all fastqc objects
        """
        temp = []
        for i in range(2):
            for qc in [row[i] for row in self.matchedLists]:
                if qc.fileName:
                    temp.append(qc.seqLength[0][0])
        return max(temp)
    
    def maxLength(self):
        """Maximal length for all fastqc objects in fastqcList
        """
        temp = []
        for i in range(2):
            for qc in [row[i] for row in self.matchedLists]:
                if qc.fileName:
                    temp.append(qc.seqLength[-1][0])
        return max(temp)

    def writeHeaderShort(self):
        """Writes the header for the short table file
        """
        s = 'File Name\t'
        s += 'Library Name\t'
        s += 'nb Sequences\t'
        s += 'nb Sequences (trimmed)\t'
        s += 'GC content(%)\t'
        s += 'GC content(%) (trimmed)\t'
        s += 'Mean quality score of first 10 bases\t'
        s += 'Mean quality score of last 10 bases\t'
        s += 'Mean quality score of first 10 bases (trimmed)\t'
        s += 'Mean quality score of last 10 bases (trimmed)\t'
        s += 'mode for average quality\t'
        s += '% of sequences with the mode\t'
        s += 'mode for average quality (trimmed)\t'
        s += '% of sequences with the mode (trimmed)\t'
        s += 'mode for sequence length\t'
        s += '% of sequences with the mode\t'
        s += 'mode for sequence length (trimmed)\t'
        s += '% of sequences with the mode (trimmed)\t'
        s += '% Duplication\t'
        s += '% Duplication (trimmed)\t'
        s += 'MAPQ >= 20 (%)\t'
        s += 'MAPQ >= 20\t'
        s += '20 > MAPQ >= 3 (%)\t'
        s += '20 > MAPQ >= 3\t'
        s += 'MAPQ < 3 or unmapped (%)\t'
        s += 'MAPQ < 3 or unmapped\n'
        
        f = open(self.path + 'tableShort.tab', 'w')
        f.write(s)
        f.close()
        
    def firstHeaderLevel(self):
        """Generates the first level for the long table's header
        
        Returns:
            the header as a string
        """
        s = 'General information' + multiTabs(12)
        for _i in range(2):
            s += 'Per base sequence quality' + multiTabs(20)
        for _i in range(2):
            s += 'Number of sequences per quality score' + multiTabs(self.maxQual())
        for _i in range(2):
            s += 'Sequence Length Distribution' + multiTabs(self.maxLength() - self.minLength()+1)
        for _i in range(2):
            s += 'Sequences duplication levels' + multiTabs(11)
        s += 'MAPQ' + multiTabs(11)
        s += '\n'
        return s
    
    def secondHeaderLevel(self):
        """Generates the second level for the long table's header
        
        Returns:
            the header as a string
        """
        #General information
        s = 'File Name\t'
        s += 'nb Sequences\t'
        s += 'nb Sequences (trimmed)\t'
        s += 'GC content(%)\t'
        s += 'GC content(%) (trimmed)\t'
        s += 'Instrument\t'
        s += 'Alias\t'
        s += 'Organism\t'
        s += 'Library Name\t'
        s += 'Library Strategy\t'
        s += 'Library Layout\t'
        s += 'Submissions\t'
        
        #Per base sequence quality
        s += 'Untrimmed'+ multiTabs(20)
        s += 'Trimmed'+ multiTabs(20)
        
        #Per sequence quality scores
        s += 'Untrimmed'+ multiTabs(self.maxQual())
        s += 'Trimmed'+ multiTabs(self.maxQual())
            
        #Sequence Length Distribution
        s += 'Untrimmed\t'+ multiTabs(self.maxLength() - self.minLength())
        s += 'Trimmed\t'+ multiTabs(self.maxLength() - self.minLength())
            
        #Sequences duplication levels
        s += 'Untrimmed\t'+ multiTabs(10)
        s += 'Trimmed\t'+ multiTabs(10)
        
        #MAPQ
        s += multiTabs(12)
        s += '\n'
        
        return s
    
    def thirdHeaderLevel(self):
        """Generates the third level for the long table's header
        
        Returns:
            the header as a string
        """
        #General information
        s = multiTabs(12)
        
        #Per base sequence quality
        for i in range(2):
            for i in range(10):
                s += str(i+1) + '\t'
            for i in range(10):
                s += str(i-10) + '\t'
            
        #Per sequence quality scores
        for i in range(2):
            for i in range(self.maxQual()):
                s += str(i) + '\t'
            
        #Sequence Length Distribution
        for i in range(2):
            for i in range(self.minLength(), self.maxLength()+1):
                s += str(i) + '\t'
            
        #Sequences duplication levels
        for i in range(2):
            s += 'Total\t'
            for i in range(9):
                s += str(i+1) + '\t'
            s += '10++\t'
        
        #MAPQ
        s += "MAPQ >= 30 (%)" + '\t'
        s += "MAPQ >= 30" + '\t'
        s += "MAPQ >= 20 (%)" + '\t'
        s += "MAPQ >= 20" + '\t'
        s += "MAPQ >= 10 (%)" + '\t'
        s += "MAPQ >= 10" + '\t'
        s += "MAPQ >= 3 (%)" + '\t'
        s += "MAPQ >= 3" + '\t'
        s += "MAPQ < 3 (%)" + '\t'
        s += "MAPQ < 3" + '\t'
        s += "Unmapped (%)" + '\t'
        s += "Unmapped" + '\t'
        s += '\n'
        
        return s
    
    def writeHeaderLong(self):
        """writes the header for the long table file
        """
        others.writeFile(self.path + 'tableLong.tab', self.firstHeaderLevel()+
                                                      self.secondHeaderLevel()+
                                                      self.thirdHeaderLevel())
    
    def searchName(self, qcbefore, qcafter, sam, meta):
        """Finds a filename to be used by tables
        
        the order is the same as the order in MatchedLists sublists
        
        Args:
            qcbefore:Fastqc object before trimming
            qcafter: Fastqc object after trimming
            sam: Sam object
            meta: Meta object
        
        Returns: a filename
        """
        if qcbefore.fileName:
            name = qcbefore.fileName
        elif qcafter.fileName:
            name = qcafter.fileName
        elif sam.name:
            name = sam.name
        elif meta.name:
            name = meta.name
        else:
            name = '-'
        return name
    
    def writeTableShort(self, qcbefore, qcafter, sam, meta):
        """Appends the content of the short table for matching before/after/sam/meta files
        
        Args:
            qcbefore:Fastqc object before trimming
            qcafter: Fastqc object after trimming
            sam: Sam object
            meta: Meta object
        """
        fastqcs = [qcbefore, qcafter]
        
        #File name
        output = self.searchName(qcbefore, qcafter, sam, meta) + '\t'
        
        #Library name
        if meta.meta:
            output += meta.meta['Library name'] + '\t'
        else:
            output += '-\t'

        #nb sequences
        for fastqc in fastqcs:
            if fastqc.nbSequences:
                output += str(fastqc.nbSequences) + '\t'
            else:
                output += '-\t'

        #GC content
        for fastqc in fastqcs:
            if fastqc.gcContent:
                output += str(fastqc.gcContent) + '\t'
            else:
                output += '-\t'

        #Average base sequence quality for first and last 10 bases
        for fastqc in fastqcs:
            if fastqc.posQuality:
                output += str(others.mean(fastqc.posQuality[:10])) + '\t'
                output += str(others.mean(fastqc.posQuality[-10:])) + '\t'
            else:
                output += '-\t-\t'

        #Mode for sequence quality scores and how many sequences ahve the mode
        #in absolute and %
        for fastqc in fastqcs:
            if fastqc.qual:
                i = others.argmax(fastqc.qual)
                output += str(i) + '\t'
                output += str(fastqc.qual[i]/sum(fastqc.qual)) + '\t'
            else:
                output += '-\t-\t'

        #Mode for sequence lengthand how many sequences ahve the mode
        #in absolute and %
        for fastqc in fastqcs:
            if fastqc.seqLength:
                i = others.argmax([row[1] for row in fastqc.seqLength])
                output += str(fastqc.seqLength[i][0]) + '\t'
                output += str(fastqc.seqLength[i][1]/sum([row[1] for row in fastqc.seqLength])) + '\t'
            else:
                output += '-\t-\t'

        #Total duplication levels (%)
        for fastqc in fastqcs:
            if fastqc.dup:
                output += str(fastqc.dup[0]) + '\t'
            else:
                output += '-\t'
        
        #MAPQ, joins each category 2 at a time
        if sam.MAPQ:
            for i in range(3):
                output += str(float(sam.MAPQ[2*i][0]) + float(sam.MAPQ[2*i+1][0])) + '\t'
                output += str(int(sam.MAPQ[2*i][1]) + int(sam.MAPQ[2*i+1][1])) + '\t'
        else:
            for i in range(6):
                output += '-\t'

        others.appendFile(self.path+'tableShort.tab', output[:-1] + '\n')

    def writeTableLong(self, qcbefore, qcafter, sam, meta):
        """Appends the content of the short table for matching before/after/sam/meta files
        
        Args:
            qcbefore:Fastqc object before trimming
            qcafter: Fastqc object after trimming
            sam: Sam object
            meta: Meta object
        """
        fastqcs = [qcbefore, qcafter]
        
        #File name
        output = self.searchName(qcbefore, qcafter, sam, meta) + '\t'
        
        #nb sequences
        for fastqc in fastqcs:
            if fastqc.nbSequences:
                output += str(fastqc.nbSequences) + '\t'
            else:
                output += '-\t'

        #GC content
        for fastqc in fastqcs:
            if fastqc.gcContent:
                output += str(fastqc.gcContent) + '\t'
            else:
                output += '-\t'

        #Meta data
        if meta.meta:
            output += meta.meta['Instrument'] + '\t'
            output += meta.meta['Alias'] + '\t'
            output += meta.meta['Organism'] + '\t'
            output += meta.meta['Library name'] + '\t'
            output += meta.meta['Library strategy'] + '\t'
            output += meta.meta['Library layout'] + '\t'
            output += meta.meta['Submissions'] + '\t'
        else:
            for i in range(7):
                output += '-\t'
            
        #Per base sequence quality
        for fastqc in fastqcs:
            if fastqc.posQuality:
                for i in fastqc.posQuality[:10]:
                    output += str(i) + '\t'
                for i in fastqc.posQuality[-10:]:
                    output += str(i) + '\t'
            else:
                for i in range(20):
                    output += '-\t'

        #Per sequence quality scores
        for fastqc in fastqcs:
            for i in range(self.maxQual()):
                if fastqc.qual:
                    if i < len(fastqc.qual):
                        output += str(fastqc.qual[i]) + '\t'
                    else:
                        output += '0\t'
                else:
                    output += '-\t'
        
        #Sequence Length Distribution
        for fastqc in fastqcs:
            for i in range(self.minLength(), self.maxLength()+1):
                if fastqc.seqLength:
                    if i in [row[0] for row in fastqc.seqLength]:
                        output += str(fastqc.seqLength[[row[0] for row in fastqc.seqLength].index(i)][1]) + '\t'
                    else:
                        output += '0\t'
                else:
                    output += '-\t'

        #Sequences duplication levels
        for fastqc in fastqcs:
            if fastqc.dup:
                for i in fastqc.dup:
                    output += str(i) + '\t'
            else:
                for i in range(11):
                    output += '-\t'

        #MAPQ >=30, >=20, >=10, >=3, < 3, Unmapped  (nb)
        if sam.MAPQ:
            for i in range(6):
                output += '%s\t%s\t' % (float(sam.MAPQ[i][0]), int(sam.MAPQ[i][1]))
        else:
            for i in range(6):
                output += '-\t-\t'

        others.appendFile(self.path + 'tableLong.tab', output[:-1] + '\n')
        
def multiTabs(i):
    """Chains multiple tabs in a row as a string
    
    Args:
        i: number of tabs
    
    Returns:
        s: result string
    """
    s=''
    for i in range(i):
        s += '\t'
    return s