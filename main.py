"""
Created on May 8, 2014

@author: Jonathan Laperle(jonathan.laperle@usherbrooke.ca)
"""
from geffq.fastqc import Fastqc
from geffq.sam import Sam
from geffq.meta import Meta
from geffq.table import Table
from geffq.graph import GraphMaker
import os
import sys
import getopt

def loadFromInputFile(inputFile):
    """Loads objects from the filepaths provided in inputFile
    
    Args:
        inputFile: the path to the inputFile
            the inputfile must be of format:
            fastqc\tfastqc\tsam\tmeta\n
            ...
            with each line representing a sequencing
            misisng files must be named N\A
            fastqc can be a zip archive but the file must than be named 
            fastqc_data.txt
            
    Returns:
        A list of lists with each inner list being of format
        (fastqc, fastqc, sam, meta), each of those being files for the same
        sequencing
    """
    matchedList = []
    with open(inputFile, 'r') as f:
        for line in f:
            line = line.split()
            fastqcs = [Fastqc(), Fastqc()]
            for i in range(2):
                if line[i] != 'N/A':
                    if line[i].endswith('.zip'):
                        fastqcs[i].loadFromZip(line[i])
                    else:
                        fastqcs[i].loadFromFile(line[i])
            sam = Sam()
            if line[2] != 'N/A':
                sam.loadFromFile(line[2])
            meta = Meta(line[3].split('/')[-1])
            if line[3] != 'N/A':
                meta.loadFromFile(line[3])
            matchedList.append([fastqcs[0], fastqcs[1], sam, meta])
    return matchedList

def hasFastqc(fastqcList):
    """Verifies that at least one fastqc object is not empty
    Args:
        fastqcList: list of fastqc objects
        
    Returns: 
        True if a fastqc was found in the list
    """
    for fastqc in fastqcList:
        if fastqc.fileName:
            return True
    return False

def prepareOutputDir(outputPath):
    """Creates directories to hold the tables and graphs
    
    Args:
        outputPath: path to a directory where to output will be written
    """
    for directory in [outputPath + 'output/', outputPath + 'output/before/', outputPath + 'output/after/']:
        if not os.path.exists(directory):
            os.makedirs(directory)

def usage():
    """Prints the proper usage for main.py
    """
    print('Usage: python main.py -i <inputfile> -o <outputpath>')
    
def noInput():
    """Warns about empty input file
    """
    print('Input file is empty.')
    
def main(argv):
    """Takes 2 arguments, -i and -o form the command line and calls the proper 
        functions to generate the output
        
        -i: the path to the inputFile
            the inputfile must be of format:
            fastqc\tfastqc\tsam\tmeta\n
            ...
            with each line representing a sequencing
            misisng files must be named N\A
            fastqc can be a zip archive but the file must than be named 
            fastqc_data.txt
        -o: path to a directory where to output will be written
        --help(-h): prints proper usage syntax
    """
    inputfile = ''
    outputPath = ''
    try:
        opts, _args = getopt.getopt(argv,"hi:o:", ["help"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ["-h", "--help"]:
            usage()
            sys.exit()
        elif opt == "-i":
            inputfile = arg
        elif opt == "-o":
            outputPath = arg
    if not inputfile and outputPath:
        usage()
        sys.exit(2)

    prepareOutputDir(outputPath)
    matchedLists = loadFromInputFile(inputfile)
    if not matchedLists:
        noInput()
        sys.exit(2)

    hasNonTrimmed = hasFastqc([row[0] for row in matchedLists])
    hasTrimmed = hasFastqc([row[1] for row in matchedLists])
    
    if hasNonTrimmed or hasTrimmed:
        Table(matchedLists, outputPath + 'output/').main()
    else:
        print('No valid fastqc file, could not produce tables')
    if hasNonTrimmed:
        GraphMaker([row[0] for row in matchedLists], outputPath + 'output/before/').generateAll()
    else:
        print('No valid untrimmed fastqc file, some graphs could not be produced')
    if hasTrimmed:
        GraphMaker([row[1] for row in matchedLists], outputPath + 'output/after/').generateAll()
    else:
        print('No valid trimmed fastqc file, some graphs could not be produced')

if __name__ == '__main__':
    main(sys.argv[1:])