"""
Created on May 8, 2014

@author: Jonathan Laperle(jonathan.laperle@usherbrooke.ca)
"""
from geecq.fastqc import Fastqc
from geecq.sam import Sam
from geecq.meta import Meta
from geecq.table import Table
from geecq.graph import GraphMaker
import os
import sys
import getopt
import csv

def load_from_input(input_file):
    r"""Loads objects from the filepaths provided in inputFile

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
    input_matrix = []
    with open(input_file, 'r') as in_file:
        for line in in_file:
            line = line.split()
            if len(line) != 4:
                print "Error: Every line should contain 4 elements"
                exit(0)
            fastqcs = [Fastqc(), Fastqc()]
            for i in range(2):
                if line[i] != 'N/A':
                    if line[i].endswith('.zip'):
                        fastqcs[i].load_from_zip(line[i])
                    else:
                        fastqcs[i].load_from_file(line[i])
            sam = Sam()
            if line[2] != 'N/A':
                sam.load_from_file(line[2])
            meta = Meta(line[3].split('/')[-1])
            if line[3] != 'N/A':
                meta.load_from_file(line[3])
            input_matrix.append([fastqcs[0], fastqcs[1], sam, meta])
    verify_input_matrix(input_matrix)
    return input_matrix

def verify_input_matrix(input_matrix):
    """Gives an error message and terminates if input has no valid file

    Args:
        input_matrix: input file in matrix format
    """
    if not input_matrix:
        no_input()
        sys.exit(2)

def has_fastqc(fastqc_list):
    """Verifies that at least one fastqc object is not empty
    Args:
        fastqcList: list of fastqc objects

    Returns:
        True if a fastqc was found in the list
    """
    for fastqc in fastqc_list:
        if fastqc.name:
            return True
    return False

def prepare_output_dir(output_path):
    """Creates directories to hold the tables and graphs

    Args:
        outputPath: path to a directory where to output will be written
    """
    for directory in [output_path + 'output/',
                      output_path + 'output/before/',
                      output_path + 'output/after/']:
        if not os.path.exists(directory):
            os.makedirs(directory)

def usage():
    """Prints the proper usage for main.py
    """
    print 'Usage: python main.py -i <inputfile> -o <outputpath>'

def no_input():
    """Warns about empty input file
    """
    print 'Input file is empty.'

def write_csv(table, path):
    csv_file = csv.writer(open(path, 'w'), dialect='excel-tab')
    for line in table:
        csv_file.writerow(line)

def launch(input_matrix, output_path):
    """Does multiple checks to ensure the list has the nessessary data
    and launches the modules to produce output

    Args:
        input_matrix: input file in matrix format
    """
    has_ntrimmed = has_fastqc([row[0] for row in input_matrix])
    has_trimmed = has_fastqc([row[1] for row in input_matrix])

    output_path = output_path + 'output/'

    if has_ntrimmed or has_trimmed:
        table_short, table_long = Table(input_matrix).make_tables()
        write_csv(table_short, output_path + 'tableShort.tab')
        write_csv(table_long, output_path + 'tableLong.tab')
    else:
        print 'No valid fastqc file, could not produce tables'
    if has_ntrimmed:
        GraphMaker([row[0] for row in input_matrix],
                   output_path + 'before/').generate_all()
    else:
        print 'No valid untrimmed fastqc file, ' \
              'some graphs will not be produced'
    if has_trimmed:
        GraphMaker([row[1] for row in input_matrix],
                   output_path + 'after/').generate_all()
    else:
        print 'No valid trimmed fastqc file, ' \
              'some graphs will not be produced'


def main(argv):
    r"""Takes 2 arguments, -i and -o form the command line and calls the proper
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
    input_file = ''
    output_path = ''
    try:
        opts, _ = getopt.getopt(argv, "hi:o:", ["help"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ["-h", "--help"]:
            usage()
            sys.exit()
        elif opt == "-i":
            input_file = arg
        elif opt == "-o":
            output_path = arg
    if not input_file and output_path:
        usage()
        sys.exit(2)

    prepare_output_dir(output_path)
    input_matrix = load_from_input(input_file)
    launch(input_matrix, output_path)

if __name__ == '__main__':
    main(sys.argv[1:])
