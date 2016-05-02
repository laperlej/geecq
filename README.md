# geecq

Takes 2 arguments, -i and -o form the command line and calls the proper
functions to generate the output

-i: the path to the inputFile
the inputfile must be of format:
fastqc\tfastqc\tsamstat\tmeta\n
...
with each line representing a sequencing
misisng files must be named N\A
fastqc can be a zip archive but the file must than be named
fastqc_data.txt
-o: path to a directory where to output will be written
--help(-h): prints proper usage syntax
