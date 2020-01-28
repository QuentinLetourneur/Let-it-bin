#! /usr/bin/env python
# -*- coding: utf8 -*-

"""
count_matrix.py : Script that merge several count table into a count matrix
                  Created count matrix contains 
                  in columns: 
                     1) reference sequence 
                     2) sequence length  
                     3) number of read mapped to the reference for sample 1
                     4) ...
                     5) ...    other samples.
"""

__author__ = "Anita Annamalé"
__version__  = "2.0"
__copyright__ = "copyleft"
__date__ = "2016/07"

#-------------------------- MODULES IMPORTATION -------------------------------#

import sys
import argparse
import re
import os.path
import glob

#----------------------------- DICTIONNARIES ----------------------------------#

MATRIX = {}
COUNT = {}

#-------------------------- FUNCTIONS DEFINITION ------------------------------#

def read(filename, nb_samples, id_sample):
    """
    Function that read a count table 

    Args: 
        filename [STR] = name of the count table
        nb_samples [INT] = total number of samples
        id_sample [INT] = id of the count table
    
    No returns.
    """
    with open(filename, "rt") as f:
        end_re = re.compile("^\*")
        for line in f:
            match = end_re.search(line)
            if match:
                continue
            # separate columns
            columns = line.split()
            # get species name
            ref = columns[0]
            # add the ref into global dictionnary MATRIX & COUNT
            if not MATRIX.has_key(ref):
                size = columns[1]
                MATRIX[ref] = size
                COUNT[ref] = ['0'] * nb_samples
            # get number of mapped reads and add to count table
            count = columns[2]
            COUNT[ref][id_sample] = count


def write(output):
    """
    Function that write the count matrix

    Args:
        output [STR] = output file name

    No returns
    """
    with open(output, "w") as o:
        o.write('ref_genome\tSize\t{0}\n'.format('\t'.join(filenames))) # header
        for key,values in COUNT.items():
            o.write('{0}\t{1}\t{2}\n'.format(key, MATRIX[key],'\t'.join(values)))


#---------------------------- PROGRAMME MAIN ----------------------------------#

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description='Create count matrix.')
    parser.add_argument("-i","--input", 
                        type=str,
                        action='store', 
                        nargs='+',
                        help='files containing count table, to be merged')
    parser.add_argument("-d","--directory",
                        type=str, 
                        action='store', 
                        help='directory containing files (count tables), to be\
 merged')
    parser.add_argument("-o", "--output",
                        type=str,
                        action='store',
                        default='count_matrix.txt',
                        help='name of the output file')

    # command-line parsing
    if not len(sys.argv) >= 2:
        sys.exit("usage: python parse_compt_genome.py [-i <file1.txt> \
<file2.txt>] [-d <directory>] [-o <path_to_directory>]")

    # parse commandline arguments
    arguments = parser.parse_args()
    arguments = dict(arguments._get_kwargs())

    # get the all count table
    if arguments['input']:
        files = os.path.abspath(arguments['input'])

    if arguments['directory']:
        path = os.path.abspath(arguments['directory'])
        if os.path.exists(path):
            files = [ file1 for file1 in glob.glob(path + "/*.txt") 
                      if os.stat(file1).st_size != 0 ] # all the file ending by
                                                       # .txt in the directory

    # merging files
    filenames=[]  # label list of samples, with the extension '.txt'
    nb_tables = len(files) # total number of file to merge

    for i in xrange(nb_tables):
        # get sample name
        table = files[i]
        base=os.path.basename(table)
        filenames.append(base[:-4])
        # get information from the table
        read(table, nb_tables, i) 

    # write the count matrix
    output = os.path.abspath(arguments['output']) # file name
    write(output) 


