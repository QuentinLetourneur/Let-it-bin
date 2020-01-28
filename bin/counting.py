#! /usr/bin/env python
# -*- coding: utf8 -*-

"""
counting.py : Script that count alignment from a bam file using three different
              methods (best, ex-aequo, shared) and output a count table.

              Mapping generates two sort of reads, reads that map to a unique
              location with the best alignment score, called "unique reads"
              and reads those map to several locations with the best alignment
              score, called as "multiple reads"
              The Three counting methods differ by the way they handle multiple
              reads :
                - Best : a multiple read provide 1 point to one location randomly
                - Ex-aequo : a multiple read provide 1 point to each location
                - Shared : a multiple read is weights according the probability
                            that the alignment is the true point of origin of
                            the read

                Created count table contains
                in columns:
                     1) reference sequence
                     2) sequence length
                     3) number of read mapped to the reference for sample
"""

__author__ = "Anita Annamalé"
__version__  = "2.0"
__copyright__ = "copyleft"
__date__ = "2016/07"

#-------------------------- MODULES IMPORTATION -------------------------------#

import sys
import re
import os
import pysam
import argparse
from collections import defaultdict
import itertools


#-------------------------- FUNCTIONS DEFINITION ------------------------------#

def is_sam(path):
    """
    Check if a path is an existing SAM file.

    Args:
        path [STR] = Path to the file

    Returns:
        abs_path [STR] = absolute path
        or quit
    """
    abs_path = os.path.abspath(path) # get absolute path

    # if not a file
    if not os.path.isfile(abs_path):
        # check if it is a path to a dir
        if os.path.isdir(abs_path):
            msg = "{0} is a directory not a file.".format(abs_path)
        # else the path doesn't not exist
        else:
            msg = "The path {0} does not exist.".format(abs_path)
        raise argparse.ArgumentTypeError(msg)

    else :
        ext = os.path.splitext(abs_path)[1] # get the file extension
        if ext != '.sam':
            msg = ("{0} isn't a sam file. "
                  "Please, provide a file with extension .sam".format(abs_path))

    return abs_path


def path_not_exists(path):
    """
    Check if a path is exist or not.

    Args:
        path [string] = Path to check

    Returns:
        abs_path [string] = absolute path if it doesn't exist
        or quit
    """
    abs_path = os.path.abspath(path)

    # if path exists
    if os.path.exists(abs_path):
        msg = "The path {0} already exist. "
        "Please give another file name".format(abs_path)
        raise argparse.ArgumentTypeError(msg)

    return abs_path


def create_bam(filename, threads=0):
    """
    Function that create a BAM file from a SAM file.

    Args :
        filename [STR] = SAM filename

    Returns:
        bamfile [STR] = BAM filename
    """
    # name of the bam file to create
    bamfile = os.path.dirname(filename)[:-3] + "bam/" + os.path.basename(filename)[:-3] + "bam"
    # convert sam to bam using pysam
    pysam.view('-@', str(threads - 1), '-S', '-b', '-o', bamfile, filename, catch_stdout=False)
    os.remove(filename)

    return bamfile


def sort_bam(filename, memory, threads=0):
    """
    Function that sort a BAM file using pysam

    Args:
        filename [STR] = unsorted BAM file

    Returns :
        sortfile [STR] = sorted BAM filename
    """
    sortfile = '{0}/sorted_{1}'.format(os.path.dirname(filename),
                                       os.path.basename(filename)[:-4])
                                   
    if threads > 0:
        mem_per_thread = memory / threads
        # print "with", mem_per_thread, "M memory per thread"
        # sort the bam file
        pysam.sort('-o', sortfile + '.bam', '-@', str(threads - 1), '-m', str(mem_per_thread) + 'M', '-O', 'bam', '-T', sortfile, filename)
    else:
        pysam.sort('-o', sortfile + '.bam', '-m', str(memory) + 'M', '-O', 'bam', '-T', sortfile, filename)
    
    return sortfile + ".bam"

    
def only_mapped(filename, threads=1):
    """
    Function that keep only mapped reads in the BAM file

    Args:
        filename [STR] = BAM file, containing all alignments

    Returns:
        mappedfile [STR] = BAM file with only the mapped reads
    """
    # new bamfile name
    mappedfile = '{0}/filtered_{1}'.format(os.path.dirname(filename),
                                       os.path.basename(filename)[7:])
    # get only mapped reads
    pysam.view('-@', str(threads - 1), '-b', '-F', '4', filename, '-o', mappedfile, catch_stdout=False)
    
    return mappedfile


def write_table(filename, outfile):
    """
    Function that create a count table using pysam. First index the BAM file,
    then count reads using the function idxstats from pysam, and output a count
    table.

    Args :
        filename [STR] : BAM file to count
        outfile [STR] : count table name

    No Returns
    """
    # index the bam file
    pysam.index(filename)
    # create count table
    table = pysam.idxstats(filename)
    # write the count table
    with open(outfile, 'wt') as out:
        for line in table:
            out.write(line)


def read_bam(filename, count, bam=False):
    """
    Function that count reads from a BAM file, using the given methods in count:
    "ex-aequo" or "shared". If bam is set to 'True', a BAM file containing only
    the used alignment for the counting is generated.

    Args:
        filename [STR] = BAM file to count
        count [STR] = counting method, either 'ex-aequo', either 'shared'
        bam [True/False] = create or not a BAM file

    Returns:
        if count = ex-aequo : new_filename [STR] = BAM file
        if count = shared :
                database [DICT] = contains length of reference genomes.
                                    key : reference genome
                                    value : size
                tmp_genomes [DICT] =
    """
    # initialise
    tmp_score = dict()
    tmp_genomes = defaultdict(list)
    # open the BAM file
    bamfile = pysam.AlignmentFile(filename, "rb")
    
    if (count == "ex-aequo") or bam:
        reads = defaultdict(list) # contains a list a alignment of each read
        # name of the filtered BAM file
        new_filename = (os.path.dirname(filename) + "/unsorted_filtered_" +
                       os.path.basename(filename))
    if count == "shared":
        # create a dictionary containing the length of reference genomes
        references = bamfile.references # get name of reference sequence
        lengths = bamfile.lengths   # get reference length
        database = dict(zip(references, lengths))
    
    for element in bamfile:
        if not element.is_unmapped: # if read mapped
            if element.has_tag('AS'):  # if an alignment score is calculated
                if element.is_read1:
                    read_id = '{0}_1'.format(element.qname) # get read direction
                else :
                    read_id = '{0}_2'.format(element.qname)
                score = element.get_tag('AS') # get alignment score
                prev_score = tmp_score.get(read_id, score) # get previous score
                                                           # for the read
                # if same score
                if prev_score == score:
                    tmp_score[read_id] = score
                    # add the genome to the list if it doesn't exist
                    if element.reference_name not in tmp_genomes[read_id]:
                        tmp_genomes[read_id].append(element.reference_name)
                        if (count == "ex-aequo") or bam:
                            reads[read_id].append(element) # append line for
                                                           # filtered BAM file
                # if previous score is lower
                elif prev_score < score:
                    tmp_score[read_id] = score # set the new score
                    # reinitialise all
                    tmp_genomes[read_id] = [element.reference_name]
                    if (count == "ex-aequo") or bam:
                        reads[read_id] = [element]
            else :
                sys.exit("[FATAL error] Parsing sam file : no optional 'AS' field.\n")

    bamfile.close()

    # writing a bam file for ex-aequo counting and if a bam=True
    if (count == "ex-aequo") or bam:
        read_list = list(itertools.chain(reads.values()))
        merged_list = list(itertools.chain.from_iterable(read_list))
        # writing the filtered BAM file
        ex_aequo_reads = pysam.AlignmentFile(new_filename, "wb", template=bamfile)
        for element in merged_list:
            ex_aequo_reads.write(element)
        ex_aequo_reads.close()
        # return the filename, once written
        if count == "ex-aequo":
            return new_filename
        # if bam, sort the BAM file
        if bam:
            sortfile = os.path.dirname(filename) + "/filtered_" + os.path.basename(filename)[:-4]
            pysam.sort(new_filename, sortfile)
    
    # returns reads and alignements informations if shared counting
    if (count =="shared") or bam:
        for (key,values) in tmp_genomes.iteritems():
            tmp_genomes[key] = list(set(values))
        return database, tmp_genomes


def uniq_from_mult(genome_dict, unique_dict):
    """
    Function that filter unique reads from all reads. Multiple reads are
    reads that map to more than one genome. And Unique reads are reads that map
    only on one genome.

    Args:
        genome_dict [DICT] = dictionary containing reads as key and a list of
                             genome where read mapped as value
        unique_dict [DICT] = contains for each reference genome the number of
                             unique reads

    Returns:
        genome_dict [DICT] = the dictionary without unique reads
        unique_dict [DICT] = nb of unique read of each reference genome
    """
    unique_reads = []

    for key in genome_dict:
        # if read mapp with best score to only one genome:
        if (len(genome_dict[key]) == 1):
            genome = ''.join(genome_dict[key]) # get genome id
            unique_dict[genome] += 1 # add to unique read dictionnary
            unique_reads.append(key)
    
    for key in unique_reads:
        del genome_dict[key] # delete the key from TMP dictionnary

    return genome_dict, unique_dict


def calculate_Co(genome_dict, unique_dict):
    """
    Calculate genome specific coefficient "Co" for each multiple read.

    Args:
        genome_dict [DICT] = Contains for each multiple read, all the reference
                             genomes name where he mapped
        unique_dict [DICT] = nb of unique read of each reference genome

    Returns:
        Co [DICT] = contains coefficient values for each couple (read, genome)
        read_dict [DICT] = contains all the multiple reads of a genome
                            key : genome
                            value : list a multiple reads
    """
    # initialise
    Co = {}
    read_dict = dict()

    for key in genome_dict:
        # for a multiple read, gets nb of unique reads from each genome
        s = [ unique_dict[genome] for genome in genome_dict[key] ]
        som = reduce(lambda x,y : x+y, s) # total number of unique reads
        if (som != 0):
            for genome in genome_dict[key]:
                nb_unique = unique_dict[genome] # get the nb of unique reads
                if(nb_unique != 0):
                    # calculate Co of multiple read for the given genome
                    Co[(key, genome)] = nb_unique/float(som)
                    read_dict.setdefault(genome, []).append(key) # append to
                                                                 # the dict

    # get all the multiple reads of a genome
    for genome, reads in read_dict.iteritems():
        read_dict[genome] = list(set(reads))

    return read_dict, Co


def calcul_AbM(unique_dict, read_dict, Co_dict, multiple_dict):
    """
    Calculates multiple reads abundance for each genome.
    Multiple reads abundance of a genome is equal to the sum of all Co
    coefficient

    Args:
        unique_dict [DICT] = nb of unique read of each reference genome
        read_dict [DICT] = list of multiple read for each reference genome
        Co_dict [DICT] = contains coefficient values for each couple
                         (read, genome)
        multiple_dict [DICT] = abundance of multiple reads for each genome

    Returns:
        multiple_dict [DICT] = abundance of multiple reads for each genome
    """
    for genome in read_dict:
        # get a list of all the Co coefficient for each reference genome
        s = [Co_dict[(read, genome)] for read in read_dict[genome]]
        som = sum(s) # sum them, (we obtain the abundance)
        multiple_dict[genome] = som

    return multiple_dict


def calcul_AbS(unique_dict, multiple_dict):
    """
    Calculates the abundance of a each genome.
        Abundance = Abundance unique reads + Abundance multiple reads

    Args:
        unique_dict [DICT] = nb of unique read of each reference genome
        multiple_dict [DICT] = abundance of multiple reads for each genome

    Returns:
        abundance_dict [DICT] = contains abundance of each reference genome
    """
    abundance_dict = {}

    for genome in unique_dict:
        abundance_dict[genome] = unique_dict[genome] + multiple_dict[genome]

    return abundance_dict


def write_stat(output, abundance_dict, database):
    """
    Write count table.
    
    Args:
        output [STRING] = output filename
        abundance_dict [DICT] = contains abundance of each reference genome
        database [DICT] = contrains length of each reference genome
    
    No Returns
    """
    with open(output, 'wt') as out:
        for genome, abundance in abundance_dict.items():
            out.write('{0}\t{1}\t{2}\n'.format(genome, database[genome], abundance))


#---------------------------- PROGRAMME MAIN ----------------------------------#

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description = " Reads counting Program ")
    
    parser.add_argument('samfile',
                        metavar = 'SAMFILE',
                        type = is_sam,
                        help = "input bam file to parse\n")
    
    parser.add_argument('output',
                        metavar = 'OUTPUT',
                        type = path_not_exists,
                        help = "Path to output the count table\n")
    
    parser.add_argument('-bam',
                        action = 'store_const',
                        const = 'bam',
                        help = "Write the BAM file\n")
                        
    parser.add_argument('-t','--threads',
                        type = int,
                        help = "number of threads to use to sort bam files\n")
                        
    parser.add_argument('-m', '--mem',
                        type = int,
                        help = "allocated memory for the job in Mb. Will be divided by the number of threads to optimise memory per threads in samtools sort\n")

    counting_method = parser.add_mutually_exclusive_group(required=True)
    
    counting_method.add_argument('--best',
                                 action = 'store_const',
                                 const = 'best',
                                 help = "Count each reads only once. "
                                         "If a read map to multiple location with"
                                         "a best score, one location is chosen "
                                         "randomly\n")
    
    counting_method.add_argument('--ex_aequo',
                                 action='store_const',
                                 const = 'ex_aequo',
                                 help = "Count all best hits. "
                                          "If a read have multiple best hits, "
                                          "all best hits are taken into "
                                          "account\n")

    counting_method.add_argument('--shared',
                                 action='store_const',
                                 const = 'shared',
                                 help = "Weights each reads according to the"
                                         "probability that the alignment does "
                                         "not correspond to the read's true "
                                         "point of origin\n")
    
    # command-line parsing
    if not len(sys.argv) > 2:
        parser.print_help()
        exit(1)

    try :
        args = parser.parse_args()
        args = vars(args)
        
    except:
        exit(1)
    
    
    # PARSING PARAMETERS -------------------------------------------------------
    
        # input & output
    samfile = args['samfile']
    output = args['output']
        
        # count mode
    
    # clean the dictionary
    for key, value in args.items():
            if value==None:
                del args[key]

    # First Step: Create BAM file
    if 'threads' in args:
        bamfile = create_bam(samfile, args['threads'])
    else:
        bamfile = create_bam(samfile)
    
    # Second Step : Count reads & create count table
    if 'best' in args:
        if 'threads' in args:
            # print "Multithreading done with", args['threads'], "threads"
            sort_file = sort_bam(bamfile, args['mem'], args['threads']) # sort bamfile
            mapped_file = only_mapped(sort_file, args['threads']) # mapped reads only
        else:
            # print "No Multithreading"
            sort_file = sort_bam(bamfile, args['mem'])
            mapped_file = only_mapped(sort_file)
        
        write_table(mapped_file, output) # write count table
    
    elif 'ex_aequo' in args:
        filtered_bam = read_bam(bamfile, "ex-aequo") # create a bamfile
        if 'threads' in args:
            sort_file = sort_bam(filtered_bam, args['mem'], args['threads'])
        else:
            sort_file = sort_bam(filtered_bam, args['mem'])
        write_table(sort_file, output)

    else : # shared
        # parsing du fichier bam
        if 'bam' in args:
            bam=True
        else:
            bam = False

        db, genomes = read_bam(bamfile, "shared", bam)

        # create dictionaries
        unique = dict.fromkeys(db, 0)
        multiple = dict.fromkeys(db, 0)
        
        # Filter unique reads from multiple reads
        genomes, unique = uniq_from_mult(genomes, unique)
        
        # For multiple reads
        reads, coef_read = calculate_Co(genomes, unique)   # calculate Co
        multiple = calcul_AbM(unique, reads, coef_read, multiple) # calculate
                                                                   # abundance
         
        # Calculate reference abundance & write count table
        abundance = calcul_AbS(unique, multiple)
        write_stat(output, abundance, db)
