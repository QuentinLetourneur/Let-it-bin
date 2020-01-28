#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Extract the annotation in case of blast on imomi database."""

from __future__ import print_function
import os
import sys
import argparse
import csv
#import requests
#from bs4 import BeautifulSoup


__author__ = "Amine Ghozlane"
__copyright__ = "Copyright 2014, INRA"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Amine Ghozlane"
__email__ = "amine.ghozlane@jouy.inra.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
       Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def isdir(path):
    """Check if path is an existing file.
       Parameters:
          path: Path to the file
    """
    if not os.path.isdir(path):
        if os.path.isfile(path):
            msg = "{0} is a file.".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


#===================
# parameters
#===================
def get_arguments():
    """
    Extract program options
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     usage="{0} -h [options] [arg]"
                                     .format(sys.argv[0]))
    parser.add_argument('-f', '--BlastResultFile', dest='blast_result_file',
                        type=isfile, required=True,
                        help='Input blast result file, in m8 mode.')
    parser.add_argument('-g', '--gi_taxid_taxonomy', dest='taxonomy_file',
                        type=isfile, required=True,
                        help='gi_taxid_taxonomy file for gi to taxid taxonomy '
                             'correspondancy')
    parser.add_argument('-nb', dest='nbest', type=int, default=0,
                        help='Number of best selected (default:0 - '
                        'based on the number of based aligned)')
    parser.add_argument('-fc', dest='filter_coverage', type=int, default=80,
                        help='Filter the coverage (default >= 80 percent).')
    parser.add_argument('-fi', dest='filter_identity', action='store_true',
                        default=False,
                        help='Remove filtering on the identity (default False, '
                        'superkingdom >= 0.0, phylum >= 65, class >= 75, '
                        'genus >= 85, specie >=95 percent ).')
    parser.add_argument('-id', dest='identity', type=str, default="",
                        help="Sample identity")
    parser.add_argument('-o', '--output_file', dest='output_file', type=str,
                        help='Output file')
    parser.add_argument('-r', dest='results', type=isdir,
                        default=os.curdir + os.sep, help='Path to result '
                        'directory.')
    return parser.parse_args()


#===========================================
# parse file
#===========================================

def parse_gi_to_taxid_taxonomy_file(taxonomy_file):
    """
    """
    gi_taxid_taxonomy_dict = {}
    try:
        with open(taxonomy_file, "rt") as taxonomy:
            taxonomy_reader = csv.reader(taxonomy, delimiter='\t')
            taxonomy_reader.next()
            for line in taxonomy_reader:
                #print(line)
                gi_taxid_taxonomy_dict[int(line[0])] = line[2]
    except IOError:
        sys.exit("Error cannot open {0}".format(taxonomy_file))
    return gi_taxid_taxonomy_dict


#===========================================
# Extract blast annotation
#===========================================

def extract_annotation(blast_result_file, gi_taxid_taxonomy_dict):
    """
    Extract blast annotation
    """
    blast_dict = {}
    try:
        # Remove query dict
        with open(blast_result_file, "rt") as blast_result:
            blast_reader = csv.reader(blast_result, delimiter="\t")
            for line in blast_reader:
                annotation_list = line[1].split('|')
                gi = int(annotation_list[1])
                #print(gi)
                if gi in gi_taxid_taxonomy_dict:
                    annotation = str(gi) + ";" + gi_taxid_taxonomy_dict[gi]
                else:
                    annotation = None
                # id identity coverage
                if line[0] in blast_dict:
                    blast_dict[line[0]] += [[annotation,
                                             float(line[10]), float(line[11]),
                                             float(line[12])]]
                else:
                    blast_dict[line[0]] = [[annotation,
                                             float(line[10]), float(line[11]), 
                                             float(line[12])]]
            assert(len(blast_dict) > 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(blast_result_file))
    except AssertionError:
        sys.exit("Error nothing read from {0}".format(blast_result_file))
    return blast_dict

#filter_identity, 
def write_annotation(blast_dict, nbest, filter_coverage, filter_identity,
                     output_file, results, identity):
    """
    Write the result
    """
    split_thres = {0:0, 65:2, 75:3, 85:6, 95:7}
    if not output_file:
        output_file = results + os.sep + 'ncbi_taxonomic_annotation.txt'
    try:
        with open(output_file, "wt") as output:
            output_writer = csv.writer(output, delimiter="\t")
            #output_writer.writerow(["ContigName", "GI", "superkingdom", "kingdom",
            #                        "phylum", "class", "order", "family",
            #                        "genus","species", "PourcID",
            #                        "Coverage", "evalue"])
            for key in blast_dict:
                #if filter_coverage > 0:
                data = []
                for element in blast_dict[key]:
                    #print(element)
                    if element[2] >= float(filter_coverage):
                        data.append(element)
                #else:
                    #data = blast_data[key]
                # Sort depending on the identity + coverage
                data.sort(key=lambda x: x[1] + x[2], reverse=True)
                if nbest > 0:
                    short_set = data[0: nbest]
                else:
                    short_set = data
                #print(short_set)
                
                matched_genome = ""
                
                for hit in short_set:
                    #print(hit)
                    if hit[0]:
                        if hit[1] >= 95.0 or filter_identity:
                            hit[0] = hit[0].split(";")
                        elif hit[1] >= 85.0:
                            hit[0] = hit[0].split(";")[0:8] + ["NA"]
                        elif hit[1] >= 75.0:
                            hit[0] = hit[0].split(";")[0:5] + ["NA"] * 4
                        elif hit[1] >= 65.0:
                            hit[0] = hit[0].split(";")[0:4] + ["NA"] * 5
                        elif hit[1] > 0.0:
                            hit[0] = hit[0].split(";")[0:3] + ["NA"] * 6
                    else:
                        hit[0] = ["NA"] * 9
                        
                    # get the genome name corresponding to the GI
                    #~ req = requests.get("https://www.ncbi.nlm.nih.gov/nuccore/" + hit[0][0])
                    #~ html_text = BeautifulSoup(req.text, "html.parser")
                    #~ matched_genome = str(html_text.find_all("h1")[1].text)
                    #~ matched_genome = matched_genome.split(",")[0]
                    #~ matched_genome = matched_genome.replace(".","")
                    #~ matched_genome = matched_genome.replace(" ","_")
                    #~ print(hit[0].split(";")[0])
                
                length = key.split("_")[3]
                
                if identity != "":
                    for element in short_set:
                        output_writer.writerow([identity + "_" + key] + [length] + element[0] + element[1:])
                else:
                    for element in short_set:
                        output_writer.writerow([key] + [length] + element[0] + element[1:])
    except IOError:
        sys.exit("Error cannot open {0}".format(output_file))


#===================
# MAIN
#===================
def main():
    """
    Main program
    """
    args = get_arguments()
    ## Parse gi to taxid
    #print('Parse gi to taxid file')
    #taxid_gi_dict = parse_gi_to_taxid_file(args.gi_taxid_file)
    #print(taxid_gi_dict)
    ## Parse taxid to taxonomy
    print('Parse taxid to taxonomy file')
    #gi_taxonomy_dict = parse_taxid_to_taxonomy_file(args.taxid_taxonomy_file,
                                                     #taxid_gi_dict)
    gi_taxid_taxonomy_dict = parse_gi_to_taxid_taxonomy_file(args.taxonomy_file)
    #print(gi_taxid_taxonomy_dict)
    # Parse blast result
    print('Parse blast results')
    blast_dict = extract_annotation(args.blast_result_file,
                                    gi_taxid_taxonomy_dict)
    #print(blast_dict)
    # Write annotation
    print('Write annotation')
    #args.filter_identity,
    write_annotation(blast_dict, args.nbest, args.filter_coverage,
                     args.filter_identity, args.output_file, args.results,
                     args.identity)


if __name__ == "__main__":
    main()
# END
