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


"""Get NCBI taxonomy with lineage."""

from __future__ import print_function
import argparse
import os
import sys
#import fileinput
try:
    from cogent.parse.ncbi_taxonomy import NcbiTaxonomyFromFiles
except ImportError:
    sys.exit("The program requires the cogent package")


__author__ = "Amine Ghozlane"
__copyright__ = "Copyright 2014, INRA"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Amine Ghozlane"
__email__ = "amine.ghozlane@jouy.inra.fr"
__status__ = "Developpement"


def isfile(path):
    """
    Check if path is an existing file.
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


def get_arguments():
    """
    Retrieves the arguments of the program.
        Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h (see also "
                                     "ftp://ftp.ncbi.nih.gov/pub/taxonomy/)"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='blast_output_file', type=isfile,
                        required=True, help="blast_output_file")
    parser.add_argument('-t', dest='ncbi_taxid_file', type=isfile,
                        required=True, help="ncbi_taxid_file")
    parser.add_argument('-n', dest='ncbi_names_file', type=isfile,
                        required=True, help="ncbi_names_file")
    parser.add_argument('-d', dest='ncbi_nodes_file', type=isfile,
                        required=True, help="ncbi_nodes_file")
    parser.add_argument('-l', dest="lineage", action="store",
                        default=['superkingdom','kingdom','phylum','class',
                                 'order','family','genus','species'],
                        type=list, help="""List of ranks to reported in
                        taxonomy lineage (default ['superkingdom','kingdom',
                        'phylum','class','order','family','genus','species'])""")
    parser.add_argument('-o', dest='taxonomy_file', type=str, required=True,
                        help="Output taxonomy_file")
    return parser.parse_args()


def extract_genbank_id(blast_output_file):
    """
    Extracting Genbank IDS from BLAST output
    """
    len_gis = 0
    gis = set()
    try:
        with open(blast_output_file, "rt") as blast_output :
            for aline in blast_output :
                if aline[0] == '#' or aline.strip() == '' :
                    pass
                elif "|" in aline:
                    gis.add(aline.split()[1].split('|')[1])
                    gis.add(aline.split('|')[1])
            assert(len(gis) > 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(blast_output_file))
    except AssertionError:
        sys.exit("No gi read from {0}".format(blast_output_file))
    except IndexError:
        sys.exit("{0}".format(aline))
    return gis


def load_taxid_nucl(ncbi_taxid_file, gis):
    """
    Parsing NCBI file 'gi_taxid_nucl.dmp'
    """
    dict_gi_taxid = {}
    try:
        with open(ncbi_taxid_file, "rt") as ncbi_taxid :
            for aline in ncbi_taxid :
                ids = aline.split()
                if ids[0] in gis :
                    dict_gi_taxid[ int(ids[0]) ] = int(ids[1])
            assert(len(dict_gi_taxid.keys()) > 0)
    except IOError:
        sys.exit("Error cannot open {0}".format(ncbi_taxid_file))
    except AssertionError:
        sys.exit("No gi read from {0}".format(ncbi_taxid_file))
    return dict_gi_taxid


def load_ncbi_tree(ncbi_nodes_file, ncbi_names_file):
    """
    Loading NCBI taxnomy tree
    """
    try:
        ncbi_tree = NcbiTaxonomyFromFiles(open(ncbi_nodes_file, "rt"),
                                          open(ncbi_names_file, "rt"))
    except IOError:
        sys.exit("Error cannot open {0} or {1}".format(ncbi_nodes_file,
                                                       ncbi_names_file))
    return ncbi_tree


def assign_taxa_lineage(dict_gi_taxid, ncbi_tree, lineage):
    """
    Generating taxa lineages from NCBI taxnomy tree
    """
    lineage_len = len(lineage)
    dict_taxid_lineage = {}
    dict_taxid_tmp = {}
    unknown = []
    for taxid in dict_gi_taxid.values() :
        # if taxid has already been processed, do nothing and move to next iteration
        names = ['NA'] * lineage_len
        if taxid in ncbi_tree.ById:
            curr = ncbi_tree.ById[taxid]
            newtaxid = taxid
            while curr.ParentId != 1 :
                # if already stored, get rank and name for current taxid
                if newtaxid in dict_taxid_tmp.keys() :
                    rank, name = dict_taxid_tmp[newtaxid]
                # else get info from tree and store
                else :
                    try :
                        rank, name = curr.Rank, curr.Name
                        dict_taxid_tmp[newtaxid] = (rank, name)
                    except :
                        # report unfound taxid
                        unknown.append(newtaxid)
                        pass
                # update taxon lineage if current rank is requested
                if rank in lineage :
                    names[lineage.index(rank)] = name
                # process parent node in next iteration
                newtaxid = curr.ParentId
                curr = ncbi_tree.ById[newtaxid]
                # store taxon lineage
        dict_taxid_lineage[taxid] = ";".join(names)
    return dict_taxid_lineage, unknown


def write_results(ranks, dict_gi_taxid, dict_taxid_lineage, taxonomy_file):
    """
    Writing results to file
    """
    try:
        with open(taxonomy_file, "wt") as output:
            output.write("gi\ttaxid\t{0}\n".format(";".join(ranks)))
            for gi, taxid in dict_gi_taxid.items():
                if dict_taxid_lineage.has_key(taxid):
                    output.write("{0}\t{1}\t{2}\n".format(
                                    gi, taxid, dict_taxid_lineage[taxid]))
    except IOError:
        sys.exit("Error cannot open {0}".format(taxonomy_file))


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    #ranks = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family',
             #'genus', 'species']
    # Get arguments
    args = get_arguments()
    # Step 1
    print("STEP 1: Extracting Genbank IDS from BLAST output...")
    gis = extract_genbank_id(args.blast_output_file)
    print("Found {0} ids !".format(len(gis)))
    # Step 2
    print("STEP 2: Parsing NCBI file 'gi_taxid_nucl.dmp'...")
    dict_gi_taxid = load_taxid_nucl(args.ncbi_taxid_file, gis)
    print("Stored {0} ids !".format(len(dict_gi_taxid.keys())))
    # Step 3
    print("STEP 3: Loading NCBI taxonomy tree...")
    #, ncbi_root
    try:
        ncbi_tree = load_ncbi_tree(args.ncbi_nodes_file, args.ncbi_names_file)
        print("Tree Loaded !")
        # Step 4
        print("STEP 4: Generating taxa lineages from NCBI taxonomy tree...")
        dict_taxid_lineage, unknown = assign_taxa_lineage(dict_gi_taxid,
                                                          ncbi_tree,
                                                          args.lineage)
        if unknown != []:
            print("WARNING: taxids not found in ncbi : {0}".format(unknown))
        print("Writing results to file '{0}'...".format(args.taxonomy_file))
        write_results(args.lineage, dict_gi_taxid, dict_taxid_lineage, args.taxonomy_file)
        print("DONE !")
    except ValueError as what:
        print(what, file=sys.stderr)


if __name__ == '__main__':
    main()
