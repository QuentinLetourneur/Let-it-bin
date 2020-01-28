#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

#std lib
import os
import argparse
# my lib
from fasta import Fasta

def bin():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_filename", help="Name of the fasta file from which to extract sequence")
    parser.add_argument("tsv_filename", help="Cluster file containing group ID and sequence ID")
    parser.add_argument("out", default = ".", help="Folder to place output. Will be created if do not exist")
    parser.add_argument("-p", "--prefix", default = "", help="Can place a prefix on the group names")
    args = parser.parse_args()

    # Parse the index file
    index = parse_tsv(args.tsv_filename)
    
    if not os.path.exists(args.out):
        os.makedirs(args.out)

    # create empty files
    create_files(index.values(), args.out, args.prefix)

    # read the fasta and fill bins
    fill_bins(args.fasta_filename, index, args.out, args.prefix)
    

def parse_tsv(tsv_filename):
    index = {}

    # Read the index file
    with open(tsv_filename) as file_reader:
        for line in file_reader:
            # Parse the line assuming group first, fasta header second
            group, header = line.strip().split("\t")

            # update the index
            if header in index:
                index[header].append(group)
            else:
                index[header] = [group]

    return index


def create_files(unsorted_groups, out, prefix):
    groups = set()

    # Create a set to filter duplicates
    for unsorted in unsorted_groups:
        groups.update(unsorted)

    # Create files
    for group in groups:
        open("{}/{}{}.fa".format(out, prefix, group), "w").close()


def fill_bins(fasta_filename, index, out, prefix):
    fasta = Fasta(fasta_filename)

    # Read sequence by sequence
    for sequence in fasta.read():
        if sequence["id"] in index:
            # get the sequence groups
            groups = index[sequence["id"]]

            # Write the sequence into the right bins
            for group in groups:
                with open("{}/{}{}.fa".format(out, prefix, group), "a") as fw:
                    fw.write(">{}\n{}\n".format(sequence["id"], sequence["value"]))


if __name__ == "__main__":
    bin()
