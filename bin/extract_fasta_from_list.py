#!/usr/bin/env python2.7

import os
import sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("fasta_file")
parser.add_argument("header_list")
parser.add_argument("output")

args=parser.parse_args()

fasta_ind = SeqIO.index(args.fasta_file, "fasta")

with open(args.header_list, 'r') as hl:
    for line in hl.readline():
        f
