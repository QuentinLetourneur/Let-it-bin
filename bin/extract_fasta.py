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

"""Find extract homologous proteins in a multifasta."""

from __future__ import print_function
import argparse
import os
import sys
import re
import textwrap


def isfile(path):
    """
    Check if path is an existing file.
    Arguments:
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
    """
    Check if path is an existing file.
    Arguments:
        path: Path to the file
    """
    if not os.path.isdir(path):
        if os.path.isfile(path):
            msg = "{0} is a file.".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def getArguments():
    """
    Retrieves the arguments of the program.
    Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h".format(sys.argv[0]))
    parser.add_argument('-q', dest='query_file', type=isfile, required=True,
                        help='Path to the query file (blast).')
    parser.add_argument('-t', dest='target_file', type=isfile, required=True,
                        help='Path to the database file.')
    parser.add_argument('-n', dest='not_in_database', action='store_true',
                        help='Select instead elements which are not in the'
                        ' list.')
    parser.add_argument('-v', dest='verbose', action='store_true',
                        help='Print operation details.')
    parser.add_argument('-o', dest='output_file', type=str,
                        help='Output file.')
    parser.add_argument('-r', dest='results', type=isdir,
                        default=os.curdir + os.sep,
                        help='Path to result directory.')
    args = parser.parse_args()
    return args


def get_query(query_file):
    """
    """
    regex_query = re.compile("^(\S+)\s")
    query_list = []
    try:
        with open(query_file, "rt") as query:
            for line in query:
                if line[0] != "#":
                    match_query = regex_query.match(line)
                    if match_query:
                        query_list += [match_query.group(1)]
    except IOError:
        sys.exit("Error cannot open {0}".format(query_file))
    return query_list


def check_reference(title, list_query):
    """
    """
    for q in list_query:
        if q in title:
            return q
    return None


def get_sequence(list_query, target_file, not_in_database):
    """
    """
    result = {}
    known = None
    try:
        with open(target_file, "rt") as target:
            for line in target:
                if known and line[0] != ">":
                    result[known] += line[0:].strip().replace("\n", "").replace("\r", "")
                if line[0] == ">":
                    if len(list_query) == 0:
                        break
                    known = None
                    title = line[1:].replace("\n", "").replace("\r", "")
                    if " " in title:
                        title = title.split(" ")[0]
                    known = check_reference(title, list_query)
                    if known and not not_in_database:
                        list_query.pop(list_query.index(known))
                        result[known] = ""
                    elif not known and not_in_database:
                        known = title
                        result[known] = ""
                    elif known and not_in_database:
                        known = None
            if not not_in_database:
                assert(len(list_query) == 0)
    except IOError:
        sys.exit("Error : cannot open {0}".format(target_file))
    except AssertionError:
        print("The program have not find every query sequence, "
            "check the following elements :", file=sys.stderr)
        print(list_query, file=sys.stderr)
    return result


def fill(text, width=80):
    """Split text"""
    return os.linesep.join(text[i:i+width] for i in xrange(0, len(text), width))


def write_sequence(results, sequence_data, output_file, verbose):
    """
    """
    if not output_file:
        output_file = results + os.sep + "extracted_sequence.fasta"
    try:
        with open(output_file, "wt") as output:
            for seq in sequence_data:
                output.write(">{1}{0}{2}{0}".format(os.linesep, seq,
                                                    fill(sequence_data[seq])))
                if verbose:
                    print("write : {0}, length : {1}".format(seq,
                        len(sequence_data[seq])))
    except IOError:
        sys.exit("Error : cannot open {0}".format(output_file))


def main():
    """
    """
    # Load the arguments
    args = getArguments()
    # Load query elements
    print("Load query list")
    list_query = get_query(args.query_file)
    # Grab query sequence in the database
    print("Load database sequence")
    sequence_data = get_sequence(list_query, args.target_file,
                                 args.not_in_database)
    # Write the new fasta file
    print("Write the new fasta")
    write_sequence(args.results, sequence_data, args.output_file, args.verbose)
    print("Done")


if __name__ == '__main__':
    main()
