#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Command Line Interface.

Examples:

Print peaks:  
  $ python cli.py test.csv
  $ cat test.csv | python cli.py -
  $ cat test.csv | python cli.py
  $ ./cli.py test.csv
  $ cat test.csv | ./cli.py -
  $ cat test.csv | ./cli.py

Print valleys:
  $ python cli.py --valleys test.csv
  $ python cli.py -v test.csv
  $ cat test.csv | python cli.py --valleys
  $ cat test.csv | python cli.py -v

Display help:
  $ python cli.py --help


Created on Fri Dec 29 18:36:00 2023.

@author: Eivind Tostesen

"""


import sys
import argparse
import csv
from peaks_and_valleys import find_peaks, find_valleys, NumSlice
from trees import Tree


def _data_from_csv(args):
    """Extract data from input."""
    data = []
    args.inputfile.reconfigure(newline="")
    with args.inputfile as f:
        reader = csv.reader(f, delimiter=args.delim, quoting=csv.QUOTE_NONE)
        for row in reader:
            data.append(float(row[args.field - 1]))
    return data


def _tree_of_peaks(data):
    """Find peaks in data and make tree."""
    return Tree.from_peaks(
        map(lambda t: NumSlice.from_start_end(*t[0:2], data), find_peaks(data))
    )


def _tree_of_valleys(data):
    """Find valleys in data and make tree."""
    return Tree.from_valleys(
        map(lambda t: NumSlice.from_start_end(*t[0:2], data), find_valleys(data))
    )


# Define CLI with arguments and options:
parser = argparse.ArgumentParser(
    description="Print a tree of peak regions. Regions are written in slice notation.",
    epilog="see also: https://github.com/eivindtostesen/hierarchical_peak_finding",
)
parser.add_argument(
    "inputfile",
    nargs="?",
    type=argparse.FileType("r"),
    default=sys.stdin,
    help="a file or stdin to read data from",
)
parser.add_argument(
    "-d",
    "--delimiter",
    default=" ",
    dest="delim",
    help="the delimiter used in inputfile",
)
parser.add_argument(
    "-f",
    "--field",
    type=int,
    default=1,
    help="column number (count from 1) with data in inputfile",
)
parser.add_argument(
    "-v",
    "--valleys",
    action="store_const",
    const=_tree_of_valleys,
    default=_tree_of_peaks,
    dest="tree",
    help="print valley regions instead of peak regions",
)


# Parse command-line arguments:
args = parser.parse_args()
# Compute and print tree:
print(args.tree(_data_from_csv(args)))
