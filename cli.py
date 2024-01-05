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
from peaks_and_valleys import find_peaks, find_valleys, NumSlice
from trees import Tree


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
    description="Print a tree of peak regions.",
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

# Extract data from input:
datalist = []
with args.inputfile as f:
    for line in f:
        datalist.append(float(line.strip()))

# Compute and print tree:
print(args.tree(datalist))
