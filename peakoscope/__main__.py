# -*- coding: utf-8 -*-
"""Command Line Interface.

The CLI is run by running the package (using the -m flag).

Examples:

Display version:
  $ python -m peakoscope --version

Print peaks:  
  $ python -m peakoscope test.csv
  $ cat test.csv | python -m peakoscope -
  $ cat test.csv | python -m peakoscope

Print valleys:
  $ python -m peakoscope --valleys test.csv
  $ python -m peakoscope -v test.csv
  $ cat test.csv | python -m peakoscope --valleys
  $ cat test.csv | python -m peakoscope -v

Display help:
  $ python -m peakoscope -h
  $ python -m peakoscope --help


Created on Fri Dec 29 18:36:00 2023.

@author: Eivind Tostesen

"""


import sys
import argparse
import csv
import peakoscope


def _data_from_csv(args):
    """Return list with data from csv input."""
    data = []
    args.inputfile.reconfigure(newline="")
    with args.inputfile as f:
        reader = csv.reader(f, delimiter=args.delim, quoting=csv.QUOTE_NONE)
        for row in reader:
            data.append(float(row[args.field - 1]))
    return data


# Define CLI with arguments and options:
parser = argparse.ArgumentParser(
    prog="python -m peakoscope",
    description="Print a tree of peak regions. Regions are written in slice notation.",
    epilog="see also: https://github.com/eivindtostesen/hierarchical_peak_finding",
)
parser.add_argument(
    "--version", action="version", version="Peakoscope " + peakoscope.__version__
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
    action="store_true",
    help="print valley regions instead of peak regions",
)

# Parse command-line arguments:
args = parser.parse_args()
# Compute and print tree:
print(peakoscope.tree(_data_from_csv(args), valleys=args.valleys))
