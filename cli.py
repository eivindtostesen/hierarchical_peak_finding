#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Command Line Interface.

Examples:

$ python cli.py test.csv
$ cat test.csv | python cli.py -
$ cat test.csv | python cli.py
$ ./cli.py test.csv
$ cat test.csv | ./cli.py -
$ cat test.csv | ./cli.py

Created on Fri Dec 29 18:36:00 2023.

@author: Eivind Tostesen

"""


import fileinput
from peaks_and_valleys import find_peaks, NumSlice
from trees import Tree


data = []
with fileinput.FileInput() as f:
    for line in f:
        data.append(float(line.strip()))

tree = Tree.from_peaks(
    map(lambda peak: NumSlice.from_start_end(*peak[0:2], data), find_peaks(data))
)
print(tree)
