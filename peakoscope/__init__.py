# -*- coding: utf-8 -*-
# This file is part of Peakoscope.
# Copyright (C) 2021-2024  Eivind Tøstesen
# Peakoscope is licensed under GPLv3.
"""Peakoscope package.

Peakoscope is a python package for hierarchical analysis of peak and valley regions in numeric data.

Usage examples:
---------------

Compute the tree of nested peak regions in a data set:

>>> data = [10, 30, 40, 30, 10, 50, 70, 70, 50, 80]
>>> print(tree(data))
0:10
├─5:10
│ ├─9:10
│ └─6:8
└─1:4
  └─2:3

From the tree, select default peak regions and print their subarrays of data:

>>> for peak in tree(data).size_filter():
...    print(peak.subarray(data))
... 
[80]
[70, 70]
[30, 40, 30]

Copyright (C) 2021-2024 Eivind Tøstesen. This software is licensed under GPL-3.

"""


__version__ = "1.0.0"


# Import names:
from peakoscope.errors import PeakyBlunder
from peakoscope.utilities import ChainedAttributes
from peakoscope.trees import tree_from_peaks, Tree, HyperTree
from peakoscope.peaks import find_peaks, find_valleys, Scope6, Region, Scope
from peakoscope.formats import TreeStrings
from peakoscope.data import example_1, example_2


# Wrapper function:
def tree(data, *, valleys=False):
    """Find all peaks (or valleys) in data and return their Tree."""
    if valleys:
        return Tree.from_valleys(
            map(
                lambda t: Scope.from_attrs(Scope6(*t), data),
                find_valleys(data),
            )
        )
    else:
        return Tree.from_peaks(
            map(
                lambda t: Scope.from_attrs(Scope6(*t), data),
                find_peaks(data),
            )
        )
