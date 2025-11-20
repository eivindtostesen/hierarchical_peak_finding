# -*- coding: utf-8 -*-
# This file is part of Peakoscope.
# Copyright (C) 2021-2025  Eivind Tøstesen
# Peakoscope is licensed under GPLv3.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Peakoscope package.

Peakoscope is a python package for hierarchical analysis of peak and valley regions in numeric data.

Usage examples:
---------------

Compute the tree of nested peak regions in a data set:

>>> data = [10, 30, 40, 30, 10, 50, 70, 70, 50, 80]
>>> peaktree = tree(data)
>>> print(peaktree)
0:10
├─5:10
│ ├─9:10
│ └─6:8
└─1:4
  └─2:3

Print the tree minus its root:

>>> peaktree = tree(data, forest=True)
>>> print(peaktree - peaktree.roots())
5:10
├─9:10
└─6:8
1:4
└─2:3

Print the leaf nodes (which are local maxima in data) and corresponding slices of data:

>>> for leaf in peaktree.leaf_nodes():
...    print(leaf, leaf.subarray(data))
...
9:10 [80]
6:8 [70, 70]
2:3 [40]

Copyright (C) 2021-2025 Eivind Tøstesen. This software is licensed under GPL-3.0-or-later.

"""


__version__ = "1.2.0"


# Import names:
from peakoscope.errors import PeakyBlunder
from peakoscope.utilities import ChainedAttributes
from peakoscope.trees import (
    tree_from_peaks,
    forest_from_peaks,
    Tree,
    HyperTree,
    Forest,
    HyperForest,
)
from peakoscope.peaks import find_peaks, find_valleys, Scope6, Region, Scope
from peakoscope.formats import TreeStrings
from peakoscope.data import example_1, example_2


# Wrapper function:
def tree(data, *, valleys=False, forest=False):
    """Return a single tree of all peaks (or valleys) in data."""
    treeclass = Forest if forest else Tree
    _pipe1 = find_peaks(data, reverse=valleys)
    _pipe2 = map(lambda t: Scope.from_attrs(Scope6(*t), data), _pipe1)
    return treeclass(
        _pipe2,
        are_valleys=valleys,
        presorted=True,
    )
