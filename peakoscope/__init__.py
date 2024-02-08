# -*- coding: utf-8 -*-
# This file is part of Peakoscope.
# Copyright (C) 2021-2024  Eivind Tøstesen
# Peakoscope is licensed under GPLv3.
"""Peakoscope package.

Peakoscope is a python package.

"""

__version__ = "0.9.0"


# Import names:
from peakoscope.errors import PeakyBlunder
from peakoscope.utilities import ChainedAttributes
from peakoscope.trees import tree_from_peaks, Tree, HyperTree
from peakoscope.peaks import find_peaks, find_valleys, Region, Scope
from peakoscope.formats import TreeStrings
from peakoscope.data import example_1, example_2


def tree(data, valleys=False):
    """Find all peaks (or valleys) in data and return their Tree."""
    if valleys:
        return Tree.from_valleys(
            map(
                lambda t: Scope.from_start_end(*t[0:2], data),
                find_valleys(data),
            )
        )
    else:
        return Tree.from_peaks(
            map(
                lambda t: Scope.from_start_end(*t[0:2], data),
                find_peaks(data),
            )
        )
