# -*- coding: utf-8 -*-
"""Peakoscope package.

Peakoscope is a python package.

Created on Wed Jan 10 23:18:00 2024.

@author: Eivind Tostesen

"""


# Import names:
from peakoscope.errors import PeakyBlunder
from peakoscope.utilities import ChainedAttributes
from peakoscope.trees import tree_from_peaks, Tree, HyperTree
from peakoscope.peaks_and_valleys import find_peaks, find_valleys, SliceStr, NumSlice
from peakoscope.strings import TreeStrings
from peakoscope.peaky_data import example_1, example_2
