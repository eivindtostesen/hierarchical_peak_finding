#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Python module for using pandas in peak analysis.

Created on Sat Jan  7 16:29:38 2023

@author: Eivind Tostesen

"""


import pandas as pd
from utilities import ChainedAttributes


# Classes:


class PeakTreeMethods(ChainedAttributes):
    """Pandas methods to be owned by a PeakTree."""

    def __init__(self, tree, attrname='pandas'):
        """Attach a pandas-aware object to a PeakTree."""
        super().__init__()
        self.setattr(obj=tree, attrname=attrname)

    def dataframe(self,
                  generate_columns=["node", "parent", "children", "high",
                                    "low", "root", "full", "mode",
                                    "is_nonroot", "has_children", "size",
                                    "height", "base_height"],
                  **kwargs,
                  ):
        """Return a dataframe with data related to the PeakTree."""
        table = {}
        if "node" in generate_columns:
            table["node"] = {n: n for n in self.rootself}
        if "parent" in generate_columns:
            table["parent"] = self.rootself._parent
        if "children" in generate_columns:
            table["children"] = self.rootself._children
        if "full" in generate_columns:
            table["full"] = self.rootself._full
        if "mode" in generate_columns:
            table["mode"] = self.rootself._mode
        if "root" in generate_columns:
            table["root"] = self.rootself._root
        if "high" in generate_columns:
            table["high"] = {n: self.rootself.high(n)
                             for n in self.rootself
                             if self.rootself.has_children(n)}
        if "low" in generate_columns:
            table["low"] = {n: self.rootself.low(n)
                            for n in self.rootself
                            if self.rootself.has_children(n)}
        if "size" in generate_columns:
            table["size"] = {n: self.rootself.size(n)
                             for n in self.rootself}
        if "base_height" in generate_columns:
            table["base_height"] = {n: self.rootself.base_height(n)
                                    for n in self.rootself}
        if "height" in generate_columns:
            table["height"] = {n: self.rootself.height(n)
                               for n in self.rootself}
        if "is_nonroot" in generate_columns:
            table["is_nonroot"] = {n: self.rootself.is_nonroot(n)
                                   for n in self.rootself}
        if "has_children" in generate_columns:
            table["has_children"] = {n: self.rootself.has_children(n)
                                     for n in self.rootself}
        if "tree_index" in generate_columns:
            table["tree_index"] = {n: self.rootself.index(n)
                                   for n in self.rootself}
        return (
            pd.DataFrame({**table, **kwargs},
                         columns=generate_columns + list(kwargs.keys()),
                         )
            .convert_dtypes()
            .sort_index(key=lambda col: col.map(self.rootself.index))
            .reset_index(drop=True)
            )

    def data_attributes(self):
        """Return a dump of the PeakTree's data attributes."""
        df = pd.DataFrame(self.rootself.as_dict_of_dicts())
        df.index.name = "node"
        return (
            df
            .reset_index()
            .convert_dtypes()
            .sort_values(by="node",
                         key=lambda col: col.map(self.rootself.index),
                         )
            )
