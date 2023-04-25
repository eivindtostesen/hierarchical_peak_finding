#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Python module for using pandas in peak analysis.

Tested with Pandas version 1.4.4

Created on Sat Jan  7 16:29:38 2023

@author: Eivind Tostesen

"""


import pandas as pd
from itertools import chain
from utilities import ChainedAttributes


# Alias when None means "use default":
_default = None


# Classes:


class PeakTreePandas(ChainedAttributes):
    """Pandas methods to be owned by a PeakTree."""

    def __init__(
        self,
        tree,
        attrname="pandas",
        location={},
        slices={},
        flanks={},
        **kwargs,
    ):
        """Attach this pandas-aware object to a PeakTree."""
        super().__init__()
        self.setattr(obj=tree, attrname=attrname)
        self.location = location
        self.slices = slices
        self.flanks = flanks
        if self.location:
            self.start = lambda n: self.location[n][0]
            self.end = lambda n: self.location[n][1]
        if self.slices:
            self.label_slice = lambda n: self.slices[n][0]
            self.value_slice = lambda n: self.slices[n][1]
        if self.flanks:
            self.left_flank = lambda n: self.flanks[n][0]
            self.right_flank = lambda n: self.flanks[n][1]
        self.node = lambda n: n
        self.root = lambda n: self.rootself.root()
        self.high = (
            lambda n: self.rootself.high(n) if self.rootself.has_children(n) else None
        )
        for name in (
            "parent children low full top is_nonroot "
            "has_children size height base_height _index"
        ).split():
            setattr(self, name, getattr(self.rootself, name))
        for name in (
            "root_path top_path subtree high_descendants "
            "low_descendants full_nodes leaf_nodes "
            "branch_nodes linear_nodes"
        ).split():
            setattr(
                self,
                name,
                lambda node, attr=name: list(getattr(self.rootself, attr)(node)),
            )
        self.set_definitions(**kwargs)

    def set_definitions(self, **kwargs):
        """Set attributes with functions, mappings or Series."""
        for name in kwargs:
            setattr(self, name, kwargs[name])

    def series(self, name="node", filter=_default, *, definitions={}):
        """Return series with one row per PeakTree node."""
        if filter is _default:
            filter = self.rootself
        return (
            pd.Series(filter)
            .map(
                definitions[name] if name in definitions else getattr(self, name),
                na_action="ignore",
            )
            .convert_dtypes()
            .rename(name)
        )

    def dataframe(self, columns="node", filter=_default, *, definitions={}):
        """Return dataframe with one row per PeakTree node."""
        if filter is _default:
            filter = self.rootself
        series = pd.Series(filter)
        return pd.concat(
            [
                self.series(name, series, definitions=definitions)
                for name in columns.split()
            ],
            axis=1,
        )

    def assign_columns(self, dataframe, columns="", *, definitions={}):
        """Assign extra columns to a given dataframe."""
        return dataframe.assign(
            **{
                name: self.series(name, dataframe["node"], definitions=definitions)
                for name in columns.split()
            }
        )

    def sort(self, dataframe, by="_index", **kwargs):
        """Return sorted dataframe."""
        return dataframe.sort_values(
            by="node",
            key=lambda col: col.map(getattr(self.rootself, by)),
            ignore_index=True,
            **kwargs,
        )

    # Out-of-the-box dataframes:

    def dump_data_attributes(self):
        """Return a dump of the PeakTree's data attributes."""
        df = pd.DataFrame(self.rootself.as_dict_of_dicts())
        df.index.name = "node"
        return df.reset_index().convert_dtypes().pipe(self.sort)

    def tree_structure(self):
        """Return dataframe with topological attributes."""
        return self.dataframe(
            "top children node parent full root",
            filter=chain.from_iterable(
                self.rootself.path(node, self.rootself.full(node), self.rootself.parent)
                for node in self.rootself.leaf_nodes()
            ),
        )

    def numeric_properties(self):
        """Return dataframe with numeric (vertical) properties."""
        if hasattr(self, "value_slice"):
            columns = "node height size base_height value_slice"
        else:
            columns = "node height size base_height"
        return (
            self.dataframe(columns)
            .pipe(self.sort, "size", ascending=False)
            .pipe(self.sort, "height", ascending=False)
        )

    def location_properties(self):
        """Return dataframe with locational (horizontal) properties."""
        df = self.dataframe(
            "node location start end",
            definitions=dict(
                location=lambda n: list(self.location[n]),
            ),
        )
        if hasattr(self, "label_slice"):
            df = df.pipe(
                self.assign_columns,
                columns="label_slice slice_length",
                definitions=dict(
                    slice_length=lambda n: len(self.slices[n][0]),
                ),
            )
        return df
