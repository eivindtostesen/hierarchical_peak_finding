#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Python module for using pandas in peak analysis.

Created on Sat Jan  7 16:29:38 2023

@author: Eivind Tostesen

"""


from itertools import chain
import pandas as pd
from utilities import ChainedAttributes


# Alias when None means "use default":
_default = None


# Classes:


class TreePandas(ChainedAttributes):
    """Pandas methods to be owned by a Tree."""

    def __init__(
        self,
        tree,
        attrname="pandas",
        X=None,
        objecttype=(
            "node root parent full top children high low"
            "root_path top_path subtree high_descendants "
            "low_descendants full_nodes leaf_nodes "
            "branch_nodes linear_nodes"
        ).split(),
        **kwargs,
    ):
        """Attach this pandas-aware object to a Tree."""
        super().__init__()
        self.setattr(obj=tree, attrname=attrname)
        if X is None:
            self.x_start = lambda n: n.start
            self.x_end = lambda n: n.end
        else:
            self.x_start = lambda n: X[n.start]
            self.x_end = lambda n: X[n.end]
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
        self.objecttype = objecttype
        self.set_definitions(**kwargs)

    def set_definitions(self, **kwargs):
        """Set attributes with functions, mappings or Series."""
        for name in kwargs:
            setattr(self, name, kwargs[name])

    def series(self, name="node", filter=_default, *, definitions={}):
        """Return series with one row per Tree node."""
        if filter is _default:
            filter = self.rootself
        dtype = "object" if name in self.objecttype else None
        s = (
            pd.Series(filter, dtype=dtype)
            .map(
                definitions[name] if name in definitions else getattr(self, name),
                na_action="ignore",
            )
            .rename(name)
        )
        if name not in self.objecttype:
            s.convert_dtypes()
        return s

    def dataframe(self, columns="node", filter=_default, *, definitions={}):
        """Return dataframe with one row per Tree node."""
        if filter is _default:
            filter = self.rootself
        nodes = list(iter(filter))
        return pd.concat(
            [
                self.series(name, nodes, definitions=definitions)
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

    def sort_by_height_and_size(self, dataframe, **kwargs):
        """Return dataframe sorted by descending height and size."""
        return dataframe.pipe(self.sort, "size", ascending=False, **kwargs).pipe(
            self.sort, "height", ascending=False, kind="stable", **kwargs
        )

    # Out-of-the-box dataframes:

    def dump_data_attributes(self):
        """Return a dump of the Tree's data attributes."""
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
        return self.dataframe("node height size base_height").pipe(
            self.sort_by_height_and_size
        )

    def location_properties(self):
        """Return dataframe with locational (horizontal) properties."""
        return self.dataframe(
            "node location x_start x_end",
            definitions=dict(
                location=lambda n: f"{self.x_start(n)}..{self.x_end(n)}",
            ),
        )
