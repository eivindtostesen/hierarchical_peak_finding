#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Python module for using polars in peak analysis.

Tested with Polars version 0.15.14

Created on Sat Jan  7 16:30:46 2023

@author: Eivind Tostesen

"""


import polars as pl
from utilities import ChainedAttributes


# Alias when None means "use default":
_default = None


# Functions:


def _siblings(tree, node, exclude_self=False):
    """Return tuple with siblings including input node."""
    if tree.is_nonroot(node):
        siblings = list(tree.children(tree.parent(node)))
        if exclude_self:
            siblings.remove(node)
        return tuple(siblings)
    else:
        return ()


def _full_path(tree, node):
    """Yield nodes on path to the full node."""
    return tree.path(node, tree.full(node), tree.parent)


# Classes:


class PeakTreePolars(ChainedAttributes):
    """Polars methods to be owned by a PeakTree."""

    def __init__(
        self,
        tree,
        attrname="polars",
        location={},
        slices={},
        flanks={},
        **kwargs,
    ):
        """Attach this polars-aware object to a PeakTree."""
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
            "parent full top is_nonroot has_children size height base_height _index"
        ).split():
            setattr(self, name, getattr(self.rootself, name))
        for name in (
            "children low root_path top_path subtree high_descendants "
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
        """Set attributes with functions."""
        for name in kwargs:
            setattr(self, name, kwargs[name])

    def series(self, name="node", filter=_default, *, definitions={}):
        """Return series with one row per PeakTree node."""
        if filter is _default:
            filter = self.rootself
        return pl.Series(
            name,
            map(
                definitions[name] if name in definitions else getattr(self, name),
                filter,
            ),
        )

    def dataframe(self, columns="node", filter=_default, *, definitions={}):
        """Return dataframe with one row per PeakTree node."""
        if filter is _default:
            filter = iter(self.rootself)
        nodes = pl.Series(filter)
        return pl.DataFrame().hstack(
            [
                self.series(name, nodes, definitions=definitions)
                for name in columns.split()
            ],
            in_place=True,
        )

    def assign_columns(self, dataframe, columns="", *, definitions={}):
        """Assign extra columns to a given dataframe."""
        return pl.concat(
            [
                dataframe,
                self.dataframe(
                    columns,
                    filter=dataframe.get_column("node"),
                    definitions=definitions,
                ),
            ],
            how="horizontal",
        )

    def sort(self, dataframe, by="_index", **kwargs):
        """Return sorted dataframe."""
        return (
            dataframe.pipe(
                self.assign_columns,
                columns="sorting_column",
                definitions=dict(sorting_column=getattr(self, by)),
            )
            .sort("sorting_column", **kwargs)
            .select([pl.exclude("sorting_column")])
        )

    # Out-of-the-box dataframes:

    def dump_data_attributes(self):
        """Return a dump of the PeakTree's data attributes."""

        def _listorscalar(data):
            if type(data) is dict:
                return [data[n] for n in self.rootself]
            else:
                return data

        return pl.concat(
            [
                self.dataframe(),
                pl.DataFrame(
                    {
                        n: _listorscalar(d)
                        for n, d in self.rootself.as_dict_of_dicts().items()
                    }
                ),
            ],
            how="horizontal",
        ).pipe(self.sort)
