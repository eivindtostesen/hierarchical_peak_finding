# -*- coding: utf-8 -*-
# This file is part of Peakoscope.
# Copyright (C) 2021-2025  Eivind Tøstesen
# Peakoscope is licensed under GPLv3.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Python module for using polars in peak/valley analysis.

Class TreePolars provides methods called dataframe, tree_structure,
numeric_properties, location_properties and node_properties
that return polars dataframes with one row per Tree node.

"""


from itertools import chain
from operator import attrgetter, methodcaller
import polars as pl
from peakoscope.utilities import ChainedAttributes


# Alias when None means "use default":
_default = None


# Classes:


class TreePolars(ChainedAttributes):
    """Polars methods to be owned by a Tree."""

    def __init__(
        self,
        tree,
        attrname="df",
        X=None,
        node_attributes=(
            "start stop istop "
            "max min argmax argmin size "
            "argext argcut extremum cutoff "
        ).split(),
        node_methods=(
            "pre post "
            "is_peak is_valley is_local_maximum is_local_minimum "
            "__len__ __repr__ __str__ "
        ).split(),
        node_methods_to_lists=("__iter__ subarray ").split(),
        tree_methods=(
            "main_child parent full tip is_nonroot has_children _index"
        ).split(),
        tree_methods_to_lists=(
            "children lateral "
            "root_path main_path size_filter "
            "subtree main_descendants lateral_descendants full_nodes "
            "leaf_nodes branch_nodes linear_nodes"
        ).split(),
        objecttype=("node root parent full tip main_child").split(),
        **kwargs,
    ):
        """Attach this polars-aware object to a Tree."""
        super().__init__()
        self.setattr(obj=tree, attrname=attrname)
        if X:
            self.x_start = lambda n: X[n.start]
            self.x_end = lambda n: X[n.istop]
        self.node = lambda n: n
        self.root = lambda _: self.rootself.root()
        for name in node_attributes:
            setattr(self, name, attrgetter(name))
        for name in node_methods:
            setattr(self, name, methodcaller(name))
        for name in node_methods_to_lists:
            setattr(self, name, lambda node, met=name: list(methodcaller(met)(node)))
        for name in tree_methods:
            setattr(self, name, getattr(self.rootself, name))
        for name in tree_methods_to_lists:
            setattr(
                self,
                name,
                lambda node, attr=name: list(getattr(self.rootself, attr)(node)),
            )
        self.objecttype = objecttype
        self.set_definitions(**kwargs)

    def set_definitions(self, **kwargs):
        """Set attributes with functions."""
        for name in kwargs:
            setattr(self, name, kwargs[name])

    def series(self, name="node", filter=_default, *, definitions={}):
        """Return series with one row per Tree node."""
        if filter is _default:
            filter = self.rootself
        dtype = pl.Object if name in self.objecttype else None
        return pl.Series(
            name,
            map(
                definitions[name] if name in definitions else getattr(self, name),
                filter,
            ),
            dtype=dtype,
        )

    def dataframe(self, columns="node", filter=_default, *, definitions={}):
        """Return dataframe with one row per Tree node."""
        if filter is _default:
            filter = iter(self.rootself)
        nodes = list(iter(filter))
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

    def sort_by_max_and_size(self, dataframe):
        """Return dataframe sorted by max and size."""
        return dataframe.pipe(self.sort, by="size", descending=True).pipe(
            self.sort, by="max", descending=True
        )

    # Out-of-the-box dataframes:

    def dump_data_attributes(self):
        """Return a dump of the Tree's data attributes."""

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

    def tree_structure(self):
        """Return dataframe with topological attributes."""
        return self.dataframe(
            "tip children node parent full root",
            filter=chain.from_iterable(
                self.rootself.path(node, self.rootself.full(node), self.rootself.parent)
                for node in self.rootself.leaf_nodes()
            ),
        )

    def numeric_properties(self):
        """Return dataframe with numeric (vertical) properties."""
        return self.dataframe("node max size min").pipe(self.sort_by_max_and_size)

    def location_properties(self):
        """Return dataframe with locational (horizontal) properties."""
        return self.dataframe(
            "node location x_start x_end",
            definitions=dict(
                location=lambda n: f"{self.x_start(n)}..{self.x_end(n)}",
            ),
        )

    def node_properties(self):
        """Return dataframe with node attributes and methods."""
        return self.dataframe(
            (
                "__str__ __repr__ __len__ "
                "start stop istop pre post "
                "argmax argmin max min size "
                "argext argcut extremum cutoff "
                "is_peak is_valley is_local_maximum is_local_minimum "
            )
        )
