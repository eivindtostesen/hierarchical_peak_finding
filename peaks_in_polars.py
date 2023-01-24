#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Python module for using polars in peak analysis.

Created on Sat Jan  7 16:30:46 2023

@author: Eivind Tostesen

"""


import polars as pl
from utilities import ChainedAttributes


# Functions:


def siblings(tree, node, exclude_self=False):
    """Return tuple with siblings including input node."""
    if tree.is_nonroot(node):
        siblings = list(tree.children(tree.parent(node)))
        if exclude_self:
            siblings.remove(node)
        return tuple(siblings)
    else:
        return ()


def full_path(tree, node):
    """Yield nodes on path to the full node."""
    return tree.path(node, tree.full(node), tree.parent)


# Classes:


class PeakTreeMethods(ChainedAttributes):
    """Polars methods to be owned by a PeakTree."""

    def __init__(self, tree, attrname='polars',
                 location={}, slices={}, flanks={},
                 ):
        """Attach a polars-aware object to a PeakTree."""
        super().__init__()
        self.setattr(obj=tree, attrname=attrname)
        self.location = location
        self.slices = slices
        self.flanks = flanks

    def dump_data_attributes(self):
        """Return a dump of the PeakTree's data attributes."""

        def listing(data):
            if type(data) is dict:
                return [data[n] for n in self.rootself]
            else:
                return data

        return self.dataframe(
            generate_columns=["node"],
            **{n: listing(d)
                for n, d in self.rootself.as_dict_of_dicts().items()
               }
        )

    def dataframe(self,
                  nodes=None,
                  generate_columns=["node", "parent", "children", "high",
                                    "low", "root", "full", "mode",
                                    "is_nonroot", "has_children", "size",
                                    "height", "base_height", "node_index"],
                  **kwargs,
                  ):
        """Return a dataframe with properties of PeakTree nodes."""
        if nodes is None:
            nodes = self.rootself
        nodes = sorted(nodes, key=self.rootself.index)

        table = {}

        def listing(method, valid=None):
            if valid is None:
                return list(map(method, nodes))
            else:
                return [(method(n) if valid(n) else None) for n in nodes]

        if "node" in generate_columns:
            table["node"] = nodes
        if "parent" in generate_columns:
            table["parent"] = listing(self.rootself.parent)
        if "children" in generate_columns:
            table["children"] = listing(self.rootself.children)
        if "full" in generate_columns:
            table["full"] = listing(self.rootself.full)
        if "mode" in generate_columns:
            table["mode"] = listing(self.rootself.mode)
        if "root" in generate_columns:
            table["root"] = self.rootself.root()
        if "high" in generate_columns:
            table["high"] = listing(self.rootself.high,
                                    self.rootself.has_children
                                    )
        if "low" in generate_columns:
            table["low"] = listing(self.rootself.low)
        if "size" in generate_columns:
            table["size"] = listing(self.rootself.size)
        if "base_height" in generate_columns:
            table["base_height"] = listing(self.rootself.base_height)
        if "height" in generate_columns:
            table["height"] = listing(self.rootself.height)
        if "is_nonroot" in generate_columns:
            table["is_nonroot"] = listing(self.rootself.is_nonroot)
        if "has_children" in generate_columns:
            table["has_children"] = listing(self.rootself.has_children)
        if "node_index" in generate_columns:
            table["node_index"] = listing(self.rootself.index)

        return (
            pl.DataFrame({**table, **kwargs})
            .select(
                pl.col(generate_columns + list(kwargs.keys()))
            )
        )

    def tree_structure(self, nodes=None):
        """Return dataframe with hierarchical properties."""
        return self.dataframe(
            nodes,
            generate_columns=[
                "mode",
                "high",
                "children",
                "low",
                "node",
                "parent",
                "full",
                "root",
                "is_nonroot",
                "has_children",
                "node_index",
            ]
        )

    def peak_properties(self, nodes=None):
        """Return dataframe with peak properties."""
        if nodes is None:
            nodes = self.rootself
        nodes = sorted(nodes, key=self.rootself.index)
        return self.dataframe(
            nodes,
            generate_columns=["node", "size", "height", "base_height", "mode"],
            start=[self.location[n][0]
                   for n in nodes
                   ],
            end=[self.location[n][1]
                 for n in nodes
                 ],
            label_slice=[self.slices[n][0]
                         for n in nodes
                         ],
            value_slice=[self.slices[n][1]
                         for n in nodes
                         ],
            slice_length=[len(self.slices[n][0])
                          for n in nodes
                          ],
            left_flank=[self.flanks[n][0]
                        for n in nodes
                        ],
            right_flank=[self.flanks[n][1]
                         for n in nodes
                         ],
        )

    def graph_theoretical_properties(self, nodes=None):
        """Return dataframe with graph theoretical properties.

        Terminology: https://en.wikipedia.org/wiki/Tree_(data_structure)
        """
        if nodes is None:
            nodes = self.rootself
        nodes = sorted(nodes, key=self.rootself.index)
        return self.dataframe(
            nodes,
            generate_columns=["node"],
            siblings=[siblings(self.rootself, n)
                      for n in nodes
                      ],
            sibling_index=[(siblings(self.rootself, n).index(n)
                            if self.rootself.is_nonroot(n)
                            else None
                            )
                           for n in nodes
                           ],
            root_distance=[len(list(self.rootself.root_path(n))) - 1
                           for n in nodes
                           ],
            full_distance=[len(list(full_path(self.rootself, n))) - 1
                           for n in nodes
                           ],
            mode_distance=[len(list(self.rootself.mode_path(n))) - 1
                           for n in nodes
                           ],
            degree=[len(self.rootself.children(n))
                    for n in nodes
                    ],
            subtree_size=[len(list(self.rootself.subtree(n)))
                          for n in nodes
                          ],
        )
