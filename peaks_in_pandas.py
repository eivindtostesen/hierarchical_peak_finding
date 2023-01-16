#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Python module for using pandas in peak analysis.

Created on Sat Jan  7 16:29:38 2023

@author: Eivind Tostesen

"""


import pandas as pd
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
    """Pandas methods to be owned by a PeakTree."""

    def __init__(self, tree, attrname='pandas',
                 location={}, slices={}, flanks={},
                 ):
        """Attach a pandas-aware object to a PeakTree."""
        super().__init__()
        self.setattr(obj=tree, attrname=attrname)
        self.location = location
        self.slices = slices
        self.flanks = flanks

    def dump_data_attributes(self):
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
        table = {}

        def mapping(method):
            return {node: method(node) for node in nodes}

        if "node" in generate_columns:
            table["node"] = {n: n for n in nodes}
        if "parent" in generate_columns:
            table["parent"] = mapping(self.rootself.parent)
        if "children" in generate_columns:
            table["children"] = mapping(self.rootself.children)
        if "full" in generate_columns:
            table["full"] = mapping(self.rootself.full)
        if "mode" in generate_columns:
            table["mode"] = mapping(self.rootself.mode)
        if "root" in generate_columns:
            table["root"] = self.rootself.root()
        if "high" in generate_columns:
            table["high"] = {n: self.rootself.high(n)
                             for n in nodes
                             if self.rootself.has_children(n)
                             }
        if "low" in generate_columns:
            table["low"] = mapping(self.rootself.low)
        if "size" in generate_columns:
            table["size"] = mapping(self.rootself.size)
        if "base_height" in generate_columns:
            table["base_height"] = mapping(self.rootself.base_height)
        if "height" in generate_columns:
            table["height"] = mapping(self.rootself.height)
        if "is_nonroot" in generate_columns:
            table["is_nonroot"] = mapping(self.rootself.is_nonroot)
        if "has_children" in generate_columns:
            table["has_children"] = mapping(self.rootself.has_children)
        if "node_index" in generate_columns:
            table["node_index"] = mapping(self.rootself.index)

        return (
            pd.DataFrame({**table, **kwargs},
                         columns=generate_columns + list(kwargs.keys()),
                         )
            .convert_dtypes()
            .sort_index(key=lambda col: col.map(self.rootself.index))
            .reset_index(drop=True)
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
        return self.dataframe(
            nodes,
            generate_columns=["node", "size", "height", "base_height", "mode"],
            start={n: self.location[n][0]
                   for n in nodes
                   },
            end={n: self.location[n][1]
                 for n in nodes
                 },
            label_slice={n: self.slices[n][0]
                         for n in nodes
                         },
            value_slice={n: self.slices[n][1]
                         for n in nodes
                         },
            slice_length={n: len(self.slices[n][0])
                          for n in nodes
                          },
            left_flank={n: self.flanks[n][0]
                        for n in nodes
                        },
            right_flank={n: self.flanks[n][1]
                         for n in nodes
                         },
            )

    def graph_theoretical_properties(self, nodes=None):
        """Return dataframe with graph theoretical properties.

        Terminology: https://en.wikipedia.org/wiki/Tree_(data_structure)
        """
        if nodes is None:
            nodes = self.rootself
        return self.dataframe(
            nodes,
            generate_columns=["node"],
            siblings={n: siblings(self.rootself, n)
                      for n in nodes
                      },
            sibling_index={n: siblings(self.rootself, n).index(n)
                           for n in nodes
                           if self.rootself.is_nonroot(n)
                           },
            root_distance={n: len(list(self.rootself.root_path(n))) - 1
                           for n in nodes
                           },
            full_distance={n: len(list(full_path(self.rootself, n))) - 1
                           for n in nodes
                           },
            mode_distance={n: len(list(self.rootself.mode_path(n))) - 1
                           for n in nodes
                           },
            degree={n: len(self.rootself.children(n))
                    for n in nodes
                    },
            subtree_size={n: len(list(self.rootself.subtree(n)))
                          for n in nodes
                          },
            )
