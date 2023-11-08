#!/usr/bin/env python3.
# -*- coding: utf-8 -*-
"""Python module for string methods.

Created on Tue Oct 24 21:30:00 2023

@author: Eivind Tostesen

"""


import pprint
from utilities import ChainedAttributes


# Classes:


class TreeStrings(ChainedAttributes):
    """String methods to be owned by a Tree or HyperTree."""

    def __init__(
        self,
        tree,
        attrname="string",
    ):
        """Attach a string methods object."""
        super().__init__()
        self.setattr(obj=tree, attrname=attrname)
        self.level = self.rootself.levels()

    def repr(self, return_string=False):
        """Prettyprint a repr that can reconstruct the tree."""
        if hasattr(self.rootself, "L"):
            string = f"HyperTree(\n{self.rootself.L.string.repr(return_string=True)}, \n{self.rootself.R.string.repr(return_string=True)}\n)"
        else:
            string = f"Tree.from_levels(\n{pprint.pformat(self.level, indent=2, sort_dicts=False)}\n)"
        return string if return_string else print(string)

    def indented_list(self, indent="| ", return_string=False):
        """Print tree in indented list notation."""
        string = "\n".join(
            [level * indent + str(node) for node, level in self.level.items()]
        )
        return string if return_string else print(string)

    def riverflow(self, localroot=None, return_string=False):
        """Print subtree in riverflow notation."""
        # defaults:
        if localroot is None:
            localroot = self.rootself.root()
        string = "# Notation: <high> /& <low>/ => <parent>\n"
        for full in self.rootself.full_nodes(localroot):
            if not self.rootself.has_children(full):
                continue
            for node in self.rootself.path(
                self.rootself.top(full), self.rootself.high(full), self.rootself.parent
            ):
                string += f"{node}"
                if len(self.rootself.children(self.rootself.parent(node))) == 2:
                    string += f" /& {self.rootself.low(self.rootself.parent(node))[0]}/"
                if len(self.rootself.children(self.rootself.parent(node))) > 2:
                    string += f" /& {self.rootself.low(self.rootself.parent(node))}/"
                string += " => "
            string += f"{full}.\n"
        return string if return_string else print(string)
