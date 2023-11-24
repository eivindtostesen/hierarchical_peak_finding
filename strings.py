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

    def repr(self, return_string=False):
        """Prettyprint a repr that can reconstruct the tree."""
        if hasattr(self.rootself, "L"):
            string = f"HyperTree(\n{self.rootself.L.string.repr(return_string=True)}, \n{self.rootself.R.string.repr(return_string=True)}\n)"
        else:
            string = f"Tree.from_levels(\n{pprint.pformat(dict(self.rootself.levels()), indent=2, sort_dicts=False)}\n)"
        return string if return_string else print(string)

    def indented_list(self, localroot=None, indent="| ", return_string=False):
        """Print tree in indented list notation."""
        # defaults:
        if localroot is None:
            localroot = self.rootself.root()
        if return_string:
            return "\n".join(
                level * indent + str(node)
                for node, level in self.rootself.levels(localroot=localroot)
            )
        else:
            for node, level in self.rootself.levels(localroot=localroot):
                print(level * indent + str(node))

    def riverflow(self, localroot=None, return_string=False):
        """Print subtree in riverflow notation."""
        # defaults:
        if localroot is None:
            localroot = self.rootself.root()
        string = "# Notation: <main> /& <lateral>/ => <parent>\n"
        for full in self.rootself.full_nodes(localroot):
            if not self.rootself.has_children(full):
                continue
            for node in self.rootself.path(
                self.rootself.core(full), self.rootself.main(full), self.rootself.parent
            ):
                string += f"{node}"
                if len(self.rootself.children(self.rootself.parent(node))) == 2:
                    string += (
                        f" /& {self.rootself.lateral(self.rootself.parent(node))[0]}/"
                    )
                if len(self.rootself.children(self.rootself.parent(node))) > 2:
                    string += (
                        " /& "
                        + ", ".join(self.rootself.lateral(self.rootself.parent(node)))
                        + "/"
                    )
                string += " => "
            string += f"{full}.\n"
        return string if return_string else print(string)
