#!/usr/bin/env python3.
# -*- coding: utf-8 -*-
"""Python module for string methods.

Created on Tue Oct 24 21:30:00 2023

@author: Eivind Tostesen

"""


from utilities import ChainedAttributes


# Classes:


class TreeStrings(ChainedAttributes):
    """String methods to be owned by a PeakTree or HyperPeakTree."""

    def __init__(
        self,
        tree,
        attrname="string",
    ):
        """Attach a string methods object."""
        super().__init__()
        self.setattr(obj=tree, attrname=attrname)
        self.level = self.rootself.levels()

    def indented_list(self, indent = "| ", return_string=False):
        """Print subtree in indented list notation."""
        string = "\n".join([level * indent + str(node) for node, level in self.level.items()])
        if return_string:
            return string
        else:
            print(string)

    def riverflow(self, localroot=None, return_string=False):
        """Print subtree in riverflow notation."""
        # defaults:
        if localroot is None:
            localroot = self.rootself.root()
        string = "# Notation: <high> /& <low>/ => <parent>\n"
        for full in self.rootself.full_nodes(localroot):
            if not self.rootself.has_children(full):
                continue
            for node in self.rootself.path(self.rootself.top(full), self.rootself.high(full), self.rootself.parent):
                string += f"{node}"
                if len(self.rootself.children(self.rootself.parent(node))) == 2:
                    string += f" /& {self.rootself.low(self.rootself.parent(node))[0]}/"
                if len(self.rootself.children(self.rootself.parent(node))) > 2:
                    string += f" /& {self.rootself.low(self.rootself.parent(node))}/"
                string += " => "
            string += f"{full}.\n"
        if return_string:
            return string
        else:
            print(string)

