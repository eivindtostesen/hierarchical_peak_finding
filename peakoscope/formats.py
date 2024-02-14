# -*- coding: utf-8 -*-
# This file is part of Peakoscope.
# Copyright (C) 2021-2024  Eivind Tøstesen
# Peakoscope is licensed under GPLv3.
"""Python module for string/text methods and formats.

"""


import pprint
from peakoscope.utilities import ChainedAttributes


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
                self.rootself.tip(full),
                self.rootself.main_child(full),
                self.rootself.parent,
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

    def boxdrawing(
        self,
        localroot=None,
        return_string=False,
        nonlast="├─",
        last="└─",
        vert="│ ",
        spaces="  ",
        margin="",
    ):
        """Print subtree using box drawing characters."""
        # defaults:
        if localroot is None:
            localroot = self.rootself.root()

        indent = [margin]
        lines = []
        for node, level in self.rootself.levels(localroot=localroot):
            if level == 0:
                # if node is localroot:
                lines.append(margin + str(node))
            elif node == self.rootself.children(self.rootself.parent(node))[-1]:
                # if node is a last child:
                del indent[level:]
                lines.append("".join([*indent, last, str(node)]))
                indent.append(spaces)
            else:
                # if node is a non-last child:
                del indent[level:]
                lines.append("".join([*indent, nonlast, str(node)]))
                indent.append(vert)
        return "\n".join(lines) if return_string else print(*lines, sep="\n")
