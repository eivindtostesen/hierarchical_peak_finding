# -*- coding: utf-8 -*-
# This file is part of Peakoscope.
# Copyright (C) 2021-2025  Eivind Tøstesen
# Peakoscope is licensed under GPLv3.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Python module for trees of peak regions or valley regions.

This module contains algorithms and classes
for building trees that represent the hierarchical nesting
of regions and subregions containing peaks or valleys
in numeric one-dimensional or higher-dimensional data.

The tree classes provide methods for
searching, sorting and selecting regions.

"""


from operator import attrgetter
from peakoscope.utilities import pairwise


# Functions:


def tree_from_peaks(
    peaks,
    *,
    presorted=True,
    reverse=False,
    getstart=attrgetter("start"),
    getistop=attrgetter("istop"),
    getcutoff=attrgetter("cutoff"),
    getextremum=attrgetter("extremum"),
):
    """Return (parent, root, children, tip) from peak regions having start, istop, cutoff, extremum."""
    parent = {}
    children = {}
    tip = {}
    in_spe = []
    if not presorted:
        peaks = list(peaks)
        # Order same as given by 'find_peaks' function:
        peaks.sort(key=getcutoff, reverse=not reverse)
        peaks.sort(key=getistop)
    for p in peaks:
        children[p] = []
        while in_spe and getstart(p) <= getstart(in_spe[-1]):
            c = in_spe.pop()
            children[p].append(c)
            parent[c] = p
        children[p].sort(key=getstart)
        children[p].sort(key=getextremum, reverse=not reverse)
        children[p] = tuple(children[p])
        tip[p] = tip[children[p][0]] if children[p] else p
        in_spe.append(p)
    root = in_spe.pop()
    parent[root] = None
    return parent, root, children, tip


def forest_from_peaks(
    peaks,
    *,
    presorted=False,
    reverse=False,
    getstart=attrgetter("start"),
    getistop=attrgetter("istop"),
    getcutoff=attrgetter("cutoff"),
    getextremum=attrgetter("extremum"),
    getargext=attrgetter("argext"),
):
    """Return (roots, children, tip) from peak regions having start, istop, cutoff, extremum, argext."""

    def sorted_tuple(alist):
        alist.sort(key=getstart)
        alist.sort(key=getextremum, reverse=not reverse)
        return tuple(alist)

    children = {}
    tip = {}
    in_spe = []
    if not presorted:
        peaks = list(peaks)
        # Order same as given by 'find_peaks' function:
        peaks.sort(key=getcutoff, reverse=not reverse)
        peaks.sort(key=getistop)
    for p in peaks:
        children[p] = []
        while in_spe and getstart(p) <= getstart(in_spe[-1]):
            c = in_spe.pop()
            children[p].append(c)
        children[p] = sorted_tuple(children[p])
        if children[p] and getargext(children[p][0]) == getargext(p):
            tip[p] = tip[children[p][0]]
        else:
            tip[p] = p
        in_spe.append(p)
    roots = sorted_tuple(in_spe)
    return roots, children, tip


# Classes:


class Tree:
    """Tree of regions in univariate data.

    A Tree represents the hierarchical nesting of
    peak or valley regions in 1D data such as a sequence of numbers,
    a time series, a function y(x) or other univariate data.

    A Tree is initialized with an iterable of regions
    that have start, istop, cutoff, extremum. The regions
    must be unique hashable objects to be used as dictionary keys.

    Notes
    -----
    Background literature for the Tree class is
    the subsection titled "1D peaks" in the article [1]_.

    References
    ----------
    .. [1] Tostesen, E. "A stitch in time: Efficient computation of
       genomic DNA melting bubbles." Algorithms for Molecular
       Biology 3, 10 (2008).
       Open access: https://doi.org/10.1186/1748-7188-3-10
    """

    # Class variables and class methods:

    getstart = attrgetter("start")
    getistop = attrgetter("istop")
    getcutoff = attrgetter("cutoff")
    getextremum = attrgetter("extremum")

    @classmethod
    def from_peaks(cls, peaks, **kwargs):
        """Return new Tree from iterable of peak regions."""
        obj = cls.__new__(cls)
        obj._parent, obj._root, obj._children, obj._tip = tree_from_peaks(
            peaks, **kwargs
        )
        obj._find_full()
        return obj

    @classmethod
    def from_valleys(cls, valleys, **kwargs):
        """Return new Tree from iterable of valley regions."""
        obj = cls.__new__(cls)
        obj._parent, obj._root, obj._children, obj._tip = tree_from_peaks(
            valleys, reverse=True, **kwargs
        )
        obj._find_full()
        return obj

    @classmethod
    def from_levels(cls, levelsdict, /):
        """Return new Tree from other tree's levels-dict."""

        def leaf_and_tip(node):
            children[node] = []
            for n in obj.path(node, obj._full[node], obj.parent):
                obj._tip[n] = node

        obj = cls.__new__(cls)
        obj._parent = {}
        children = {}
        obj._tip = {}
        obj._full = {}
        for (A, a), (B, b) in pairwise(levelsdict.items()):
            if a == 0:  # first item is the root:
                obj._parent[A] = None
                obj._full[A] = A
                obj._root = A
                stack = [A]
            if b == a + 1:  # a subtree grows:
                obj._full[B] = obj._full[A]
                children[stack[-1]] = [B]  # B is the main child (of A)
            else:  # (then a >= b) a subtree finishes:
                del stack[b:]  # pop a slice
                obj._full[B] = B
                children[stack[-1]].append(B)  # B is a lateral child
                leaf_and_tip(A)  # A is a leaf and tip
            obj._parent[B] = stack[-1]
            stack.append(B)
        leaf_and_tip(B)  # the last B is a leaf and tip
        obj._children = {p: tuple(c) for p, c in children.items()}
        return obj

    def __init__(self, peaks, *, are_valleys=False, presorted=False):
        """Initialize Tree from iterable of peak (or valley) regions."""
        self._parent, self._root, self._children, self._tip = tree_from_peaks(
            peaks,
            presorted=presorted,
            reverse=are_valleys,
            getstart=Tree.getstart,
            getistop=Tree.getistop,
            getcutoff=Tree.getcutoff,
            getextremum=Tree.getextremum,
        )
        self._find_full()

    def __contains__(self, node):
        """Return True if the input is a node in the Tree."""
        return node in self._full

    def __iter__(self):
        """Iterate over nodes in the Tree."""
        return iter(self._full)

    def __len__(self):
        """Return number of nodes in the Tree."""
        return len(self._full)

    def __matmul__(self, other):
        """Return product self @ other (other is a Tree or HyperTree)."""
        return HyperTree(self, other)

    def __repr__(self) -> str:
        """Return string that can reconstruct the Tree."""
        return f"Tree.from_levels({repr(dict(self.levels()))})"

    def __str__(self):
        """Return tree as string using box drawing characters."""
        indent = [""]
        lines = []
        for node, level in self.levels():
            if level == 0:
                # if node is root:
                lines.append(str(node))
            elif node == self.children(self.parent(node))[-1]:
                # if node is a last child:
                del indent[level:]
                lines.append("".join([*indent, "└─", str(node)]))
                indent.append("  ")
            else:
                # if node is a non-last child:
                del indent[level:]
                lines.append("".join([*indent, "├─", str(node)]))
                indent.append("│ ")
        return "\n".join(lines)

    def as_dict_of_dicts(self):
        """Return data attributes as a dict of dicts."""
        return {
            "_parent": self._parent,
            "_children": self._children,
            "_tip": self._tip,
            "_full": self._full,
            "_root": self._root,
        }

    def set_nodes(self, changes={}):
        """Replace Tree nodes by using given mapping."""

        def new(node):
            return changes[node] if node in changes else node

        self._tip = dict((new(x), new(y)) for (x, y) in self._tip.items())
        self._full = dict((new(x), new(y)) for (x, y) in self._full.items())
        self._parent = dict((new(x), new(y)) for (x, y) in self._parent.items())
        self._children = dict(
            (new(x), tuple(new(z) for z in y)) for (x, y) in self._children.items()
        )
        self._root = new(self._root)
        return None

    def levels(self, localroot=None, level=0):
        """Yield ordered sequence of (node, level) tuples (root is zero level)."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        yield (localroot, level)
        for child in self.children(localroot):
            yield from self.levels(child, level + 1)

    def root(self):
        """Return the root node of the Tree."""
        return self._root

    def roots(self):
        """Return tuple with the root node of the Tree."""
        return (self.root(),)

    def is_nonroot(self, node):
        """Return True if the given node has a parent."""
        return node != self._root

    def tip(self, node):
        """Return the smallest region with same argext as the given node."""
        return self._tip[node]

    def has_children(self, node):
        """Return True if the given node has sub regions."""
        return node != self._tip[node]

    def size(self, node):
        """Return vertical size of given node (max minus min)."""
        return abs(Tree.getextremum(node) - Tree.getcutoff(node))

    def parent(self, node):
        """Return the parent (containing region) or None."""
        return self._parent[node]

    def children(self, node):
        """Return ordered tuple of children of given node."""
        return self._children[node]

    def main_child(self, node):
        """Return child of given node with same argext, or None."""
        if self.has_children(node):
            return self.children(node)[0]
        else:
            return None

    def lateral(self, node):
        """Return ordered tuple of the given node's lateral children."""
        return self.children(node)[1:]

    def full(self, node):
        """Return the largest region with same argext as the given node."""
        return self._full[node]

    def _index(self, node):
        """Return the zero-based index of the given node."""
        # accessing the preserved insertion order:
        return list(self._full).index(node)

    # public recursive algorithms:

    def path(self, start, istop, step):
        """Yield nodes on a path in the tree."""
        climber = start
        yield climber
        while climber != istop:
            climber = step(climber)
            yield climber

    def root_path(self, node):
        """Yield nodes on the parent path from given node to its root."""
        yield from self.path(node, self.root(), self.parent)

    def main_path(self, node):
        """Yield nodes on the main_child path from given node to its tip."""
        yield from self.path(node, self.tip(node), self.main_child)

    def subtree(self, localroot=None):
        """Yield all nodes in the given subtree."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        yield localroot
        for child in self.children(localroot):
            yield from self.subtree(child)

    def main_descendants(self, localroot=None):
        """Yield main child nodes in the given subtree."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        if self.has_children(localroot):
            yield self.main_child(localroot)
            for child in self.children(localroot):
                yield from self.main_descendants(child)

    def lateral_descendants(self, localroot=None):
        """Yield lateral child nodes in the given subtree."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        if self.has_children(localroot):
            yield from self.lateral(localroot)
            for child in self.children(localroot):
                yield from self.lateral_descendants(child)

    def full_nodes(self, localroot=None):
        """Yield full nodes in the given subtree."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        if localroot == self.full(localroot):
            yield localroot
        yield from self.lateral_descendants(localroot)

    def leaf_nodes(self, localroot=None):
        """Yield subtree nodes that have no children."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        return (
            node for node in self.subtree(localroot) if len(self.children(node)) == 0
        )

    def branch_nodes(self, localroot=None):
        """Yield subtree nodes that have two or more children."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        return (
            node for node in self.subtree(localroot) if len(self.children(node)) > 1
        )

    def linear_nodes(self, localroot=None):
        """Yield subtree nodes that have one child."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        return (
            node for node in self.subtree(localroot) if len(self.children(node)) == 1
        )

    def size_filter(self, localroot=None, *, maxsize=None):
        """Yield subtree nodes filtered by size."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        if maxsize is None:
            maxsize = 0.2 * self.size(self.root())
        # The 'MAXDEEP algorithm' in reverse:
        for climber in self.main_path(localroot):
            if self.size(climber) >= maxsize:
                for child in self.lateral(climber):
                    yield from self.size_filter(maxsize=maxsize, localroot=child)
            elif climber == self.root() or self.size(self.parent(climber)) >= maxsize:
                yield climber
                break

    def innermost(self, nodes, localroot=None):
        """Yield innermost nodes of the given nodes."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        filter = list(nodes)
        countdown = {n: len(self.children(n)) for n in self.subtree(localroot)}

        # Recursive "bottom-up" search via parents:
        def _yield_or_propagate(node):
            if node in filter:
                yield node
            elif node != localroot:
                parent = self.parent(node)
                countdown[parent] -= 1
                if countdown[parent] == 0:
                    yield from _yield_or_propagate(parent)

        for node in self.leaf_nodes(localroot):
            yield from _yield_or_propagate(node)

    def outermost(self, nodes, localroot=None):
        """Yield outermost nodes of the given nodes."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        filter = list(nodes)

        # Recursive "top-down" search via children:
        def _yield_or_branch(node):
            if node in filter:
                yield node
            elif self.has_children(node):
                for child in self.children(node):
                    yield from _yield_or_branch(child)

        yield from _yield_or_branch(localroot)

    # Initialization algorithms:

    def _find_full(self):
        """Compute attribute: self._full."""

        def fullnodes():
            yield self.root()
            yield from self.lateral_descendants()

        self._full = {
            node: full for full in fullnodes() for node in self.main_path(full)
        }


class Forest:
    """Trees of regions in univariate data.

    A Forest represents the hierarchical nesting of
    peak or valley regions in 1D data such as a sequence of numbers,
    a time series, a function y(x) or other univariate data.

    A Forest is initialized with an iterable of regions
    that have start, istop, cutoff, extremum, argext. The regions
    must be unique hashable objects to be used as dictionary keys
    and they become the tree nodes.

    A Forest consists of zero, one or more trees with a root each.
    But a Forest is a container of tree nodes, not trees.

    Notes
    -----
    Background literature for the Forest class is
    the subsection titled "1D peaks" in the article [1]_.

    References
    ----------
    .. [1] Tostesen, E. "A stitch in time: Efficient computation of
       genomic DNA melting bubbles." Algorithms for Molecular
       Biology 3, 10 (2008).
       Open access: https://doi.org/10.1186/1748-7188-3-10
    """

    # Class variables and class methods:

    getstart = attrgetter("start")
    getistop = attrgetter("istop")
    getcutoff = attrgetter("cutoff")
    getextremum = attrgetter("extremum")
    getargext = attrgetter("argext")

    @classmethod
    def from_peaks(cls, peaks, **kwargs):
        """Return new Forest from iterable of peak regions."""
        obj = cls.__new__(cls)
        obj.are_valleys = False
        obj._roots, obj._children, obj._tip = forest_from_peaks(peaks, **kwargs)
        obj._find_full_parent()
        return obj

    @classmethod
    def from_valleys(cls, valleys, **kwargs):
        """Return new Forest from iterable of valley regions."""
        obj = cls.__new__(cls)
        obj.are_valleys = True
        obj._roots, obj._children, obj._tip = forest_from_peaks(
            valleys, reverse=True, **kwargs
        )
        obj._find_full_parent()
        return obj

    @classmethod
    def from_levels(cls, levelsdict, /, *, are_valleys=False):
        """Return new Forest from other forest's levels-dict."""

        def _make_tip(node):
            # full path shares this tip
            for n in obj.path(node, obj._full[node], obj.parent):
                obj._tip[n] = node

        def _make_leaf(node):
            # has no children and is tip
            children[node] = []
            _make_tip(node)

        def _make_root(node):
            # has no parent and is full
            obj._parent[node] = None
            obj._full[node] = node
            obj._roots.append(node)

        obj = cls.__new__(cls)
        obj.are_valleys = are_valleys
        obj._parent = {}
        children = {}
        obj._tip = {}
        obj._full = {}
        obj._roots = []
        stack = []
        for (A, a), (B, b) in pairwise(levelsdict.items()):
            if not stack:
                _make_root(A)  # the first element is a root
                stack = [A]
            if b == 0:  # B is a root
                _make_leaf(A)
                _make_root(B)
                stack = [B]
            elif b == a + 1:
                children[A] = [B]  # B is the first child of A
                obj._parent[B] = A
                if Forest.getargext(A) == Forest.getargext(
                    stack[-1]
                ):  # main path continues
                    obj._full[B] = obj._full[A]
                else:  # main path ends at A
                    _make_tip(A)
                    obj._full[B] = B
                stack.append(B)
            else:  # then a >= b > 0:
                _make_leaf(A)  # A is a leaf and tip
                obj._full[B] = B
                del stack[b:]  # remove finished nodes
                children[stack[-1]].append(B)  # B is a lateral child
                obj._parent[B] = stack[-1]
                stack.append(B)
        _make_leaf(B)  # the last element is a leaf and tip
        obj._children = {p: tuple(c) for p, c in children.items()}
        obj._roots = tuple(obj._roots)
        return obj

    # dunder methods:

    def __init__(self, peaks, *, are_valleys=False, presorted=False):
        """Initialize Forest from iterable of peak (or valley) regions."""
        self.are_valleys = are_valleys
        self._roots, self._children, self._tip = forest_from_peaks(
            peaks,
            presorted=presorted,
            reverse=are_valleys,
            getstart=Forest.getstart,
            getistop=Forest.getistop,
            getcutoff=Forest.getcutoff,
            getextremum=Forest.getextremum,
            getargext=Forest.getargext,
        )
        self._find_full_parent()

    def __contains__(self, node):
        """Return True if the input is a node in the Forest."""
        return node in self._full

    def __iter__(self):
        """Iterate over nodes in the Forest."""
        return iter(self._full)

    def __len__(self):
        """Return number of nodes in the Forest."""
        return len(self._full)

    def __repr__(self) -> str:
        """Return string that can reconstruct the Forest."""
        return f"Forest.from_levels({repr(dict(self.levels()))}, are_valleys={self.are_valleys})"

    def __str__(self):
        """Return Forest as string using box drawing characters."""
        indent = [""]
        lines = []
        for node, level in self.levels():
            if level == 0:
                # if node is root:
                lines.append(str(node))
            elif node == self.children(self.parent(node))[-1]:
                # if node is a last child:
                del indent[level:]
                lines.append("".join([*indent, "└─", str(node)]))
                indent.append("  ")
            else:
                # if node is a non-last child:
                del indent[level:]
                lines.append("".join([*indent, "├─", str(node)]))
                indent.append("│ ")
        return "\n".join(lines)

    def __matmul__(self, other):
        """Return product self @ other (other is a Forest or HyperForest)."""
        return HyperForest(self, other)

    def __sub__(self, other):
        """Return new Forest of nodes in self not other (difference)."""
        return Forest(
            set(self) - set(other), are_valleys=self.are_valleys, presorted=False
        )

    def __and__(self, other):
        """Return new Forest of nodes in self and other (intersection)."""
        return Forest(
            set(self) & set(other), are_valleys=self.are_valleys, presorted=False
        )

    def __or__(self, other):
        """Return new Forest of nodes in self or other (union)."""
        return Forest(
            set(self) | set(other), are_valleys=self.are_valleys, presorted=False
        )

    def __xor__(self, other):
        """Return new Forest of nodes in either self or other (symmetric difference)."""
        return Forest(
            set(self) ^ set(other), are_valleys=self.are_valleys, presorted=False
        )

    # whole forest methods:

    def as_dict_of_dicts(self):
        """Return data attributes as a dict of dicts."""
        rootsdict = {node: root for root in self._roots for node in self.subtree(root)}
        return {
            "_parent": self._parent,
            "_children": self._children,
            "_tip": self._tip,
            "_full": self._full,
            "_roots": rootsdict,
        }

    def set_nodes(self, changes={}):
        """Replace Forest nodes by using given mapping."""

        def _new(node):
            return changes[node] if node in changes else node

        self._tip = dict((_new(x), _new(y)) for (x, y) in self._tip.items())
        self._full = dict((_new(x), _new(y)) for (x, y) in self._full.items())
        self._parent = dict((_new(x), _new(y)) for (x, y) in self._parent.items())
        self._children = dict(
            (_new(x), tuple(_new(z) for z in y)) for (x, y) in self._children.items()
        )
        self._roots = tuple(_new(root) for root in self._roots)
        return None

    # node methods:

    def root(self, node=None):
        """Return root of single tree, root of given node, or None."""
        if len(self.roots()) == 1:
            return self.roots()[0]
        elif node is not None and len(self.roots()) > 1:
            climber = node
            while self.is_nonroot(climber):
                climber = self.parent(climber)
                climber = self.full(climber)
            return climber
        else:
            return None

    def roots(self):
        """Return tuple of roots of trees in forest."""
        return self._roots

    def is_nonroot(self, node):
        """Return True if the given node has a parent."""
        return node not in self._roots

    def tip(self, node):
        """Return the smallest region with same argext as the given node."""
        return self._tip[node]

    def has_children(self, node):
        """Return True if the given node has sub regions."""
        return bool(self._children[node])

    def size(self, node):
        """Return vertical size of given node (max minus min)."""
        return abs(Forest.getextremum(node) - Forest.getcutoff(node))

    def parent(self, node):
        """Return the parent (containing region) or None."""
        return self._parent[node]

    def children(self, node):
        """Return ordered tuple of children of given node."""
        return self._children[node]

    def main_child(self, node):
        """Return child of given node with same argext, or None."""
        if node == self._tip[node]:
            return None
        else:
            return self._children[node][0]

    def lateral(self, node):
        """Return ordered tuple of the given node's lateral children."""
        if node == self.tip(node):
            return self.children(node)
        else:
            return self.children(node)[1:]

    def full(self, node):
        """Return the largest region with same argext as the given node."""
        return self._full[node]

    # generator methods (yielding nodes):

    def path(self, start, istop, step):
        """Yield nodes on a path in the tree."""
        climber = start
        yield climber
        while climber != istop:
            climber = step(climber)
            yield climber

    def root_path(self, node):
        """Yield nodes on the parent path from given node to its root."""
        yield from self.path(node, self.root(node), self.parent)

    def main_path(self, node):
        """Yield nodes on the main_child path from given node to its tip."""
        yield from self.path(node, self.tip(node), self.main_child)

    def subtree(self, localroot=None):
        """Yield all nodes in forest or given subtree."""
        # defaults:
        if localroot is None:
            roots = self.roots()
        else:
            roots = (localroot,)
        for root in roots:
            yield root
            for child in self.children(root):
                yield from self.subtree(child)

    def levels(self, localroot=None, level=0):
        """Yield ordered sequence of (node, level) tuples (roots are zero level)."""
        # defaults:
        if localroot is None:
            roots = self.roots()
        else:
            roots = (localroot,)
        for root in roots:
            yield (root, level)
            for child in self.children(root):
                yield from self.levels(child, level + 1)

    def main_descendants(self, localroot=None):
        """Yield main child nodes in forest or given subtree."""
        # defaults:
        if localroot is None:
            roots = self.roots()
        else:
            roots = (localroot,)
        for root in roots:
            if (mc := self.main_child(root)) is not None:
                yield mc
            for child in self.children(root):
                yield from self.main_descendants(child)

    def lateral_descendants(self, localroot=None):
        """Yield lateral child nodes in forest or given subtree."""
        # defaults:
        if localroot is None:
            roots = self.roots()
        else:
            roots = (localroot,)
        for root in roots:
            yield from self.lateral(root)
            for child in self.children(root):
                yield from self.lateral_descendants(child)

    def full_nodes(self, localroot=None):
        """Yield full nodes in forest or given subtree."""
        # defaults:
        if localroot is None:
            yield from self.roots()
        elif localroot == self.full(localroot):
            yield localroot
        yield from self.lateral_descendants(localroot)

    def leaf_nodes(self, localroot=None):
        """Yield forest or subtree nodes that have no children."""
        return (
            node for node in self.subtree(localroot) if len(self.children(node)) == 0
        )

    def branch_nodes(self, localroot=None):
        """Yield forest or subtree nodes that have two or more children."""
        return (
            node for node in self.subtree(localroot) if len(self.children(node)) > 1
        )

    def linear_nodes(self, localroot=None):
        """Yield forest or subtree nodes that have one child."""
        return (
            node for node in self.subtree(localroot) if len(self.children(node)) == 1
        )

    def size_filter(self, localroot=None, *, maxsize=None):
        """Yield forest or subtree nodes filtered by size."""
        # defaults:
        if localroot is None:
            roots = self.roots()
        else:
            roots = (localroot,)
        if maxsize is None:
            maxsize = 0.2 * max(self.size(root) for root in self.roots())
        # The 'MAXDEEP algorithm' in reverse:
        for root in roots:
            for climber in self.main_path(root):
                if self.size(climber) >= maxsize:
                    for child in self.lateral(climber):
                        yield from self.size_filter(maxsize=maxsize, localroot=child)
                elif (
                    not self.is_nonroot(climber)
                    or self.size(self.parent(climber)) >= maxsize
                ):
                    yield climber
                    break

    def innermost(self, nodes, localroot=None):
        """Yield innermost nodes of the given nodes."""
        # defaults:
        if localroot is None:
            roots = self.roots()
        else:
            roots = (localroot,)
        filter = list(nodes)
        countdown = {
            n: len(self.children(n)) for root in roots for n in self.subtree(root)
        }

        # Recursive "bottom-up" search via parents:
        def _yield_or_propagate(node):
            if node in filter:
                yield node
            elif node not in roots:
                parent = self.parent(node)
                countdown[parent] -= 1
                if countdown[parent] == 0:
                    yield from _yield_or_propagate(parent)

        for root in roots:
            for node in self.leaf_nodes(root):
                yield from _yield_or_propagate(node)

    def outermost(self, nodes, localroot=None):
        """Yield outermost nodes of the given nodes."""
        # defaults:
        if localroot is None:
            roots = self.roots()
        else:
            roots = (localroot,)
        filter = list(nodes)

        # Recursive "top-down" search via children:
        def _yield_or_branch(node):
            if node in filter:
                yield node
            else:
                for child in self.children(node):
                    yield from _yield_or_branch(child)

        for root in roots:
            yield from _yield_or_branch(root)

    # Implementation details (may change):

    def _index(self, node):
        """Return the zero-based index of the given node."""
        # accessing the preserved insertion order:
        return list(self._full).index(node)

    def _find_full_parent(self):
        """Compute attributes: self._full self._parent."""

        def _make_parent(node):
            for child in self.children(node):
                self._parent[child] = node
                if self.tip(node) == self.tip(child):
                    self._full[child] = self._full[node]
                else:
                    self._full[child] = child
                _make_parent(child)

        self._parent = {}
        self._full = {}
        for root in self.roots():
            self._parent[root] = None
            self._full[root] = root
            _make_parent(root)


class HyperTree(Tree):
    """Tree of higher-dimensional regions.

    A HyperTree represents the hierarchical nesting of peak regions
    or valley regions in more dimensions, such as a mountain landscape
    with height z as a function of x and y.

    A HyperTree assumes dimensional decoupling, i.e. the landscape is a
    sum z(x,y) = f(x) + g(y) or a product z(x,y) = f(x) * g(y).

    A HyperTree is constructed as a pair of trees that are
    of type Tree or HyperTree.

    A HyperTree is a kind of product tree, but it is not the Cartesian
    product.

    Notes
    -----
    Background literature for the HyperTree class is
    the subsection titled "2D peaks" in the article [1]_.

    References
    ----------
    .. [1] Tostesen, E. "A stitch in time: Efficient computation of
       genomic DNA melting bubbles." Algorithms for Molecular
       Biology 3, 10 (2008).
       Open access: https://doi.org/10.1186/1748-7188-3-10
    """

    def __init__(self, left_tree, right_tree):
        self.L = left_tree
        self.R = right_tree

    def __contains__(self, pair):
        """Return True if the input is a node in the HyperTree."""
        a, b = pair
        # test if (a, b) is 'parent-above':
        return (
            a == self.L.root() or self.L.size(self.L.parent(a)) > self.R.size(b)
        ) and (b == self.R.root() or self.R.size(self.R.parent(b)) > self.L.size(a))

    def __iter__(self):
        """Iterate over nodes in the HyperTree."""
        yield from self.subtree()

    def __len__(self):
        """Return number of nodes in the HyperTree."""
        return len(list(self.__iter__()))

    def __repr__(self) -> str:
        """Return string that can reconstruct the HyperTree."""
        return f"HyperTree({repr(self.L)}, {repr(self.R)})"

    def root(self):
        """Return the root node of the HyperTree."""
        return (self.L.root(), self.R.root())

    def is_nonroot(self, node):
        """Return True if given node has a parent."""
        a, b = node
        return self.L.is_nonroot(a) or self.R.is_nonroot(b)

    def tip(self, node):
        """Return the given node's tip node."""
        a, b = node
        return (self.L.tip(a), self.R.tip(b))

    def has_children(self, node):
        """Return True if given node has children."""
        a, b = node
        return self.L.has_children(a) or self.R.has_children(b)

    def size(self, node):
        """Return the given node's size."""
        a, b = node
        return max(self.L.size(a), self.R.size(b))

    def parent(self, node):
        """Return the given node's parent or None."""
        a, b = node
        if self.L.is_nonroot(a) and self.R.is_nonroot(b):
            pa, pb = self.L.parent(a), self.R.parent(b)
            if self.L.size(pa) > self.R.size(pb):
                return (a, pb)
            elif self.L.size(pa) < self.R.size(pb):
                return (pa, b)
            elif self.L.size(pa) == self.R.size(pb):
                return (pa, pb)
        elif not self.L.is_nonroot(a) and not self.R.is_nonroot(b):
            return None
        elif self.L.is_nonroot(a):
            # then b is root
            return (self.L.parent(a), b)
        else:
            # then a is root and b nonroot
            return (a, self.R.parent(b))

    def children(self, node):
        """Return the given node's children."""
        a, b = node
        if not self.has_children(node):
            return ()
        elif self.L.size(a) > self.R.size(b):
            return tuple((ca, b) for ca in self.L.children(a))
        elif self.L.size(a) < self.R.size(b):
            return tuple((a, cb) for cb in self.R.children(b))
        elif self.L.size(a) == self.R.size(b):
            return tuple(
                (ca, cb) for ca in self.L.children(a) for cb in self.R.children(b)
            )

    def main_child(self, node):
        """Return the main child (the child that has the same tip)."""
        a, b = node
        if not self.has_children(node):
            return None
        elif self.L.size(a) > self.R.size(b):
            return (self.L.main_child(a), b)
        elif self.L.size(a) < self.R.size(b):
            return (a, self.R.main_child(b))
        elif self.L.size(a) == self.R.size(b):
            return (self.L.main_child(a), self.R.main_child(b))

    def full(self, node):
        """Return the largest node with same tip as given node."""
        climber = node
        while self.is_nonroot(climber) and self.tip(climber) == self.tip(
            nextstep := self.parent(climber)
        ):
            climber = nextstep
        return climber

    def _index(self, node):
        """Return a tuple of (nested) indices for given node."""
        a, b = node
        return self.L._index(a), self.R._index(b)

    def leaf_nodes(self, localroot=None):
        """Yield leaf nodes."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        ra, rb = localroot
        for a in self.L.leaf_nodes(localroot=ra):
            for b in self.R.leaf_nodes(localroot=rb):
                yield a, b

    def size_filter(self, localroot=None, *, maxsize=None):
        """Yield grid nodes in given subtree."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        if maxsize is None:
            maxsize = 0.2 * max(self.L.size(self.L.root()), self.R.size(self.R.root()))
        ra, rb = localroot
        for a in self.L.size_filter(maxsize=maxsize, localroot=ra):
            for b in self.R.size_filter(maxsize=maxsize, localroot=rb):
                yield a, b

    def from_peaks(self):
        """Return that it is NotImplemented."""
        return NotImplemented

    def from_valleys(self):
        """Return that it is NotImplemented."""
        return NotImplemented

    def from_levels(self):
        """Return that it is NotImplemented."""
        return NotImplemented

    def as_dict_of_dicts(self):
        """Return that it is NotImplemented."""
        return NotImplemented

    def set_nodes(self):
        """Return that it is NotImplemented."""
        return NotImplemented

    def _find_full(self):
        """Return that it is NotImplemented."""
        return NotImplemented


class HyperForest(Forest):
    """Forest of trees of higher-dimensional regions.

    A HyperForest represents the hierarchical nesting of peak regions
    or valley regions in more dimensions, for example, a mountain landscape
    with height z as a function of x and y.

    A HyperForest assumes dimensional decoupling, for example, a
    landscape z(x, y) that is a sum f(x) + g(y) or product f(x) * g(y).

    A HyperForest is constructed as a pair of forests that are
    of type Forest or HyperForest.

    A HyperForest is a kind of product tree, but it is not a Cartesian
    product.

    Notes
    -----
    Background literature for the HyperForest class is
    the subsection titled "2D peaks" in the article [1]_.

    References
    ----------
    .. [1] Tostesen, E. "A stitch in time: Efficient computation of
       genomic DNA melting bubbles." Algorithms for Molecular
       Biology 3, 10 (2008).
       Open access: https://doi.org/10.1186/1748-7188-3-10
    """

    # dunder methods:

    def __init__(self, left_forest, right_forest):
        """Initialize HyperForest from pair of Forest or HyperForest objects."""
        self.L = left_forest
        self.R = right_forest

    def __contains__(self, pair):
        """Return True if the input is a node in the HyperForest."""
        a, b = pair
        # test if (a, b) is 'parent-above' or 'below-leaf':
        return (
            a in self.L.roots()
            or self.L.size(self.L.parent(a)) > self.R.size(b)
            or (not self.R.has_children(b) and self.R.size(b) >= self.L.size(a))
        ) and (
            b in self.R.roots()
            or self.R.size(self.R.parent(b)) > self.L.size(a)
            or (not self.L.has_children(a) and self.L.size(a) >= self.R.size(b))
        )

    def __iter__(self):
        """Iterate over nodes in the HyperForest."""
        yield from self.subtree()

    def __len__(self):
        """Return number of nodes in the HyperForest."""
        return len(list(self.__iter__()))

    def __repr__(self) -> str:
        """Return string that can reconstruct the HyperForest."""
        return f"HyperForest({repr(self.L)}, {repr(self.R)})"

    # node methods:

    def roots(self):
        """Return tuple of roots of trees in HyperForest."""
        return tuple((a, b) for a in self.L.roots() for b in self.R.roots())

    def is_nonroot(self, node):
        """Return True if given node has a parent."""
        a, b = node
        return self.L.is_nonroot(a) or self.R.is_nonroot(b)

    def tip(self, node):
        """Return the given node's tip node."""
        a, b = node
        climber = (self.L.tip(a), self.R.tip(b))
        while climber not in self:
            climber = self.parent(climber)
        return climber

    def has_children(self, node):
        """Return True if given node has children."""
        a, b = node
        return self.L.has_children(a) or self.R.has_children(b)

    def size(self, node):
        """Return the given node's size."""
        a, b = node
        if (not self.L.has_children(a) and self.L.size(a) > self.R.size(b)) or (
            not self.R.has_children(b) and self.L.size(a) < self.R.size(b)
        ):
            return min(self.L.size(a), self.R.size(b))
        else:
            return max(self.L.size(a), self.R.size(b))

    def parent(self, node):
        """Return the given node's parent or None."""
        a, b = node
        if self.L.is_nonroot(a) and self.R.is_nonroot(b):
            pa, pb = self.L.parent(a), self.R.parent(b)
            if self.L.size(pa) > self.R.size(pb):
                return (a, pb)
            elif self.L.size(pa) < self.R.size(pb):
                return (pa, b)
            elif self.L.size(pa) == self.R.size(pb):
                return (pa, pb)
        elif not self.L.is_nonroot(a) and not self.R.is_nonroot(b):
            return None
        elif self.L.is_nonroot(a):
            # then b is root
            return (self.L.parent(a), b)
        else:
            # then a is root and b nonroot
            return (a, self.R.parent(b))

    def children(self, node):
        """Return the given node's children."""
        a, b = node
        if self.L.has_children(a) and self.R.has_children(b):
            if self.L.size(a) > self.R.size(b):
                return tuple((ca, b) for ca in self.L.children(a))
            elif self.L.size(a) < self.R.size(b):
                return tuple((a, cb) for cb in self.R.children(b))
            elif self.L.size(a) == self.R.size(b):
                return tuple(
                    (ca, cb) for ca in self.L.children(a) for cb in self.R.children(b)
                )
        elif not self.L.has_children(a) and not self.R.has_children(b):
            return ()
        elif self.L.has_children(a):
            return tuple((ca, b) for ca in self.L.children(a))
        else:
            return tuple((a, cb) for cb in self.R.children(b))

    def main_child(self, node):
        """Return the main child (the child that has the same tip)."""
        a, b = node
        if (children := self.children(node)) == ():
            return None
        else:
            a1, b1 = children[0]
            if self.L.tip(a1) == self.L.tip(a) and self.R.tip(b1) == self.R.tip(b):
                return a1, b1
            else:
                return None

    def full(self, node):
        """Return the largest node with same tip as given node."""
        a, b = node
        climber = (self.L.full(a), self.R.full(b))
        while climber not in self:
            climber = self.main_child(climber)
        return climber

    # generator methods (yielding nodes):

    def leaf_nodes(self, localroot=None):
        """Yield leaf nodes."""
        # defaults:
        if localroot is None:
            roots = self.roots()
        else:
            roots = (localroot,)
        for ra, rb in roots:
            for a in self.L.leaf_nodes(localroot=ra):
                for b in self.R.leaf_nodes(localroot=rb):
                    yield a, b

    def size_filter(self, localroot=None, *, maxsize=None):
        """Yield grid nodes in given subtree."""
        # defaults:
        if localroot is None:
            roots = self.roots()
        else:
            roots = (localroot,)
        if maxsize is None:
            maxsize = 0.2 * max(self.size(root) for root in self.roots())
        for ra, rb in roots:
            for a in self.L.size_filter(maxsize=maxsize, localroot=ra):
                for b in self.R.size_filter(maxsize=maxsize, localroot=rb):
                    yield a, b

    # Implementation details (may change):

    def _index(self, node):
        """Return a tuple of (nested) indices for given node."""
        a, b = node
        return self.L._index(a), self.R._index(b)

    # Forest methods that are not HyperForest methods:

    def from_peaks(self):
        """Return that it is NotImplemented."""
        return NotImplemented

    def from_valleys(self):
        """Return that it is NotImplemented."""
        return NotImplemented

    def from_levels(self):
        """Return that it is NotImplemented."""
        return NotImplemented

    def as_dict_of_dicts(self):
        """Return that it is NotImplemented."""
        return NotImplemented

    def set_nodes(self):
        """Return that it is NotImplemented."""
        return NotImplemented

    def _find_full_parent(self):
        """Return that it is NotImplemented."""
        return NotImplemented

    def __sub__(self):
        """Return that it is NotImplemented."""
        return NotImplemented

    def __and__(self):
        """Return that it is NotImplemented."""
        return NotImplemented

    def __or__(self):
        """Return that it is NotImplemented."""
        return NotImplemented

    def __xor__(self):
        """Return that it is NotImplemented."""
        return NotImplemented
