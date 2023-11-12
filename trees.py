#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Python module for trees of peak regions or valley regions.

This module contains algorithms for building trees that
represent the hierarchical nesting of regions and subregions
containing peaks or valleys in numeric
one-dimensional or higher-dimensional data.
The tree classes provide methods for searching and selecting regions.

Requires Python 3.7+

Created on Wed Mar 24 14:49:04 2021.

@author: Eivind Tostesen

"""


from operator import attrgetter
from utilities import pairwise


# Functions:


def tree_from_peaks(
    peaks,
    presorted=True,
    getstart=attrgetter("start"),
    getend=attrgetter("end"),
    getmin=attrgetter("min"),
    getmax=attrgetter("max"),
):
    """Return (parent, root, children, core) from peak regions having start, end, min, max."""
    parent = {}
    children = {}
    core = {}
    in_spe = []
    if not presorted:
        peaks = list(peaks)
        # Order same as given by 'find_peaks' function:
        peaks.sort(key=getmin, reverse=True)
        peaks.sort(key=getend)
    for p in peaks:
        children[p] = []
        while in_spe and getstart(p) <= getstart(in_spe[-1]):
            c = in_spe.pop()
            children[p].append(c)
            parent[c] = p
        children[p].sort(key=getstart)
        children[p].sort(key=getmax, reverse=True)
        children[p] = tuple(children[p])
        core[p] = core[children[p][0]] if children[p] else p
        in_spe.append(p)
    root = in_spe.pop()
    parent[root] = None
    return parent, root, children, core


# Classes:


class Tree:
    """Tree of regions in univariate data.

    A Tree represents the hierarchical nesting of
    peak or valley regions in 1D data such as a sequence of numbers,
    a time series, a function y(x) or other univariate data.

    A Tree is initialized with an iterable of regions
    that have start, end, min, max. The regions
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

    @classmethod
    def from_peaks(cls, peaks, **kwargs):
        """Return new Tree from iterable of peak regions."""
        obj = cls.__new__(cls)
        obj._parent, obj._root, obj._children, obj._core = tree_from_peaks(
            peaks, **kwargs
        )
        obj._find_full()
        return obj

    @classmethod
    def from_valleys(
        cls,
        valleys,
        presorted=True,
        getstart=attrgetter("start"),
        getend=attrgetter("end"),
        getmin=lambda p: -p.max,
        getmax=lambda p: -p.min,
    ):
        """Return new Tree from iterable of valley regions."""
        obj = cls.__new__(cls)
        obj._parent, obj._root, obj._children, obj._core = tree_from_peaks(
            valleys,
            presorted=presorted,
            getstart=getstart,
            getend=getend,
            getmin=getmin,
            getmax=getmax,
        )
        obj._find_full()
        return obj

    @classmethod
    def from_levels(cls, levels):
        """Return new Tree from other tree's levels() output."""

        def leaf_and_core(node):
            children[node] = []
            for n in obj.path(node, obj._full[node], obj.parent):
                obj._core[n] = node

        obj = cls.__new__(cls)
        obj._parent = {}
        children = {}
        obj._core = {}
        obj._full = {}
        for (A, a), (B, b) in pairwise(levels.items()):
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
                leaf_and_core(A)  # A is a leaf and core
            obj._parent[B] = stack[-1]
            stack.append(B)
        leaf_and_core(B)  # the last B is a leaf and core
        obj._children = {p: tuple(c) for p, c in children.items()}
        return obj

    def __init__(self, data):
        pass  # TODO

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
        return f"Tree.from_levels({repr(self.levels())})"

    def __str__(self):
        """Return string with tree in indented list notation."""
        indent = "| "
        return "\n".join(
            [level * indent + str(node) for node, level in self.levels().items()]
        )

    def as_dict_of_dicts(self):
        """Return data attributes as a dict of dicts."""
        return {
            "_parent": self._parent,
            "_children": self._children,
            "_core": self._core,
            "_full": self._full,
            "_root": self._root,
        }

    def set_nodes(self, changes={}):
        """Replace Tree nodes by using given mapping."""

        def new(node):
            return changes[node] if node in changes else node

        self._core = dict((new(x), new(y)) for (x, y) in self._core.items())
        self._full = dict((new(x), new(y)) for (x, y) in self._full.items())
        self._parent = dict((new(x), new(y)) for (x, y) in self._parent.items())
        self._children = dict(
            (new(x), tuple(new(z) for z in y)) for (x, y) in self._children.items()
        )
        self._root = new(self._root)
        return None

    def levels(self):
        """Return ordered dict of node:level pairs (root is zero level)."""
        levels = {}
        for node in self.subtree():
            if self.is_nonroot(node):
                levels[node] = 1 + levels[self.parent(node)]
            else:
                levels[node] = 0
        return levels

    def root(self):
        """Return the root node of the Tree."""
        return self._root

    def is_nonroot(self, node):
        """Return True if the given node has a parent."""
        return node != self._root

    def core(self, node):
        """Return the given node's core region (main leaf in subtree)."""
        return self._core[node]

    def has_children(self, node):
        """Return True if the given node has sub regions."""
        return node != self._core[node]

    def max(self, node):
        """Return the max value in the given region."""
        return node.max

    def min(self, node):
        """Return the min value in the given region."""
        return node.min

    def size(self, node):
        """Return vertical size of given node (max minus min)."""
        return node.size

    def parent(self, node):
        """Return the parent (containing region) or None."""
        return self._parent[node]

    def children(self, node):
        """Return ordered tuple of children."""
        return self._children[node]

    def main(self, node):
        """Return the given node's main child (first in tuple)."""
        return self.children(node)[0]

    def lateral(self, node):
        """Return the given node's lateral children (tuple except first)."""
        return self.children(node)[1:]

    def full(self, node):
        """Return the largest region with same core as the given node."""
        return self._full[node]

    def _index(self, node):
        """Return the zero-based index of the given node."""
        # accessing the preserved insertion order:
        return list(self._full).index(node)

    # public recursive algorithms:

    def path(self, start, stop, step):
        """Yield nodes on a path in the tree."""
        climber = start
        yield climber
        while climber != stop:
            climber = step(climber)
            yield climber

    def root_path(self, node):
        """Yield nodes on the path from given node to its root."""
        yield from self.path(node, self.root(), self.parent)

    def core_path(self, node):
        """Yield nodes on the path from given node to its core."""
        yield from self.path(node, self.core(node), self.main)

    def subtree(self, localroot=None):
        """Yield all nodes in the given subtree."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        yield localroot
        for child in self.children(localroot):
            yield from self.subtree(child)

    def main_descendants(self, localroot=None):
        """Yield main nodes in the given subtree."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        if self.has_children(localroot):
            yield self.main(localroot)
            for child in self.children(localroot):
                yield from self.main_descendants(child)

    def lateral_descendants(self, localroot=None):
        """Yield lateral nodes in the given subtree."""
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

    def filter(self, *, maxsize=None, localroot=None):
        """Yield subtree nodes filtered by size."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        if maxsize is None:
            maxsize = 0.2 * self.size(self.root())
        # The 'maxdeep' algorithm:
        if self.size(localroot) < maxsize:
            yield localroot
        else:
            climber = self.core(localroot)
            while self.size(self.parent(climber)) < maxsize:
                climber = self.parent(climber)
            yield climber
            while climber != localroot:
                climber = self.parent(climber)
                for child in self.lateral(climber):
                    yield from self.filter(maxsize=maxsize, localroot=child)

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
            node: full for full in fullnodes() for node in self.core_path(full)
        }


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
        # test if (a, b) is 'sigma-above':
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

    def core(self, node):
        """Return the given node's core (main leaf node)."""
        a, b = node
        return (self.L.core(a), self.R.core(b))

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
        if self.L.size(a) > self.R.size(b):
            return tuple((ca, b) for ca in self.L.children(a))
        elif self.L.size(a) < self.R.size(b):
            return tuple((a, cb) for cb in self.R.children(b))
        elif self.L.size(a) == self.R.size(b):
            return tuple(
                (ca, cb) for ca in self.L.children(a) for cb in self.R.children(b)
            )

    def main(self, node):
        """Return the main child (the child that has the same core)."""
        a, b = node
        if self.L.size(a) > self.R.size(b):
            return (self.L.main(a), b)
        elif self.L.size(a) < self.R.size(b):
            return (a, self.R.main(b))
        elif self.L.size(a) == self.R.size(b):
            return (self.L.main(a), self.R.main(b))

    def full(self, node):
        """Return the largest node with same core as given node."""
        climber = node
        while self.is_nonroot(climber):
            nextstep = self.parent(climber)
            if self.core(nextstep) == self.core(climber):
                climber = nextstep
            else:
                break
        return climber

    def _index(self, node):
        """Return a tuple of (nested) indices for given node."""
        a, b = node
        return self.L._index(a), self.R._index(b)

    def leaf_nodes(self):
        """Yield leaf nodes."""
        for a in self.L.leaf_nodes():
            for b in self.R.leaf_nodes():
                yield a, b

    def filter(self, *, maxsize=None, localroot=None):
        """Yield grid nodes in given subtree."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        if maxsize is None:
            maxsize = 0.2 * max(self.L.size(self.L.root()), self.R.size(self.R.root()))
        ra, rb = localroot
        for a in self.L.filter(maxsize=maxsize, localroot=ra):
            for b in self.R.filter(maxsize=maxsize, localroot=rb):
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

    def max(self, frame):
        """Return that max is NotImplemented."""
        return NotImplemented

    def min(self, node):
        """Return that min is NotImplemented."""
        return NotImplemented

    def _find_full(self):
        """Return that it is NotImplemented."""
        return NotImplemented
