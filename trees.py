#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Python module for trees of peaks in numeric data.

This module contains algorithms for building trees
representing the hierarchical nesting of peak and subpeak regions
in one-dimensional or higher-dimensional data.
It also contains methods for searching trees and selecting peaks.

Requires Python 3.7+

Created on Wed Mar 24 14:49:04 2021.

@author: Eivind Tostesen

"""


from operator import attrgetter


# Functions:


def tree_from_peak_objects(peaks, presorted=True):
    """Return (parent, root, children, top) from peaks having attrs: start, end, min, max."""
    parent = {}
    children = {}
    top = {}
    in_spe = []
    if not presorted:
        peaks = list(peaks)
        # Order same as given by 'peaks' function:
        peaks.sort(key=attrgetter("min"), reverse=True)
        peaks.sort(key=attrgetter("end"))
    for p in peaks:
        children[p] = []
        while in_spe and p.start <= in_spe[-1].start:
            c = in_spe.pop()
            children[p].append(c)
            parent[c] = p
        children[p].sort(key=attrgetter('start'))
        children[p].sort(key=attrgetter('max'), reverse=True)
        children[p] = tuple(children[p])
        top[p] = top[children[p][0]] if children[p] else p
        in_spe.append(p)
    root = in_spe.pop()
    parent[root] = None
    return parent, root, children, top


# Classes:


class PeakTree:
    """Tree of peaks in univariate data.

    A PeakTree represents the hierarchical nesting of
    peak and subpeak regions in 1D data such as a sequence of numbers,
    a time series, a function y(x) or other univariate data.

    A PeakTree is initialized with an iterable of peak objects
    that have attributes start, end, min, max. The peak objects
    must be unique hashable objects to be used as dictionary keys.

    Notes
    -----
    Background literature for the PeakTree class is
    the subsection titled "1D peaks" in the article [1]_.

    References
    ----------
    .. [1] Tostesen, E. "A stitch in time: Efficient computation of
       genomic DNA melting bubbles." Algorithms for Molecular
       Biology 3, 10 (2008).
       Open access: https://doi.org/10.1186/1748-7188-3-10
    """

    @classmethod
    def from_peak_objects(cls, peaks, presorted=True):
        """Build a PeakTree and compute its data attributes."""
        obj = cls.__new__(cls)
        obj._parent, obj._root, obj._children, obj._top = tree_from_peak_objects(peaks, presorted=presorted)
        obj._find_full()
        return obj

    def __init__(self, data):
        pass  # TODO

    def __contains__(self, node):
        """Return True if the input is a node in the PeakTree."""
        return node in self._full

    def __iter__(self):
        """Iterate over nodes in the PeakTree."""
        return iter(self._full)

    def __len__(self):
        """Return number of nodes in the PeakTree."""
        return len(self._full)

    def __matmul__(self, other):
        """Return product with other tree (PeakTree or HyperPeakTree)."""
        return HyperPeakTree(self, other)

    def __str__(self):
        """Return printable tree structure."""
        indent = "| "
        return "\n".join([level * indent + str(node) for node, level in self.levels().items()])

    def as_string(self, localroot):
        """Return printable subtree structure."""
        str = "# Notation: <high> /& <low>/ => <parent>\n"
        for full in self.full_nodes(localroot):
            if not self.has_children(full):
                continue
            for node in self.path(self.top(full), self.high(full), self.parent):
                str += f"{node}"
                if len(self.children(self.parent(node))) == 2:
                    str += f" /& {self.low(self.parent(node))[0]}/"
                if len(self.children(self.parent(node))) > 2:
                    str += f" /& {self.low(self.parent(node))}/"
                str += " => "
            str += f"{full}.\n"
        return str

    def as_dict_of_dicts(self):
        """Return data attributes as a dict of dicts."""
        return {
            "_parent": self._parent,
            "_children": self._children,
            "_top": self._top,
            "_full": self._full,
            "_root": self._root,
        }

    def set_nodes(self, changes={}):
        """Replace PeakTree nodes by using given mapping."""

        def new(node):
            return changes[node] if node in changes else node

        self._top = dict((new(x), new(y)) for (x, y) in self._top.items())
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
        """Return the root node of the PeakTree."""
        return self._root

    def is_nonroot(self, node):
        """Return True if the input node has a parent."""
        return node != self._root

    def top(self, node):
        """Return the highest leaf node in the node's subtree."""
        return self._top[node]

    def has_children(self, node):
        """Return True if the input node has subpeaks."""
        return node != self._top[node]

    def height(self, node):
        """Return the height at the top of the input peak."""
        return node.max

    def base_height(self, node):
        """Return the height at the base of the input peak."""
        return node.min

    def size(self, node):
        """Return vertical distance between top and base."""
        return node.size

    def parent(self, node):
        """Return the input peak merged with its neighboring peak."""
        return self._parent[node]

    def children(self, node):
        """Return tuple of children (subpeaks) ordered after height."""
        return self._children[node]

    def high(self, node):
        """Return the highest of the input peak's subpeaks."""
        return self.children(node)[0]

    def low(self, node):
        """Return tuple of the input peak's lower subpeaks."""
        return self.children(node)[1:]

    def full(self, node):
        """Return the largest peak with same top as the input peak."""
        return self._full[node]

    def _index(self, node):
        """Return the zero-based index of the input node."""
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
        """Yield nodes on path from a node to its root."""
        yield from self.path(node, self.root(), self.parent)

    def top_path(self, node):
        """Yield nodes on path from a node to its top."""
        yield from self.path(node, self.top(node), self.high)

    def subtree(self, localroot=None):
        """Yield all nodes in the input node's subtree."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        yield localroot
        for child in self.children(localroot):
            yield from self.subtree(child)

    def high_descendants(self, localroot=None):
        """Yield high children in the input node's subtree."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        if self.has_children(localroot):
            yield self.high(localroot)
            for child in self.children(localroot):
                yield from self.high_descendants(child)

    def low_descendants(self, localroot=None):
        """Yield lower children in the input node's subtree."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        if self.has_children(localroot):
            yield from self.low(localroot)
            for child in self.children(localroot):
                yield from self.low_descendants(child)

    def full_nodes(self, localroot=None):
        """Yield full nodes in the input node's subtree."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        if localroot == self.full(localroot):
            yield localroot
        yield from self.low_descendants(localroot)

    def leaf_nodes(self, localroot=None):
        """Yield leaf nodes (local maxima) in (sub)tree."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        return (
            node for node in self.subtree(localroot) if len(self.children(node)) == 0
        )

    def branch_nodes(self, localroot=None):
        """Yield nodes in (sub)tree with two or more children."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        return (
            node for node in self.subtree(localroot) if len(self.children(node)) > 1
        )

    def linear_nodes(self, localroot=None):
        """Yield nodes in (sub)tree with one child."""
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
            climber = self.top(localroot)
            while self.size(self.parent(climber)) < maxsize:
                climber = self.parent(climber)
            yield climber
            while climber != localroot:
                climber = self.parent(climber)
                for child in self.low(climber):
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
            yield from self.low_descendants()

        self._full = {
            node: full for full in fullnodes() for node in self.top_path(full)
        }


class HyperPeakTree(PeakTree):
    """Tree of higher-dimensional peaks.

    A HyperPeakTree represents the hierarchical nesting of peaks
    and subpeaks in more dimensions, such as a mountain landscape
    with height z as a function of x and y.

    A HyperPeakTree assumes dimensional decoupling, i.e. the landscape is a
    sum z(x,y) = f(x) + g(y) or a product z(x,y) = f(x) * g(y).

    A HyperPeakTree is constructed as a pair of lower-dimensional trees
    of type PeakTree or HyperPeakTree.

    A HyperPeakTree is a kind of product tree, but it is not the Cartesian
    product.

    Notes
    -----
    Background literature for the HyperPeakTree class is
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
        """Return True if the input is a node in the HyperPeakTree."""
        a, b = pair
        # test if (a, b) is 'sigma-above':
        return (
            a == self.L.root() or self.L.size(self.L.parent(a)) > self.R.size(b)
        ) and (b == self.R.root() or self.R.size(self.R.parent(b)) > self.L.size(a))

    def __iter__(self):
        """Iterate over nodes in the HyperPeakTree."""
        yield from self.subtree()

    def __len__(self):
        """Return number of nodes in the HyperPeakTree."""
        return len(list(self.__iter__()))

    def root(self):
        """Return the root node of the HyperPeakTree."""
        return (self.L.root(), self.R.root())

    def is_nonroot(self, node):
        """Return True if given node has a parent."""
        a, b = node
        return self.L.is_nonroot(a) or self.R.is_nonroot(b)

    def top(self, node):
        """Return the given node's top (a leaf node)."""
        a, b = node
        return (self.L.top(a), self.R.top(b))

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

    def high(self, node):
        """Return the child that has the same top."""
        a, b = node
        if self.L.size(a) > self.R.size(b):
            return (self.L.high(a), b)
        elif self.L.size(a) < self.R.size(b):
            return (a, self.R.high(b))
        elif self.L.size(a) == self.R.size(b):
            return (self.L.high(a), self.R.high(b))

    def full(self, node):
        """Return the largest node with same top as given node."""
        climber = node
        while self.is_nonroot(climber):
            nextstep = self.parent(climber)
            if self.top(nextstep) == self.top(climber):
                climber = nextstep
            else:
                break
        return climber

    def _index(self, node):
        """Return a tuple of (nested) indices for given node."""
        a, b = node
        return self.L._index(a), self.R._index(b)

    def leaf_nodes(self):
        """Yield leaf nodes (local maxima) in (sub)tree."""
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

    def from_peak_objects(self):
        """Return that it is NotImplemented."""
        return NotImplemented
    
    def as_dict_of_dicts(self):
        """Return that it is NotImplemented."""
        return NotImplemented

    def set_nodes(self):
        """Return that it is NotImplemented."""
        return NotImplemented

    def height(self, frame):
        """Return that height is NotImplemented."""
        return NotImplemented

    def base_height(self, node):
        """Return that base_height is NotImplemented."""
        return NotImplemented

    def _find_full(self):
        """Return that it is NotImplemented."""
        return NotImplemented
