#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Python module for finding peaks in numeric data.

This module contains algorithms for building binary trees
representing the hierarchical structure of peaks and subpeaks
in one-dimensional or higher-dimensional data.
It also contains methods for searching the tree and selecting peaks.

Created on Wed Mar 24 14:49:04 2021.

@author: Eivind Tostesen

"""


from itertools import chain


# Iteration tools:


def _forward_backward(iterable):
    """Return pair of iterators: (forward, backward)."""
    data = list(iterable)
    return iter(data), reversed(data)


def _pairwise(iterable):
    """Yield nearest neighbor pairs."""
    it = iter(iterable)
    a = next(it)
    for b in it:
        yield a, b
        a = b


def _tripletwise(iterable):
    """Yield nearest neighbor triplets."""
    it = iter(iterable)
    a, b = next(it), next(it)
    for c in it:
        yield a, b, c
        a, b = b, c


# Functions:


def filter_local_extrema(datapoints):
    """Let only maxima/minima pass from a stream of points."""
    previous = None
    for (x1, e1), (x2, e2) in _pairwise(datapoints):
        if e1 == e2:
            continue
        is_uphill = e2 > e1
        if is_uphill == previous:
            continue
        yield x1, e1
        previous = is_uphill
    else:
        yield x2, e2


def peak_locations(peakpoints, curvepoints, revpeaks=None, revcurve=None):
    """Return dict of (start, end)-locations for the input points."""
    # defaults:
    if revpeaks is None:
        peakpoints, revpeaks = _forward_backward(peakpoints)
    if revcurve is None:
        curvepoints, revcurve = _forward_backward(curvepoints)

    def trace(curve, peaks):
        # flow: input -> peaklist -> level -> outputdict
        peaklist = list(peaks)
        level = {}
        outputdict = {}
        previousx = None
        for x, e in curve:
            for p in [p for p in level if level[p] > e]:
                outputdict[p] = previousx
                del level[p]
            if peaklist and peaklist[-1][0] == x:
                level.update([peaklist.pop()])
            previousx = x
        for p in level:
            outputdict[p] = previousx
        return outputdict

    # compute starts and ends:
    startdict = trace(revcurve, peakpoints)
    enddict = trace(curvepoints, revpeaks)
    # return location intervals:
    return {x: (startdict[x], enddict[x]) for x in enddict}


# Classes:


class PeakTree:
    """Binary tree of peaks in univariate data.

    A PeakTree represents the hierarchical structure of
    peaks and subpeaks in data such as a sequence of numbers,
    a time series, a function y(x) or other univariate data.

    A PeakTree is initialized with an iterable of (label, value)-tuples,
    representing the local extrema in the data set
    (i.e. the alternating maximum and minimum points).

    The values must be numeric and the labels must be
    unique hashable objects for use as dictionary keys.

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

    def __init__(self, data):
        """Build a PeakTree and compute its data attributes."""
        self.data = dict(data)
        # requires preserved insertion order:
        nodes = list(self.data)
        # remove possible minimum at the beginning:
        if self.data[nodes[0]] < self.data[nodes[1]]:
            del self.data[nodes[0]]
            del nodes[0]
        # remove possible minimum at the end:
        if self.data[nodes[-1]] < self.data[nodes[-2]]:
            del self.data[nodes[-1]]
            del nodes[-1]
        # compute data attributes:
        self._find_parent_and_root()
        self._find_mode_high_low()
        self._find_full()

    def __contains__(self, node):
        """Return True if the input is a node in the PeakTree."""
        return node in self.data

    def __iter__(self):
        """Iterate over nodes in the PeakTree."""
        return iter(self.data)

    def __len__(self):
        """Return number of nodes in the PeakTree."""
        return len(self.data)

    def __matmul__(self, other):
        """Return product with other tree (PeakTree or FrameTree)."""
        return FrameTree(self, other)

    def __str__(self):
        """Return printable tree structure."""
        return self.as_string(self.root())

    def as_string(self, localroot):
        """Return printable subtree structure."""
        str = "# Notation: <high> /& <low>/ => <parent>\n"
        for full in chain([localroot],
                          self.low_descendants(localroot)
                          ):
            if not self.has_children(full):
                continue
            for node in self.path(self.mode(full),
                                  self.high(full),
                                  self.parent
                                  ):
                str += f'{node} /& {self.low(self.parent(node))}/ => '
            str += f'{full}\n'
        return str

    def as_dict_of_dicts(self):
        """Return data attributes as a dict of dicts."""
        return {
            "elevation": self.data,
            "parent": self._parent,
            "high": self._high,
            "low": self._low,
            "mode": self._mode,
            "full": self._full,
        }

    def maxima(self):
        """Return ordered list of leaf nodes (local maxima)."""
        return list(self.data)[0::2]

    def minima(self):
        """Return ordered list of branch nodes (local minima)."""
        return list(self.data)[1::2]

    def root(self):
        """Return the root node of the PeakTree."""
        return self._root

    def is_nonroot(self, node):
        """Return True if the input peak has a parent."""
        return node != self._root

    def mode(self, node):
        """Return the input peak's top (highest leaf node)."""
        return self._mode[node]

    def has_children(self, node):
        """Return True if the input peak has subpeaks."""
        return node != self._mode[node]

    def height(self, node):
        """Return the height at the mode of the input peak."""
        return self.data[self.mode(node)]

    def base_height(self, node):
        """Return the height at the base of the input peak."""
        return self.data[node]

    def size(self, node):
        """Return vertical distance between mode and base."""
        return self.height(node) - self.base_height(node)

    def parent(self, node):
        """Return the input peak merged with its neighboring peak."""
        return self._parent[node]

    def high(self, node):
        """Return the highest of the input peak's two subpeaks."""
        return self._high[node]

    def low(self, node):
        """Return the lowest of the input peak's two subpeaks."""
        return self._low[node]

    def full(self, node):
        """Return the largest peak with same mode as the input peak."""
        return self._full[self.mode(node)]

    def index(self, node):
        """Return the zero-based index of the input peak."""
        # accessing the preserved insertion order:
        return list(self.data).index(node)

    # public recursive algorithms:

    def path(self, start, stop, step):
        """Yield nodes on a path in the tree."""
        climber = start
        yield climber
        while (climber != stop):
            climber = step(climber)
            yield climber

    def root_path(self, node):
        """Yield nodes on the path from the input node to the root."""
        yield from self.path(node, self.root(), self.parent)

    def subtree(self, localroot=None):
        """Yield all nodes in the input node's subtree."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        yield localroot
        yield from self.high_descendants(localroot)
        yield from self.low_descendants(localroot)

    def high_descendants(self, localroot=None):
        """Yield high children in the input node's subtree."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        if self.has_children(localroot):
            yield self.high(localroot)
            yield from self.high_descendants(self.high(localroot))
            yield from self.high_descendants(self.low(localroot))

    def low_descendants(self, localroot=None):
        """Yield low children in the input node's subtree."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        if self.has_children(localroot):
            yield self.low(localroot)
            yield from self.low_descendants(self.low(localroot))
            yield from self.low_descendants(self.high(localroot))

    def full_to_mode_path(self, node):
        """Yield path nodes from full to mode."""
        yield from self.path(self.full(node), self.mode(node), self.high)

    def filter(self, *, maxsize, localroot=None):
        """Yield subtree nodes filtered by size."""
        # defaults:
        if localroot is None:
            localroot = self.root()

        def height_above_base(x):
            while self.is_nonroot(x) \
                    and self.data[x] == self.data[self.parent(x)]:
                x = self.parent(x)
            return self.size(x)

        # The 'maxdeep' algorithm:
        if height_above_base(localroot) < maxsize:
            yield localroot
        else:
            climber = self.mode(localroot)
            while height_above_base(self.parent(climber)) < maxsize:
                climber = self.parent(climber)
            yield climber
            while climber != localroot:
                climber = self.parent(climber)
                yield from self.filter(maxsize=maxsize, localroot=self.low(climber))

    # Initialization algorithms:

    def _find_parent_and_root(self):
        """Compute attributes: self._parent and self._root."""
        self._parent = {}
        children_in_spe = []
        children = []
        for peak in self:
            while children_in_spe:
                if self.data[children_in_spe[-1]] < self.data[peak]:
                    break
                else:
                    children.append(children_in_spe.pop())
            if children:
                child = children.pop()
                self._parent[child] = peak
                while children:
                    grandchild = children.pop()
                    self._parent[grandchild] = child
                    child = grandchild
            children_in_spe.append(peak)
        child = children_in_spe.pop()
        while children_in_spe:
            parent = children_in_spe.pop()
            self._parent[child] = parent
            child = parent
        self._root = parent

    def _find_mode_high_low(self):
        """Compute attributes: self._mode, self._high and self._low."""
        self._high = {}
        self._low = {}

        def propagate_mode(peak):
            parent = self.parent(peak)
            if parent not in self._high:
                self._high[parent] = peak
            else:
                self._low[parent] = peak
                hf = self.height(self.high(parent))
                hm = self.height(self.low(parent))
                if hf < hm:
                    self._low[parent] = self._high[parent]
                    self._high[parent] = peak
                elif hf == hm:
                    # low goes left, high goes right:
                    if self.index(peak) > self.index(self.high(parent)):
                        self._low[parent] = self._high[parent]
                        self._high[parent] = peak
                self._mode[parent] = self.mode(self.high(parent))
                if self.is_nonroot(parent):
                    propagate_mode(parent)

        self._mode = {peak: peak for peak in self.maxima()}
        for peak in self.maxima():
            propagate_mode(peak)

    def _find_full(self):
        """Compute attribute: self._full."""
        self._full = {self.mode(peak): peak for peak in self._low.values()}
        self._full[self.mode(self._root)] = self._root


class FrameTree(PeakTree):
    """Binary tree of higher-dimensional peaks.

    A FrameTree represents the hierarchical structure of peaks
    and subpeaks in more dimensions, like a "mountain" landscape
    with height z as a function of x and y. It assumes dimensional
    decoupling, i.e. the height is a sum z(x,y) = f(x) + g(y)
    or a product z(x,y) = f(x) * g(y).

    A Frametree is constructed as a pair of lower-dimensional
    trees (of type PeakTree or Frametree). It is a kind of
    product tree, but it is not the Cartesian product.

    Notes
    -----
    Background literature for the FrameTree class is
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

    def __contains__(self, frame):
        """Return True if the input is a node in the FrameTree."""
        a, b = frame
        # test if frame is 'sigma-above':
        return (
            (a == self.L.root()
             or
             self.L.size(self.L.parent(a)) > self.R.size(b)
             )
            and
            (b == self.R.root()
             or
             self.R.size(self.R.parent(b)) > self.L.size(a)
             )
        )

    def __iter__(self):
        """Iterate over nodes in the FrameTree."""
        yield from self.subtree()

    def __len__(self):
        """Return number of nodes in the FrameTree."""
        return len(list(self.subtree()))

    def maxima(self):
        """Yield leaf nodes (local maxima)."""
        for a in self.L.maxima():
            for b in self.R.maxima():
                yield a, b

    def minima(self):
        """Yield local minima."""
        for a in self.L.minima():
            for b in self.R.minima():
                yield a, b

    def root(self):
        """Return the root node of the FrameTree."""
        return (self.L.root(), self.R.root())

    def is_nonroot(self, frame):
        """Return True if the input frame has a parent."""
        a, b = frame
        return self.L.is_nonroot(a) or self.R.is_nonroot(b)

    def mode(self, frame):
        """Return the input frame's top (a leaf node)."""
        a, b = frame
        return (self.L.mode(a), self.R.mode(b))

    def has_children(self, frame):
        """Return True if the input frame has subframes."""
        a, b = frame
        return self.L.has_children(a) or self.R.has_children(b)

    def height(self, frame):
        """Return that height is NotImplemented."""
        return NotImplemented

    def base_height(self, frame):
        """Return that base_height is NotImplemented."""
        return NotImplemented

    def size(self, frame):
        """Return the input frame's size."""
        a, b = frame
        return max(self.L.size(a), self.R.size(b))

    def parent(self, frame):
        """Return the input frame merged with its neighboring frame."""
        a, b = frame
        if self.L.is_nonroot(a) and self.R.is_nonroot(b):
            pa, pb = self.L.parent(a), self.R.parent(b)
            if self.L.size(pa) > self.R.size(pb):
                return (a, pb)
            else:
                return (pa, b)
        elif not self.L.is_nonroot(a) and not self.R.is_nonroot(b):
            return None
        elif self.L.is_nonroot(a):
            # then b is root
            return (self.L.parent(a), b)
        else:
            # then a is root and b nonroot
            return (a, self.R.parent(b))

    def high(self, frame):
        """Return the input frame's high subframe."""
        a, b = frame
        if self.L.size(a) > self.R.size(b):
            return (self.L.high(a), b)
        else:
            return (a, self.R.high(b))

    def low(self, frame):
        """Return the input frame's low subframe."""
        a, b = frame
        if self.L.size(a) > self.R.size(b):
            return (self.L.low(a), b)
        else:
            return (a, self.R.low(b))

    def full(self, frame):
        """Return the largest frame with same mode as the input frame."""
        climber = frame
        while self.is_nonroot(climber):
            nextstep = self.parent(climber)
            if self.mode(nextstep) == self.mode(climber):
                climber = nextstep
            else:
                break
        return climber

    def index(self, frame):
        """Return a tuple of (nested) indices for the input frame."""
        a, b = frame
        return self.L.index(a), self.R.index(b)

    def filter(self, *, maxsize, localroot=None):
        """Yield grid frames contained in the input frame."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        ra, rb = localroot
        for a in self.L.filter(maxsize=maxsize, localroot=ra):
            for b in self.R.filter(maxsize=maxsize, localroot=rb):
                yield a, b

    def _find_parent_and_root(self):
        """Return that it is NotImplemented."""
        return NotImplemented

    def _find_mode_high_low(self):
        """Return that it is NotImplemented."""
        return NotImplemented

    def _find_full(self):
        """Return that it is NotImplemented."""
        return NotImplemented
