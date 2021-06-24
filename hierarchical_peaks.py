#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Python module for analysis of peaks in numeric data.

This module contains algorithms for building a binary tree
representing the hierarchical structure of peaks in data.
It also contains methods for searching the tree and selecting peaks.

Created on Wed Mar 24 14:49:04 2021.

@author: Eivind Tostesen

"""


from itertools import chain


# Iteration tools:


def _forward_backward(iterable):
    """Return pair of iterators: (forward, backward)."""
    mylist = list(iterable)
    return iter(mylist), reversed(mylist)


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
    """Return dict of (start, end)-locations for the input peaks."""
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
    """Binary tree of hierarchical peaks in data.

    A PeakTree represents the hierarchical structure of
    all peaks in a univariate data set (such as a time series,
    a numeric sequence or a function x -> h).

    A PeakTree is initialized with an iterable of (x,h)-tuples,
    representing the local extrema in the data set
    (the alternating maximum and minimum points).

    The h-values must be numeric (elevation values).
    The x-values must be unique hashable objects for use as keys.

    Notes
    -----
    Recommended literature for the PeakTree class is
    the subsection titled "1D peaks" in the article [1]_.

    References
    ----------
    .. [1] Tostesen, E. "A stitch in time: Efficient computation of
       genomic DNA melting bubbles." Algorithms for Molecular
       Biology 3, 10 (2008).
       Open access: https://doi.org/10.1186/1748-7188-3-10
    """

    def __init__(self, maxima_and_minima):
        """Build a PeakTree and compute its data attributes."""
        self.elevation = dict(maxima_and_minima)
        # requires preserved insertion order:
        nodes = list(self.elevation)
        # remove possible minimum at the beginning:
        if self.elevation[nodes[0]] < self.elevation[nodes[1]]:
            del self.elevation[nodes[0]]
            del nodes[0]
        # remove possible minimum at the end:
        if self.elevation[nodes[-1]] < self.elevation[nodes[-2]]:
            del self.elevation[nodes[-1]]
            del nodes[-1]
        # compute data attributes:
        self._find_successor_and_root()
        self._find_mode_father_mother()
        self._find_full()

    def __contains__(self, peak):
        """Return True if the input is a node in the PeakTree."""
        return peak in self.elevation

    def __iter__(self):
        """Iterate over nodes in the PeakTree."""
        return iter(self.elevation)

    def __len__(self):
        """Return number of nodes in the PeakTree."""
        return len(self.elevation)

    def __mul__(self, other):
        """Return product with other tree (PeakTree or FrameTree)."""
        return FrameTree(self, other)

    def __str__(self):
        """Return printable tree structure."""
        return self.as_string(self.root())

    def as_string(self, localroot):
        """Return printable subtree structure."""
        str = "# Notation: <father> /& <mother>/ => <successor>\n"
        for full in chain([localroot],
                          self.foremothers(localroot)
                          ):
            if not self.has_parents(full):
                continue
            for node in self.path(self.mode(full),
                                  self.father(full),
                                  self.successor
                                  ):
                str += f'{node} /& {self.mother(self.successor(node))}/ => '
            str += f'{full}\n'
        return str

    def maxima(self):
        """Return ordered list of leaf nodes (local maxima)."""
        return list(self.elevation)[0::2]

    def minima(self):
        """Return ordered list of branch nodes (local minima)."""
        return list(self.elevation)[1::2]

    def root(self):
        """Return the root node of the PeakTree."""
        return self._root

    def is_nonroot(self, peak):
        """Return True if the input peak has a successor."""
        return peak != self._root

    def mode(self, peak):
        """Return the input peak's top (a leaf node)."""
        return self._mode[peak]

    def has_parents(self, peak):
        """Return True if the input peak has a father and mother."""
        return peak != self._mode[peak]

    def height(self, peak):
        """Return the input peak's top elevation."""
        return self.elevation[self.mode(peak)]

    def depth(self, peak):
        """Return the input peak's depth."""
        return self.height(peak) - self.elevation[peak]

    def successor(self, peak):
        """Return the input peak merged with its neighboring peak."""
        return self._successor[peak]

    def father(self, peak):
        """Return the highest of the input peak's two subpeaks."""
        return self._father[peak]

    def mother(self, peak):
        """Return the lowest of the input peak's two subpeaks."""
        return self._mother[peak]

    def full(self, peak):
        """Return the deepest peak with same mode as the input peak."""
        return self._full[self.mode(peak)]

    def index(self, peak):
        """Return the index of the input peak."""
        # accessing the preserved insertion order:
        return list(self.elevation).index(peak)

    # public recursive algorithms:

    def path(self, start, stop, step):
        """Yield nodes on a path in the tree."""
        climber = start
        yield climber
        while (climber != stop):
            climber = step(climber)
            yield climber

    def successors(self, node):
        """Yield nodes on the path from the input node to the root."""
        yield from self.path(node, self.root(), self.successor)

    def ancestors(self, localroot=None):
        """Yield all nodes in the input node's subtree."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        yield localroot
        yield from self.forefathers(localroot)
        yield from self.foremothers(localroot)

    def forefathers(self, localroot=None):
        """Yield fathers in the input node's subtree."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        if self.has_parents(localroot):
            yield self.father(localroot)
            yield from self.forefathers(self.father(localroot))
            yield from self.forefathers(self.mother(localroot))

    def foremothers(self, localroot=None):
        """Yield mothers in the input node's subtree."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        if self.has_parents(localroot):
            yield self.mother(localroot)
            yield from self.foremothers(self.mother(localroot))
            yield from self.foremothers(self.father(localroot))

    def paternal_line(self, node):
        """Yield nodes on the input node's paternal line."""
        yield from self.path(self.full(node), self.mode(node), self.father)

    def maxdeep(self, Dmax=3.0, localroot=None):
        """Yield nodes in subtree according to a maximum depth."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        # maxdeep algorithm:
        if self.depth(localroot) < Dmax:
            yield localroot
        else:
            climber = self.mode(localroot)
            while self.depth(self.successor(climber)) < Dmax:
                climber = self.successor(climber)
            yield climber
            while climber != localroot:
                climber = self.successor(climber)
                yield from self.maxdeep(Dmax, self.mother(climber))

    # Initialization algorithms:

    def _find_successor_and_root(self):
        """Compute attributes: self._successor and self._root."""
        self._successor = {}
        parents_in_spe = []
        parents = []
        for peak in self:
            while parents_in_spe:
                if self.elevation[parents_in_spe[-1]] < self.elevation[peak]:
                    break
                else:
                    parents.append(parents_in_spe.pop())
            if parents:
                parent = parents.pop()
                self._successor[parent] = peak
                while parents:
                    grandparent = parents.pop()
                    self._successor[grandparent] = parent
                    parent = grandparent
            parents_in_spe.append(peak)
        parent = parents_in_spe.pop()
        while parents_in_spe:
            child = parents_in_spe.pop()
            self._successor[parent] = child
            parent = child
        self._root = child

    def _find_mode_father_mother(self):
        """Compute attributes: self._mode, self._father and self._mother."""
        self._father = {}
        self._mother = {}

        def propagate_mode(peak):
            successor = self.successor(peak)
            if successor not in self._father:
                self._father[successor] = peak
            else:
                self._mother[successor] = peak
                hf = self.height(self.father(successor))
                hm = self.height(self.mother(successor))
                if hf < hm:
                    self._mother[successor] = self._father[successor]
                    self._father[successor] = peak
                elif hf == hm:
                    # mother goes left, father goes right:
                    if self.index(peak) > self.index(self.father(successor)):
                        self._mother[successor] = self._father[successor]
                        self._father[successor] = peak
                self._mode[successor] = self.mode(self.father(successor))
                if self.is_nonroot(successor):
                    propagate_mode(successor)

        self._mode = {peak: peak for peak in self.maxima()}
        for peak in self.maxima():
            propagate_mode(peak)

    def _find_full(self):
        """Compute attribute: self._full."""
        self._full = {self.mode(peak): peak for peak in self._mother.values()}
        self._full[self.mode(self._root)] = self._root


class FrameTree(PeakTree):
    """Binary tree of higher-dimensional hierarchical peaks.

    A FrameTree represents the hierarchical structure of
    peaks in more dimensions (such as a "mountain" landscape
    with height z as a function of x and y). It is based on
    an assumption of decoupling: The height function z(x,y)
    is a sum f(x) + g(y) or a product f(x) * g(y).

    A Frametree is constructed as a pair of lower-dimensional
    trees (of type PeakTree or Frametree). It is a kind of
    product tree, but it is not the Cartesian product.

    Subclass
    --------
    From PeakTree, the FrameTree class
    - overrides: __init__, __contains__, __iter__, 'is_nonroot',
      'mode', 'mother', 'root', 'successor',
    - copies: __str__, as_string, __mul__,
    - removes (NotImplemented): height, '_find_full',
      '_find_mode_father_mother', '_find_successor_and_root',

    Notes
    -----
    Recommended literature for the FrameTree class is
    the subsection titled "2D peaks" in the article [1]_.

    References
    ----------
    .. [1] Tostesen, E. "A stitch in time: Efficient computation of
       genomic DNA melting bubbles." Algorithms for Molecular
       Biology 3, 10 (2008).
       Open access: https://doi.org/10.1186/1748-7188-3-10
    """

    def __init__(self, x_tree, y_tree):
        self.x = x_tree
        self.y = y_tree

    def __contains__(self, frame):
        """Return True if the input is a node in the FrameTree."""
        a, b = frame
        # test if frame is sigma-above:
        return (
                (a == self.x.root()
                 or
                 self.x.depth(self.x.successor(a)) > self.y.depth(b)
                 )
                and
                (b == self.y.root()
                 or
                 self.y.depth(self.y.successor(b)) > self.x.depth(a)
                 )
                )

    def __iter__(self):
        """Iterate over nodes in the FrameTree."""
        yield from self.ancestors()

    def __len__(self):
        """Return number of nodes in the FrameTree."""
        return len(list(self.ancestors()))

    def maxima(self):
        """Yield leaf nodes (local maxima)."""
        for a in self.x.maxima():
            for b in self.y.maxima():
                yield a, b

    def minima(self):
        """Yield local minima."""
        for a in self.x.minima():
            for b in self.y.minima():
                yield a, b

    def root(self):
        """Return the root node of the FrameTree."""
        return (self.x.root(), self.y.root())

    def is_nonroot(self, frame):
        """Return True if the input frame has a successor."""
        a, b = frame
        return self.x.is_nonroot(a) or self.y.is_nonroot(b)

    def mode(self, frame):
        """Return the input frame's top (a leaf node)."""
        a, b = frame
        return (self.x.mode(a), self.y.mode(b))

    def has_parents(self, frame):
        """Return True if the input frame has a father and mother."""
        a, b = frame
        return self.x.has_parents(a) or self.y.has_parents(b)

    def height(self, frame):
        """Return that height is NotImplemented."""
        # todo ? sum of heights
        return NotImplemented

    def depth(self, frame):
        """Return the input frame's depth."""
        a, b = frame
        d1, d2 = self.x.depth(a), self.y.depth(b)
        return d1 if d1 > d2 else d2

    def successor(self, frame):
        """Return the input frame merged with its neighboring frame."""
        a, b = frame
        if self.x.is_nonroot(a) and self.y.is_nonroot(b):
            sa, sb = self.x.successor(a), self.y.successor(b)
            if self.x.depth(sa) > self.y.depth(sb):
                return (a, sb)
            else:
                return (sa, b)
        elif not self.x.is_nonroot(a) and not self.y.is_nonroot(b):
            return None
        elif self.x.is_nonroot(a):
            # then b is root
            return (self.x.successor(a), b)
        else:
            # then a is root and b nonroot
            return (a, self.y.successor(b))

    def father(self, frame):
        """Return the input frame's father subframe."""
        a, b = frame
        if self.x.depth(a) > self.y.depth(b):
            return (self.x.father(a), b)
        else:
            return (a, self.y.father(b))

    def mother(self, frame):
        """Return the input frame's mother subframe."""
        a, b = frame
        if self.x.depth(a) > self.y.depth(b):
            return (self.x.mother(a), b)
        else:
            return (a, self.y.mother(b))

    def full(self, frame):
        """Return the largest frame with same mode as the input frame."""
        climber = frame
        while self.is_nonroot(climber):
            nextstep = self.successor(climber)
            if self.mode(nextstep) == self.mode(climber):
                climber = nextstep
            else:
                break
        return climber

    def index(self, frame):
        """Return a tuple of (nested) indices for the input frame."""
        a, b = frame
        return self.x.index(a), self.y.index(b)

    def maxdeep(self, Dmax=3.0, localroot=None):
        """Yield grid frames contained in the input frame."""
        # defaults:
        if localroot is None:
            localroot = self.root()
        rx, ry = localroot
        for a in self.x.maxdeep(Dmax, rx):
            for b in self.y.maxdeep(Dmax, ry):
                yield a, b

    def _find_successor_and_root(self):
        """Return that _find_successor_and_root is NotImplemented."""
        return NotImplemented

    def _find_mode_father_mother(self):
        """Return that _find_mode_father_mother is NotImplemented."""
        return NotImplemented

    def _find_full(self):
        """Return that _find_full is NotImplemented."""
        return NotImplemented
