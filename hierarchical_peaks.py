#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Python module for analysis of peaks in numeric data.

This module contains algorithms for building a binary tree
representing the hierarchical structure of all peaks in data.
It also contains methods for searching the tree and selecting peaks.

Created on Wed Mar 24 14:49:04 2021.

@author: Eivind Tostesen

"""


from itertools import chain


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


def _pairwise(iterable):
    it = iter(iterable)
    a = next(it)
    for b in it:
        yield a, b
        a = b


def find_peak_locations(peaks):
    pass


def _surf():
    pass


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
        # make use of the preserved insertion order:
        nodes = list(self.elevation)
        # remove minimum at the beginning:
        if self.elevation[nodes[0]] < self.elevation[nodes[1]]:
            del self.elevation[nodes[0]]
            del nodes[0]
        # remove minimum at the end:
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

    def as_string(self, peak):
        """Return printable subtree structure."""
        str = "# Notation: <father> /& <mother>/ => <successor>\n"
        for full in chain([peak],
                          self.foremothers(peak)
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
        return self.elevation[self.mode(peak)] - self.elevation[peak]

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

    # public recursive algorithms:

    def path(self, start, stop, step):
        """Yield nodes on a path in the tree."""
        climber = start
        yield climber
        while (climber != stop):
            climber = step(climber)
            yield climber

    def successors(self, peak):
        """Yield nodes on the path from the input node to the root."""
        yield from self.path(peak, self.root(), self.successor)

    def ancestors(self, peak):
        """Yield all nodes in the input node's subtree."""
        yield peak
        yield from self.forefathers(peak)
        yield from self.foremothers(peak)

    def forefathers(self, peak):
        """Yield fathers in the input node's subtree."""
        if self.has_parents(peak):
            yield self.father(peak)
            yield from self.forefathers(self.father(peak))
            yield from self.forefathers(self.mother(peak))

    def foremothers(self, peak):
        """Yield mothers in the input node's subtree."""
        if self.has_parents(peak):
            yield self.mother(peak)
            yield from self.foremothers(self.mother(peak))
            yield from self.foremothers(self.father(peak))

    def paternal_line(self, peak):
        """Yield nodes on the input node's paternal line."""
        yield from self.path(self.full(peak), self.mode(peak), self.father)

    def maxdeep(self, peak, Dmax=3.0):
        """Yield nodes in subtree according to a maximum depth."""
        if self.depth(peak) < Dmax:
            yield peak
        else:
            climber = self.mode(peak)
            while self.depth(self.successor(climber)) < Dmax:
                climber = self.successor(climber)
            yield climber
            while climber != peak:
                climber = self.successor(climber)
                yield from self.maxdeep(self.mother(climber), Dmax)

    # Initialization algorithms:

    def _find_successor_and_root(self):
        """Compute attributes: self._successor and self._root."""
        self._successor = {}
        parents_in_spe = []
        parents = []
        for peak in self.elevation.keys():
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
        maxima = list(self.elevation)[0::2]
        self._mode = {peak: peak for peak in maxima}
        for peak in maxima:
            self._propagate_mode(peak)

    def _propagate_mode(self, peak):
        """Help _find_mode_father_mother()."""
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
                # to access the preserved insertion order:
                i = list(self.elevation).index(peak)
                j = list(self.elevation).index(self.father(successor))
                if i > j:
                    self._mother[successor] = self._father[successor]
                    self._father[successor] = peak
            self._mode[successor] = self.mode(self.father(successor))
            if self.is_nonroot(successor):
                self._propagate_mode(successor)

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
    - overwrites: __init__, __contains__, __iter__, 'is_nonroot',
      'mode', 'mother', 'root', 'successor',
    - copies: __str__, as_string, __mul__, 
    - removes (NotImplemented): height, '_find_full',
      '_find_mode_father_mother', '_find_successor_and_root',
      '_propagate_mode'

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
        yield from self.ancestors(self.root())

    def __len__(self):
        """Return number of nodes in the FrameTree."""
        return len(self.ancestors(self.root()))

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
        return NotImplemented

    def depth(self, frame):
        """Return the input frame's depth."""
        a, b = frame
        d1, d2 = self.x.depth(a), self.y.depth(b)
        return d1 if d1 > d2 else d2

    def successor(self, frame):
        """Return the input frame merged with its neighboring frame."""
        a, b = frame
        sa = self.x.successor(a) if self.x.is_nonroot(a) else None
        sb = self.y.successor(b) if self.y.is_nonroot(b) else None
        if (sa, sb) == (None, None):
            return None
        elif (sa is None or self.x.depth(sa) > self.y.depth(sb)):
            return (a, sb)
        elif (sb is None or self.x.depth(sa) < self.y.depth(sb)):
            return (sa, b)
        else:
            return (sa, b)

    def father(self, frame):
        a, b = frame
        if self.x.depth(a) > self.y.depth(b):
            return (self.x.father(a), b)
        else:
            return (a, self.y.father(b))

    def mother(self, frame):
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

    def grid(self, frame, Dmax):
        a, b = frame
        return (list(self.x.maxdeep(a, Dmax)), list(self.y.maxdeep(b, Dmax)))
