#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Docstring for the hierarchical_peaks.py module.

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
    """Tree of hierarchical peaks in 1D data.

    Peak finding based on tree algorithms.

    Notes
    -----
    The subsection titled "1D peaks" in the article [1]_
    is relevant literature for the PeakTree class.

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
        nodes = sorted(self.elevation)
        # remove minimum at the beginning:
        if self.elevation[nodes[0]] < self.elevation[nodes[1]]:
            del self.elevation[nodes[0]]
            del nodes[0]
        # remove minimum at the end:
        if self.elevation[nodes[-1]] < self.elevation[nodes[-2]]:
            del self.elevation[nodes[-1]]
            del nodes[-1]
        # compute data attributes:
        self._find_successor_and_root(nodes)
        self._find_mode_father_mother(nodes)
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

    def __str__(self):
        """Return printable tree structure."""
        str = "# Notation: <father> /& <mother>/ => <successor>\n"
        for full in chain([self.root()],
                          self.foremothers(self.root())
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

    def _find_successor_and_root(self, nodes):
        """Compute attributes: self._successor and self._root."""
        self._successor = {}
        parents_in_spe = []
        parents = []
        for peak in nodes:
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

    def _find_mode_father_mother(self, nodes):
        """Compute attributes: self._mode, self._father and self._mother."""
        self._father = {}
        self._mother = {}
        maxima = nodes[0::2]
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
            if hf < hm or (hf == hm and peak > self.father(successor)):
                self._mother[successor] = self._father[successor]
                self._father[successor] = peak
            self._mode[successor] = self.mode(self.father(successor))
            if self.is_nonroot(successor):
                self._propagate_mode(successor)

    def _find_full(self):
        """Compute attribute: self._full."""
        self._full = {self.mode(peak): peak for peak in self._mother.values()}
        self._full[self.mode(self._root)] = self._root
