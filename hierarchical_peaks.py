#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 14:49:04 2021.

@author: Eivind Tostesen

"""

import itertools as it


def filter_local_extrema(sequence):
    previous = None
    for (x1, e1), (x2, e2) in _pairwise(sequence):
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


class PeakTree:

    def __init__(self, maxima_and_minima):
        self.elevation = dict(maxima_and_minima)
        nodes = sorted(self.elevation)
        if self.elevation[nodes[0]] < self.elevation[nodes[1]]:
            del self.elevation[nodes[0]]
            del nodes[0]
        if self.elevation[nodes[-1]] < self.elevation[nodes[-2]]:
            del self.elevation[nodes[-1]]
            del nodes[-1]
        self._find_successor_and_root(nodes)
        self._find_mode_father_mother(nodes)
        self._find_full()

    def __contains__(self, peak):
        return peak in self.elevation

    def __iter__(self):
        return iter(self.elevation)

    def __len__(self):
        return len(self.elevation)

    def __str__(self):
        str = "# Notation: <father> /& <mother>/ => <successor>\n"
        for full in it.chain([self.root()], self.foremothers(self.root())):
            if not self.has_parents(full):
                continue
            for peak in self.path(self.mode(full),
                                  self.father(full),
                                  self.successor):
                str += f'{peak} /& {self.mother(self.successor(peak))}/ => '
            str += f'{full}\n'

        return str

    def root(self):
        return self._root

    def is_nonroot(self, peak):
        return peak != self._root

    def mode(self, peak):
        return self._mode[peak]

    def has_parents(self, peak):
        return peak != self._mode[peak]

    def height(self, peak):
        return self.elevation[self.mode(peak)]

    def depth(self, peak):
        return self.elevation[self.mode(peak)] - self.elevation[peak]

    def successor(self, peak):
        return self._successor[peak]

    def father(self, peak):
        return self._father[peak]

    def mother(self, peak):
        return self._mother[peak]

    def full(self, peak):
        return self._full[self.mode(peak)]

    # public recursive algorithms:

    def maxdeep(self, peak, Dmax=3.0):
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

    def forefathers(self, peak):
        if self.has_parents(peak):
            yield self.father(peak)
            yield from self.forefathers(self.father(peak))
            yield from self.forefathers(self.mother(peak))

    def foremothers(self, peak):
        if self.has_parents(peak):
            yield self.mother(peak)
            yield from self.foremothers(self.mother(peak))
            yield from self.foremothers(self.father(peak))

    def ancestors(self, peak):
        yield peak
        yield from self.forefathers(peak)
        yield from self.foremothers(peak)

    def path(self, start, stop, step):
        climber = start
        yield climber
        while (climber != stop):
            climber = step(climber)
            yield climber

    def paternal_line(self, peak):
        yield from self.path(self.full(peak), self.mode(peak), self.father)

    # Initialization algorithms:

    def _find_successor_and_root(self, nodes):
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
        self._father = {}
        self._mother = {}
        maxima = nodes[0::2]
        self._mode = {peak: peak for peak in maxima}
        for peak in maxima:
            self._propagate_mode(peak)

    def _propagate_mode(self, peak):
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
        self._full = {self.mode(peak): peak for peak in self._mother.values()}
        self._full[self.mode(self._root)] = self._root
