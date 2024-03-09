# -*- coding: utf-8 -*-
# This file is part of Peakoscope.
# Copyright (C) 2021-2024  Eivind Tøstesen
# Peakoscope is licensed under GPLv3.
"""Python module for generating data sets.

Synthetic data sets based on random walks:

The function randomwalk() takes as argument
a list of numeric steps and accumulates them
into a walk. Random steps are returned by the
functions discrete_steps() and continuous_steps().

The function alternating_steps() zips together
two lists of steps. Can be used to generate
zig-zagging random walks.

Functions example_1() and example_2() return data sets
intended as mini examples for peakoscope analysis.

Usage examples:
---------------

Generate a random walk with equally probable up and down steps:

>>> randomwalk(start=5, steps=discrete_steps(length=20, moves=[1, -1]))
[5, 4, 5, 4, 5, 4, 3, 2, 3, 4, 3, 4, 3, 2, 1, 2, 1, 2, 1, 0, 1]

A sequence with probabilities 0.85 and 0.15 of adding 1 and 10, respectively:

>>> randomwalk(steps=discrete_steps(length=20, moves=[1, 10], weights=[85, 15]))
[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 19, 20, 30, 31, 32, 33, 43, 44, 45, 55, 56]

A continuous random walk using the default parameter settings:

>>> randomwalk(steps=continuous_steps())
[0, -0.6659222888958249, 0.2856116398771149, -0.141155356156913, 0.17382464266411013]

"""

import random
import string
from itertools import accumulate, chain
from peakoscope.errors import PeakyBlunder


# Random walks:


def discrete_steps(
    *,
    length=10,
    moves=[1, -1],
    weights=None,
    randomseed="It's...",
):
    """Return list of discrete random steps."""
    random.seed(a=randomseed)
    return random.choices(population=moves, weights=weights, k=length)


def continuous_steps(
    *,
    length=4,
    moves=[1, -1],
    weights=None,
    randomseed="It's...",
):
    """Return list of continuous random steps."""
    if weights is None:
        # weights can also be set to random.triangular
        weights = random.uniform
    random.seed(a=randomseed)
    return [weights(*moves) for _ in range(length)]


def alternating_steps(
    steps1=discrete_steps(length=5, moves=[-1, -2, -3]),
    steps2=discrete_steps(length=5, moves=[1, 2, 3]),
):
    """Return alternating list of steps."""
    return list(chain.from_iterable(zip(steps1, steps2)))


def randomwalk(
    *,
    start=0,
    steps=discrete_steps(length=100, moves=[2, 1, 0, -1, -2]),
):
    """Return random walk as a list."""
    return list(accumulate(steps, initial=start))


def example_1(length=52, *, randomseed="it's....."):
    """Return mini example data set as two lists X, Y."""
    if length > 52:
        raise PeakyBlunder(f"Maximum length is 52.")
    return list(string.ascii_letters)[0:length], randomwalk(
        start=5.0,
        steps=discrete_steps(
            length=length - 1,
            moves=[0.2, 0.1, 0, -0.1, -0.2],
            randomseed=randomseed,
        ),
    )


def example_2(length=201, *, randomseed="It's..."):
    """Return mini example data set as two lists X, Y."""
    return list(range(1900, 1900 + length)), randomwalk(
        start=0.0,
        steps=discrete_steps(
            length=length - 1,
            moves=[-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0],
            randomseed=randomseed,
        ),
    )
