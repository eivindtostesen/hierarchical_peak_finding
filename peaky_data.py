#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Python module for generating data sets.

Requires Python 3.8+

Created on Mon Oct 11 09:34:12 2021

@author: Eivind Tostesen
"""

import random
import string
from itertools import accumulate, chain
from errors import PeakyBlunder


# Random walks:


def discrete_steps(
    length=10,
    moves=[1, -1],
    weights=None,
    randomseed="It's...",
):
    """Return list of discrete random steps."""
    random.seed(a=randomseed)
    return random.choices(population=moves, weights=weights, k=length)


def continuous_steps(
    length=10,
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
