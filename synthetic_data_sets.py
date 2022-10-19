#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Python module for generating synthetic data.

Requires Python 3.8+

Created on Mon Oct 11 09:34:12 2021

@author: Eivind Tostesen
"""

import random
from itertools import accumulate, chain


# Random walks:


def discrete(length=10, moves=[1, -1], weights=None, randomseed="It's...",):
    """Return list of discrete random steps."""
    random.seed(a=randomseed)
    return random.choices(population=moves, weights=weights, k=length)


def continuous(length=10, moves=[1, -1], weights=None, randomseed="It's...",):
    """Return list of continuous random steps."""
    if weights is None:
         # weights can also be set to random.triangular
        weights = random.uniform
    random.seed(a=randomseed)
    return [weights(*moves) for _ in range(length)]


def alternating(steps1=discrete(length=5, moves=[-1, -2, -3]),
                steps2=discrete(length=5, moves=[1, 2, 3]),
               ):
    """Return alternating list of steps."""
    return list(chain.from_iterable(zip(steps1, steps2)))


def randomwalk(start=0, steps=discrete(length=100, moves=[2, 1, 0, -1, -2])):
    """Return random walk as a list."""
    return list(accumulate(steps, initial=start))


