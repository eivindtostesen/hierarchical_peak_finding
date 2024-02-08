# -*- coding: utf-8 -*-
# This file is part of Peakoscope.
# Copyright (C) 2021-2024  Eivind Tøstesen
# Peakoscope is licensed under GPLv3.
"""Python module containing tools and utilities.

"""


from peakoscope.errors import PeakyBlunder


# Iteration tools:


def forward_backward(iterable):
    """Return pair of iterators: (forward, backward)."""
    data = list(iterable)
    return iter(data), reversed(data)


def pairwise(iterable):
    """Yield nearest neighbor pairs."""
    it = iter(iterable)
    a = next(it)
    for b in it:
        yield a, b
        a = b


def tripletwise(iterable):
    """Yield nearest neighbor triplets."""
    it = iter(iterable)
    a, b = next(it), next(it)
    for c in it:
        yield a, b, c
        a, b = b, c


def filter_local_extrema(datapoints):
    """Let only maxima/minima pass from a stream of points."""
    previous = None
    for (x1, e1), (x2, e2) in pairwise(datapoints):
        if e1 == e2:
            continue
        is_uphill = e2 > e1
        if is_uphill == previous:
            continue
        yield x1, e1
        previous = is_uphill
    else:
        yield x2, e2


# Classes:


class ChainedAttributes:
    """Base class for adding chained attributes to naive objects."""

    def __init__(self):
        """Create initially empty object."""
        self.attrname = None
        self.parentself = None
        self.rootself = None

    def setattr(self, obj, attrname=None, rootself=None):
        """Attach self object to a parent object as an attribute."""
        if attrname is None:
            attrname = str.lower(type(self).__name__)

        if self.parentself is not None:
            raise PeakyBlunder(
                f"Attributed as '{self.attrname}'. Call method 'delattr' to break."
            )
        elif hasattr(obj, attrname):
            raise PeakyBlunder(f"The {attrname=} already exists.")

        self.attrname = attrname
        setattr(obj, self.attrname, self)
        self.parentself = obj

        if rootself is not None:
            self.rootself = rootself
        elif isinstance(obj, type(self)):
            self.rootself = obj.rootself
        else:
            self.rootself = obj

        return self

    def delattr(self):
        """Detach self object from parent object."""
        delattr(self.parentself, self.attrname)
        self.attrname = None
        self.parentself = None
        self.rootself = None
