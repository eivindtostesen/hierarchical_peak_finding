# -*- coding: utf-8 -*-
# This file is part of Peakoscope.
# Copyright (C) 2021-2024  Eivind Tøstesen
# Peakoscope is licensed under GPLv3.
"""Python module containing tools and utilities.

"""


from peakoscope.errors import PeakyBlunder


# Iteration tools:


def pairwise(iterable):
    """Yield nearest neighbor pairs."""
    it = iter(iterable)
    a = next(it)
    for b in it:
        yield a, b
        a = b


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
