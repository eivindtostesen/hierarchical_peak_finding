# -*- coding: utf-8 -*-
# This file is part of Peakoscope.
# Copyright (C) 2021-2025  Eivind Tøstesen
# Peakoscope is licensed under GPLv3.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Python module containing tools and utilities.

The class ChainedAttributes is a base class for writing classes
that can add functionality to objects at runtime,
which is then called via a new attribute.
(A kind of decorator design pattern.)

Usage examples:
---------------

Create a naive object:

>>> class Student:
...     pass
...
>>> student1 = Student()

Give it an external object as a new attribute:

>>> class Degree(ChainedAttributes):
...     pass
...
>>> degree1 = Degree().setattr(student1)
>>> student1.degree == degree1
True

Extend with a chained attribute:

>>> degree2 = Degree().setattr(student1.degree, attrname="higher_degree")
>>> student1.degree.higher_degree == degree2
True

The external objects access the original instance as 'rootself':

>>> student1.degree.rootself == student1
True
>>> student1.degree.higher_degree.rootself == student1
True

Make the object naive again:

>>> student1.degree.delattr()
>>> hasattr(student1, "degree")
False

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
                f"Attribute has already been set. Call method delattr() to break."
            )
        elif hasattr(obj, attrname):
            raise PeakyBlunder(
                f"The attribute name '{attrname=}' exists. Give a different attrname."
            )

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
