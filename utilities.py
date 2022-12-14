#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 15:48:10 2022

@author: Eivind Tostesen

"""


# Classes:


class ChainedAttributes():
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

        if hasattr(obj, attrname) or self.parentself is not None:
            pass

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
