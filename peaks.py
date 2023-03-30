#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Python module defining peaks as slices.

Created on Thu Mar 30 16:31:00 2023

@author: Eivind Tostesen

"""


# Classes:


class SliceStr(str):
    """A slice class with rich properties."""

    sep = ":"

    def __getattr__(self, name):
        if name == "slice":
            return slice(*map(int, self.split(self.sep)))
        elif name == "tuple":
            return tuple(map(int, self.split(self.sep)))
        elif name == "left":
            return self.slice.start - 1
        elif name == "right":
            return self.slice.stop
        elif name == "start":
            return self.slice.start
        elif name == "stop":
            return self.slice.stop
        elif name == "end":
            return self.stop - 1
        else:
            raise AttributeError()

    def __le__(self, other) -> bool:
        return other.start <= self.start and self.stop <= other.stop

    def __ge__(self, other) -> bool:
        return other.start >= self.start and self.stop >= other.stop

    def __lt__(self, other) -> bool:
        return self <= other and self != other

    def __gt__(self, other) -> bool:
        return self >= other and self != other

    def __contains__(self, i):
        return self.start <= i < self.stop

    def __iter__(self):
        yield from range(*self.tuple)

    def __len__(self):
        return self.stop - self.start

    def __repr__(self) -> str:
        return f'SliceStr("{self}")'

    def __dir__(self):
        return ("slice", "tuple", "start", "stop", "end", "right", "left")


class Peak(SliceStr):
    """A Peak class that extends a slice class."""

    def __new__(cls, slicestr, values):
        self = super().__new__(cls, slicestr)
        self.values = values
        return self

    def height(self):
        return max(self.values[self.slice])

    def base_height(self):
        return min(self.values[self.slice])

    def size(self):
        return self.height() - self.base_height()

    def argmax(self):
        return self.values[self.slice].index(self.height()) + self.start

    def argmin(self):
        return self.values[self.slice].index(self.base_height()) + self.start

    def __repr__(self) -> str:
        return f'Peak("{self}")'
