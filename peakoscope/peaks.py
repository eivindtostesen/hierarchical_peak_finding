# -*- coding: utf-8 -*-
# This file is part of Peakoscope.
# Copyright (C) 2021-2024  Eivind Tøstesen
# Peakoscope is licensed under GPLv3.
"""Python module for regions in a numeric sequence.

The function find_peaks goes through an iterable of numbers
in one pass yielding all peak regions as they are found.
The function find_valleys is the reverse mode of find_peaks.

A Region object represents any region in a sequence of numbers.
It is stored as start:stop positions and the sequence reference.

A Scope object represents a peak or valley region in a numeric sequence.
It is derived from the Region class, but in addition to start:stop,
it stores the positions argext and argcut.

'argext' is the (first) position of the 'extremum', defined as the
maximum in a peak region or the minimum in a valley region.

'argcut' is the (first) position of the 'cutoff', defined as the
minimum in a peak region or the maximum in a valley region.

Usage examples:
---------------

Find all peak regions in a data set:

>>> data = [10, 30, 40, 30, 10, 50, 70, 70, 50, 80]
>>> for p in find_peaks(data):
...     print(Scope6(*p))
... 
Scope6(start=2, istop=2, argext=2, argcut=2, extremum=40, cutoff=40)
Scope6(start=1, istop=3, argext=2, argcut=1, extremum=40, cutoff=30)
Scope6(start=6, istop=7, argext=6, argcut=6, extremum=70, cutoff=70)
Scope6(start=9, istop=9, argext=9, argcut=9, extremum=80, cutoff=80)
Scope6(start=5, istop=9, argext=9, argcut=5, extremum=80, cutoff=50)
Scope6(start=0, istop=9, argext=9, argcut=0, extremum=80, cutoff=10)


"""


from collections import namedtuple
from operator import lt, gt, le, ge
from peakoscope.utilities import pairwise


# Functions:


def find_peaks(values, reverse=False):
    """Yield peak regions as tuples: (start, istop, argext, argcut, extremum, cutoff)."""
    if not reverse:
        lessthan, greaterthan, greater_equal = lt, gt, ge
    else:
        lessthan, greaterthan, greater_equal = gt, lt, le
    regions = []
    for i, (y1, y2) in enumerate(pairwise(values)):
        if i == 0:  # first data point:
            regions.append([i, None, i, i, y1, y1])
        if greaterthan(y2, y1):  # if uphill:
            for r in reversed(regions):
                if greater_equal(r[4], y2):
                    break
                r[2], r[4] = i + 1, y2  # update argext and extremum value
            regions.append([i + 1, None, i + 1, i + 1, y2, y2])
        elif y2 == y1:
            pass  # region already created
        elif lessthan(y2, y1):  # if downhill:
            while regions and lessthan(y2, regions[-1][5]):
                popped = regions.pop()
                popped[1] = i  # update istop value
                yield tuple(popped)
            if not (regions and y2 == regions[-1][5]):
                regions.append([popped[0], None, popped[2], i + 1, popped[4], y2])
    for r in reversed(regions):
        r[1] = i + 1  # use last i value
        yield tuple(r)


def find_valleys(values):
    """Yield valley regions as tuples: (start, istop, argext, argcut, extremum, cutoff)."""
    yield from find_peaks(values, reverse=True)


# Classes:


Scope6 = namedtuple(
    "Scope6",
    # A Scope6 encodes the position and orientation of a peak or valley region:
    "start, istop, argext, argcut, extremum, cutoff",
    defaults=(None, None, None, None),
)


class Region(str):
    """String representing a region in a numeric sequence."""

    sep = ":"  # slice notation

    def __new__(cls, slicestr, values):
        """Create region with a reference to a numeric sequence."""
        self = super().__new__(cls, slicestr)
        self.values = values
        return self

    @classmethod
    def from_attrs(cls, obj, values):
        """Return Region from object's attributes and the sequence."""
        return Region(f"{obj.start}:{obj.istop + 1}", values)

    def __getattr__(self, name):
        """Get attribute."""
        if name == "start":
            # The start integer position:
            return int(self.split(self.sep)[0])
        elif name == "stop":
            # The (exclusive) stop integer position:
            return int(self.split(self.sep)[1])
        elif name == "istop":
            # The inclusive stop (index of the last item):
            return self.stop - 1
        elif name == "max":
            # Maximum value in the region:
            return max(self.subarray())
        elif name == "min":
            # Minimum value in the region:
            return min(self.subarray())
        elif name == "argmax":
            # Index of (the first) maximum value in the region:
            return self.values.index(self.max, self.start, self.stop)
        elif name == "argmin":
            # Index of (the first) minimum value in the region:
            return self.values.index(self.min, self.start, self.stop)
        elif name == "size":
            # Distance between maximum and minimum value:
            return self.max - self.min
        else:
            raise AttributeError()

    def __dir__(self):
        """Return list of attribute/method names."""
        return [
            "start",
            "stop",
            "istop",
            "values",
            "max",
            "min",
            "argmax",
            "argmin",
            "size",
            "pre",
            "post",
            "is_peak",
            "is_valley",
            "is_local_maximum",
            "is_local_minimum",
            "subarray",
        ]

    def __repr__(self) -> str:
        """Return string that can reconstruct the object."""
        return f'Region("{self}", _)'

    def __contains__(self, index):
        """Return True if index is in the region."""
        return self.start <= index < self.stop

    def __len__(self):
        """Return number of items in the region."""
        return self.stop - self.start

    def __iter__(self):
        """Yield indices in the region."""
        yield from range(self.start, self.stop)

    def __le__(self, other) -> bool:
        """Return True if set(self) <= set(other)."""
        return other.start <= self.start and self.stop <= other.stop

    def __ge__(self, other) -> bool:
        """Return True if set(self) >= set(other)."""
        return other.start >= self.start and self.stop >= other.stop

    def __lt__(self, other) -> bool:
        """Return True if set(self) < set(other)."""
        return self <= other and (other.start != self.start or self.stop != other.stop)

    def __gt__(self, other) -> bool:
        """Return True if set(self) > set(other)."""
        return self >= other and (other.start != self.start or self.stop != other.stop)

    def pre(self):
        """Return index of left neighbor if exists else None."""
        if self.start > 0:
            return self.start - 1
        else:
            return None

    def post(self):
        """Return index of right neighbor if exists else None."""
        if self.stop < len(self.values):
            return self.istop + 1
        else:
            return None

    def is_peak(self):
        """Return True if region is a peak."""
        return (self.pre() is None or self.values[self.pre()] < self.min) and (
            self.post() is None or self.values[self.post()] < self.min
        )

    def is_valley(self):
        """Return True if region is a valley."""
        return (self.pre() is None or self.values[self.pre()] > self.max) and (
            self.post() is None or self.values[self.post()] > self.max
        )

    def is_local_maximum(self):
        """Return True if region is a local maximum."""
        return self.size == 0 and self.is_peak()

    def is_local_minimum(self):
        """Return True if region is a local minimum."""
        return self.size == 0 and self.is_valley()

    def subarray(self, array=None):
        """Return region of (other) array as a new sequence."""
        if array is None:
            array = self.values
        return array[slice(self.start, self.stop)]


class Scope(Region):
    """String representing a peak region or valley region."""

    def __new__(cls, slicestr, argext, argcut, values):
        """Create region and add argext, argcut attributes."""
        self = super().__new__(cls, slicestr, values)
        self.argext = argext
        self.argcut = argcut
        return self

    @classmethod
    def from_attrs(cls, obj, values):
        """Return Scope from object's attributes and the sequence."""
        return Scope(f"{obj.start}:{obj.istop + 1}", obj.argext, obj.argcut, values)

    def __getattr__(self, name):
        """Get attribute."""
        if name == "extremum":
            # Return max of peak or min of valley:
            return self.values[self.argext]
        elif name == "cutoff":
            # Return min of peak or max of valley:
            return self.values[self.argcut]
        elif name == "max":
            # Maximum value in the region:
            return max(self.cutoff, self.extremum)
        elif name == "min":
            # Minimum value in the region:
            return min(self.cutoff, self.extremum)
        elif name == "argmax":
            # Index of (the first) maximum value in the region:
            return self.argcut if self.cutoff > self.extremum else self.argext
        elif name == "argmin":
            # Index of (the first) minimum value in the region:
            return self.argcut if self.cutoff < self.extremum else self.argext
        else:
            # Attribute of Region:
            return super().__getattr__(name)

    def __dir__(self):
        """Return list of attribute/method names."""
        return super().__dir__() + [
            "argext",
            "argcut",
            "extremum",
            "cutoff",
        ]

    def __repr__(self) -> str:
        """Return string that can reconstruct the object."""
        return f'Scope("{self}", {self.argext}, {self.argcut}, _)'
