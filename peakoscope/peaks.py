# -*- coding: utf-8 -*-
# This file is part of Peakoscope.
# Copyright (C) 2021-2024  Eivind Tøstesen
# Peakoscope is licensed under GPLv3.
"""Python module for peak and valley regions in a numeric sequence.

"""


from collections import namedtuple
from operator import lt, gt, le, ge
from peakoscope.utilities import pairwise


# Functions:


def find_peaks(values, reverse=False):
    """Yield peak regions as (start, istop, argext, argcut, extremum, cutoff) Pose tuples."""
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
                yield Pose(*popped)
            if not (regions and y2 == regions[-1][5]):
                regions.append([popped[0], None, popped[2], i + 1, popped[4], y2])
    for r in reversed(regions):
        r[1] = i + 1  # use last i value
        yield Pose(*r)


def find_valleys(values):
    """Yield valley regions as (start, istop, argext, argcut, extremum, cutoff) Pose tuples."""
    yield from find_peaks(values, reverse=True)


# Classes:


Pose = namedtuple(
    "Pose",
    # A pose encodes the position and orientation of a peak or valley region:
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
    def from_start_stop(cls, start, stop, values):
        """Return Region from (start, stop) integers and the sequence."""
        self = super().__new__(cls, f"{start}:{stop}")
        self.values = values
        return self

    @classmethod
    def from_start_istop(cls, start, istop, values):
        """Return Region from (start, istop) integers and the sequence."""
        self = super().__new__(cls, f"{start}:{istop + 1}")
        self.values = values
        return self

    def __getattr__(self, name):
        """Get attribute."""
        if name == "slice":
            # Convert to python's slice object:
            return slice(*map(int, self.split(self.sep)))
        elif name == "tuple":
            # Tuple of integers (start, stop):
            return tuple(map(int, self.split(self.sep)))
        elif name == "start":
            # The start integer position:
            return self.slice.start
        elif name == "stop":
            # The (exclusive) stop integer position:
            return self.slice.stop
        elif name == "istop":
            # The index of the last item:
            return self.stop - 1
        if name == "max":
            # Maximum value in the region:
            return max(self.values[self.slice])
        if name == "min":
            # Minimum value in the region:
            return min(self.values[self.slice])
        if name == "argmax":
            # Index of (the first) maximum value in the region:
            return self.values.index(self.max, *self.tuple)
        if name == "argmin":
            # Index of (the first) minimum value in the region:
            return self.values.index(self.min, *self.tuple)
        if name == "size":
            # Distance between maximum and minimum value:
            return self.max - self.min
        else:
            raise AttributeError()

    def __dir__(self):
        """Return list of attribute/method names."""
        return [
            "slice",
            "tuple",
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
        yield from range(*self.tuple)

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


class Scope(Region):
    """String representing a peak region or valley region."""

    def __dir__(self):
        """Return list of attribute/method names."""
        return super().__dir__() + [
            "boundary_value",
        ]

    def __repr__(self) -> str:
        """Return string that can reconstruct the object."""
        return f'Scope("{self}", _)'

    def boundary_value(self):
        """Return peak's min or valley's max or None."""
        if self.is_peak():
            return self.min
        elif self.is_valley():
            return self.max
        else:
            return None
