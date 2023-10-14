#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Python module for peaks and valleys in a numeric sequence.

Created on Thu Mar 30 16:31:00 2023

@author: Eivind Tostesen

"""


from utilities import pairwise


# Functions:


def peaks(values):
    """Yield peak regions as (start, end, min, max) tuples."""
    regions = []
    for i, (y1, y2) in enumerate(pairwise(values)):
        if i == 0:  # first data point:
            regions.append([i, None, y1, y1])
        if y2 > y1:  # if uphill:
            for r in reversed(regions):
                if r[3] >= y2:
                    break
                r[3] = y2  # update max value
            regions.append([i + 1, None, y2, y2])
        elif y2 == y1:
            pass  # region already created
        elif y2 < y1:  # if downhill:
            while regions and y2 < regions[-1][2]:
                popped = regions.pop()
                popped[1] = i  # update end value
                yield tuple(popped)
            if not (regions and y2 == regions[-1][2]):
                regions.append([popped[0], None, y2, popped[3]])
    for r in reversed(regions):
        r[1] = i + 1  # use last i value
        yield tuple(r)


# Classes:


class SliceStr(str):

    sep = ":"  # slice notation

    def __getattr__(self, name):
        """Get attribute."""
        if name == "slice":
            # Convert to python's slice object:
            return slice(*map(int, self.split(self.sep)))
        elif name == "tuple":
            # Tuple of integers (start, stop):
            return tuple(map(int, self.split(self.sep)))
        elif name == "start":
            # The start as integer:
            return self.slice.start
        elif name == "stop":
            # The stop as integer:
            return self.slice.stop
        elif name == "end":
            # The index of the last item:
            return self.stop - 1
        else:
            raise AttributeError()

    def __dir__(self):
        """Return list of attribute names."""
        return ["slice", "tuple", "start", "stop", "end"]

    def __repr__(self) -> str:
        """Return string that can reconstruct the object."""
        return f'SliceStr("{self}")'

    def __contains__(self, index):
        """Return True if index is in the slice."""
        return self.start <= index < self.stop

    def __len__(self):
        """Return number of items in the slice."""
        return self.stop - self.start

    def __iter__(self):
        """Yield indices in the slice."""
        yield from range(*self.tuple)

    def __le__(self, other) -> bool:
        """Return True if set(self) <= set(other)."""
        return other.start <= self.start and self.stop <= other.stop

    def __ge__(self, other) -> bool:
        """Return True if set(self) >= set(other)."""
        return other.start >= self.start and self.stop >= other.stop

    def __lt__(self, other) -> bool:
        """Return True if set(self) < set(other)."""
        return self <= other and self != other

    def __gt__(self, other) -> bool:
        """Return True if set(self) > set(other)."""
        return self >= other and self != other


class NumSlice(SliceStr):
    def __new__(cls, slicestr, values):
        """Create slice string with a reference to a numeric sequence."""
        self = super().__new__(cls, slicestr)
        self.values = values
        return self

    @classmethod
    def from_start_stop(cls, start, stop, values):
        self = super().__new__(cls, f"{start}:{stop}")
        self.values = values
        return self

    @classmethod
    def from_start_end(cls, start, end, values):
        self = super().__new__(cls, f"{start}:{end + 1}")
        self.values = values
        return self

    def __getattr__(self, name):
        """Get attribute."""
        if name == "max":
            # Maximum value in the slice of values:
            return max(self.values[self.slice])
        if name == "min":
            # Minimum value in the slice of values:
            return min(self.values[self.slice])
        if name == "argmax":
            # Index of (the first) maximum value in the slice:
            return self.values.index(self.max, *self.tuple)
        if name == "argmin":
            # Index of (the first) minimum value in the slice:
            return self.values.index(self.min, *self.tuple)
        if name == "size":
            # Distance between maximum and minimum value:
            return self.max - self.min
        else:
            # Attribute of SliceStr:
            return super().__getattr__(name)

    def __dir__(self):
        """Return list of attribute names."""
        return super().__dir__() + [
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
        """Return "start:stop" - string."""
        return str(self)

    def pre(self):
        """Return index of left neighbor if exists else None."""
        if self.start > 0:
            return self.start - 1
        else:
            return None

    def post(self):
        """Return index of right neighbor if exists else None."""
        if self.stop < len(self.values):
            return self.end + 1
        else:
            return None

    def is_peak(self):
        """Return True if slice of values is a peak."""
        return (self.pre() is None or self.values[self.pre()] < self.min) and (
            self.post() is None or self.values[self.post()] < self.min
        )

    def is_valley(self):
        """Return True if slice of values is a valley."""
        return (self.pre() is None or self.values[self.pre()] > self.max) and (
            self.post() is None or self.values[self.post()] > self.max
        )

    def is_local_maximum(self):
        """Return True if slice of values is a local maximum."""
        return self.size == 0 and self.is_peak()

    def is_local_minimum(self):
        """Return True if slice of values is a local minimum."""
        return self.size == 0 and self.is_valley()
