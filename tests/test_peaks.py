# -*- coding: utf-8 -*-
# This file is part of Peakoscope.
# Copyright (C) 2021-2024  Eivind Tøstesen
# Peakoscope is licensed under GPLv3.


import pytest
from itertools import chain
import peakoscope
import peakoscope.data
from peakoscope import find_peaks, find_valleys, Region, Scope, Scope6


def _scope6(t):
    return Scope6(*t)


def _scope(t, values):
    return Scope.from_attrs(Scope6(*t), values)


def _region(t, values):
    return Region.from_attrs(Scope6(*t), values)


def _scope_equality(self, other):
    """For Scopes, the == is not enough, because Scope inherits Region.__eq__."""
    return self == other and self.argext == other.argext and self.argcut == other.argcut


@pytest.fixture(
    params=[
        # Discrete random walk data set of length 100:
        (
            123,
            peakoscope.data.discrete_steps(
                length=99, moves=[3, 1, 0, -1, -3], randomseed=None
            ),
        ),
        # Continuous random walk data set of length 100:
        (
            123.0,
            peakoscope.data.continuous_steps(length=99, moves=[3, -3], randomseed=None),
        ),
        # Pathological data sets of lengths 2 and 5:
        (123, [0]),
        (123, [1]),
        (123, [-1]),
        (123.0, [0.0, 0.0, 0.0, 0.0]),
    ],
)
def data(request):
    start, steps = request.param
    return peakoscope.data.randomwalk(start=start, steps=steps)


def test_version():
    assert peakoscope.__version__ == "1.1.0.dev2"


def test_mirror_symmetry(data):
    """Test that negating data swaps peaks and valleys."""
    assert all(
        v == (p[0], p[1], p[2], p[3], -p[4], -p[5])
        for v, p in zip(find_valleys(data), find_peaks(-y for y in data))
    )
    assert all(
        v == (p[0], p[1], p[2], p[3], -p[4], -p[5])
        for p, v in zip(find_peaks(data), find_valleys(-y for y in data))
    )


def test_corner_case(data):
    """Test that whole sequence is both peak and valley."""
    whole = Region(f"0:{len(data)}", data)
    assert whole.is_peak() and whole.is_valley()


def test_Region(data):
    """Test methods and attributes of Region objects.

    Where: all possible regions
    Objects: Region
    Attributes: size, argmax, argmin
    Methods: from_attrs, subarray, is_peak, is_valley
    Special methods: new, len, iter, contains
    """
    _peak_regions = {_region(t, data) for t in find_peaks(data)}
    _valley_regions = {_region(t, data) for t in find_valleys(data)}
    for start in range(0, len(data)):
        for istop in range(start, len(data)):
            region = Region(f"{start}:{istop + 1}", data)
            assert region.size == 0 if len(region) == 1 else True
            assert region.argmax in region
            assert region.argmin in region
            assert region.subarray() == [data[i] for i in region]
            assert region.is_peak() == (region in _peak_regions)
            assert region.is_valley() == (region in _valley_regions)
            assert all(i in region for i in region)


def test_local_extrema(data):
    """Test local maxima and minima.

    Where: all possible regions
    Objects: Region
    Methods: from_attrs, is_local_maximum, is_local_minimum
    Special methods: new
    """
    assert {
        t[0:2] for t in find_peaks(data) if _region(t, data).is_local_maximum()
    } == {
        (start, istop)
        for start in range(0, len(data))
        for istop in range(start, len(data))
        if Region(f"{start}:{istop + 1}", data).is_local_maximum()
    }
    assert {
        t[0:2] for t in find_valleys(data) if _region(t, data).is_local_minimum()
    } == {
        (start, istop)
        for start in range(0, len(data))
        for istop in range(start, len(data))
        if Region(f"{start}:{istop + 1}", data).is_local_minimum()
    }


def test_eval_repr(data):
    """Test that repr is readable by eval.

    Where: peak regions
    Objects: Scope6, Region, Scope
    Methods: from_attrs
    Special methods: repr
    """
    for t in find_peaks(data):
        _ = data
        assert _scope6(t) == eval(repr(_scope6(t))) == t
        assert _region(t, data) == eval(repr(_region(t, data)))
        assert _scope_equality(_scope(t, data), eval(repr(_scope(t, data))))


def test_horizontally(data):
    """Test attributes start, stop, istop.

    Where: peak and valley regions
    Objects: Scope6, Region, Scope
    Attributes: start, stop, istop
    Methods: from_attrs
    """
    assert all(
        _scope(t, data).start
        == _region(t, data).start
        == _scope6(t).start
        <= _scope(t, data).istop
        == _region(t, data).istop
        == _scope6(t).istop
        < _scope(t, data).stop
        == _region(t, data).stop
        == _scope6(t).istop + 1
        for t in chain(find_peaks(data), find_valleys(data))
    )


def test_vertically(data):
    """Test attributes max, min, argmax, argmin, extremum, cutoff, argext, argcut.

    Where: peak and valley regions
    Objects: Scope6, Region, Scope
    Attributes: max, min, argmax, argmin, extremum, cutoff, argext, argcut
    Methods: from_attrs
    """
    assert all(
        _scope(t, data).max
        == _region(t, data).max
        == _scope(t, data).extremum
        == _scope6(t).extremum
        and _scope(t, data).min
        == _region(t, data).min
        == _scope(t, data).cutoff
        == _scope6(t).cutoff
        and _scope(t, data).argmax == _region(t, data).argmax == _scope(t, data).argext
        and _scope(t, data).argmin == _region(t, data).argmin == _scope(t, data).argcut
        for t in find_peaks(data)
    )
    assert all(
        _scope(t, data).max
        == _region(t, data).max
        == _scope(t, data).cutoff
        == _scope6(t).cutoff
        and _scope(t, data).min
        == _region(t, data).min
        == _scope(t, data).extremum
        == _scope6(t).extremum
        and _scope(t, data).argmax == _region(t, data).argmax == _scope(t, data).argcut
        and _scope(t, data).argmin == _region(t, data).argmin == _scope(t, data).argext
        for t in find_valleys(data)
    )
