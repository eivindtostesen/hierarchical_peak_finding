# -*- coding: utf-8 -*-
# This file is part of Peakoscope.
# Copyright (C) 2021-2025  Eivind Tøstesen
# Peakoscope is licensed under GPLv3.
# SPDX-License-Identifier: GPL-3.0-or-later


import pytest
from itertools import chain
import peakoscope
import peakoscope.data
import peakoscope.testing as testing
from peakoscope import find_peaks, find_valleys, Region, Scope, Scope6


def _scope6(t):
    """Convert tuple to Scope6."""
    return Scope6(*t)


def _scope(t, values):
    """Convert tuple to Scope."""
    return Scope.from_attrs(Scope6(*t), values)


def _region(t, values):
    """Convert tuple to Region."""
    return Region.from_attrs(Scope6(*t), values)


def _all_regions(values):
    """Yield all possible regions."""
    for start in range(0, len(values)):
        for istop in range(start, len(values)):
            yield Region(f"{start}:{istop + 1}", values)


def _peak_regions(values):
    """Yield all peak regions."""
    for t in find_peaks(values):
        yield _region(t, values)


def _valley_regions(values):
    """Yield all valley regions."""
    for t in find_valleys(values):
        yield _region(t, values)


def test_mirror_symmetry(data0):
    """Test that negating data swaps peaks and valleys."""
    assert all(
        v == (p[0], p[1], p[2], p[3], -p[4], -p[5])
        for v, p in zip(find_valleys(data0), find_peaks(-y for y in data0))
    )
    assert all(
        v == (p[0], p[1], p[2], p[3], -p[4], -p[5])
        for p, v in zip(find_peaks(data0), find_valleys(-y for y in data0))
    )


def test_corner_case(data0):
    """Test that whole sequence is both peak and valley."""
    whole = Region(f"0:{len(data0)}", data0)
    assert whole.is_peak() and whole.is_valley()


def test_Region(data0):
    """Test some methods and attributes of Region objects.

    Attributes: size, max, min, argmax, argmin
    Methods: subarray
    Special methods: new, len, iter, contains
    """
    for region in _all_regions(data0):
        assert region.size == 0 if len(region) == 1 else True
        assert region.argmax in region
        assert region.argmin in region
        assert region.max == region.values[region.argmax]
        assert region.min == region.values[region.argmin]
        assert region.subarray() == [data0[i] for i in region]
        assert all(i in region for i in region)


def test_all_peaks_and_valleys_are_found(data0):
    """Test that find_peaks and find_valleys do not miss anything."""
    assert set(_peak_regions(data0)) == {r for r in _all_regions(data0) if r.is_peak()}
    assert set(_valley_regions(data0)) == {
        r for r in _all_regions(data0) if r.is_valley()
    }


def test_all_local_extrema_are_found(data0):
    """Test that find_peaks and find_valleys include all local extrema."""
    assert {r for r in _all_regions(data0) if r.is_local_maximum()} == {
        r for r in _peak_regions(data0) if r.is_local_maximum()
    }
    assert {r for r in _all_regions(data0) if r.is_local_minimum()} == {
        r for r in _valley_regions(data0) if r.is_local_minimum()
    }


def test_eval_repr(data0):
    """Test that repr is readable by eval."""
    Region.default_data = data0
    Scope.default_data = data0
    for t in find_peaks(data0):
        assert _scope6(t) == eval(repr(_scope6(t))) == t
        testing.assert_region_equal(_region(t, data0), eval(repr(_region(t, data0))))
        testing.assert_scope_equal(_scope(t, data0), eval(repr(_scope(t, data0))))
    Region.default_data = None
    Scope.default_data = None


def test_horizontal_attributes(data0):
    """Test attributes start, stop, istop."""
    assert all(
        _scope(t, data0).start
        == _region(t, data0).start
        == _scope6(t).start
        <= _scope(t, data0).istop
        == _region(t, data0).istop
        == _scope6(t).istop
        < _scope(t, data0).stop
        == _region(t, data0).stop
        == _scope6(t).istop + 1
        for t in chain(find_peaks(data0), find_valleys(data0))
    )


def test_vertical_attributes(data0):
    """Test attributes max, min, argmax, argmin, extremum, cutoff, argext, argcut."""
    assert all(
        _scope(t, data0).max
        == _region(t, data0).max
        == _scope(t, data0).extremum
        == _scope6(t).extremum
        and _scope(t, data0).min
        == _region(t, data0).min
        == _scope(t, data0).cutoff
        == _scope6(t).cutoff
        and _scope(t, data0).argmax
        == _region(t, data0).argmax
        == _scope(t, data0).argext
        and _scope(t, data0).argmin
        == _region(t, data0).argmin
        == _scope(t, data0).argcut
        for t in find_peaks(data0)
    )
    assert all(
        _scope(t, data0).max
        == _region(t, data0).max
        == _scope(t, data0).cutoff
        == _scope6(t).cutoff
        and _scope(t, data0).min
        == _region(t, data0).min
        == _scope(t, data0).extremum
        == _scope6(t).extremum
        and _scope(t, data0).argmax
        == _region(t, data0).argmax
        == _scope(t, data0).argcut
        and _scope(t, data0).argmin
        == _region(t, data0).argmin
        == _scope(t, data0).argext
        for t in find_valleys(data0)
    )
