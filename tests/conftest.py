# -*- coding: utf-8 -*-
# This file is part of Peakoscope.
# Copyright (C) 2021-2024  Eivind Tøstesen
# Peakoscope is licensed under GPLv3.


import pytest
import peakoscope.data


# Data sets:


X1, Y1 = peakoscope.example_1()


X2, Y2 = peakoscope.example_2()


discrete_randomwalk = peakoscope.data.randomwalk(
    start=123,
    steps=peakoscope.data.discrete_steps(
        length=49, moves=[3, 1, 0, -1, -3], randomseed=None
    ),
)


continuous_randomwalk = peakoscope.data.randomwalk(
    start=123.0,
    steps=peakoscope.data.continuous_steps(length=49, moves=[3, -3], randomseed=None),
)


zigzag_randomwalk = peakoscope.data.randomwalk(
    start=123,
    steps=peakoscope.data.alternating_steps(
        steps1=peakoscope.data.discrete_steps(length=30, moves=[-1, -2, -3]),
        steps2=peakoscope.data.discrete_steps(length=30, moves=[1, 2, 3]),
    ),
)


# Pathological data sets:
empty = []
singleton = [123]
flat = [123, 123]
up = [123, 124]
down = [123, 122]
nosignal = [123.0, 123.0, 123.0, 123.0, 123.0, 123.0, 123.0, 123.0, 123.0]


# Test data:


@pytest.fixture(
    params=[
        discrete_randomwalk,
        continuous_randomwalk,
        zigzag_randomwalk,
        Y1,
        Y2,
        flat,
        up,
        down,
        nosignal,
    ]
)
def data0(request):
    data = request.param
    return data


@pytest.fixture(params=[zigzag_randomwalk])
def zigzag_data(request):
    data = request.param
    return data


@pytest.fixture(params=[discrete_randomwalk, continuous_randomwalk, zigzag_randomwalk])
def data1(request):
    data = request.param
    return data


@pytest.fixture(params=[discrete_randomwalk, continuous_randomwalk, zigzag_randomwalk])
def data2(request):
    data = request.param
    return data
