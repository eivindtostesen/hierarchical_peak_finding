# -*- coding: utf-8 -*-
# This file is part of Peakoscope.
# Copyright (C) 2021-2025  Eivind Tøstesen
# Peakoscope is licensed under GPLv3.


import pytest
import peakoscope


def test_version():
    assert peakoscope.__version__ == "1.1.0.dev8"
