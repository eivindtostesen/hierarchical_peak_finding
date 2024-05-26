# Peakoscope
[![PyPI version](https://badge.fury.io/py/peakoscope.svg)](https://badge.fury.io/py/peakoscope)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Peakoscope is a python package for hierarchical analysis of peak and valley regions in numeric data.

![peak plot](https://github.com/eivindtostesen/hierarchical_peak_finding/blob/v1.0.0/output.png?raw=true "fig, (peaks.plot.ax, valleys.plot.ax) = plt.subplots(2, 1, sharex=True, figsize=(4, 4));
peaks.plot.crowns(peaks.size_filter(maxsize=7))
peaks.plot.bounding_boxes(peaks.size_filter(maxsize=7))
valleys.plot.crowns(valleys.size_filter(maxsize=7), facecolor='C9')
valleys.plot.bounding_boxes(valleys.size_filter(maxsize=7), edgecolor='C1')
peaks.plot.ax.set_title('Peak regions')
valleys.plot.ax.set_title('Valley regions')
peaks.plot.ax.plot(X, Y, linewidth=2, color='black')
valleys.plot.ax.plot(X, Y, linewidth=2, color='black');")

* Peak and valley regions can be nested, for example, when a large peak region contains smaller subpeak regions.
* Based on a one-pass algorithm that finds all peak regions and orders them into a tree.
* Classes for peak/valley objects and tree objects.
* Optional interfaces to matplotlib, pandas and polars.

## Usage examples
Compute the tree of nested peak regions in a data set:
```python
>>> import peakoscope
>>> data = [10, 30, 40, 30, 10, 50, 70, 70, 50, 80]
>>> print(peakoscope.tree(data))
0:10
├─5:10
│ ├─9:10
│ └─6:8
└─1:4
  └─2:3
```
From the tree, select default peak regions and print their subarrays of data:
```python
>>> for peak in peakoscope.tree(data).size_filter():
...    print(peak.subarray(data))
... 
[30, 40, 30]
[70, 70]
[80]
```

## Documentation
The github repo contains tutorials and a glossary:
* [plotting_tutorial.ipynb](https://nbviewer.org/github/eivindtostesen/hierarchical_peak_finding/blob/v1.0.0/plotting_tutorial.ipynb)
* [dataframes_tutorial.ipynb](https://nbviewer.org/github/eivindtostesen/hierarchical_peak_finding/blob/v1.0.0/dataframes_tutorial.ipynb)
* [glossary.rst](https://github.com/eivindtostesen/hierarchical_peak_finding/blob/v1.0.0/glossary.rst)

## Authors
* Eivind Tøstesen, <contact@tostesen.no>

## License
Copyright (C) 2021-2024 Eivind Tøstesen. This software is licensed under [GPLv3](https://github.com/eivindtostesen/hierarchical_peak_finding/blob/v1.0.0/LICENSE?raw=true "included LICENSE file")

## Citation
Citation can include one or more of:

* Peakoscope + version
* Github URL: https://github.com/eivindtostesen/hierarchical_peak_finding
* PyPI URL: https://pypi.org/project/peakoscope/
* The open-access article:

    >Tøstesen, E.
    >A stitch in time: Efficient computation of genomic DNA melting bubbles.
    >*Algorithms for Molecular Biology*, 3, 10 (2008).
    >[DOI: 10.1186/1748-7188-3-10](http://dx.doi.org/10.1186/1748-7188-3-10)
