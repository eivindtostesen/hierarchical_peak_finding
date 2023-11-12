# Hierarchical Peak and Valley Finding

Python library for hierarchical analysis of peaks and valleys in numeric data.
* Peak and valley regions can be nested, for example, when a large peak region contains smaller subpeak regions.
* Based on a one-pass algorithm that finds all peak regions and orders them into a tree.
* Classes for peak/valley objects and tree objects, with optional interfaces to matplotlib, pandas or polars.

![peak plot](output.png "fig, (peaks.plot.ax, valleys.plot.ax) = plt.subplots(2, 1, sharex=True, figsize=(4, 4));
peaks.plot.crowns(peaks.filter(maxsize=7))
peaks.plot.bounding_boxes(peaks.filter(maxsize=7))
valleys.plot.crowns(valleys.filter(maxsize=7), facecolor="C9")
valleys.plot.bounding_boxes(valleys.filter(maxsize=7), edgecolor="C1")
peaks.plot.ax.set_title('Peak regions')
valleys.plot.ax.set_title('Valley regions')
peaks.plot.ax.plot(X, Y, linewidth=2, color="black")
valleys.plot.ax.plot(X, Y, linewidth=2, color="black");")

## Howto files
For tutorials and documentation of the usage, jupyter notebooks are (or will be) included:
* [plotting_tutorial.ipynb](plotting_tutorial.ipynb)
* [dataframes_tutorial.ipynb](dataframes_tutorial.ipynb)

## Citation
In publications, please cite the current homepage of this software,

* https://github.com/eivindtostesen/hierarchical_peak_finding

and this open-access article:

>Tøstesen, E.
>A stitch in time: Efficient computation of genomic DNA melting bubbles.
>*Algorithms for Molecular Biology*, 3, 10 (2008).
>[DOI: 10.1186/1748-7188-3-10](http://dx.doi.org/10.1186/1748-7188-3-10)


## Authors
* Eivind Tøstesen, <contact@tostesen.no>

## Copying and license
Copyright 2021 Eivind Tøstesen.

License: GPL v3
