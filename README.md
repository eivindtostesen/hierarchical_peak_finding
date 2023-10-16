# Hierarchical Peak Finding

Python tools for hierarchical analysis of peaks (and valleys) in numeric data. Based on a one-pass algorithm that finds all peak regions and orders them into a tree. Peak regions can be nested, for example, when a large peak region contains smaller subpeak regions. This is represented by classes for peak objects and tree objects, with optional interfaces to matplotlib, pandas or polars.

![peak plot](output.png "tree.plot.new(figsize=(10, 3)).ax.set_title('Bounding Boxes around peaks of maxsize=8')
tree.plot.crowns(tree.filter(maxsize=8), alpha=0.4)
tree.plot.bounding_boxes(tree.filter(maxsize=8), edgecolor='C4', linewidth=5)
tree.plot.ax.plot(X, Y, linewidth=2, color='C0');")

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
