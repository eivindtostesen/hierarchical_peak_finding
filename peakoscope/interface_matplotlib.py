# -*- coding: utf-8 -*-
# This file is part of Peakoscope.
# Copyright (C) 2021-2024  Eivind Tøstesen
# Peakoscope is licensed under GPLv3.
"""Python module for visualizing peaks and valleys.

Matplotlib is used in the object-oriented way to plot peaks and valleys.

For each peak or valley, the functions add_L_arrow, add_bounding_box,
add_crown and add_bar can be called to add elements to an existing plot.

For a Tree of peaks or valleys, the class TreeMatPlotLib adds the methods
arrows, bounding_boxes, crowns and pyramids for plotting nodes in the Tree.

TreeMatPlotLib stores dicts containing x,y coordinates, with either default
or customized values, as an intermediary in the data flow from a Tree
to matplotlib, which provides more adaptability.

"""

import matplotlib
import matplotlib.pyplot as plt
from peakoscope.utilities import ChainedAttributes


# Functions:


def add_L_arrow(
    axes, tail_x, tail_y, head_x, head_y, *, color="C5", linewidth=1, **kwargs
):
    """Plot an L-shaped arrow to indicate branch in Tree."""
    axes.add_patch(
        matplotlib.patches.FancyArrowPatch(
            (tail_x, tail_y),
            (head_x, tail_y),
            arrowstyle="-",
            shrinkA=2,
            shrinkB=0,
            linewidth=linewidth,
            color=color,
            **kwargs,
        )
    )
    axes.add_patch(
        matplotlib.patches.FancyArrowPatch(
            (head_x, tail_y),
            (head_x, head_y),
            arrowstyle="->,head_length=2, head_width=1.5",
            shrinkA=0,
            shrinkB=1,
            linewidth=linewidth + 0.5,
            color=color,
            **kwargs,
        )
    )


def add_bounding_box(
    ax, x1, x2, y1, y2, *, edgecolor="C4", fill=False, linewidth=3, **kwargs
):
    """Plot bounding box around a peak or valley."""
    ax.add_patch(
        matplotlib.patches.Rectangle(
            xy=(x1, y1),
            width=x2 - x1,
            height=y2 - y1,
            fill=fill,
            linewidth=linewidth,
            edgecolor=edgecolor,
            **kwargs,
        )
    )


def add_crown(ax, xslice, yslice, y1, *, facecolor="gold", alpha=0.9, **kwargs):
    """Color the area of a peak or valley."""
    ax.fill_between(
        xslice,
        yslice,
        y1,
        facecolor=facecolor,
        alpha=alpha,
        **kwargs,
    )


def add_bar(axes, x1, x2, y1, *, height=0.5, color="C7", fill=True, **kwargs):
    """Plot a bar to indicate peak or valley location."""
    axes.add_patch(
        matplotlib.patches.Rectangle(
            xy=(x1, y1),
            width=x2 - x1,
            height=height,
            color=color,
            fill=fill,
            **kwargs,
        )
    )


# Classes:


class TreeMatPlotLib(ChainedAttributes):
    """Plotting methods to be owned by a Tree."""

    def __init__(
        self,
        tree,
        attrname="plot",
        ax=None,
        fig=None,
        X=None,
        Y=None,
        xlim=None,
        ylim=None,
        xy={},
        xinterval={},
        yinterval={},
        cutoff={},
        slice_of_x={},
        slice_of_y={},
    ):
        """Attach a plotting object to a Tree."""
        super().__init__()
        self.setattr(obj=tree, attrname=attrname)
        self.ax = ax
        self.fig = fig
        self.X = X
        self.Y = Y
        self.xlim = xlim
        self.ylim = ylim
        self.xy = xy
        self.xinterval = xinterval
        self.yinterval = yinterval
        self.cutoff = cutoff
        self.slice_of_x = slice_of_x
        self.slice_of_y = slice_of_y
        self.level = dict(self.rootself.levels())
        if not xlim:
            self.xlim = (
                self.X[self.rootself.root().start],
                self.X[self.rootself.root().istop],
            )
        if not ylim:
            self.ylim = (self.rootself.root().min, self.rootself.root().max)
        if not xy:
            self.xy = {n: (self.X[n.argext], n.cutoff) for n in self.rootself}
        if not xinterval:
            self.xinterval = {
                n: (self.X[n.start], self.X[n.istop]) for n in self.rootself
            }
        if not yinterval:
            self.yinterval = {n: (n.min, n.max) for n in self.rootself}
        if not cutoff:
            self.cutoff = {n: n.cutoff for n in self.rootself}
        plt.style.use("fast")

    def new(self, *, figsize=(10.0, 4.0)):
        """Initialize new figure and axes."""
        self.fig = plt.figure(figsize=figsize)
        self.ax = self.fig.add_axes([0.1, 0.1, 1, 1])
        self.ax.set_xlim(self.xlim)
        self.ax.set_ylim(self.ylim)
        self.ax.set_xlabel("Label")
        self.ax.set_ylabel("Value")
        return self

    def arrows(self, nodes=None, **kwargs):
        """Plot L arrows to nodes from their parent nodes."""
        if nodes is None:
            nodes = self.rootself
        if self.ax is None:
            self.new()
        for n in nodes:
            if self.rootself.is_nonroot(n):
                add_L_arrow(
                    self.ax,
                    *self.xy[self.rootself.parent(n)],
                    *self.xy[n],
                    **kwargs,
                )
        return self

    def bounding_boxes(self, nodes=None, **kwargs):
        """Plot bounding boxes around peaks or valleys."""
        if nodes is None:
            nodes = self.rootself
        if self.ax is None:
            self.new()
        for n in nodes:
            add_bounding_box(
                self.ax,
                *self.xinterval[n],
                *self.yinterval[n],
                **kwargs,
            )
        return self

    def crowns(self, nodes=None, **kwargs):
        """Fill color inside peaks or valleys."""
        if nodes is None:
            nodes = self.rootself
        if self.ax is None:
            self.new()
        for n in nodes:
            add_crown(
                self.ax,
                self.slice_of_x[n],
                self.slice_of_y[n],
                self.cutoff[n],
                **kwargs,
            )
        return self

    def pyramids(self, nodes=None, **kwargs):
        """Plot pyramids of stacked locations of peaks/valleys."""
        if nodes is None:
            nodes = self.rootself
        if self.ax is None:
            self.new()
            self.ax.set_ylabel("Level")
            self.ax.set_ylim([0, max(self.level.values()) + 1])
        for n in nodes:
            add_bar(
                self.ax,
                *self.xinterval[n],
                self.level[n],
                height=0.8,
                **kwargs,
            )
        return self
