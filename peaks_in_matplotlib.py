#!/usr/bin/env python3.
# -*- coding: utf-8 -*-
"""Python module for visualizing peaks in matplotlib.

Created on Thu Dec  1 16:13:14 2022

@author: Eivind Tostesen

"""

import matplotlib
import matplotlib.pyplot as plt
from utilities import ChainedAttributes


# Functions:


def add_L_arrow(axes, tail_x, tail_y, head_x, head_y, *,
                color="C2", linewidth=1, **kwargs):
    """Plot an L-shaped arrow to indicate branch in PeakTree."""
    axes.add_patch(
        matplotlib.patches.FancyArrowPatch(
            (tail_x, tail_y),
            (head_x, tail_y),
            arrowstyle='-',
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
            arrowstyle='->,head_length=2, head_width=1.5',
            shrinkA=0,
            shrinkB=1,
            linewidth=linewidth + 0.5,
            color=color,
            **kwargs,
        )
    )


def add_bounding_box(ax, x1, x2, base_height, size, *,
                     edgecolor="C2", fill=False, linewidth=3, **kwargs):
    """Plot bounding box around a peak."""
    ax.add_patch(
        matplotlib.patches.Rectangle(
            xy=(x1, base_height),
            width=x2 - x1,
            height=size,
            fill=fill,
            linewidth=linewidth,
            edgecolor=edgecolor,
            **kwargs,
        )
    )


def add_pedestal(ax, x1, x2, base_height, *,
                 fill=True, linewidth=1, edgecolor="C2",
                 facecolor="gold", alpha=0.6, **kwargs):
    """Plot a pedestal (below a bounding box)."""
    ax.add_patch(
        matplotlib.patches.Rectangle(
            xy=(x1, 0),
            width=x2 - x1,
            height=base_height,
            fill=fill,
            linewidth=linewidth,
            edgecolor=edgecolor,
            facecolor=facecolor,
            alpha=alpha,
            **kwargs,
        )
    )


def add_crown(ax, xslice, yslice, base_height, *,
              facecolor='gold', alpha=0.9, **kwargs):
    """Color the area under a peak."""
    ax.fill_between(
        xslice, yslice, base_height,
        facecolor=facecolor,
        alpha=alpha,
        **kwargs,
    )


def add_bar(axes, x1, x2, y, *,
            height=0.5, color="C3", fill=True, **kwargs):
    """Plot a bar to indicate peak location."""
    axes.add_patch(
        matplotlib.patches.Rectangle(
            xy=(x1, y),
            width=x2 - x1,
            height=height,
            color=color,
            fill=fill,
            **kwargs,
        )
    )


# Classes:


class PeakTreePlotter(ChainedAttributes):
    """Plotting methods to be owned by a PeakTree."""

    def __init__(self, tree, attrname='plot',
                 ax=None, fig=None,
                 xy={}, location={}, base_height={}, size={}, slices={},
                 ):
        """Attach a plotting object to a PeakTree."""
        super().__init__()
        self.setattr(obj=tree, attrname=attrname)
        self.ax = None
        self.fig = None
        self.xy = {n: (self.rootself.mode(n), self.rootself.base_height(n))
                   for n in self.rootself}
        self.base_height = {n: self.rootself.base_height(n)
                            for n in self.rootself}
        self.size = {n: self.rootself.size(n) for n in self.rootself}
        if ax:
            self.ax = ax
        if fig:
            self.fig = fig
        if xy:
            self.xy = xy
        if location:
            self.location = location
        if base_height:
            self.base_height = base_height
        if size:
            self.size = size
        if slices:
            self.slices = slices
        plt.style.use('seaborn')

    def new(self):
        """Initialize new figure and axes."""
        self.fig = plt.figure(figsize=(10.0, 4.0))
        self.ax = self.fig.add_axes([.1, .1, 1, 1])
        self.ax.set_xlim([min(self.rootself), max(self.rootself)])
        self.ax.set_ylim([min(self.rootself.data.values()),
                          max(self.rootself.data.values())])
        self.ax.set_xlabel('Label')
        self.ax.set_ylabel('Value')
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
        """Plot bounding boxes around peaks."""
        if nodes is None:
            nodes = self.rootself
        if self.ax is None:
            self.new()
        for n in nodes:
            add_bounding_box(
                self.ax,
                *self.location[n],
                self.base_height[n],
                self.size[n],
                **kwargs,
            )
        return self

    def crowns(self, nodes=None, **kwargs):
        """Fill color in the area below peaks."""
        if nodes is None:
            nodes = self.rootself
        if self.ax is None:
            self.new()
        for n in nodes:
            add_crown(
                self.ax,
                *self.slices[n],
                self.base_height[n],
                **kwargs,
            )

    def bars(self, nodes=None, **kwargs):
        """Plot bars that indicate peak sizes and locations."""
        if nodes is None:
            nodes = self.rootself
        # height = [height]
        if self.ax is None:
            self.new()
            self.ax.set_ylabel('Size')
            self.ax.set_ylim([0, self.rootself.size(self.rootself.root())])
        for n in nodes:
            add_bar(
                self.ax,
                *self.location[n],
                self.size[n],
                **kwargs,
            )
        return self
