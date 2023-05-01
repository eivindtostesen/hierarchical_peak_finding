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


def add_L_arrow(
    axes, tail_x, tail_y, head_x, head_y, *, color="C2", linewidth=1, **kwargs
):
    """Plot an L-shaped arrow to indicate branch in PeakTree."""
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
    ax, x1, x2, y1, y2, *, edgecolor="C2", fill=False, linewidth=3, **kwargs
):
    """Plot bounding box around a peak."""
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


def add_pedestal(
    ax,
    x1,
    x2,
    y2,
    *,
    y1=0,
    fill=True,
    linewidth=1,
    edgecolor="C2",
    facecolor="gold",
    alpha=0.6,
    **kwargs
):
    """Plot a pedestal (below a bounding box)."""
    ax.add_patch(
        matplotlib.patches.Rectangle(
            xy=(x1, y1),
            width=x2 - x1,
            height=y2 - y1,
            fill=fill,
            linewidth=linewidth,
            edgecolor=edgecolor,
            facecolor=facecolor,
            alpha=alpha,
            **kwargs,
        )
    )


def add_crown(ax, xslice, yslice, y1, *, facecolor="gold", alpha=0.9, **kwargs):
    """Color the area under a peak."""
    ax.fill_between(
        xslice,
        yslice,
        y1,
        facecolor=facecolor,
        alpha=alpha,
        **kwargs,
    )


def add_bar(axes, x1, x2, y1, *, height=0.5, color="C3", fill=True, **kwargs):
    """Plot a bar to indicate peak location."""
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


class PeakTreeMatPlotLib(ChainedAttributes):
    """Plotting methods to be owned by a PeakTree."""

    def __init__(
        self,
        tree,
        attrname="plot",
        ax=None,
        fig=None,
        xy={},
        xinterval={},
        yinterval={},
        slices={},
    ):
        """Attach a plotting object to a PeakTree."""
        super().__init__()
        self.setattr(obj=tree, attrname=attrname)
        self.ax = ax
        self.fig = fig
        self.xy = xy
        self.xinterval = xinterval
        self.yinterval = yinterval
        self.slices = slices
        self.level = {
            n: len(list(self.rootself.root_path(n))) - 1 for n in self.rootself
        }
        if not xy:
            self.xy = {
                n: (self.rootself.top(n), self.rootself.base_height(n))
                for n in self.rootself
            }
        if not yinterval:
            self.yinterval = {
                n: (self.rootself.base_height(n), self.rootself.height(n))
                for n in self.rootself
            }
        plt.style.use("seaborn")

    def new(self):
        """Initialize new figure and axes."""
        self.fig = plt.figure(figsize=(10.0, 4.0))
        self.ax = self.fig.add_axes([0.1, 0.1, 1, 1])
        self.ax.set_xlim([min(self.rootself), max(self.rootself)])
        self.ax.set_ylim(
            [min(self.rootself._data.values()), max(self.rootself._data.values())]
        )
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
        """Plot bounding boxes around peaks."""
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
        """Fill color inside peaks above base lines."""
        if nodes is None:
            nodes = self.rootself
        if self.ax is None:
            self.new()
        for n in nodes:
            add_crown(
                self.ax,
                *self.slices[n],
                self.yinterval[n][0],
                **kwargs,
            )
        return self

    def pyramids(self, nodes=None, **kwargs):
        """Plot pyramids of stacked locations of peaks."""
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
                height=0.9,
                **kwargs,
            )
        return self
