# -*- coding: utf-8 -*-
# This file is part of Peakoscope.
# Copyright (C) 2021-2025  Eivind Tøstesen
# Peakoscope is licensed under GPLv3.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Python module for testing.

Collection of assertion functions for use in testing.

Two of the functions assert equality between objects.
The other 19 functions assert expected tree properties.
Input to the each of these is either a Tree or a HyperTree object.

Usage examples:
---------------

>>> import peakoscope
>>> data = [10, 30, 40, 30, 10, 50, 70, 70, 50, 80]
>>> tree = peakoscope.tree(data)
>>> assert_iteration_produces_members(tree)
>>> assert_leafs_have_no_children_and_root_has_no_parent(tree)
>>> assert_parent_and_children_are_inverse_of_each_other(tree)
>>> assert_level_is_length_of_root_path(tree)
>>> assert_root_is_outermost_and_leafs_are_innermost(tree)

"""


from itertools import chain
from peakoscope.errors import PeakyBlunder


# equality assertions:


def assert_region_equal(self, other):
    """Assert equality between two Region objects."""
    if self.values != other.values:
        raise PeakyBlunder("Can not compare regions in different sequences.")
    assert other.start == self.start and self.stop == other.stop


def assert_scope_equal(self, other):
    """Assert equality between two Scope objects."""
    assert_region_equal(self, other)
    assert self.argext == other.argext and self.argcut == other.argcut


# Generic tree consistency assertions:


def assert_iteration_produces_members(tree):
    """Assert iterating over a tree gives nodes in the tree."""
    assert all(node in tree for node in tree)


def assert_leafs_have_no_children_and_root_has_no_parent(tree):
    """Assert leaf nodes have no children, root has no parent."""
    assert all(tree.has_children(x) is False for x in tree.leaf_nodes())
    assert all(tree.children(x) == () for x in tree.leaf_nodes())
    assert tree.is_nonroot(tree.root()) is False
    assert tree.parent(tree.root()) is None


def assert_parent_and_children_are_inverse_of_each_other(tree):
    """Assert nodes are parents of children and children of parents.

    See also Proposition 5.
    """
    assert all(x in tree.children(tree.parent(x)) for x in tree if tree.is_nonroot(x))
    assert all(x == tree.parent(y) for x in tree for y in tree.children(x))


def assert_level_is_length_of_root_path(tree):
    """Assert level of a node is distance to root."""
    assert all(
        level == len(list(tree.root_path(node))) - 1 for node, level in tree.levels()
    )


def assert_root_is_outermost_and_leafs_are_innermost(tree):
    """Assert root is outermost and leaf nodes are innermost."""
    assert {tree.root()} == set(tree.outermost(tree))
    assert set(tree.leaf_nodes()) == set(tree.innermost(tree))


# tree partitioning assertions:


def assert_tree_consists_of_children_and_root(tree):
    """Assert all nodes are children or root."""
    assert set(tree) == (
        set(chain.from_iterable(tree.children(x) for x in tree)) | {tree.root()}
    )


def assert_tree_consists_of_leafs_and_linears_and_branches(tree):
    """Assert all nodes are leaf, linear or branch nodes."""
    assert set(tree) == (
        set(tree.leaf_nodes()) | set(tree.linear_nodes()) | set(tree.branch_nodes())
    )


def assert_tree_consists_of_main_paths(tree):
    """Assert all nodes are on a main path."""
    assert set(tree) == (
        set(chain.from_iterable(tree.main_path(x) for x in tree.full_nodes()))
    )


def assert_tree_consists_of_full_nodes_and_main_descendants(tree):
    """Assert all nodes are full nodes or main descendants."""
    assert set(tree) == set(tree.full_nodes()) | set(tree.main_descendants())


def assert_full_nodes_consists_of_lateral_descendants_plus_root(tree):
    """Assert all full nodes are lateral descendants or root."""
    assert set(tree.full_nodes()) == set(tree.lateral_descendants()) | {tree.root()}


def assert_children_consist_of_main_child_plus_lateral_children(tree):
    """Assert children is list of lateral extended with main if exists."""
    assert all(
        list(tree.children(x)) == [tree.main_child(x)] + list(tree.lateral(x))
        for x in tree
        if tree.main_child(x) is not None
    )
    assert all(
        list(tree.children(x)) == list(tree.lateral(x))
        for x in tree
        if tree.main_child(x) is None
    )


# tip and full assertions:


def assert_being_tip_is_having_no_main_child(tree):
    """Assert that nodes are tip nodes iff they have no main child."""
    assert all((x == tree.tip(x)) == (tree.main_child(x) is None) for x in tree)


def assert_main_child_keeps_tip_and_lateral_changes_tip(tree):
    """Assert main/lateral child has same/different tip as parent."""
    assert all(tree.tip(x) == tree.tip(tree.parent(x)) for x in tree.main_descendants())
    assert all(
        tree.tip(x) != tree.tip(tree.parent(x)) for x in tree.lateral_descendants()
    )


def assert_main_path_shares_tip_and_full(tree):
    """Assert nodes on main path have same tip node and full node."""
    assert all(
        tree.full(x) == tree.full(y) and tree.tip(x) == tree.tip(y)
        for x in tree.full_nodes()
        for y in tree.main_path(x)
    )


def assert_root_is_full_and_leafs_are_tips(tree):
    """Assert leaf nodes are tip nodes and root is a full node."""
    assert all(tree.tip(x) == x for x in tree.leaf_nodes())
    assert tree.full(tree.root()) == tree.root()


# node size assertions:


def assert_parent_size_is_strictly_greater(tree):
    """Assert each node is smaller than its parent.

    See also Proposition 1.
    """
    assert all(
        tree.size(tree.parent(x)) > tree.size(x) for x in tree if tree.is_nonroot(x)
    )


def assert_if_local_extremum_then_leaf(tree):
    """Assert zero-size nodes are included in leaf nodes."""
    assert set(x for x in tree if tree.size(x) == 0) <= set(tree.leaf_nodes())


# size_filter() assertions:


def assert_size_filter_equals_definition(tree, fractions=None):
    """Assert size_filter produces same nodes as a set comprehension.

    See also Lemma 1.
    """
    if fractions is None:
        fractions = (-0.1, 0.0, 0.2, 0.6, 2.0)
    rootsize = tree.size(tree.root())
    for maxsize in (fraction * rootsize for fraction in fractions):
        assert set(tree.size_filter(maxsize=maxsize)) == set(
            x
            for x in tree
            if tree.size(x) < maxsize
            and (x == tree.root() or tree.size(tree.parent(x)) >= maxsize)
        )


def assert_size_filter_equals_outermost_of_below_maxsize(tree, fractions=None):
    """Assert size_filter produces the outermost of all nodes below maxsize."""
    if fractions is None:
        fractions = (-0.1, 0.0, 0.2, 0.6, 2.0)
    rootsize = tree.size(tree.root())
    for maxsize in (fraction * rootsize for fraction in fractions):
        assert set(tree.size_filter(maxsize=maxsize)) == (
            set(tree.outermost(x for x in tree if tree.size(x) < maxsize))
        )
