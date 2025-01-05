# -*- coding: utf-8 -*-
# This file is part of Peakoscope.
# Copyright (C) 2021-2025  Eivind Tøstesen
# Peakoscope is licensed under GPLv3.
"""Python module for testing.

Contains functions that assert expected tree properties.
Input to the each of these functions
can be either a Tree or HyperTree object.

"""


from itertools import chain


# Generic tree consistency assertions:


def assert_iteration_produces_members(tree):
    assert all(node in tree for node in tree)


def assert_leafs_have_no_children_and_root_has_no_parent(tree):
    assert all(tree.has_children(x) is False for x in tree.leaf_nodes())
    assert all(tree.children(x) == () for x in tree.leaf_nodes())
    assert tree.is_nonroot(tree.root()) is False
    assert tree.parent(tree.root()) is None


def assert_parent_and_children_are_inverse_of_each_other(tree):
    """Proposition 5."""
    assert all(x in tree.children(tree.parent(x)) for x in tree if tree.is_nonroot(x))
    assert all(x == tree.parent(y) for x in tree for y in tree.children(x))


def assert_level_is_length_of_root_path(tree):
    assert all(
        level == len(list(tree.root_path(node))) - 1 for node, level in tree.levels()
    )


def assert_root_is_outermost_and_leafs_are_innermost(tree):
    assert {tree.root()} == set(tree.outermost(tree))
    assert set(tree.leaf_nodes()) == set(tree.innermost(tree))


# tree partitioning assertions:


def assert_tree_consists_of_children_and_root(tree):
    assert set(tree) == (
        set(chain.from_iterable(tree.children(x) for x in tree)) | {tree.root()}
    )


def assert_tree_consists_of_leafs_and_linears_and_branches(tree):
    assert set(tree) == (
        set(tree.leaf_nodes()) | set(tree.linear_nodes()) | set(tree.branch_nodes())
    )


def assert_tree_consists_of_main_paths(tree):
    assert set(tree) == (
        set(chain.from_iterable(tree.main_path(x) for x in tree.full_nodes()))
    )


def assert_tree_consists_of_full_nodes_and_main_descendants(tree):
    assert set(tree) == set(tree.full_nodes()) | set(tree.main_descendants())


def assert_full_nodes_consists_of_lateral_descendants_plus_root(tree):
    assert set(tree.full_nodes()) == set(tree.lateral_descendants()) | {tree.root()}


def assert_children_consist_of_main_child_plus_lateral_children(tree):
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
    assert all((x == tree.tip(x)) == (tree.main_child(x) is None) for x in tree)


def assert_main_child_keeps_tip_and_lateral_changes_tip(tree):
    assert all(tree.tip(x) == tree.tip(tree.parent(x)) for x in tree.main_descendants())
    assert all(
        tree.tip(x) != tree.tip(tree.parent(x)) for x in tree.lateral_descendants()
    )


def assert_main_path_shares_tip_and_full(tree):
    assert all(
        tree.full(x) == tree.full(y) and tree.tip(x) == tree.tip(y)
        for x in tree.full_nodes()
        for y in tree.main_path(x)
    )


def assert_root_is_full_and_leafs_are_tips(tree):
    assert all(tree.tip(x) == x for x in tree.leaf_nodes())
    assert tree.full(tree.root()) == tree.root()


# node size assertions:


def assert_parent_size_is_strictly_greater(tree):
    """Proposition 1."""
    assert all(
        tree.size(tree.parent(x)) > tree.size(x) for x in tree if tree.is_nonroot(x)
    )


def assert_if_local_extremum_then_leaf(tree):
    assert set(x for x in tree if tree.size(x) == 0) <= set(tree.leaf_nodes())


# size_filter() assertions:


def assert_size_filter_equals_definition(tree, fractions=None):
    """Lemma 1."""
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
    if fractions is None:
        fractions = (-0.1, 0.0, 0.2, 0.6, 2.0)
    rootsize = tree.size(tree.root())
    for maxsize in (fraction * rootsize for fraction in fractions):
        assert set(tree.size_filter(maxsize=maxsize)) == (
            set(tree.outermost(x for x in tree if tree.size(x) < maxsize))
        )
