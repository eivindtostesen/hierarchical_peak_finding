# -*- coding: utf-8 -*-
# This file is part of Peakoscope.
# Copyright (C) 2021-2024  Eivind Tøstesen
# Peakoscope is licensed under GPLv3.


import pytest
from itertools import chain
from operator import attrgetter
import peakoscope
import peakoscope.data
from peakoscope import Tree, HyperTree, Scope, Scope6, find_peaks


# test Tree objects (one-dimensional):


@pytest.fixture(params=[True, False])
def tree(data0, request):
    valleys = request.param
    return peakoscope.tree(data0, valleys=valleys)


def test_tree(tree):
    """Test that Tree objects have all known behaviors."""
    assert_basic_tree_consistency(tree)
    # for 1D only:
    assert_tree_nodes_are_nested_scopes(tree)
    assert_main_keeps_argext_and_lateral_changes_argext(tree)
    assert_dict_of_dicts_of_same_length(tree)
    # tree partitions:
    assert_tree_consists_of_full_nodes_and_main_descendants(tree)
    assert_full_nodes_consists_of_lateral_descendants_plus_root(tree)
    assert_children_consist_of_main_child_plus_lateral_children(tree)
    assert_tree_consists_of_main_paths(tree)
    # tip and full:
    assert_being_tip_is_having_no_main_child(tree)
    assert_main_child_keeps_tip_and_lateral_changes_tip(tree)
    assert_main_path_shares_tip_and_full(tree)
    assert_root_is_full_and_leafs_are_tips(tree)
    # node size:
    assert_parent_size_is_strictly_greater(tree)
    assert_if_local_extremum_then_leaf(tree)
    # size filter:
    assert_size_filter_equals_outermost_of_smaller_nodes(tree)
    assert_size_filter_equals_definition(tree)


# test HyperTree objects (pairs of 1D trees):


@pytest.fixture(params=[True, False])
def pair_of_trees(data1, data2, request):
    valleys = request.param
    return (
        peakoscope.tree(data1, valleys=valleys),
        peakoscope.tree(data2, valleys=valleys),
    )


def test_hypertree(pair_of_trees):
    """Test that HyperTree objects have all known behaviors."""
    t1, t2 = pair_of_trees
    assert_basic_tree_consistency(t1 @ t2)
    # for hypertree only:
    assert_parents_and_children_belong_to_tree(t1 @ t2)
    assert_grid_nodes_belong_to_tree(t1 @ t2)
    assert_leafs_belong_to_tree(t1 @ t2)
    assert_parent_size_equals_minimum(t1 @ t2)
    # tree partitions:
    assert_tree_consists_of_full_nodes_and_main_descendants(t1 @ t2)
    assert_full_nodes_consists_of_lateral_descendants_plus_root(t1 @ t2)
    assert_children_consist_of_main_child_plus_lateral_children(t1 @ t2)
    assert_tree_consists_of_main_paths(t1 @ t2)
    # tip and full:
    assert_being_tip_is_having_no_main_child(t1 @ t2)
    assert_main_child_keeps_tip_and_lateral_changes_tip(t1 @ t2)
    assert_main_path_shares_tip_and_full(t1 @ t2)
    assert_root_is_full_and_leafs_are_tips(t1 @ t2)
    # node size:
    assert_parent_size_is_strictly_greater(t1 @ t2)
    assert_if_local_extremum_then_leaf(t1 @ t2)
    # size filter:
    assert_size_filter_equals_outermost_of_smaller_nodes(t1 @ t2)
    assert_size_filter_equals_definition(t1 @ t2)
    assert_size_filter_recursion_equals_cartesian_product(t1 @ t2)


# special tests:


def test_non_linear_tree(zigzag_data):
    """Test that zigzag data gives tree without linear nodes."""
    assert len(list(peakoscope.tree(zigzag_data).linear_nodes())) == 0


def test_eval_repr(data1):
    """Test that repr is readable by eval."""
    _ = data1
    tree = peakoscope.tree(_)
    # Tree:
    assert repr(tree) == repr(eval(repr(tree)))
    # HyperTree:
    assert repr(tree @ tree) == repr(eval(repr(tree @ tree)))
    assert repr(tree @ tree @ tree) == repr(eval(repr(tree @ tree @ tree)))


def test_children_are_sorted_by_extremum(data1):
    peaks = peakoscope.tree(data1)
    valleys = peakoscope.tree(data1, valleys=True)
    assert all(
        list(peaks.children(x))
        == sorted(peaks.children(x), key=attrgetter("extremum"), reverse=True)
        for x in peaks
    )
    assert all(
        list(valleys.children(x))
        == sorted(valleys.children(x), key=attrgetter("extremum"), reverse=False)
        for x in valleys
    )


def test_changing_the_node_type(data1):
    # create tree with nodes of type Scope:
    scope_tree = peakoscope.tree(data1)
    # create tree with nodes of type Scope6:
    scope6_tree = Tree.from_peaks(map(lambda t: Scope6(*t), find_peaks(data1)))
    # lists of nodes are NOT equal:
    assert list(scope_tree) != list(scope6_tree)
    # change from Scope6 to Scope:
    scope6_tree.set_nodes({s6: Scope.from_attrs(s6, data1) for s6 in scope6_tree})
    # now lists of nodes are equal:
    assert list(scope_tree) == list(scope6_tree)


########### Consistency assertions: #################


# for 1D only:


def assert_tree_nodes_are_nested_scopes(tree):
    assert all(x < tree.parent(x) for x in tree if tree.is_nonroot(x))
    assert all(y < x for x in tree for y in tree.children(x))
    assert all(tree.tip(x) <= x <= tree.full(x) <= tree.root() for x in tree)


def assert_main_keeps_argext_and_lateral_changes_argext(tree):
    assert all(x.argext == tree.parent(x).argext for x in tree.main_descendants())
    assert all(x.argext != tree.parent(x).argext for x in tree.lateral_descendants())


def assert_dict_of_dicts_of_same_length(tree):
    assert all(
        len(tree) == len(d)
        for (_, d) in tree.as_dict_of_dicts().items()
        if type(d) == dict
    )


# for HyperTree only:


def assert_parent_size_equals_minimum(tree):
    """Proposition 2."""
    assert all(
        tree.size(tree.parent((a, b)))
        == min(tree.L.size(tree.L.parent(a)), tree.R.size(tree.R.parent(b)))
        for (a, b) in tree
        if tree.L.is_nonroot(a) and tree.R.is_nonroot(b)
    )


def assert_parents_and_children_belong_to_tree(tree):
    """Proposition 4."""
    # parent in tree:
    assert all(tree.parent(x) in tree for x in tree if tree.is_nonroot(x))
    # children in tree:
    assert all(x in tree for y in tree for x in tree.children(y))


def assert_grid_nodes_belong_to_tree(tree):
    """Proposition 7."""
    assert all(x in tree for x in tree.size_filter())


def assert_leafs_belong_to_tree(tree):
    assert all(x in tree for x in tree.leaf_nodes())


# tree partitioning assertions:


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
    rsize = tree.size(tree.root())
    for maxsize in (fraction * rsize for fraction in fractions):
        assert set(tree.size_filter(maxsize=maxsize)) == set(
            x
            for x in tree
            if tree.size(x) < maxsize
            and (x == tree.root() or tree.size(tree.parent(x)) >= maxsize)
        )


def assert_size_filter_equals_outermost_of_smaller_nodes(tree, fractions=None):
    if fractions is None:
        fractions = (-0.1, 0.0, 0.2, 0.6, 2.0)
    rsize = tree.size(tree.root())
    for maxsize in (fraction * rsize for fraction in fractions):
        assert set(tree.size_filter(maxsize=maxsize)) == (
            set(tree.outermost(x for x in tree if tree.size(x) < maxsize))
        )


def assert_size_filter_recursion_equals_cartesian_product(tree, fractions=None):
    if fractions is None:
        fractions = (0.0, 0.01, 0.2, 0.6, 2.0)
    rsize = tree.size(tree.root())
    base_class_method = Tree.size_filter
    overriding_method = HyperTree.size_filter
    for maxsize in (fraction * rsize for fraction in fractions):
        HyperTree.size_filter = base_class_method
        by_recursion = set(tree.size_filter(maxsize=maxsize))
        HyperTree.size_filter = overriding_method
        by_cartesian_product = set(tree.size_filter(maxsize=maxsize))
        assert by_recursion == by_cartesian_product


# Generic tree consistency assertions:


def assert_basic_tree_consistency(tree):
    assert_iteration_produces_members(tree)
    assert_leafs_have_no_children_and_root_has_no_parent(tree)
    assert_parent_and_children_are_inverse_of_each_other(tree)
    assert_tree_consists_of_children_and_root(tree)
    assert_tree_consists_of_leafs_and_linears_and_branches(tree)
    assert_level_is_length_of_root_path(tree)
    assert_root_is_outermost_and_leafs_are_innermost(tree)


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


def assert_tree_consists_of_children_and_root(tree):
    assert set(tree) == (
        set(chain.from_iterable(tree.children(x) for x in tree)) | {tree.root()}
    )


def assert_tree_consists_of_leafs_and_linears_and_branches(tree):
    assert set(tree) == (
        set(tree.leaf_nodes()) | set(tree.linear_nodes()) | set(tree.branch_nodes())
    )


def assert_level_is_length_of_root_path(tree):
    assert all(
        level == len(list(tree.root_path(node))) - 1 for node, level in tree.levels()
    )


def assert_root_is_outermost_and_leafs_are_innermost(tree):
    assert {tree.root()} == set(tree.outermost(tree))
    assert set(tree.leaf_nodes()) == set(tree.innermost(tree))
