# -*- coding: utf-8 -*-
# This file is part of Peakoscope.
# Copyright (C) 2021-2025  Eivind Tøstesen
# Peakoscope is licensed under GPLv3.
# SPDX-License-Identifier: GPL-3.0-or-later


import pytest
from operator import attrgetter
import peakoscope
import peakoscope.testing as testing
from peakoscope import Tree, Forest, HyperTree, HyperForest, Scope, Scope6, find_peaks


# Assert functions:


def assert_all_assertions(tree):
    """Assert all tree assertions in peakoscope.testing."""
    testing.assert_iteration_produces_members(tree)
    testing.assert_leafs_have_no_children_and_root_has_no_parent(tree)
    testing.assert_parent_and_children_are_inverse_of_each_other(tree)
    testing.assert_level_is_length_of_root_path(tree)
    testing.assert_root_is_outermost_and_leafs_are_innermost(tree)
    # tree partitions:
    testing.assert_tree_consists_of_children_and_root(tree)
    testing.assert_tree_consists_of_leafs_and_linears_and_branches(tree)
    testing.assert_tree_consists_of_full_nodes_and_main_descendants(tree)
    testing.assert_full_nodes_consists_of_lateral_descendants_plus_root(tree)
    testing.assert_children_consist_of_main_child_plus_lateral_children(tree)
    testing.assert_tree_consists_of_main_paths(tree)
    # tip and full:
    testing.assert_being_tip_is_having_no_main_child(tree)
    testing.assert_main_child_keeps_tip_and_lateral_changes_tip(tree)
    testing.assert_main_path_shares_tip_and_full(tree)
    testing.assert_root_is_full_and_leafs_are_tips(tree)
    # node size:
    testing.assert_parent_size_is_strictly_greater(tree)
    testing.assert_if_local_extremum_then_leaf(tree)
    # size filter:
    testing.assert_size_filter_equals_outermost_of_below_maxsize(tree)
    testing.assert_size_filter_equals_definition(tree)


def assert_tree_nodes_are_nested_scopes(tree):
    """Assert nesting of regions: parent > child and root > full > tip."""
    assert all(x < tree.parent(x) for x in tree if tree.is_nonroot(x))
    assert all(y < x for x in tree for y in tree.children(x))
    assert all(tree.tip(x) <= x <= tree.full(x) <= tree.root() for x in tree)


def assert_main_children_keep_argext_lateral_children_move_away(tree):
    """Assert that main/lateral child has same/different argext as parent."""
    assert all(x.argext == tree.parent(x).argext for x in tree.main_descendants())
    assert all(x.argext != tree.parent(x).argext for x in tree.lateral_descendants())


def assert_dict_of_dicts_of_same_length(tree):
    """Assert tree or forest length equals dict sizes."""
    assert all(
        len(tree) == len(d) for d in tree.as_dict_of_dicts().values() if type(d) == dict
    )


def assert_parents_and_children_belong_to_tree(tree):
    """Assert that hyper parent and children methods produce tree or forest nodes.

    See also Proposition 4.
    """
    # parent in tree:
    assert all(tree.parent(x) in tree for x in tree if tree.is_nonroot(x))
    # children in tree:
    assert all(x in tree for y in tree for x in tree.children(y))


def assert_grid_nodes_belong_to_tree(tree):
    """Assert that hyper size_filter produces tree or forest nodes.

    See also Proposition 7.
    """
    assert all(x in tree for x in tree.size_filter())


def assert_leafs_belong_to_tree(tree):
    """Assert that hyper leaf_nodes produces tree or forest nodes."""
    assert all(x in tree for x in tree.leaf_nodes())


def assert_parent_size_equals_minimum(tree):
    """Assert that hyper parent size is the smallest of L/R parent sizes.

    See also Proposition 2.
    """
    assert all(
        tree.size(tree.parent((a, b)))
        == min(tree.L.size(tree.L.parent(a)), tree.R.size(tree.R.parent(b)))
        for (a, b) in tree
        if tree.L.is_nonroot(a) and tree.R.is_nonroot(b)
    )


# parameter fixtures:


@pytest.fixture(params=[0.0, 0.01, 0.2, 0.6, 2.0])
def fraction(request):
    return request.param


@pytest.fixture(params=[False, True])
def are_valleys(request):
    return request.param


@pytest.fixture(params=[False, True])
def is_forest(request):
    return request.param


def prune1(f):
    """Strip input forest of roots and some leafs."""
    return f - f.roots() - (set(f.leaf_nodes()) - set(f.full_nodes()))


@pytest.fixture(params=[lambda f: f, prune1])
def prune(request):
    return request.param


# test Tree objects (one-dimensional):


@pytest.fixture()
def tree(data0, are_valleys):
    return peakoscope.tree(data0, valleys=are_valleys)


def test_tree(tree):
    """Test Tree objects."""
    assert_all_assertions(tree)
    assert_tree_nodes_are_nested_scopes(tree)
    assert_main_children_keep_argext_lateral_children_move_away(tree)
    assert_dict_of_dicts_of_same_length(tree)


# test Forest objects (one-dimensional):


@pytest.fixture()
def forest(all_data, are_valleys, prune):
    return prune(peakoscope.tree(all_data, valleys=are_valleys, forest=True))


def test_forest(forest):
    """Test Forest objects."""
    assert_all_assertions(forest)
    assert_main_children_keep_argext_lateral_children_move_away(forest)
    assert_dict_of_dicts_of_same_length(forest)


def test_tree_nodes_are_nested_scopes(forest):
    """Assert nesting of regions: parent > child and root > full > tip."""
    tree = forest
    assert all(x < tree.parent(x) for x in tree if tree.is_nonroot(x))
    assert all(y < x for x in tree for y in tree.children(x))
    assert all(tree.tip(x) <= x <= tree.full(x) <= tree.root(x) for x in tree)


def test_forest_operations(forest):
    """Test set operations on forests."""
    assert repr(forest) == repr(Forest(are_valleys=forest.are_valleys) | forest)
    assert repr(forest & forest.leaf_nodes()) == repr(
        forest - forest.linear_nodes() - forest.branch_nodes()
    )


# test HyperTree objects (pairs of 1D trees):


@pytest.fixture()
def hypertree(data1, data2, are_valleys):
    return HyperTree(
        peakoscope.tree(data1, valleys=are_valleys),
        peakoscope.tree(data2, valleys=are_valleys),
    )


def test_hypertree(hypertree):
    """Test HyperTree objects."""
    assert_all_assertions(hypertree)
    assert_parents_and_children_belong_to_tree(hypertree)
    assert_grid_nodes_belong_to_tree(hypertree)
    assert_leafs_belong_to_tree(hypertree)
    assert_parent_size_equals_minimum(hypertree)


def test_recursion_equals_cartesian_product(hypertree, monkeypatch, fraction):
    """Assert that two algorithms are equivalent."""
    tree = hypertree
    rootsize = tree.size(tree.root())
    # output of size_filter:
    by_cartesian_product = set(tree.size_filter(maxsize=fraction * rootsize))
    # now change the algorithm behind size_filter:
    monkeypatch.setattr(HyperTree, "size_filter", Tree.size_filter)
    # output of size_filter reloaded:
    by_recursion = set(tree.size_filter(maxsize=fraction * rootsize))
    # output are equal as sets of nodes:
    assert by_recursion == by_cartesian_product


# test HyperForest objects (pairs of 1D forests):


@pytest.fixture()
def hyperforest(data1, data2, are_valleys, prune):
    return HyperForest(
        peakoscope.tree(data1, valleys=are_valleys, forest=True),
        prune(peakoscope.tree(data2, valleys=are_valleys, forest=True)),
    )


def test_hyperforest(hyperforest):
    """Test HyperForest objects."""
    assert_all_assertions(hyperforest)
    assert_parents_and_children_belong_to_tree(hyperforest)
    assert_grid_nodes_belong_to_tree(hyperforest)
    assert_parent_size_equals_minimum(hyperforest)


def test_recursion_equals_cartesian_product2(hyperforest, monkeypatch, fraction):
    """Assert that two algorithms are equivalent."""
    tree = hyperforest
    rootsize = max((tree.size(root) for root in tree.roots()), default=0)
    # output of size_filter:
    by_cartesian_product = set(tree.size_filter(maxsize=fraction * rootsize))
    # now change the algorithm behind size_filter:
    monkeypatch.setattr(HyperForest, "size_filter", Forest.size_filter)
    # output of size_filter reloaded:
    by_recursion = set(tree.size_filter(maxsize=fraction * rootsize))
    # output are equal as sets of nodes:
    assert by_recursion == by_cartesian_product


# special tests with special data:


def test_identity_element(all_data, flat_data):
    """Test that a Forest of flat data is an identity element."""
    forest = peakoscope.tree(all_data, forest=True)
    unity = peakoscope.tree(flat_data, forest=True)
    assert all(
        f == fu[0] == uf[1] for f, fu, uf in zip(forest, forest @ unity, unity @ forest)
    )
    assert all(
        forest.size(f) == (forest @ unity).size(fu) == (unity @ forest).size(uf)
        for f, fu, uf in zip(forest, forest @ unity, unity @ forest)
    )


def test_zero_element(all_data):
    """Test that an empty Forest is a zero element."""
    forest = peakoscope.tree(all_data, forest=True)
    zero = Forest()
    assert len(zero) == len(forest @ zero) == len(zero @ forest) == 0


def test_non_linear_tree(zigzag_data, is_forest):
    """Assert that zigzag data gives tree or forest without linear nodes."""
    assert len(list(peakoscope.tree(zigzag_data, forest=is_forest).linear_nodes())) == 0


def test_str_forest(data1):
    """Assert that forest str is a concatenation of tree str."""
    forest1 = prune1(peakoscope.tree(data1, forest=True))
    assert str(forest1) == "\n".join(
        str(Tree.from_peaks(forest1.subtree(root), presorted=False))
        for root in forest1.roots()
    )


def test_eval_repr_tree(data1):
    """Assert that tree repr is readable by eval."""
    Scope.default_data = data1
    tree = peakoscope.tree(data1)
    # Tree:
    assert repr(tree) == repr(eval(repr(tree)))
    # HyperTree:
    assert repr(tree @ tree) == repr(eval(repr(tree @ tree)))
    assert repr(tree @ tree @ tree) == repr(eval(repr(tree @ tree @ tree)))
    Scope.default_data = None


def test_eval_repr_forest(all_data, prune):
    """Assert that forest repr is readable by eval."""
    Scope.default_data = all_data
    forest = prune(peakoscope.tree(all_data, forest=True))
    # Forest:
    assert repr(forest) == repr(eval(repr(forest)))
    # HyperForest:
    assert repr(forest @ forest) == repr(eval(repr(forest @ forest)))
    Scope.default_data = None


def test_children_are_sorted_by_extremum(data1, is_forest, are_valleys):
    """Assert peak children sorted by max, valley children sorted by min."""
    tree = peakoscope.tree(data1, valleys=are_valleys, forest=is_forest)
    assert all(
        list(tree.children(x))
        == sorted(tree.children(x), key=attrgetter("extremum"), reverse=not are_valleys)
        for x in tree
    )


def test_tree_set_nodes(data1):
    """Test the method Tree.set_nodes."""
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


def test_forest_set_nodes(data1):
    """Test the method Forest.set_nodes."""
    # create forest with nodes of type Scope:
    scope_tree = peakoscope.tree(data1, forest=True)
    # create forest with nodes of type Scope6:
    scope6_tree = Forest.from_peaks(map(lambda t: Scope6(*t), find_peaks(data1)))
    # lists of nodes are NOT equal:
    assert list(scope_tree) != list(scope6_tree)
    # change from Scope6 to Scope:
    scope6_tree.set_nodes({s6: Scope.from_attrs(s6, data1) for s6 in scope6_tree})
    # now lists of nodes are equal:
    assert list(scope_tree) == list(scope6_tree)
