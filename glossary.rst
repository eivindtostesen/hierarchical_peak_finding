Glossary
--------

argcut
  The argcut of a peak (valley) region is its argmin (argmax).
  See also cutoff.

argext
  The argext of a peak (valley) region is its argmax (argmin).
  See also extremum.

argmax
  The argmax of a region ``i:j`` is the index of the first occurrence of the region's maximum
  (at or after index ``i`` and before index ``j``).

argmin
  The argmin of a region ``i:j`` is the index of the first occurrence of the region's minimum
  (at or after index ``i`` and before index ``j``).

branching node
  A node that has two or more children.
  Implemented in method ``Tree.branch_nodes()``.
  See also leaf node and linear node.

cutoff
  The cutoff of a peak (valley) region is its minimum (maximum).
  See also argcut.

descendant node
  Node ``x`` is a descendant (child, grandchild, etc.)
  of node ``y`` if there is a path of parents from ``x`` to ``y``.

extremum
  The extremum of a peak (valley) region is its maximum (minimum).
  See also argext.

full node
  The outermost of all nodes having the same argext.
  A node is full if and only if it is not a main child.
  Implemented in the methods ``Tree.full()`` and ``Tree.full_nodes()``.
  See also tip node.

``i:j``
  Slice notation of a region with ``start = i`` (inclusive) and ``stop = j`` (exclusive).
  See also region.

innermost node
  An innermost node of a set of nodes is a
  node in the set that has no descendants in the set.
  Implemented in method ``Tree.innermost()``.
  See also outermost node.

``istop``
  The inclusive stop position of a region.
  Equal to the index of the last item in the region (``istop = stop - 1``).
  See also ``stop``.

lateral child
  A child node that has a different argext than its parent.
  Implemented in methods ``Tree.lateral()`` and ``Tree.lateral_descendants()``.
  See also main child.

leaf node
  A node that has no children.
  Implemented in method ``Tree.leaf_nodes()``.
  See also linear node and branching node.

level
  Number of steps to the root, in other words, a root node has level 0,
  its children have level 1, grandchildren level 2, etc.
  Implemented in method ``Tree.levels()``.

linear node
  A node that has one child.
  Implemented in method ``Tree.linear_nodes()``.
  See also leaf node and branching node.

local maximum
  A region containing equal values where all surrounding (``pre`` and ``post``) values
  are less than the value in the region.
  The definition is implemented in the method ``Region.is_local_maximum()``.
  See also peak.

local minimum
  A region containing equal values where all surrounding (``pre`` and ``post``) values
  are greater than the value in the region.
  The definition is implemented in the method ``Region.is_local_minimum()``.
  See also valley.

main child
  A child node that has the same argext as its parent.
  Implemented in methods ``Tree.main_child()`` and ``Tree.main_descendants()``.
  See also lateral child.

``max``
  Largest value in a region or numeric sequence.

``min``
  Smallest value in a region or numeric sequence.

nested region
  ``i:j`` is nested inside ``k:l``, written as ``i:j <= k:l``,
  if ``k<=i`` and ``j<=l``.
  Implemented in comparison methods ``Region.__le__()``, ``Region.__ge__()``,
  ``Region.__lt__()`` and ``Region.__gt__()``.  

node
  A ``Tree`` node is a ``Scope`` or other hashable object representing a peak or valley.
  A ``HyperTree`` node is a tuple of nodes.

outermost node
  An outermost node of a set of nodes is a
  node in the set that is not descendant of any node in the set.
  Implemented in method ``Tree.outermost()``.
  See also innermost node.

parent node
  Next node on path to root.

peak
  A region where all surrounding (``pre`` and ``post``) values
  are less than the region's minimum (cutoff) value.
  Implemented in the function ``find_peaks()`` and the method ``Region.is_peak()``.
  See also local maximum.

pose
  Six values specifying the position and orientation of a peak or valley:
  ``(start, istop, argext, argcut, extremum, cutoff)``.
  Implemented as namedtuple ``Pose``.

post
  Position(s) coming after/to the right of/succeeding a region.
  Implemented in method ``Region.post()``.

pre
  Position(s) coming before/to the left of/preceeding a region.
  Implemented in method ``Region.pre()``.

region
  A contiguous part of a numeric sequence.
  Implemented in classes ``Region`` and ``Scope``.
  See also ``i:j``.

root node
  A node that has no parent.

scope
  A region containing a peak or a valley.
  Implemented in the ``Scope`` class.
  See also peak and valley.

size
  The size of a region ``r`` is the difference between its maximum and minimum
  or equivalently: ``abs(r.extremum - r.cutoff)``.

``start``
  The start position of a region.
  Equal to the index of the first item in the region.

``stop``
  The exclusive stop position of a region.
  Equal to 1 + the index of the last item in the region (``stop = istop + 1``).
  See also ``istop``.

subarray
  A subarray is a new sequence corresponding to a region.
  Implemented in the method ``Region.subarray()``.

subtree
  A node and all its descendants.
  Implemented in the method ``Tree.subtree()``.

tip node
  The innermost of all nodes having the same argext.
  A node is a tip if and only if it has no main child.
  Implemented in the method ``Tree.tip()``.
  See also full node.

tree
  A graph of nested peaks or valleys.
  Implemented in the function ``tree()`` and classes ``Tree`` and ``HyperTree``.
  See also nested region.

valley
  A region where all surrounding (``pre`` and ``post``) values
  are greater than the region's maximum (cutoff) value.
  Implemented in the function ``find_valleys()`` and the method ``Region.is_valley()``.
  See also local minimum.

