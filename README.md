Archimedean Tilings
===================

This is a Python implementation of the ideas described in Chapter 19 (_Archimedean Tilings_) of the book
_Symmetries of Things_ by John H. Conway, Heidi Burgiel, and Chaim Goodman-Strauss.  All page number
references below are for this book.

First Example
-------------
Code to be input to the Python interpreter is preceded by `>>>` below, and the result is shown on the next line.

    >>> from archimedean_tilings import PermutationSymbol

Our first example is the permutation symbol `[0](1,2)[3,4].` which appears on page 254.

Note that the notation differs slightly from the text: When there are two numbers within `()`, `[]`, or `<>`, we
separate them with a comma.  We use a period `.` or asterisk `*` to denote the local symmetry at the vertex (cyclic or dihedral,
respectively).

    >>> ps = PermutationSymbol('[0](1,2)[3,4].')

### Face Code

    >>> ps.face_code
    [(8, 0), (1, 1), (8, 0), (8, 0), (8, 0)]

In the text (bottom of page 254), the face code is written as `(8b,a,8b,8b,8b)^n`. This differs from our notation above
slightly.

First, we do not show the `^n` portion at all.  It is implicitly understood that all face codes should be repeated `n` times, where `n`
is the parameter at the vertex.

Next, the product `8b` has become `(8, 0)`. The 0 is a face index (see details below).  The pair `(8, 0)` indicates that
we should multiply 8 by the parameter (gyration or kaleidoscopic point) associated with face 0.  So face 0 corresponds
to the face touching the kaleidoscopic point `b` in the diagram at the bottom of page 254.  Similarly, face 1 corresponds
to the interior face with gyration point `a`.

### Orbifold

Next, we display the signature of the orbifold associated to this permutation symbol.

    >>> ps.orbifold
    ([1],n)*([0])x

This is similar to the orbifold notation in the text discussed in Chapter 2 (see page 27 for a summary).

The first part, `([1],n)` is the list of gyrations. The number enclosed in square brackets `[1]` means the parameter associated
to face index 1 (which is called `a` in the text).  The symbol `n` always denotes the local parameter at the vertex.
Because our original permutation symbol ends in a dot, indicating cyclic symmetry, the local parameter is a gyration point.

The next part, `*([0])` indicates a kaleidoscope.  As in the text, `*` indicates the start of a kaleidoscope.  Here, the
kaleidoscope has only one value associated with it, `[0]`, meaning the parameter associated to face index 0 (which is called
`b` in the text).

Finally, the `x` indicates a crosscap (miracle), as in the text.

This example did not have any handles (wonders), but if it did, an appropriate number of `o` symbols would precede the gyration points.

### Boundary components

Next, we show the boundary components (kaleidoscopes) of the orbifold.

    >>> ps.boundary_components
    [[0: [0U, 1L, 2U, 3L, 4L, 3U, 4U, 0L]]]

The value of `ps.boundary_components` is a list of boundary components. Each boundary component consists of a list of
faces, and each face consists of a face index (preceding the colon) and a list of edge lines (discussed below).  In
this case, we have only one kaleidoscope, and there is only one face touching the kaleidoscope, which is face 0.

The list `[0U, 1L, 2U, 3L, 4L, 3U, 4U, 0L]` is a list of _edge lines_ describing face 0.  In an edge line like `2U`,
the number 2 is an edge index, the same index as in the original permutation symbol.  The symbols `L` and `U` denote
the two sides of the edge ("lower" and "upper").  Each edgeline corresponds to half of a side in the original tiling.

### Interior faces

    >>> ps.interior_faces
    [1: [1U, 2L]]

The value of `ps.interior_faces` is a list of interior faces.  The notation is the same as for boundary faces.

### Edgelines of a face

We can programmatically access the list of edgelines of a face.

    >>> component = ps.boundary_components[0]
    >>> face_0 = component.faces[0]
    >>> face_0
    0: [0U, 1L, 2U, 3L, 4L, 3U, 4U, 0L]
    >>> face_0.index
    0
    >>> face_0.edgelines
    [0U, 1L, 2U, 3L, 4L, 3U, 4U, 0L]
    >>> face_0.edgelines[3].side
    0

0 denotes the lower side, 1 denotes the upper side.

    >>> face_0.edgelines[3].edge.index
    3
    >>> face_0.edgelines[3].edge.feature
    4

4 denotes a twisted band.  The mapping between features and integers is defined in the `DoilyFeature` class in the source code.

### Features of the orbifold

We can programmatically access the features of the orbifold.

    >>> ps.orbifold.gyrations
    ['[1]', 'n']
    >>> ps.orbifold.kaleidoscopes
    [['[0]']]
    >>> ps.orbifold.handles
    0
    >>> ps.orbifold.crosscaps
    1

Second Example
--------------

The next example is `<0,8>(1)[2,6](3,4)(5,7)*` which appears on pages 255-256.

    >>> ps2 = PermutationSymbol('<0,8>(1)[2,6](3,4)(5,7)*')
    >>> ps2.orbifold
    (2,[0],[1])*(n)xx

Notice that the `2` in the gyration component above is not enclosed in square brackets, so it denotes a fixed integer
`2`.  This `2` gyration point comes from the rotary arm `(1)` in the permutation symbol.

Because the local symmetry is dihedral, the vertex parameter `n` occurs in a kaleidoscope.

    >>> ps2.face_code
    [(7, 0), (7, 0), (7, 0), (1, 1), (7, 0), (7, 0), (7, 0), (7, 0), (7, 0), (7, 0), (7, 0), (7, 0), (1, 1), (7, 0), (7, 0), (7, 0)]

    >>> ps2.boundary_components
    [[]]

We have one boundary component, but it _touches no faces_. Such a boundary component comes from the
half-band `<0,8>`.

    >>> ps2.interior_faces
    [0: [0U, 8L, 7U, 5L, 4U, 3L, 2U, 6U, 7L, 5U, 6L, 2L, 1U, 1L], 1: [3U, 4L]]

    >>> ps2.orbifold.gyrations
    [2, '[0]', '[1]']
    >>> ps2.orbifold.kaleidoscopes
    [['n']]
    >>> ps2.orbifold.handles
    0
    >>> ps2.orbifold.crosscaps
    2

