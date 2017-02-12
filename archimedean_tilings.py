"""
This is a Python implementation of the ideas described in Chapter 19 (Archimedean Tilings) of the book
"Symmetries of Things" by John H. Conway, Heidi Burgiel, and Chaim Goodman-Strauss.
"""

import re


class DoilyFeature(object):
    ROTARY_ARM = 0
    FOLDED_BAND = 1
    HALF_ARM = 2
    UNTWISTED_BAND = 3
    TWISTED_BAND = 4
    HALF_BAND = 5


class Location(object):
    INTERIOR = 0
    BOUNDARY = 1


class EdgeLineSide(object):
    LOWER = 0
    UPPER = 1


class Edge:
    def __init__(self, index, endpoint, feature):
        self.index = index
        self.endpoint = endpoint
        self.feature = feature

    def __eq__(self, other):
        return other is not None and self.__dict__ == other.__dict__


class EdgeLine:
    def __init__(self, edge_index, side, permutation_symbol):
        if edge_index is not None and not isinstance(edge_index, int):
            raise ValueError("edge_index parameter to EdgeLine.__init__ must be None or int")
        if permutation_symbol is None:
            raise ValueError('Permutation_symbol=None in EdgeLine __init__')

        self.permutation_symbol = permutation_symbol
        if edge_index is not None:
            self.edge = permutation_symbol.edges[edge_index]
        else:
            self.edge = None
        self.side = side

    def __eq__(self, other):
        return other is not None and self.__dict__ == other.__dict__

    def __repr__(self):
        if self.edge:
            return str(self.edge.index) + ('U' if self.side == EdgeLineSide.UPPER else 'L')
        else:
            return 'B' + ('U' if self.side == EdgeLineSide.UPPER else 'L')

    def __hash__(self):
        return hash((None if not self.edge else self.edge.index, self.side))

    def touches_boundary(self):
        if self.edge is None:
            return True
        if self.edge.feature in (DoilyFeature.FOLDED_BAND, DoilyFeature.HALF_ARM):
            return True
        return False

    def connected_edgeline(self):
        edge = self.edge
        side = self.side
        ps = self.permutation_symbol

        if edge is None:  # boundary edgeline
            return None

        if edge.feature == DoilyFeature.ROTARY_ARM:
            return EdgeLine(edge.index, 1 - side, ps)
        elif edge.feature == DoilyFeature.UNTWISTED_BAND:
            return EdgeLine(edge.endpoint, 1 - side, ps)
        elif edge.feature == DoilyFeature.TWISTED_BAND:
            return EdgeLine(edge.endpoint, side, ps)
        elif edge.feature == DoilyFeature.HALF_BAND:
            assert (edge.index == 0 and side == EdgeLineSide.UPPER) or \
                   (edge.index == ps.num_edges - 1 and side == EdgeLineSide.LOWER)
            if edge.index == 0:
                return EdgeLine(ps.num_edges - 1, EdgeLineSide.LOWER, ps)
            elif edge.index == ps.num_edges - 1:
                return EdgeLine(0, EdgeLineSide.UPPER, ps)
        else:  # folded band or half arm
            return None

    def adjacent_edgeline(self):
        edge = self.edge
        side = self.side
        ps = self.permutation_symbol

        if edge is None and side == EdgeLineSide.LOWER:
            return EdgeLine(0, EdgeLineSide.LOWER, ps)
        elif edge is None and side == EdgeLineSide.UPPER:
            return EdgeLine(ps.num_edges - 1, EdgeLineSide.UPPER, ps)
        elif edge.index == 0 and side == EdgeLineSide.LOWER and ps.has_lower_boundary_edgeline():
            return EdgeLine(None, EdgeLineSide.LOWER, ps)
        elif edge.index == ps.num_edges - 1 and side == EdgeLineSide.UPPER and ps.has_upper_boundary_edgeline():
            return EdgeLine(None, EdgeLineSide.UPPER, ps)
        elif side == EdgeLineSide.UPPER:
            return EdgeLine((edge.index + 1) % ps.num_edges, EdgeLineSide.LOWER, ps)
        elif side == EdgeLineSide.LOWER:
            if edge.index > 0:
                return EdgeLine(edge.index - 1, EdgeLineSide.UPPER, ps)
            else:
                return EdgeLine(ps.num_edges - 1, EdgeLineSide.UPPER, ps)
        else:
            assert False, "We shouldn't get here in adjacent_edgeline"

    def touching_boundary_edgeline(self):
        assert self.touches_boundary()

        edge = self.edge
        side = self.side
        ps = self.permutation_symbol

        if edge is not None and edge.feature == DoilyFeature.FOLDED_BAND:
            return EdgeLine(edge.index, 1 - side, ps)
        if ps.num_half_arms == 0:
            assert edge is None
            return EdgeLine(None, 1 - side, ps)
        elif ps.num_half_arms == 1:
            assert side == EdgeLineSide.UPPER and (edge is None or edge.index == 0)
            if side == EdgeLineSide.UPPER and edge is None:
                return EdgeLine(0, EdgeLineSide.UPPER, ps)
            elif side == EdgeLineSide.UPPER and edge.index == 0:
                return EdgeLine(None, EdgeLineSide.UPPER, ps)
        elif ps.num_half_arms == 2:
            assert (edge.index == 0 and side == EdgeLineSide.UPPER) or \
                   (edge.index == ps.num_edges - 1 and side == EdgeLineSide.LOWER)
            if edge.index == 0 and side == EdgeLineSide.UPPER:
                return EdgeLine(ps.num_edges - 1, EdgeLineSide.LOWER, ps)
            elif edge.index == ps.num_edges - 1 and side == EdgeLineSide.LOWER:
                return EdgeLine(0, EdgeLineSide.UPPER, ps)

        assert False, "We shouldn't get here in touching_boundary_edgeline"


class Face:
    def __init__(self, index, location, edgelines):
        self.index = index
        self.location = location
        self.edgelines = edgelines

    def __eq__(self, other):
        return other is not None and self.__dict__ == other.__dict__

    def __repr__(self):
        return str(self.index) + ': ' + str(self.edgelines)

    def add_edgeline(self, edgeline, all_edgelines):
        if not isinstance(edgeline, EdgeLine):
            raise ValueError("edgeline parameter to Face.add_edgeline must be of class EdgeLine")
        self.edgelines.append(edgeline)
        del all_edgelines[all_edgelines.index(edgeline)]

    def number_of_sides(self):
        half_sides = 0
        for el in self.edgelines:
            if el.edge is not None:
                half_sides += 1
        if self.location == Location.BOUNDARY:
            sides = half_sides  # reflection multiplies by 2
        else:
            assert half_sides % 2 == 0
            sides = half_sides / 2
        return sides


class BoundaryComponent:
    def __init__(self):
        self.faces = []

    def __eq__(self, other):
        return other is not None and self.__dict__ == other.__dict__

    def __repr__(self):
        return str(self.faces)

    def add_face(self, face):
        self.faces.append(face)


class Orbifold:
    def __init__(self, gyrations, kaleidoscopes, handles, crosscaps):
        """
            Orbifolds are stored with the following normalizations.

            Gyrations: sorted and all 1's removed.
            Kaleidoscopes: all 1's removed.
            Handles and crosscaps: Crosscaps = 0,1,or 2; and number of handles increased accordingly.
        """
        self.gyrations = sorted([x for x in gyrations if x != 1])
        self.kaleidoscopes = [[x for x in k if x != 1] for k in kaleidoscopes]
        self.handles = handles
        self.crosscaps = crosscaps
        if self.crosscaps > 2:
            if self.crosscaps % 2 == 1:
                self.handles += (self.crosscaps - 1) / 2
                self.crosscaps = 1
            else:
                self.handles += (self.crosscaps - 2) / 2
                self.crosscaps = 2

    def __repr__(self):
        # display the signature in the order described in SoT, p. 27.
        signature = ''
        signature += 'o' * self.handles
        if self.gyrations:
            signature += '(' + ','.join(map(str, self.gyrations)) + ')'
        for kaleidoscope in self.kaleidoscopes:
            signature += '*(' + ','.join(map(str, kaleidoscope)) + ')'
        signature += 'x' * self.crosscaps
        return signature

    @staticmethod
    def cyclic_shift(x):
        if len(x) == 0:
            return x
        return x[1:] + [x[0]]

    @staticmethod
    def normalized_kaleidoscopes_equal(k1, k2):
        """

        :param k1: one kaleidoscope
        :param k2: another kaleidoscope
        :return:
            0 if the kaleidoscopes are not equivalent
            1 if the kaleidoscopes are equivalent with the same orientation, but not the opposite orientation
           -1 if the kaleidoscopes are equivalent with the opposite orientation, but not the same orientation
            2 if the kaleidoscopes are equivalent with both orientations
        """

        if len(k1) != len(k2):
            return 0
        if k1 == []:
            return 2

        same_orientation = False
        for i in xrange(len(k2)):
            if k1 == k2:
                same_orientation = True
                break
            k2 = Orbifold.cyclic_shift(k2)

        k2_reversed = [x for x in reversed(k2)]
        opposite_orientation = False
        for i in xrange(len(k2_reversed)):
            if k1 == k2_reversed:
                opposite_orientation = True
                break
            k2_reversed = Orbifold.cyclic_shift(k2_reversed)

        if same_orientation:
            if opposite_orientation:
                return 2
            return 1
        if opposite_orientation:
            return -1
        return 0

    def __eq__(self, other):
        if not isinstance(other, Orbifold):
            return False
        if self.handles != other.handles:
            return False
        if self.crosscaps != other.crosscaps:
            return False
        if self.gyrations != other.gyrations:
            return False
        if len(self.kaleidoscopes) != len(other.kaleidoscopes):
            return False

        matched_other_indices = range(len(other.kaleidoscopes))
        self_index = 0
        matching_orientation = 0
        while matched_other_indices:
            self_kaleidoscope = self.kaleidoscopes[self_index]
            match_found = False
            for i, m in enumerate(matched_other_indices):
                other_kaleidoscope = other.kaleidoscopes[m]
                match_result = Orbifold.normalized_kaleidoscopes_equal(self_kaleidoscope, other_kaleidoscope)
                if match_result == 0:
                    continue
                elif match_result == 2:
                    match_found = True
                elif self.crosscaps > 0:  # non-orientable orbifold, we don't care about orientation
                    match_found = True
                elif matching_orientation == 0:  # we don't yet have a defined orientation, set it to be this one
                    matching_orientation = match_result
                    match_found = True
                elif matching_orientation == match_result:
                    match_found = True
                if match_found:
                    break
            if not match_found:
                return False
            del matched_other_indices[i]
            self_index += 1
        return True


class PermutationSymbol:
    def __init__(self, symbol_string):
        self.symbol_string = symbol_string
        numeric_re = "(?: *[0-9]+ *(?:, *[0-9]+ *)?)"
        feature_re = "((?:\\(" + numeric_re + "\\))|(?:\\[" + numeric_re + "\\])|(?:<" + numeric_re + ">))"

        if not re.match("^" + feature_re + "+[.*]$", symbol_string):
            raise ValueError('Invalid Input')

        feature_strings = re.findall(feature_re, symbol_string)
        max_number = -1
        edge_dict = {}

        self.is_orientable = True
        self.has_half_band = False
        self.num_half_arms = 0
        self.num_rotary_arms = 0
        half_arms = []
        half_band = []
        for fs in feature_strings:
            numbers = fs[1:-1].split(',')
            n1 = int(numbers[0].strip())
            n2 = None
            if len(numbers) == 2:
                n2 = int(numbers[1].strip())

            if n1 in edge_dict:
                raise ValueError('Duplicate value %d in permutation symbol %s' % (n1, symbol_string))
            if n2 in edge_dict:
                raise ValueError('Duplicate value %d in permutation symbol %s' % (n2, symbol_string))

            if n2 is None:
                if fs[0] == '(':
                    edge_dict[n1] = Edge(n1, None, DoilyFeature.ROTARY_ARM)
                    self.num_rotary_arms += 1
                elif fs[0] == '[':
                    edge_dict[n1] = Edge(n1, None, DoilyFeature.FOLDED_BAND)
                elif fs[0] == '<':
                    self.num_half_arms += 1
                    if self.num_half_arms > 2:
                        raise ValueError('Cannot have more than 2 half arms in permutation symbol %s' % symbol_string)
                    edge_dict[n1] = Edge(n1, None, DoilyFeature.HALF_ARM)
                    half_arms.append(n1)
            else:
                feature = None
                if fs[0] == '(':
                    feature = DoilyFeature.UNTWISTED_BAND
                elif fs[0] == '[':
                    feature = DoilyFeature.TWISTED_BAND
                    self.is_orientable = False
                elif fs[0] == '<':
                    if self.has_half_band:
                        raise ValueError('Cannot have more than 1 half band in permutation symbol %s' % symbol_string)
                    self.has_half_band = True
                    feature = DoilyFeature.HALF_BAND
                    half_band = (min(n1, n2), max(n1, n2))

                edge_dict[n1] = Edge(n1, n2, feature)
                edge_dict[n2] = Edge(n2, n1, feature)

            if n1 > max_number:
                max_number = n1
            if n2 is not None and n2 > max_number:
                max_number = n2

        # Constraints on half arms:
        #   Cannot be more than 2.
        #   Must be the first or last edge.  If only one, must be the first (normalization).
        #
        # Constraints on half bands:
        #   Can only be one.
        #   Must join the first and last edges.
        if self.has_half_band:
            if not (half_band[0] == 0 and half_band[1] == max_number):
                raise ValueError('Half band must join first and last edges in permutation symbol %s' % symbol_string)
        elif self.num_half_arms == 1:
            if half_arms[0] != 0:
                raise ValueError('A single half arm must be the first edge in permutation symbol %s' % symbol_string)
        elif self.num_half_arms == 2:
            half_arms = sorted(half_arms)
            if not (half_arms[0] == 0 and half_arms[1] == max_number):
                raise ValueError('Two half arms must be first and last edges in permutation symbol %s' % symbol_string)

        assert symbol_string[-1] in ('.', '*')
        if symbol_string[-1] == '.':
            if self.has_half_band or self.num_half_arms > 0:
                raise ValueError('Local symmetry in %s must be * when there is a half-band or half-arm' % symbol_string)
            self.local_symmetry = Location.INTERIOR
        elif symbol_string[-1] == '*':
            self.local_symmetry = Location.BOUNDARY

        self.edges = []
        for i in xrange(max_number + 1):
            if i not in edge_dict:
                raise ValueError('Edge number %d missing from permutation symbol %s' % (i, symbol_string))
            self.edges.append(edge_dict[i])
        self.num_edges = len(self.edges)

        self.set_face_info()

    def has_lower_boundary_edgeline(self):
        if self.num_half_arms == 0 and not self.has_half_band and self.local_symmetry == Location.BOUNDARY:
            return True
        return False

    def has_upper_boundary_edgeline(self):
        if self.num_half_arms == 1:
            assert self.edges[0].feature == DoilyFeature.HALF_ARM
            return True
        if self.num_half_arms == 0 and not self.has_half_band and self.local_symmetry == Location.BOUNDARY:
            return True
        return False

    @staticmethod
    def update_edgeline_to_face(edgeline_to_face, edgeline, adjacent_edgeline, current_face):
        if edgeline.edge is None and edgeline.side == EdgeLineSide.LOWER:
            key = edgeline
        elif adjacent_edgeline.edge is None and adjacent_edgeline.side == EdgeLineSide.LOWER:
            key = adjacent_edgeline
        elif edgeline.edge is not None and edgeline.side == EdgeLineSide.UPPER:
            key = edgeline
        elif adjacent_edgeline.edge is not None and adjacent_edgeline.side == EdgeLineSide.UPPER:
            key = adjacent_edgeline
        else:
            assert False, 'Could not find a suitable key in update_edgeline_to_face'

        edgeline_to_face[key] = current_face

    @staticmethod
    def fundamental_iteration(edgeline, current_face, all_edgelines, edgeline_to_face):
        while True:
            connected_edgeline = edgeline.connected_edgeline()
            if not connected_edgeline or connected_edgeline not in all_edgelines:
                break  # this face is done
            current_face.add_edgeline(connected_edgeline, all_edgelines)

            edgeline = connected_edgeline
            adjacent_edgeline = edgeline.adjacent_edgeline()
            if not adjacent_edgeline:  # we hit a boundary wall
                break

            # we returned to where we started.  This still needs an update to edgeline_to_face before we break
            PermutationSymbol.update_edgeline_to_face(edgeline_to_face, edgeline, adjacent_edgeline, current_face)
            if adjacent_edgeline not in all_edgelines:
                break  # this face is done
            current_face.add_edgeline(adjacent_edgeline, all_edgelines)

            edgeline = adjacent_edgeline

    # Return all edgelines in order.
    def get_all_edgelines(self):
        all_edgelines = []

        if self.has_lower_boundary_edgeline():
            all_edgelines.append(EdgeLine(None, EdgeLineSide.LOWER, self))
        for i in xrange(self.num_edges):
            if 0 < i < self.num_edges - 1:
                all_edgelines.append(EdgeLine(i, EdgeLineSide.LOWER, self))
                all_edgelines.append(EdgeLine(i, EdgeLineSide.UPPER, self))
            elif i == 0:
                if not (self.has_half_band or self.edges[0].feature == DoilyFeature.HALF_ARM):
                    all_edgelines.append(EdgeLine(i, EdgeLineSide.LOWER, self))
                all_edgelines.append(EdgeLine(i, EdgeLineSide.UPPER, self))
            elif i == self.num_edges - 1:
                all_edgelines.append(EdgeLine(i, EdgeLineSide.LOWER, self))
                if not (self.has_half_band or self.edges[self.num_edges - 1].feature == DoilyFeature.HALF_ARM):
                    all_edgelines.append(EdgeLine(i, EdgeLineSide.UPPER, self))
        if self.has_upper_boundary_edgeline():
            all_edgelines.append(EdgeLine(None, EdgeLineSide.UPPER, self))

        return all_edgelines

    @staticmethod
    def get_next_boundary_edgeline(all_edgelines):
        # We must start each component with either a lower boundary edge, or an upper non-boundary edge.
        # This ensures (in the orientable case) that we traverse all boundary components with the same orientation
        # (clockwise or counterclockwise).  This matters in the orientable case when we're
        # calculating the orbifold.
        edgeline = None
        for el in all_edgelines:
            if el.touches_boundary():
                if (el.edge is None and el.side == EdgeLineSide.LOWER) or \
                        (el.edge is not None and el.side == EdgeLineSide.UPPER):
                    edgeline = el
                    break
        return edgeline

    def get_boundary_components(self, all_edgelines, edgeline_to_face):
        boundary_components = []
        num_faces = 0

        edgeline = self.get_next_boundary_edgeline(all_edgelines)
        while edgeline:
            current_boundary_component = BoundaryComponent()
            while edgeline in all_edgelines:  # loop until we cycle around to the beginning of this boundary component
                current_face = Face(num_faces, Location.BOUNDARY, [])
                num_faces += 1
                current_face.add_edgeline(edgeline, all_edgelines)

                adjacent_edgeline = edgeline.adjacent_edgeline()
                current_face.add_edgeline(adjacent_edgeline, all_edgelines)
                self.update_edgeline_to_face(edgeline_to_face, edgeline, adjacent_edgeline, current_face)
                edgeline = adjacent_edgeline

                self.fundamental_iteration(edgeline, current_face, all_edgelines, edgeline_to_face)

                current_boundary_component.add_face(current_face)

                last_edgeline = current_face.edgelines[-1]
                edgeline = last_edgeline.touching_boundary_edgeline()
            boundary_components.append(current_boundary_component)

            # Find a starting boundary edgeline for the next boundary component
            edgeline = self.get_next_boundary_edgeline(all_edgelines)

        if self.has_half_band:  # the edge of the half-band is a boundary component that touches no face.
            boundary_components.append(BoundaryComponent())
        return boundary_components

    def get_interior_faces(self, all_edgelines, edgeline_to_face, num_faces):
        interior_faces = []

        # Find interior faces
        # We assume that all boundary components have already been found, so all_edgelines no longer contains
        # any edges that touch the boundary.  We start with an upper edge and repeatedly apply the fundamental
        # iteration.  Starting with an upper edge guarantees (in the orientable case) that we traverse all faces
        # with the same orientation.
        while all_edgelines:
            for edgeline in all_edgelines:
                assert edgeline.edge is not None
                if edgeline.side == EdgeLineSide.UPPER:
                    break
            current_face = Face(num_faces, Location.INTERIOR, [])
            num_faces += 1
            current_face.add_edgeline(edgeline, all_edgelines)
            self.fundamental_iteration(edgeline, current_face, all_edgelines, edgeline_to_face)
            interior_faces.append(current_face)
        return interior_faces

    def get_face_code(self, all_edgelines, edgeline_to_face):
        face_code = []
        for el in all_edgelines:  # sorted
            # edgeline_to_face should only be on the interior upper edges, except for the lower boundary edge
            if el.side == EdgeLineSide.LOWER and el.edge is not None:  # skip all non-boundary lower edges
                continue
            if el.side == EdgeLineSide.UPPER and el.edge is None:  # skip upper boundary edge
                continue

            face = edgeline_to_face[el]
            face_code.append((face.number_of_sides(), face.index))

        if self.has_half_band or self.num_half_arms == 2:
            for fc in reversed(face_code):
                face_code.append(fc)
        else:
            has_boundary = False
            interior_face_code = face_code
            if self.has_lower_boundary_edgeline():
                has_boundary = True
                interior_face_code = interior_face_code[1:]
            if self.has_upper_boundary_edgeline():
                has_boundary = True
                interior_face_code = interior_face_code[:-1]
            if has_boundary:
                for fc in reversed(interior_face_code):
                    face_code.append(fc)

        if self.has_lower_boundary_edgeline():
            # in this case, we have to rotate face_code one position to the left, so that the
            # face between 0 and 1 is listed first, per the normalization defined in SoT, p. 252.
            first = face_code[0]
            face_code = face_code[1:]
            face_code.append(first)

        return face_code

    def get_orbifold(self):
        gyrations = []
        kaleidoscopes = []
        vertex_parameter_added = False

        # gyrations:
        #   each rotary arm gets a 2
        #   each interior faces gets a parameter
        #   if local symmetry is '.', then local vertex parameter 'n' is a gyration
        gyrations.extend([2]*self.num_rotary_arms)
        for f in self.interior_faces:
            gyrations.append('[' + str(f.index) + ']')
        if self.local_symmetry == Location.INTERIOR:
            gyrations.append('n')
            vertex_parameter_added = True

        # kaleidoscopes:
        #    each boundary component is a kaleidoscope
        #    each face is a parameter
        #    half arm at 0 gets an additional 2 prepended
        #    half arm at m-1 gets an additional 2 appended
        #    if local symmetry is '*', then component containing m-1 (if not half-band) gets local vertex parameter 'n'
        #       appended

        for bc in self.boundary_components:
            kaleidoscope = []
            has_lower_half_arm = False
            has_upper_half_arm = False
            has_upper_boundary = False
            if len(bc.faces) == 0:  # half-band
                has_upper_boundary = True
            for f in bc.faces:
                kaleidoscope.append('[' + str(f.index) + ']')
                edge_indices = [el.edge.index for el in f.edgelines if el.edge is not None]
                if 0 in edge_indices and self.num_half_arms > 0:
                    has_lower_half_arm = True
                if self.num_edges - 1 in edge_indices and not self.has_half_band:
                    has_upper_boundary = True
                    if self.num_half_arms == 2:
                        has_upper_half_arm = True
            if has_lower_half_arm:
                kaleidoscope.insert(0, 2)
            if has_upper_half_arm and self.num_edges > 1:
                kaleidoscope.append(2)
            if self.local_symmetry == Location.BOUNDARY and has_upper_boundary:
                assert not vertex_parameter_added, 'Attempted to add vertex parameter twice'
                kaleidoscope.append('n')
                vertex_parameter_added = True
            kaleidoscopes.append(kaleidoscope)
        assert vertex_parameter_added, 'Vertex parameter was never added'

        # handles and crosscaps:
        #    Compute s = 2 - Euler characteristic - number of boundary components.
        #    If orientable, then s must be even, and we have s/2 handles and 0 crosscaps
        #    If non-orientable, and s is even, then we have (s-2)/2 handles and 2 crosscaps
        #    If non-orientable, and s is odd, then we have (s-1)/2 handles and 1 crosscap.

        if self.local_symmetry == Location.INTERIOR:
            orbifold_half_vertices = 2
        else:
            orbifold_half_vertices = 1
        orbifold_half_edges = 0
        for e in self.edges:
            if e.feature == DoilyFeature.FOLDED_BAND:
                orbifold_half_edges += 2
            elif e.feature in (DoilyFeature.ROTARY_ARM, DoilyFeature.HALF_ARM, DoilyFeature.HALF_BAND):
                continue
            else:  # twisted or untwisted bands
                orbifold_half_edges += 1  # these bands will be added twice and so contribute 2
        if self.local_symmetry == Location.BOUNDARY:
            orbifold_half_edges += 1  # the boundary component containing the half-vertex is another half-edge
        orbifold_faces = self.num_faces
        orbifold_boundaries = len(self.boundary_components)

        double_euler_characteristic = orbifold_half_vertices - orbifold_half_edges + 2 * orbifold_faces
        assert double_euler_characteristic % 2 == 0
        euler_characteristic = double_euler_characteristic / 2
        assert euler_characteristic <= 2
        s = 2 - euler_characteristic - orbifold_boundaries
        if self.is_orientable:
            assert s % 2 == 0 and s >= 0
            handles = s/2
            crosscaps = 0
        else:  # non-orientable
            assert s > 0
            if s % 2 == 0:
                handles = (s-2)/2
                crosscaps = 2
            else:
                handles = (s-1)/2
                crosscaps = 1

        return Orbifold(gyrations, kaleidoscopes, handles, crosscaps)

    def set_face_info(self):
        all_edgelines = self.get_all_edgelines()
        copy_all_edgelines = [x for x in all_edgelines]  # we'll be destroying all_edgelines throughout the algorithm

        edgeline_to_face = {}
        self.boundary_components = self.get_boundary_components(all_edgelines, edgeline_to_face)
        num_faces = sum([len(x.faces) for x in self.boundary_components])
        self.interior_faces = self.get_interior_faces(all_edgelines, edgeline_to_face, num_faces)
        num_faces += len(self.interior_faces)
        self.num_faces = num_faces

        self.face_code = self.get_face_code(copy_all_edgelines, edgeline_to_face)

        self.orbifold = self.get_orbifold()