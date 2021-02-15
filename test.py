import unittest
import re
from archimedean_tilings import PermutationSymbol, Orbifold


class TestBadPermutationSymbols(unittest.TestCase):
    def test_invalid_half_band_1(self):
        with self.assertRaises(ValueError):
            PermutationSymbol('<1,3>[0][2]*')

    def test_invalid_half_band_2(self):
        with self.assertRaises(ValueError):
            PermutationSymbol('<0,2>[1][3]*')

    def test_invalid_half_arms(self):
        with self.assertRaises(ValueError):
            PermutationSymbol('[0][1]<2>*')

    def test_too_many_half_bands(self):
        with self.assertRaises(ValueError):
            PermutationSymbol('<0,2><1,3>*')

    def test_too_many_half_arms(self):
        with self.assertRaises(ValueError):
            PermutationSymbol('<0><1><2>*')


class TestKaleidoscopeEquality(unittest.TestCase):
    def test_exactly_same(self):
        self.assertEqual(Orbifold.normalized_kaleidoscopes_equal([2, 3, 4, 5], [2, 3, 4, 5]), 1)

    def test_empty_kaleidoscopes(self):
        self.assertEqual(Orbifold.normalized_kaleidoscopes_equal([], []), 2)

    def test_different_lengths(self):
        self.assertEqual(Orbifold.normalized_kaleidoscopes_equal([2, 3, 4], [2, 3]), 0)

    def test_equal_lengths_totally_different(self):
        self.assertEqual(Orbifold.normalized_kaleidoscopes_equal([2, 3, 4], [5, 6, 7]), 0)

    def test_cyclic_shift(self):
        self.assertEqual(Orbifold.normalized_kaleidoscopes_equal([2, 3, 4, 5, 6], [4, 5, 6, 2, 3]), 1)

    def test_reverse_cyclic_shift(self):
        self.assertEqual(Orbifold.normalized_kaleidoscopes_equal([2, 3, 4, 5, 6], [4, 3, 2, 6, 5]), -1)

    def test_permutation_unequal(self):
        self.assertEqual(Orbifold.normalized_kaleidoscopes_equal([2, 3, 4, 5, 6], [4, 3, 5, 6, 2]), 0)

    def test_both_orientations_two_elements(self):
        self.assertEqual(Orbifold.normalized_kaleidoscopes_equal([2, 3], [2, 3]), 2)

    def test_both_orientations_two_runs(self):
        self.assertEqual(Orbifold.normalized_kaleidoscopes_equal([2, 2, 2, 3, 3], [2, 2, 3, 3, 2]), 2)


class TestOrbifoldNormalization(unittest.TestCase):
    def test_gyration_normalization(self):
        o1 = Orbifold(gyrations=[3, 1, 1, 2], kaleidoscopes=[], handles=0, crosscaps=0)
        self.assertEqual(o1.gyrations, [2, 3])

    def test_kaleidoscope_normalization(self):
        o1 = Orbifold(gyrations=[], kaleidoscopes=[[1, 1, 2, 3], [4, 1, 3]], handles=0, crosscaps=0)
        self.assertEqual(o1.kaleidoscopes, [[2, 3], [4, 3]])

    def test_odd_crosscaps_normalization(self):
        o1 = Orbifold(gyrations=[], kaleidoscopes=[], handles=3, crosscaps=5)
        self.assertEqual(o1.handles, 5)
        self.assertEqual(o1.crosscaps, 1)

    def test_even_crosscaps_normalization(self):
        o1 = Orbifold(gyrations=[], kaleidoscopes=[], handles=2, crosscaps=6)
        self.assertEqual(o1.handles, 4)
        self.assertEqual(o1.crosscaps, 2)

    def test_no_crosscaps_normalization(self):
        o1 = Orbifold(gyrations=[], kaleidoscopes=[], handles=3, crosscaps=0)
        self.assertEqual(o1.handles, 3)
        self.assertEqual(o1.crosscaps, 0)


class TestOrbifoldsEqual(unittest.TestCase):
    def test_different_handles(self):
        o1 = Orbifold(gyrations=[], kaleidoscopes=[], handles=1, crosscaps=1)
        o2 = Orbifold(gyrations=[], kaleidoscopes=[], handles=0, crosscaps=1)
        self.assertTrue(o1 != o2)

    def test_different_crosscaps(self):
        o1 = Orbifold(gyrations=[], kaleidoscopes=[], handles=1, crosscaps=1)
        o2 = Orbifold(gyrations=[], kaleidoscopes=[], handles=1, crosscaps=2)
        self.assertTrue(o1 != o2)

    def test_gyrations_equal(self):
        o1 = Orbifold(gyrations=[2, 4, 5, 1, 6], kaleidoscopes=[], handles=0, crosscaps=0)
        o2 = Orbifold(gyrations=[5, 6, 1, 1, 4, 2], kaleidoscopes=[], handles=0, crosscaps=0)
        self.assertTrue(o1 == o2)

    def test_gyrations_unequal(self):
        o1 = Orbifold(gyrations=[2, 4, 5, 1, 6, 6], kaleidoscopes=[], handles=0, crosscaps=0)
        o2 = Orbifold(gyrations=[5, 6, 1, 1, 4, 2], kaleidoscopes=[], handles=0, crosscaps=0)
        self.assertTrue(o1 != o2)

    def test_num_kaleidoscopes_different(self):
        o1 = Orbifold(gyrations=[], kaleidoscopes=[[]], handles=0, crosscaps=0)
        o2 = Orbifold(gyrations=[], kaleidoscopes=[[], []], handles=0, crosscaps=0)
        self.assertTrue(o1 != o2)

    def test_same_kaleidoscope_ordering(self):
        o1 = Orbifold(gyrations=[], kaleidoscopes=[[1, 2, 3], [4, 5, 6, 7]], handles=0, crosscaps=0)
        o2 = Orbifold(gyrations=[], kaleidoscopes=[[2, 3], [4, 5, 6, 7]], handles=0, crosscaps=0)
        self.assertTrue(o1 == o2)

    def test_totally_different_kaleidoscopes(self):
        o1 = Orbifold(gyrations=[], kaleidoscopes=[[2, 3], [4, 5, 6]], handles=0, crosscaps=0)
        o2 = Orbifold(gyrations=[], kaleidoscopes=[[4, 5], [7, 8, 9]], handles=0, crosscaps=0)
        self.assertTrue(o1 != o2)

    def test_same_oriented_kaleidoscopes(self):
        o1 = Orbifold(gyrations=[], kaleidoscopes=[[2, 3, 4], [5, 6, 7], [8, 9, 10, 11]], handles=0, crosscaps=0)
        o2 = Orbifold(gyrations=[], kaleidoscopes=[[9, 10, 11, 8], [3, 4, 2], [7, 5, 6]], handles=0, crosscaps=0)
        self.assertTrue(o1 == o2)

    def test_different_oriented_kaleidoscopes(self):
        o1 = Orbifold(gyrations=[], kaleidoscopes=[[2, 3, 4], [5, 6, 7], [8, 9, 10, 11]], handles=0, crosscaps=0)
        o2 = Orbifold(gyrations=[], kaleidoscopes=[[3, 4, 2], [9, 10, 11, 8], [7, 6, 5]], handles=0, crosscaps=0)
        self.assertTrue(o1 != o2)

    def test_different_oriented_kaleidoscopes_nonorientable(self):
        o1 = Orbifold(gyrations=[], kaleidoscopes=[[2, 3, 4], [5, 6, 7], [8, 9, 10, 11]], handles=0, crosscaps=1)
        o2 = Orbifold(gyrations=[], kaleidoscopes=[[3, 4, 2], [9, 10, 11, 8], [7, 6, 5]], handles=0, crosscaps=1)
        self.assertTrue(o1 == o2)

    def test_both_oriented_kaleidoscopes_1(self):
        o1 = Orbifold(gyrations=[], kaleidoscopes=[[2, 2, 3, 3], [5, 6, 7]], handles=0, crosscaps=0)
        o2 = Orbifold(gyrations=[], kaleidoscopes=[[2, 2, 3, 3], [7, 5, 6]], handles=0, crosscaps=0)
        self.assertTrue(o1 == o2)

    def test_both_oriented_kaleidoscopes_2(self):
        o1 = Orbifold(gyrations=[], kaleidoscopes=[[2, 2, 3, 3], [5, 6, 7]], handles=0, crosscaps=0)
        o2 = Orbifold(gyrations=[], kaleidoscopes=[[2, 2, 3, 3], [7, 6, 5]], handles=0, crosscaps=0)
        self.assertTrue(o1 == o2)

    def test_two_oppositely_oriented_kaleidoscopes_1(self):
        o1 = Orbifold(gyrations=[], kaleidoscopes=[[2, 3, 4], [5, 6, 7], [7, 6, 5]], handles=0, crosscaps=0)
        o2 = Orbifold(gyrations=[], kaleidoscopes=[[2, 3, 4], [6, 5, 7], [6, 7, 5]], handles=0, crosscaps=0)
        self.assertTrue(o1 == o2)

    def test_two_oppositely_oriented_kaleidoscopes_2(self):
        o1 = Orbifold(gyrations=[], kaleidoscopes=[[2, 3, 4], [5, 6, 7], [7, 6, 5]], handles=0, crosscaps=0)
        o2 = Orbifold(gyrations=[], kaleidoscopes=[[2, 3, 4], [6, 7, 5], [6, 5, 7]], handles=0, crosscaps=0)
        self.assertTrue(o1 == o2)

    def test_empty_kaleidoscope_1(self):
        o1 = Orbifold(gyrations=[], kaleidoscopes=[[], [2, 3, 4]], handles=0, crosscaps=0)
        o2 = Orbifold(gyrations=[], kaleidoscopes=[[1, 1], [3, 4, 2]], handles=0, crosscaps=0)
        self.assertTrue(o1 == o2)

    def test_empty_kaleidoscope_2(self):
        o1 = Orbifold(gyrations=[], kaleidoscopes=[[1, 1, 1], [2, 3, 4]], handles=0, crosscaps=0)
        o2 = Orbifold(gyrations=[], kaleidoscopes=[[], [4, 3, 2]], handles=0, crosscaps=0)
        self.assertTrue(o1 == o2)


class TestFaceCodes(unittest.TestCase):
    @staticmethod
    def get_normalization_map(fc):
        """

        :param fc: a parametrized face code as a list of pairs (number_of_sides, parameter)
        :return: a dictionary from the parameter in the face code to normalized parameters which are
           strings of the form [0], [1], ..., where [0] is the first distinct parameter in the original
           face code, [1] is the second distinct parameter, etc.

        For example:
          [(3,1), (3,1), (4,0)]
        returns the dictionary
          1 -> '[0]'
          2 -> '[1]'

        We use '[0]' and not just an integer 0, so that when the map is applied to orbifolds, we don't collide with
        fixed integers in the orbifold.
        """
        normalization_map = {}
        next_index = 0
        for (a, b) in fc:
            if b not in normalization_map:
                normalization_map[b] = '[%s]' % next_index
                next_index += 1
        return normalization_map

    @staticmethod
    def normalize_face_parameters(fc, normalization_map):
        return [(a, normalization_map[b]) for (a, b) in fc]

    @staticmethod
    def apply_normalization_map(normalization_map, parameter_list):
        return [a if type(a) == int or a == 'n' else normalization_map[a] for a in parameter_list]

    @staticmethod
    def normalize_orbifold_parameters(orb, normalization_map):
        """

        :param orb:  an Orbifold
        :param normalization_map:  a normalization map as returned by get_normalization_map on a face code
        :return: a new Orbifold with all parameters except 'n' in gyrations and kaleidoscopes normalized
          with the parameter names according to normalization_map.  We assume that 'n' is the vertex parameter
          so it wouldn't be in the normalization_map (which is on faces).
        """
        normalized_gyrations = TestFaceCodes.apply_normalization_map(normalization_map, orb.gyrations)
        normalized_kaleidoscopes = [TestFaceCodes.apply_normalization_map(normalization_map, k)
                                    for k in orb.kaleidoscopes]
        return Orbifold(normalized_gyrations, normalized_kaleidoscopes, orb.handles, orb.crosscaps)

    @staticmethod
    def parse_facecode(fc):
        """
            Parse a comma-separated string like 8b,a,8b,8b,8b into a face code
            [(8, 1), (1, 0), (8, 1), (8, 1), (8, 1)]
        """
        faces = fc.split(',')
        result = []
        for fg in faces:
            result.extend(TestFaceCodes.parse_face_group(fg))
        return result

    @staticmethod
    def parse_face(m):
        if m.group(1):
            sides = int(m.group(1))
        else:
            sides = 1
        variable = m.group(2)
        return sides, variable

    @staticmethod
    def parse_face_group(face_group):
        face_group = face_group.strip()
        m = re.match('^([0-9]+)?([a-z])$', face_group)  # match something like 7a or b
        if m is not None:
            sides, variable = TestFaceCodes.parse_face(m)
            count = 1
        else:
            m = re.match('\\(([0-9]+)?([a-z])\\)\\^([0-9]+)', face_group)  # match something like (7a)^3
            if m is None:
                raise ValueError('Invalid face %s' % face_group)
            sides, variable = TestFaceCodes.parse_face(m)
            count = int(m.group(3))
        return [(sides, variable)]*count

    @staticmethod
    def parse_orbifold(orb_string):
        """
            Parse an orbifold string like o234*2an*3bcxx into an Orbifold object.
            We assume that:
             - there are only single digit gyration and kaleidoscope numbers;
             - any symbol other than o, *, x, 0, 1, ..., 9 is a variable representing a parameter.
        """
        orb_regex = 'o+|\\*[^ox*]*|x+|[^ox*]+'
        if not re.match('^(%s)+$' % orb_regex, orb_string):
            raise ValueError('Invalid orbifold string %s' % orb_string)
        feature_strings = re.findall(orb_regex, orb_string)
        gyrations = []
        kaleidoscopes = []
        handles = 0
        crosscaps = 0
        digits = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
        for fs in feature_strings:
            if fs[0] == 'o':
                handles = len(fs)
            elif fs[0] == 'x':
                crosscaps = len(fs)
            elif fs[0] == '*':
                kaleidoscopes.append([int(x) if x in digits else x for x in list(fs[1:])])
            else:
                gyrations = list([int(x) if x in digits else x for x in list(fs)])

        return Orbifold(gyrations, kaleidoscopes, handles, crosscaps)


def generate_permutation_symbol_test(symbol_string, correct_facecode_string, correct_orbifold_string=None):
    def test(self):  # self is of type TestFaceCode
        ps = PermutationSymbol(symbol_string)
        correct_facecode = self.parse_facecode(correct_facecode_string)

        normalization_map_ps = self.get_normalization_map(ps.face_code)
        normalization_map_correct = self.get_normalization_map(correct_facecode)

        normalized_ps_face_indices = self.normalize_face_parameters(ps.face_code, normalization_map_ps)
        normalized_correct_face_indices = self.normalize_face_parameters(correct_facecode, normalization_map_correct)

        self.assertEqual(normalized_ps_face_indices, normalized_correct_face_indices)

        if correct_orbifold_string:
            correct_orbifold = self.parse_orbifold(correct_orbifold_string)
            orbifold_normalization_map_ps = {}
            for k, v in normalization_map_ps.items():  # because (3,1) in the face code becomes '[1]' in orbifold
                orbifold_normalization_map_ps['[%s]' % k] = v
            normalized_orbifold_ps = self.normalize_orbifold_parameters(ps.orbifold, orbifold_normalization_map_ps)
            normalized_orbifold_correct = self.normalize_orbifold_parameters(correct_orbifold, normalization_map_correct)
            self.assertEqual(normalized_orbifold_ps, normalized_orbifold_correct)
    return test


if __name__ == '__main__':
    f = open('test_cases.txt')
    count = 1
    for testcase_string in f:
        testcase_string = testcase_string.strip()
        try:
            testcase_string = testcase_string[:testcase_string.index('#')]
        except ValueError:
            pass
        if testcase_string == '':
            continue
        test_name = 'test_example_%s' % count
        testcase = testcase_string.split(None, 3)
        test = generate_permutation_symbol_test(testcase[0], testcase[1], testcase[2] if len(testcase) == 3 else None)
        setattr(TestFaceCodes, test_name, test)
        count += 1
    f.close()
    unittest.main()
