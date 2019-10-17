from .bond import Bond

class Angle(object):
    """ An angle between three Atom objects.

        Currently, there is no reference to the actual Atom objects but instead
        they are indexed in the parent molecule class.

        An angle between three atoms:

        A
        |
        |
        B----C

        is defined in the Angle class as

        >>> a = Angle(B.get_idx(), A.get_idx(), C.get_idx())

        where B is the vertex of the angle.
    """
    def __init__(self, vertex, id1, id2):
        self._vertex = vertex
        self._id1 = id1
        self._id2 = id2
        s = set([vertex, id1, id2])
        if len(s) < 3:
            raise ValueError("Indices cannot be the same. Initialized as {}".format(repr(self)))
        if vertex < 0:
            raise ValueError("Argument vertex supplied negative value.")
        if id1 < 0:
            raise ValueError("Argument id1 supplied negative value.")
        if id2 < 0:
            raise ValueError("Argument id2 supplied negative value.")

    def __getitem__(self, key):
        if not isinstance(key, int):
            raise TypeError("provided key must be integer.")

        if key < 0:
            raise IndexError("negative indices not supported")

        if key > 2:
            raise IndexError("Angle only has three indices.")

        # kinda bad that Angle is not derived from a simple list
        if key == 0: return self._vertex
        if key == 1: return self._id1
        if key == 2: return self._id2

    def __eq__(self, other):
        if self._vertex != other._vertex:
            return False

        # make bonds extending from the vertex
        b1_self = Bond(self._vertex, self._id1)
        b2_self = Bond(self._vertex, self._id2)
        b1_other = Bond(other._vertex, other._id1)
        b2_other = Bond(other._vertex, other._id2)

        # and compare bonds
        value1 = (b1_self == b1_other) or (b1_self == b2_other)
        value2 = (b2_self == b1_other) or (b2_self == b2_other)

        return value1 and value2

    def __repr__(self):
        return("Angle({0:d}, {1:d}, {2:d})".format(self._vertex, self._id1, self._id2))
