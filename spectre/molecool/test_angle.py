import pytest

from angle import Angle


def test_angle():
    with pytest.raises(ValueError):
        Angle(0, 0, 1)

    with pytest.raises(ValueError):
        Angle(-1, 1, 2)

    with pytest.raises(ValueError):
        Angle(1, -1, 2)

    with pytest.raises(ValueError):
        Angle(1, 0, -2)

    a1 = Angle(0, 1, 2)
    a2 = Angle(0, 2, 1)
    assert a1 == a2

    a3 = Angle(1, 2, 0)
    assert a3 != a1
    assert a3 != a2

    a4 = Angle(0, 1, 3)
    assert a4 != a1

    # a funky setup to trigger the first return statement in __eq__ of angle
    a1 = Angle(0, 1, 2)
    a2 = Angle(1, 2, 3)
    the_same = a1 == a2
    assert not the_same

def test_angle_indices():
    a = Angle(1, 4, 2)
    assert a[0] == 1
    assert a[1] == 4
    assert a[2] == 2

    with pytest.raises(IndexError):
        a[-1]

    with pytest.raises(IndexError):
        a[3]

    with pytest.raises(TypeError):
        a[-1.0]
