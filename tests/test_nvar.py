from math import inf
from fmsl.nvar import limit_n

V4 = ["x", "y", "z", "w"]


def test_range_enclosed():
    r = limit_n("x*w/(x^2+y^2+z^2+w^2)", V4)
    assert r["proven"] is True and r["limit_exists"] is False
    (a, b), = r["components"]
    assert abs(a + 0.5) < 1e-6 and abs(b - 0.5) < 1e-6
    lo, hi = r["enclosure"][0][0], r["enclosure"][1][1]
    assert lo <= -0.5 and hi >= 0.5


def test_limits_that_exist():
    assert limit_n("x^3/(x^2+y^2+z^2+w^2)", V4)["limit"] == "0"
    assert limit_n("(x^2+y^2+z^2+w^2)/(x^2+y^2+z^2+w^2+x^3)", V4)["limit"] == "1"
    assert limit_n("1/(x^2+y^4+z^2+w^2)", V4)["limit"] == "+oo"


def test_whole_circle():
    assert limit_n("x/(x^2+y^2+z^2+w^2)", V4)["components"] == [(-inf, inf)]


def test_falls_back_to_inner_bound():
    r = limit_n("x*y/(x^2+y^2-z*w)", V4)          # denominator vanishes on a cone: not covered
    assert r["method"] == "witness curves"
