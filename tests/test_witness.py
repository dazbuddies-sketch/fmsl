from fractions import Fraction as Fr
from math import inf
from fmsl.witness import inner_bound, parse_n, _ev


def test_four_variables_known_answers():
    r = inner_bound("x*w/(x^2+y^2+z^2+w^2)", ["x", "y", "z", "w"])
    assert r["components"] == [(-0.5, 0.5)] and r["limit_exists"] is False
    r = inner_bound("x*y*z*w/(x^4+y^4+z^4+w^4)", ["x", "y", "z", "w"])
    (lo, hi), = r["components"]
    assert -0.25 <= lo < -0.2499 and 0.2499 < hi <= 0.25


def test_whole_circle_is_exact():
    r = inner_bound("x/y", ["x", "y"])
    assert r["whole_circle"] and r["components"] == [(-inf, inf)]
    # the Part 4 stress-test case: the denominator's zero hugs the z-axis at order z^3
    f = "((-2)*x*y*(1+z)+(-2)*z*x^2*x)/((-2)*x^2*z^2+(-1)*x*y*(2-z)+(-1)*z^3*x)"
    assert inner_bound(f, ["x", "y", "z"])["whole_circle"]


def test_limit_exists_is_never_claimed():
    r = inner_bound("(x^2+y^2+z^2+w^2)/(1+x)", ["x", "y", "z", "w"])
    assert r["limit_exists"] is None and r["components"] == [(0.0, 0.0)]


def test_witness_curves():
    names = ["x", "y", "z", "w"]
    f = "(x^2*y-z*w^2)/(x^4+y^2+z^2+w^4)"
    P, Q = parse_n(f, names)
    for val, w, c in inner_bound(f, names)["witnesses"]:
        if val in ("segment", "inf", "-inf"): continue
        t = Fr(1, 10 ** 12)
        x = tuple(ci * t ** wi for ci, wi in zip(c, w))
        assert abs(_ev(P, x) / _ev(Q, x) - Fr(val)) < Fr(1, 10 ** 6)


def test_zoom_reaches_two_term_curves():
    from fmsl.witness import curve_point
    f = "x*y*z/(x^2+z*y^2)"                         # exact answer: the whole circle
    r = inner_bound(f, ["x", "y", "z"])
    assert r["components"][0][0] == -inf and r["components"][-1][1] == inf
    P, Q = parse_n(f, ["x", "y", "z"])
    zoomed = [w for w in r["witnesses"] if len(w) == 4 and w[0] not in ("segment", "inf", "-inf")]
    assert zoomed
    for wt in zoomed[:50]:
        x = curve_point(wt, Fr(1, 10 ** 60))
        assert abs(_ev(P, x) / _ev(Q, x) - Fr(wt[0])) < Fr(1, 10 ** 6) * (1 + abs(Fr(wt[0])))
