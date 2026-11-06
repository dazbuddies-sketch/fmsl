from fractions import Fraction as Fr
from fmsl.jelonek_n import analyse, arc_point, _polys
from fmsl.witness import _ev

V = ["x", "y", "z"]


def _in_outer(r, p):
    return all(any(_ev(L["lead"], p) == 0 for L in basis) for basis in r["outer"])


def test_plane():
    r = analyse(["x", "x*y", "z"], V)
    assert r["inner_points"] and all(p[0] == 0 for p in r["inner_points"])
    assert r["surface_lines"][0] > 100                 # S_F is the plane u1 = 0
    assert all(_in_outer(r, p) for p in r["inner_points"])


def test_automorphism_is_proper():
    r = analyse(["x", "y", "z+x*y"], V)
    assert not r["inner_points"] and r["surface_lines"][0] == 0


def test_codimension_two_line():
    m = ["x", "y", "(x^2+y^2)*z^2+z"]                    # S_F is the z-axis
    r = analyse(m, V)
    assert r["inner_points"] and all(p[0] == 0 and p[1] == 0 for p in r["inner_points"])
    assert r["surface_lines"][0] == 0
    assert all(_in_outer(r, p) for p in r["inner_points"])
    F = _polys(m, V)
    for pt, arc in r["inner_arcs"][:20]:               # each point is approached along its arc
        x = arc_point(arc, Fr(1, 10 ** 30))
        assert max(abs(c) for c in x) > 10 ** 9
        assert all(abs(_ev(P, x) - v) < Fr(1, 10 ** 6) for P, v in zip(F, pt))
