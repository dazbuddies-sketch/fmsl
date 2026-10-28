"""Horizon or true singularity? FMSL on the formulas of black holes (python examples/black_holes.py).

Units with G = c = 1 and mass M = 1, so a non-spinning (Schwarzschild) black hole has its horizon
at r = 2. The metric component g_rr = 1/(1 - 2/r) blows up there, and flips sign across it. Is
space itself singular there? A curvature invariant answers that: the Kretschmann scalar
K = 48 M^2 / r^6 does not depend on the coordinates, and it is finite at r = 2 (the horizon is a
coordinate singularity) but infinite at r = 0 (the true singularity).

For a spinning (Kerr) black hole with M = 1 and spin a = 1, the Kretschmann scalar is
    K = 48 (r^6 - 15 r^4 c^2 + 15 r^2 c^4 - c^6) / (r^2 + c^2)^6,     c = cos(theta),
and the singularity is the ring r = 0, theta = pi/2, that is (r, c) = (0, 0). There the limit of
K depends on the direction of approach: its fluid state is the whole extended line.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from fmsl.verify import limit_at, limit_at_origin


def table():
    rows = [
        ("Schwarzschild g_rr = 1/(1-2/r), horizon r = 2, from outside", limit_at("1/(1-2/r)", 2, side="+")["one_sided_limit"]),
        ("Schwarzschild g_rr = 1/(1-2/r), horizon r = 2, from inside", limit_at("1/(1-2/r)", 2, side="-")["one_sided_limit"]),
        ("Schwarzschild g_tt = -(1-2/r), horizon r = 2", limit_at("-(1-2/r)", 2)["set"]),
        ("Kretschmann 48/r^6, horizon r = 2", limit_at("48/r^6", 2)["set"]),
        ("Kretschmann 48/r^6, centre r = 0, from outside", limit_at("48/r^6", 0, side="+")["one_sided_limit"]),
        ("Kerr Kretschmann at the ring (r, cos theta) -> (0, 0)",
         limit_at_origin("48*(r^6-15*r^4*c^2+15*r^2*c^4-c^6)/(r^2+c^2)^6", ["r", "c"])["set"]),
    ]
    return rows


if __name__ == "__main__":
    for name, value in table():
        print(f"{name:66s} {value}")
