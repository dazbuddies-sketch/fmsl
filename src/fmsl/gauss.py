"""Exact Gaussian rationals a + b*i with a, b in Q (the field Q(i))."""
from fractions import Fraction as Fr
from numbers import Rational

_NUM = (int, Rational)


def _f(v):
    return v if isinstance(v, Fr) else Fr(v)


class G:
    __slots__ = ("re", "im")

    def __init__(s, re=0, im=0):
        if isinstance(re, G):
            s.re, s.im = re.re, re.im
        else:
            s.re, s.im = _f(re), _f(im)

    # ---- arithmetic ----
    def __add__(s, o):
        if isinstance(o, G): return G(s.re + o.re, s.im + o.im)
        if not isinstance(o, _NUM): return NotImplemented
        return G(s.re + o, s.im)
    __radd__ = __add__

    def __sub__(s, o):
        if isinstance(o, G): return G(s.re - o.re, s.im - o.im)
        if not isinstance(o, _NUM): return NotImplemented
        return G(s.re - o, s.im)

    def __rsub__(s, o):
        if not isinstance(o, _NUM): return NotImplemented
        return G(o - s.re, -s.im)
    def __neg__(s): return G(-s.re, -s.im)

    def __mul__(s, o):
        if isinstance(o, G):
            if not o.im: return G(s.re * o.re, s.im * o.re)
            if not s.im: return G(s.re * o.re, s.re * o.im)
            return G(s.re * o.re - s.im * o.im, s.re * o.im + s.im * o.re)
        if not isinstance(o, _NUM): return NotImplemented
        return G(s.re * o, s.im * o)
    __rmul__ = __mul__

    def __truediv__(s, o):
        if isinstance(o, G):
            if not o.im: return G(s.re / o.re, s.im / o.re)
            d = o.re * o.re + o.im * o.im
            return G((s.re * o.re + s.im * o.im) / d, (s.im * o.re - s.re * o.im) / d)
        if not isinstance(o, _NUM): return NotImplemented
        return G(s.re / o, s.im / o)

    def __rtruediv__(s, o):
        if not isinstance(o, _NUM): return NotImplemented
        return G(o) / s

    def __pow__(s, n):
        if n < 0: return G(1) / (s ** -n)
        out, b = G(1), s
        while n:
            if n & 1: out = out * b
            b = b * b; n >>= 1
        return out

    # ---- comparison ----
    def __eq__(s, o):
        if isinstance(o, G): return s.re == o.re and s.im == o.im
        if isinstance(o, tuple): return (s.re, s.im) == o
        if not isinstance(o, _NUM): return NotImplemented
        return s.im == 0 and s.re == o

    def __ne__(s, o):
        r = s.__eq__(o)
        return r if r is NotImplemented else not r
    def __hash__(s): return hash((s.re, s.im))
    def __bool__(s): return bool(s.re) or bool(s.im)

    def _real(s, op):
        if s.im: raise TypeError(f"ordering ({op}) is undefined for complex value {s}")
        return s.re

    def __gt__(s, o): return s._real(">") > (o.re if isinstance(o, G) else o)
    def __lt__(s, o): return s._real("<") < (o.re if isinstance(o, G) else o)
    def __ge__(s, o): return s._real(">=") >= (o.re if isinstance(o, G) else o)
    def __le__(s, o): return s._real("<=") <= (o.re if isinstance(o, G) else o)

    # ---- helpers ----
    def is_real(s): return s.im == 0
    def conj(s): return G(s.re, -s.im)
    def abs_upper(s): return abs(s.re) + abs(s.im)          # >= |z|
    def abs_lower(s): return max(abs(s.re), abs(s.im))      # <= |z|
    def __float__(s): return float(s._real("float"))
    def __complex__(s): return complex(float(s.re), float(s.im))
    def pair(s): return (s.re, s.im)

    def __repr__(s):
        if not s.im: return str(s.re)
        if not s.re: return f"{s.im}i" if s.im not in (1, -1) else ("i" if s.im == 1 else "-i")
        sign = "+" if s.im > 0 else "-"
        im = abs(s.im)
        return f"{s.re}{sign}{'' if im == 1 else im}i"
    __str__ = __repr__


I = G(0, 1)
