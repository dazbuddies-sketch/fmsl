"""FMSL for SymPy users: exact multivariable limits with SymPy expressions in and SymPy objects out.

SymPy's `limit` handles one variable. This bridge adds two and three variables (rational
functions), with the complete set of limit values when the limit does not exist:

    >>> from sympy import symbols, Rational
    >>> from fmsl.sympy_bridge import limit_set, multivariate_limit
    >>> x, y = symbols("x y", real=True)
    >>> limit_set(x**3*y/(x**6 + y**2), (x, y), (0, 0))
    Interval(-1/2, 1/2)
    >>> multivariate_limit(x**2*y**2/(x**2 + y**2), (x, y), (0, 0))
    0
    >>> limit_set((x**2 + y**2)/(x*y), (x, y))           # infinities are included as points
    Union(FiniteSet(-oo, oo), Interval(-oo, -2), Interval(2, oo))

Endpoints are exact: rationals, or real algebraic numbers returned as SymPy expressions
(radicals when SymPy finds them, otherwise CRootOf). The point may be any point with rational
coordinates. `pip install "fmsl[sympy]"`.
"""
from fractions import Fraction


class NoLimit(ValueError):
    """the limit does not exist; .values is the SymPy set of all limit values"""
    def __init__(self, values):
        super().__init__(f"the limit does not exist; the set of limit values is {values}")
        self.values = values


def _sympy():
    try:
        import sympy
    except ImportError:
        raise ImportError('fmsl.sympy_bridge needs SymPy: pip install "fmsl[sympy]"')
    return sympy


def _to_fmsl_text(expr, variables, point):
    sp = _sympy()
    variables = list(variables)
    if not 1 <= len(variables) <= 3:
        raise ValueError("one, two or three variables")
    point = [0] * len(variables) if point is None else list(point)
    names = ["x", "y", "z"][:len(variables)]
    tmp = [sp.Symbol(f"_fmsl_{n}", real=True) for n in names]
    sub = {}
    for v, t, p in zip(variables, tmp, point):
        p = sp.nsimplify(p)
        if not p.is_Rational:
            raise ValueError(f"the point must have rational coordinates (got {p})")
        sub[v] = t + p                                   # move the point to the origin
    e = sp.together(sp.sympify(expr).subs(sub, simultaneous=True))
    num, den = sp.fraction(e)
    for part in (num, den):
        if not part.is_polynomial(*tmp):
            raise ValueError("two or three variables: rational functions only")
    text = f"({sp.expand(num)})/({sp.expand(den)})"
    for t, n in zip(tmp, names):
        text = text.replace(str(t), n)
    return text.replace("**", "^"), names


def _to_sympy_number(v):
    sp = _sympy()
    from .cluster import _Inf
    if isinstance(v, _Inf): return sp.oo if v.sgn > 0 else -sp.oo
    if isinstance(v, Fraction): return sp.Rational(v.numerator, v.denominator)
    if isinstance(v, (int,)): return sp.Integer(v)
    if isinstance(v, float):
        if v in (float("inf"), float("-inf")): return sp.oo if v > 0 else -sp.oo
        return sp.nsimplify(v)
    # a real algebraic number g(alpha), alpha the root of F.p isolated in (F.lo, F.hi]
    from .algebraic import count_roots, root_bound
    F = v.F
    p = [Fraction(c) for c in F.p]
    t = sp.Symbol("t")
    P = sp.Poly(sum(sp.Rational(c.numerator, c.denominator) * t ** i for i, c in enumerate(p)), t)
    index = count_roots(p, -root_bound(p) - 1, F.lo)     # real roots of p below alpha
    alpha = sp.CRootOf(P, index)
    c = F.reduce(v.c)
    val = sum(sp.Rational(Fraction(ci).numerator, Fraction(ci).denominator) * alpha ** i for i, ci in enumerate(c))
    try:
        s = sp.nsimplify(sp.radsimp(sp.expand(val)))
        if abs(float(s) - float(v)) < 1e-12 * max(1.0, abs(float(v))): return s
    except (TypeError, ValueError, NotImplementedError):
        pass
    return val


def _components(expr, variables, point):
    from .verify import limit_at_origin
    text, names = _to_fmsl_text(expr, variables, point)
    if len(names) == 1:
        raise ValueError("one variable: use sympy.limit (both one-sided limits)")
    if len(names) == 2:
        from .multivar import profile
        cl = profile(text)["cluster"]
        proven = bool(cl.get("exact")) and bool(cl.get("complete", True))
        comps = cl["signed_components"]
    else:
        from .space import profile3
        r = profile3(text)
        proven = bool(r.get("exact", True)) and bool(r.get("complete"))
        comps = r["signed_components"]
    if not proven:
        raise ArithmeticError("FMSL could not prove the set of limit values for this expression; nothing is claimed")
    return comps


def limit_set(expr, variables, point=None):
    """The complete set of limit values of `expr` as `variables` -> `point` (default: the origin),
       as a SymPy set. +-oo are included as points when f approaches them."""
    sp = _sympy()
    parts = []
    for lo, hi in _components(expr, variables, point):
        a, b = _to_sympy_number(lo), _to_sympy_number(hi)
        if a == b:
            parts.append(sp.FiniteSet(a))
        else:
            parts.append(sp.Interval(a, b))
            for e in (a, b):
                if e in (sp.oo, -sp.oo): parts.append(sp.FiniteSet(e))
    return sp.Union(*parts)


def multivariate_limit(expr, variables, point=None):
    """The limit of `expr` as `variables` -> `point`, proven. Raises NoLimit (with the set of all
       limit values in .values) when the limit does not exist."""
    sp = _sympy()
    comps = _components(expr, variables, point)
    from .cluster import cmp
    if len(comps) == 1 and cmp(comps[0][0], comps[0][1]) == 0:
        return _to_sympy_number(comps[0][0])            # a real number, or +-oo (diverges)
    raise NoLimit(limit_set(expr, variables, point))
