"""FMSL as an MCP server: lets an AI assistant compute proven limits and check its own answers
instead of guessing.

    pip install "fmsl[mcp]"
    fmsl-mcp                      # speaks MCP over stdio

Claude Desktop / Claude Code configuration:
    {"mcpServers": {"fmsl": {"command": "fmsl-mcp"}}}

Tools:
    fmsl_limit(expression, variables?)          the complete set of limit values at the origin, proven
    fmsl_check_limit(expression, claim, ...)    is a claimed limit or set right? what is missing or extra?
    fmsl_rate(expression, variables?)           how fast the limit is approached: |f - L| <= C r^theta
    fmsl_extremum(expression, kind)             exact global min/max over the plane; attained or not
    fmsl_nonproperness(f, g)                    the Jelonek set of the polynomial map (f, g) of the plane
"""
import multiprocessing as mp
from typing import List, Optional

from . import __version__
from .verify import check_claim, limit_at, limit_at_origin, rate_at_origin

TIMEOUT = 60          # seconds per call; hard cases are reported, not left hanging

INSTRUCTIONS = (
    "FMSL computes, with proof, the complete set of values a function approaches at the origin "
    "(the limit if it exists, otherwise every limit value). Use it whenever an answer depends on a "
    "limit at a 0/0 point, a multivariable limit, or whether a range of values is ever approached: "
    "these are easy to get wrong by testing a few paths. One variable: expressions with + - * / ^, "
    "sin, cos, exp, abs, sign. Two or three variables: rational functions. Move a point (a, b) to the "
    "origin first by substituting x -> x + a.")


def _worker(q, fn, args):
    try:
        q.put(("ok", fn(*args)))
    except Exception as exc:                            # report the engine's message to the caller
        q.put(("error", f"{type(exc).__name__}: {exc}"))


def _run(fn, *args):
    ctx = mp.get_context("spawn")
    q = ctx.Queue()
    p = ctx.Process(target=_worker, args=(q, fn, args), daemon=True)
    p.start(); p.join(TIMEOUT)
    if p.is_alive():
        p.terminate()
        return {"error": f"no answer within {TIMEOUT} s; FMSL claims nothing for this expression"}
    status, val = q.get() if not q.empty() else ("error", "the engine stopped without an answer")
    return val if status == "ok" else {"error": val}


def _limit(expression, variables, point, side):
    if point is None and side is None: return limit_at_origin(expression, variables)
    from fractions import Fraction
    pts = [Fraction(str(p)) for p in (point or [0])]
    return limit_at(expression, pts, variables, side)


def fmsl_limit(expression: str, variables: Optional[List[str]] = None,
               point: Optional[List[str]] = None, side: Optional[str] = None) -> dict:
    """The complete set of values `expression` approaches as its variables tend to `point`
    (default: the origin; rational coordinates such as ["2"] or ["1/3", "0"]), proven.
    One variable: side "-" (from below) or "+" (from above) also gives the one-sided limit.

    expression: e.g. "x^3*y/(x^6+y^2)", "sin(x)/x", "(x*y+y*z+z*x)/(x^2+y^2+z^2)", "x*y/sqrt(x^2+y^2)" (use ^ for powers; a Euclidean
    distance to the point, such as sqrt(x^2+y^2+z^2), and sin, cos, exp, log, ... are supported
    in two and three variables, e.g. "sin(x*y)/(x^2+y^2)").
    variables: optional, e.g. ["x", "y"]; by default every name in the expression except functions.
    Returns limit_exists, limit, set (e.g. "[-1/2, 1/2]" or "[-oo, -2] u [2, +oo]"), proven, an
    explanation, and in two variables witnesses: curves along which the end values are approached.
    """
    return _run(_limit, expression, variables, point, side)


def fmsl_check_limit(expression: str, claim: str, variables: Optional[List[str]] = None) -> dict:
    """Check a claimed answer for the limit of `expression` at the origin against FMSL's proven answer.

    claim: a value ("0", "1/2", "+oo", "LIMIT 3") or a set ("[-1/2, 1/2]", "[-oo, -2] u [2, +oo]",
    "{0} u [1, 2]"). Returns correct, the proven set, missing (values approached but left out of the
    claim), extra (claimed values that are never approached) and a one-line verdict.
    """
    return _run(check_claim, expression, claim, variables)


def fmsl_rate(expression: str, variables: Optional[List[str]] = None) -> dict:
    """How fast `expression` approaches its limit L at the origin: the best exponent theta with
    |f - L| <= C * r^theta near the origin (r = distance to the origin), proven and attained along
    the reported curve. One or two variables; for example "x^2*y/(x^2+y^2)" gives theta = 1.
    """
    return _run(rate_at_origin, expression, variables)


def _extremum(expression, kind):
    from .optimize import global_min, global_max
    from .verify import _variables, _rename
    vs = _variables(expression)
    if len(vs) != 2: raise ValueError("global minimum/maximum: functions of exactly two variables")
    r = (global_max if kind == "max" else global_min)(_rename(expression, vs, ["x", "y"]))
    r.pop("value_exact", None)
    back = dict(zip(["x", "y"], vs))
    r["text"] = r["text"].replace("(x, y)", f"({vs[0]}, {vs[1]})")
    r["variables"] = vs
    return r


def fmsl_extremum(expression: str, kind: str = "min") -> dict:
    """The exact global minimum (kind="min") or maximum (kind="max") of a function of two variables
    over the whole plane, and whether it is ATTAINED (and where) or only approached at infinity.
    Polynomials, or rational functions whose denominator never vanishes. Example: "x^2+(x*y-1)^2"
    has infimum 0, which is not attained.
    """
    return _run(_extremum, expression, kind)


def _jelonek(f, g):
    from .jelonek import jelonek_set
    r = jelonek_set(f, g)
    return {"text": r["text"], "empty": r["empty"], "complete": r["complete"]}


def fmsl_nonproperness(f: str, g: str) -> dict:
    """The exact non-properness (Jelonek) set of the polynomial map F = (f, g) of the real plane: the
    points that F approaches along sequences going to infinity. F is proper exactly off this set; a
    map with nonvanishing Jacobian is a bijection iff the set is empty (it is then proper). Polynomials in x and y.
    """
    return _run(_jelonek, f, g)


def build_server():
    try:
        from mcp.server.mcpserver import MCPServer as Server      # mcp >= 2
    except ImportError:
        from mcp.server.fastmcp import FastMCP as Server          # mcp 1.x
    server = Server("fmsl", instructions=INSTRUCTIONS)
    server.tool()(fmsl_limit)
    server.tool()(fmsl_check_limit)
    server.tool()(fmsl_rate)
    server.tool()(fmsl_extremum)
    server.tool()(fmsl_nonproperness)
    return server


def main():
    try:
        server = build_server()
    except ImportError:
        raise SystemExit('The MCP server needs the mcp package: pip install "fmsl[mcp]"')
    server.run()


if __name__ == "__main__":
    main()
