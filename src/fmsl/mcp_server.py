"""FMSL as an MCP server: lets an AI assistant compute proven limits and check its own answers
instead of guessing.

    pip install "fmsl[mcp]"
    fmsl-mcp                      # speaks MCP over stdio

Claude Desktop / Claude Code configuration:
    {"mcpServers": {"fmsl": {"command": "fmsl-mcp"}}}

Tools:
    fmsl_limit(expression, variables?)          the complete set of limit values at the origin, proven
    fmsl_check_limit(expression, claim, ...)    is a claimed limit or set right? what is missing or extra?
"""
import multiprocessing as mp
from typing import List, Optional

from . import __version__
from .verify import check_claim, limit_at_origin

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


def fmsl_limit(expression: str, variables: Optional[List[str]] = None) -> dict:
    """The complete set of values `expression` approaches as all its variables tend to 0, proven.

    expression: e.g. "x^3*y/(x^6+y^2)", "sin(x)/x", "(x*y+y*z+z*x)/(x^2+y^2+z^2)" (use ^ for powers).
    variables: optional, e.g. ["x", "y"]; by default every name in the expression except functions.
    Returns limit_exists, limit, set (e.g. "[-1/2, 1/2]" or "[-oo, -2] u [2, +oo]"), proven, an
    explanation, and in two variables witnesses: curves along which the end values are approached.
    """
    return _run(limit_at_origin, expression, variables)


def fmsl_check_limit(expression: str, claim: str, variables: Optional[List[str]] = None) -> dict:
    """Check a claimed answer for the limit of `expression` at the origin against FMSL's proven answer.

    claim: a value ("0", "1/2", "+oo", "LIMIT 3") or a set ("[-1/2, 1/2]", "[-oo, -2] u [2, +oo]",
    "{0} u [1, 2]"). Returns correct, the proven set, missing (values approached but left out of the
    claim), extra (claimed values that are never approached) and a one-line verdict.
    """
    return _run(check_claim, expression, claim, variables)


def build_server():
    try:
        from mcp.server.mcpserver import MCPServer as Server      # mcp >= 2
    except ImportError:
        from mcp.server.fastmcp import FastMCP as Server          # mcp 1.x
    server = Server("fmsl", instructions=INSTRUCTIONS)
    server.tool()(fmsl_limit)
    server.tool()(fmsl_check_limit)
    return server


def main():
    try:
        server = build_server()
    except ImportError:
        raise SystemExit('The MCP server needs the mcp package: pip install "fmsl[mcp]"')
    server.run()


if __name__ == "__main__":
    main()
