"""FMSL: certified analysis and safe code generation for formulas with 0/0 traps."""
from .engine import V, x, Exp, Sin, Cos, Essential
from .gauss import G, I
from .compiler import run, compile_text, analyze
from .siletz import proven_series_error

__version__ = "0.28.0"
# Matrix functions (requires numpy):  from fmsl import matrix
# Two variables / exact arcs / algebraic numbers / iteration:
#   from fmsl import multivar, cluster, algebraic, bifurcation
# Three variables / exact two-variable ranges:  from fmsl import space, plane
__all__ = ["V", "x", "Exp", "Sin", "Cos", "Essential", "G", "I",
           "run", "compile_text", "analyze", "proven_series_error"]
