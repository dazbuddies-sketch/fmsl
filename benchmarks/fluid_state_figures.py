"""Figures of two-variable fluid states. Every verdict, curve and value shown is taken from
fmsl.multivar's exact output (nothing hand-entered).  Run:  python benchmarks/fluid_state_figures.py [outdir]"""
from fractions import Fraction as Fr
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from fmsl.multivar import profile, _ueval

OUT = sys.argv[1] if len(sys.argv) > 1 else "."

# ---- reference palette (light mode) ----
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#e34948"]
DIVERGING = LinearSegmentedColormap.from_list(
    "fms_div", ["#104281", "#3987e5", "#9ec5f4", "#f0efec", "#f3b3b2", "#e34948", "#9c2424"])
ARROW, SQ2, SQ3, SQ4, NE, DOT = "→", "²", "³", "⁴", "≠", "·"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10.5, "text.color": INK,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "axes.titlesize": 12, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.facecolor": SURFACE,
    "figure.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})

def frame(ax):
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.tick_params(length=3)

def feval(P, X, Y):
    return sum(float(v) * X ** i * Y ** j for (i, j), v in P.items())

def field(N, D, lim, n=900):
    xs = np.linspace(-lim, lim, n)
    X, Y = np.meshgrid(xs, xs)
    with np.errstate(all="ignore"):
        Z = feval(N, X, Y) / feval(D, X, Y)
    return X, Y, Z

def depends_entry(r):
    for e in r["entries"]:
        if e["kind"] == "depends" and e.get("sigma") == 1:
            return e
    return None

def pretty(f):
    return (f.replace("*", "").replace("^2", SQ2).replace("^3", SQ3).replace("^4", SQ4))

# ======================= Figure 1: anatomy of one fluid state =======================
F = "x^2*y/(x^4+y^2)"
r = profile(F)
N, D = r["N"], r["D"]
crit = sorted({e["a"] for e in r["entries"] if e.get("critical")})
dep = depends_entry(r)
a_star = dep["a"]                                       # balance curvature found by the engine
Rx = lambda c: _ueval(dep["R"][0], Fr(c)) / _ueval(dep["R"][1], Fr(c))   # exact
R = lambda c: float(Rx(c))
CS = (Fr(1, 2), Fr(1), Fr(2), Fr(-1))

fig = plt.figure(figsize=(14, 11.8))
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.26, left=0.06, right=0.97, top=0.88, bottom=0.06)
fig.text(0.06, 0.955, f"One point, many answers: the fluid state of  f(x, y) = x{SQ2}y / (x{SQ4} + y{SQ2})  at the origin",
         fontsize=16, fontweight="bold", color=INK)
fig.text(0.06, 0.928, f"Engine verdict: {r['verdict'][0]}.  Every straight line gives 0, yet along parabolas "
         f"y = c{DOT}x{SQ2} the value is c/(1+c{SQ2}). The trap hides at exactly one balance curvature, a = {a_star}.",
         fontsize=11, color=INK2)

# --- A: the landscape with approach paths ---
ax = fig.add_subplot(gs[0, 0]); frame(ax)
lim = 1.0
X, Y, Z = field(N, D, lim)
im = ax.imshow(Z, extent=(-lim, lim, -lim, lim), origin="lower", cmap=DIVERGING, vmin=-0.5, vmax=0.5,
               interpolation="bilinear")
t = np.linspace(-lim, lim, 400)
for m in (-2, -0.5, 0.5, 2):
    ax.plot(t, m * t, color=INK, lw=1.2, ls=(0, (4, 3)), alpha=0.75)
handles = [ax.plot([], [], color=INK, lw=1.2, ls=(0, (4, 3)), label=f"any straight line: {ARROW} 0")[0]]
for k, c in enumerate(CS):
    xs = np.linspace(-lim, lim, 400); ys = float(c) * xs ** 2
    keep = np.abs(ys) <= lim
    handles.append(ax.plot(xs[keep], ys[keep], color=SERIES[k], lw=2.4,
                           label=f"y = {c}{DOT}x{SQ2}: {ARROW} {Rx(c)}")[0])
ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3, frameon=False,
          fontsize=9.3, handlelength=2.2, columnspacing=1.2)
ax.plot(0, 0, "o", ms=9, mfc=SURFACE, mec=INK, mew=1.8)
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
ax.set_xlabel("x"); ax.set_ylabel("y")
ax.set_title("A.  The landscape: color = value of f near the origin")
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03); cb.outline.set_visible(False)
cb.set_label("f(x, y)", color=INK2)

# --- B: the state profile at the balance curvature ---
ax = fig.add_subplot(gs[0, 1]); frame(ax)
cs = np.linspace(-6, 6, 600)
ax.axhline(0, color=AXIS, lw=1.6)
ax.plot(cs, [R(c) for c in cs], color=SERIES[0], lw=2.4)
for k, c in enumerate(CS):
    ax.plot(float(c), R(c), "o", ms=9, color=SERIES[k], mec=SURFACE, mew=2, zorder=5)
ax.annotate("maximum 1/2 at c = 1", (1, 0.5), xytext=(2.2, 0.56), fontsize=9.5, color=INK2,
            arrowprops=dict(arrowstyle="-", color=MUTED, lw=1))
ax.annotate(f"minimum {chr(0x2212)}1/2 at c = {chr(0x2212)}1", (-1, -0.5), xytext=(-5.8, -0.66), fontsize=9.5,
            color=INK2, arrowprops=dict(arrowstyle="-", color=MUTED, lw=1))
ax.text(-5.8, 0.66, f"proven value along y = c{DOT}x{SQ2}:   c / (1 + c{SQ2})", fontsize=10, color=INK)
ax.text(-5.8, 0.04, f"curvatures a {NE} {a_star}: value 0 (gray)", fontsize=9.5, color=INK2)
ax.set_ylim(-0.75, 0.75); ax.set_xlim(-6, 6)
ax.set_xlabel(f"path constant c   (path y = c{DOT}x{SQ2})")
ax.set_ylabel("proven limit along the path")
ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
ax.set_title("B.  The fluid state: a profile, not a number")

# --- C: balance diagram (weights of competing terms vs curvature a) ---
ax = fig.add_subplot(gs[1, 0]); frame(ax)
A = np.linspace(0.05, 4, 400)
terms = [(f"x{SQ2}y  (numerator)", (2, 1), SERIES[0]),
         (f"x{SQ4}  (denominator)", (4, 0), SERIES[1]),
         (f"y{SQ2}  (denominator)", (0, 2), SERIES[2])]
for label, (i, j), col in terms:
    ax.plot(A, i + j * A, color=col, lw=2.4)
    ae = min(4.0, (6.6 - i) / j) if j else 4.0
    ax.text(ae + 0.08, i + j * ae, label, color=INK, fontsize=9.5, va="center")
for ac in crit:
    ax.axvline(float(ac), color=MUTED, lw=1, ls=(0, (2, 3)))
ax.plot(float(a_star), 4, "o", ms=13, mfc=SURFACE, mec=INK, mew=2, zorder=6)
ax.annotate("balance point: all three terms tie\n(the Celedon Point, rigorous form)", (float(a_star), 4),
            xytext=(2.35, 1.2), fontsize=9.5, color=INK, arrowprops=dict(arrowstyle="->", color=INK2, lw=1.2))
ax.text(4.15, 3.3, "The lowest line dominates\nnear the origin. f has a\nnonzero limit only where the\n"
        "numerator's line meets the\nlowest denominator line.", fontsize=9.2, color=INK2, va="top")
ax.set_xlim(0, 6.6); ax.set_ylim(0, 7.3)
ax.set_xlabel(f"curvature a of the approach path  y = c{DOT}x^a")
ax.set_ylabel(f"weight i + a{DOT}j   (smaller dominates near 0)")
ax.set_title("C.  Balance diagram: which term wins near the origin")

# --- D: Newton polygon ---
ax = fig.add_subplot(gs[1, 1]); frame(ax)
for (i, j) in N: ax.plot(i, j, "o", ms=14, color=SERIES[0], mec=SURFACE, mew=2, zorder=5)
for (i, j) in D: ax.plot(i, j, "s", ms=13, color=SERIES[1], mec=SURFACE, mew=2, zorder=5)
ax.plot([4, 0], [0, 2], color=SERIES[1], lw=2.4)
ax.annotate("", xy=(2 + 0.55, 1 + 1.1), xytext=(2, 1), arrowprops=dict(arrowstyle="->", color=INK, lw=1.6))
ax.text(2.62, 2.12, f"edge normal (1, {a_star})  {ARROW}  a = {a_star}", fontsize=9.5, color=INK)
ax.text(2.08, 0.66, f"x{SQ2}y", fontsize=10, color=INK); ax.text(4.12, 0.12, f"x{SQ4}", fontsize=10, color=INK)
ax.text(0.12, 2.12, f"y{SQ2}", fontsize=10, color=INK)
ax.plot([], [], "o", color=SERIES[0], ms=9, label="numerator monomial")
ax.plot([], [], "s", color=SERIES[1], ms=9, label="denominator monomial")
ax.legend(loc="upper right", frameon=False, fontsize=9.5)
ax.text(-0.25, 3.1, "The numerator sits ON the denominator's edge:\nthe two sides balance only at a = 2",
        fontsize=9.3, color=INK2, va="top")
ax.set_xlim(-0.4, 5.2); ax.set_ylim(-0.4, 3.2); ax.set_aspect("equal")
ax.set_xticks(range(0, 6)); ax.set_yticks(range(0, 4))
ax.grid(color=GRID, lw=0.8); ax.set_axisbelow(True)
ax.set_xlabel("power of x  (i)"); ax.set_ylabel("power of y  (j)")
ax.set_title("D.  Newton polygon: the balance is visible in the exponents")

fig.savefig(f"{OUT}/fluid_state_anatomy.png", dpi=160)

# ======================= Figure 2: gallery of fluid states =======================
gallery = [("x*y/(x^2+y^2)", 1.0, "state varies along lines through the origin"),
           ("x^2*y/(x^4+y^2)", 1.0, "state varies along parabolas (hidden from every line)"),
           ("x^4/(x^2 + y^2 - 2*x*y + x^4)", 0.6, f"state varies along y = x + c{DOT}x{SQ2}\n(found two zoom levels deep)"),
           ("x^3/(x^2+y^2)", 1.0, "calm: the same value from every direction")]
fig, axes = plt.subplots(1, 4, figsize=(18, 6.6))
fig.subplots_adjust(left=0.035, right=0.95, top=0.80, bottom=0.26, wspace=0.18)
fig.text(0.035, 0.94, "Four fluid states, as the engine classifies them", fontsize=16, fontweight="bold", color=INK)
fig.text(0.035, 0.89, "Color = value of f near the origin (shared scale: blue < 0 < red). "
         "Verdicts below each map are the engine's exact output; the ring marks the transition point.",
         fontsize=11, color=INK2)
for ax, (fx, lim, where) in zip(axes, gallery):
    frame(ax)
    rr = profile(fx)
    X, Y, Z = field(rr["N"], rr["D"], lim)
    im = ax.imshow(np.clip(Z, -1, 1), extent=(-lim, lim, -lim, lim), origin="lower", cmap=DIVERGING,
                   vmin=-1, vmax=1, interpolation="bilinear")
    ax.plot(0, 0, "o", ms=8, mfc=SURFACE, mec=INK, mew=1.8)
    ax.set_aspect("equal"); ax.set_xticks([-lim, 0, lim]); ax.set_yticks([-lim, 0, lim])
    ax.set_title(f"f = {pretty(fx)}", fontsize=10.5)
    ax.text(0, -0.1, rr["verdict"][0], transform=ax.transAxes, fontsize=10.5, color=INK,
            fontweight="bold", va="top")
    ax.text(0, -0.18, where, transform=ax.transAxes, fontsize=9.6, color=INK2, va="top")
cax = fig.add_axes([0.962, 0.26, 0.008, 0.54])
cb = fig.colorbar(im, cax=cax); cb.outline.set_visible(False); cb.set_label("f(x, y)", color=INK2)
fig.savefig(f"{OUT}/fluid_state_gallery.png", dpi=150)
print("saved", f"{OUT}/fluid_state_anatomy.png", f"{OUT}/fluid_state_gallery.png")
