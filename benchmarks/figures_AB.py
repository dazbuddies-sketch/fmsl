"""Figures for Part A (ensemble vs exact) and Part B (bifurcations). Reads the JSON written by
ensemble_vs_exact.py and bifurcation_demo.py.  Run:  python benchmarks/figures_AB.py [dir]"""
import json, sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = sys.argv[1] if len(sys.argv) > 1 else "."
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#e34948"]
GOOD, BAD = "#0ca30c", "#d03b3b"
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

import os
# ================================ Part A ================================
def part_a():
  E = json.load(open(f"{D}/ensemble_results.json"))
  fig = plt.figure(figsize=(17, 6.6))
  gs = fig.add_gridspec(1, 3, wspace=0.32, left=0.09, right=0.98, top=0.78, bottom=0.14)
  rnd = E["random"]
  fig.text(0.05, 0.94, "Repeating the question: random sampling vs exact analysis", fontsize=16, fontweight="bold")
  fig.text(0.05, 0.895, f"Monte Carlo = 1,000 random approach paths per function, judged by the spread of results. "
           f"On {rnd['n']} random functions it matched the exact verdict {rnd['agree']} times "
           f"({100 * rnd['agree'] / rnd['n']:.0f}%). Where it fails is the interesting part.", fontsize=11, color=INK2)

  # A1: agreement by exact category
  ax = fig.add_subplot(gs[0, 0]); frame(ax)
  cats = ["LIMIT", "NO LIMIT", "DIVERGES"]
  conf = rnd["confusion"]
  tot = [sum(v for k, v in conf.items() if k.split("|")[0] == c) for c in cats]
  ok = [conf.get(f"{c}|{c}", 0) for c in cats]
  y = np.arange(len(cats))
  ax.barh(y, ok, color=GOOD, height=0.55, label="sampling agrees")
  ax.barh(y, [t - o for t, o in zip(tot, ok)], left=ok, color=BAD, height=0.55, label="sampling wrong")
  for i, (t, o) in enumerate(zip(tot, ok)):
      ax.text(t + max(tot) * 0.02, i, f"{o}/{t}", va="center", fontsize=10, color=INK)
  ax.set_yticks(y); ax.set_yticklabels([f"exact: {c}" for c in cats]); ax.invert_yaxis()
  ax.set_xlabel("number of random functions"); ax.set_xlim(0, max(tot) * 1.18)
  ax.legend(frameon=False, loc="lower right", fontsize=9.5)
  ax.set_title("A1.  Does sampling reach the exact verdict?")

  # A2: detection probability vs repetitions
  ax = fig.add_subplot(gs[0, 1]); frame(ax)
  Ks = E["Ks"]
  for k, (name, probs) in enumerate(E["detect"].items()):
      short = name.split("  ")[0]
      if probs[-1] == 0: short += "  (never found)"
      ax.plot(Ks, probs, marker="o", color=SERIES[k], lw=2.2, ms=6 + 3 * (k == 1), mec=SURFACE, mew=1.5, label=short)
  ax.legend(frameon=False, loc="center right", fontsize=9.5)
  ax.set_xscale("log"); ax.set_ylim(-0.04, 1.06); ax.set_xlim(0.8, Ks[-1] * 1.5)
  ax.set_xlabel("number of random paths K (repetitions)")
  ax.set_ylabel("chance the trap is detected")
  ax.grid(color=GRID, lw=0.8); ax.set_axisbelow(True)
  ax.set_title("A2.  More repetitions help only some traps")

  # A3: what the samples show for the parabola trap
  ax = fig.add_subplot(gs[0, 2]); frame(ax)
  v = np.array(E["parabola_samples"])
  ax.axvspan(-0.5, 0.5, color=SERIES[0], alpha=0.08)
  bins = np.linspace(-0.55, 0.55, 111)
  ax.hist(v, bins=bins, color=SERIES[0], alpha=0.95)
  ax.set_yscale("log")
  ax.axvline(-0.5, color=SERIES[0], lw=1.2, ls=(0, (3, 3))); ax.axvline(0.5, color=SERIES[0], lw=1.2, ls=(0, (3, 3)))
  ax.text(-0.47, ax.get_ylim()[1] * 0.35, "exact profile range\n[-1/2, 1/2]", ha="left", fontsize=9.5, color=INK2)
  ax.set_xlabel("value of f along a sampled path (x near 0)")
  ax.set_ylabel("number of paths (log scale)")
  ax.set_title("A3.  Parabola trap: samples sketch the profile")
  fig.savefig(f"{D}/ensemble_vs_exact.png", dpi=150)

if os.path.exists(f"{D}/ensemble_results.json"): part_a()

# ================================ Part B ================================
B = json.load(open(f"{D}/bifurcation_results.json"))
fig = plt.figure(figsize=(17, 6.8))
gs = fig.add_gridspec(1, 3, wspace=0.3, left=0.05, right=0.98, top=0.78, bottom=0.13, width_ratios=[1.35, 1, 0.85])
fig.text(0.05, 0.94, "Repeating the function: bifurcations are balance points", fontsize=16, fontweight="bold")
fig.text(0.05, 0.895, "Iterating x → f(x) millions of times (dots) versus the exact Newton-polygon predictions (lines). "
         "Each new orbit is born along a branch e = c·|u|^a found exactly at a 0/0 point.", fontsize=11, color=INK2)

# B1: logistic bifurcation diagram with exact branches
ax = fig.add_subplot(gs[0, 0]); frame(ax)
ax.scatter(B["diagram"]["r"], B["diagram"]["x"], s=0.25, color=INK2, alpha=0.45, linewidths=0, rasterized=True)
r = np.linspace(2.85, 3.0, 50); ax.plot(r, 1 - 1 / r, color=INK, lw=1.4)
cases = {c["name"]: c for c in B["cases"]}
c1 = cases["period doubling 1 (logistic, r0 = 3)"]["c"]
e = np.linspace(0, 0.3, 200)
ax.plot(3 + e, 2 / 3 + np.sqrt(e / c1), color=SERIES[0], lw=2.2)
ax.plot(3 + e, 2 / 3 - np.sqrt(e / c1), color=SERIES[0], lw=2.2)
c2 = cases["period doubling 2 (logistic, r0 = 1+sqrt6)"]["c"]
x4, r4 = B["x4"], B["r4"]
e = np.linspace(0, 0.09, 200)
ax.plot(r4 + e, x4 + np.sqrt(e / c2), color=SERIES[1], lw=2.2)
ax.plot(r4 + e, x4 - np.sqrt(e / c2), color=SERIES[1], lw=2.2)
for (rr, xx) in ((3, 2 / 3), (r4, x4)):
    ax.plot(rr, xx, "o", ms=9, mfc=SURFACE, mec=INK, mew=1.8, zorder=5)
ax.annotate("r = 3: e = 9u²\n(balance curvature 2)", (3, 2 / 3), xytext=(2.88, 0.80), fontsize=9.5,
            arrowprops=dict(arrowstyle="-", color=MUTED))
ax.annotate(f"r = 1+√6: e = {c2:.2f}·u²\n(c exact in ℚ(√2, √3))", (r4, x4),
            xytext=(3.235, 0.64), fontsize=9.5, arrowprops=dict(arrowstyle="-", color=MUTED))
ax.set_xlim(2.85, 3.57); ax.set_ylim(0.3, 0.95)
ax.set_xlabel("parameter r"); ax.set_ylabel("long-run values of x")
ax.set_title("B1.  Logistic map x → r·x(1−x): iteration vs exact branches")

# B2: measured orbit size vs distance from the bifurcation
ax = fig.add_subplot(gs[0, 1]); frame(ax)
for k, c in enumerate(B["cases"]):
    es = np.array([row["e"] for row in c["rows"]]); pred = np.array([row["pred"] for row in c["rows"]])
    meas = np.array([row["meas"] for row in c["rows"]])
    hollow = c["name"].startswith("pitchfork")             # identical to saddle-node (both sqrt(e))
    ax.plot(es, pred, color=SERIES[k], lw=1.8, ls=(0, (4, 2)) if hollow else "-")
    lab = c["name"].split(" (")[0] + f"  (measured slope {c['slopes'][-1]:.3f})"
    ax.plot(es, meas, "o", color=SERIES[k], ms=11 if hollow else 8, mfc=SURFACE if hollow else SERIES[k],
            mec=SERIES[k] if hollow else SURFACE, mew=2 if hollow else 1.5, label=lab)
ax.legend(loc="lower right", frameon=False, fontsize=9, handletextpad=0.4)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(5e-6, 3e-2); ax.set_ylim(3e-6, 0.3)
ax.set_xlabel("distance from the bifurcation  e = r − r₀")
ax.set_ylabel("size of the new orbit |u|")
ax.grid(color=GRID, lw=0.8, which="major"); ax.set_axisbelow(True)
ax.set_title("B2.  Dots: iteration.  Lines: exact prediction")

# B3: Newton polygon of Q2 at the first period doubling
ax = fig.add_subplot(gs[0, 2]); frame(ax)
pts = B["newton_Q2"]
for i, j, v in pts:
    ax.plot(i, j, "o", ms=12, color=SERIES[0], mec=SURFACE, mew=2, zorder=5)
ax.plot([2, 0], [0, 1], color=SERIES[1], lw=2.4)
ax.text(2.25, 0.3, "balance edge: u² against e\n⇒ e ∝ u² (curvature 2)", fontsize=9.5, color=INK)
ax.set_xlim(-0.4, 4.4); ax.set_ylim(-0.4, 2.6); ax.set_aspect("equal")
ax.set_xticks(range(0, 5)); ax.set_yticks(range(0, 3))
ax.grid(color=GRID, lw=0.8); ax.set_axisbelow(True)
ax.set_xlabel("power of u = x − 2/3"); ax.set_ylabel("power of e = r − 3")
ax.set_title("B3.  Newton polygon of Q₂ at r = 3")
fig.savefig(f"{D}/bifurcations.png", dpi=150)
print("saved", f"{D}/ensemble_vs_exact.png", f"{D}/bifurcations.png")
