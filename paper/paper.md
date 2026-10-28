---
title: 'FMSL: exact fluid states of 0/0 points, in several variables and for moving targets'
tags:
  - Python
  - multivariate limits
  - cluster sets
  - Newton polygon
  - Puiseux series
  - real algebraic numbers
  - parametric families
authors:
  - name: MJ Celedon
    affiliation: 1
affiliations:
  - name: Independent researcher
    index: 1
date: 28 October 2026
bibliography: paper.bib
---

# Summary

At a point where a function takes the form 0/0, a single number is often the wrong answer: in
several variables the value can depend on the direction, and even the curvature, of the path of
approach. For a rational function $f = N/D$ of two or three variables, FMSL returns the limit
with a proof when it exists, and otherwise the complete set of values $f$ approaches, which we
call its *fluid state* (the cluster set of $f$ at the point). Viewing $f$ as a map into the
circle $\mathbb{R}\cup\{\infty\}$, this set is a point, a closed arc or the whole circle; FMSL
returns it exactly, with rational or real algebraic endpoints.

The same exact machinery answers related questions. For a family $f(x,y,t)$ FMSL gives the
fluid state at each time and the exact times where it changes character, which we call
*Celedon points*. It computes how fast a limit is approached (the best $\theta$ with
$|f-L|\le C r^{\theta}$, a Łojasiewicz-type exponent, with a curve attaining it); the global
minimum of a polynomial in two variables and whether it is attained or only approached at
infinity; and the non-properness set of a polynomial map of the plane [@jelonek1993]. A
one-variable engine handles formulas such as $\sin t/t$ and generates code with proven error
bounds.

# Statement of need

Multivariate limits are a standard topic in calculus and a recurring source of errors in
computer algebra. SymPy [@sympy] handles one-variable limits [@gruntz1996] but has no
multivariable limit routine. @zengxiao2020 decide whether limits of bivariate rational functions
exist; @dgnp2022 compute the range of limit values when the denominator has an isolated zero,
report a case where Maple 2021 returned a wrong range, and leave other denominators to future
work. Cylindrical algebraic decomposition [@collins1975] decides such questions in principle
but does not return the structure of the answer.

FMSL is free, open source and dependency-free. It reproduces all 21 examples of @dgnp2022,
including the one Maple got wrong, computes the case they left open, and extends the
computation to three variables and to parametric families. Intended users are educators,
researchers, authors of numerical software, and AI assistants, which can call FMSL through a
Model Context Protocol server to check their answers. A SymPy bridge returns SymPy objects, and inclusion in SymPy has been proposed.

# Software design

The engines are pure Python with exact arithmetic: real algebraic numbers are square-free
polynomials with isolating intervals, compared by Sturm sequences, also over number fields. The two-variable engine follows the Newton polygon with Newton–Puiseux zooming; simple
roots are settled by the implicit function and intermediate value theorems, and completeness
rests on the curve selection lemma [@milnor1968]. The three-variable engine works on the
arrangement of tie lines in the weight plane, zooms at rational and algebraic points and along
curves, and computes exact ranges of two-variable forms from critical values, values at infinity
and one exact test value per gap. The rate and the non-properness set reuse the zooming, the latter along arcs to infinity
$x=\pm s^{-1}$, $y=\sum c_i s^{e_i}$. Anything that cannot be decided
is reported as UNPROVEN, never guessed.

# Verification

Every result type is checked by an algorithm that shares no logic with the engine:

- **Fluid states.** Exact images on small squares agree on 113 of 113 two-variable functions,
  and on small cube surfaces on 102 of 105 three-variable functions (the other three are
  confirmed on a plane or by families of paths). Deliberately damaged answers were caught 53 of
  53 times.
- **Moving targets.** 75 of 75 snapshots and 36 of 36 moving states agree with the exact images.
- **Rate.** On 139 functions, exact rational points near the origin never contradict the bound,
  and on every reported curve the measured order matches $\theta$.
- **Global extrema.** On 67 polynomials, an exact real-points test confirms that no value lies
  below the infimum, that it is approached, and whether it is attained.
- **Non-properness set.** On 153 maps and 612 generic lines, candidates found by elimination and
  decided by exact solving agree with FMSL in all 304 cases.
- **Experiment A.** 1,120 expressions with ground truth fixed by construction, following the
  author's earlier Fluid Mathematics System research brief, run against SymPy 1.13.3: FMSL 1,120
  correct; SymPy 963 correct, 147 wrong, 10 unresolved. In several variables, limits along
  standard paths found the exact set of values in 67 of 224 cases.

# Limitations

In three variables FMSL reports UNPROVEN when it would have to zoom along a curve $c^pd^q=u_0$
with $|p|,|q|>1$, use a number field of degree above 16, or zoom deeper than four levels.
Celedon points that arise only deep inside the zoom are found by sampling and reported as
approximate. The multivariable engines handle rational functions; with exponentials and
absolute values, equality can be undecidable [@richardson1968].

# Research impact statement

LimitBench, 716 multivariable limit problems with proven answers, tests AI models. In October
2026 Claude Sonnet 5.5 answered 707 exactly, and all 7 of its errors on the subset built around
gaps in the set of limit values claimed that every value is approached. Told to check its
answers with FMSL, it answered all 225 problems of that subset exactly. The non-properness set
is the central object of the real Jacobian conjecture, whose counterexample by @pinchuk1994 has
a nonempty one. The software is archived on Zenodo (doi:10.5281/zenodo.23682015).

# AI usage disclosure

Generative AI was used extensively. Anthropic's Claude (model Claude Opus 5.5, via Claude Code)
wrote most of the source code, tests, benchmarks, documentation and the first draft of this
paper, working under the author's direction. The author set the goals and the research
direction, named the concepts, reviewed the results, and decided what to keep. All results are
validated by the automated tests and by the independent cross-check algorithms in the repository.

# Acknowledgements

FMSL is the rigorous continuation of the author's Fluid Mathematics System; the Celedon points
and the scaling map are rigorous forms of its Celedon Point and Dimensional Unit Scaling.

# References
