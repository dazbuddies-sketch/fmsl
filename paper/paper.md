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
approach. For a function $f$ of two or three variables, FMSL returns the limit with a proof when
it exists, and otherwise the complete set of values $f$ approaches, which we call its *fluid
state* (the cluster set of $f$ at the point). Viewing $f$ as a map into the circle
$\mathbb{R}\cup\{\infty\}$, this set is a point, a closed arc or the whole circle; FMSL returns
it exactly, with rational or real algebraic endpoints. Besides rational functions, $f$ may
contain distances such as $\sqrt{x^2+2y^2}$ and the functions $\sin$, $\cos$, $\exp$, $\log$ and
their relatives.

The same machinery gives, for a family $f(x,y,t)$, the exact times where the fluid state changes
character (*Celedon points*); the rate $\theta$ in $|f-L|\le C r^{\theta}$; whether a
polynomial's global minimum is attained; and the non-properness set of a polynomial map of the
plane [@jelonek1993]. A one-variable engine generates code with proven error bounds.

# Statement of need

Multivariate limits are a standard topic in calculus and a recurring source of errors in
computer algebra. SymPy [@sympy] handles one-variable limits [@gruntz1996] but has no
multivariable limit routine. @zengxiao2020 decide whether limits of bivariate rational functions
exist; @dgnp2022 compute the range of limit values when the denominator has an isolated zero,
report a case where Maple 2021 returned a wrong range, and leave other denominators to future
work; @strzebonski2021 reduces limits of analytic quotients to rational ones when the
denominator's zero is isolated. Cylindrical algebraic decomposition [@collins1975] decides such
questions in principle but does not return the structure of the answer.

FMSL is free, open source and dependency-free. It reproduces all 21 examples of @dgnp2022,
including the one Maple got wrong, computes the case they left open, and extends the computation
to three variables, distances, analytic functions and parametric families. Intended users are
educators, researchers, numerical-software authors, and AI assistants, which can call FMSL
through a Model Context Protocol server. A SymPy bridge returns SymPy
objects.

# Software design

The engines are pure Python with exact arithmetic; real algebraic numbers are isolated and
compared by Sturm sequences. The two-variable engine follows the Newton polygon with
Newton–Puiseux zooming; simple roots are settled by the implicit function and intermediate value
theorems, and completeness rests on the curve selection lemma [@milnor1968]. The three-variable
engine works on the arrangement of tie lines in the weight plane, zooms at rational and
algebraic points and along curves, and computes exact ranges of two-variable forms from critical
values, values at infinity and one exact test value per gap. The rate and the non-properness set
reuse the zooming, the latter along arcs to infinity $x=\pm s^{-1}$, $y=\sum c_i s^{e_i}$. A
square root of a positive-definite quadratic form is removed by a squaring map,
$X_1=u^2-D_2v^2-D_3w^2$, $X_2=2uv$, $X_3=2uw$, under which it becomes $s\,(u^2+D_2v^2+D_3w^2)$;
the map is onto and preserves approach to the point, so limit sets agree. For analytic $f=A/B$,
exact Taylor polynomials replace $A$ and $B$ once the engine proves $|B|\ge c\,r^L$, which makes
the remainders irrelevant; common factors along rational lines are first cancelled, with their
orders bounded by exact restriction to each line. Anything that cannot be decided is reported as
UNPROVEN, never guessed.

# Verification

Independent algorithms check each result type:

- **Fluid states.** Exact images on small squares agree on 113 of 113 two-variable functions,
  and on small cube surfaces on 102 of 105 three-variable functions (the other three are
  confirmed on a plane or by families of paths). Deliberately damaged answers were caught 53 of
  53 times.
- **Other results.** Moving targets (111 of 111 states), rates (139 functions), global extrema
  (67 polynomials) and non-properness sets (153 maps, 304 candidates) agree with exact
  independent methods.
- **Distances.** In two variables the squaring map and an unrelated polar blow-up agree on 74 of
  74 functions; in three, direction-only functions match numerical sphere ranges (36 of 36).
- **Analytic functions and weighted distances.** Direct 320-digit evaluation of the original
  formula, near the point, never leaves the proven set (126 and 38 decided functions); 32 of 32
  damaged answers were rejected.
- **Experiment A.** On 1,120 expressions with ground truth fixed by construction, FMSL is
  correct on all; SymPy 1.13.3 on 963, and wrong on 147.

# Limitations

In three variables FMSL reports UNPROVEN when it would have to zoom along a curve $c^pd^q=u_0$
with $|p|,|q|>1$, use a number field of degree above 16, or zoom deeper than four levels.
Celedon points that arise only deep inside the zoom are found by sampling and reported as
approximate. Analytic functions are proven when the denominator's zero becomes isolated after
cancelling common factors along rational lines, or the denominator is a polynomial whose zero
curves the numerator avoids; square roots only of quadratic forms that
take a rational square value. Zero-testing is undecidable in general [@richardson1968].

# Research impact statement

LimitBench, 716 multivariable limit problems with proven answers, tests AI models. Claude Sonnet
5.5 answered 707 exactly; checking with FMSL, it solved all 225 problems built around gaps. The
non-properness set is the central object of the real Jacobian conjecture, whose counterexample
by @pinchuk1994 has a nonempty one. The software is archived on Zenodo
(doi:10.5281/zenodo.23682015).

# AI usage disclosure

Generative AI was used extensively. Anthropic's Claude (model Claude Opus 5.5, via Claude Code)
wrote most of the source code, tests, benchmarks, documentation and the first draft of this
paper, working under the author's direction. The author set the goals and the research
direction, named the concepts, reviewed the results, and decided what to keep. All results are
validated by the automated tests and by the independent cross-check algorithms in the
repository.

# Acknowledgements

FMSL is the rigorous continuation of the author's Fluid Mathematics System; the Celedon points
and the scaling map are rigorous forms of its Celedon Point and Dimensional Unit Scaling.

# References
