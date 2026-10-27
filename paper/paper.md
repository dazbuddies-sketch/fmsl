---
title: 'FMSL: exact fluid states of 0/0 points in two and three variables'
tags:
  - Python
  - multivariate limits
  - cluster sets
  - Newton polygon
  - Puiseux series
  - real algebraic numbers
authors:
  - name: MJ Celedon
    affiliation: 1
affiliations:
  - name: Independent researcher
    index: 1
date: 27 October 2026
bibliography: paper.bib
---

# Summary

At a point where a function takes the form 0/0, a single number is often the wrong answer.
For a function of several variables the value can depend on the direction, and even on the
curvature, of the path of approach. FMSL computes exactly what happens there. For a rational
function $f = N/D$ of two or three variables it returns the limit with a proof when one exists.
Otherwise it returns the complete set of values $f$ approaches, which we call its *fluid
state* (the cluster set of $f$ at the point).

Treating $f$ as a continuous map into the circle $\mathbb{R}\cup\{\infty\}$ shows that this set
is always closed and connected: a single point, a closed arc, or the whole circle. FMSL returns
it exactly. Endpoints are rationals or real algebraic numbers, decided with Sturm sequences, and
a signed variant keeps $+\infty$ and $-\infty$ apart.

Hidden paths are found from the Newton polygon (in three variables, the arrangement of weight
lines in the $(a,b)$ plane), with Newton–Puiseux zooming where initial forms have repeated or
shared roots. Simple roots are settled by the implicit function and intermediate value
theorems, which also covers pole curves whose Puiseux series never terminates. Completeness rests
on the curve selection lemma [@milnor1968]. Every reported path value is re-checked by direct
exact evaluation, and every fluid state by an independent algorithm: the exact image of $f$ on
small squares (two variables) or cube surfaces (three variables), which uses no Newton polygons.

FMSL also contains a one-variable engine for 0/0 traps in formulas such as $\sin t/t$, with
generated code that carries a proven error bound, and an exact analysis of bifurcations of
iterated maps as 0/0 balance points.

# Statement of need

Multivariate limits are a standard topic in calculus and a recurring source of errors in
computer algebra. SymPy [@sympy] handles one-variable limits [@gruntz1996]. Đinh, Guo, Nguyễn
and Phạm [@dgnp2022] gave an algorithm, implemented in Maple, that computes the range of limit
values of a real bivariate rational function when the denominator has an isolated zero. They
also report a case where Maple 2021's `limit` returned a wrong range, and they leave the case of
a denominator vanishing along curves to future work. Cylindrical algebraic decomposition
[@collins1975] decides such questions in principle, but at a cost that limits practical use.

FMSL is free, open source, and dependency-free for its exact core. It reproduces all 21
examples of @dgnp2022 (including the one Maple got wrong), and it computes ranges when the
denominator vanishes along curves (for example, $[-\infty,-3]\cup[3,+\infty]$ for their
Example 5). It extends the computation to three variables, where the vertex forms are genuine
two-variable polynomials. Their exact range over the plane is found from critical values,
values at infinity, and an exact cylindrical test of one value per gap. Intended users are
educators showing why a limit fails, researchers in real algebraic geometry who need exact
experiments, and authors of numerical software who must handle removable singularities safely.

# Verification

The repository's tests and benchmarks reproduce the following:

- **Path values.** Thousands of reported path values agree with direct exact evaluation.
- **Fluid states.** Two variables: 113 of 113 agree with the independent loop computation.
  Three variables: 51 of 51 agree with the cube computation.
- **Deliberately broken answers.** The independent checks caught 42 of 42 in two variables and
  11 of 11 in three.
- **The two-variable range algorithm** agrees with dense sampling on 120 of 120 random cases.

# Limitations

In three variables, version 1 reports UNPROVEN, and never guesses, when the analysis would need
to zoom along a whole curve of special points or at irrational points (about 1 in 200 random
functions). The exact multivariate engines handle rational functions; for
expressions mixing exponentials and absolute values, equality can be undecidable
[@richardson1968].

# Acknowledgements

FMSL is the rigorous continuation of the author's Fluid Mathematics System. Development was
assisted by an AI system (Anthropic's Claude), which wrote much of the code and tests under the
author's direction; all results are verified by the tests and independent cross-checks in the
repository.

# References
