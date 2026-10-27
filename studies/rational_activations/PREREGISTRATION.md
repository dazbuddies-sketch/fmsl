# Pre-registration: do exact pole events predict training instability in rational activations?

Written and committed BEFORE any training run. The git timestamp of this file is the evidence
that the hypothesis and criteria were fixed in advance. Results go in RESULTS.md, whatever they are.

## Unit studied
The rational activation of Boullé, Nakatsukasa & Townsend, "Rational neural networks"
(NeurIPS 2020, arXiv:2004.01902): R(x) = P(x)/Q(x), type (3,2),
P = a0 + a1 x + a2 x^2 + a3 x^3, Q = b0 + b1 x + b2 x^2, one trainable R per layer, initialised
at their Table 1 values a = (1.1915, 1.5957, 0.5000, 0.0218), b = (2.3830, 0.0000, 1.0000)
(best rational approximation of ReLU on [-1, 1]). No safeguard on Q (as in the paper).

## Setup
MLP 1 -> 32 -> 32 -> 1, rational activation after each hidden layer (2 activations, 7
coefficients each), Glorot-uniform weights, exact gradients (NumPy), Adam, mini-batch 128.
Task: regression of y = sin(3x) + 0.3 x^2 on x uniform in [-2, 2] (2000 points, fixed).
Grid: learning rates {1e-3, 3e-3, 1e-2, 3e-2, 1e-1} x 8 seeds = 40 runs, 3000 steps each.

## Exact event (FMSL), computed at every step for each activation
The float coefficients are converted exactly to rationals. Let [lo, hi] be the range of the
pre-activations entering that activation at that step (on the full training set).
POLE-IN-RANGE at step t: Q has a real root r in [lo, hi] that is a genuine pole, i.e. not
cancelled by P (r is a root of Q of higher multiplicity than of P). Decided exactly with Sturm
sequences and polynomial gcds (fmsl.algebraic); a root cancelled by P is a removable 0/0 point
and does NOT count. Also recorded: REAL-POLE (a genuine real pole anywhere on the real line;
the coefficient-space boundary is the discriminant b1^2 - 4 b0 b2 = 0, or b2 = 0).

## Instability event
First step t >= 100 where the training loss is not finite, or exceeds 10 x the median of the
previous 50 steps. A run with no such step is STABLE.

## Hypothesis H and criteria (fixed now)
H1 (necessity): every unstable run has a POLE-IN-RANGE event at or before its instability step,
    at most 25 steps earlier.
H2 (sufficiency): every run with a POLE-IN-RANGE event becomes unstable within 50 steps after it.
H is SUPPORTED if both hold in at least 90% of the relevant runs and at least 5 runs are
unstable. H is REFUTED if either holds in fewer than 60%. Otherwise: PARTIAL.

## Baselines (to judge what exactness adds)
B1: a float check, min |Q| over 1000 grid points of [lo, hi] below 1e-3.
B2: gradient-norm spike (norm > 10 x median of previous 50 steps).
For each, the same H1/H2 scores are reported. If B1 does as well as the exact test, that is
reported plainly: exactness then adds certainty, not predictive power, on this task.

## What would be worth writing up
H supported AND the exact test separates cases a float check gets wrong (narrow pole pairs,
or removable 0/0 points that are not poles). If H is refuted or partial, that is written up
honestly as a negative or mixed result.
