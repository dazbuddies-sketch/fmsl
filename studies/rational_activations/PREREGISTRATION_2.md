# Pre-registration 2: unsafe Padé activation units

Written and committed BEFORE any run of this second experiment. Results go in RESULTS_2.md.

## Why
Molina, Schramowski & Kersting, "Padé Activation Units" (ICLR 2020, arXiv:1907.06732) use a
"safe" denominator 1 + |b1 x + ... + b4 x^4| and report no experiments on the unsafe version
(denominator 1 + b1 x + ... + b4 x^4), leaving it to future work. Experiment 1 (Boullé et al.'s
unit) found no pole ever reaching the data; this tests the case its authors flagged as risky.

## Unit
Unsafe PAU, degrees (5, 4): R(x) = (a0 + ... + a5 x^5) / (1 + b1 x + ... + b4 x^4), one trainable
unit per hidden layer, initialised at the paper's Table 6 LeakyReLU (slope 0.2) coefficients
a = (0.02557776, 0.66182815, 1.58182975, 2.94478759, 0.95287794, 0.23319681),
b = (0.50962605, 4.18376890, 0.37832090, 0.32407314).
Checked before any run: the safe form approximates LeakyReLU(0.2) on [-3, 3] with max error
0.28; the unsafe denominator has no real root at initialisation (min 0.984 on [-3, 3]).

## Setup (as in experiment 1, so the two are comparable)
MLP 1 -> 32 -> 32 -> 1, same task (y = sin(3x) + 0.3 x^2, x in [-2, 2], 2000 points), Glorot
weights, Adam, batch 128, 3000 steps. Learning rates {0.002 (the paper's Adam setting), 0.01,
0.03, 0.1, 0.3} x 8 seeds = 40 runs. Exact gradients, checked by finite differences first.

## Events
POLE-IN-RANGE (exact, FMSL): as in experiment 1, with Q of degree 4: a root of
Q / gcd(P, Q) in the range of the activation's inputs on the full training set, decided with
Sturm sequences on the exact rational values of the float coefficients.
INSTABILITY (absolute, stricter than experiment 1): first step t >= 100 where the full-set loss
is not finite or exceeds its value at step 0 (the untrained network).

## Hypothesis and criteria
H1: every unstable run has a POLE-IN-RANGE event at or before its instability step, at most 25
    steps earlier.
H2: every run with a POLE-IN-RANGE event becomes unstable within 50 steps after it.
SUPPORTED if both hold in at least 90% of the relevant runs, with at least 5 unstable runs and
at least 5 runs with a pole event. REFUTED if either holds in fewer than 60% of at least 5
relevant runs. NOT TESTABLE if fewer than 5 runs are unstable or fewer than 5 have pole events.
Otherwise PARTIAL.
Baselines as before: B1 float grid check (min |Q| < 1e-3), B2 gradient-norm spike.
