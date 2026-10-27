# Results: exact pole events and training instability in rational activations

Pre-registration: `PREREGISTRATION.md` (committed before any run). Code: `train.py` (gradients
checked against finite differences, worst relative difference 4e-8), `analyze.py` (the
pre-registered analysis; exact pole tests with `fmsl.algebraic`). Data: `results.json`.

## Pre-registered result: H is REFUTED

40 runs (5 learning rates x 8 seeds, 3,000 steps each) of a 1-32-32-1 network with the
type-(3,2) rational activation of Boullé, Nakatsukasa & Townsend (2020), initialised as in
their Table 1.

| | runs |
|---|---|
| unstable (pre-registered definition) | 30 of 40 |
| a genuine pole ever inside the activation's input range (exact test) | **0 of 40** |
| the same, float check (min \|Q\| < 1e-3 on a 1,000-point grid) | 0 of 40 |

- **H1** (every unstable run has a pole in range shortly before): **0 of 30** -> below the
  pre-registered 60% threshold, so **H is refuted**.
- **H2** (pole events are followed by instability): not testable; no pole event occurred.
- The exact test and the float check agreed at every one of the 240,000 checked steps, so in
  this setting exactness added certainty but no predictive power.
- Gradient-norm spikes (baseline B2) matched the instability step in 1 of 30 runs.

Note: `analyze.py` prints the verdict "NOT TESTABLE" because no run was flagged; the
pre-registration's rule ("REFUTED if either holds in fewer than 60%") applies to H1 at 0%, so
the verdict is REFUTED. The pre-registration takes precedence.

## Exploratory findings (not pre-registered)

1. **The pre-registered instability definition was too sensitive.** The typical "instability" is
   a spike from a loss of about 2e-5 to 3e-4 late in training: ordinary optimiser noise, not a
   breakdown. No run produced a non-finite loss.

   | learning rate | runs flagged | median spike loss (from) | smallest min\|Q\| on the inputs during training |
   |---|---|---|---|
   | 0.001 | 0/8 | - | 0.992 |
   | 0.003 | 6/8 | 3.0e-4 (2.4e-5) | 0.679 |
   | 0.01 | 8/8 | 5.0e-4 (3.7e-5) | 0.405 |
   | 0.03 | 8/8 | 3.7e-3 (2.3e-4) | 2.159 |
   | 0.1 | 8/8 | 1.6e-2 (1.3e-3) | 2.281 |

2. **Real poles appeared during training, but stayed away from the data.** In 5 runs (learning
   rates 0.003 and 0.01) the first layer's denominator acquired real roots, i.e. the activation
   developed genuine poles on the real line. In every case they sat far outside the inputs the
   activation ever received:

   | run | real poles of layer 1 | appeared at step | input range at that step |
   |---|---|---|---|
   | lr 0.003, seed 2 | 6.73, 8.52 | 363 | [-1.83, 1.84] |
   | lr 0.01, seed 1 | 5.86, 10.05 | 83 | [-2.04, 2.60] |
   | lr 0.01, seed 2 | -79.77, 414.85 | 110 | [-2.17, 2.32] |
   | lr 0.01, seed 5 | 5.15, 8.19 | 126 | [-2.43, 2.94] |
   | lr 0.01, seed 7 | 4.91, 5.74 | 1159 | [-2.01, 1.75] |

   This is consistent with the authors' statement that a good initialisation keeps the network
   away from poles, and FMSL makes it exact: the pole never entered the data range at any step.

## What this means

In this setting, with the published initialisation, training instability of rational
activations is **not** caused by poles reaching the data. A write-up claiming that exact pole
events predict instability would be false, so none is made. The negative result itself, and the
observation that poles do form but stay outside the data, are reported here as they are.

## The open question that remains

Molina, Schramowski & Kersting (ICLR 2020) avoid poles with a "safe" denominator and state that
they ran no experiments on the unsafe version, leaving it to future work. The experiment that
could still show a pole-instability link is: unsafe Padé activation units (degrees 5/4,
denominator 1 + b1 x + ... + b4 x^4), initialised as in that paper, with instability defined in
absolute terms (a non-finite loss, or a loss above its initial value). That needs its own
pre-registration before it is run.
