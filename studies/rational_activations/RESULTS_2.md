# Results 2: unsafe Padé activation units

Pre-registration: `PREREGISTRATION_2.md` (committed before any run). Code: `train2.py`
(gradients checked, worst relative difference 1e-6), `analyze2.py`. Data: `results2.json`.

## Pre-registered result: NOT TESTABLE

40 runs (learning rates 0.002, 0.01, 0.03, 0.1, 0.3 x 8 seeds) of the same network as
experiment 1 with unsafe Padé activation units (degrees 5/4, denominator 1 + b1 x + ... + b4 x^4),
initialised at the paper's LeakyReLU(0.2) coefficients.

| | runs |
|---|---|
| unstable (loss not finite, or above its untrained value, at step 100 or later) | 3 of 40 |
| a genuine pole inside the activation's input range (exact test) | 1 of 40 |
| a genuine real pole anywhere (exact test) | 21 of 40 |

The pre-registration requires at least 5 unstable runs and 5 runs with a pole event, so the
verdict is NOT TESTABLE.

## Exploratory observations (not pre-registered; weaker evidence)

1. **Unsafe Padé units do form real poles**: in 21 of 40 runs, against 5 of 40 for the unit of
   experiment 1. In all but one run they stayed outside the data.
2. **The one pole that reached the data coincided with a real blow-up.** Learning rate 0.3,
   seed 7: at step 2 the second activation had genuine poles at -3.82 and 5.42 inside its input
   range [-37.5, 24.3], and at that step the loss rose from 0.63 to 9.3e8. The pre-registered
   analysis could not count this, because it only looked for instability from step 100 on.
3. **The exact test caught a pole the pre-registered float check missed.** Near that pole the
   denominator is steep (|Q'| of 27 to 117), so no point of a 1,000-point grid came within 1e-3
   of zero. A float check that looks for sign changes would have caught it, so the advantage is
   over the naive check, not over every float check.
4. **Not every blow-up is a pole.** Learning rate 0.3, seed 5: the loss reached 11,777 at step 1
   while the poles were at -20.3 and -7.7, outside the inputs [-4.4, 4.5]. Very large learning
   rates break training on their own.

## Conclusion

Across both experiments, on this task, poles of trainable rational activations rarely reach the
data, even when the safeguard is removed. One case shows that when a pole does reach the data,
the loss explodes at that step, and exact detection catches it where a naive grid check does not.
One case is an anecdote, not a result. Establishing a link would need a pre-registered study
designed to produce many pole events (larger learning rates, harder tasks, counting from step 1),
and is left there.
