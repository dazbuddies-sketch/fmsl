# Experiment A: singularity classification, head to head

96 one-variable and 24 several-variable expressions (120 total); ground truth by construction; time limit 20 s each.
Unresolved = no answer, undecided, timed out, or not supported. Several variables: 'correct' means the exact set of limit values; 'verdict' means limit / no limit only.

| category | method | n | correct | wrong | unresolved | verdict right | time (s) |
|---|---|---|---|---|---|---|---|
| 1 removable | FMSL | 12 | 12 | 0 | 0 | - | 0.0 |
| 1 removable | float | 12 | 9 | 3 | 0 | - | 0.0 |
| 10 multivariable | FMSL | 12 | 12 | 0 | 0 | 12 | 0.2 |
| 2 pole (even order) | FMSL | 12 | 12 | 0 | 0 | - | 0.0 |
| 2 pole (even order) | float | 12 | 12 | 0 | 0 | - | 0.0 |
| 3 one-sided divergence | FMSL | 12 | 12 | 0 | 0 | - | 0.0 |
| 3 one-sided divergence | float | 12 | 12 | 0 | 0 | - | 0.0 |
| 4 oscillatory | FMSL | 12 | 12 | 0 | 0 | - | 0.0 |
| 4 oscillatory | float | 12 | 9 | 3 | 0 | - | 0.0 |
| 5 piecewise | FMSL | 12 | 12 | 0 | 0 | - | 0.0 |
| 5 piecewise | float | 12 | 12 | 0 | 0 | - | 0.0 |
| 6 catastrophic cancellation | FMSL | 12 | 12 | 0 | 0 | - | 0.0 |
| 6 catastrophic cancellation | float | 12 | 0 | 12 | 0 | - | 0.0 |
| 7 ill-conditioned | FMSL | 12 | 12 | 0 | 0 | - | 0.0 |
| 7 ill-conditioned | float | 12 | 0 | 12 | 0 | - | 0.0 |
| 8 parameter-dependent | FMSL | 12 | 12 | 0 | 0 | 12 | 0.1 |
| 9 essential | FMSL | 12 | 12 | 0 | 0 | - | 0.0 |
| 9 essential | float | 12 | 4 | 0 | 8 | - | 0.0 |

| method | n | correct | wrong | unresolved | time (s) |
|---|---|---|---|---|---|
| FMSL | 120 | 120 | 0 | 0 | 0.4 |
| float | 96 | 58 | 30 | 8 | 0.0 |
