# Experiment A: singularity classification, head to head

896 one-variable and 224 several-variable expressions (1120 total); ground truth by construction; time limit 20 s each.
Unresolved = no answer, undecided, timed out, or not supported. Several variables: 'correct' means the exact set of limit values; 'verdict' means limit / no limit only.

| category | method | n | correct | wrong | unresolved | verdict right | time (s) |
|---|---|---|---|---|---|---|---|
| 1 removable | FMSL | 112 | 112 | 0 | 0 | - | 0.0 |
| 1 removable | SymPy | 112 | 112 | 0 | 0 | - | 0.4 |
| 1 removable | float | 112 | 92 | 20 | 0 | - | 0.0 |
| 10 multivariable | FMSL | 112 | 112 | 0 | 0 | 112 | 0.9 |
| 10 multivariable | SymPy | 112 | 63 | 49 | 0 | 103 | 0.4 |
| 2 pole (even order) | FMSL | 112 | 112 | 0 | 0 | - | 0.0 |
| 2 pole (even order) | SymPy | 112 | 112 | 0 | 0 | - | 0.3 |
| 2 pole (even order) | float | 112 | 112 | 0 | 0 | - | 0.0 |
| 3 one-sided divergence | FMSL | 112 | 112 | 0 | 0 | - | 0.0 |
| 3 one-sided divergence | SymPy | 112 | 112 | 0 | 0 | - | 0.2 |
| 3 one-sided divergence | float | 112 | 112 | 0 | 0 | - | 0.0 |
| 4 oscillatory | FMSL | 112 | 112 | 0 | 0 | - | 0.0 |
| 4 oscillatory | SymPy | 112 | 112 | 0 | 0 | - | 0.1 |
| 4 oscillatory | float | 112 | 77 | 35 | 0 | - | 0.0 |
| 5 piecewise | FMSL | 112 | 112 | 0 | 0 | - | 0.0 |
| 5 piecewise | SymPy | 112 | 112 | 0 | 0 | - | 0.1 |
| 5 piecewise | float | 112 | 112 | 0 | 0 | - | 0.0 |
| 6 catastrophic cancellation | FMSL | 112 | 112 | 0 | 0 | - | 0.1 |
| 6 catastrophic cancellation | SymPy | 112 | 112 | 0 | 0 | - | 0.4 |
| 6 catastrophic cancellation | float | 112 | 0 | 112 | 0 | - | 0.0 |
| 7 ill-conditioned | FMSL | 112 | 112 | 0 | 0 | - | 0.2 |
| 7 ill-conditioned | SymPy | 112 | 112 | 0 | 0 | - | 0.4 |
| 7 ill-conditioned | float | 112 | 0 | 112 | 0 | - | 0.0 |
| 8 parameter-dependent | FMSL | 112 | 112 | 0 | 0 | 112 | 0.5 |
| 8 parameter-dependent | SymPy | 112 | 4 | 98 | 10 | 102 | 0.3 |
| 9 essential | FMSL | 112 | 112 | 0 | 0 | - | 0.0 |
| 9 essential | SymPy | 112 | 112 | 0 | 0 | - | 0.4 |
| 9 essential | float | 112 | 52 | 0 | 60 | - | 0.0 |

| method | n | correct | wrong | unresolved | time (s) |
|---|---|---|---|---|---|
| FMSL | 1120 | 1120 | 0 | 0 | 1.8 |
| SymPy | 1120 | 963 | 147 | 10 | 2.9 |
| float | 896 | 557 | 279 | 60 | 0.1 |
