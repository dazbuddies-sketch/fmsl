# Celedon point error catalog

67 distinct families, logs: m_101.jsonl

## Counts by mode and class

| mode | class | count |
|---|---|---|
| deep | PHANTOM | 1 |
| deep | WRONG_END_FORMULA | 8 |
| default | MISSED_CHANGE | 6 |
| default | MISSED_FLAGGED | 13 |
| default | MISSED_SILENT | 5 |
| default | PHANTOM | 15 |
| default | WRONG_END_FORMULA | 6 |

## By source

- **dgnp**: deep: honest 'cannot prove' branch 4
- **hidden1**: default: MISSED_CHANGE 1, default: MISSED_FLAGGED 5, default: MISSED_SILENT 2, default: PHANTOM 5
- **hidden2**: default: MISSED_FLAGGED 8, default: PHANTOM 6
- **irrat**: deep: honest 'cannot prove' branch 21
- **merge**: deep: PHANTOM 1, deep: WRONG_END_FORMULA 8, deep: honest 'cannot prove' branch 5, default: MISSED_CHANGE 5, default: MISSED_SILENT 3, default: PHANTOM 4, default: WRONG_END_FORMULA 6
- **path**: deep: honest 'cannot prove' branch 2
- **randA**: deep: honest 'cannot prove' branch 9

## Every error

| source | mode | class | t | family | detail | loop-confirmed |
|---|---|---|---|---|---|---|
| hidden1 | default | MISSED_CHANGE | 1/13 | `(x^3*y)/((y-(1-t)*x)^2 + (t-(1/2)+1*(1-t)^3)*x^4 - 1*x*y^3)` | type differs between t=1/13 and t=7/13, inside reported interval (None, 1) |  |
| merge | default | MISSED_CHANGE | -18/13 | `x^4/((y-(1)*x)^2-1*x^4) + 2*t*x^2/(x^2+y^2)` | type differs between t=-18/13 and t=-12/13, inside reported interval (None, 0) |  |
| merge | default | MISSED_CHANGE | 25/26 | `x^4/((y-(1)*x)^2-1*x^4) + 2*t*x^2/(x^2+y^2)` | type differs between t=25/26 and t=37/26, inside reported interval (1/2, None) |  |
| merge | default | MISSED_CHANGE | 6/13 | `x^4/((y-(1)*x)^2-2*x^4) + t*x*y/(x^2+y^2)` | type differs between t=6/13 and t=12/13, inside reported interval (0, None) |  |
| merge | default | MISSED_CHANGE | -66/13 | `x^4/((y-(2)*x)^2-2*x^4) + t*x*y/(x^2+y^2)` | type differs between t=-66/13 and t=-60/13, inside reported interval (None, 0) |  |
| merge | default | MISSED_CHANGE | 6/13 | `x^4/((y-(2)*x)^2-2*x^4) + t*x*y/(x^2+y^2)` | type differs between t=6/13 and t=12/13, inside reported interval (0, None) |  |
| hidden1 | default | MISSED_FLAGGED | 1 | `(x^4)/((y-(t)*x)^2 + (t-(1)+1*(t)^3)*x^4 - 1*x*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden1 | default | MISSED_FLAGGED | -3/2 | `(x^4)/((y-(2*t)*x)^2 + (t-(-3/2)+2*(2*t)^3)*x^4 - 2*x*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden1 | default | MISSED_FLAGGED | 1 | `(x^3*y)/((y-(t+1)*x)^2 + (t-(1)+3*(t+1)^3)*x^4 - 3*x*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden1 | default | MISSED_FLAGGED | 1 | `(x^5)/((y-(1)*x)^2 + (t-(1)+1*(1)^3)*x^4 - 1*x*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden1 | default | MISSED_FLAGGED | -3/2 | `(x^3*y)/((y-(1-t)*x)^2 + (t-(-3/2)+3*(1-t)^3)*x^4 - 3*x*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | 1/2 | `(x^6)/((y-x-(-1)*x^2)^2 + (t-(1/2)+1)*x^6 - 1*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | 1 | `(x^5*y)/((y-x-(2)*x^2)^2 + (t-(1)+1)*x^6 - 1*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | 1 | `(x^5*y)/((y-x-(2)*x^2)^2 + (t-(1)+2)*x^6 - 2*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | 0 | `(x^7)/((y-x-(2)*x^2)^2 + (t-(0)+2)*x^6 - 2*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | -1 | `(x^7)/((y-x-(2)*x^2)^2 + (t-(-1)+1)*x^6 - 1*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | -1 | `(x^6)/((y-x-(-1)*x^2)^2 + (t-(-1)+2)*x^6 - 2*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | -1 | `(x^7)/((y-x-(1)*x^2)^2 + (t-(-1)+2)*x^6 - 2*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | 0 | `(x^6)/((y-x-(2)*x^2)^2 + (t-(0)+1)*x^6 - 1*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden1 | default | MISSED_SILENT | 1/2 | `(x^3*y)/((y-(1-t)*x)^2 + (t-(1/2)+1*(1-t)^3)*x^4 - 1*x*y^3)` | real change not reported and no warning given | True |
| hidden1 | default | MISSED_SILENT | 1/2 | `(x^4+x^2*y^2)/((y-(t)*x)^2 + (t-(1/2)+3*(t)^3)*x^4 - 3*x*y^3)` | real change not reported and no warning given | True |
| merge | default | MISSED_SILENT | -1 | `x^4/((y-(1)*x)^2-1*x^4) + 2*t*x^2/(x^2+y^2)` | real change not reported and no warning given | True |
| merge | default | MISSED_SILENT | 1 | `x^4/((y-(1)*x)^2-1*x^4) + 2*t*x^2/(x^2+y^2)` | real change not reported and no warning given | True |
| merge | default | MISSED_SILENT | 1/2 | `x^4/((y-(1)*x)^2-2*x^4) + t*x*y/(x^2+y^2)` | real change not reported and no warning given | True |
| merge | deep | PHANTOM | 5/8 | `x^4/((y-(2)*x)^2-2*x^4) + t*x*y/(x^2+y^2)` | deep candidate; no real change anywhere reported | True |
| hidden1 | default | PHANTOM | ~0.6823278 | `(x^4)/((y-(t)*x)^2 + (t-(1)+1*(t)^3)*x^4 - 1*x*y^3)` | top-level candidate; nearest real change at 1, offset -0.317672 | True |
| hidden1 | default | PHANTOM | ~-0.40859132 | `(x^4)/((y-(2*t)*x)^2 + (t-(-3/2)+2*(2*t)^3)*x^4 - 2*x*y^3)` | top-level candidate; nearest real change at -1.5, offset +1.09141 | True |
| hidden1 | default | PHANTOM | ~-0.25258475 | `(x^3*y)/((y-(t+1)*x)^2 + (t-(1)+3*(t+1)^3)*x^4 - 3*x*y^3)` | top-level candidate; nearest real change at -1, offset +0.747415 | True |
| hidden1 | default | PHANTOM | ~0.36001477 | `(x^4+x^2*y^2)/((y-(t)*x)^2 + (t-(1/2)+3*(t)^3)*x^4 - 3*x*y^3)` | top-level candidate; nearest real change at 0.5, offset -0.139985 | True |
| hidden1 | default | PHANTOM | 0 | `(x^5)/((y-(1)*x)^2 + (t-(1)+1*(1)^3)*x^4 - 1*x*y^3)` | top-level candidate; nearest real change at 1, offset -1 | True |
| hidden2 | default | PHANTOM | -1/2 | `(x^6)/((y-x-(-1)*x^2)^2 + (t-(1/2)+1)*x^6 - 1*x^3*y^3)` | top-level candidate; nearest real change at 0.5, offset -1 | True |
| hidden2 | default | PHANTOM | 0 | `(x^5*y)/((y-x-(2)*x^2)^2 + (t-(1)+1)*x^6 - 1*x^3*y^3)` | top-level candidate; nearest real change at 1, offset -1 | True |
| hidden2 | default | PHANTOM | -1 | `(x^5*y)/((y-x-(2)*x^2)^2 + (t-(1)+2)*x^6 - 2*x^3*y^3)` | top-level candidate; nearest real change at 1, offset -2 | True |
| hidden2 | default | PHANTOM | -2 | `(x^7)/((y-x-(2)*x^2)^2 + (t-(-1)+1)*x^6 - 1*x^3*y^3)` | top-level candidate; nearest real change at -1, offset -1 | True |
| hidden2 | default | PHANTOM | -3 | `(x^6)/((y-x-(-1)*x^2)^2 + (t-(-1)+2)*x^6 - 2*x^3*y^3)` | top-level candidate; nearest real change at -1, offset -2 | True |
| hidden2 | default | PHANTOM | -1 | `(x^6)/((y-x-(2)*x^2)^2 + (t-(0)+1)*x^6 - 1*x^3*y^3)` | top-level candidate; nearest real change at 0, offset -1 | True |
| merge | default | PHANTOM | 0 | `x^4/((y-(1)*x)^2-1*x^4) + 2*t*x^2/(x^2+y^2)` | top-level candidate; nearest real change at -1, offset +1 | True |
| merge | default | PHANTOM | 1/2 | `x^4/((y-(1)*x)^2-1*x^4) + 2*t*x^2/(x^2+y^2)` | top-level candidate; nearest real change at 1, offset -0.5 | True |
| merge | default | PHANTOM | 0 | `x^4/((y-(1)*x)^2-2*x^4) + t*x*y/(x^2+y^2)` | top-level candidate; nearest real change at 0.5, offset -0.5 | True |
| merge | default | PHANTOM | 0 | `x^4/((y-(2)*x)^2-2*x^4) + t*x*y/(x^2+y^2)` | top-level candidate; no real change anywhere reported | True |
| merge | deep | WRONG_END_FORMULA | -6/5 | `x^4/((y-(-1)*x)^2+1*x^4) + t*y^2/(x^2+y^2)` | end #1 = 2/5 does not satisfy {(0, 1): Fraction(1, 1)} (formula from interval (None, None)) |  |
| merge | deep | WRONG_END_FORMULA | 6/5 | `x^4/((y-(-1)*x)^2+1*x^4) + t*y^2/(x^2+y^2)` | end #0 = 0 does not satisfy {(0, 1): Fraction(1, 1), (1, 0): Fraction(-1, 1)} (formula from interval (None, None)) |  |
| merge | deep | WRONG_END_FORMULA | 1/5 | `x^4/((y-(1)*x)^2-1*x^4) + 2*t*x^2/(x^2+y^2)` | end #1 = 0 does not satisfy {(0, 1): Fraction(1, 1), (1, 0): Fraction(-2, 1)} (formula from interval (-1, 1)) |  |
| merge | deep | WRONG_END_FORMULA | 6/5 | `x^4/((y-(2)*x)^2+2*x^4) + t*y^2/(x^2+y^2)` | end #0 = 0 does not satisfy {(0, 1): Fraction(1, 1), (1, 0): Fraction(-1, 1)} (formula from interval (None, None)) |  |
| merge | deep | WRONG_END_FORMULA | 6/5 | `x^4/((y-(2)*x)^2+2*x^4) + t*y^2/(x^2+y^2)` | end #1 = 73/50 does not satisfy {(0, 1): Fraction(1, 1)} (formula from interval (None, None)) |  |
| merge | deep | WRONG_END_FORMULA | 6/5 | `x^4/((y-(1)*x)^2+1*x^4) + (1-t)*x^2/(x^2+y^2)` | end #0 = -1/5 does not satisfy {(0, 1): Fraction(1, 1)} (formula from interval (None, None)) |  |
| merge | deep | WRONG_END_FORMULA | 6/5 | `x^4/((y-(1)*x)^2+1*x^4) + (1-t)*x^2/(x^2+y^2)` | end #1 = 9/10 does not satisfy {(0, 1): Fraction(1, 1), (0, 0): Fraction(-1, 1), (1, 0): Fraction(1, 1)} (formula from interval (None, None)) |  |
| merge | deep | WRONG_END_FORMULA | -6/5 | `x^4/((y-(1)*x)^2+2*x^4) + (t^2-1)*y^2/(x^2+y^2)` | end #1 = 18/25 does not satisfy {(0, 1): Fraction(1, 1), (0, 0): Fraction(1, 1), (2, 0): Fraction(-1, 1)} (formula from interval (None, None)) |  |
| merge | default | WRONG_END_FORMULA | -6/5 | `x^4/((y-(-1)*x)^2+1*x^4) + t*y^2/(x^2+y^2)` | end #1 = 2/5 does not satisfy {(0, 1): Fraction(1, 1)} (formula from interval (None, None)) |  |
| merge | default | WRONG_END_FORMULA | 6/5 | `x^4/((y-(-1)*x)^2+1*x^4) + t*y^2/(x^2+y^2)` | end #0 = 0 does not satisfy {(0, 1): Fraction(1, 1), (1, 0): Fraction(-1, 1)} (formula from interval (None, None)) |  |
| merge | default | WRONG_END_FORMULA | 6/5 | `x^4/((y-(2)*x)^2+2*x^4) + t*y^2/(x^2+y^2)` | end #0 = 0 does not satisfy {(0, 1): Fraction(1, 1), (1, 0): Fraction(-1, 1)} (formula from interval (None, None)) |  |
| merge | default | WRONG_END_FORMULA | 6/5 | `x^4/((y-(2)*x)^2+2*x^4) + t*y^2/(x^2+y^2)` | end #1 = 73/50 does not satisfy {(0, 1): Fraction(1, 1)} (formula from interval (None, None)) |  |
| merge | default | WRONG_END_FORMULA | 6/5 | `x^4/((y-(1)*x)^2+1*x^4) + (1-t)*x^2/(x^2+y^2)` | end #0 = -1/5 does not satisfy {(0, 1): Fraction(1, 1)} (formula from interval (None, None)) |  |
| merge | default | WRONG_END_FORMULA | -6/5 | `x^4/((y-(1)*x)^2+2*x^4) + (t^2-1)*y^2/(x^2+y^2)` | end #1 = 18/25 does not satisfy {(0, 1): Fraction(1, 1), (0, 0): Fraction(1, 1), (2, 0): Fraction(-1, 1)} (formula from interval (None, None)) |  |
