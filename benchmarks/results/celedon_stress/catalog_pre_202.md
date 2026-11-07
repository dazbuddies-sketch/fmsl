# Celedon point error catalog

57 distinct families, logs: pre_202.jsonl

## Counts by mode and class

| mode | class | count |
|---|---|---|
| deep | PHANTOM | 7 |
| deep | TIMEOUT | 2 |
| default | MISSED_CHANGE | 1 |
| default | MISSED_FLAGGED | 13 |
| default | MISSED_SILENT | 1 |
| default | PHANTOM | 13 |
| default | TIMEOUT | 2 |

## By source

- **dgnp**: deep: honest 'cannot prove' branch 2
- **hidden1**: deep: PHANTOM 7, deep: TIMEOUT 1, default: MISSED_CHANGE 1, default: MISSED_FLAGGED 5, default: MISSED_SILENT 1, default: PHANTOM 7, default: TIMEOUT 1
- **hidden2**: default: MISSED_FLAGGED 8, default: PHANTOM 6
- **irrat**: deep: honest 'cannot prove' branch 8
- **randB**: deep: TIMEOUT 1, default: TIMEOUT 1

## Every error

| source | mode | class | t | family | detail | loop-confirmed |
|---|---|---|---|---|---|---|
| hidden1 | default | MISSED_CHANGE | -30/13 | `(x^4-x*y^3)/((y-(-1)*x)^2 + (t-(-3/2)+3*(-1)^3)*x^4 - 3*x*y^3)` | type differs between t=-30/13 and t=-18/13, inside reported interval (None, None) |  |
| hidden1 | default | MISSED_FLAGGED | 1 | `(x^4-x*y^3)/((y-(t+1)*x)^2 + (t-(1)+3*(t+1)^3)*x^4 - 3*x*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden1 | default | MISSED_FLAGGED | 2 | `(x^3*y)/((y-(2*t)*x)^2 + (t-(2)+2*(2*t)^3)*x^4 - 2*x*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden1 | default | MISSED_FLAGGED | -1 | `(x^4-x*y^3)/((y-(1)*x)^2 + (t-(-1)+2*(1)^3)*x^4 - 2*x*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden1 | default | MISSED_FLAGGED | 2 | `(x^4)/((y-(2*t)*x)^2 + (t-(2)+2*(2*t)^3)*x^4 - 2*x*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden1 | default | MISSED_FLAGGED | -3/2 | `(x^4)/((y-(2*t)*x)^2 + (t-(-3/2)+3*(2*t)^3)*x^4 - 3*x*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | 0 | `(x^7)/((y-x-(-1)*x^2)^2 + (t-(0)+2)*x^6 - 2*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | 1 | `(x^5*y)/((y-x-(-1)*x^2)^2 + (t-(1)+1)*x^6 - 1*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | 1/2 | `(x^5*y)/((y-x-(-1)*x^2)^2 + (t-(1/2)+1)*x^6 - 1*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | 0 | `(x^6)/((y-x-(-1)*x^2)^2 + (t-(0)+2)*x^6 - 2*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | 0 | `(x^5*y)/((y-x-(2)*x^2)^2 + (t-(0)+1)*x^6 - 1*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | -1 | `(x^7)/((y-x-(2)*x^2)^2 + (t-(-1)+2)*x^6 - 2*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | -1 | `(x^5*y)/((y-x-(1)*x^2)^2 + (t-(-1)+1)*x^6 - 1*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | 1 | `(x^6)/((y-x-(1)*x^2)^2 + (t-(1)+2)*x^6 - 2*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden1 | default | MISSED_SILENT | -3/2 | `(x^4-x*y^3)/((y-(-1)*x)^2 + (t-(-3/2)+3*(-1)^3)*x^4 - 3*x*y^3)` | real change not reported and no warning given | True |
| hidden1 | deep | PHANTOM | ~-0.20629947 | `(x^4-x*y^3)/((y-(t+1)*x)^2 + (t-(1)+3*(t+1)^3)*x^4 - 3*x*y^3)` | top-level candidate; nearest real change at 0, offset -0.206299 | True |
| hidden1 | deep | PHANTOM | ~-0.25258475 | `(x^4-x*y^3)/((y-(t+1)*x)^2 + (t-(1)+3*(t+1)^3)*x^4 - 3*x*y^3)` | top-level candidate; nearest real change at 0, offset -0.252585 | True |
| hidden1 | deep | PHANTOM | ~1.541834 | `(x^4-x*y^3)/((y-(t/2)*x)^2 + (t-(2)+1*(t/2)^3)*x^4 - 1*x*y^3)` | top-level candidate; nearest real change at 2, offset -0.458166 | True |
| hidden1 | deep | PHANTOM | ~1.5874011 | `(x^4-x*y^3)/((y-(t/2)*x)^2 + (t-(2)+1*(t/2)^3)*x^4 - 1*x*y^3)` | top-level candidate; nearest real change at 2, offset -0.412599 | True |
| hidden1 | deep | PHANTOM | ~0.45843774 | `(x^3*y)/((y-(2*t)*x)^2 + (t-(2)+2*(2*t)^3)*x^4 - 2*x*y^3)` | top-level candidate; nearest real change at 0, offset +0.458438 | True |
| hidden1 | deep | PHANTOM | ~0.45843774 | `(x^4)/((y-(2*t)*x)^2 + (t-(2)+2*(2*t)^3)*x^4 - 2*x*y^3)` | top-level candidate; nearest real change at 2, offset -1.54156 | True |
| hidden1 | deep | PHANTOM | ~-0.36195109 | `(x^4)/((y-(2*t)*x)^2 + (t-(-3/2)+3*(2*t)^3)*x^4 - 3*x*y^3)` | top-level candidate; nearest real change at -1.5, offset +1.13805 | True |
| hidden1 | default | PHANTOM | ~-0.20629947 | `(x^4-x*y^3)/((y-(t+1)*x)^2 + (t-(1)+3*(t+1)^3)*x^4 - 3*x*y^3)` | top-level candidate; nearest real change at 0, offset -0.206299 | True |
| hidden1 | default | PHANTOM | ~-0.25258475 | `(x^4-x*y^3)/((y-(t+1)*x)^2 + (t-(1)+3*(t+1)^3)*x^4 - 3*x*y^3)` | top-level candidate; nearest real change at 0, offset -0.252585 | True |
| hidden1 | default | PHANTOM | ~1.541834 | `(x^4-x*y^3)/((y-(t/2)*x)^2 + (t-(2)+1*(t/2)^3)*x^4 - 1*x*y^3)` | top-level candidate; nearest real change at 2, offset -0.458166 | True |
| hidden1 | default | PHANTOM | ~1.5874011 | `(x^4-x*y^3)/((y-(t/2)*x)^2 + (t-(2)+1*(t/2)^3)*x^4 - 1*x*y^3)` | top-level candidate; nearest real change at 2, offset -0.412599 | True |
| hidden1 | default | PHANTOM | ~0.45843774 | `(x^3*y)/((y-(2*t)*x)^2 + (t-(2)+2*(2*t)^3)*x^4 - 2*x*y^3)` | top-level candidate; nearest real change at 0, offset +0.458438 | True |
| hidden1 | default | PHANTOM | ~0.45843774 | `(x^4)/((y-(2*t)*x)^2 + (t-(2)+2*(2*t)^3)*x^4 - 2*x*y^3)` | top-level candidate; nearest real change at 2, offset -1.54156 | True |
| hidden1 | default | PHANTOM | ~-0.36195109 | `(x^4)/((y-(2*t)*x)^2 + (t-(-3/2)+3*(2*t)^3)*x^4 - 3*x*y^3)` | top-level candidate; nearest real change at -1.5, offset +1.13805 | True |
| hidden2 | default | PHANTOM | 0 | `(x^5*y)/((y-x-(-1)*x^2)^2 + (t-(1)+1)*x^6 - 1*x^3*y^3)` | top-level candidate; nearest real change at 1, offset -1 | True |
| hidden2 | default | PHANTOM | -1/2 | `(x^5*y)/((y-x-(-1)*x^2)^2 + (t-(1/2)+1)*x^6 - 1*x^3*y^3)` | top-level candidate; nearest real change at 0.5, offset -1 | True |
| hidden2 | default | PHANTOM | -2 | `(x^6)/((y-x-(-1)*x^2)^2 + (t-(0)+2)*x^6 - 2*x^3*y^3)` | top-level candidate; nearest real change at 0, offset -2 | True |
| hidden2 | default | PHANTOM | -1 | `(x^5*y)/((y-x-(2)*x^2)^2 + (t-(0)+1)*x^6 - 1*x^3*y^3)` | top-level candidate; nearest real change at 0, offset -1 | True |
| hidden2 | default | PHANTOM | -2 | `(x^5*y)/((y-x-(1)*x^2)^2 + (t-(-1)+1)*x^6 - 1*x^3*y^3)` | top-level candidate; nearest real change at -1, offset -1 | True |
| hidden2 | default | PHANTOM | -1 | `(x^6)/((y-x-(1)*x^2)^2 + (t-(1)+2)*x^6 - 2*x^3*y^3)` | top-level candidate; nearest real change at 1, offset -2 | True |
| hidden1 | deep | TIMEOUT | None | `(x^4-x*y^3)/((y-(t/2)*x)^2 + (t-(0)+3*(t/2)^3)*x^4 - 3*x*y^3)` | >240s |  |
| randB | deep | TIMEOUT | None | `((t^2+3*t+-2)*x^1*y^3 + (1)*x^1*y^2 + (-1*t+-2)*x^2*y^2)/((t^2+-3*t+-3)*x^0*y^5 + (-2)*x^2*y^3 + (2)*t*x^0*y^5)` | >240s |  |
| hidden1 | default | TIMEOUT | None | `(x^4-x*y^3)/((y-(t/2)*x)^2 + (t-(0)+3*(t/2)^3)*x^4 - 3*x*y^3)` | >240s |  |
| randB | default | TIMEOUT | None | `((t^2+3*t+-2)*x^1*y^3 + (1)*x^1*y^2 + (-1*t+-2)*x^2*y^2)/((t^2+-3*t+-3)*x^0*y^5 + (-2)*x^2*y^3 + (2)*t*x^0*y^5)` | >240s |  |
