# Celedon point error catalog

60 distinct families, logs: pre_101.jsonl

## Counts by mode and class

| mode | class | count |
|---|---|---|
| deep | PHANTOM | 7 |
| default | MISSED_CHANGE | 1 |
| default | MISSED_FLAGGED | 13 |
| default | MISSED_SILENT | 2 |
| default | PHANTOM | 14 |

## By source

- **dgnp**: deep: honest 'cannot prove' branch 4
- **hidden1**: deep: PHANTOM 7, default: MISSED_CHANGE 1, default: MISSED_FLAGGED 5, default: MISSED_SILENT 2, default: PHANTOM 8
- **hidden2**: default: MISSED_FLAGGED 8, default: PHANTOM 6
- **irrat**: deep: honest 'cannot prove' branch 12

## Every error

| source | mode | class | t | family | detail | loop-confirmed |
|---|---|---|---|---|---|---|
| hidden1 | default | MISSED_CHANGE | 1/13 | `(x^3*y)/((y-(1-t)*x)^2 + (t-(1/2)+1*(1-t)^3)*x^4 - 1*x*y^3)` | type differs between t=1/13 and t=7/13, inside reported interval (None, 1) |  |
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
| hidden1 | deep | PHANTOM | ~0.6823278 | `(x^4)/((y-(t)*x)^2 + (t-(1)+1*(t)^3)*x^4 - 1*x*y^3)` | top-level candidate; nearest real change at 1, offset -0.317672 | True |
| hidden1 | deep | PHANTOM | ~-0.40859132 | `(x^4)/((y-(2*t)*x)^2 + (t-(-3/2)+2*(2*t)^3)*x^4 - 2*x*y^3)` | top-level candidate; nearest real change at -1.5, offset +1.09141 | True |
| hidden1 | deep | PHANTOM | ~-0.25258475 | `(x^3*y)/((y-(t+1)*x)^2 + (t-(1)+3*(t+1)^3)*x^4 - 3*x*y^3)` | top-level candidate; nearest real change at -1, offset +0.747415 | True |
| hidden1 | deep | PHANTOM | ~0.39685026 | `(x^4-x*y^3)/((y-(2*t)*x)^2 + (t-(0)+3*(2*t)^3)*x^4 - 3*x*y^3)` | top-level candidate; nearest real change at 0.5, offset -0.10315 | True |
| hidden1 | deep | PHANTOM | ~2.1914879 | `(x^3*y)/((y-(1-t)*x)^2 + (t-(1/2)+1*(1-t)^3)*x^4 - 1*x*y^3)` | top-level candidate; nearest real change at 1, offset +1.19149 | True |
| hidden1 | deep | PHANTOM | ~0.36001477 | `(x^4+x^2*y^2)/((y-(t)*x)^2 + (t-(1/2)+3*(t)^3)*x^4 - 3*x*y^3)` | top-level candidate; nearest real change at 0.5, offset -0.139985 | True |
| hidden1 | deep | PHANTOM | ~2.058566 | `(x^3*y)/((y-(1-t)*x)^2 + (t-(-3/2)+3*(1-t)^3)*x^4 - 3*x*y^3)` | top-level candidate; nearest real change at 1, offset +1.05857 | True |
| hidden1 | default | PHANTOM | ~0.6823278 | `(x^4)/((y-(t)*x)^2 + (t-(1)+1*(t)^3)*x^4 - 1*x*y^3)` | top-level candidate; nearest real change at 1, offset -0.317672 | True |
| hidden1 | default | PHANTOM | ~-0.40859132 | `(x^4)/((y-(2*t)*x)^2 + (t-(-3/2)+2*(2*t)^3)*x^4 - 2*x*y^3)` | top-level candidate; nearest real change at -1.5, offset +1.09141 | True |
| hidden1 | default | PHANTOM | ~-0.25258475 | `(x^3*y)/((y-(t+1)*x)^2 + (t-(1)+3*(t+1)^3)*x^4 - 3*x*y^3)` | top-level candidate; nearest real change at -1, offset +0.747415 | True |
| hidden1 | default | PHANTOM | ~0.39685026 | `(x^4-x*y^3)/((y-(2*t)*x)^2 + (t-(0)+3*(2*t)^3)*x^4 - 3*x*y^3)` | top-level candidate; nearest real change at 0.5, offset -0.10315 | True |
| hidden1 | default | PHANTOM | ~2.1914879 | `(x^3*y)/((y-(1-t)*x)^2 + (t-(1/2)+1*(1-t)^3)*x^4 - 1*x*y^3)` | top-level candidate; nearest real change at 1, offset +1.19149 | True |
| hidden1 | default | PHANTOM | ~0.36001477 | `(x^4+x^2*y^2)/((y-(t)*x)^2 + (t-(1/2)+3*(t)^3)*x^4 - 3*x*y^3)` | top-level candidate; nearest real change at 0.5, offset -0.139985 | True |
| hidden1 | default | PHANTOM | 0 | `(x^5)/((y-(1)*x)^2 + (t-(1)+1*(1)^3)*x^4 - 1*x*y^3)` | top-level candidate; nearest real change at 1, offset -1 | True |
| hidden1 | default | PHANTOM | ~2.058566 | `(x^3*y)/((y-(1-t)*x)^2 + (t-(-3/2)+3*(1-t)^3)*x^4 - 3*x*y^3)` | top-level candidate; nearest real change at 1, offset +1.05857 | True |
| hidden2 | default | PHANTOM | -1/2 | `(x^6)/((y-x-(-1)*x^2)^2 + (t-(1/2)+1)*x^6 - 1*x^3*y^3)` | top-level candidate; nearest real change at 0.5, offset -1 | True |
| hidden2 | default | PHANTOM | 0 | `(x^5*y)/((y-x-(2)*x^2)^2 + (t-(1)+1)*x^6 - 1*x^3*y^3)` | top-level candidate; nearest real change at 1, offset -1 | True |
| hidden2 | default | PHANTOM | -1 | `(x^5*y)/((y-x-(2)*x^2)^2 + (t-(1)+2)*x^6 - 2*x^3*y^3)` | top-level candidate; nearest real change at 1, offset -2 | True |
| hidden2 | default | PHANTOM | -2 | `(x^7)/((y-x-(2)*x^2)^2 + (t-(-1)+1)*x^6 - 1*x^3*y^3)` | top-level candidate; nearest real change at -1, offset -1 | True |
| hidden2 | default | PHANTOM | -3 | `(x^6)/((y-x-(-1)*x^2)^2 + (t-(-1)+2)*x^6 - 2*x^3*y^3)` | top-level candidate; nearest real change at -1, offset -2 | True |
| hidden2 | default | PHANTOM | -1 | `(x^6)/((y-x-(2)*x^2)^2 + (t-(0)+1)*x^6 - 1*x^3*y^3)` | top-level candidate; nearest real change at 0, offset -1 | True |
