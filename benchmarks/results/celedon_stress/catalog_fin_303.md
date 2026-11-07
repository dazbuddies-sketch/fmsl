# Celedon point error catalog

60 distinct families, logs: fin_303.jsonl

## Counts by mode and class

| mode | class | count |
|---|---|---|
| deep | TIMEOUT | 2 |
| default | MISSED_CHANGE | 1 |
| default | MISSED_FLAGGED | 11 |
| default | MISSED_SILENT | 2 |
| default | PHANTOM | 11 |
| default | TIMEOUT | 2 |

## By source

- **dgnp**: deep: honest 'cannot prove' branch 4
- **hidden1**: deep: TIMEOUT 2, default: MISSED_CHANGE 1, default: MISSED_FLAGGED 3, default: MISSED_SILENT 2, default: PHANTOM 4, default: TIMEOUT 2
- **hidden2**: default: MISSED_FLAGGED 8, default: PHANTOM 7
- **irrat**: deep: honest 'cannot prove' branch 12

## Every error

| source | mode | class | t | family | detail | loop-confirmed |
|---|---|---|---|---|---|---|
| hidden1 | default | MISSED_CHANGE | -6/13 | `(x^3*y)/((y-(2)*x)^2 + (t-(0)+1*(2)^3)*x^4 - 1*x*y^3)` | type differs between t=-6/13 and t=6/13, inside reported interval (None, None) |  |
| hidden1 | default | MISSED_FLAGGED | 2 | `(x^4-x*y^3)/((y-(1-t)*x)^2 + (t-(2)+1*(1-t)^3)*x^4 - 1*x*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden1 | default | MISSED_FLAGGED | 1 | `(x^3*y)/((y-(t+1)*x)^2 + (t-(1)+1*(t+1)^3)*x^4 - 1*x*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden1 | default | MISSED_FLAGGED | -3/2 | `(x^4)/((y-(t)*x)^2 + (t-(-3/2)+3*(t)^3)*x^4 - 3*x*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | 1 | `(x^5*y)/((y-x-(-1)*x^2)^2 + (t-(1)+2)*x^6 - 2*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | 1/2 | `(x^7)/((y-x-(-1)*x^2)^2 + (t-(1/2)+1)*x^6 - 1*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | 0 | `(x^7)/((y-x-(1)*x^2)^2 + (t-(0)+2)*x^6 - 2*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | 1 | `(x^5*y)/((y-x-(2)*x^2)^2 + (t-(1)+2)*x^6 - 2*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | -1 | `(x^7)/((y-x-(2)*x^2)^2 + (t-(-1)+1)*x^6 - 1*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | 1 | `(x^6)/((y-x-(1)*x^2)^2 + (t-(1)+2)*x^6 - 2*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | 1 | `(x^5*y)/((y-x-(1)*x^2)^2 + (t-(1)+1)*x^6 - 1*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | 1/2 | `(x^6)/((y-x-(-1)*x^2)^2 + (t-(1/2)+1)*x^6 - 1*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden1 | default | MISSED_SILENT | 0 | `(x^3*y)/((y-(2)*x)^2 + (t-(0)+1*(2)^3)*x^4 - 1*x*y^3)` | real change not reported and no warning given | True |
| hidden1 | default | MISSED_SILENT | -3/2 | `(x^4-x*y^3)/((y-(t+1)*x)^2 + (t-(-3/2)+1*(t+1)^3)*x^4 - 1*x*y^3)` | real change not reported and no warning given | True |
| hidden1 | default | PHANTOM | ~-1.4238538 | `(x^4-x*y^3)/((y-(t+1)*x)^2 + (t-(-3/2)+1*(t+1)^3)*x^4 - 1*x*y^3)` | top-level candidate; nearest real change at -1.5, offset +0.0761462 | True |
| hidden1 | default | PHANTOM | 1 | `(x^4-x*y^3)/((y-(1-t)*x)^2 + (t-(2)+1*(1-t)^3)*x^4 - 1*x*y^3)` | top-level candidate; nearest real change at 0, offset +1 | True |
| hidden1 | default | PHANTOM | 0 | `(x^3*y)/((y-(t+1)*x)^2 + (t-(1)+1*(t+1)^3)*x^4 - 1*x*y^3)` | top-level candidate; nearest real change at -1, offset +1 | True |
| hidden1 | default | PHANTOM | ~-0.6554042 | `(x^4)/((y-(t)*x)^2 + (t-(-3/2)+3*(t)^3)*x^4 - 3*x*y^3)` | top-level candidate; nearest real change at -1.5, offset +0.844596 | True |
| hidden2 | default | PHANTOM | -1 | `(x^5*y)/((y-x-(-1)*x^2)^2 + (t-(1)+2)*x^6 - 2*x^3*y^3)` | top-level candidate; nearest real change at 1, offset -2 | True |
| hidden2 | default | PHANTOM | -1/2 | `(x^7)/((y-x-(-1)*x^2)^2 + (t-(1/2)+1)*x^6 - 1*x^3*y^3)` | top-level candidate; nearest real change at 0.5, offset -1 | True |
| hidden2 | default | PHANTOM | -1 | `(x^5*y)/((y-x-(2)*x^2)^2 + (t-(1)+2)*x^6 - 2*x^3*y^3)` | top-level candidate; nearest real change at 1, offset -2 | True |
| hidden2 | default | PHANTOM | -2 | `(x^7)/((y-x-(2)*x^2)^2 + (t-(-1)+1)*x^6 - 1*x^3*y^3)` | top-level candidate; nearest real change at -1, offset -1 | True |
| hidden2 | default | PHANTOM | -1 | `(x^6)/((y-x-(1)*x^2)^2 + (t-(1)+2)*x^6 - 2*x^3*y^3)` | top-level candidate; nearest real change at 1, offset -2 | True |
| hidden2 | default | PHANTOM | 0 | `(x^5*y)/((y-x-(1)*x^2)^2 + (t-(1)+1)*x^6 - 1*x^3*y^3)` | top-level candidate; nearest real change at 1, offset -1 | True |
| hidden2 | default | PHANTOM | -1/2 | `(x^6)/((y-x-(-1)*x^2)^2 + (t-(1/2)+1)*x^6 - 1*x^3*y^3)` | top-level candidate; nearest real change at 0.5, offset -1 | True |
| hidden1 | deep | TIMEOUT | None | `(x^4-x*y^3)/((y-(t/2)*x)^2 + (t-(0)+3*(t/2)^3)*x^4 - 3*x*y^3)` | >240s |  |
| hidden1 | deep | TIMEOUT | None | `(x^4-x*y^3)/((y-(t/2)*x)^2 + (t-(1/2)+2*(t/2)^3)*x^4 - 2*x*y^3)` | >240s |  |
| hidden1 | default | TIMEOUT | None | `(x^4-x*y^3)/((y-(t/2)*x)^2 + (t-(0)+3*(t/2)^3)*x^4 - 3*x*y^3)` | >240s |  |
| hidden1 | default | TIMEOUT | None | `(x^4-x*y^3)/((y-(t/2)*x)^2 + (t-(1/2)+2*(t/2)^3)*x^4 - 2*x*y^3)` | >240s |  |
