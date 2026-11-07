# Celedon point error catalog

70 distinct families, logs: v_404.jsonl

## Counts by mode and class

| mode | class | count |
|---|---|---|
| default | MISSED_CHANGE | 7 |
| default | MISSED_FLAGGED | 20 |
| default | MISSED_SILENT | 11 |
| default | PHANTOM | 28 |

## By source

- **hidden1**: default: MISSED_CHANGE 2, default: MISSED_FLAGGED 5, default: MISSED_SILENT 2, default: PHANTOM 5
- **hidden2**: default: MISSED_FLAGGED 8, default: PHANTOM 8
- **irrat**: default: MISSED_CHANGE 1, default: MISSED_FLAGGED 2, default: MISSED_SILENT 3, default: PHANTOM 4
- **merge**: default: MISSED_CHANGE 4, default: MISSED_FLAGGED 5, default: MISSED_SILENT 6, default: PHANTOM 11

## Every error

| source | mode | class | t | family | detail | loop-confirmed |
|---|---|---|---|---|---|---|
| hidden1 | default | MISSED_CHANGE | 6/13 | `(x^3*y)/((y-(2)*x)^2 + (t-(1)+1*(2)^3)*x^4 - 1*x*y^3)` | type differs between t=6/13 and t=18/13, inside reported interval (None, None) |  |
| hidden1 | default | MISSED_CHANGE | 6/13 | `(x^5)/((y-(2)*x)^2 + (t-(1/2)+1*(2)^3)*x^4 - 1*x*y^3)` | type differs between t=6/13 and t=18/13, inside reported interval (None, None) |  |
| irrat | default | MISSED_CHANGE | 6/13 | `(x^8-x^6*y^2)/((y^2-3*x^2)^2 + (t-(1)+3^2)*x^8 - x^4*y^4)` | type differs between t=6/13 and t=18/13, inside reported interval (None, None) |  |
| merge | default | MISSED_CHANGE | -18/13 | `x^4/((y-(-1)*x)^2-1*x^4) + 2*t*x^2/(x^2+y^2)` | type differs between t=-18/13 and t=-12/13, inside reported interval (None, 0) |  |
| merge | default | MISSED_CHANGE | 25/26 | `x^4/((y-(-1)*x)^2-1*x^4) + 2*t*x^2/(x^2+y^2)` | type differs between t=25/26 and t=37/26, inside reported interval (1/2, None) |  |
| merge | default | MISSED_CHANGE | -5/13 | `x^4/((y-(1)*x)^2-2*x^4) + (t-1)*y^2/(x^2+y^2)` | type differs between t=-5/13 and t=1/13, inside reported interval (None, 1) |  |
| merge | default | MISSED_CHANGE | 51/26 | `x^4/((y-(1)*x)^2-2*x^4) + (t-1)*y^2/(x^2+y^2)` | type differs between t=51/26 and t=63/26, inside reported interval (3/2, None) |  |
| hidden1 | default | MISSED_FLAGGED | 1 | `(x^5)/((y-(t/2)*x)^2 + (t-(1)+3*(t/2)^3)*x^4 - 3*x*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden1 | default | MISSED_FLAGGED | 1/2 | `(x^4-x*y^3)/((y-(1-t)*x)^2 + (t-(1/2)+3*(1-t)^3)*x^4 - 3*x*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden1 | default | MISSED_FLAGGED | 0 | `(x^4+x^2*y^2)/((y-(1-t)*x)^2 + (t-(0)+1*(1-t)^3)*x^4 - 1*x*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden1 | default | MISSED_FLAGGED | -3/2 | `(x^4)/((y-(1-t)*x)^2 + (t-(-3/2)+1*(1-t)^3)*x^4 - 1*x*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden1 | default | MISSED_FLAGGED | 1 | `(x^4+x^2*y^2)/((y-(1)*x)^2 + (t-(1)+1*(1)^3)*x^4 - 1*x*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | -1 | `(x^5*y)/((y-x-(-1)*x^2)^2 + (t-(-1)+2)*x^6 - 2*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | 0 | `(x^5*y)/((y-x-(1)*x^2)^2 + (t-(0)+2)*x^6 - 2*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | 0 | `(x^7)/((y-x-(2)*x^2)^2 + (t-(0)+1)*x^6 - 1*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | 1 | `(x^6)/((y-x-(1)*x^2)^2 + (t-(1)+2)*x^6 - 2*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | 1 | `(x^6)/((y-x-(1)*x^2)^2 + (t-(1)+1)*x^6 - 1*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | -1 | `(x^7)/((y-x-(2)*x^2)^2 + (t-(-1)+1)*x^6 - 1*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | 1 | `(x^7)/((y-x-(-1)*x^2)^2 + (t-(1)+1)*x^6 - 1*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| hidden2 | default | MISSED_FLAGGED | 1 | `(x^5*y)/((y-x-(1)*x^2)^2 + (t-(1)+1)*x^6 - 1*x^3*y^3)` | real change not reported (a 'change near' warning was given) | True |
| irrat | default | MISSED_FLAGGED | 1 | `(x^8-x^6*y^2)/((y^2-(t+2)*x^2)^2 + (t-(1)+(t+2)^2)*x^8 - x^4*y^4)` | real change not reported (a 'change near' warning was given) | True |
| irrat | default | MISSED_FLAGGED | 1 | `(x^8)/((y^2-(2*t)*x^2)^2 + (t-(1)+(2*t)^2)*x^8 - x^4*y^4)` | real change not reported (a 'change near' warning was given) | True |
| merge | default | MISSED_FLAGGED | ~-1.4142136 | `x^4/((y-(-1)*x)^2-2*x^4) + (t^2-1)*x^2/(x^2+y^2)` | real change not reported (a 'change near' warning was given) | True |
| merge | default | MISSED_FLAGGED | 0 | `x^4/((y-(-1)*x)^2-2*x^4) + (t^2-1)*x^2/(x^2+y^2)` | real change not reported (a 'change near' warning was given) | True |
| merge | default | MISSED_FLAGGED | ~1.4142136 | `x^4/((y-(-1)*x)^2-2*x^4) + (t^2-1)*x^2/(x^2+y^2)` | real change not reported (a 'change near' warning was given) | True |
| merge | default | MISSED_FLAGGED | -5/2 | `x^4/((y-(2)*x)^2-2*x^4) + 2*t*x*y/(x^2+y^2)` | real change not reported (a 'change near' warning was given) | True |
| merge | default | MISSED_FLAGGED | 5/18 | `x^4/((y-(2)*x)^2-2*x^4) + 2*t*x*y/(x^2+y^2)` | real change not reported (a 'change near' warning was given) | True |
| hidden1 | default | MISSED_SILENT | 1 | `(x^3*y)/((y-(2)*x)^2 + (t-(1)+1*(2)^3)*x^4 - 1*x*y^3)` | real change not reported and no warning given | True |
| hidden1 | default | MISSED_SILENT | 1/2 | `(x^5)/((y-(2)*x)^2 + (t-(1/2)+1*(2)^3)*x^4 - 1*x*y^3)` | real change not reported and no warning given | True |
| irrat | default | MISSED_SILENT | 1/2 | `(x^8+t*x^4*y^4)/((y^2-t*x^2)^2 + (t-(1/2)+t^2)*x^8 - x^4*y^4)` | real change not reported and no warning given | True |
| irrat | default | MISSED_SILENT | 1 | `(x^8-x^6*y^2)/((y^2-3*x^2)^2 + (t-(1)+3^2)*x^8 - x^4*y^4)` | real change not reported and no warning given | True |
| irrat | default | MISSED_SILENT | 1/2 | `(x^8)/((y^2-t*x^2)^2 + (t-(1/2)+t^2)*x^8 - x^4*y^4)` | real change not reported and no warning given | True |
| merge | default | MISSED_SILENT | -1 | `x^4/((y-(-1)*x)^2-1*x^4) + 2*t*x^2/(x^2+y^2)` | real change not reported and no warning given | True |
| merge | default | MISSED_SILENT | 1 | `x^4/((y-(-1)*x)^2-1*x^4) + 2*t*x^2/(x^2+y^2)` | real change not reported and no warning given | True |
| merge | default | MISSED_SILENT | 0 | `x^4/((y-(1)*x)^2-2*x^4) + (t-1)*y^2/(x^2+y^2)` | real change not reported and no warning given | True |
| merge | default | MISSED_SILENT | 2 | `x^4/((y-(1)*x)^2-2*x^4) + (t-1)*y^2/(x^2+y^2)` | real change not reported and no warning given | True |
| merge | default | MISSED_SILENT | ~-1.2747549 | `x^4/((y-(2)*x)^2-2*x^4) + (t^2-1)*y^2/(x^2+y^2)` | real change not reported and no warning given | True |
| merge | default | MISSED_SILENT | ~1.2747549 | `x^4/((y-(2)*x)^2-2*x^4) + (t^2-1)*y^2/(x^2+y^2)` | real change not reported and no warning given | True |
| hidden1 | default | PHANTOM | ~0.80463988 | `(x^5)/((y-(t/2)*x)^2 + (t-(1)+3*(t/2)^3)*x^4 - 3*x*y^3)` | top-level candidate; nearest real change at 1, offset -0.19536 | True |
| hidden1 | default | PHANTOM | 1 | `(x^4-x*y^3)/((y-(1-t)*x)^2 + (t-(1/2)+3*(1-t)^3)*x^4 - 3*x*y^3)` | top-level candidate; nearest real change at 0.5, offset +0.5 | True |
| hidden1 | default | PHANTOM | ~0.20629947 | `(x^4-x*y^3)/((y-(1-t)*x)^2 + (t-(1/2)+3*(1-t)^3)*x^4 - 3*x*y^3)` | top-level candidate; nearest real change at 0, offset +0.206299 | True |
| hidden1 | default | PHANTOM | 1 | `(x^4+x^2*y^2)/((y-(1-t)*x)^2 + (t-(0)+1*(1-t)^3)*x^4 - 1*x*y^3)` | top-level candidate; nearest real change at 0, offset +1 | True |
| hidden1 | default | PHANTOM | 0 | `(x^4+x^2*y^2)/((y-(1)*x)^2 + (t-(1)+1*(1)^3)*x^4 - 1*x*y^3)` | top-level candidate; nearest real change at 1, offset -1 | True |
| hidden2 | default | PHANTOM | -3 | `(x^5*y)/((y-x-(-1)*x^2)^2 + (t-(-1)+2)*x^6 - 2*x^3*y^3)` | top-level candidate; nearest real change at -1, offset -2 | True |
| hidden2 | default | PHANTOM | -2 | `(x^5*y)/((y-x-(1)*x^2)^2 + (t-(0)+2)*x^6 - 2*x^3*y^3)` | top-level candidate; nearest real change at 0, offset -2 | True |
| hidden2 | default | PHANTOM | -1 | `(x^7)/((y-x-(2)*x^2)^2 + (t-(0)+1)*x^6 - 1*x^3*y^3)` | top-level candidate; nearest real change at 0, offset -1 | True |
| hidden2 | default | PHANTOM | -1 | `(x^6)/((y-x-(1)*x^2)^2 + (t-(1)+2)*x^6 - 2*x^3*y^3)` | top-level candidate; nearest real change at 1, offset -2 | True |
| hidden2 | default | PHANTOM | 0 | `(x^6)/((y-x-(1)*x^2)^2 + (t-(1)+1)*x^6 - 1*x^3*y^3)` | top-level candidate; nearest real change at 1, offset -1 | True |
| hidden2 | default | PHANTOM | -2 | `(x^7)/((y-x-(2)*x^2)^2 + (t-(-1)+1)*x^6 - 1*x^3*y^3)` | top-level candidate; nearest real change at -1, offset -1 | True |
| hidden2 | default | PHANTOM | 0 | `(x^7)/((y-x-(-1)*x^2)^2 + (t-(1)+1)*x^6 - 1*x^3*y^3)` | top-level candidate; nearest real change at 1, offset -1 | True |
| hidden2 | default | PHANTOM | 0 | `(x^5*y)/((y-x-(1)*x^2)^2 + (t-(1)+1)*x^6 - 1*x^3*y^3)` | top-level candidate; nearest real change at 1, offset -1 | True |
| irrat | default | PHANTOM | 0 | `(x^8-x^6*y^2)/((y^2-(t+2)*x^2)^2 + (t-(1)+(t+2)^2)*x^8 - x^4*y^4)` | top-level candidate; nearest real change at -1, offset +1 | True |
| irrat | default | PHANTOM | ~0.3660254 | `(x^8+t*x^4*y^4)/((y^2-t*x^2)^2 + (t-(1/2)+t^2)*x^8 - x^4*y^4)` | top-level candidate; nearest real change at 0.5, offset -0.133975 | True |
| irrat | default | PHANTOM | ~0.3903882 | `(x^8)/((y^2-(2*t)*x^2)^2 + (t-(1)+(2*t)^2)*x^8 - x^4*y^4)` | top-level candidate; nearest real change at 0, offset +0.390388 | True |
| irrat | default | PHANTOM | ~0.3660254 | `(x^8)/((y^2-t*x^2)^2 + (t-(1/2)+t^2)*x^8 - x^4*y^4)` | top-level candidate; nearest real change at 0.5, offset -0.133975 | True |
| merge | default | PHANTOM | 1 | `x^4/((y-(-1)*x)^2-2*x^4) + (t^2-1)*x^2/(x^2+y^2)` | top-level candidate; nearest real change at 1.41421, offset -0.414214 | True |
| merge | default | PHANTOM | ~-1.2247449 | `x^4/((y-(-1)*x)^2-2*x^4) + (t^2-1)*x^2/(x^2+y^2)` | top-level candidate; nearest real change at -1.41421, offset +0.189469 | True |
| merge | default | PHANTOM | ~1.2247449 | `x^4/((y-(-1)*x)^2-2*x^4) + (t^2-1)*x^2/(x^2+y^2)` | top-level candidate; nearest real change at 1.41421, offset -0.189469 | True |
| merge | default | PHANTOM | -1 | `x^4/((y-(-1)*x)^2-2*x^4) + (t^2-1)*x^2/(x^2+y^2)` | top-level candidate; nearest real change at -1.41421, offset +0.414214 | True |
| merge | default | PHANTOM | 0 | `x^4/((y-(-1)*x)^2-1*x^4) + 2*t*x^2/(x^2+y^2)` | top-level candidate; nearest real change at -1, offset +1 | True |
| merge | default | PHANTOM | 1/2 | `x^4/((y-(-1)*x)^2-1*x^4) + 2*t*x^2/(x^2+y^2)` | top-level candidate; nearest real change at 1, offset -0.5 | True |
| merge | default | PHANTOM | 0 | `x^4/((y-(2)*x)^2-2*x^4) + 2*t*x*y/(x^2+y^2)` | top-level candidate; nearest real change at 0.277778, offset -0.277778 | True |
| merge | default | PHANTOM | 1 | `x^4/((y-(1)*x)^2-2*x^4) + (t-1)*y^2/(x^2+y^2)` | top-level candidate; nearest real change at 0, offset +1 | True |
| merge | default | PHANTOM | 3/2 | `x^4/((y-(1)*x)^2-2*x^4) + (t-1)*y^2/(x^2+y^2)` | top-level candidate; nearest real change at 2, offset -0.5 | True |
| merge | default | PHANTOM | ~-1.2247449 | `x^4/((y-(2)*x)^2-2*x^4) + (t^2-1)*y^2/(x^2+y^2)` | top-level candidate; nearest real change at -1.27475, offset +0.05001 | True |
| merge | default | PHANTOM | ~1.2247449 | `x^4/((y-(2)*x)^2-2*x^4) + (t^2-1)*y^2/(x^2+y^2)` | top-level candidate; nearest real change at 1.27475, offset -0.05001 | True |
