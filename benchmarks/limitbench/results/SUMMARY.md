### fmsl  (491 items)

| subset | n | verdict right | exact answer | unparsed |
|---|---|---|---|---|
| all | 491 | 100.0% | 100.0% | 0 |
| tier 1 | 184 | 100.0% | 100.0% | 0 |
| tier 2 | 58 | 100.0% | 100.0% | 0 |
| tier 3 | 208 | 100.0% | 100.0% | 0 |
| tier 4 | 41 | 100.0% | 100.0% | 0 |
| path-test traps | 132 | 100.0% | 100.0% | 0 |

### paths  (450 items)

| subset | n | verdict right | exact answer | unparsed |
|---|---|---|---|---|
| all | 450 | 70.7% | 23.1% | 0 |
| tier 1 | 184 | 93.5% | 32.6% | 0 |
| tier 2 | 58 | 82.8% | 63.8% | 0 |
| tier 3 | 208 | 47.1% | 3.4% | 0 |
| path-test traps | 132 | 0.0% | 0.0% | 0 |

### anthropic_claude-haiku-5-5  (50 items)

| subset | n | verdict right | exact answer | unparsed |
|---|---|---|---|---|
| all | 50 | 94.0% | 92.0% | 3 |
| tier 1 | 21 | 85.7% | 81.0% | 3 |
| tier 2 | 4 | 100.0% | 100.0% | 0 |
| tier 3 | 21 | 100.0% | 100.0% | 0 |
| tier 4 | 4 | 100.0% | 100.0% | 0 |
| path-test traps | 17 | 100.0% | 100.0% | 0 |

Tokens used: 16,427 input, 86,271 output.


**Regraded** (grader fix: endpoints written as expressions such as `(1-sqrt(10))/2` were left unparsed;
all three were correct). anthropic_claude-haiku-5-5, same 50 items:

| subset | n | verdict right | exact answer | unparsed |
|---|---|---|---|---|
| all | 50 | 100.0% | 98.0% | 0 |
| tier 1 | 21 | 100.0% | 95.2% | 0 |
| tier 2 | 4 | 100.0% | 100.0% | 0 |
| tier 3 | 21 | 100.0% | 100.0% | 0 |
| tier 4 | 4 | 100.0% | 100.0% | 0 |
| path-test traps | 17 | 100.0% | 100.0% | 0 |

For comparison, the path method on the 46 two-variable items of this subset: verdict 29/46 (63%), exact 10/46 (22%).
