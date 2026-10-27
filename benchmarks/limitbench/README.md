# LimitBench

Limits of functions of several variables, with **proven** answers: a test of whether a solver,
human or AI, reasons about a limit or pattern-matches.

Each item gives a rational function f of two or three variables and asks for the set of all
values f approaches at the origin: the limit if it exists, and otherwise the complete closed set
of limit values (for example [-1/2, 1/2], or [-oo, -3] u [3, +oo]).

**Ground truth.** Every answer is confirmed by two independent methods that must agree: either a
construction with a known closed form (AM-GM bounds, eigenvalues of quadratic forms) or FMSL's
exact engine, and FMSL's independent exact image of f on small squares (r = 10^-16). Items where
they disagree, or where FMSL is incomplete, are not included.

**Traps.** Each two-variable item records what the common "try a few paths" method concludes
(limits along y = 0, x = 0, y = ±x, y = 2x, y = ±x², x = y²) and whether that conclusion is
wrong. Example: x³y/(x⁶ + y²), where all eight paths give 0 but the values fill [-1/2, 1/2].

**Tiers.** 1: classic (AM-GM ratios, quadratic forms). 2: hidden curves (values found only along
curves like y = c·x^k). 3: random rational functions. 4: three variables.

## Files

- `limitbench.jsonl`: the items (`f`, `variables`, `limit_exists`, `limit`, `fluid_state`,
  `fluid_state_float`, `path_test`, `tier`, `category`, `verified_by`).
- `generate.py`: rebuilds the dataset deterministically.
- `grade.py`: the prompt and the grader. Answers end with `ANSWER: LIMIT <value>` or
  `ANSWER: SET <set>`; scored for the right verdict (does the limit exist?) and the exact answer.
- `evaluate.py`: runs the path method, FMSL (must score 100%), or a Claude model through the API.
- `results/`: scores; `.github/workflows/limitbench.yml` runs models in GitHub Actions.

## Running models

Add an `ANTHROPIC_API_KEY` repository secret (Settings -> Secrets and variables -> Actions), then
commit `benchmarks/limitbench/RUN.json`, for example `{"models": ["claude-sonnet-5-5"], "limit": null}`.
Model runs use the key's account and are billed to it.
