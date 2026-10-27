# Contributing to FMSL

Thank you for your interest. FMSL is developed in the open, and contributions of every size help.

## Reporting a problem or a wrong answer
Open an issue with the function you tried (for example `x^2*y/(x^4+y^2)`), what FMSL returned,
and what you expected. A wrong verdict or a wrong fluid state is the most valuable bug report we
can get: every result is supposed to be exact, so a counterexample is always a real bug.

## Proposing a change
1. Fork the repository and create a branch.
2. Add a test in `tests/` that fails before your change and passes after it.
3. Run the full suite: `python -m unittest discover tests`.
4. Open a pull request describing the mathematics behind the change.

## Standards
- Exactness first: no floating-point decisions in the engines. Anything that cannot be decided
  exactly must be reported as UNPROVEN, never guessed.
- New capabilities need an independent cross-check (see `benchmarks/loop_crosscheck.py` and
  `benchmarks/cube_crosscheck.py`).

## Questions
Open a GitHub Discussion or an issue labelled `question`.

## Code of conduct
Be respectful and constructive. Harassment of any kind is not tolerated.
