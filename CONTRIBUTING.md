# Contributing to xastk

Thanks for helping. This project aims for analysis code that a referee can trust, so the bar is
correctness and traceability first.

## Workflow

1. Open an issue describing the problem or method, with a reference if it is a published method.
2. Branch from `main`, add or change code **with tests**. A new numerical method needs a test
   against an analytical result, a synthetic dataset with known truth, or a published value.
3. Run `ruff check .` and `pytest` locally.
4. Open a pull request. CI runs on Linux, macOS and Windows with Python 3.11–3.13.

## Code rules

- The core (`src/xastk`) is MIT. Do not import GPL packages from the core. GPL integrations
  (for example pyFitIt) belong in a separately distributed plugin.
- Do not copy code from repositories without a license. Re-implement from the paper and cite it.
- Heavy dependencies (Larch, emcee, python-docx, originpro) are imported inside the functions that
  need them, never at module import time.
- Any use of a private Larch function must be pinned by a golden-value test in
  `tests/test_larch_integration.py`.
- Numbers quoted in `docs/` must be produced by a test or by a script in `examples/`.

## AI-assistance disclosure

AI tools may be used, under these rules (modelled on the
[Larch AI policy](https://github.com/xraypy/xraylarch/blob/master/AI_POLICY.md)):

- **Disclose** in the pull request which tool was used and for what (code, tests, docs).
- **Understand and own** every line you submit; you must be able to explain it in review.
- **No autonomous pull requests.** A human reviews and submits.
- **Upstream contributions** to Larch or other projects must follow *their* policies. Larch
  requires disclosure and does not accept autonomous AI pull requests.

The initial v0.1.0.dev0 scaffold, code, tests and technical report were drafted with Claude
(Anthropic) at the maintainer's request and are under maintainer review.

## Commit messages

Imperative mood, short subject line, body explaining *why*. Reference issues with `#123`.
