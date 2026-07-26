# Interactive risk simulator — withdrawn July 2026

**Taken down, not patched. Retained as the record of what was served publicly.**

`clinical_sim.html` and `index.html` are byte-identical copies of the same page.
`index.html` was the GitHub Pages landing page at arnava25.github.io, so this
calculator is what a visitor got by default. Both are kept here because the
question "what did the public page actually compute" needs an answer that
survives.

A banner was briefly added above the calculator on 26 July 2026. That was
judged an insufficient mitigation: a warning above a working tool that still
renders a clinical tier badge is not a withdrawal. The page was replaced with a
static notice the same day.

## What it served

The page performed **no external data load** — no `fetch`, no `XMLHttpRequest`,
no `.json` or `.csv` read. Every value was a JavaScript constant compiled into
the file, roughly 270 lines across six blocks:

| Constant | Line | Contents |
|---|---|---|
| `DQTC` | 387 | Bazett ΔQTc for all 66 drug pairs |
| `DAPD` | 456 | genuine ΔAPD90 for all pairs |
| `IKR` | 524 | hERG block percent per pair |
| `TRIPLES` | 629 | 18 three-drug combinations with genuine, Bazett, and tier |
| `PEDIATRIC_CMAX` | 634 | 12 pediatric free Cmax values, hand-annotated with citations |
| `ADULT_REF_PK` | 649 | 12 adult free Cmax values plus dose-scaling parameters |

This constituted **a full second copy of the withdrawn parameter set**, entirely
independent of `data/herg_master_params.csv`. Correcting the CSV would not have
changed a single number on this page. `PEDIATRIC_CMAX` is the same failure mode
as `src/cyp2d6_ito.py`: hand-typed free Cmax values carrying source comments,
living in a file that never reads the table it claims to cite.

## The composite score

Line 720 computed and displayed a clinical tier badge:

```js
let score = 0.50*c_dqtc + 0.20*c_ikr + 0.30*c_faers;
if (hasCyp) score += 15; if (hasCond) score += 10;
```

`scoreLabel()` mapped this to HIGH / MODERATE / LOW-MOD / LOW.

This is the composite metric removed from the manuscript in June 2026 as
circular — `composite_score.py` and `composite_scores.csv` have been in
`archive/src/` and `archive/results/` since then. It nonetheless remained the
driver of the tier shown to any visitor of the live page for a further month.
That gap between "removed from the paper" and "removed from the artifact" is the
substantive lesson here, and is recorded in `params/rebuild_record.md` § 4.6.

## Regeneration

`src/generate_sim_data.py` writes the `DQTC`, `DAPD`, `IKR`, and `TRIPLES`
blocks into this page from `results/risk_grid_results.csv`. It is now marked
DO-NOT-RUN. Any future simulator should read its parameters at runtime from
`params/herg_params_v2.csv` rather than having them compiled in, and must not
reintroduce a composite tier.
