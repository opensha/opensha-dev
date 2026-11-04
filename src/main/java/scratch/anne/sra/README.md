# scratch.anne.sra — Seismic Risk Analysis

Regional portfolio seismic-risk assessment within the OpenSHA ecosystem,
with direct integration into the full hazard logic tree and support for
time-dependence.

The library incorporates OpenSHA hazard with structural-response models
(vulnerability or fragility) to produce expected loss and loss-exceedance
results, from single portfolio runs up to full logic-tree ensembles.

## Entry points

- `engine.lec_pipeline.PipelineMain` — the full multi-stage true-mean
  logic-tree LEC workflow (per-rupture EALs → branch contributions →
  pseudo loss curves → LEC summary). Each stage is gated by output
  existence, so re-runs only recompute what is missing.
- `examples.risk_convolution.PortfolioExpectedLoss` — a single
  expected-loss portfolio run (FULL_HCURVE or PER_RUPTURE), driving
  `engine.portfolio_workflow`.

See the `examples` package for smaller, single-concept usage samples.

## Package map

Indentation shows package nesting: an indented entry is a subpackage of
the entry above it.

```
sra/
├── domain/                     immutable data model
│   ├── asset/                  portfolio assets
│   │   ├── vulnerability/
│   │   └── fragility/
│   ├── hazard/                 hazard parameters, curves, ERF building, loss-exceedance curves
│   ├── portfolio/              portfolios
│   │   └── portfolio_wrappers/ risk-convolution wrappers around a portfolio
│   └── structural_response/    response-model libraries
│       ├── fragilities/
│       └── vulnerabilities/
├── engine/                     calculators that operate on the domain model
│   ├── convolution/            the core risk-convolution integral (Riemann / closed-form)
│   ├── portfolio_workflow/     portfolio-level expected-loss calculators
│   ├── per_rupture/            per-rupture loss and probability accumulation
│   ├── logic_tree/             per-branch EAL and LEC-distribution building
│   └── lec_pipeline/           orchestrated 4-stage LEC workflow (and PipelineMain)
├── io/                         reading inputs and writing results
│   ├── readers/                CSV, JSON, Parquet, XML inputs
│   └── writers/                curves, per-rupture/per-branch results, summaries
├── util/                       numeric, parsing, and grouping helpers
│   └── enums/                  shared enums (IM type, limit state, fragility type, COV distribution)
├── examples/                   runnable usage samples and research drivers
│   ├── io/                     single-concept reader examples
│   ├── risk_convolution/       single-concept calculation examples
│   ├── PBR/                    precarious-balanced-rock research drivers
│   └── workflow_tests/         sensitivity, timing, and round-trip harnesses
└── gui/                        interactive UI (PBRAppUI)
```
