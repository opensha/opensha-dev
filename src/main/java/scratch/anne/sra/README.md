# scratch.anne.sra — Seismic Risk Analysis

Regional portfolio seismic-risk assessment within the OpenSHA ecosystem,
with direct integration into the full hazard logic tree and support for
time-dependence.

The library incorporates OpenSHA hazard with structural-response models
(vulnerability or fragility) to produce expected loss and loss-exceedance
results, from single portfolio runs up to full logic-tree ensembles.

## Entry points

- `engine.lecPipeline.PipelineMain` — the full multi-stage true-mean
  logic-tree LEC workflow (per-rupture EALs → branch contributions →
  pseudo loss curves → LEC summary). Each stage is gated by output
  existence, so re-runs only recompute what is missing.
- `examples.riskConvolution.PortfolioExpectedLossExample` — a single
  expected-loss portfolio run (FULL_HCURVE or PER_RUPTURE), driving
  `engine.portfolioWorkflow`.

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
│   │   └── portfolioWrappers/ risk-convolution wrappers around a portfolio
│   └── structuralResponse/    response-model libraries
│       ├── fragilities/
│       └── vulnerabilities/
├── engine/                     calculators that operate on the domain model
│   ├── convolution/            the core risk-convolution integral (Riemann / closed-form)
│   ├── portfolioWorkflow/     portfolio-level expected-loss calculators
│   ├── perRupture/            per-rupture loss and probability accumulation
│   ├── logicTree/             per-branch EAL and LEC-distribution building
│   └── lecPipeline/           orchestrated 4-stage LEC workflow (and PipelineMain)
├── io/                         reading inputs and writing results
│   ├── readers/                CSV, JSON, Parquet, XML inputs
│   └── writers/                curves, per-rupture/per-branch results, summaries
├── util/                       numeric, parsing, and grouping helpers
│   └── enums/                  shared enums (IM type, limit state, fragility type, COV distribution)
├── examples/                   runnable usage samples and research drivers
│   ├── io/                     single-concept reader examples
│   ├── riskConvolution/       single-concept calculation examples
│   ├── PBR/                    precarious-balanced-rock research drivers
│   └── workflowTests/         sensitivity, timing, and round-trip harnesses
└── gui/                        interactive UI (PBRAppUI)
```
