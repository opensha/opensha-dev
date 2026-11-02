/**
 * Calculators that operate on the {@link scratch.anne.sra.domain} model.
 *
 * <p>The core is {@code convolution}: the risk-convolution integral that
 * combines a hazard curve with a structural-response model to produce
 * expected loss or damage-state exceedance. Higher-level subpackages build
 * on it — {@code portfolio_workflow} runs convolution across a whole
 * portfolio, {@code per_rupture} accumulates loss and probability per
 * rupture, {@code logic_tree} builds per-branch EAL and loss-exceedance
 * distributions, and {@code lec_pipeline} orchestrates the full multi-stage
 * true-mean LEC workflow (see
 * {@link scratch.anne.sra.engine.lec_pipeline.PipelineMain}).</p>
 */
package scratch.anne.sra.engine;
