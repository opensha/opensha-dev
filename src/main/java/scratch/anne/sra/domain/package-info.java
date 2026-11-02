/**
 * Immutable data model for seismic risk analysis.
 *
 * <p>These types describe <em>what</em> is being analyzed, independent of
 * any calculation. Subpackages cover the portfolio of exposed assets
 * ({@code asset}, {@code portfolio}), the seismic hazard they are exposed
 * to ({@code hazard}), and the structural-response models — vulnerability
 * or fragility — that map intensity to loss or damage-state exceedance
 * ({@code structural_response}).</p>
 *
 * <p>Model objects are immutable value types; the calculators that consume
 * them live in {@link scratch.anne.sra.engine}.</p>
 */
package scratch.anne.sra.domain;
