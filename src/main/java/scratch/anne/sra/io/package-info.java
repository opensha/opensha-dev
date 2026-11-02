/**
 * Reading inputs and writing results.
 *
 * <p>{@code readers} parses portfolios, vulnerability and fragility
 * libraries, run-parameter configuration, and intermediate rupture-loss
 * data from CSV, JSON, Parquet, and XML. {@code writers} exports hazard
 * curves, per-rupture and per-branch results, loss-exceedance curves, and
 * ensemble summaries. File-level provenance metadata is propagated through
 * these readers/writers so results remain traceable to the run that
 * produced them.</p>
 */
package scratch.anne.sra.io;
