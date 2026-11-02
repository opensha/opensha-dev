package scratch.anne.sra.engine.logic_tree;

import java.nio.file.Path;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import scratch.anne.sra.domain.hazard.LossExceedanceCurve;
import scratch.anne.sra.domain.hazard.RuptureKey;

/**
 * DuckDB-backed twin of StreamingBranchLecBuilder. Same contract, same
 * contiguity enforcement -- only the row source differs (a vectorized
 * SQL read via DuckDB-JDBC instead of AvroParquetReader/GenericRecord).
 *
 * <p>Requires an explicit ORDER BY branchIndex in the query so rows
 * arrive grouped, matching the same contiguous-by-branch assumption
 * StreamingBranchLecBuilder relies on -- DuckDB does not guarantee row
 * order from read_parquet() without one, even though the file itself
 * was written in branch order.</p>
 */
public final class StreamingBranchLecBuilder {
	
	@FunctionalInterface
	public interface BranchCurveConsumer {
	    void accept(int branchIndex, LossExceedanceCurve curve);
	}

    public static final class NonContiguousBranchException extends RuntimeException {
        public NonContiguousBranchException(int branchIndex) {
            super("branchIndex " + branchIndex + " reappeared after already being flushed -- "
                    + "contributions for this branch are not contiguous in the query result.");
        }
    }

    private StreamingBranchLecBuilder() {}

    public static void streamDiscreteRankCurves(Path contributionsParquet,
                                                  Map<RuptureKey, Double> fixedLossLookup,
                                                  StreamingBranchLecBuilder.BranchCurveConsumer consumer) throws SQLException {

        String sql = "SELECT branchIndex, sourceId, ruptureId, rate "
                + "FROM read_parquet('" + contributionsParquet.toAbsolutePath().normalize().toString().replace("\\", "/") + "') ";

        System.out.println("Starting DuckDB streaming query...");
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
        	     Statement stmt = conn.createStatement();
        	     ResultSet rs = stmt.executeQuery(sql)) {

        	    int currentBranch = -1;
        	    Map<RuptureKey, Double> currentBranchRates = new HashMap<>();
        	    Set<Integer> closedBranches = new HashSet<>();

        	    long startTime = System.nanoTime();
        	    long rowCount = 0;
        	    int branchesCompleted = 0;

        	    while (rs.next()) {
        	        rowCount++;
        	        int branchIndex = rs.getInt("branchIndex");
        	        int sourceId = rs.getInt("sourceId");
        	        int ruptureId = rs.getInt("ruptureId");
        	        double rate = rs.getDouble("rate");

        	        if (branchIndex != currentBranch) {
        	            if (closedBranches.contains(branchIndex))
        	                throw new NonContiguousBranchException(branchIndex);

        	            if (currentBranch >= 0 && !currentBranchRates.isEmpty()) {
        	                LossExceedanceCurve curve = LecAccumulation.discreteRank(currentBranchRates, fixedLossLookup);
        	                consumer.accept(currentBranch, curve);
        	                closedBranches.add(currentBranch);
        	                branchesCompleted++;

        	                if (branchesCompleted % 50 == 0) {
        	                    double elapsedSec = (System.nanoTime() - startTime) / 1e9;
        	                    double rowsPerSec = rowCount / elapsedSec;
        	                    System.out.printf("Processed %d branches, %,d rows, %.1f sec elapsed, %.0f rows/sec%n",
        	                            branchesCompleted, rowCount, elapsedSec, rowsPerSec);
        	                }
        	            }
        	            currentBranch = branchIndex;
        	            currentBranchRates = new HashMap<>();
        	        }

        	        currentBranchRates.put(new RuptureKey(sourceId, ruptureId), rate);
        	    }

        	    if (currentBranch >= 0 && !currentBranchRates.isEmpty()) {
        	        LossExceedanceCurve curve = LecAccumulation.discreteRank(currentBranchRates, fixedLossLookup);
        	        consumer.accept(currentBranch, curve);
        	        closedBranches.add(currentBranch);
        	        branchesCompleted++;
        	    }

        	    double totalElapsedSec = (System.nanoTime() - startTime) / 1e9;
        	    System.out.printf("DONE: %d branches, %,d rows, %.1f sec total, %.0f rows/sec average%n",
        	            branchesCompleted, rowCount, totalElapsedSec, rowCount / totalElapsedSec);
        	}
    }
}