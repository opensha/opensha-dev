package scratch.anne.sra.examples.workflowTests;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import scratch.anne.sra.domain.hazard.RuptureKey;
import scratch.anne.sra.io.readers.BranchRuptureContributionReader;
import scratch.anne.sra.io.readers.BranchRuptureContributionReader.BranchRuptureContributionRow;

/**
 * One-time sanity check: verifies BranchContributionExporter's Parquet
 * output reproduces the same per-branch totals that the inline
 * BranchEALCalculator accumulation produced during the original run.
 *
 * <p>Uses the full row type (branchIndex, sourceId, ruptureId, loss, rate)
 * rather than the identity-stripped LossRateContribution shape, so that
 * any mismatch or duplicate can be traced to a specific rupture rather
 * than just flagged at the branch level.</p>
 *
 * <p>Not a permanent part of the pipeline -- run once after a
 * BranchSolutionFromRuptureResults export to confirm the Parquet file
 * is trustworthy before building any accumulation strategy on top of
 * it.</p>
 */
public class BranchRuptureContributionRoundTripCheck {

    public static void main(String[] args) throws IOException {

        Path baseFolder = Path.of("C:\\Users\\ahulsey\\OneDrive - DOI\\Desktop\\Research\\BERM\\results\\EAL\\true_mean_solution");
        Path inputFolder = Path.of("RES1-W1-MC_Porter-vuln_vs-360");
        String fileTag = inputFolder.getFileName().toString();
        
		// Cap the number of branches processed, for quick tests. Set to Integer.MAX_VALUE to run the full tree.
		int maxBranches = 10;
//		int maxBranches = Integer.MAX_VALUE;
		String fileSuffix = maxBranches != Integer.MAX_VALUE ? "_subset-n-" + maxBranches : "";
		
		Path branchEALsCsv = baseFolder.resolve(inputFolder).resolve(
		        fileTag + "_branch-eals" + fileSuffix + ".csv");
		
		Path branchContributionsParquet = baseFolder.resolve(inputFolder).resolve(fileTag + "_branch-rupture-contributions" + fileSuffix + ".parquet");

        System.out.println("Reading contributions from: " + branchContributionsParquet);
        List<BranchRuptureContributionRow> allRows = BranchRuptureContributionReader.readAll(branchContributionsParquet);

        Map<Integer, List<BranchRuptureContributionRow>> byBranch = new HashMap<>();
        for (BranchRuptureContributionRow row : allRows)
            byBranch.computeIfAbsent(row.branchIndex(), k -> new ArrayList<>()).add(row);
        System.out.println("Loaded " + allRows.size() + " contribution rows across " + byBranch.size() + " branches");

        System.out.println("Reading reference EALs from: " + branchEALsCsv);
        List<String> lines = Files.readAllLines(branchEALsCsv);
        String header = lines.get(0);
        int ealCol = List.of(header.split(",")).indexOf("eal");
        if (ealCol < 0)
            throw new IllegalStateException("Could not find 'eal' column in header: " + header);

        double maxAbsDiff = 0d;
        double maxRelDiff = 0d;
        int checked = 0;
        int branchesWithDuplicates = 0;
        int totalDuplicateRows = 0;

        for (int row = 1; row < lines.size(); row++) {
            String[] tokens = lines.get(row).split(",");
            int branchIndex = Integer.parseInt(tokens[0].trim());
            double referenceEal = Double.parseDouble(tokens[ealCol].trim());

            List<BranchRuptureContributionRow> contributions = byBranch.get(branchIndex);
            double recomputedEal = 0d;

            Set<RuptureKey> seenKeys = new HashSet<>();
            List<RuptureKey> duplicateKeys = new ArrayList<>();

            if (contributions != null) {
                for (BranchRuptureContributionRow c : contributions) {
                    recomputedEal += c.loss() * c.rate();

                    RuptureKey key = new RuptureKey(c.sourceId(), c.ruptureId());
                    if (!seenKeys.add(key))
                        duplicateKeys.add(key);
                }
            }

            double absDiff = Math.abs(recomputedEal - referenceEal);
            double relDiff = referenceEal == 0d ? (recomputedEal == 0d ? 0d : Double.POSITIVE_INFINITY)
                    : absDiff / Math.abs(referenceEal);

            maxAbsDiff = Math.max(maxAbsDiff, absDiff);
            maxRelDiff = Math.max(maxRelDiff, relDiff);
            checked++;

            if (relDiff > 1e-9) {
                System.out.println("MISMATCH branch " + branchIndex
                        + ": reference eal=" + referenceEal
                        + ", recomputed from parquet=" + recomputedEal
                        + " (contributions found: " + (contributions == null ? 0 : contributions.size()) + ")");
            }

            if (!duplicateKeys.isEmpty()) {
                branchesWithDuplicates++;
                totalDuplicateRows += duplicateKeys.size();
                System.out.println("WARNING branch " + branchIndex + ": " + duplicateKeys.size()
                        + " duplicate rupture key(s) in contributions, e.g. " + duplicateKeys.get(0)
                        + " (sourceId=" + duplicateKeys.get(0).sourceId() + ", ruptureId=" + duplicateKeys.get(0).ruptureId() + ")");
            }
        }

        System.out.println("Checked " + checked + " branches");
        System.out.println("Max absolute diff: " + maxAbsDiff);
        System.out.println("Max relative diff: " + maxRelDiff);
        System.out.println("Branches with duplicate rupture keys: " + branchesWithDuplicates
                + " (" + totalDuplicateRows + " duplicate rows total)");
        System.out.println(maxRelDiff < 1e-9 && branchesWithDuplicates == 0
                ? "PASS: parquet contributions reproduce branch EALs exactly, no duplicate rupture keys."
                : "FAIL: see mismatches/warnings above.");
    }
}