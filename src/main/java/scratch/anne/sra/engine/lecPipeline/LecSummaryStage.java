package scratch.anne.sra.engine.lecPipeline;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import scratch.anne.sra.domain.hazard.RuptureKey;
import scratch.anne.sra.domain.risk.LossExceedanceCurve;
import scratch.anne.sra.engine.logicTree.BranchLecDistributionBuilder;
import scratch.anne.sra.engine.logicTree.StreamingBranchLecBuilder;
import scratch.anne.sra.io.readers.ParquetMetadataReader;
import scratch.anne.sra.io.readers.RuptureLossLookupReader;
import scratch.anne.sra.io.writers.LogicTreeBranchLecExporter;
import scratch.anne.sra.io.writers.LogicTreeBranchPercentileSummaryExporter;

/**
 * Pipeline stage 4: per-branch LEC summary CSV + ensemble-level
 * percentile summary CSV, streamed from stage 2's branch contributions.
 *
 * <p>Replicates {@code LogicTreeLecSummary}'s body. No existing class is
 * modified; this class only calls the same public readers/writers the
 * original script already used.</p>
 *
 * <p><b>COV mode is a placeholder only.</b> {@link PipelineConfig#getCovMode()}
 * currently has no effect on behavior here beyond the output file names
 * (via {@link PipelineConfig#covModeTag()}) -- this stage, like the
 * original {@code LogicTreeLecSummary}, is hardcoded to
 * {@link StreamingBranchLecBuilder#streamDiscreteRankCurves}, i.e.
 * {@code LecAccumulation.discreteRank}. There is currently no
 * COV-producing counterpart wired into {@code StreamingBranchLecBuilder}
 * at all ({@code LecAccumulation.logNormalCOV} exists but nothing calls
 * it). If/when COV-mode support is added, THIS is the one place that
 * would need to change -- the call to {@code streamDiscreteRankCurves}
 * below would need to branch on {@code config.getCovMode()} and call a
 * COV-producing stream method instead. Nothing else in this pipeline
 * (stages 1-3) would need to change or be rerun: they are entirely
 * COV-agnostic, and stage 4's own output file names already fold
 * {@code covMode} in via {@link PipelineConfig#covModeTag()}, so merely
 * changing {@code covMode} naturally produces new output file names and
 * re-triggers only this stage (via the existence check below), without
 * needing to touch any other stage's gating.</p>
 *
 * <p>Gated by existence: if both
 * {@link PipelineConfig#getBranchEalsWithLecCsv} and
 * {@link PipelineConfig#getBranchPercentilesCsv} already exist, this
 * stage is skipped -- no SolutionLogicTree load or Parquet streaming
 * needed just to check.</p>
 */
public final class LecSummaryStage {

    private LecSummaryStage() {}

    public static void run(PipelineConfig config) throws IOException, SQLException {

        Path branchSummaryCsv = config.getBranchEalsWithLecCsv();
        Path branchPercentilesCsv = config.getBranchPercentilesCsv();

        if (Files.exists(branchSummaryCsv) && Files.exists(branchPercentilesCsv)) {
            System.out.println("[Stage 4] Skipping -- already exist: " + branchSummaryCsv + ", " + branchPercentilesCsv);
            return;
        }

        System.out.println("[Stage 4] Not found -- running LEC summary (covMode=" + config.getCovMode() + ")...");

        Path exportedCsv = config.getPerRuptureEalCsv();
        Path branchEALsCsv = config.getBranchEalsCsv();
        Path contributionsParquet = config.getBranchContributionsParquet();

        // fixed loss lookup: same one built in stage 2, read fresh here
        // since this stage can run independently of that one
        System.out.println("[Stage 4] Loading fixed loss lookup from: " + exportedCsv);
        Map<RuptureKey, Double> fixedLossLookup = RuptureLossLookupReader.read(exportedCsv);
        System.out.println("[Stage 4] Loaded fixed loss lookup for " + fixedLossLookup.size() + " ruptures");

        // exceedance-rate levels to evaluate each branch's curve at
        double[] rates = StandardRateLikelihoods.get();

        // weighted percentiles to report across the branch ensemble, per level
        double[] percentilesOfInterest = config.getPercentilesOfInterest();

        // branch weights, needed for the percentile pass
        Map<Integer, Double> weightByBranch = readBranchWeights(branchEALsCsv);
        List<BranchLecDistributionBuilder.BranchLevelSample> ealSamples = readEalSamples(branchEALsCsv);

        // ------------------------------------------------------------
        // SINGLE streaming pass over the Parquet file, feeding both the
        // per-branch output and the ensemble percentile output.
        //
        // COV-MODE PLACEHOLDER: this call is hardcoded to
        // streamDiscreteRankCurves (discrete-rank accumulation). See
        // this class's javadoc -- if COV-mode support is ever added,
        // this is the call site that would need to branch on
        // config.getCovMode() and call a COV-producing stream method
        // instead.
        // ------------------------------------------------------------
        Map<Integer, double[]> lossByBranchAndLevel = new HashMap<>();

        Map<Double, List<BranchLecDistributionBuilder.BranchLevelSample>> samplesByLikelihood = new LinkedHashMap<>();
        for (double level : rates)
            samplesByLikelihood.put(level, new ArrayList<>());

        System.out.println("[Stage 4] Streaming branch contributions from: " + contributionsParquet);
        StreamingBranchLecBuilder.streamDiscreteRankCurves(contributionsParquet, fixedLossLookup, (branchIndex, curve) -> {
            List<LossExceedanceCurve.LossAtLikelihood> results = curve.sampleAtLikelihoods(rates);

            Double weight = weightByBranch.get(branchIndex);
            if (weight == null)
                throw new IllegalStateException("No branch weight found for branchIndex " + branchIndex);

            double[] losses = new double[rates.length];
            for (int i = 0; i < results.size(); i++) {
                double loss = results.get(i).loss();
                losses[i] = loss;
                samplesByLikelihood.get(rates[i]).add(
                        new BranchLecDistributionBuilder.BranchLevelSample(branchIndex, weight, loss));
            }
            lossByBranchAndLevel.put(branchIndex, losses);
        });
        System.out.println("[Stage 4] Built LEC-derived losses for " + lossByBranchAndLevel.size() + " branches");

        // build the EAL distribution and every likelihood-level distribution
        BranchLecDistributionBuilder.LevelDistribution ealDistribution =
                BranchLecDistributionBuilder.summarize(BranchLecDistributionBuilder.EAL_LEVEL, ealSamples, percentilesOfInterest);
        Map<Double, BranchLecDistributionBuilder.LevelDistribution> levelDistributions =
                BranchLecDistributionBuilder.summarizeAll(samplesByLikelihood, percentilesOfInterest);

        // ------------------------------------------------------------
        // Write output 1: per-branch CSV
        // ------------------------------------------------------------
        LogicTreeBranchLecExporter.write(branchEALsCsv, rates, lossByBranchAndLevel, branchSummaryCsv);

        // ------------------------------------------------------------
        // Write output 2: ensemble percentile summary
        // ------------------------------------------------------------
        // branch-percentiles.csv has no per-branch rows of its own (one
        // row per METRIC, not per branch), so unlike branch-eals[-with-lec].csv
        // there is no other way to tell how many branches this summary
        // was built from, which solution logic tree zip they came from,
        // or which HazardParameters produced the true-mean ERF they're
        // mapped against -- propagate the SAME provenance
        // BranchContributionsStage already wrote into
        // contributionsParquet's own file-level metadata (see
        // ParquetMetadataReader/LogicTreeBranchRuptureContributionExporter),
        // rather than independently reconstructing a smaller subset of
        // it from branchEALsCsv alone (this stage never loads the
        // SolutionLogicTree -- see this class's own javadoc).
        Map<String, String> contribMetadata = ParquetMetadataReader.readKeyValueMetadata(contributionsParquet);
        String metadataComment = contribMetadata.isEmpty() ? null : formatMetadataComment(contribMetadata);
        LogicTreeBranchPercentileSummaryExporter.write(ealDistribution, levelDistributions, percentilesOfInterest,
                metadataComment, branchPercentilesCsv);

        System.out.println("[Stage 4] DONE");
    }

    /**
     * Flattens a key-value map into one {@code key=val,key2=val2,...}
     * line, for the CSV comment line -- same helper as
     * {@code PseudoLossCurveStage}'s own.
     */
    private static String formatMetadataComment(Map<String, String> metadata) {
        StringBuilder sb = new StringBuilder();
        for (Map.Entry<String, String> e : metadata.entrySet()) {
            if (sb.length() > 0)
                sb.append(',');
            sb.append(e.getKey()).append('=').append(e.getValue());
        }
        return sb.toString();
    }

    private static Map<Integer, Double> readBranchWeights(Path branchEALsCsv) throws IOException {
        List<String> lines = Files.readAllLines(branchEALsCsv);
        int headerIndex = (!lines.isEmpty() && lines.get(0).startsWith("#")) ? 1 : 0;
        String[] header = lines.get(headerIndex).split(",");
        int weightCol = List.of(header).indexOf("branchWeight");
        if (weightCol < 0)
            throw new IllegalStateException("Could not find 'branchWeight' column in header: " + lines.get(headerIndex));

        Map<Integer, Double> weights = new HashMap<>();
        for (int i = headerIndex + 1; i < lines.size(); i++) {
            String[] tokens = lines.get(i).split(",");
            weights.put(Integer.parseInt(tokens[0].trim()), Double.parseDouble(tokens[weightCol].trim()));
        }
        return weights;
    }

    private static List<BranchLecDistributionBuilder.BranchLevelSample> readEalSamples(Path branchEALsCsv) throws IOException {
        List<String> lines = Files.readAllLines(branchEALsCsv);
        int headerIndex = (!lines.isEmpty() && lines.get(0).startsWith("#")) ? 1 : 0;
        String[] header = lines.get(headerIndex).split(",");
        int weightCol = List.of(header).indexOf("branchWeight");
        int ealCol = List.of(header).indexOf("eal");
        if (weightCol < 0 || ealCol < 0)
            throw new IllegalStateException("Could not find 'branchWeight'/'eal' column in header: " + lines.get(headerIndex));

        List<BranchLecDistributionBuilder.BranchLevelSample> samples = new ArrayList<>();
        for (int i = headerIndex + 1; i < lines.size(); i++) {
            String[] tokens = lines.get(i).split(",");
            int branchIndex = Integer.parseInt(tokens[0].trim());
            double weight = Double.parseDouble(tokens[weightCol].trim());
            double eal = Double.parseDouble(tokens[ealCol].trim());
            samples.add(new BranchLecDistributionBuilder.BranchLevelSample(branchIndex, weight, eal));
        }
        return samples;
    }
}
