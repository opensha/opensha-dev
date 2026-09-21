package scratch.kevin.sampling;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.junit.Test;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.logicTree.sampling.SamplingMethod;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;

import scratch.kevin.sampling.HazardConvergenceCalcs.*;

public class HazardConvergenceCalcsTest {
	@Test public void spanReferencesUseFixedSizeReserveReplacement() {
		GriddedRegion grid = new GriddedRegion(new Location(0, 0), new Location(0.1, 0.1),
				1d, new Location(0, 0));
		assertEquals(1, grid.getNodeCount());
		RunPeriodData first = run("first", 0, 7);
		RunPeriodData second = run("second", 7, 5);
		RunPeriodData reserve = run("reserve", 100, 4);
		double[][] primaryMaps = java.util.stream.Stream.of(first, second)
				.flatMap(data -> Arrays.stream(data.branchMaps())).toArray(double[][]::new);
		ReferenceStatistics mcs = new ReferenceStatistics(HazardConvergenceCalcs.MCS_REFERENCE_NAME,
				null, 12, HazardConvergenceCalcs.calcHazardStatistics(primaryMaps, 12, new double[] {1}));
		ReferenceStatistics sobol = new ReferenceStatistics(HazardConvergenceCalcs.POOLED_SOBOL_REFERENCE_NAME,
				null, 7, HazardConvergenceCalcs.calcHazardStatistics(first.branchMaps(), 7, new double[] { 1 }));
		List<ReferenceComparison> rows = new ArrayList<>();
		List<Realization> realizations = HazardConvergenceCalcs.appendMCSSpanComparisons(rows,
				List.of(first, second), reserve, new int[] {2, 4}, mcs, sobol, grid, ReturnPeriods.TWO_IN_50);
		assertEquals(11, realizations.size());
		// Seven size-2 realizations and four size-4 realizations, including one reserve prefix of each size.
		assertEquals(135, HazardConvergenceCalcs.buildRealizationPairComparisons(realizations, grid).size());
		// Each primary span has replacement-MCS and Sobol comparisons; each reserve prefix additionally has both.
		assertEquals(110, rows.size());
		for (ReferenceComparison row : rows) {
			assertEquals(0, row.startIndex() % row.sampleCount());
			assertTrue(row.startIndex()+row.sampleCount() <= row.run().maxSamples());
			if (!row.referenceName().equals(HazardConvergenceCalcs.REPLACED_MCS_REFERENCE_NAME))
				continue;
			assertEquals(12, row.referenceSampleCount());
			int from = row.startIndex();
			int to = from+row.sampleCount();
			if (row.metric() == ConvergenceMetric.STANDARD_DEVIATION) {
				double[] included = java.util.stream.IntStream.concat(
						java.util.stream.IntStream.range(0, 12).filter(i -> i < from || i >= to)
								.map(i -> i+1),
						java.util.stream.IntStream.range(0, row.sampleCount()).map(i -> 101+i))
						.asDoubleStream().toArray();
				double mean = Arrays.stream(included).average().orElseThrow();
				double variance = Arrays.stream(included).map(v -> (v-mean)*(v-mean)).average().orElseThrow();
				double spanVariance = (row.sampleCount()*row.sampleCount()-1d)/12d;
				assertEquals(100d*(Math.sqrt(spanVariance/variance)-1), row.comparison().meanPercentChange(), 1e-10);
			}
			if (row.metric() == ConvergenceMetric.MEAN_HAZARD) {
				double spanScale = 0, replacementScale = 0;
				for (int i=0; i<12; i++) {
					if (i >= from && i < to) spanScale += i+1;
					else replacementScale += i+1;
				}
				for (int i=0; i<row.sampleCount(); i++)
					replacementScale += 101+i;
				double test = HazardConvergenceCalcs.buildCurveMeanMap(curves(spanScale), first.curveX(),
						row.sampleCount(), ReturnPeriods.TWO_IN_50)[0];
				double ref = HazardConvergenceCalcs.buildCurveMeanMap(curves(replacementScale), first.curveX(),
						row.referenceSampleCount(), ReturnPeriods.TWO_IN_50)[0];
				assertEquals(100d*(test/ref-1), row.comparison().meanPercentChange(), 1e-10);
			}
		}
	}

	private static RunPeriodData run(String name, int offset, int size) {
		double[][] maps = new double[size][1];
		Map<Integer, double[][]> boundaries = new TreeMap<>();
		double sum = 0;
		for (int i=0; i<size; i++) {
			maps[i][0] = offset+i+1;
			sum += maps[i][0];
			boundaries.put(i+1, curves(sum));
		}
		return new RunPeriodData(new RunSpec(name, null, null, SamplingMethod.MONTE_CARLO, offset, size),
				maps, new double[] {0.1, 1, 10}, curves(sum), Map.of(), boundaries);
	}

	private static double[][] curves(double scale) {
		return new double[][] {{0.01*scale, 0.0001*scale, 0.000001*scale}};
	}

	@Test public void replacementsShareUntouchedRowsAndUseReplacementRows() {
		double[][] rows = {{1}, {2}, {3}, {4}, {5}};
		double[][] replacements = {{8}, {9}};
		double[][] replaced = HazardConvergenceCalcs.replaceSpan(rows, replacements, 1, 3);
		assertSame(rows[0], replaced[0]);
		assertSame(replacements[0], replaced[1]);
		assertSame(replacements[1], replaced[2]);
		assertSame(rows[3], replaced[3]);
		assertSame(rows[4], replaced[4]);
	}

	@Test public void weightedBootstrapStatisticsMatchExpandedSample() {
		double[][] source = {{1}, {2}, {3}, {4}};
		HazardStatistics weighted = HazardConvergenceCalcs.calcBootstrapHazardStatistics(
				source, new int[] {2, 0, 1, 1}, 4);
		double[][] expanded = {{1}, {1}, {3}, {4}};
		HazardStatistics direct = HazardConvergenceCalcs.calcHazardStatistics(expanded, 4, new double[] {0});
		for (ConvergenceMetric metric : ConvergenceMetric.values()) {
			if (metric != ConvergenceMetric.MEAN_HAZARD)
				assertArrayEquals(direct.values(metric), weighted.values(metric), 0d);
		}
	}
}
