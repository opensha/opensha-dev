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
	@Test public void spanReferencesRetainRemaindersAndOtherRuns() {
		GriddedRegion grid = new GriddedRegion(new Location(0, 0), new Location(0.1, 0.1),
				1d, new Location(0, 0));
		assertEquals(1, grid.getNodeCount());
		RunPeriodData first = run("first", 0, 7);
		RunPeriodData second = run("second", 7, 5);
		ReferenceStatistics sobol = new ReferenceStatistics(HazardConvergenceCalcs.POOLED_SOBOL_REFERENCE_NAME,
				null, 7, HazardConvergenceCalcs.calcHazardStatistics(first.branchMaps(), 7, new double[] { 1 }));
		List<ReferenceComparison> rows = new ArrayList<>();
		List<Realization> spans = HazardConvergenceCalcs.appendMCSSpanComparisons(rows, List.of(first, second), new int[] {2, 4},
				sobol, grid, ReturnPeriods.TWO_IN_50);
		assertEquals(7, spans.size());
		// C(5,2) at size 2 plus C(2,2) at size 4, each with five metrics.
		assertEquals(55, HazardConvergenceCalcs.buildRealizationPairComparisons(spans, grid).size());
		// Five two-sample spans and two four-sample spans, each with two references and five metrics.
		assertEquals(70, rows.size());
		for (ReferenceComparison row : rows) {
			assertEquals(0, row.startIndex() % row.sampleCount());
			assertTrue(row.startIndex()+row.sampleCount() <= row.run().maxSamples());
			if (!row.referenceName().equals(HazardConvergenceCalcs.LOO_MCS_REFERENCE_NAME))
				continue;
			assertEquals(12-row.sampleCount(), row.referenceSampleCount());
			int from = (row.run().id().equals("first") ? 0 : 7)+row.startIndex();
			int to = from+row.sampleCount();
			if (row.metric() == ConvergenceMetric.STANDARD_DEVIATION) {
				double[] included = java.util.stream.IntStream.range(0, 12)
						.filter(i -> i < from || i >= to).mapToDouble(i -> i+1).toArray();
				double mean = Arrays.stream(included).average().orElseThrow();
				double variance = Arrays.stream(included).map(v -> (v-mean)*(v-mean)).average().orElseThrow();
				double spanVariance = (row.sampleCount()*row.sampleCount()-1d)/12d;
				assertEquals(100d*(Math.sqrt(spanVariance/variance)-1), row.comparison().meanPercentChange(), 1e-10);
			}
			if (row.metric() == ConvergenceMetric.MEAN_HAZARD) {
				double spanScale = 0, remainingScale = 0;
				for (int i=0; i<12; i++) {
					if (i >= from && i < to) spanScale += i+1;
					else remainingScale += i+1;
				}
				double test = HazardConvergenceCalcs.buildCurveMeanMap(curves(spanScale), first.curveX(),
						row.sampleCount(), ReturnPeriods.TWO_IN_50)[0];
				double ref = HazardConvergenceCalcs.buildCurveMeanMap(curves(remainingScale), first.curveX(),
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

	@Test public void exclusionsShareRowsAndHandleEnds() {
		double[][] rows = {{1}, {2}, {3}, {4}, {5}};
		assertSame(rows[0], HazardConvergenceCalcs.excludeSpan(rows, 1, 3)[0]);
		assertSame(rows[3], HazardConvergenceCalcs.excludeSpan(rows, 1, 3)[1]);
		assertSame(rows[2], HazardConvergenceCalcs.excludeSpan(rows, 0, 2)[0]);
		assertEquals(3, HazardConvergenceCalcs.excludeSpan(rows, 3, 5).length);
	}
}
