package scratch.kevin.nshm23;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.SectBValuePlot;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SparseGutenbergRichterSolver;
import org.opensha.sha.magdist.SparseGutenbergRichterSolver.SpreadingMethod;

import com.google.common.base.Preconditions;

public class SparseGRTests {
	
	public static void main(String[] args) throws IOException {
		//FaultSystemRupSet rupSet = FaultSystemRupSet.load(new File("/tmp/rupture_set.zip"));
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(
				new File("/home/kevin/markdown/inversions/fm3_1_u3ref_uniform_reproduce_ucerf3.zip"));

		boolean preserveRates = false;
		SpreadingMethod method = SpreadingMethod.NEAREST_GROUP;
		double sampleDiscr = 0.0;
//		double sampleDiscr = 0.01;

		DefaultXY_DataSet scatter1 = new DefaultXY_DataSet();
		DefaultXY_DataSet scatter2 = new DefaultXY_DataSet();

		DefaultXY_DataSet totRateScatter = new DefaultXY_DataSet();
		DefaultXY_DataSet totMoRateScatter = new DefaultXY_DataSet();

		DefaultXY_DataSet scatter1Avg = new DefaultXY_DataSet();
		DefaultXY_DataSet scatter2Avg = new DefaultXY_DataSet();

		int debugIndex = 100;

		//for (double bValue : new double[] {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4}) {
		for (double bValue : new double[] {0.0, 0.4, 0.8, 1.2}) {
			System.out.println("Doing b="+bValue);
			SupraSeisBValInversionTargetMFDs targets = new SupraSeisBValInversionTargetMFDs.Builder(rupSet, bValue)
					.clearTargetAdjustments().sparseGR(false).build();

			MinMaxAveTracker equivBTrack = new MinMaxAveTracker();
			MinMaxAveTracker equivBTrack2 = new MinMaxAveTracker();

			for (int s=0; s<rupSet.getNumSections(); s++) {
				IncrementalMagFreqDist orig = targets.getOnFaultSupraSeisNucleationMFDs().get(s);

				int minNonZero = -1;
				int maxNonZero = 0;
				for (int i=0; i<orig.size(); i++) {
					if (orig.getY(i) > 0d) {
						if (minNonZero < 0)
							minNonZero = i;
						maxNonZero = i;
					}
				}
				if (minNonZero == maxNonZero)
					continue;
				GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(
						orig.getX(minNonZero), 1+maxNonZero-minNonZero, orig.getDelta(), orig.getTotalMomentRate(), bValue);

				if ((float)gr.getTotalIncrRate() != (float)orig.getTotalIncrRate()) {
					System.err.println("WARNING: rate mismatch on input for "+s+": "
							+(float)gr.getTotalIncrRate()+" != "+(float)orig.getTotalIncrRate());
					for (int i=0; i<orig.size(); i++) {
						double x = orig.getX(i);
						double y1 = orig.getY(i);
						if (gr.hasX(x)) {
							double y2 = gr.getY(i);
							System.out.println(x+"\t"+y1+"\t"+y2);
						} else {
							Preconditions.checkState(y1 == 0d);
						}
					}
				}

				double minMag = gr.getMinX() - 0.5*gr.getDelta();

				double minEncounter = Double.POSITIVE_INFINITY;
				double minAdded = Double.POSITIVE_INFINITY;
				double maxAdded = 0d;
				List<Double> mags = new ArrayList<>();
				for (int r : rupSet.getRupturesForSection(s)) {
					double mag = rupSet.getMagForRup(r);
					minEncounter = Math.min(minEncounter, mag);
					if ((float)mag >= (float)minMag) {
						mags.add(mag);
						minAdded = Math.min(minAdded, mag);
					}
					maxAdded = Math.max(maxAdded, mag);
				}
				//		System.out.println(s+". minEncountered="+(float)minEncounter+"\tminAdded="
				//				+(float)minAdded+"\tfirstBin="+(float)gr.getMinX()+"\tlowerEdge="+(float)minMag);

				IncrementalMagFreqDist sparseGR = SparseGutenbergRichterSolver.getEquivGR(orig, mags, orig.getTotalMomentRate(),
						bValue, sampleDiscr, method, preserveRates);
				//		IncrementalMagFreqDist sparseGR = invertEquivGR(orig, mags, orig.getTotalMomentRate(), bValue);

				double equivB = SectBValuePlot.estBValue(gr.getMinX(), gr.getMaxX(),
						sparseGR.getTotalIncrRate(), sparseGR.getTotalMomentRate());
				scatter1.set(bValue, equivB);
				equivBTrack.addValue(equivB);

				//		equivB = SectBValuePlot.estBValue(minAdded, maxAdded,
				//				sparseGR.getTotalIncrRate(), sparseGR.getTotalMomentRate());
				GutenbergRichterMagFreqDist testGR = new GutenbergRichterMagFreqDist(minAdded, maxAdded, 100);
				testGR.setAllButBvalue(minAdded, maxAdded, sparseGR.getTotalMomentRate(), sparseGR.getTotalIncrRate());
				//		GutenbergRichterMagFreqDist testGR = new GutenbergRichterMagFreqDist(gr.getMinX(), gr.getMaxX(), gr.size());
				//		testGR.setAllButBvalue(gr.getMinX(), gr.getMaxX(), gr.getTotalMomentRate(), gr.getTotalIncrRate());

				double equivB2 = testGR.get_bValue();
				scatter2.set(bValue, equivB2);
				equivBTrack2.addValue(equivB2);

				totRateScatter.set(orig.getTotalIncrRate(), sparseGR.getTotalIncrRate());
				totMoRateScatter.set(orig.getTotalMomentRate(), sparseGR.getTotalMomentRate());

				if (s == debugIndex) {
					List<XY_DataSet> funcs = new ArrayList<>();
					List<PlotCurveCharacterstics> chars = new ArrayList<>();

					funcs.add(gr);
					chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));

					funcs.add(sparseGR);
					chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_SQUARE, 5f, Color.RED));

					MinMaxAveTracker nonZeroTrack = new MinMaxAveTracker();
					for (XY_DataSet func : funcs)
						for (Point2D pt : func)
							if (pt.getY() > 0)
								nonZeroTrack.addValue(pt.getY());

					Range yRange = new Range(Math.min(1e-8, nonZeroTrack.getMin()*0.9), nonZeroTrack.getMax()*1.1);

					DefaultXY_DataSet magFunc = new DefaultXY_DataSet();
					for (double mag : mags)
						magFunc.set(mag, yRange.getLowerBound());
					funcs.add(magFunc);
					chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_SQUARE, 5f, Color.GREEN.darker()));

					double origRate = sampleDiscr == 0d ? gr.calcSumOfY_Vals() : testGR.calcSumOfY_Vals();

					String title = rupSet.getFaultSectionData(s).getSectionName()+" , b="+(float)bValue
							+" , bEst="+(float)equivB+", GR rate="+(float)origRate
							+", new rate="+(float)sparseGR.calcSumOfY_Vals();
					GraphWindow gw = new GraphWindow(funcs, title, chars);
					gw.setYLog(true);
					gw.setX_AxisRange(gr.getMinX()-0.5*gr.getDelta(), gr.getMaxX()+0.5*gr.getDelta());
					gw.setY_AxisRange(yRange);
					gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
				}
			}
			System.out.println("Equivalent b-value distributions (w/ snapped min/max) for b="+(float)bValue+": "+equivBTrack);
			System.out.println("Equivalent b-value distributions (w/ actual min/max) for b="+(float)bValue+": "+equivBTrack2);

			scatter1Avg.set(bValue, equivBTrack.getAverage());
			scatter2Avg.set(bValue, equivBTrack2.getAverage());
		}
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();

		funcs.add(scatter1);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, Color.BLACK));
		funcs.add(scatter1Avg);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 4f, Color.GREEN.darker()));
		GraphWindow gw = new GraphWindow(funcs, "B-Value Scatter 1", chars);
		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);

		funcs.clear();
		chars.clear();
		funcs.add(scatter2);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, Color.BLUE));
		funcs.add(scatter2Avg);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 4f, Color.GREEN.darker()));
		gw = new GraphWindow(funcs, "B-Value Scatter 2", chars);
		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);

		funcs.clear();
		chars.clear();
		funcs.add(totRateScatter);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, Color.GREEN.darker()));
		Range range = calcBoundedLogRange(totRateScatter);
		funcs.add(oneToOne(range));
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
		gw = new GraphWindow(funcs, "Total Rate Scatter", chars);
		gw.setXLog(true);
		gw.setYLog(true);
		gw.setAxisRange(range, range);
		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);

		funcs.clear();
		chars.clear();
		funcs.add(totMoRateScatter);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, Color.RED.darker()));
		range = calcBoundedLogRange(totMoRateScatter);
		funcs.add(oneToOne(range));
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
		gw = new GraphWindow(funcs, "Moment Rate Scatter", chars);
		gw.setXLog(true);
		gw.setYLog(true);
		gw.setAxisRange(range, range);
		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
	}

	private static Range calcBoundedLogRange(XY_DataSet scatter) {
		double min = Math.min(scatter.getMinX(), scatter.getMinY());
		Preconditions.checkState(min > 0);
		double max = Math.max(scatter.getMaxX(), scatter.getMaxY());
		Preconditions.checkState(max > 0);

		return new Range(Math.pow(10, Math.floor(Math.log10(min))), Math.pow(10, Math.ceil(Math.log10(max))));
	}

	private static XY_DataSet oneToOne(Range range) {
		DefaultXY_DataSet line = new DefaultXY_DataSet();
		line.set(range.getLowerBound(), range.getLowerBound());
		line.set(range.getUpperBound(), range.getUpperBound());
		return line;
	}

}
