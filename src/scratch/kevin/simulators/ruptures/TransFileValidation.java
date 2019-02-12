package scratch.kevin.simulators.ruptures;

import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jfree.data.Range;
import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.util.DataUtils;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.srf.RSQSimEventSlipTimeFunc;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class TransFileValidation {

	public static void main(String[] args) throws IOException {
		RSQSimCatalog catalog;
		if (args.length == 1) {
			File dir = new File(args[0]);
			catalog = new RSQSimCatalog(dir, "Temp", "None", null, "Metadata", null, null);
		} else {
			catalog = Catalogs.BRUCE_2829.instance(new File("/home/kevin/Simulators/catalogs"));
		}
		
		DefaultXY_DataSet momPDiffScatter = new DefaultXY_DataSet();
		DefaultXY_DataSet magPDiffScatter = new DefaultXY_DataSet();
		DefaultXY_DataSet magDiffScatter = new DefaultXY_DataSet();
		
		double minMag = 6.5;
		double printThreshold = 5;
		FileWriter debugFW = new FileWriter(new File(catalog.getCatalogDir(), "trans_file_validation.txt"));
		
		int numZeroInTrans = 0;
		int numZeroInList = 0;
		SummaryStatistics stats = new SummaryStatistics();
		SummaryStatistics absDiffStats = new SummaryStatistics();
		int numAbove0p1 = 0;
		int numAbove1 = 0;
		long testCount = 0;
		
		List<SimulatorElement> allElems = catalog.getElements();
		Map<Integer, Double> patchAreas = new HashMap<>();
		for (SimulatorElement elem : allElems)
			patchAreas.put(elem.getID(), elem.getArea());
		
		SummaryStatistics momentStats = new SummaryStatistics();
		SummaryStatistics momentAbsDiffStats = new SummaryStatistics();
		
		int errCount = 0;
		
		Iterable<RSQSimEvent> eventsIt = catalog.loader().minMag(minMag).iterable();
		int count = 0;
		for (RSQSimEvent event : eventsIt) {
			try {
				if (count % 1000 == 0)
					System.out.println("Processing event: "+count+" (ID="+event.getID()+")");
				count++;
				RSQSimEventSlipTimeFunc slipTimeFunc = catalog.getSlipTimeFunc(event);
				int[] elems = event.getAllElementIDs();
				double[] slips = event.getAllElementSlips();
				
				double listMoment = 0d;
				for (EventRecord rec : event)
					listMoment += rec.getMoment();
				
				double transMoment = 0d;
				
				double time = slipTimeFunc.getEndTime();
				for (int i=0; i<elems.length; i++) {
					testCount++;
					double transSlip = slipTimeFunc.getCumulativeEventSlip(elems[i], time);
					transMoment += FaultMomentCalc.getMoment(patchAreas.get(elems[i]), transSlip);
					if (transSlip == 0 && slips[i] > 0) {
						numZeroInTrans++;
					} else if (transSlip > 0 && slips[i] == 0) {
						numZeroInList++;
					} else {
						double pDiff = DataUtils.getPercentDiff(transSlip, slips[i]);
						stats.addValue(pDiff);
						double diff = Math.abs(transSlip - slips[i]);
						absDiffStats.addValue(diff);
						if (pDiff >= 1d)
							numAbove1++;
						if (pDiff >= 0.1d)
							numAbove0p1++;
//					if (pDiff > printThreshold) {
//						System.out.println("Bad slip for event "+event.getID()+" (M="+(float)event.getMagnitude()+") at t="
//								+event.getTime()+" ("+(float)event.getTimeInYears()+" yrs)");
//						System.out.println("\tList File Slip:  "+(float)slips[i]);
//						System.out.println("\tTrans File Slip: "+(float)transSlip);
//						System.out.println("\t% Diff: "+(float)pDiff);
//						System.out.println("\tAbs Diff: "+(float)diff);
//					}
					}
				}
				double pDiff = DataUtils.getPercentDiff(transMoment, listMoment);
				momentStats.addValue(pDiff);
				double absDiff = Math.abs(transMoment - listMoment);
				momentAbsDiffStats.addValue(absDiff);
				
				double timeYears = event.getTimeInYears();
				if (transMoment > listMoment)
					momPDiffScatter.set(timeYears, pDiff);
				else
					momPDiffScatter.set(timeYears, -pDiff);
				double transMag = MagUtils.momentToMag(transMoment);
				double magDiff = transMag - event.getMagnitude();
				double magPDiff = DataUtils.getPercentDiff(transMag, event.getMagnitude());
				magDiffScatter.set(timeYears, magDiff);
				if (magDiff > 0)
					magPDiffScatter.set(timeYears, magPDiff);
				else
					magPDiffScatter.set(timeYears, -magPDiff);
				
				if (pDiff > printThreshold) {
					List<String> lines = new ArrayList<>();
					
					lines.add("Bad moment for event "+event.getID()+" (M="+(float)event.getMagnitude()+") at t="
							+event.getTime()+" ("+(float)event.getTimeInYears()+" yrs)");
					lines.add("\tList File Moment:  "+(float)listMoment+", M="+(float)MagUtils.momentToMag(listMoment));
					lines.add("\tTrans File Moment:  "+(float)transMoment+", M="+(float)MagUtils.momentToMag(transMoment));
					lines.add("\t% Diff: "+(float)pDiff);
					lines.add("\tAbs Diff: "+(float)absDiff+", M="+(float)MagUtils.momentToMag(absDiff));
					lines.add("");
					for (String line : lines) {
						System.out.println(line);
						debugFW.write(line+"\n");
					}
					debugFW.flush();
				}
				
				if (count % 10000 == 0) {
					System.out.println("=== SLIP STATS ===");
					System.out.println("Average % diff: "+(float)stats.getMean());
					System.out.println("Max % diff: "+(float)stats.getMax());
					System.out.println("Average diff: "+(float)absDiffStats.getMean());
					System.out.println("Max diff: "+(float)absDiffStats.getMax());
					System.out.println("=== MOMENT STATS ===");
					System.out.println("Average % diff: "+(float)momentStats.getMean());
					System.out.println("Max % diff: "+(float)momentStats.getMax());
					System.out.println("Average diff: "+(float)momentAbsDiffStats.getMean());
					System.out.println("Max diff: "+(float)momentAbsDiffStats.getMax()
						+", M="+((float)MagUtils.momentToMag(momentAbsDiffStats.getMax())));
					System.out.println("====================");
				}
			} catch (Exception e) {
				System.err.println("Errr with event!");
				errCount++;
				e.printStackTrace();
			}
		}
		if (errCount > 0) {
			System.out.println("WARNING!!!! "+errCount+" events errored!");
		}
		List<String> lines = new ArrayList<>();
		lines.add("");
		lines.add("*** SUMMARY INFO ****");
		double percentZeroTrans = 100d*(double)numZeroInTrans/(double)testCount;
		lines.add(numZeroInTrans+"/"+testCount+" ("+(float)percentZeroTrans+" %) zero in trans, nonzero in list");
		double percentZeroList = 100d*(double)numZeroInList/(double)testCount;
		lines.add(numZeroInList+"/"+testCount+" ("+(float)percentZeroList+" %) nonzero in trans, zero in list");
		lines.add("Average % diff: "+(float)stats.getMean());
		lines.add("Max % diff: "+(float)stats.getMax());
		lines.add("% diff std. dev.: "+(float)stats.getStandardDeviation());
		lines.add("Average diff: "+(float)absDiffStats.getMean());
		lines.add("Max diff: "+(float)absDiffStats.getMax());
		double percentAbove1 = 100d*(double)numAbove1/(double)testCount;
		lines.add(numAbove1+"/"+testCount+" ("+(float)percentAbove1+" %) above 1 %");
		double percentAbove0p1 = 100d*(double)numAbove0p1/(double)testCount;
		lines.add(numAbove0p1+"/"+testCount+" ("+(float)percentAbove0p1+" %) above 0.1 %");
		lines.add("");
		lines.add("Trans moment ave % diff: "+(float)momentStats.getMean());
		lines.add("Trans moment max % diff: "+(float)momentStats.getMax());
		lines.add("Trans moment ave diff: "+(float)momentAbsDiffStats.getMean());
		lines.add("Trans moment max diff: "+(float)momentAbsDiffStats.getMax());
		
		for (String line : lines) {
			System.out.println(line);
			debugFW.write(line+"\n");
		}
		
		debugFW.close();
		
		List<PlotSpec> specs = new ArrayList<>();
		
		String title = "Transition vs List File Descrapancies";
		specs.add(buildScatterSpec(momPDiffScatter, title, "Trans Moment vs List Moment (% diff)"));
		specs.add(buildScatterSpec(magPDiffScatter, title, "Trans Mag vs List Mag (% diff)"));
		specs.add(buildScatterSpec(magDiffScatter, title, "Trans Mag - List Mag"));
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		gp.setLegendFontSize(18);
		gp.setBackgroundColor(Color.WHITE);
		
		List<Range> xRanges = new ArrayList<>();
		xRanges.add(new Range(momPDiffScatter.getMinX(), momPDiffScatter.getMaxX()));
		
		gp.drawGraphPanel(specs, false, false, xRanges, null);
		gp.getChartPanel().setSize(5000, 1000);
		gp.saveAsPNG(new File(catalog.getCatalogDir(), "trans_scatter_comparisons.png").getAbsolutePath());
	}
	
	private static PlotSpec buildScatterSpec(XY_DataSet xy, String title, String yAxisLabel) {
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		funcs.add(xy);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, new Color(0, 0, 0, 80)));
		
		DefaultXY_DataSet line = new DefaultXY_DataSet();
		line.set(0d, 0d);
		line.set(xy.getMinX(), 0d);
		line.set(xy.getMaxX(), 0d);
		
		funcs.add(line);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		
		return new PlotSpec(funcs, chars, title, "Time (years)", yAxisLabel);
	}

}
