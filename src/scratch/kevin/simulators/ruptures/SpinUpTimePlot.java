package scratch.kevin.simulators.ruptures;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class SpinUpTimePlot {

	public static void main(String[] args) throws IOException {
		double maxTime = 150000;
		RSQSimCatalog catalog = Catalogs.BRUCE_4983_STITCHED.instance();
		
		double[] minMags = { 0d };
		double[] printTimes = { 1000d, 5000d, 10000d, 20000d, 50000d,
				100000d, 150000d };
		
		EvenlyDiscretizedFunc timeDiscr = new EvenlyDiscretizedFunc(0d, maxTime, 1000);
		
		int[][][] elemsRuptured = new int[minMags.length][catalog.getElements().size()][timeDiscr.size()];
		Table<Integer, Double, int[]> subSectsRuptured = HashBasedTable.create();
		RSQSimSubSectionMapper mapper = catalog.getSubSectMapper();
		Table<Integer, Double, int[]> sectsRuptured = HashBasedTable.create();
		for (SimulatorElement elem : catalog.getElements()) {
			FaultSection sect = mapper.getMappedSection(elem);
			if (!subSectsRuptured.containsRow(sect.getSectionId()))
				for (double minMag : minMags)
					subSectsRuptured.put(sect.getSectionId(), minMag, new int[timeDiscr.size()]);
			if (!sectsRuptured.containsRow(sect.getParentSectionId()))
				for (double minMag : minMags)
					sectsRuptured.put(sect.getParentSectionId(), minMag, new int[timeDiscr.size()]);
		}
		int minElemID = catalog.getElements().get(0).getID();
		System.out.println(subSectsRuptured.size()+" sub sects");
		System.out.println(sectsRuptured.size()+" sects");
		double firstTime = Double.NaN;
		HashSet<Integer> sectsRupturedSet = new HashSet<>();
		for (RSQSimEvent event : catalog.loader().maxDuration(maxTime).iterable()) {
			double time = event.getTimeInYears();
			if (Double.isNaN(firstTime))
				firstTime = time;
			double timeDelta = time - firstTime;
			int timeIndex = timeDiscr.getClosestXIndex(timeDelta);
			double mag = event.getMagnitude();
			for (SimulatorElement elem : event.getAllElements()) {
				FaultSection sect = mapper.getMappedSection(elem);
				int elemIndex = elem.getID()-minElemID;
				for (int m=0; m<minMags.length; m++) {
					if (minMags[m] > 0 && mag < minMags[m])
						break;
					elemsRuptured[m][elemIndex][timeIndex]++;
					subSectsRuptured.get(sect.getSectionId(), minMags[m])[timeIndex]++;
					sectsRuptured.get(sect.getParentSectionId(), minMags[m])[timeIndex]++;
					if (!sectsRupturedSet.contains(sect.getParentSectionId())) {
						sectsRupturedSet.add(sect.getParentSectionId());
						int numLeft = sectsRuptured.rowKeySet().size() - sectsRupturedSet.size();
						System.out.println("A fault ("+numLeft+" left) had it's "
								+ "first rupture at "+(float)timeDelta+": "+sect.getParentSectionName());
					}
				}
			}
			
			
		}
		
		System.out.println("Converting to cumulative...");
		
		// convert to cumulative
		for (int m=0; m<minMags.length; m++)
			for (int i=0; i<elemsRuptured[m].length; i++)
				toCumulative(elemsRuptured[m][i]);
		for (int[] timeCounts : subSectsRuptured.values())
			toCumulative(timeCounts);
		for (int[] timeCounts : sectsRuptured.values())
			toCumulative(timeCounts);
		
		EvenlyDiscretizedFunc[] fractElemsRuptured = new EvenlyDiscretizedFunc[minMags.length];
		EvenlyDiscretizedFunc[] fractSubSectsRuptured = new EvenlyDiscretizedFunc[minMags.length];
		EvenlyDiscretizedFunc[] fractSectsRuptured = new EvenlyDiscretizedFunc[minMags.length];
		
		List<PlotSpec> specs = new ArrayList<>();
		Range xRange = new Range(0d, timeDiscr.getMaxX());
		Range yRange = new Range(0d, 1d);
		List<Range> xRanges = new ArrayList<>();
		xRanges.add(xRange);
		List<Range> yRanges = new ArrayList<>();
		
		for (int m=0; m<minMags.length; m++) {
			System.out.println("Plotting M>="+minMags[m]);
			
			fractElemsRuptured[m] = timeDiscr.deepClone();
			fractSubSectsRuptured[m] = timeDiscr.deepClone();
			fractSectsRuptured[m] = timeDiscr.deepClone();
			
			for (int i=0; i<timeDiscr.size(); i++) {
				int elemCount = 0;
				for (int e=0; e<elemsRuptured[m].length; e++)
					if (elemsRuptured[m][e][i] > 0)
						elemCount++;
				fractElemsRuptured[m].set(i, (double)elemCount/(double)elemsRuptured[m].length);
				fractSubSectsRuptured[m].set(i, count(subSectsRuptured, minMags[m], i));
				fractSectsRuptured[m].set(i, count(sectsRuptured, minMags[m], i));
			}
			
			List<DiscretizedFunc> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			fractSectsRuptured[m].setName("Sections");
			funcs.add(fractSectsRuptured[m]);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
			
			fractSubSectsRuptured[m].setName("Subsections");
			funcs.add(fractSubSectsRuptured[m]);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
			
			fractElemsRuptured[m].setName("Elements");
			funcs.add(fractElemsRuptured[m]);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
			
			for (double time : printTimes) {
				System.out.println("\t"+(int)time);
				System.out.println("\t\t"+(float)fractElemsRuptured[m].getInterpolatedY(time)+" elems");
				System.out.println("\t\t"+(float)fractSubSectsRuptured[m].getInterpolatedY(time)+" subsects");
				System.out.println("\t\t"+(float)fractSectsRuptured[m].getInterpolatedY(time)+" sects");
			}
			
			PlotSpec spec = new PlotSpec(funcs, chars, "Spin-Up Time Comparison", "Time (years)",
					"Fract Ruptured");
			spec.setLegendVisible(specs.isEmpty());
			
			String text = minMags[m] > 0 ? "M≥"+(float)minMags[m] : "Any Mag";
			XYTextAnnotation magAnn = new XYTextAnnotation(text, 0.975*timeDiscr.getMaxX(), 0.025);
			magAnn.setTextAnchor(TextAnchor.BOTTOM_RIGHT);
			magAnn.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 24));
			spec.addPlotAnnotation(magAnn);
			
			specs.add(spec);
			yRanges.add(yRange);
		}
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(28);
		gp.setBackgroundColor(Color.WHITE);
		
		gp.drawGraphPanel(specs, false, false, xRanges, yRanges);
		
		File file = new File(catalog.getCatalogDir(), "spin_up_time");
		gp.getChartPanel().setSize(1000, 1200);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
	}
	
	private static void toCumulative(int[] timeCounts) {
		int curTot = 0;
		for (int i=0; i<timeCounts.length; i++) {
			curTot += timeCounts[i];
			timeCounts[i] = curTot;
		}
	}
	
	private static double count(Table<Integer, Double, int[]> table, double minMag, int timeIndex) {
		int rows = table.rowKeySet().size();
		int rowsWith = 0;
		for (int id : table.rowKeySet()) {
			int[] counts = table.get(id, minMag);
			if (counts[timeIndex] > 0)
				rowsWith++;
		}
		return (double)rowsWith/(double)rows;
	}

}
