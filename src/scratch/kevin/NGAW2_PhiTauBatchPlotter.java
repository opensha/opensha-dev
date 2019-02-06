package scratch.kevin;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.sha.imr.attenRelImpl.ngaw2.ASK_2014;
import org.opensha.sha.imr.attenRelImpl.ngaw2.BSSA_2014;
import org.opensha.sha.imr.attenRelImpl.ngaw2.CB_2014;
import org.opensha.sha.imr.attenRelImpl.ngaw2.CY_2014;
import org.opensha.sha.imr.attenRelImpl.ngaw2.FaultStyle;
import org.opensha.sha.imr.attenRelImpl.ngaw2.IMT;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_GMM;
import org.opensha.sha.imr.attenRelImpl.ngaw2.ScalarGroundMotion;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.collect.Table.Cell;

public class NGAW2_PhiTauBatchPlotter {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/home/kevin/git/misc-research/nga_w2_phi_tau");
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		NGAW2_GMM[] gmpes = { new ASK_2014(), new BSSA_2014(), new CB_2014(), new CY_2014() };
		
		String[] names = {"ASK", "BSSA", "CB", "CY"};
		String phi_html = "&phi;";
		String tau_html = "&tau;";
		
		double minDist = 10d;
		double maxDist = 200d;
		double deltaDist = 10d;
		double[] fixedDists = {20, 50, 100, 200};
		
		double minMag = 5d;
		double maxMag = 8d;
		double deltaMag = 0.1;
		double[] fixedMags = { 6, 6.5, 7, 7.5 };
		
		double[] periods = { 1d, 2d, 3d, 4, 5d, 7.5d, 10d };
		double[] highlightPeriods = { 3d, 5d, 7.5d, 10d };
		IMT[] imts = new IMT[periods.length];
		for (int p=0; p<periods.length; p++)
			imts[p]= IMT.getSA(periods[p]);
		
		EvenlyDiscretizedFunc distFunc = new EvenlyDiscretizedFunc(minDist, (int)((maxDist-minDist)/deltaDist)+1, deltaDist);
		EvenlyDiscretizedFunc magFunc = new EvenlyDiscretizedFunc(minMag, (int)((maxMag-minMag)/deltaMag)+1, deltaMag);
		
		Table<Double, Double, Table<NGAW2_GMM, Double, ScalarGroundMotion>> resultsTable = HashBasedTable.create();
		
		List<String> lines = new ArrayList<>();
		lines.add("# NGA-West2 Phi/Tau Dependence");
		lines.add("");
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		lines.add("## Constant parameters");
		lines.add(topLink); lines.add("");
		lines.add("All calculations hold the following parameters constant:");
		lines.add("");
		lines.add("* Dip: 90");
		lines.add("* Fault Type: Strike-Slip");
		lines.add("* Vs30: 863");
		lines.add("* Vs30 Inferred: false");
		lines.add("* Width: 10");
		lines.add("* zHyp: 5");
		lines.add("* zTop: 0");
		lines.add("* rJB = rRup = rX");
		
		for (NGAW2_GMM gmm : gmpes) {
			gmm.set_dip(90);
			gmm.set_fault(FaultStyle.STRIKE_SLIP);
			gmm.set_vs30(863);
			gmm.set_vsInf(false);
			gmm.set_width(10);
			gmm.set_zHyp(5d);
			gmm.set_zTop(0);
		}
		
		for (int i=0; i<distFunc.size(); i++) {
			double distance = distFunc.getX(i);
			for (int j=0; j<magFunc.size(); j++) {
				double mag = magFunc.getX(j);
				
				Table<NGAW2_GMM, Double, ScalarGroundMotion> table = HashBasedTable.create();
				resultsTable.put(distance, mag, table);
				
				for (NGAW2_GMM gmm : gmpes) {
					gmm.set_Mw(mag);
					gmm.set_rJB(distance);
					gmm.set_rX(distance);
					gmm.set_rRup(distance);
					for (int p=0; p<periods.length; p++) {
						gmm.set_IMT(imts[p]);
						table.put(gmm, periods[p], gmm.calc());
					}
				}
			}
		}
		
		lines.add("## Summary Table");
		lines.add(topLink); lines.add("");
		List<String> tableHeader = new ArrayList<>();
		tableHeader.add("**Distance (km)**");
		tableHeader.add("**Mw**");
		tableHeader.add("**Period (s)**");
		for (String name : names) {
			tableHeader.add("**"+name+" "+phi_html+"**");
			tableHeader.add("**"+name+" "+tau_html+"**");
		}
		double minY = Double.POSITIVE_INFINITY;
		double maxY = Double.NEGATIVE_INFINITY;
		TableBuilder table = MarkdownUtils.tableBuilder();
		for (double distance : fixedDists) {
			table.addLine(tableHeader);
			for (double mag : fixedMags) {
				Table<NGAW2_GMM, Double, ScalarGroundMotion> results = resultsTable.get(distance, mag);
				for (double period : highlightPeriods) {
					table.initNewLine();
					table.addColumn(optionalDigitDF.format(distance));
					table.addColumn(optionalDigitDF.format(mag));
					table.addColumn(optionalDigitDF.format(period));
					for (NGAW2_GMM gmm : gmpes) {
						ScalarGroundMotion gm = results.get(gmm, period);
						minY = Double.min(minY, gm.phi());
						minY = Double.min(minY, gm.tau());
						maxY = Double.max(maxY, gm.phi());
						maxY = Double.max(maxY, gm.tau());
						table.addColumn(optionalDigitDF.format(gm.phi()));
						table.addColumn(optionalDigitDF.format(gm.tau()));
					}
					table.finalizeLine();
				}
			}
		}
		lines.addAll(table.build());
		lines.add("");
		
		Range yRange = new Range(Math.floor(minY*10d)/10d, Math.ceil(maxY*10d)/10d);
		
		lines.add("## Distance Dependence");
		lines.add(topLink); lines.add("");
		
		table = MarkdownUtils.tableBuilder();
		table.initNewLine();
		table.addColumn("Mw, Period");
		for (double period : highlightPeriods)
			table.addColumn(optionalDigitDF.format(period)+" s");
		table.finalizeLine();
		
		for (double mag : fixedMags) {
			table.initNewLine();
			table.addColumn("**"+optionalDigitDF.format(mag)+"**");
			for (double period : highlightPeriods) {
				File plotFile = plotForFixed(resourcesDir, resultsTable, gmpes, names, yRange, null, mag, period);
				table.addColumn("![Plot](resources/"+plotFile.getName()+")");
			}
			table.finalizeLine();
		}
		lines.addAll(table.build());
		lines.add("");
		
		lines.add("## Magnitude Dependence");
		lines.add(topLink); lines.add("");
		
		table = MarkdownUtils.tableBuilder();
		table.initNewLine();
		table.addColumn("Distance, Period");
		for (double period : highlightPeriods)
			table.addColumn(optionalDigitDF.format(period)+" s");
		table.finalizeLine();
		
		for (double distance : fixedDists) {
			table.initNewLine();
			table.addColumn("**"+optionalDigitDF.format(distance)+" km**");
			for (double period : highlightPeriods) {
				File plotFile = plotForFixed(resourcesDir, resultsTable, gmpes, names, yRange, distance, null, period);
				table.addColumn("![Plot](resources/"+plotFile.getName()+")");
			}
			table.finalizeLine();
		}
		lines.addAll(table.build());
		lines.add("");
		
		lines.add("## Period Dependence");
		lines.add(topLink); lines.add("");
		
		table = MarkdownUtils.tableBuilder();
		table.initNewLine();
		table.addColumn("Mw, Distance");
		for (double distance : fixedDists)
			table.addColumn(optionalDigitDF.format(distance)+" km");
		table.finalizeLine();
		
		for (double mag : fixedMags) {
			table.initNewLine();
			table.addColumn("**"+optionalDigitDF.format(mag)+"**");
			for (double distance : fixedDists) {
				File plotFile = plotForFixed(resourcesDir, resultsTable, gmpes, names, yRange, distance, mag, null);
				table.addColumn("![Plot](resources/"+plotFile.getName()+")");
			}
			table.finalizeLine();
		}
		lines.addAll(table.build());
		lines.add("");
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2, 4));
		lines.add(tocIndex, "## Table Of Contents");
		
		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	private static Color[] colors = { Color.BLACK, Color.BLUE, Color.RED, Color.GREEN };
	private static PlotLineType phi_line_type = PlotLineType.SOLID;
	private static PlotLineType tau_line_type = PlotLineType.DASHED;
	
	private static File plotForFixed(File resourcesDir, Table<Double, Double, Table<NGAW2_GMM, Double, ScalarGroundMotion>> resultsTable,
			NGAW2_GMM[] gmms, String[] names, Range yRange, Double fixedDist, Double fixedMag, Double fixedPeriod) throws IOException {
		String xAxisLabel = null;
		String prefix;
		if (fixedDist == null) {
			xAxisLabel = "Distance (km)";
			prefix = "distVar";
		} else {
			prefix = "dist"+optionalDigitDF.format(fixedDist);
		}
		
		if (fixedMag == null) {
			Preconditions.checkState(xAxisLabel == null, "Multiple variable quantities!");
			xAxisLabel = "Magnitude";
			prefix += "_magVar";
		} else {
			prefix += "_mag"+optionalDigitDF.format(fixedMag);
		}
		
		if (fixedPeriod == null) {
			Preconditions.checkState(xAxisLabel == null, "Multiple variable quantities!");
			xAxisLabel = "Period (s)";
			prefix += "_periodVar";
		} else {
			prefix += "_period"+optionalDigitDF.format(fixedPeriod);
		}
		
		Preconditions.checkNotNull(xAxisLabel, "No variable quantities!");
		
		DiscretizedFunc[] phiFuncs = new DiscretizedFunc[gmms.length];
		DiscretizedFunc[] tauFuncs = new DiscretizedFunc[gmms.length];
		for (int i=0; i<gmms.length; i++) {
			phiFuncs[i] = new ArbitrarilyDiscretizedFunc(names[i]+" ϕ");
			tauFuncs[i] = new ArbitrarilyDiscretizedFunc(names[i]+" τ");
		}
		
		for (Cell<Double, Double, Table<NGAW2_GMM, Double, ScalarGroundMotion>> cell : resultsTable.cellSet()) {
			Double x = null;
			
			double distance = cell.getRowKey();
			if (fixedDist == null)
				x = distance;
			else if ((float)distance != fixedDist.floatValue())
				continue;
			
			double mag = cell.getColumnKey();
			if (fixedMag == null)
				x = mag;
			else if ((float)mag != fixedMag.floatValue())
				continue;
			
			Table<NGAW2_GMM, Double, ScalarGroundMotion> results = resultsTable.get(distance, mag);
			for (Double period : results.columnKeySet()) {
				if (fixedPeriod == null)
					x = period;
				else if (period.floatValue() != fixedPeriod.floatValue())
					continue;
				for (int i=0; i<gmms.length; i++) {
					ScalarGroundMotion gm = results.get(gmms[i], period);
					phiFuncs[i].set(x, gm.phi());
					tauFuncs[i].set(x, gm.tau());
				}
			}
		}
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		for (int i=0; i<gmms.length; i++) {
			funcs.add(phiFuncs[i]);
			chars.add(new PlotCurveCharacterstics(phi_line_type, 3f, colors[i % colors.length]));
			funcs.add(tauFuncs[i]);
			chars.add(new PlotCurveCharacterstics(tau_line_type, 3f, colors[i % colors.length]));
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, null, xAxisLabel, "Standard Deviation");
		spec.setLegendVisible(true);

		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(18);
		plotPrefs.setAxisLabelFontSize(20);
		plotPrefs.setPlotLabelFontSize(21);
		plotPrefs.setLegendFontSize(24);
		plotPrefs.setBackgroundColor(Color.WHITE);

		HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
		
		Range xRange = null;

		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		gp.getChartPanel().setSize(800, 600);
		File pngFile = new File(resourcesDir, prefix+".png");
		gp.saveAsPNG(pngFile.getAbsolutePath());
		return pngFile;
	}
	
	static final DecimalFormat optionalDigitDF = new DecimalFormat("0.##");

}
