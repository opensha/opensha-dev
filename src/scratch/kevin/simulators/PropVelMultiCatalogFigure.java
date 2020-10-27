package scratch.kevin.simulators;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYBoxAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.data.Range;
import org.jfree.chart.ui.TextAnchor;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.IDPairing;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class PropVelMultiCatalogFigure {

	public static void main(String[] args) throws IOException {
		boolean logX = true;
		boolean median = false;
		
		EvenlyDiscretizedFunc xVals;
		double maxDist = 200d;
		if (logX)
			xVals = new EvenlyDiscretizedFunc(0d, (int)(Math.log10(maxDist)/0.1+1), 0.1);
		else
			xVals = new EvenlyDiscretizedFunc(5d, 20, 10d);
		
		File outputDir = new File("/home/kevin/Documents/papers/2020_RSQSim_PSHA");
		String prefix = "figure_5";
		
		CSVFile<String> csv = new CSVFile<>(true);
		
		double xStart = xVals.getMinX()-0.5*xVals.getDelta();
		double xEnd = xVals.getMaxX()+0.5*xVals.getDelta();
		
		RSQSimCatalog[] catalogs = {
			Catalogs.BRUCE_2585.instance(),
			Catalogs.BRUCE_2740.instance(),
			Catalogs.BRUCE_4317.instance(),
//			Catalogs.BRUCE_4860_10X.instance()
			Catalogs.BRUCE_4983_STITCHED.instance()
		};
		
		String[] names = {
				"Shaw et al. (2018)",
				"+ Finite Receiver",
				"+ Variable Slip Speed",
				"+ Time-Delay"
		};
		
		PlotCurveCharacterstics[] pChars = {
				new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.BLACK),
				new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, Color.BLACK),
				new PlotCurveCharacterstics(PlotLineType.DOTTED_AND_DASHED, 3f, Color.BLACK),
				new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK)
		};
		
		double minMag = 6.5;
		double skipYears = 5000;
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		Range xRange = logX ? new Range(Math.pow(10, xStart), Math.pow(10, xEnd)) : new Range(xStart, xEnd);
		Range yRange = new Range(0d, 6d);
		
		List<XYAnnotation> anns = new ArrayList<>();
		Font annFont = new Font(Font.SANS_SERIF, Font.BOLD, 18);
		
		boolean plotCatVal = false;
		double refAnnCenterFract = plotCatVal ? 0.77 : 0.92;
		
		double annLeftFract = refAnnCenterFract-0.035;
		double annRightFract = refAnnCenterFract+0.035;
		
		double annRightJustify, annLeftJustify;
		if (logX) {
			annLeftJustify = Math.pow(10, xStart + (xEnd-xStart)*annLeftFract);
			annRightJustify = Math.pow(10, xStart + (xEnd-xStart)*annRightFract);
		} else {
			annLeftJustify = xRange.getLowerBound() + xRange.getLength()*annLeftFract;
			annRightJustify = xRange.getLowerBound() + xRange.getLength()*annRightFract;
		}
		double annDeltaY = 0.3;
		double annY = yRange.getUpperBound()-annDeltaY;
		
		for (int c=0; c<catalogs.length; c++) {
			RSQSimCatalog catalog = catalogs[c];
			List<List<Double>> binVals = new ArrayList<>();
			for (int i=0; i<xVals.size(); i++)
				binVals.add(new ArrayList<>());
			
			Map<IDPairing, Double> distsMap = new HashMap<>();
			
			System.out.println("Processing "+catalog.getName());
			
			List<Double> rupVals = new ArrayList<>();
			
			for (RSQSimEvent e : catalog.loader().minMag(minMag).skipYears(skipYears).iterable()) {
				List<SimulatorElement> elems = e.getAllElements();
				double[] timeFirsts = e.getAllElementTimes();
				
				SimulatorElement hypoEl = null;
				Location hypoLoc = null;
				double tFirst = Double.POSITIVE_INFINITY;
				
				for (int i=0; i<timeFirsts.length; i++) {
					if (timeFirsts[i] < tFirst) {
						tFirst = timeFirsts[i];
						hypoEl = elems.get(i);
						hypoLoc = hypoEl.getCenterLocation();
					}
				}
				
				int hypoID = hypoEl.getID();
				
				List<Double> myVals = new ArrayList<>();
				
				for (int i=0; i<timeFirsts.length; i++) {
					double time = timeFirsts[i] - tFirst;
					SimulatorElement elem = elems.get(i);
					
					int id = elem.getID();
					
					if (id == hypoID || time == 0d)
						continue;
					
					IDPairing pair = id < hypoID ? new IDPairing(id, hypoID) : new IDPairing(hypoID, id);
					
					Double dist = distsMap.get(pair);
					if (dist == null) {
						dist = LocationUtils.linearDistanceFast(hypoLoc, elem.getCenterLocation());
						Preconditions.checkState(dist > 0, "Bad dist between %s and %s: ", id, hypoID, dist);
						distsMap.put(pair, dist);
					}
					
					
					double vel = dist/time;
					if (logX)
						dist = Math.log10(dist);
					if (dist < xStart || dist > xEnd)
						continue;
					int index = xVals.getClosestXIndex(dist);
					binVals.get(index).add(vel);
					myVals.add(vel);
				}
				
				double[] array = Doubles.toArray(myVals);
				double myVal = median ? DataUtils.median(array) : StatUtils.mean(array);
				rupVals.add(myVal);
			}
			
			ArbitrarilyDiscretizedFunc func = new ArbitrarilyDiscretizedFunc(names[c]);
			
			List<String> header = null;
			if (csv.getNumRows() == 0) {
				header = new ArrayList<>();
				header.add("Model");
			}
			List<String> line = new ArrayList<>();
			line.add(names[c]);
			
			for (int i=0; i<xVals.size(); i++) {
				double x = xVals.getX(i);
				if (logX)
					x = Math.pow(10, x);
				
				if (header != null)
					header.add(new DecimalFormat("0.#").format(x));
				
				List<Double> vals = binVals.get(i);
				if (vals.isEmpty())
					continue;
				double[] array = Doubles.toArray(vals);
				double y = median ? DataUtils.median(array) : StatUtils.mean(array);
				Preconditions.checkState(Double.isFinite(y));
				
				line.add(new DecimalFormat("0.00").format(y));
				
				func.set(x, y);
			}
			
			if (header != null)
				csv.addLine(header);
			csv.addLine(line);
			
			double[] array = Doubles.toArray(rupVals);
			double catVal = median ? DataUtils.median(array) : StatUtils.mean(array);
			
			XYTextAnnotation nameAnn = new XYTextAnnotation(names[c]+"  ", annLeftJustify, annY);
			nameAnn.setFont(annFont);
			nameAnn.setTextAnchor(TextAnchor.CENTER_RIGHT);
			anns.add(nameAnn);
			if (plotCatVal) {
				String valStr = median ? "mdn=" : "mean=";
				valStr += new DecimalFormat("0.0").format(catVal);
				XYTextAnnotation valAnn = new XYTextAnnotation("  "+valStr, annRightJustify, annY);
				valAnn.setFont(annFont);
				valAnn.setTextAnchor(TextAnchor.CENTER_LEFT);
				anns.add(valAnn);
			}
			
			DiscretizedFunc legend = new ArbitrarilyDiscretizedFunc();
			legend.set(annLeftJustify, annY);
			legend.set(annRightJustify, annY);
			
			annY -= annDeltaY;
			
			funcs.add(legend);
			chars.add(pChars[c]);
			
			funcs.add(func);
			chars.add(pChars[c]);
		}
		
		String yAxisLabel = median ? "Median" : "Mean";
		yAxisLabel += " Propagation Velocity (km/s)";
		
		PlotSpec spec = new PlotSpec(funcs, chars, " ", "Hypocentral Distance (km)", yAxisLabel);
//		spec.setLegendVisible(true);
		spec.setPlotAnnotations(anns);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(22);
		gp.setBackgroundColor(Color.WHITE);
		gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
		
		gp.drawGraphPanel(spec, logX, false, xRange, yRange);
		
		File file = new File(outputDir, prefix);
		gp.getChartPanel().setSize(800, 600);
		gp.saveAsPNG(file.getAbsolutePath()+".png");
		gp.saveAsPDF(file.getAbsolutePath()+".pdf");
		
		csv.writeToFile(new File(outputDir, prefix+".csv"));
	}

}
