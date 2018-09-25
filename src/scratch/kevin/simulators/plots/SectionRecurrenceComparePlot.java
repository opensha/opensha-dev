package scratch.kevin.simulators.plots;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.IDPairing;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.utils.RSQSimUtils;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Lists;
import com.google.common.collect.Table;
import com.google.common.primitives.Doubles;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;

public class SectionRecurrenceComparePlot extends AbstractPlot {
	
	private Map<Integer, SimulatorElement> elemsMap;
	private String compName;
	private FaultSystemSolution compSol;
	private SectType sectType;
	private double minFractForInclusion;
	private double[] minMags;
	private double overallMinMag;
	
	public enum SectType {
		ELEMENT("Element", "Subsection", "elem"),
		SUBSECTION("Subsection", "Subsection", "subsect"),
		PARENT("Section", "Section", "parent");
		
		private String simType;
		private String compType;
		private String prefix;

		private SectType(String simType, String compType, String prefix) {
			this.simType = simType;
			this.compType = compType;
			this.prefix = prefix;
		}

		public String getSimType() {
			return simType;
		}

		public String getCompType() {
			return compType;
		}

		public String getPrefix() {
			return prefix;
		}
	}
	
	private Table<Integer, Double, Integer> sectCounts;
	
	private List<FaultSectionPrefData> compSects;
	private Map<IDPairing, Double> subSectDistsCache;
	private Map<Integer, Double> subSectAreas;
	private int minElemSectID;
	
	public SectionRecurrenceComparePlot(List<SimulatorElement> elems, FaultSystemSolution compSol, String compName,
			SectType sectType, double minFractForInclusion, double... minMags) {
		this.compSol = compSol;
		this.compName = compName;
		this.sectType = sectType;
		this.minFractForInclusion = minFractForInclusion;
		if (minMags == null || minMags.length == 0)
			minMags = new double[] { 0d };
		this.minMags = minMags;
		this.overallMinMag = StatUtils.min(minMags);

		elemsMap = new HashMap<>();
		for (SimulatorElement e : elems)
			elemsMap.put(e.getID(), e);
		
		sectCounts = HashBasedTable.create();
		subSectDistsCache = new HashMap<>();
		compSects = compSol.getRupSet().getFaultSectionDataList();
		subSectAreas = RSQSimUtils.calcSubSectAreas(elems, compSects);
		minElemSectID = RSQSimUtils.getSubSectIndexOffset(elems, compSects);
	}

	@Override
	protected void doProcessEvent(SimulatorEvent e) {
		double mag = e.getMagnitude();
		if (mag < overallMinMag)
			return;
		
		HashSet<Integer> ids = new HashSet<>();
		
		if (sectType == SectType.ELEMENT) {
			for (int id : e.getAllElementIDs())
				ids.add(id);
		} else {
			List<List<FaultSectionPrefData>> sects = RSQSimUtils.getSectionsForRupture(e, minElemSectID, compSects,
					subSectDistsCache, minFractForInclusion, subSectAreas);
			Preconditions.checkState(!sects.isEmpty());
			for (List<FaultSectionPrefData> sectsForParent : sects) {
				Preconditions.checkState(!sectsForParent.isEmpty());
				for (FaultSectionPrefData sect : sectsForParent) {
					if (sectType == SectType.PARENT)
						ids.add(sect.getParentSectionId());
					else
						ids.add(sect.getSectionId());
				}
			}
		}
		
		for (double minMag : minMags) {
			if (mag >= minMag) {
				for (Integer id : ids) {
					Integer prev = sectCounts.get(id, minMag);
					if (prev == null)
						prev = 0;
					sectCounts.put(id, minMag, prev + 1);
				}
			}
		}
	}

	@Override
	public void finalizePlot() throws IOException {
		double durationYears = getCurrentDurationYears();
		
		FaultSystemRupSet rupSet = compSol.getRupSet();
		
		File outputDir = getOutputDir();
		String prefix = getOutputPrefix();
		
		for (double minMag : minMags) {
			DefaultXY_DataSet scatter = new DefaultXY_DataSet();
			
			Map<Integer, Double> compVals = new HashMap<>();
			
			for (int id : sectCounts.rowKeySet()) {
				Integer count = sectCounts.get(id, minMag);
				if (count != null) {
					double simRI = durationYears/count.doubleValue();
					
					Integer compID;
					if (sectType == SectType.ELEMENT)
						compID = elemsMap.get(id).getSectionID()-minElemSectID;
					else
						compID = id;
					
					Double compRI = compVals.get(compID);
					if (compRI == null) {
						double compRate = 0d;
						List<Integer> compRups;
						if (sectType == SectType.PARENT)
							compRups = rupSet.getRupturesForParentSection(compID);
						else
							compRups = rupSet.getRupturesForSection(compID);
						for (int rup : compRups)
							if (rupSet.getMagForRup(rup) >= minMag)
								compRate += compSol.getRateForRup(rup);
						compRI = 1d/compRate;
						compVals.put(compID, compRI);
					}
					if (!Double.isFinite(compRI))
						continue;
					
//					System.out.println(u3RI+" "+simRI);
					scatter.set(compRI, simRI);
				}
			}
			
			if (scatter.size() == 0)
				continue;
			boolean hasFinite = false;
			for (int i=0; i<scatter.size(); i++) {
				if (Double.isFinite(scatter.getY(i))) {
					hasFinite = true;
					break;
				}
			}
			if (!hasFinite)
				continue;
			
			// Scatter plot
			System.out.println("Plotting scatter with "+scatter.size()+" values");
			String magLabel = getCleanMagLabel(minMag);
			String myPrefix = prefix+"_m"+magLabel;
			
			List<XY_DataSet> funcs = Lists.newArrayList();
			List<PlotCurveCharacterstics> chars = Lists.newArrayList();
			
			double minVal = minNonZero(scatter, true);
			minVal = Math.min(minVal, minNonZero(scatter, false));
			double maxVal = Math.max(scatter.getMaxX(), scatter.getMaxY());
			Range range = calcEncompassingLog10Range(minVal, maxVal);
			minVal = range.getLowerBound();
			maxVal = range.getUpperBound();
			
			funcs.add(scatter);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, Color.BLACK));
			
			funcs.add(getLine(null, minVal, minVal, maxVal, maxVal));
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GRAY));
			
			String title = "M≥"+magLabel+" "+sectType.simType+" Interevent Times";
			PlotSpec plot = new PlotSpec(funcs, chars, title, compName+" (years)", getCatalogName()+" (years)");
			plot.setLegendVisible(true);
			
			HeadlessGraphPanel gp = getGraphPanel();
			gp.drawGraphPanel(plot, true, true, range, range);
			gp.getChartPanel().setSize(getPlotWidth(), getPlotHeight());
//			System.out.println("About to plot "+minVal+" "+maxVal);
			gp.saveAsPNG(new File(getOutputDir(), myPrefix+"_scatter.png").getAbsolutePath());
//			System.out.println("Done PNG");
			gp.saveAsPDF(new File(getOutputDir(), myPrefix+"_scatter.pdf").getAbsolutePath());
//			System.out.println("Done PDF");
			
			// XYZ plot (2D hist)
			int numBins = 51;
			
			double logMin = Math.log10(range.getLowerBound());
			double logMax = Math.log10(range.getUpperBound());
			double gridSpacing = (logMax - logMin)/(numBins-1);
			
			EvenlyDiscrXYZ_DataSet xyz = new EvenlyDiscrXYZ_DataSet(numBins, numBins, logMin, logMin, gridSpacing);
			
			for (Point2D pt : scatter) {
				int index = xyz.indexOf(Math.log10(pt.getX()), Math.log10(pt.getY()));
				if (index < 0 || index >= xyz.size())
					throw new IllegalStateException("Scatter point not in XYZ range. x: "
								+pt.getX()+" ["+xyz.getMinX()+" "+xyz.getMaxX()
							+"], y: "+pt.getY()+" ["+xyz.getMinY()+" "+xyz.getMaxY()+"]");
				xyz.set(index, xyz.get(index)+1);
			}
//			// convert to density
//			for (int i=0; i<xyz.size(); i++) {
//				// convert to density
//				Point2D pt = xyz.getPoint(i);
//				double x = pt.getX();
//				double y = pt.getY();
//				double binWidth = Math.pow(10, x + 0.5*gridSpacing) - Math.pow(10, x - 0.5*gridSpacing);
//				double binHeight = Math.pow(10, y + 0.5*gridSpacing) - Math.pow(10, y - 0.5*gridSpacing);
//				double area = binWidth * binHeight;
//				xyz.set(i, xyz.get(i)*area);
//			}
//			xyz.scale(1d/xyz.getSumZ());
			
			// set all zero to NaN so that it will plot white
			for (int i=0; i<xyz.size(); i++) {
				if (xyz.get(i) == 0)
					xyz.set(i, Double.NaN);
			}
			xyz.log10();
			
			double minZ = Double.POSITIVE_INFINITY;
			double maxZ = Double.NEGATIVE_INFINITY;
			for (int i=0; i<xyz.size(); i++) {
				double val = xyz.get(i);
				if (!Doubles.isFinite(val))
					continue;
				if (val < minZ)
					minZ = val;
				if (val > maxZ)
					maxZ = val;
			}
			
//			System.out.println("MinZ: "+minZ);
//			System.out.println("MaxZ: "+maxZ);
			
			CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, maxZ);
//			if ((float)minZ == (float)maxZ)
//				cpt = cpt.rescale(minZ, minZ*2);
//			else if (!Doubles.isFinite(minZ))
//				cpt = cpt.rescale(0d, 1d);
//			else
//				cpt = cpt.rescale(minZ, maxZ);
			cpt.setNanColor(Color.WHITE);
			
			String zAxisLabel = "Log10 Num "+sectType.simType+"s";
			XYZPlotSpec xyzSpec = new XYZPlotSpec(xyz, cpt, title, "Log10 "+compName+" (years)",
					"Log10 "+getCatalogName()+" (years)", zAxisLabel);
			funcs = Lists.newArrayList();
			chars = Lists.newArrayList();
			funcs.add(getLine(null, logMin, logMin, logMax, logMax));
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GRAY));
			xyzSpec.setXYElems(funcs);
			xyzSpec.setXYChars(chars);
			
			XYZGraphPanel xyzGP = getXYZGraphPanel();
			Range logRange = new Range(logMin, logMax);
			xyzGP.drawPlot(xyzSpec, false, false, logRange, logRange);
			// write plot
			xyzGP.getChartPanel().setSize(getPlotWidth(), getPlotHeight());
			xyzGP.saveAsPNG(new File(outputDir, myPrefix+"_hist2D.png").getAbsolutePath());
			xyzGP.saveAsPDF(new File(outputDir, myPrefix+"_hist2D.pdf").getAbsolutePath());
		}
	}

	@Override
	public Collection<SimulatorElement> getApplicableElements() {
		return null;
	}

}
