package scratch.kevin.simulators.plots;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SubSectionMapping;

import scratch.UCERF3.utils.paleoRateConstraints.PaleoProbabilityModel;
import scratch.UCERF3.utils.paleoRateConstraints.PaleoRateConstraint;
import scratch.UCERF3.utils.paleoRateConstraints.UCERF3_PaleoProbabilityModel;
import scratch.UCERF3.utils.paleoRateConstraints.UCERF3_PaleoRateConstraintFetcher;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class PaleoRecurrencePlot extends AbstractPlot {
	
	private static boolean D = false;
	
	private List<PaleoResult> results;
	private Map<SimulatorElement, List<PaleoResult>> elemPaleoMappings;
	private Map<FaultSection, List<PaleoResult>> sectPaleoMappings;
	
	private PaleoProbabilityModel paleoProbModel;
	
	private static DiscretizedFunc slipProbFunc = PaleoOpenIntervalPlot.slipProbFunc;

	private RSQSimSubSectionMapper mapper;
	private HashSet<SimulatorElement> applicableElements;
	
	public PaleoRecurrencePlot(List<SimulatorElement> elements, RSQSimSubSectionMapper mapper) throws IOException {
		this.mapper = mapper;
		List<? extends FaultSection> subSects = mapper.getSubSections();
		ArrayList<PaleoRateConstraint> constraints = UCERF3_PaleoRateConstraintFetcher.getConstraints(subSects);
		
		paleoProbModel = UCERF3_PaleoProbabilityModel.load();
		
		sectPaleoMappings = new HashMap<>();
		elemPaleoMappings = new HashMap<>();
		results = new ArrayList<>();
		
		applicableElements = new HashSet<>();
		
		for (PaleoRateConstraint constr : constraints) {
			if (D) System.out.println("Paleo constraint: "+constr.getPaleoSiteName());
			int sectIndex = constr.getSectionIndex();
			FaultSection subSect = subSects.get(sectIndex);
			if (D) System.out.println("\tMapped section: "+subSect.getSectionName());
			if (!mapper.isMapped(subSect)) {
				if (D) System.out.println("\tSkipping section (no mapped elements)");
				continue;
			}
			applicableElements.addAll(mapper.getElementsForSection(subSect));
			PaleoResult result = new PaleoResult(constr);
			results.add(result);
			List<PaleoResult> sectResults = sectPaleoMappings.get(subSect);
			if (sectResults == null) {
				sectResults = new ArrayList<>();
				sectPaleoMappings.put(subSect, sectResults);
			}
			sectResults.add(result);
			// look for surface rupturing element
			SimulatorElement closestElem = null;
			double closestDist = Double.POSITIVE_INFINITY;
			for (SimulatorElement elem : mapper.getElementsForSection(subSect)) {
				int numSurface = 0;
				for (Location loc : elem.getVertices())
					if (loc.getDepth() == 0d)
						numSurface++;
				if (numSurface >= 2) {
					double dist = LocationUtils.horzDistanceFast(elem.getCenterLocation(), constr.getPaleoSiteLoction());
					if (dist < closestDist) {
						closestDist = dist;
						closestElem = elem;
					}
				}
			}
			if (closestElem == null) {
				if (D) System.out.println("\tNo surface rupturing elements");
			} else {
				if (D) System.out.println("\tClosest surface rupturing element is "+(float)closestDist+" km away");
				if (closestDist < 5d) {
					List<PaleoResult> elemResults = elemPaleoMappings.get(closestElem);
					if (elemResults == null) {
						elemResults = new ArrayList<>();
						elemPaleoMappings.put(closestElem, elemResults);
					}
					elemResults.add(result);
				}
			}
		}
	}
	
	private class PaleoResult {
		private PaleoRateConstraint constraint;
		
		private int sectEventCount = 0;
		private int elemEventCount = 0;
		private double sectPaleoProbCount = 0d;
		private double elemPaleoProbCount = 0d;
		
		public PaleoResult(PaleoRateConstraint constraint) {
			this.constraint = constraint;
		}
		
		public void addSubSectMatch(double mag, List<FaultSection> rupSections, int sectionIndex) {
			sectEventCount++;
			sectPaleoProbCount += paleoProbModel.getProbPaleoVisible(mag, rupSections, sectionIndex);
		}
		
		public void addElementMatch(double mag, double slip) {
			elemEventCount++;
			elemPaleoProbCount += slip >= slipProbFunc.getMaxX() ? slipProbFunc.getMaxY() : slipProbFunc.getInterpolatedY(slip);
		}
	}

	@Override
	protected void doProcessEvent(SimulatorEvent e) {
		RSQSimEvent event = (RSQSimEvent)e;
		List<List<SubSectionMapping>> mappedSects = mapper.getFilteredSubSectionMappings(event);
		List<FaultSection> sects = new ArrayList<>();
		for (List<SubSectionMapping> bundle : mappedSects)
			for (SubSectionMapping mapping : bundle)
				sects.add(mapping.getSubSect());
		
		for (FaultSection sect : sects)
			if (sectPaleoMappings.containsKey(sect))
				for (PaleoResult result : sectPaleoMappings.get(sect))
					result.addSubSectMatch(event.getMagnitude(), sects, sect.getSectionId());
		
		ArrayList<SimulatorElement> allElements = event.getAllElements();
		double[] slips = event.getAllElementSlips();
		for (int i=0; i<allElements.size(); i++) {
			SimulatorElement elem = allElements.get(i);
			if (elemPaleoMappings.containsKey(elem))
				for (PaleoResult result : elemPaleoMappings.get(elem))
					result.addElementMatch(event.getMagnitude(), slips[i]);
		}
	}

	@Override
	public void finalizePlot() throws IOException {
		double totalDuration = getCurrentDurationYears();
		
		Map<PaleoRateConstraint, Double> rawSectRates = new HashMap<>();
		Map<PaleoRateConstraint, Double> paleoProbSectRates = new HashMap<>();
		Map<PaleoRateConstraint, Double> rawElemRates = new HashMap<>();
		Map<PaleoRateConstraint, Double> paleoProbElemRates = new HashMap<>();
		
		for (PaleoResult result : results) {
			rawSectRates.put(result.constraint, result.sectEventCount/totalDuration);
			paleoProbSectRates.put(result.constraint, result.sectPaleoProbCount/totalDuration);
			rawElemRates.put(result.constraint, result.elemEventCount/totalDuration);
			paleoProbElemRates.put(result.constraint, result.elemPaleoProbCount/totalDuration);
		}
		
		File outputDir = getOutputDir();
		String prefix = getOutputPrefix();
		plot(outputDir, prefix+"_raw_sect_rate", rawSectRates, "Raw Simulated Subsection Rates",
				"Simulated Reucrrence Rate (1/yr)");
		plot(outputDir, prefix+"_paleo_sect_rate", paleoProbSectRates, "Paleo-Detectable Simulated Subsection Rates",
				"Simulated Reucrrence Rate (1/yr)");
		plot(outputDir, prefix+"_raw_elem_rate", rawElemRates, "Raw Simulated Element Rates",
				"Simulated Reucrrence Rate (1/yr)");
		plot(outputDir, prefix+"_paleo_elem_rate", paleoProbElemRates, "Paleo-Detectable Simulated Element Rates",
				"Simulated Reucrrence Rate (1/yr)");
		
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("Paleoseismic Site Name", "UCERF3 Rate", "UCERF3 95% Conf", "UCERF3 68% Conf", "Sim Subsection Rate",
				"Sim Paleo-Detectable Subsection Rate", "Sim Element Rate", "Sim Paleo-Detectable Element Rate");
		for (PaleoResult result : results) {
			List<String> line = new ArrayList<>();
			PaleoRateConstraint constraint = result.constraint;
			line.add(constraint.getPaleoSiteName());
			line.add((float)constraint.getMeanRate()+"");
			line.add("["+(float)constraint.getLower95ConfOfRate()+" "+(float)result.constraint.getUpper95ConfOfRate()+"]");
			line.add("["+(float)constraint.getLower68ConfOfRate()+" "+(float)result.constraint.getUpper68ConfOfRate()+"]");
			line.add(rawSectRates.get(constraint).floatValue()+"");
			line.add(paleoProbSectRates.get(constraint).floatValue()+"");
			line.add(rawElemRates.get(constraint).floatValue()+"");
			line.add(paleoProbElemRates.get(constraint).floatValue()+"");
			csv.addLine(line);
		}
		csv.writeToFile(new File(outputDir, prefix+".csv"));
	}
	
	private void plot(File outputDir, String prefix, Map<PaleoRateConstraint, Double> rateMap, String title, String yAxisLabel)
			throws IOException {
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();

		DefaultXY_DataSet scatterInside68 = new DefaultXY_DataSet();
		DefaultXY_DataSet scatterInside95 = new DefaultXY_DataSet();
		DefaultXY_DataSet scatterOutside = new DefaultXY_DataSet();
		
		double logWhiskerDelta95 = 0.03;
		double logWhiskerDelta68 = 0.04;
		
		for (PaleoRateConstraint constr : rateMap.keySet()) {
			double rate = rateMap.get(constr);
			if (rate == 0d)
				continue;
			double paleoRate = constr.getMeanRate();
			
			double lower68 = constr.getLower68ConfOfRate();
			double upper68 = constr.getUpper68ConfOfRate();
			double lower95 = constr.getLower95ConfOfRate();
			double upper95 = constr.getUpper95ConfOfRate();
			if (rate >= lower68 && rate <= upper68)
				scatterInside68.set(paleoRate, rate);
			else if (rate >= lower95 && rate <= upper95)
				scatterInside95.set(paleoRate, rate);
			else
				scatterOutside.set(paleoRate, rate);
			
			double whiskerAbove95 = Math.pow(10, Math.log10(rate)+logWhiskerDelta95);
			double whiskerBelow95 = Math.pow(10, Math.log10(rate)-logWhiskerDelta95);
			
			DefaultXY_DataSet confRange95 = new DefaultXY_DataSet();
			confRange95.set(lower95, whiskerAbove95);
			confRange95.set(lower95, whiskerBelow95);
			confRange95.set(lower95, rate);
			confRange95.set(upper95, rate);
			confRange95.set(upper95, whiskerAbove95);
			confRange95.set(upper95, whiskerBelow95);
			
			funcs.add(confRange95);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.LIGHT_GRAY));
			
			double whiskerAbove68 = Math.pow(10, Math.log10(rate)+logWhiskerDelta68);
			double whiskerBelow68 = Math.pow(10, Math.log10(rate)-logWhiskerDelta68);
			
			DefaultXY_DataSet confRange68 = new DefaultXY_DataSet();
			confRange68.set(lower68, whiskerAbove68);
			confRange68.set(lower68, whiskerBelow68);
			confRange68.set(lower68, rate);
			confRange68.set(upper68, rate);
			confRange68.set(upper68, whiskerAbove68);
			confRange68.set(upper68, whiskerBelow68);
			
			funcs.add(confRange68);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
		}
		
		if (scatterInside68.size() > 0) {
			scatterInside68.setName("Inside 68% Conf");
			funcs.add(scatterInside68);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 7f, Color.GREEN.darker()));
		}
		
		if (scatterInside95.size() > 0) {
			scatterInside95.setName("Inside 95% Conf");
			funcs.add(scatterInside95);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 7f, Color.BLUE.darker()));
		}

		if (scatterOutside.size() > 0) {
			scatterOutside.setName("Outside 95% Conf");
			funcs.add(scatterOutside);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 7f, Color.BLACK));
		}
		
		double minVal = Double.POSITIVE_INFINITY;
		double maxVal = Double.NEGATIVE_INFINITY;
		if (scatterInside68.size() > 0) {
			minVal = Math.min(minVal, scatterInside68.getMinX());
			minVal = Math.min(minVal, scatterInside68.getMinY());
			maxVal = Math.max(maxVal, scatterInside68.getMaxX());
			maxVal = Math.max(maxVal, scatterInside68.getMaxY());
		}
		if (scatterInside95.size() > 0) {
			minVal = Math.min(minVal, scatterInside95.getMinX());
			minVal = Math.min(minVal, scatterInside95.getMinY());
			maxVal = Math.max(maxVal, scatterInside95.getMaxX());
			maxVal = Math.max(maxVal, scatterInside95.getMaxY());
		}
		if (scatterOutside.size() > 0) {
			minVal = Math.min(minVal, scatterOutside.getMinX());
			minVal = Math.min(minVal, scatterOutside.getMinY());
			maxVal = Math.max(maxVal, scatterOutside.getMaxX());
			maxVal = Math.max(maxVal, scatterOutside.getMaxY());
		}
		
		Range range = Double.isInfinite(minVal) ? new Range(1e-10, 1) :
			calcEncompassingLog10Range(minVal, maxVal);
		
		DefaultXY_DataSet oneToOne = new DefaultXY_DataSet();
		oneToOne.set(range.getLowerBound(), range.getLowerBound());
		oneToOne.set(range.getUpperBound(), range.getUpperBound());
		
		if (funcs.isEmpty()) {
			funcs.add(oneToOne);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		} else {
			funcs.set(0, oneToOne);
			chars.set(0, new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		}
		
		PlotSpec plot = new PlotSpec(funcs, chars, title, "UCERF3 Paleoseismic Reucrrence Rate (1/yr)", yAxisLabel);
		plot.setLegendVisible(true);
		
		HeadlessGraphPanel gp = getGraphPanel();
		gp.drawGraphPanel(plot, true, true, range, range);
		gp.getChartPanel().setSize(getPlotWidth(), getPlotHeight());
		gp.saveAsPNG(new File(getOutputDir(), prefix+".png").getAbsolutePath());
		gp.saveAsPDF(new File(getOutputDir(), prefix+".pdf").getAbsolutePath());
	}

	@Override
	public Collection<SimulatorElement> getApplicableElements() {
		return applicableElements;
	}
	
	public static void main(String[] args) throws IOException {
		File baseDir = new File("/home/kevin/Simulators/catalogs");
		
		D = true;
		
		double skipYears = 5000;
		
		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance(baseDir);
		
		File outputDir = new File("/tmp");
		
		PaleoRecurrencePlot plot = new PaleoRecurrencePlot(catalog.getElements(), catalog.getSubSectMapper());
		plot.initialize(catalog.getName(), outputDir, "paleo_recurrence");
		
		for (RSQSimEvent e : catalog.loader().skipYears(skipYears).iterable())
			plot.processEvent(e);
		
		System.out.println("Finalizing plot...");
		
		plot.finalizePlot();

		System.out.println("DONE");
	}

}
