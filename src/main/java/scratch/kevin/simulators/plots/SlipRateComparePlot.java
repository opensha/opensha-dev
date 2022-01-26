package scratch.kevin.simulators.plots;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;
import org.dom4j.DocumentException;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.data.Range;
import org.jfree.chart.ui.TextAnchor;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import scratch.UCERF3.U3SlipEnabledSolution;
import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.inversion.InversionFaultSystemRupSet;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;
import scratch.UCERF3.utils.U3FaultSystemIO;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class SlipRateComparePlot extends AbstractPlot {
	
	private RSQSimSubSectionMapper mapper;
	private FaultModels fm;
	private DeformationModels dm;
	
	private Map<SimulatorElement, Double> elemTotalSlipsMap;
	private U3SlipEnabledSolution compSol;

	public SlipRateComparePlot(RSQSimSubSectionMapper mapper, FaultModels fm, DeformationModels dm) {
		this(mapper, fm, dm, null);
	}

	public SlipRateComparePlot(RSQSimSubSectionMapper mapper, FaultModels fm, DeformationModels dm, FaultSystemSolution compSol) {
		this.mapper = mapper;
		this.fm = fm;
		this.dm = dm;
		if (compSol != null) {
			Preconditions.checkState(mapper.getSubSections().size() == compSol.getRupSet().getFaultSectionDataList().size());
			if (compSol instanceof U3SlipEnabledSolution) {
				this.compSol = (U3SlipEnabledSolution)compSol;
			} else {
				// make it slip enabled, assuming mean UCERF3
				U3LogicTreeBranch branch = U3LogicTreeBranch.getMEAN_UCERF3(fm, dm);
				InversionFaultSystemRupSet iRupSet = new InversionFaultSystemRupSet(compSol.getRupSet(),
						branch, null, null, null, null, null);
				this.compSol = new InversionFaultSystemSolution(iRupSet, compSol.getRateForAllRups());
			}
		}
		
		mapper.trackSlipOnSections();
		
		elemTotalSlipsMap = new HashMap<>();
	}

	@Override
	protected void doProcessEvent(SimulatorEvent e) {
		ArrayList<SimulatorElement> elems = e.getAllElements();
		double[] slips = e.getAllElementSlips();
		
		for (int i=0; i<elems.size(); i++) {
			SimulatorElement elem = elems.get(i);
			double slip = slips[i];
			double totSlip = elemTotalSlipsMap.containsKey(elem) ? elemTotalSlipsMap.get(elem) : 0d;
			totSlip += slip;
			elemTotalSlipsMap.put(elem, totSlip);
		}
	}

	@Override
	public void finalizePlot() throws IOException {
		List<? extends FaultSection> subSects = mapper.getSubSections();
		double[] simSlipRates = new double[subSects.size()];
		double[] simTargetSlipRates = new double[subSects.size()];
		
		double[] u3TargetSlipRates = new double[subSects.size()];
		double[] u3SolSlipRates = compSol == null ? null : new double[subSects.size()];
		
		double[] ratioSimToTarget = new double[subSects.size()];
		double[] ratioSimToU3 = new double[subSects.size()];
		double[] ratioU3SolToTarget = compSol == null ? null : new double[subSects.size()];
		
		double durationYears = getCurrentDurationYears();
		List<LocationList> faults = new ArrayList<>();
		for (int i=0; i<subSects.size(); i++) {
			FaultSection subSect = subSects.get(i);
			if (u3SolSlipRates != null) {
				u3SolSlipRates[i] = 1e3 * compSol.calcSlipRateForSect(i);
				u3TargetSlipRates[i] = 1e3 * compSol.getRupSet().getSlipRateForSection(i);
			} else {
				u3TargetSlipRates[i] = subSect.getReducedAveSlipRate();
			}
			
			if (mapper.isMapped(subSect)) {
				double totArea = 0d;
				double totAreaWeightedSlip = 0d;
				double totAreaWeightedTargets = 0d;
				for (SimulatorElement elem : mapper.getSlipSectionsForSect(subSect)) {
					double totSlip = elemTotalSlipsMap.containsKey(elem) ? elemTotalSlipsMap.get(elem) : 0d;
					double slipRate = totSlip/durationYears;
					double area = elem.getArea();
					totAreaWeightedTargets += elem.getSlipRate()*area;
					totAreaWeightedSlip += slipRate*area;
					totArea += area;
				}
				// convert these to mm/yr
				simSlipRates[i] = 1e3 * totAreaWeightedSlip / totArea;
				simTargetSlipRates[i] = 1e3 * totAreaWeightedTargets / totArea;
				if (u3SolSlipRates != null)
					ratioU3SolToTarget[i] = u3SolSlipRates[i] / u3TargetSlipRates[i];
				ratioSimToTarget[i] = simSlipRates[i] / simTargetSlipRates[i];
				ratioSimToU3[i] = simSlipRates[i] / u3TargetSlipRates[i];
			} else {
				if (u3SolSlipRates != null)
					ratioU3SolToTarget[i] = Double.NaN;
				ratioSimToTarget[i] = Double.NaN;
				ratioSimToU3[i] = Double.NaN;
			}
			faults.add(subSect.getFaultTrace());
		}

		plotMaps(simSlipRates, u3TargetSlipRates, u3SolSlipRates, ratioSimToTarget, ratioSimToU3, ratioU3SolToTarget, faults);
		plotFaults(simSlipRates, simTargetSlipRates, u3TargetSlipRates, u3SolSlipRates);
	}

	private void plotMaps(double[] simSlipRates, double[] u3TargetSlipRates, double[] u3SolutionSlipRates, double[] ratioSimToTarget,
			double[] ratioSimToU3, double[] ratioU3SolToTarget, List<LocationList> faults) throws IOException {
		CPT slipCPT = FaultBasedMapGen.getLog10_SlipRateCPT();
		slipCPT.setNanColor(Color.GRAY);
//		CPT ratioCPT = GMT_CPT_Files.GMT_POLAR.instance().rescale(-2d, 2d);
//		CPT newRatioCPT = new CPT();
//		for (CPTVal val : ratioCPT)
//			newRatioCPT.add(new CPTVal(val.start, val.minColor.darker(), val.end, val.maxColor.darker()));
//		newRatioCPT.setNanColor(Color.GRAY);
//		ratioCPT = newRatioCPT;
		CPT ratioCPT = new CPT(-2d, 2d, Color.BLUE.darker().darker(), Color.BLUE.brighter(), new Color(230, 230, 230),
				Color.RED.brighter(), Color.RED.darker().darker());
		ratioCPT.setNanColor(Color.GRAY);
		
		Region region = new CaliforniaRegions.RELM_TESTING();
		
//		System.out.println("Building maps");
		GMT_Map simSlipMap = FaultBasedMapGen.buildMap(slipCPT, faults, FaultBasedMapGen.log10(simSlipRates),
				null, 1d, region, false, "Log@-10@- Simulator Output Slip Rates (mm/yr)");
		GMT_Map u3TargetSlipMap = FaultBasedMapGen.buildMap(slipCPT, faults, FaultBasedMapGen.log10(u3TargetSlipRates),
				null, 1d, region, false,
				"Log@-10@- UCERF3 "+fm.getName().replace("Fault Model", "FM")+" "+dm.getShortName()+" Target Slip Rates (mm/yr)");
		GMT_Map u3SolutionSlipMap = u3SolutionSlipRates == null ? null : FaultBasedMapGen.buildMap(
				slipCPT, faults, FaultBasedMapGen.log10(u3SolutionSlipRates), null, 1d, region, false,
				"Log@-10@- UCERF3 "+fm.getName().replace("Fault Model", "FM")+" "+dm.getShortName()+" Solution Slip Rates (mm/yr)");
		GMT_Map simRatioMap = FaultBasedMapGen.buildMap(ratioCPT, faults, FaultBasedMapGen.log10(ratioSimToTarget),
				null, 1d, region, false, "Log@-10@- Simulator Slip Rate / Target");
		GMT_Map simRatioU3Map = FaultBasedMapGen.buildMap(ratioCPT, faults, FaultBasedMapGen.log10(ratioSimToU3),
				null, 1d, region, false, "Log@-10@- Simulator Slip Rate / UCERF3");
		GMT_Map u3RatioMap = ratioU3SolToTarget == null ? null : FaultBasedMapGen.buildMap(
				ratioCPT, faults, FaultBasedMapGen.log10(ratioU3SolToTarget),
				null, 1d, region, false, "Log@-10@- UCERF3 Solution / Target");
		
		String prefix = getOutputPrefix();
		
//		System.out.println("Plotting");
		try {
			FaultBasedMapGen.plotMap(getOutputDir(), prefix+"_sim_map", false, simSlipMap);
			FaultBasedMapGen.plotMap(getOutputDir(), prefix+"_u3_target_map", false, u3TargetSlipMap);
			if (u3SolutionSlipMap != null)
				FaultBasedMapGen.plotMap(getOutputDir(), prefix+"_u3_sol_map", false, u3SolutionSlipMap);
			FaultBasedMapGen.plotMap(getOutputDir(), prefix+"_sim_ratio_map", false, simRatioMap);
			FaultBasedMapGen.plotMap(getOutputDir(), prefix+"_sim_u3_ratio_map", false, simRatioU3Map);
			if (u3RatioMap != null)
				FaultBasedMapGen.plotMap(getOutputDir(), prefix+"_u3_ratio_map", false, u3RatioMap);
		} catch (GMT_MapException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
	}
	
	private void plotFaults(double[] simSlipRates, double[] simSlipTargets, double[] u3TargetSlipRates, double[] u3SolSlipRates)
			throws IOException {
		Map<String, List<Integer>> namedFaults = fm.getNamedFaultsMapAlt();
		
		List<? extends FaultSection> allSubSects = mapper.getSubSections();
		
		for (String fault : namedFaults.keySet()) {
			HashSet<Integer> parentIDs = new HashSet<>(namedFaults.get(fault));
			
			MinMaxAveTracker latTrack = new MinMaxAveTracker();
			MinMaxAveTracker lonTrack = new MinMaxAveTracker();
			
			Map<Integer, List<FaultSection>> sectionsForFault = new HashMap<>();
			
			for (int i=0; i<allSubSects.size(); i++) {
				FaultSection subSect= allSubSects.get(i);
				if (parentIDs.contains(subSect.getParentSectionId())) {
					Location first = subSect.getFaultTrace().first();
					Location last = subSect.getFaultTrace().last();
					latTrack.addValue(first.getLatitude());
					latTrack.addValue(last.getLatitude());
					lonTrack.addValue(first.getLongitude());
					lonTrack.addValue(last.getLongitude());
					
					List<FaultSection> subSectsForFault = sectionsForFault.get(subSect.getParentSectionId());
					if (subSectsForFault == null) {
						subSectsForFault = new ArrayList<>();
						sectionsForFault.put(subSect.getParentSectionId(), subSectsForFault);
					}
					subSectsForFault.add(subSect);
				}
			}
			
			if (sectionsForFault.size() < 3)
				continue;
			
			// find common prefix if any
			List<String> parentNames = Lists.newArrayList();
			for (Integer parentID : namedFaults.get(fault)) {
				List<FaultSection> sectionsForParent = sectionsForFault.get(parentID);
				if (sectionsForParent == null)
					continue;
				String parentName = sectionsForParent.get(0).getParentSectionName();
				parentNames.add(parentName);
			}
			String[] parentNamesArray = parentNames.toArray(new String[0]);
			String commonPrefix = StringUtils.getCommonPrefix(parentNamesArray);
			
			double deltaLat = latTrack.getMax() - latTrack.getMin();
			double deltaLon = lonTrack.getMax() - lonTrack.getMin();
			boolean latitudeX = deltaLat > 0.5*deltaLon; // heavily favor latitude x
			Range xRange;
			if (latitudeX)
				xRange = new Range(latTrack.getMin(), latTrack.getMax());
			else
				xRange = new Range(lonTrack.getMin(), lonTrack.getMax());
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			PlotCurveCharacterstics simTargetChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK);
			PlotCurveCharacterstics simSolChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE);
			PlotCurveCharacterstics u3TargetChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED.darker());
			PlotCurveCharacterstics u3SolChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GREEN.darker());
			
			for (List<FaultSection> parentSections : sectionsForFault.values()) {
				DefaultXY_DataSet simFunc = new DefaultXY_DataSet();
				DefaultXY_DataSet simTargetFunc = new DefaultXY_DataSet();
				DefaultXY_DataSet u3TargetFunc = new DefaultXY_DataSet();
				DefaultXY_DataSet u3SolFunc = u3SolSlipRates == null ? null : new DefaultXY_DataSet();
				for (FaultSection sect : parentSections) {
					int index = sect.getSectionId(); // 0-indexed
					
					Location first = sect.getFaultTrace().first();
					Location last = sect.getFaultTrace().last();
					double x0 = latitudeX ? first.getLatitude() : first.getLongitude();
					double x1 = latitudeX ? last.getLatitude() : last.getLongitude();
					
					double simVal = simSlipRates[index];
					double simTargetVal = simSlipTargets[index];
					double u3TargetVal = u3TargetSlipRates[index];
					
					simFunc.set(x0, simVal);
					simFunc.set(x1, simVal);
					
					simTargetFunc.set(x0, simTargetVal);
					simTargetFunc.set(x1, simTargetVal);
					
					u3TargetFunc.set(x0, u3TargetVal);
					u3TargetFunc.set(x1, u3TargetVal);
					
					if (u3SolFunc != null) {
						double u3SolVal = u3SolSlipRates[index];
						
						u3SolFunc.set(x0, u3SolVal);
						u3SolFunc.set(x1, u3SolVal);
					}
				}
				
				if (funcs.isEmpty()) {
					simFunc.setName("Simulator Output");
					simTargetFunc.setName("Simulator Target");
					u3TargetFunc.setName("UCERF3 "+fm.getName().replace("Fault Model", "FM")+" "+dm.getShortName()+" Target");
					if (u3SolFunc != null)
						u3SolFunc.setName("UCERF3 Solution");
				}
				
				funcs.add(simFunc);
				chars.add(simSolChar);
				
				funcs.add(simTargetFunc);
				chars.add(simTargetChar);
				
				funcs.add(u3TargetFunc);
				chars.add(u3TargetChar);
				
				if (u3SolFunc != null) {
					funcs.add(u3SolFunc);
					chars.add(u3SolChar);
				}
			}
			
			double minY = Double.POSITIVE_INFINITY;
			double maxY = 0d;
			for (XY_DataSet func : funcs) {
				for (Point2D pt : func) {
					if (Double.isFinite(pt.getY()) && pt.getY() > 0d) {
						minY = Math.min(minY, pt.getY());
						maxY = Math.max(maxY, pt.getY());
					}
				}
			}
			Range yRange = calcEncompassingLog10Range(minY, maxY);
			
			// now add separators and labels for parent sections
			PlotCurveCharacterstics sepChar = new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.GRAY);
			Font annFont = new Font(Font.SANS_SERIF, Font.PLAIN, 20);
			List<XYTextAnnotation> anns = new ArrayList<>();
			for (List<FaultSection> parentSects : sectionsForFault.values()) {
				String annotationName = parentSects.get(0).getParentSectionName();
				if (!commonPrefix.isEmpty())
					annotationName = annotationName.substring(commonPrefix.length());
				annotationName = annotationName.replaceAll("San Andreas", "");
				annotationName = annotationName.replaceAll("Elsinore", "");
				annotationName = annotationName.replaceAll("\\(", "").replaceAll("\\)", "");
				annotationName = annotationName.replaceAll("2011 CFM", "");
				annotationName = annotationName.replaceAll(" rev", "");
				annotationName = annotationName.trim();
				
				Location firstLoc = parentSects.get(0).getFaultTrace().first();
				Location lastLoc = parentSects.get(parentSects.size()-1).getFaultTrace().last();
				
				double startX = latitudeX ? firstLoc.getLatitude() : firstLoc.getLongitude();
				double endX = latitudeX ? lastLoc.getLatitude() : lastLoc.getLongitude();
				double midX = 0.5*(startX + endX);
				
				DefaultXY_DataSet xy = new DefaultXY_DataSet();
				xy.set(startX, yRange.getLowerBound());
				xy.set(startX, yRange.getUpperBound());
				funcs.add(xy);
				chars.add(sepChar);
				xy = new DefaultXY_DataSet();
				xy.set(endX, yRange.getLowerBound());
				xy.set(endX, yRange.getUpperBound());
				funcs.add(xy);
				chars.add(sepChar);
				
				// note - Y location will get reset in the plotting code
				if (!annotationName.isEmpty()) {
					XYTextAnnotation a = new XYTextAnnotation(" "+annotationName, midX, yRange.getLowerBound());
					a.setFont(annFont);
					a.setRotationAnchor(TextAnchor.CENTER_LEFT);
					a.setTextAnchor(TextAnchor.CENTER_LEFT);
					a.setRotationAngle(-0.5*Math.PI);
					anns.add(a);
				}
			}
			
			PlotSpec spec = new PlotSpec(funcs, chars, fault+" Slip Rates", latitudeX ? "Latitude (degrees)" : "Longitude (degrees)",
					"Slip Rate (mm/yr)");
			spec.setPlotAnnotations(anns);
			spec.setLegendVisible(true);
			
			HeadlessGraphPanel gp = getGraphPanel();
			gp.drawGraphPanel(spec, false, true, xRange, yRange);
			gp.getPlot().setDatasetRenderingOrder(DatasetRenderingOrder.REVERSE);
			gp.getChartPanel().setSize(1200, 800);
			gp.getPlot().getDomainAxis().setInverted(true);
			File outputDir = getOutputDir();
			String prefix = getOutputPrefix()+"_fault_"+fault.replaceAll("\\W+", "_");
			gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
			gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
		}
	}
	
	private static DefaultXY_DataSet buildSectXY(FaultSectionPrefData subSect, boolean latitudeX, double yVal) {
		DefaultXY_DataSet xy = new DefaultXY_DataSet();
		Location first = subSect.getFaultTrace().first();
		Location last = subSect.getFaultTrace().last();
		double x0 = latitudeX ? first.getLatitude() : first.getLongitude();
		double x1 = latitudeX ? last.getLatitude() : last.getLongitude();
		xy.set(x0, yVal);
		xy.set(x1, yVal);
		return xy;
	}

	@Override
	public Collection<SimulatorElement> getApplicableElements() {
		return null;
	}

	public static void main(String[] args) throws IOException, DocumentException {
		File baseDir = new File("/home/kevin/Simulators/catalogs");
		
		File compSolFile = new File("/home/kevin/.opensha/ucerf3_fm_dm_sols/FM3_1_GEOL_MEAN_BRANCH_AVG_SOL.zip");
		FaultSystemSolution compSol = U3FaultSystemIO.loadSol(compSolFile);
//		SlipEnabledSolution compSol = null;
		
		double skipYears = 5000;
		
		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance(baseDir);
		
		File outputDir = new File("/tmp");
		
		SlipRateComparePlot plot = new SlipRateComparePlot(catalog.getSubSectMapper(),
				catalog.getFaultModel(), catalog.getDeformationModel(), compSol);
		plot.initialize(catalog.getName(), outputDir, "slip_rate");
		
		for (RSQSimEvent e : catalog.loader().skipYears(skipYears).iterable())
			plot.processEvent(e);
		
		System.out.println("Finalizing plot...");
		
		plot.finalizePlot();

		System.out.println("DONE");
	}

}
