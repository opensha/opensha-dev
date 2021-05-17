package scratch.kevin.simulators.multiFault;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.dom4j.DocumentException;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.axis.TickUnit;
import org.jfree.chart.axis.TickUnits;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.IDPairing;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.Vertex;
import org.opensha.sha.simulators.iden.LogicalAndRupIden;
import org.opensha.sha.simulators.iden.MagRangeRuptureIdentifier;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.iden.SkipYearsLoadIden;
import org.opensha.sha.simulators.parsers.RSQSimFileReader;
import org.opensha.sha.simulators.utils.RSQSimUtils;

import com.google.common.base.Preconditions;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.inversion.CommandLineInversionRunner;
import scratch.UCERF3.utils.DeformationModelFetcher;
import scratch.UCERF3.utils.FaultSystemIO;

public class RSQSimRupJumpCompare {
	
	private static final Color u3Color = Color.GRAY;
	private static final Color rsColor = Color.RED;

	public static void main(String[] args) throws IOException, DocumentException {
		File outputDir = new File("/home/kevin/Simulators/multiFault/2017_am");
		
		File u3SolFile = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/FM3_1_GEOL_MEAN_BRANCH_AVG_SOL.zip");
		FaultSystemSolution u3Sol = FaultSystemIO.loadSol(u3SolFile);
		
		boolean calcFixedJumpHist = false;
		boolean calcLengthHist = false;
		boolean calcMaxJumpHist = true;
		
		double jumpDist = 1d;
		double minMag = 6.5;
		FaultModels fm = FaultModels.FM3_1;
		DeformationModels dm = DeformationModels.GEOLOGIC;
		
//		List<File> rsGeomFiles = new ArrayList<>();
//		List<List<? extends SimulatorElement>> rsGeoms = new ArrayList<>();
//		List<SimElemDistCalculator> rsElemDists = new ArrayList<>();
//		List<List<? extends SimulatorEvent>> rsEvents = new ArrayList<>();
//		List<File> rsDistFiles = new ArrayList<>();
//		List<Integer> rsDistSizes = new ArrayList<>();
		
		int skipYears = 5000;
		double minFractForInclusion = 0.5;
		
		FaultSystemSolution rsSol = FaultSystemIO.loadSol(new File("/data/kevin/simulators/catalogs/rundir2194_long/laugh_test/"
				+ "rsqsim_sol_m6.5_skip5000_sectArea0.2.zip"));
		String rsName = "Shaw 2194";
		
//		if (calcMaxJumpHist || calcFixedJumpHist) {
//			// need to load events
//			for (int i=0; i<rsGeomFiles.size(); i++) {
//				File geomFile = rsGeomFiles.get(i);
//				System.out.println("Loading geometry from: "+geomFile.getAbsolutePath());
//				List<SimulatorElement> geom = RSQSimFileReader.readGeometryFile(geomFile, 11, 's');
//				if (geom.get(0).getFaultID() < 0)
//					RSQSimUtils.populateFaultIDWithParentIDs(geom, rsSols.get(i).getRupSet().getFaultSectionDataList());
//				RuptureIdentifier loadIden = new LogicalAndRupIden(new MagRangeRuptureIdentifier(minMag, 10),
//						new SkipYearsLoadIden(skipYears));
//				List<RuptureIdentifier> loadIdens = new ArrayList<>();
//				loadIdens.add(loadIden);
//				System.out.println("Loading events");
//				List<RSQSimEvent> events = RSQSimFileReader.readEventsFile(geomFile.getParentFile(), geom, loadIdens);
//				rsGeoms.add(geom);
//				rsEvents.add(events);
//				SimElemDistCalculator elemDistances;
//				File distFile = new File(geomFile.getAbsolutePath()+".distances");
//				rsDistFiles.add(distFile);
//				if (distFile.exists()) {
//					System.out.println("Loading distances");
//					elemDistances = new SimElemDistCalculator(DeformationModelFetcher.readMapFile(distFile), geom);
//				} else {
//					System.out.println("Calculating distances");
//					elemDistances = precalcElemDistances(geom);
//					DeformationModelFetcher.writeMapFile(elemDistances.distances, distFile);
//				}
//				rsElemDists.add(elemDistances);
//				rsDistSizes.add(elemDistances.distances.size());
//				System.out.println("DONE");
//			}
//		}
		
		File scratchDir = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/");
		DeformationModelFetcher fetch = new DeformationModelFetcher(fm, dm, scratchDir, 0.1);
		
		Map<IDPairing, Double> distances = fetch.getSubSectionDistanceMap(1000d);
		
		if (calcFixedJumpHist) {
			plotFixedJumpDist(u3Sol, distances, rsSol, rsName, minMag, jumpDist, outputDir);
		}
		
		if (calcLengthHist) {
			// now length hists
		}
		
		if (calcMaxJumpHist) {
			plotMaxJumpHist(u3Sol, distances, rsSol, rsName, minMag, outputDir);
		}
		
//		if (!rsDistFiles.isEmpty()) {
//			for (int i=0; i<rsDistFiles.size(); i++) {
//				SimElemDistCalculator dists = rsElemDists.get(i);
//				File file = rsDistFiles.get(i);
//				int prev = rsDistSizes.get(i);
//				if (prev >= dists.distances.size())
//					continue;
//				System.out.println("Writing "+dists.distances.size()+" distances (was "+prev+")to "+file.getAbsolutePath());
//				DeformationModelFetcher.writeMapFile(dists.distances, file);
//			}
//		}
	}
	
//	private static ArbitrarilyDiscretizedFunc calcActualFixedJumpHist(List<? extends SimulatorEvent> events,
//			Map<IDPairing, Double> elemDistances, double jumpDist) {
//		ArbitrarilyDiscretizedFunc distFunc = new ArbitrarilyDiscretizedFunc();
//		
//		for (SimulatorEvent e : events) {
//			
//		}
//		
////		hist.normalizeBySumOfY_Vals();
//		
//		return distFunc;
//	}
	
//	private static ArbitrarilyDiscretizedFunc calcNumJumpsFunc(FaultSystemSolution sol, )
	
	public static void plotFixedJumpDist(FaultSystemSolution u3Sol, Map<IDPairing, Double> distances,
			FaultSystemSolution rsSol, String rsName, double minMag, double jumpDist, File outputDir) throws IOException {
		DiscretizedFunc u3Func = CommandLineInversionRunner.getJumpFuncs(u3Sol, distances, jumpDist, minMag, null)[0];
		u3Func.scale(1d/u3Func.calcSumOfY_Vals());
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		funcs.add(u3Func);
		u3Func.setName("UCERF3");
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, u3Color));
		
		DiscretizedFunc rsFunc = CommandLineInversionRunner.getJumpFuncs(rsSol, distances, jumpDist, minMag, null)[0];
		rsFunc.scale(1d/rsFunc.calcSumOfY_Vals());
		rsFunc.setName(rsName);
		funcs.add(rsFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, rsColor));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "M≥"+(float)minMag+" Jump Comparison",
				"Num Jumps ≥"+(float)jumpDist+"km", "Fraction (Rate-Weighted)");
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setBackgroundColor(Color.WHITE);
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		
		String prefix = new File(outputDir, "jumps_"+(float)jumpDist+"km").getAbsolutePath();
		
		gp.drawGraphPanel(spec, false, false, null, new Range(0d, 1d));
		TickUnits tus = new TickUnits();
		TickUnit tu = new NumberTickUnit(1d);
		tus.add(tu);
		gp.getXAxis().setStandardTickUnits(tus);
		gp.getChartPanel().setSize(1000, 500);
		gp.saveAsPNG(prefix+".png");
		gp.saveAsPDF(prefix+".pdf");
		gp.saveAsTXT(prefix+".txt");
	}
	
	public static void plotMaxLengthHist(FaultSystemSolution u3Sol, FaultSystemSolution rsSol, String rsName,
			double minMag, File outputDir) throws IOException {
		double lenBinWidth = 50;
		double lenMax = 1250;
		HistogramFunction u3LengthHist = calcLengthHist(u3Sol, minMag, lenBinWidth, lenMax);
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		u3LengthHist.setName("UCERF3");
		funcs.add(u3LengthHist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, u3Color));
		
		HistogramFunction rsLengthHist = calcLengthHist(rsSol, minMag, lenBinWidth, lenMax);
		rsLengthHist.setName(rsName);
		funcs.add(rsLengthHist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, rsColor));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "M≥"+(float)minMag+" Length Comparison",
				"Rupture Length (km)", "Fraction (Rate-Weighted)");
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setBackgroundColor(Color.WHITE);
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		
		String prefix = new File(outputDir, "lenths").getAbsolutePath();
		
		gp.drawGraphPanel(spec, false, true, new Range(0d, lenMax), new Range(1e-6, 1));
		gp.getChartPanel().setSize(1000, 500);
		gp.saveAsPNG(prefix+".png");
		gp.saveAsPDF(prefix+".pdf");
		gp.saveAsTXT(prefix+".txt");
	}
	
	public static void plotMaxJumpHist(FaultSystemSolution u3Sol, Map<IDPairing, Double> distances,
			FaultSystemSolution rsSol, String rsName, double minMag, File outputDir) throws IOException {
		double deltaJump = 0.5;
		double maxJump = 30d;
		HistogramFunction u3Hist = calcMaxJumpDistHist(u3Sol, distances, maxJump, deltaJump);
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		u3Hist.setName("UCERF3");
		funcs.add(u3Hist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, u3Color));
		
		HistogramFunction rsHist = calcMaxJumpDistHist(rsSol, distances, maxJump, deltaJump);
		rsHist.setName(rsName+" Sub Sects");
		funcs.add(rsHist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, rsColor));
		
//		for (int i=0; i<rsSols.size(); i++) {
//			HistogramFunction rsHist = calcMaxJumpDistHist(rsSols.get(i), distances, maxJump, deltaJump);
//			rsHist.setName(rsNames.get(i)+" Sub Sects");
//			funcs.add(rsHist);
//			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, rsColors.get(i)));
//			
////			rsHist = calcActualSimulatorJumpDistHist(rsEvents.get(i), rsElemDists.get(i), maxJump, deltaJump);
////			rsHist.setName(rsNames.get(i)+" Actual");
////			funcs.add(rsHist);
////			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, rsColors.get(i)));
//		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, "M≥"+(float)minMag+" Max Jump Histogram",
				"Max Jump Dist Per Rupture (km)", "Fraction (Rate-Weighted)");
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setBackgroundColor(Color.WHITE);
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		
		String prefix = new File(outputDir, "max_jumps").getAbsolutePath();
		
		gp.drawGraphPanel(spec, false, true, new Range(0d, maxJump), new Range(1e-6, 1));
		gp.getChartPanel().setSize(1000, 500);
		gp.saveAsPNG(prefix+".png");
		gp.saveAsPDF(prefix+".pdf");
		gp.saveAsTXT(prefix+".txt");
	}
	
	private static HistogramFunction calcLengthHist(FaultSystemSolution sol, double minMag, double binWidth, double maxLen) {
		HistogramFunction hist = HistogramFunction.getEncompassingHistogram(0d, maxLen, binWidth);
		
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		for (int r=0; r<rupSet.getNumRuptures(); r++) {
			if (rupSet.getMagForRup(r) >= minMag) {
				double len = 0d;
				for (FaultSection sect : rupSet.getFaultSectionDataForRupture(r))
					len += sect.getTraceLength();
				hist.add(hist.getClosestXIndex(len), sol.getRateForRup(r));
			}
		}
		
		hist.normalizeBySumOfY_Vals();
		
		return hist;
	}
	
	private static HistogramFunction calcMaxJumpDistHist(FaultSystemSolution sol, Map<IDPairing, Double> distances,
			double maxJump, double deltaJump) {
		HistogramFunction hist = HistogramFunction.getEncompassingHistogram(0, maxJump, deltaJump);
		
		FaultSystemRupSet rupSet = sol.getRupSet();
		for (int r=0; r<rupSet.getNumRuptures(); r++) {
			double rupMaxJump = 0d;
			for (double jump : getJumps(rupSet.getFaultSectionDataForRupture(r), distances))
				rupMaxJump = Math.max(rupMaxJump, jump);
			hist.add(hist.getClosestXIndex(rupMaxJump), sol.getRateForRup(r));
		}
		
		hist.normalizeBySumOfY_Vals();
		
		return hist;
	}
	
//	private static List<Double> getJumps(List<FaultSectionPrefData> rupture, Map<IDPairing, Double> distances) {
//		List<Double> jumps = new ArrayList<>();
//		
//		for (int i=1; i<rupture.size(); i++) {
//			if (rupture.get(i-1).getParentSectionId() != rupture.get(i).getParentSectionId()) {
//				// it's a jump
//				double jumpDist = Double.POSITIVE_INFINITY;
//				int myParent = rupture.get(i-1).getParentSectionId();
//				int s1 = rupture.get(i-1).getSectionId();
//				for (FaultSectionPrefData sect : rupture) {
//					if (myParent != sect.getParentSectionId()) {
//						Double dist = distances.get(new IDPairing(s1, sect.getSectionId()));
//						if (dist == null)
//							// longer than 1000km rupture
//							continue;
//						jumpDist = Math.min(jumpDist, dist);
//					}
//				}
//				Preconditions.checkState(Double.isFinite(jumpDist));
//				jumps.add(jumpDist);
//			}
//		}
//		return jumps;
//	}
	
	private static List<Double> getJumps(List<? extends FaultSection> rupture, Map<IDPairing, Double> distances) {
		List<Double> jumps = new ArrayList<>();
		
		// bin by parent section
		Map<Integer, List<FaultSection>> parentSects = new HashMap<>();
		for (FaultSection sect : rupture) {
			List<FaultSection> sects = parentSects.get(sect.getParentSectionId());
			if (sects == null) {
				sects = new ArrayList<>();
				parentSects.put(sect.getParentSectionId(), sects);
			}
			sects.add(sect);
		}
		
		if (parentSects.size() == 1)
			return jumps;
		
//		boolean debug = Math.random() < 0.001;
		boolean debug = false;
		if (debug) System.out.println("Debugging rupture with "+parentSects.size()+" parent sects");
		for (int parentID : parentSects.keySet()) {
			if (debug) System.out.println("\tDebugging Parent: "+parentSects.get(parentID).get(0).getParentSectionName());
			double minDist = Double.POSITIVE_INFINITY;
			for (int oParentID : parentSects.keySet()) {
				if (oParentID == parentID)
					continue;
				double myMinDist = Double.POSITIVE_INFINITY;
				for (FaultSection s1 : parentSects.get(parentID)) {
					for (FaultSection s2 : parentSects.get(oParentID)) {
						Double dist = distances.get(new IDPairing(s1.getSectionId(), s2.getSectionId()));
						if (dist == null)
							// longer than 1000km rupture
							continue;
						myMinDist = Math.min(myMinDist, dist);
					}
				}
				if (debug) System.out.println("\t\t"+parentSects.get(oParentID).get(0).getParentSectionName()+": "+myMinDist);
				minDist = Math.min(myMinDist, minDist);
				if ((float)minDist == 0f)
					break;
			}
			if (debug) System.out.println("\tMin dist: "+minDist);
			jumps.add(minDist);
		}
		
		return jumps;
	}
	
//	private static SimElemDistCalculator precalcElemDistances(List<SimulatorElement> elems) {
//		Map<IDPairing, Future<Double>> futures = new HashMap<>();
//		
//		ExecutorService exec = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors()-1);
//		
//		Map<Integer, List<SimulatorElement>> elemsByParent = new HashMap<>();
//		for (SimulatorElement elem : elems) {
//			List<SimulatorElement> parentElems = elemsByParent.get(elem.getFaultID());
//			if (parentElems == null) {
//				parentElems = new ArrayList<>();
//				elemsByParent.put(elem.getFaultID(), parentElems);
//			}
//			parentElems.add(elem);
//		}
//		
//		System.out.println("Precalculating distances within "+elemsByParent.size()+" parent sects");
//		
//		for (int parentID : elemsByParent.keySet()) {
//			Map<IDPairing, Future<Double>> parentFutures = new HashMap<>();
//			List<SimulatorElement> parentElems = elemsByParent.get(parentID);
//			for (int i=0; i<parentElems.size(); i++) {
//				SimulatorElement e1 = parentElems.get(i);
//				for (int j=i+1; j<parentElems.size(); j++) {
//					SimulatorElement e2 = parentElems.get(j);
//					
//					if (Math.abs(e1.getSectionID() - e2.getSectionID()) > 1.5)
//						// only neighboring sections
//						continue;
//					parentFutures.put(new IDPairing(e1.getID(), e2.getID()), exec.submit(new SimDistCallable(e1, e2)));
//				}
//			}
//			System.out.println("Submitted "+parentFutures.size()+" futures for parent "+parentID);
//			futures.putAll(parentFutures);
//		}
//		
//		System.out.println("Submitted "+futures.size()+" futures");
//		
//		Map<IDPairing, Double> distances = new HashMap<>();
//		
//		for (IDPairing pair : futures.keySet()) {
//			if (distances.size() % 1000 == 0)
//				System.out.println("Calculated "+distances.size()+"/"+futures.size()
//					+" distances ("+(100d*distances.size()/futures.size())+" %)");
//			try {
//				distances.put(pair, futures.get(pair).get());
//			} catch (Exception e) {
//				throw ExceptionUtils.asRuntimeException(e);
//			}
//		}
//		
//		return new SimElemDistCalculator(distances, elems);
//	}
//	
//	private static class SimElemDistCalculator {
//		
//		private Map<IDPairing, Double> distances;
//		private Map<Integer, SimulatorElement> elemsMap;
//
//		public SimElemDistCalculator(Map<IDPairing, Double> distances, List<SimulatorElement> elems) {
//			this.distances = distances;
//			elemsMap = new HashMap<>();
//			for (SimulatorElement elem : elems)
//				elemsMap.put(elem.getID(), elem);
//		}
//		
//		public double get(int id1, int id2) {
//			return get(new IDPairing(id1, id2));
//		}
//		
//		public synchronized double get(IDPairing pair) {
//			Double dist = distances.get(pair);
//			if (dist == null) {
//				try {
//					dist = new SimDistCallable(elemsMap.get(pair.getID1()), elemsMap.get(pair.getID2())).call();
//				} catch (Exception e) {
//					ExceptionUtils.throwAsRuntimeException(e);
//				}
//				distances.put(pair, dist);
//			}
//			return dist;
//		}
//	}
//	
//	private static class SimDistCallable implements Callable<Double> {
//		
//		private SimulatorElement e1, e2;
//		
//		public SimDistCallable(SimulatorElement e1, SimulatorElement e2) {
//			super();
//			this.e1 = e1;
//			this.e2 = e2;
//		}
//
//		@Override
//		public Double call() throws Exception {
//			double dist = Double.POSITIVE_INFINITY;
//			for (Vertex v1 : e1.getVertices())
//				for (Vertex v2 : e2.getVertices())
//					dist = Double.min(dist, LocationUtils.linearDistanceFast(v1, v2));
//			return dist;
//		}
//		
//	}
//	
//	private static HistogramFunction calcActualSimulatorJumpDistHist(List<? extends SimulatorEvent> events,
//			SimElemDistCalculator elemDistances, double maxJump, double deltaJump) {
//		HistogramFunction hist = HistogramFunction.getEncompassingHistogram(0, maxJump, deltaJump);
//		
//		for (SimulatorEvent e : events) {
//			double rupMaxJump = 0d;
//			int[] elemIDs = e.getAllElementIDs();
//			for (int i=0; i<elemIDs.length; i++)
//				rupMaxJump = Math.max(rupMaxJump, calcClosestOtherElemDist(elemDistances, elemIDs[i], elemIDs));
//			hist.add(hist.getClosestXIndex(rupMaxJump), 1d);
//		}
//		
//		hist.normalizeBySumOfY_Vals();
//		
//		return hist;
//	}
//	
//	private static double calcClosestOtherElemDist(SimElemDistCalculator distCalc, int elemID, int[] elemIDs) {
//		if (elemIDs.length == 1)
//			return 0d;
//		double minDist = Double.POSITIVE_INFINITY;
//		for (int other : elemIDs) {
//			if (other == elemID)
//				continue;
//			minDist = Double.min(minDist, distCalc.get(elemID, other));
//			if ((float)minDist == 0f)
//				return 0;
//		}
//		return minDist;
//	}

}
