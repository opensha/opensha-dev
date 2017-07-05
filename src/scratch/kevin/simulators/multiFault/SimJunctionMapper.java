package scratch.kevin.simulators.multiFault;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.dom4j.DocumentException;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotElement;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.util.FaultUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.parsers.EQSIMv06FileReader;
import org.opensha.sha.simulators.utils.General_EQSIM_Tools;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.inversion.CommandLineInversionRunner;
import scratch.UCERF3.inversion.InversionFaultSystemRupSet;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.inversion.coulomb.CoulombRates;
import scratch.UCERF3.inversion.coulomb.CoulombRatesRecord;
import scratch.UCERF3.inversion.coulomb.CoulombRatesTester;
import scratch.UCERF3.inversion.laughTest.LaughTestFilter;
import scratch.UCERF3.utils.DeformationModelFetcher;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.IDPairing;
import scratch.UCERF3.utils.MatrixIO;
import scratch.UCERF3.utils.UCERF3_DataUtils;
import scratch.kevin.simulators.erf.SimulatorFaultSystemSolution;

public class SimJunctionMapper {
	
	private InversionFaultSystemSolution ucerfSol;
	private FaultSystemSolution simSol;
	
	private Map<IDPairing, Double> simDistances;
	private Map<IDPairing, Double> ucerfDistances;
	
	private Map<Integer, Integer> simToUCERF_parentsMap;
	private File resultsDir;
	private File plotsDir;
	
	// lowest always first
	private Map<IDPairing, Double> simParentTogetherRates;
	private Map<IDPairing, Double> normSimParentTogetherRates;
	private Map<IDPairing, Double> ucerfParentTogetherRates;
	private Map<IDPairing, Double> normUCERFParentTogetherRates;
	
	private List<double[]> simSlipsMean, simSlipsMax;
	
	private CoulombRates coulombRates;
	
	private Map<Integer, Double> ucerfParentSlipsMap;
	
	// mapped with UCERF parent IDs
	private Map<IDPairing, CoulombRatesRecord[]> parentCoulombs;
	
	private CoulombRatesTester tester;

	public SimJunctionMapper(InversionFaultSystemSolution ucerfSol, FaultSystemSolution simSol,
			List<double[]> simSlipsMean, List<double[]> simSlipsMax, CoulombRates coulombRates, File resultsDir)
			throws IOException {
		this.ucerfSol = ucerfSol;
		this.simSol = simSol;
		this.resultsDir = resultsDir;
		plotsDir = new File(resultsDir, "plots");
		if (!plotsDir.exists())
			plotsDir.mkdir();
		
		this.simSlipsMax = simSlipsMax;
		this.simSlipsMean = simSlipsMean;
		
		simDistances = loadDistances(simSol.getRupSet(), "rsqsim");
		// fill out simulator distances as they can get really screwy with "X" ruptures
		FaultSystemRupSet simRupSet = simSol.getRupSet();
		for (int i=0; i<simRupSet.getNumSections(); i++) {
			for (int j=0; j<simRupSet.getNumSections(); j++) {
				IDPairing pair = new IDPairing(i, j);
				if (!simDistances.containsKey(pair))
					simDistances.put(pair, Double.POSITIVE_INFINITY);
			}
		}
		ucerfDistances = loadDistances(ucerfSol.getRupSet(), "ucerf");
		
		mapParentSects();
		
		// only for sucessfully mapped pairs
		simParentTogetherRates = calcParentTogetherRates(simSol, new HashSet<Integer>(simToUCERF_parentsMap.keySet()));
		normSimParentTogetherRates = calcNormParentTogetherRates(simParentTogetherRates, simSol);
		ucerfParentTogetherRates = calcParentTogetherRates(ucerfSol, new HashSet<Integer>(simToUCERF_parentsMap.values()));
		normUCERFParentTogetherRates = calcNormParentTogetherRates(ucerfParentTogetherRates, ucerfSol);
		
		// now remap sim parents to ucerf parents
		simParentTogetherRates = remapSimParentTogetherRates(simParentTogetherRates);
		normSimParentTogetherRates = remapSimParentTogetherRates(normSimParentTogetherRates);
		
		this.coulombRates = coulombRates;
		parentCoulombs = Maps.newHashMap();
		
		FaultSystemRupSet ucerfRupSet = ucerfSol.getRupSet();
		for (CoulombRatesRecord rec : coulombRates.values()) {
			IDPairing subSectPair = rec.getPairing();
			int parent1 = ucerfRupSet.getFaultSectionData(subSectPair.getID1()).getParentSectionId();
			int parent2 = ucerfRupSet.getFaultSectionData(subSectPair.getID2()).getParentSectionId();
			if (parent1 == parent2)
				continue;
			IDPairing parentPair = getPairing(parent1, parent2);
			if (parentCoulombs.containsKey(parentPair))
				// already exists, put in 2nd slot
				parentCoulombs.get(parentPair)[1] = rec;
			else
				// first one, put in first slot
				parentCoulombs.put(parentPair, new CoulombRatesRecord[] { rec, null });
		}
		for (CoulombRatesRecord[] recs : parentCoulombs.values())
			Preconditions.checkState(recs[0] != null && recs[1] != null);
		System.out.println("Sim has "+simParentTogetherRates.size()+" parent connections");
		System.out.println("UCERF has "+ucerfParentTogetherRates.size()+" parent connections");
		System.out.println("Coulomb has "+parentCoulombs.size()+" parent connections");
		
		tester = LaughTestFilter.getDefault().getCoulombFilter();
		
		ucerfParentSlipsMap = calcParentSectMeanSlipRates(ucerfRupSet);
	}
	
	private Map<IDPairing, Double> loadDistances(FaultSystemRupSet rupSet, String prefix) throws IOException {
		File distFile = new File(resultsDir, "distances_"+prefix+"_"+rupSet.getNumRuptures()+".txt");
		if (distFile.exists())
			return DeformationModelFetcher.readMapFile(distFile);
		Map<IDPairing, Double> distances = DeformationModelFetcher.calculateDistances(50d, rupSet.getFaultSectionDataList());
		
		HashMap<IDPairing, Double> reversed = new HashMap<IDPairing, Double>();

		// now add the reverse distance
		for (IDPairing pair : distances.keySet()) {
			IDPairing reverse = pair.getReversed();
			reversed.put(reverse, distances.get(pair));
		}
		distances.putAll(reversed);
		
		DeformationModelFetcher.writeMapFile(distances, distFile);
		return distances;
	}
	
//	private static List<IDPairing> loadJunctions() {
//		return null;
//	}
	
	private void mapParentSects() throws IOException {
		simToUCERF_parentsMap = Maps.newHashMap();
		
		File matchCSVFile = new File(resultsDir, "parent_matches.csv");
		if (matchCSVFile.exists()) {
			// load it
			System.out.println("Loading precomputed parent sects...");
			CSVFile<String> csv = CSVFile.readFile(matchCSVFile, true);
			for (int row=1; row<csv.getNumRows(); row++) {
				List<String> line = csv.getLine(row);
				Integer simParent = Integer.parseInt(line.get(1));
				Integer ucerfParent = Integer.parseInt(line.get(3));
				simToUCERF_parentsMap.put(simParent, ucerfParent);
			}
			return;
		}
		File noMatchCSVFile = new File(resultsDir, "parent_nomatches.csv");
		System.out.println("Mapping parent sects...");
		
		List<ParentSectInfo> simParents = getParentSects(simSol.getRupSet());
		List<ParentSectInfo> ucerfParents = getParentSects(ucerfSol.getRupSet());
		
		double maxEndDistTol = 20d;
		double maxClosestDistTol = 1; // make sure that the trace is close. can be overridden if names match
		double maxDipDiff = 10;
		double maxRakeDiff = 30;
		
		CSVFile<String> matchCSV = new CSVFile<String>(true);
		CSVFile<String> noMatchCSV = new CSVFile<String>(true);
		List<String> header = Lists.newArrayList("Sim Name", "Sim ID", "UCERF3 Name",
				"UCERF3 ID", "Lengh Discrep", "End Loc Discrep", "Min Trace Dist", "Rake Diff", "Dip Diff");
		matchCSV.addLine(header);
		noMatchCSV.addLine(header);
		
		for (ParentSectInfo simParent : simParents) {
			ParentSectInfo closest = null;
			double closestDist = Double.MAX_VALUE;
			double traceDist = Double.NaN;
			for (ParentSectInfo ucerfParent : ucerfParents) {
				if (ucerfParent.parentID == 230)
					// hack for pathological case
					continue;
				double dist = simParent.calcMaxDist(ucerfParent);
				double myTraceDist = simParent.calcMinDiscrTraceDist(ucerfParent);
				double myRakeDiff = simParent.getRakeDiff(ucerfParent);
				double myDipDiff = simParent.getDipDiff(ucerfParent);
				boolean passTrace = myTraceDist <= maxClosestDistTol;
				if (!passTrace && myTraceDist < 5) {
					// lets check name
					passTrace = StringUtils.getCommonPrefix(simParent.parentName, ucerfParent.parentName).length() > 4;
					if (!passTrace)
						passTrace = StringUtils.getLevenshteinDistance(simParent.parentName, ucerfParent.parentName) < 4
								&& simParent.parentName.length() > 4;
				}
				if (dist < closestDist && passTrace
						&& myRakeDiff < maxRakeDiff && myDipDiff < maxDipDiff) {
					closest = ucerfParent;
					closestDist = dist;
					traceDist = myTraceDist;
				}
			}
			double distDiscrep;
			List<String> line;
			if (closest == null) {
				distDiscrep = Double.NaN;
				line = Lists.newArrayList(simParent.parentName, simParent.parentID+"",
						"", "", distDiscrep+"", Double.NaN+"", Double.NaN+"", Double.NaN+"", Double.NaN+"");
			} else {
				distDiscrep = Math.abs(closest.length - simParent.length);
				line = Lists.newArrayList(simParent.parentName, simParent.parentID+"",
						closest.parentName, closest.parentID+"", (float)distDiscrep+"",
						(float)closestDist+"", (float)traceDist+"",
						(float)simParent.getRakeDiff(closest)+"", (float)simParent.getDipDiff(closest)+"");
			}
			if (closestDist <= maxEndDistTol) {
				// it's a match
				simToUCERF_parentsMap.put(simParent.parentID, closest.parentID);
				matchCSV.addLine(line);
			} else {
				noMatchCSV.addLine(line);
			}
		}
		
		System.out.println("Matched "+simToUCERF_parentsMap.size()+"/"+simParents.size()+" parents.");
		matchCSV.writeToFile(matchCSVFile);
		noMatchCSV.writeToFile(noMatchCSVFile);
	}
	
	private List<ParentSectInfo> getParentSects(FaultSystemRupSet rupSet) {
		Map<Integer, List<FaultSectionPrefData>> map = Maps.newHashMap();
		for (FaultSectionPrefData sect : rupSet.getFaultSectionDataList()) {
			Integer parentID = sect.getParentSectionId();
			List<FaultSectionPrefData> sects = map.get(parentID);
			if (sects == null) {
				sects = Lists.newArrayList();
				map.put(parentID, sects);
			}
			sects.add(sect);
		}
		
		List<ParentSectInfo> parents = Lists.newArrayList();
		for (List<FaultSectionPrefData> sects : map.values())
			parents.add(new ParentSectInfo(sects));
		
		return parents;
	}
	
	private static Map<IDPairing, Double> calcParentTogetherRates(FaultSystemSolution sol, HashSet<Integer> parentsToConsider) {
		Map<IDPairing, Double> map = Maps.newHashMap();
		FaultSystemRupSet rupSet = sol.getRupSet();
		for (int r=0; r<rupSet.getNumRuptures(); r++) {
			List<Integer> parents = rupSet.getParentSectionsForRup(r);
			double rate = sol.getRateForRup(r);
			if (rate == 0)
				continue;
			for (int i=1; i<parents.size(); i++) {
				IDPairing pair = getPairing(parents.get(i-1), parents.get(i));
				if (!parentsToConsider.contains(pair.getID1()) || !parentsToConsider.contains(pair.getID2()))
					continue;
				Double val = map.get(pair);
				if (val == null)
					val = 0d;
				val += rate;
				map.put(pair, val);
			}
		}
		return map;
	}
	
	private static Map<IDPairing, Double> calcNormParentTogetherRates(Map<IDPairing, Double> togetherRates, FaultSystemSolution sol) {
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		Map<IDPairing, Double> normTogetherRates = Maps.newHashMap();
		
		MinMaxAveTracker sectsAtConnectTrack = new MinMaxAveTracker();
		
		for (IDPairing pair : togetherRates.keySet()) {
			double togetherRate = togetherRates.get(pair);
			// sections of parent 1 that connect to parent 2
			HashSet<Integer> sects1 = new HashSet<Integer>();
			// sections of parent 2 that connect to parent 1
			HashSet<Integer> sects2 = new HashSet<Integer>();
			for (int rupIndex : rupSet.getRupturesForParentSection(pair.getID1())) {
				List<Integer> parentsForRup = rupSet.getParentSectionsForRup(rupIndex);
				if (parentsForRup.contains(pair.getID2())) {
					// this is a multi fault with both parents
					List<FaultSectionPrefData> sects = rupSet.getFaultSectionDataForRupture(rupIndex);
					// find the connecting sub section
					for (int i=1; i<sects.size(); i++) {
						FaultSectionPrefData s1 = sects.get(i-1);
						FaultSectionPrefData s2 = sects.get(i);
						if (s1.getParentSectionId() == pair.getID1() && s2.getParentSectionId() == pair.getID2()) {
							sects1.add(s1.getSectionId());
							sects2.add(s2.getSectionId());
						} else if (s1.getParentSectionId() == pair.getID2() && s2.getParentSectionId() == pair.getID1()) {
							sects1.add(s2.getSectionId());
							sects2.add(s1.getSectionId());
						}
					}
				}
			}
			Preconditions.checkState(!sects1.isEmpty());
			Preconditions.checkState(!sects2.isEmpty());
			sectsAtConnectTrack.addValue(sects1.size());
			sectsAtConnectTrack.addValue(sects2.size());
			
			// total rate of section on parent 1 that connects to parent 2
			double s1Rate = getTotSectRates(sects1, sol);
			// total rate of section on parent 2 that connects to parent 1
			double s2Rate = getTotSectRates(sects2, sol);
			double avg = 0.5*(s1Rate+s2Rate);
			normTogetherRates.put(pair, togetherRate/avg);
		}
		System.out.println("Sections that connect: "+sectsAtConnectTrack);
		return normTogetherRates;
	}
	
	private static double getTotSectRates(HashSet<Integer> sects, FaultSystemSolution sol) {
		FaultSystemRupSet rupSet = sol.getRupSet();
		HashSet<Integer> rups = new HashSet<Integer>();
		for (int sect : sects)
			rups.addAll(rupSet.getRupturesForSection(sect));
		double rate = 0;
		for (int rupIndex : rups)
			rate += sol.getRateForRup(rupIndex);
		return rate;
	}
	
	private Map<IDPairing, Double> remapSimParentTogetherRates(Map<IDPairing, Double> togetherRates) {
		Map<IDPairing, Double> simRemappedParentTogetherRates = Maps.newHashMap();
		for (IDPairing pair : togetherRates.keySet()) {
			Integer mapped1 = simToUCERF_parentsMap.get(pair.getID1());
			Integer mapped2 = simToUCERF_parentsMap.get(pair.getID2());
			Preconditions.checkNotNull(mapped1);
			Preconditions.checkNotNull(mapped2);
			simRemappedParentTogetherRates.put(getPairing(mapped1, mapped2), togetherRates.get(pair));
		}
		return simRemappedParentTogetherRates;
	}
	
	private static IDPairing getPairing(int id1, int id2) {
		if (id1 < id2)
			return new IDPairing(id1, id2);
		return new IDPairing(id2, id1);
	}
	
	private static Map<Integer, Double> calcParentSectMeanSlipRates(FaultSystemRupSet rupSet) {
		Map<Integer, List<Double>> parentAllSlipRates = Maps.newHashMap();
		for (int sectIndex=0; sectIndex<rupSet.getNumSections(); sectIndex++) {
			FaultSectionPrefData sect = rupSet.getFaultSectionData(sectIndex);
			double slip = sect.getOrigAveSlipRate();
			List<Double> slips = parentAllSlipRates.get(sect.getParentSectionId());
			if (slips == null) {
				slips = Lists.newArrayList();
				parentAllSlipRates.put(sect.getParentSectionId(), slips);
			}
			slips.add(slip);
		}
		Map<Integer, Double> parentSlipRates = Maps.newHashMap();
		for (Integer parentID : parentAllSlipRates.keySet()) {
			List<Double> parentRates = parentAllSlipRates.get(parentID);
			double meanRate = 0;
			for (double rate : parentRates)
				meanRate += rate;
			meanRate /= parentRates.size();
			parentSlipRates.put(parentID, meanRate);
		}
		return parentSlipRates;
	}
	
	private static final double trace_discr_km = 1d;
	
	private class ParentSectInfo {
		private int parentID;
		private String parentName;
		private double length;
		private Location loc1; // always the northermonst
		private Location loc2; // always the southernmost
		List<Location> discrTrace = Lists.newArrayList();
		
		private double avgRake;
		private double avgDip;

		public ParentSectInfo(List<FaultSectionPrefData> sectsForParent) {
			this.parentID = sectsForParent.get(0).getParentSectionId();
			this.parentName = sectsForParent.get(0).getParentSectionName();
			
			length = 0;
			List<Double> rakes = Lists.newArrayList();
			for (FaultSectionPrefData sect: sectsForParent) {
				FaultTrace trace = sect.getFaultTrace();
				double subLen = trace.getTraceLength();
				length += subLen;
				int num = (int)Math.round(subLen/trace_discr_km);
				if (num < 2)
					num = 2;
				discrTrace.addAll(FaultUtils.resampleTrace(trace, num));
				avgDip += sect.getAveDip();
				rakes.add(sect.getAveRake());
			}
			avgDip /= (double)sectsForParent.size();
			avgRake = FaultUtils.getInRakeRange(FaultUtils.getAngleAverage(rakes));
			
			loc1 = sectsForParent.get(0).getFaultTrace().first();
			loc2 = sectsForParent.get(sectsForParent.size()-1).getFaultTrace().last();
			if (loc1.getLatitude() < loc2.getLatitude()) {
				// swap em
				Location tempLoc = loc1;
				loc1 = loc2;
				loc2 = tempLoc;
			}
		}
		
		public double calcMaxDist(ParentSectInfo o) {
			return Math.max(LocationUtils.horzDistanceFast(loc1, o.loc1), LocationUtils.horzDistanceFast(loc2, o.loc2));
		}
		
		public double calcMinDiscrTraceDist(ParentSectInfo o) {
			double val = Double.MAX_VALUE;
			for (Location loc1 : discrTrace)
				for (Location loc2 : o.discrTrace)
					val = Math.min(val, LocationUtils.horzDistance(loc1, loc2));
			return val;
		}
		
		public double getDipDiff(ParentSectInfo o) {
			return Math.abs(avgDip - o.avgDip);
		}
		
		public double getRakeDiff(ParentSectInfo o) {
			double r1 = avgRake;
			double r2 = o.avgRake;
			if (r1 < -90 && r2 > 90)
				r1 += 360;
			return Math.abs(r1 - r2);
		}
	}
	
	public void plotMultiRateVsCoulomb(Map<IDPairing, Double> parentTogetherRates, String name, boolean norm)
			throws IOException {
		// first make sure that every parent together rate has a coulomb mapping
		// nevermind, RSQ sim can connect in a few extra places, ignore check
//		for (IDPairing pair : parentTogetherRates.keySet())
//			Preconditions.checkNotNull(parentCoulombs.get(pair), "No coulomb for "+pair);
		DefaultXY_DataSet passXY = new DefaultXY_DataSet();
		DefaultXY_DataSet failXY = new DefaultXY_DataSet();
		
		for (IDPairing pair : parentCoulombs.keySet()) {
			CoulombRatesRecord[] recs = parentCoulombs.get(pair);
			double dcff = Math.max(recs[0].getCoulombStressChange(), recs[1].getCoulombStressChange());
			Double rate = parentTogetherRates.get(pair);
			if (rate == null)
				continue;
//				rate = 0d;
			if (tester.doesRupturePass(Lists.newArrayList(recs[0]), Lists.newArrayList(recs[1])))
				passXY.set(dcff, rate);
			else
				failXY.set(dcff, rate);
		}
		
		List<PlotElement> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		if (passXY.size() > 0) {
			funcs.add(passXY);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 6f, Color.GREEN));
		}
		
		if (failXY.size() > 0) {
			funcs.add(failXY);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 6f, Color.RED));
		}
		
		// now do regression
		SimpleRegression reg = new SimpleRegression(true);
		MinMaxAveTracker xTrack = new MinMaxAveTracker();
		for (Point2D pt : passXY) {
			if (pt.getX() > 0 && pt.getY() > 0) {
				reg.addData(Math.log10(pt.getX()), Math.log10(pt.getY()));
				xTrack.addValue(pt.getX());
			}
		}
		for (Point2D pt : failXY) {
			if (pt.getX() > 0 && pt.getY() > 0) {
				reg.addData(Math.log10(pt.getX()), Math.log10(pt.getY()));
				xTrack.addValue(pt.getX());
			}
		}
		ArbitrarilyDiscretizedFunc regFunc = new ArbitrarilyDiscretizedFunc();
		regFunc.setName("Best Fit Line. MSE="+reg.getMeanSquareError());
		regFunc.set(xTrack.getMin(), Math.pow(10, reg.predict(Math.log10(xTrack.getMin()))));
		regFunc.set(xTrack.getMax(), Math.pow(10, reg.predict(Math.log10(xTrack.getMax()))));
		
		funcs.add(0, regFunc);
		chars.add(0, new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		
		String title = name+" Junction Rate vs DCFF";
		if (norm)
			title = "Normalized "+title;
		
//		GraphWindow gw = new GraphWindow(funcs, title, chars);
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		String xAxisLabel = "Coulomb DCFF";
		String yAxisLabel = "Junction Corupture Rate";
		String fName = "junction_corupture_vs_coulomb_"+name.toLowerCase();
		if (norm) {
			yAxisLabel = "Normalized "+yAxisLabel;
			fName += "_norm";
		}
		PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
		gp.setBackgroundColor(Color.WHITE);
		CommandLineInversionRunner.setFontSizes(gp);
		gp.drawGraphPanel(spec, true, true);
		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPNG(new File(plotsDir, fName+".png").getAbsolutePath());
		gp.saveAsPDF(new File(plotsDir, fName+".pdf").getAbsolutePath());
	}
	
	public void plotSimVsUCERFParentRates(boolean norm) throws IOException {
		Map<IDPairing, Double> simRates;
		Map<IDPairing, Double> ucerfRates;
		if (norm) {
			simRates = normSimParentTogetherRates;
			ucerfRates = normUCERFParentTogetherRates;
		} else {
			simRates = simParentTogetherRates;
			ucerfRates = ucerfParentTogetherRates;
		}
		
		DefaultXY_DataSet xy = new DefaultXY_DataSet();
		
		HashSet<IDPairing> pairs = new HashSet<IDPairing>();
		pairs.addAll(simRates.keySet());
		pairs.addAll(ucerfRates.keySet());
		
		for (IDPairing pair : pairs) {
			Double simRate = simRates.get(pair);
			if (simRate == null)
//				simRate = 0d;
				continue;
			Double ucerfRate = ucerfRates.get(pair);
			if (ucerfRate == null)
//				ucerfRate = 0d;
				continue;
			xy.set(simRate, ucerfRate);
		}
		
		List<PlotElement> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		funcs.add(xy);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.BLACK));
		
		// now do regression
		MinMaxAveTracker xTrack = new MinMaxAveTracker();
		for (Point2D pt : xy)
			if (pt.getX() > 0 && pt.getY() > 0)
				xTrack.addValue(pt.getX());
		
		ArbitrarilyDiscretizedFunc oneToOne = new ArbitrarilyDiscretizedFunc();
		oneToOne.set(xTrack.getMin(), xTrack.getMin());
		oneToOne.set(xTrack.getMax(), xTrack.getMax());
		
		funcs.add(0, oneToOne);
		chars.add(0, new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		
		String title = " UCERF vs RSQSim Junction Rates";
		if (norm)
			title = "Normalized "+title;
		
//		GraphWindow gw = new GraphWindow(funcs, title, chars);
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		CommandLineInversionRunner.setFontSizes(gp);
		String xAxisLabel, yAxisLabel, fName;
		if (norm) {
			xAxisLabel = "Normalized RSQSim Junction Corupture Rate";
			yAxisLabel = "Normalized UCERF Junction Corupture Rate";
			fName = "sim_vs_ucerf_rates_norm";
		} else {
			xAxisLabel = "RSQSim Junction Corupture Rate";
			yAxisLabel = "UCERF Junction Corupture Rate";
			fName = "sim_vs_ucerf_rates";
		}
		PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
		gp.setBackgroundColor(Color.WHITE);
		gp.drawGraphPanel(spec, true, true);
		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPNG(new File(plotsDir, fName+".png").getAbsolutePath());
		gp.saveAsPDF(new File(plotsDir, fName+".pdf").getAbsolutePath());
	}
	
	public void plotSlipVsRate() throws IOException {
		DefaultXY_DataSet func = new DefaultXY_DataSet();
		
		FaultSystemRupSet rupSet = simSol.getRupSet();
		double[] aveSlips = new double[rupSet.getNumSections()];
		int[] sectCounts = new int[rupSet.getNumSections()];
		for (int rupIndex=0; rupIndex<rupSet.getNumRuptures(); rupIndex++) {
			double[] slips = simSlipsMax.get(rupIndex);
			List<Integer> sectIndexes = rupSet.getSectionsIndicesForRup(rupIndex);
			for (int i=0; i<sectIndexes.size(); i++) {
				double slip = slips[i];
				int sectIndex = sectIndexes.get(i);
				aveSlips[sectIndex] += slip;
				sectCounts[sectIndex]++;
			}
		}
		for (int sectIndex=0; sectIndex<rupSet.getNumSections(); sectIndex++) {
			double slipRate = rupSet.getSlipRateForSection(sectIndex);
			double slip = aveSlips[sectIndex]/(double)sectCounts[sectIndex];
			func.set(slipRate, slip);
		}
		
		List<PlotElement> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		funcs.add(func);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.BLACK));
		
		String xLabel = "Slip Rate (m/yr)";
		String yLabel = "Average Slip (m)";
		String fName = "slips_vs_rate";
		String title = "RSQSim Ave Slips vs Slip Rate";
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		CommandLineInversionRunner.setFontSizes(gp);
		PlotSpec spec = new PlotSpec(funcs, chars, title, xLabel, yLabel);
		gp.setBackgroundColor(Color.WHITE);
		gp.drawGraphPanel(spec);
		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPNG(new File(plotsDir, fName+".png").getAbsolutePath());
		gp.saveAsPDF(new File(plotsDir, fName+".pdf").getAbsolutePath());
	}
	
	public void plotJumps(double jumpDist, double minMag) throws IOException {
		EvenlyDiscretizedFunc ucerfFunc = CommandLineInversionRunner.getJumpFuncs(
				ucerfSol, ucerfDistances, jumpDist, minMag, null)[0];
		EvenlyDiscretizedFunc simFunc = CommandLineInversionRunner.getJumpFuncs(
				simSol, simDistances, jumpDist, minMag, null)[0];
		ucerfFunc.scale(1d/ucerfFunc.calcSumOfY_Vals());
		simFunc.scale(1d/simFunc.calcSumOfY_Vals());
		EvenlyDiscretizedFunc empiricalFunc = new EvenlyDiscretizedFunc(ucerfFunc.getMinX(), ucerfFunc.getMaxX(), ucerfFunc.size());
		empiricalFunc.set(0, 0.405);
		empiricalFunc.set(1, 0.32);
		empiricalFunc.set(2, 0.1);
		empiricalFunc.set(3, 0.19);
		
		simFunc.setName("RSQSim");
		ucerfFunc.setName("UCERF");
		empiricalFunc.setName("Empirical");
		
		List<PlotElement> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		funcs.add(simFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		funcs.add(ucerfFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		funcs.add(empiricalFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GREEN));
		
		String xLabel = "Number of Jumps > "+(float)jumpDist+" km";
		String yLabel = "Fraction";
		String fName = "jump_rates_m"+(float)minMag;
		String title = "Multi Fault M"+(float)minMag+"+ Rupture Jumps > "+(float)jumpDist+" km";
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		CommandLineInversionRunner.setFontSizes(gp);
		PlotSpec spec = new PlotSpec(funcs, chars, title, xLabel, yLabel);
		gp.setBackgroundColor(Color.WHITE);
		gp.drawGraphPanel(spec);
		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPNG(new File(plotsDir, fName+".png").getAbsolutePath());
		gp.saveAsPDF(new File(plotsDir, fName+".pdf").getAbsolutePath());
	}
	
	public void plotJunctionSlipDists(double minDist, double maxDist, double minDiscrepFract,
			double maxDiscrepFract, boolean useMaxSlips, boolean plotJunctions) throws IOException {
		System.out.println("Plotting for dist: "+(float)minDist+"=>"+(float)maxDist+", discrep: "
			+(float)minDiscrepFract+"=>"+(float)maxDiscrepFract+", useMax="+useMaxSlips);
		// first find junctions of interest
		HashSet<IDPairing> ucerfPairs = new HashSet<IDPairing>(ucerfParentTogetherRates.keySet());
		HashSet<IDPairing> simPairs = new HashSet<IDPairing>(simParentTogetherRates.keySet());
		
		File aggregatedDir = new File(plotsDir, "junc_slip_dist_agg");
		if (!aggregatedDir.exists())
			aggregatedDir.mkdir();
		File junctDir = new File(plotsDir, "junc_slip_dist_junctions");
		if (!junctDir.exists())
			junctDir.mkdir();
		
		List<IDPairing> allPairings = Lists.newArrayList();
		// we only want parents in both solutions
		for (IDPairing ucerfPair : ucerfPairs)
			if (simPairs.contains(ucerfPair))
				allPairings.add(ucerfPair);
		
		InversionFaultSystemRupSet ucerfRupSet = ucerfSol.getRupSet();
		
		int minNumSectsPerParent = 4;
		
		int numPairingsUsed = 0;
		int numUCERFRupsUsed = 0;
		int numSimRupsUsed = 0;
		
		List<double[]> simSlipAlongs;
		if (useMaxSlips)
			simSlipAlongs = simSlipsMax;
		else
			simSlipAlongs = simSlipsMean;
		List<double[]> ucerfSlipAlongs = ucerfRupSet.getSlipOnSectionsForAllRups();
		
		double[] totalSimSlips = new double[minNumSectsPerParent*2];
		double[] totalSimSlipRates = new double[minNumSectsPerParent*2];
		double[] totalUCERFSlips = new double[minNumSectsPerParent*2];
		
		for (IDPairing pair : allPairings) {
//			System.out.println("slipLeft="+slipLeft+", slipRight="+slipRight+", discrepFract="+discrepFract);
			
			HashSet<Integer> simLeftParents = new HashSet<Integer>();
			HashSet<Integer> simRightParents = new HashSet<Integer>();
			for (Integer simParent : simToUCERF_parentsMap.keySet()) {
				int ucerfParent = simToUCERF_parentsMap.get(simParent);
				if (ucerfParent == pair.getID1())
					simLeftParents.add(simParent);
				if (ucerfParent == pair.getID2())
					simRightParents.add(simParent);
			}
			HashSet<Integer> ucerfLeftParents = new HashSet<Integer>();
			HashSet<Integer> ucerfRightParents = new HashSet<Integer>();
			ucerfLeftParents.add(pair.getID1());
			ucerfRightParents.add(pair.getID2());
			
			// ucerf slips
			SlipsAlongResult ucerfSlipsResult = getSlipsAlongAtJunction(ucerfSol, ucerfSlipAlongs,
					minNumSectsPerParent, ucerfLeftParents, ucerfRightParents,
					ucerfDistances, minDist, maxDist);
			List<double[]> ucerfSlipsList = ucerfSlipsResult.slips;
			if (ucerfSlipsList.isEmpty())
				continue;
			// sim slips
			SlipsAlongResult simSlipsResult = getSlipsAlongAtJunction(simSol, simSlipAlongs,
					minNumSectsPerParent, simLeftParents, simRightParents,
					simDistances, minDist, maxDist);
			List<double[]> simSlipsList = simSlipsResult.slips;
			List<double[]> simSlipRatesList = simSlipsResult.slipRates;
			if (simSlipsList.isEmpty())
				continue;
			
			double[] simSlips = averageSlips(simSlipsList, simSlipsResult.rates);
			double[] ucerfSlips = averageSlips(ucerfSlipsList, ucerfSlipsResult.rates);
			double[] simSlipRates = averageSlips(simSlipRatesList, simSlipsResult.rates);
			
			// now calculate the slip discrepancy
			double slipLeft = simSlipRates[simSlipRates.length/2-1];
			double slipRight = simSlipRates[simSlipRates.length/2];
			Preconditions.checkState(!Double.isNaN(slipLeft) && !Double.isNaN(slipRight));
			
			// make sure left is higher
			if (slipLeft < slipRight) {
				pair = pair.getReversed();
				double slipTemp = slipLeft;
				slipLeft = slipRight;
				slipRight = slipTemp;
				double[] newSimSlips = new double[simSlips.length];
				double[] newUCERFSlips = new double[simSlips.length];
				double[] newSimSlipRates = new double[simSlips.length];
				for (int i=0; i<simSlips.length; i++) {
					int j = simSlips.length-i-1;
					newSimSlips[i] = simSlips[j];
					newUCERFSlips[i] = ucerfSlips[j];
					newSimSlipRates[i] = simSlipRates[j];
				}
				simSlips = newSimSlips;
				simSlipRates = newSimSlipRates;
				ucerfSlips = newUCERFSlips;
			}
			double discrep = slipLeft - slipRight;
			double discrepFract = discrep / slipLeft;
			
			if (discrepFract < minDiscrepFract || discrepFract > maxDiscrepFract)
				continue;
			
			// we have a match!
			numPairingsUsed++;
			numSimRupsUsed += simSlipsList.size();
			numUCERFRupsUsed += ucerfSlipsList.size();
			
			// divide by 2 here because of bundling in slip rates
			for (int i=0; i<simSlips.length; i++)
				totalSimSlips[i] += simSlips[i];
			for (int i=0; i<simSlipRates.length; i++)
				totalSimSlipRates[i] += simSlipRates[i];
			for (int i=0; i<ucerfSlips.length; i++)
				totalUCERFSlips[i] += ucerfSlips[i];
			
			if (plotJunctions) {
				String nameLeft = null;
				String nameRight = null;
				for (FaultSectionPrefData sect : ucerfRupSet.getFaultSectionDataList()) {
					if (sect.getParentSectionId() == pair.getID1())
						nameLeft = sect.getParentSectionName();
					if (sect.getParentSectionId() == pair.getID2())
						nameRight = sect.getParentSectionName();
					if (nameLeft != null && nameRight != null)
						break;
				}
				String fName = nameLeft.replaceAll(" ", "_")+"_to_"+nameRight.replaceAll(" ", "_");
				if (useMaxSlips)
					fName += "_maxslips";
				String title = "Slip Distribution for "+nameLeft+" to "+nameRight;
				
//				double[] simSlipRates = Arrays.copyOfRange(simSlips, simSlips.length/2, simSlips.length);
				
				// normalize slip rates
//				double maxSlipRate = StatUtils.max(simSlipRates);
//				for (int i=0; i<simSlipRates.length; i++)
//					simSlipRates[i] /= maxSlipRate;
//				Preconditions.checkState((float)StatUtils.max(simSlipRates) == 1f);
				
				writeJunctionSlipPlot(simSlips, ucerfSlips, minNumSectsPerParent, numPairingsUsed,
						numSimRupsUsed, numUCERFRupsUsed, simSlipRates, nameLeft, nameRight, fName, title, junctDir);
			}
		}
		
		if (plotJunctions)
			// this was a junction plot
			return;
		
		System.out.println("Found "+numPairingsUsed+" junctions with "+numSimRupsUsed
				+" sim rups and "+numUCERFRupsUsed+" UCERF rups");
		if (numPairingsUsed == 0) {
			System.out.println("Skipping");
			return;
		}
		double maxSim = StatUtils.max(totalSimSlips);
		double maxSimRate = StatUtils.max(totalSimSlipRates);
		double maxUCERF = StatUtils.max(totalUCERFSlips);
		for (int i=0; i<totalSimSlips.length; i++) {
			totalSimSlips[i] /= maxSim;
			totalSimSlipRates[i] /= maxSimRate;
			totalUCERFSlips[i] /= maxUCERF;
		}
		
		String fName = "junc_slip_dist_"+minDist+"_to_"+maxDist+"km_fracts"
				+(float)minDiscrepFract+"_"+(float)maxDiscrepFract;
		if (useMaxSlips)
			fName += "_maxslips";
		String title = "Slip Distribution at Junctions Dists "+(float)minDist+"=>"+(float)maxDist
				+"km Slip Discrep"+(float)minDiscrepFract+"=>"+(float)maxDiscrepFract;
		
		writeJunctionSlipPlot(totalSimSlips, totalUCERFSlips, minNumSectsPerParent, numPairingsUsed,
				numSimRupsUsed, numUCERFRupsUsed, totalSimSlipRates, null, null, fName, title, aggregatedDir);
	}
	
	private static void writeJunctionSlipPlot(double[] simSlips, double[] ucerfSlips, int minNumSectsPerParent,
			int numPairingsUsed, int numSimRupsUsed, int numUCERFRupsUsed, double[] simSlipRates,
			String leftText, String rightText, String fName, String title, File plotsDir) throws IOException {
		DefaultXY_DataSet simFunc = new DefaultXY_DataSet();
		DefaultXY_DataSet ucerfFunc = new DefaultXY_DataSet();
		DefaultXY_DataSet slipRateFunc = new DefaultXY_DataSet();
		for (int i=0; i<simSlips.length; i++) {
			double startX = i-minNumSectsPerParent;
			double endX = startX + 1;
			
			simFunc.set(startX, simSlips[i]);
			simFunc.set(endX, simSlips[i]);
			
			ucerfFunc.set(startX, ucerfSlips[i]);
			ucerfFunc.set(endX, ucerfSlips[i]);
			if (simSlipRates != null) {
				slipRateFunc.set(startX, simSlipRates[i]);
				slipRateFunc.set(endX, simSlipRates[i]);
			}
		}
		DefaultXY_DataSet divider = new DefaultXY_DataSet();
		divider.set(0d, 0d);
		divider.set(0d, 1d);
		
		simFunc.setName("RSQSim");
		ucerfFunc.setName("UCERF");
		
		List<PlotElement> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		if (simSlipRates != null) {
			funcs.add(slipRateFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GREEN));
		}
		funcs.add(simFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		funcs.add(ucerfFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		funcs.add(divider);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		
		String xLabel = "Distance";
		String yLabel = "Normalized Slip On Section";
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		CommandLineInversionRunner.setFontSizes(gp);
		PlotSpec spec = new PlotSpec(funcs, chars, title, xLabel, yLabel);
		List<XYTextAnnotation> annotations = Lists.newArrayList();
		Font font = new Font(Font.SERIF, Font.PLAIN, 24);
		if (leftText != null) {
			XYTextAnnotation leftAnn = new XYTextAnnotation(leftText, -4, 0.3);
			leftAnn.setFont(font);
			leftAnn.setPaint(Color.BLACK);
			leftAnn.setTextAnchor(TextAnchor.CENTER_LEFT);
			annotations.add(leftAnn);
			XYTextAnnotation rightAnn = new XYTextAnnotation(rightText, 4, 0.3);
			rightAnn.setFont(font);
			rightAnn.setPaint(Color.BLACK);
			rightAnn.setTextAnchor(TextAnchor.CENTER_RIGHT);
			annotations.add(rightAnn);
		} else {
			XYTextAnnotation junctsAnn = new XYTextAnnotation(numPairingsUsed+" Junctions", -4, 0.3);
			junctsAnn.setFont(font);
			junctsAnn.setPaint(Color.BLACK);
			junctsAnn.setTextAnchor(TextAnchor.CENTER_LEFT);
			annotations.add(junctsAnn);
		}
		XYTextAnnotation simAnn = new XYTextAnnotation(numSimRupsUsed+" RSQSim Rups", -4, 0.2);
		simAnn.setFont(font);
		simAnn.setPaint(Color.RED);
		simAnn.setTextAnchor(TextAnchor.CENTER_LEFT);
		annotations.add(simAnn);
		XYTextAnnotation ucerfAnn = new XYTextAnnotation(numUCERFRupsUsed+" UCERF Rups", -4, 0.1);
		ucerfAnn.setFont(font);
		ucerfAnn.setPaint(Color.BLUE);
		ucerfAnn.setTextAnchor(TextAnchor.CENTER_LEFT);
		annotations.add(ucerfAnn);
		spec.setPlotAnnotations(annotations);
		gp.setBackgroundColor(Color.WHITE);
		gp.drawGraphPanel(spec);
		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPNG(new File(plotsDir, fName+".png").getAbsolutePath());
		gp.saveAsPDF(new File(plotsDir, fName+".pdf").getAbsolutePath());
	}
	
	private double[] averageSlips(List<double[]> slipsList, List<Double> rates) {
		double totWeight = 0d;
		for (double rate : rates)
			totWeight += rate;
		// this will average the NORMALIZED slips
		double[] slips = new double[slipsList.get(0).length];
		for (int j = 0; j < slipsList.size(); j++) {
			double weight = rates.get(j)/totWeight;
			double[] subSlips = slipsList.get(j);
			// just the first half
			double max = StatUtils.max(subSlips);
			for (int i=0; i<subSlips.length; i++) {
				slips[i] += weight*(subSlips[i]/max);
			}
		}
		return slips;
	}
	
	private static class SlipsAlongResult {
		List<double[]> slips;
		List<double[]> slipRates;
		List<Double> rates;
		
		public SlipsAlongResult(List<double[]> slips, List<double[]> slipRates,
				List<Double> rates) {
			super();
			this.slips = slips;
			this.slipRates = slipRates;
			this.rates = rates;
		}
	}
	private SlipsAlongResult getSlipsAlongAtJunction(FaultSystemSolution sol, List<double[]> slipAlongs,
			int minNumSectsPerParent, HashSet<Integer> leftParents, HashSet<Integer> rightParents,
			Map<IDPairing, Double> distsMap, double minDist, double maxDist) {
		FaultSystemRupSet rupSet = sol.getRupSet();
		HashSet<Integer> rupIndexes = new HashSet<Integer>();
		for (int id1 : leftParents) {
			for (int rupIndex : rupSet.getRupturesForParentSection(id1)) {
				List<Integer> parents = rupSet.getParentSectionsForRup(rupIndex);
				for (int id2 : rightParents) {
					if (parents.contains(id2)) {
						rupIndexes.add(rupIndex);
						break;
					}
				}
			}
		}
		
		List<double[]> slipVals = Lists.newArrayList();
		List<double[]> slipRates = Lists.newArrayList();
		List<Double> rates = Lists.newArrayList();
		
		for (int rupIndex : rupIndexes) {
			List<FaultSectionPrefData> sects = rupSet.getFaultSectionDataForRupture(rupIndex);
			int junctIndex = -1;
			boolean reversed = false;
			for (int i=1; i<sects.size(); i++) {
				int parentLeft = sects.get(i-1).getParentSectionId();
				int parentRight = sects.get(i).getParentSectionId();
				if (leftParents.contains(parentLeft) && rightParents.contains(parentRight)) {
					junctIndex = i;
					reversed = false;
					break;
				} else if (leftParents.contains(parentRight) && rightParents.contains(parentLeft)) {
					junctIndex = i;
					reversed = true;
					break;
				}
			}
			if (junctIndex < 0)
				continue;
			// now check distances
			double dist = distsMap.get(new IDPairing(sects.get(junctIndex-1).getSectionId(),
					sects.get(junctIndex).getSectionId()));
			if (dist < minDist || dist > maxDist)
				continue;
			int buffStartIndex = junctIndex - minNumSectsPerParent + 1;
			int buffEndIndex = buffStartIndex + 2*minNumSectsPerParent - 1;
			// make sure there's a big enough buffer
			if (buffStartIndex < 0 || buffEndIndex >= sects.size())
				continue;
			double[] rupSlips = slipAlongs.get(rupIndex);
			double[] sectSlipRates = new double[rupSlips.length];
			for (int i=0; i<sectSlipRates.length; i++)
				sectSlipRates[i] = sects.get(i).getOrigAveSlipRate();
			double[] slips = new double[2*minNumSectsPerParent];
			double[] slipRatesVals = new double[2*minNumSectsPerParent];
			int cnt = 0;
			if (reversed)
				for (int i=buffEndIndex; i>=buffStartIndex; i--)
					slips[cnt++] = rupSlips[i];
			else
				for (int i=buffStartIndex; i<=buffEndIndex; i++)
					slips[cnt++] = rupSlips[i];
			// now slip rates
			cnt = 0;
			if (reversed)
				for (int i=buffEndIndex; i>=buffStartIndex; i--)
					slipRatesVals[cnt++] = sectSlipRates[i];
			else
				for (int i=buffStartIndex; i<=buffEndIndex; i++)
					slipRatesVals[cnt++] = sectSlipRates[i];
			Preconditions.checkState(cnt == slips.length);
			slipVals.add(slips);
			slipRates.add(slipRatesVals);
			rates.add(sol.getRateForRup(rupIndex));
//			double slipLeft = sects.get(junctIndex-1).getOrigAveSlipRate();
//			double slipRight = sects.get(junctIndex).getOrigAveSlipRate();
//			if (reversed) {
//				double slipTemp = slipLeft;
//				slipLeft = slipRight;
//				slipRight = slipTemp;
//			}
//			if (slipLeft < slipRight) {
//				System.out.println("VIOLATION!!!!! slipLeft="+slipLeft+", slipRight="+slipRight+", reversed="+reversed);
//				if (!(rupSet instanceof InversionFaultSystemRupSet)) {
//					// this means it's a sim one
//					int u3parentLeft = simToUCERF_parentsMap.get(sects.get(junctIndex-1).getParentSectionId());
//					int u3parentRight = simToUCERF_parentsMap.get(sects.get(junctIndex).getParentSectionId());
//					if (reversed) {
//						int temp = u3parentLeft;
//						u3parentLeft = u3parentRight;
//						u3parentRight = temp;
//					}
//					System.out.println("U3 slips: slipLeft="+ucerfParentSlipsMap.get(u3parentLeft)
//							+", slipRight="+ucerfParentSlipsMap.get(u3parentRight));
//				}
//			}
		}
		return new SlipsAlongResult(slipVals, slipRates, rates);
	}
	
	public void plotJumpDistHists() throws IOException {
		HistogramFunction simFunc = getJumpDistHist(simSol, simDistances);
		HistogramFunction ucerfFunc = getJumpDistHist(ucerfSol, ucerfDistances);
		
		simFunc.setName("RSQSim");
		ucerfFunc.setName("UCERF");
		
		List<PlotElement> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		funcs.add(simFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		funcs.add(ucerfFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		
		String xLabel = "Max Rupture Jump Distance (km)";
		String yLabel = "Rate";
		String fName = "jump_dist_rates";
		String title = "Rupture Jump Distance Distributions";
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		CommandLineInversionRunner.setFontSizes(gp);
		PlotSpec spec = new PlotSpec(funcs, chars, title, xLabel, yLabel);
		gp.setBackgroundColor(Color.WHITE);
		gp.setUserBounds(new Range(0d, 10d), new Range(1e-6, 1e1));
		gp.drawGraphPanel(spec, false, true);
		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPNG(new File(plotsDir, fName+".png").getAbsolutePath());
		gp.saveAsPDF(new File(plotsDir, fName+".pdf").getAbsolutePath());
	}
	
	private static HistogramFunction getJumpDistHist(FaultSystemSolution sol, Map<IDPairing, Double> distances) {
		double min = 0.1;
		double max = 9.9;
		double delta = 0.2;
		int num = (int)((max - min) / delta) + 2;
		HistogramFunction func = new HistogramFunction(min, num, delta);
		
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		for (int rupIndex=0; rupIndex<rupSet.getNumRuptures(); rupIndex++) {
			double maxJump = 0;
			List<Integer> sects = rupSet.getSectionsIndicesForRup(rupIndex);
			for (int i=0; i<sects.size()-1; i++) {
				// first the next section
				// jump from i to i+1
				double minDist = distances.get(new IDPairing(sects.get(i), sects.get(i+1)));
				// now see if there's a non neighboring section that's closer (important for weird RSQSim rups
				for (int j=0; j<sects.size(); j++) {
					// don't check against neighbors or the section itself
					if (j == i-1 || j == i || j == i+1)
						continue;
					IDPairing pair = new IDPairing(sects.get(i), sects.get(j));
					Double dist = distances.get(pair);
					if (dist == null)
						continue;
					if (dist < minDist)
						minDist = dist;
				}
				if (minDist > maxJump)
					maxJump = minDist;
			}
			if (maxJump > max + 0.5*delta)
				continue;
//			System.out.println("MAX="+func.getMaxX());
//			System.out.println("Dist: "+maxJump);
			func.add(maxJump, sol.getRateForRup(rupIndex));
		}
		
		return func;
	}

	public static void main(String[] args) throws IOException, DocumentException {
		File invSolDir = new File(UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, "InversionSolutions");
		System.out.println("Loading UCERF sol");
		InversionFaultSystemSolution ucerfSol = FaultSystemIO.loadInvSol(new File(invSolDir,
						"2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		CoulombRates coulombRates = CoulombRates.loadUCERF3CoulombRates(FaultModels.FM3_1);
		
		File simMainDir = new File("/home/kevin/Simulators");
		File resultsDir = new File(simMainDir, "multiFault");
		if (!resultsDir.exists())
			resultsDir.mkdir();
		
		File simSolFile = new File(simMainDir, "rsqsim_long_sol.zip");
		File simSlipMaxFile = new File(resultsDir, "rsqsim_slips_max.bin");
		File simSlipMeanFile = new File(resultsDir, "rsqsim_slips_mean.bin");
		
		List<double[]> simSlipsMean, simSlipsMax;
		
		FaultSystemSolution simSol;
		if (simSolFile.exists() && simSlipMaxFile.exists() && simSlipMeanFile.exists()) {
			System.out.println("Loading Sim sol");
			simSol = FaultSystemIO.loadSol(simSolFile);
			simSlipsMean = MatrixIO.doubleArraysListFromFile(simSlipMeanFile);
			simSlipsMax = MatrixIO.doubleArraysListFromFile(simSlipMaxFile);
		} else {
			File geomFile = new File(simMainDir, "ALLCAL2_1-7-11_Geometry.dat");
			System.out.println("Loading geometry from "+geomFile.getAbsolutePath());
			General_EQSIM_Tools tools = new General_EQSIM_Tools(geomFile);
//			File eventFile = new File(dir, "eqs.ALLCAL2_RSQSim_sigma0.5-5_b=0.015.barall");
			File eventFile = new File(simMainDir, "eqs.ALLCAL2_RSQSim_sigma0.5-5_b=0.015.long.barall");
			System.out.println("Loading events...");
			List<? extends SimulatorEvent> events = EQSIMv06FileReader.readEventsFile(eventFile, tools.getElementsList());
			
			SimulatorFaultSystemSolution theSol = SimulatorFaultSystemSolution.build(tools.getElementsList(), events,
					tools.getSimulationDurationYears());
			simSlipsMean = theSol.getSlipAlongRupVals(false);
			simSlipsMax = theSol.getSlipAlongRupVals(true);
			MatrixIO.doubleArraysListToFile(simSlipsMean, simSlipMeanFile);
			MatrixIO.doubleArraysListToFile(simSlipsMax, simSlipMaxFile);
			simSol = theSol;
			FaultSystemIO.writeSol(simSol, simSolFile);
		}
		
		SimJunctionMapper mapper = new SimJunctionMapper(ucerfSol, simSol,
				simSlipsMean, simSlipsMax, coulombRates, resultsDir);

		mapper.plotMultiRateVsCoulomb(mapper.ucerfParentTogetherRates, "UCERF3", false);
		mapper.plotMultiRateVsCoulomb(mapper.normUCERFParentTogetherRates, "UCERF3", true);
		mapper.plotMultiRateVsCoulomb(mapper.simParentTogetherRates, "RSQSim", false);
		mapper.plotMultiRateVsCoulomb(mapper.normSimParentTogetherRates, "RSQSim", true);
		
		mapper.plotSimVsUCERFParentRates(false);
		mapper.plotSimVsUCERFParentRates(true);
		
		mapper.plotJumps(0.1d, 6d);
		mapper.plotJumps(1d, 6d);
		
		mapper.plotJumps(0.1d, 7d);
		mapper.plotJumps(1d, 7d);
		
		mapper.plotSlipVsRate();
		
		mapper.plotJumpDistHists();
		
		List<double[]> distRanges = Lists.newArrayList(new double[] {0d, 0.1},
				new double[] {0.1d, 1d}, new double[] {1d, 10d});
		List<double[]> discrepRanges = Lists.newArrayList(new double[] {0d, 0.2},
				new double[] {0.2d, 1d});
//		boolean[] maxes = { false, true };
		boolean[] maxes = { true };
		
		for (boolean useMaxSlips : maxes) {
			for (double[] distRange : distRanges) {
				for (double[] discrepRange : discrepRanges) {
					mapper.plotJunctionSlipDists(distRange[0], distRange[1],
							discrepRange[0], discrepRange[1], useMaxSlips, false);
				}
			}
			mapper.plotJunctionSlipDists(0, 100d, 0d, 1d, useMaxSlips, true);
		}
	}

}
