package scratch.kevin.simulators;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.dom4j.DocumentException;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityOptions;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.iden.ElementIden;
import org.opensha.sha.simulators.iden.LogicalAndRupIden;
import org.opensha.sha.simulators.iden.LogicalOrRupIden;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.iden.SupraSeisRupIden;
import org.opensha.sha.simulators.parsers.EQSIMv06FileReader;
import org.opensha.sha.simulators.utils.General_EQSIM_Tools;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.analysis.FaultSysSolutionERF_Calc;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.LastEventData;
import scratch.UCERF3.utils.paleoRateConstraints.PaleoRateConstraint;
import scratch.UCERF3.utils.paleoRateConstraints.UCERF3_PaleoRateConstraintFetcher;
import scratch.kevin.simulators.catBuild.RandomCatalogBuilder;
import scratch.kevin.simulators.dists.RandomDistType;

public class UC3PaleoOpenCalc {
	
	private static final double calcEndYear = 2014;
	
	private static class PaleoOpenIden extends ElementIden {
		
		private PaleoRateConstraint constr;
		private long dateOfLastEventMillis;
		
		public PaleoOpenIden(PaleoRateConstraint constr, long dateOfLastEventMillis, List<SimulatorElement> geom) {
			super(constr.getPaleoSiteName(), getForLoc(constr.getPaleoSiteName(), constr.getPaleoSiteLoction(), geom));
			this.constr = constr;
			this.dateOfLastEventMillis = dateOfLastEventMillis;
		}
		
		public PaleoRateConstraint getPaleoConstr() {
			return constr;
		}
		
		public long getDateLastEvent() {
			return dateOfLastEventMillis;
		}
		
		public double getOpenInterval() {
			return LastEventData.OPEN_INTERVAL_BASIS_YEAR - millisToYear(dateOfLastEventMillis);
		}
		
	}
	
	private static double millisToYear(long millis) {
		double year = ((double)millis / ProbabilityModelsCalc.MILLISEC_PER_YEAR) + 1970d;
		return year;
	}
	
	private static long yearToMillis(double year) {
		return Math.round((year-1970.0)*ProbabilityModelsCalc.MILLISEC_PER_YEAR);
	}
	
	private static int getForLoc(String name, Location loc, List<SimulatorElement> geom) {
		SimulatorElement closest = null;
		double closestDist = Double.POSITIVE_INFINITY;
		
		for (SimulatorElement el : geom) {
			Location centerLocation = el.getCenterLocation();
			double hDist = LocationUtils.horzDistanceFast(loc, centerLocation);
			if (hDist < closestDist) {
				closestDist = hDist;
				closest = el;
			}
		}
		
		System.out.println("Mapped "+name+" to "+closest.getName()+" (hDist="+(float)closestDist+" km)"
				+"\t(depth="+(float)closest.getCenterLocation().getDepth()+" km)");
		return closest.getID();
	}
	
	private static List<PaleoOpenIden> getUC3SitesWithOpenIntervals(
			FaultSystemSolution sol, List<SimulatorElement> geom) throws IOException {
		List<FaultSectionPrefData> fsd = sol.getRupSet().getFaultSectionDataList();
		ArrayList<PaleoRateConstraint> constraints = UCERF3_PaleoRateConstraintFetcher.getConstraints(fsd);
		
		List<PaleoOpenIden> idens = Lists.newArrayList();
		
		for (PaleoRateConstraint constr : constraints) {
			if (constr.getPaleoSiteName().toLowerCase().startsWith("compton"))
				// skip compton, not in RSQSim
				continue;
			long dateLastEvent = fsd.get(constr.getSectionIndex()).getDateOfLastEvent();
			if (dateLastEvent > Long.MIN_VALUE)
				idens.add(new PaleoOpenIden(constr, dateLastEvent, geom));
		}
		System.out.println(idens.size()+"/"+constraints.size()+" constraints had open intervals");
		
		return idens;
	}
	
	private static EvenlyDiscretizedFunc createFunc(double startYear, double endYear, double deltaYears) {
		Preconditions.checkArgument(endYear > startYear);
		Preconditions.checkState((endYear - startYear) % deltaYears == 0d);
		double numDouble = (endYear - startYear)/deltaYears + 1d;
		Preconditions.checkState(numDouble == (double)((int)numDouble));
		return new EvenlyDiscretizedFunc(startYear, (int)numDouble, deltaYears);
	}
	
	private static EvenlyDiscretizedFunc calcUCERF3(double startYear, double endYear, double deltaYears,
			FaultSystemSolution sol, List<PaleoOpenIden> idens, ProbabilityModelOptions probModel,
			MagDependentAperiodicityOptions cov) throws IOException {
		EvenlyDiscretizedFunc func = createFunc(startYear, endYear, deltaYears);
		
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		// reset last event data
		Map<Integer, List<LastEventData>> datas = LastEventData.load();
		LastEventData.populateSubSects(rupSet.getFaultSectionDataList(), datas);
		
		boolean constrainAtPaleoSiteOnly = true;
		
		String name = "UCERF3";
		if (cov == null)
			name += " ("+probModel.name()+")";
		else
			name += " ("+cov.name()+")";
		
		for (int i=func.size(); --i>=0;) {
			double year = func.getX(i);
			double duration = calcEndYear - year;
			
			long startYearInMillis = yearToMillis(year);
			
			// clear any last open intervals that are too short
			List<Integer> constrIndexes = Lists.newArrayList();
			for (FaultSectionPrefData sect : rupSet.getFaultSectionDataList()) {
				if (sect.getDateOfLastEvent() > startYearInMillis) {
//					System.out.println("Clearing date last event on "+sect.getParentSectionName());
					sect.setDateOfLastEvent(Long.MIN_VALUE);
				} else if (sect.getDateOfLastEvent() > Long.MIN_VALUE && !constrainAtPaleoSiteOnly) {
					constrIndexes.add(sect.getSectionId());
				}
			}
			if (constrainAtPaleoSiteOnly) {
				// only at paleo sites
				for (PaleoOpenIden iden : idens) {
					if (millisToYear(iden.getDateLastEvent()) < year)
						constrIndexes.add(iden.getPaleoConstr().getSectionIndex());
				}
			}
			
			System.out.println("Testing "+name+" from "+year
					+" to "+calcEndYear+" with "+constrIndexes.size()+" sites");
			
			FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
			erf.setParameter(ProbabilityModelParam.NAME, probModel);
			if (cov != null)
				erf.setParameter(MagDependentAperiodicityParam.NAME, cov);
			if (probModel != ProbabilityModelOptions.POISSON)
				erf.getTimeSpan().setStartTime((int)year);
			erf.getTimeSpan().setDuration(duration);
			if (probModel != ProbabilityModelOptions.POISSON)
				Preconditions.checkState(erf.getTimeSpan().getStartTimeYear() == (int)year);
			erf.updateForecast();
			
			double sumRate = 0d;
			Map<Integer, List<Integer>> sectsForRups = Maps.newHashMap();
			HashSet<Integer> rups = new HashSet<Integer>();
			int numDuplicates = 0;
			int total = 0;
			for (int sectIndex : constrIndexes) {
				List<Integer> myRups = rupSet.getRupturesForSection(sectIndex);
				for (int rup : myRups) {
					List<Integer> sects = sectsForRups.get(rup);
					if (sects == null) {
						sects = Lists.newArrayList();
						sectsForRups.put(rup, sects);
					}
					sects.add(sectIndex);
					total++;
					if (rups.contains(rup))
						numDuplicates++;
					rups.add(rup);
				}
			}
//			System.out.println(numDuplicates+"/"+total+" duplicates!");
			List<Double> rupProbs = Lists.newArrayList();
			for (int rup : rups) {
				int sourceID = erf.getSrcIndexForFltSysRup(rup);
				if (sourceID < 0) {
//					Preconditions.checkState(fss.getRateForRup(rup) == 0);
					continue;
				}
				double rate = erf.getSource(sourceID).computeTotalEquivMeanAnnualRate(duration)*duration;
				rupProbs.add(erf.getSource(sourceID).computeTotalProb());
				sumRate += rate;
			}
			double probNot = (Math.exp(-sumRate));
			System.out.println(name+" "+year+" survival prob: "+probNot);
			
			func.set(i, probNot);
		}
		
		func.setName(name);
		
		return func;
	}
	
	private static EvenlyDiscretizedFunc calcSimulator(double startYear, double endYear, double deltaYears,
			List<? extends SimulatorEvent> events, List<PaleoOpenIden> idens, RandomDistType randType) throws IOException {
		EvenlyDiscretizedFunc func = createFunc(startYear, endYear, deltaYears);
		
		if (randType != null)
			events = RandomCatalogBuilder.getRandomResampledCatalog(events, idens, randType, false);
		
		List<double[]> eventTimes = calcIdenEventTimes(events, idens);
		
		String name = "Simulator";
		if (randType != null)
			name += " ("+randType.getName()+")";
		
		for (int i=func.size(); --i>=0;) {
			double year = func.getX(i);
			double duration = calcEndYear - year;
			
			List<Double> ois = Lists.newArrayList();
			List<double[]> includedEventTimes = Lists.newArrayList();
			for (int j=0; j<idens.size(); j++) {
				PaleoOpenIden iden = idens.get(j);
				if (millisToYear(iden.getDateLastEvent()) < year) {
					includedEventTimes.add(eventTimes.get(j));
					ois.add(duration);
				}
			}
			
			System.out.println("Testing "+name+" from "+year
					+" to "+calcEndYear+" with "+ois.size()+" sites");
			
			double fract = calcFractTimeInWindow(events.get(0).getTimeInYears(),
					events.get(events.size()-1).getTimeInYears(), includedEventTimes, ois);
			
			System.out.println(name+" survival prob: "+fract);
			func.set(i, fract);
		}
		
		func.setName(name);
		
		return func;
	}
	
	private static List<double[]> calcIdenEventTimes(List<? extends SimulatorEvent> events, List<? extends RuptureIdentifier> idens) {
		List<double[]> idenEventTimes = Lists.newArrayList();
		for (int i=0; i<idens.size(); i++) {
			List<? extends SimulatorEvent> matches = idens.get(i).getMatches(events);
			double[] times = new double[matches.size()];
			for (int j=0; j<matches.size(); j++)
				times[j] = matches.get(j).getTimeInYears();
			idenEventTimes.add(times);
		}
		return idenEventTimes;
	}
	
	private static double calcFractTimeInWindow(double catStartTime, double catEndTime,
			List<double[]> idenEventTimes, List<Double> idenOIs) {
		int timesInWindows = 0;
		int timeSteps = 0;
		
		timeLoop:
		for (double time=catStartTime+100d; time<catEndTime; time += 1d) {
//			if (loopCnt % 10000 == 0)
//				System.out.println("Loop "+loopCnt+"/"+expected);
//			loopCnt++;
			
			boolean match = true;
			for (int i=0; i<idenEventTimes.size(); i++) {
				double[] times = idenEventTimes.get(i);
				int bin = Arrays.binarySearch(times, time);
				double eventTime;
				if (bin >= 0) {
					eventTime = times[bin];
				} else {
					int insertionPoint = -(bin + 1);
					if (insertionPoint == 0)
						// this means we're before the first event
						continue timeLoop;
					eventTime = times[insertionPoint-1];
				}
				double oi = time - eventTime;
//				if (loopCnt % 10000 == 0) 
//					System.out.println("Debug! time="+time+"\teTime="+eventTime
//							+"\toi="+oi+"\tdataOI="+openIntervals.get(i));
				if (oi < idenOIs.get(i)) {
					match = false;
					break;
				}
			}
			timeSteps++;
			
			if (match)
				timesInWindows++;
		}
		
		double fract = (double)timesInWindows/(double)timeSteps;
		
		System.out.println(timesInWindows+"/"+timeSteps+" ("+(float)(fract*100)+" %) match data");
		
		return fract;
	}
	
	private static void plot(EvenlyDiscretizedFunc... funcs) {
		List<DiscretizedFunc> elems = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		List<Color> colors = GraphWindow.generateDefaultColors();
		
		double minNonZero = 1d;
		
		for (int i=0; i<funcs.length; i++) {
			EvenlyDiscretizedFunc func = funcs[i];
			for (Point2D pt : func)
				if (pt.getY() < minNonZero && pt.getY() > 0)
					minNonZero = pt.getY();
			elems.add(func);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, colors.get(i)));
		}
		
		PlotSpec spec = new PlotSpec(elems, chars, "Survival Functions", "Year", "Survival Probability");
		spec.setLegendVisible(true);
		
		GraphWindow gw = new GraphWindow(spec, null, false, false, null, new Range(minNonZero*0.5, 1d));
		gw.setYLog(true);
		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
	}

	public static void main(String[] args) throws IOException, DocumentException {
		File simDir = new File("/home/kevin/Simulators");
		File geomFile = new File(simDir, "ALLCAL2_1-7-11_Geometry.dat");
		System.out.println("Loading geometry...");
		General_EQSIM_Tools tools = new General_EQSIM_Tools(geomFile);
		List<SimulatorElement> geoms = tools.getElementsList();
		
		FaultSystemSolution sol = FaultSystemIO.loadSol(new File("/home/kevin/workspace/OpenSHA/dev/"
				+ "scratch/UCERF3/data/scratch/InversionSolutions/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		List<PaleoOpenIden> idens = getUC3SitesWithOpenIntervals(sol, geoms);
		
		Collections.sort(idens, new Comparator<PaleoOpenIden>() {

			@Override
			public int compare(PaleoOpenIden o1, PaleoOpenIden o2) {
				return new Long(o2.getDateLastEvent()).compareTo(o1.getDateLastEvent());
			}
		});
		
		System.out.println("Paleoseismic Sites:");
		for (PaleoOpenIden iden : idens) {
			System.out.println("\t"+iden.getPaleoConstr().getPaleoSiteName().trim()
					+", Rupture Year: "+(int)(millisToYear(iden.getDateLastEvent())+0.5d));
		}
		
		double startYear = 1800d;
		double endYear = 2010d;
		double deltaYears = 10d;
//		double deltaYears = 30d;

		EvenlyDiscretizedFunc uc3Poisson = calcUCERF3(
				startYear, endYear, deltaYears, sol, idens, ProbabilityModelOptions.POISSON, null);
		EvenlyDiscretizedFunc uc3PrefBlend = calcUCERF3(
				startYear, endYear, deltaYears, sol, idens, ProbabilityModelOptions.U3_PREF_BLEND, null);
		EvenlyDiscretizedFunc uc3High = calcUCERF3(
				startYear, endYear, deltaYears, sol, idens, ProbabilityModelOptions.U3_BPT,
				MagDependentAperiodicityOptions.HIGH_VALUES);
		EvenlyDiscretizedFunc uc3Mid = calcUCERF3(
				startYear, endYear, deltaYears, sol, idens, ProbabilityModelOptions.U3_BPT,
				MagDependentAperiodicityOptions.MID_VALUES);
		EvenlyDiscretizedFunc uc3Low = calcUCERF3(
				startYear, endYear, deltaYears, sol, idens, ProbabilityModelOptions.U3_BPT,
				MagDependentAperiodicityOptions.LOW_VALUES);
		
		boolean supraSeismo = true;
		File eventFile = new File(simDir, "eqs.ALLCAL2_RSQSim_sigma0.5-5_b=0.015.long.barall");
		System.out.println("Loading events...");
		List<? extends SimulatorEvent> events;
		if (supraSeismo) {
			RuptureIdentifier supraSeismoIden = new SupraSeisRupIden(tools);
			RuptureIdentifier andIden = new LogicalAndRupIden(supraSeismoIden, new LogicalOrRupIden(idens));
			events = EQSIMv06FileReader.readEventsFile(eventFile, geoms, Lists.newArrayList(andIden));
//			events = EQSIMv06FileReader.readEventsFile(eventFile, geoms, Lists.newArrayList(supraSeismoIden));
		} else {
			events = EQSIMv06FileReader.readEventsFile(eventFile, geoms, idens);
		}
		
		EvenlyDiscretizedFunc simActual = calcSimulator(startYear, endYear, deltaYears, events, idens, null);
		EvenlyDiscretizedFunc simUnsynch = calcSimulator(startYear, endYear, deltaYears, events, idens, RandomDistType.ACTUAL);
		EvenlyDiscretizedFunc simPoisson = calcSimulator(startYear, endYear, deltaYears, events, idens, RandomDistType.PROBABILISTIC_SHUFFLE);
		
		plot(uc3Poisson, uc3PrefBlend, uc3High, uc3Mid, uc3Low, simActual, simUnsynch, simPoisson);
	}

}
