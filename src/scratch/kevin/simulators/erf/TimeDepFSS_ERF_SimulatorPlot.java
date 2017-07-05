package scratch.kevin.simulators.erf;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.iden.LogicalAndRupIden;
import org.opensha.sha.simulators.iden.MagRangeRuptureIdentifier;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.iden.SupraSeisRupIden;
import org.opensha.sha.simulators.parsers.EQSIMv06FileReader;
import org.opensha.sha.simulators.utils.General_EQSIM_Tools;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.utils.MatrixIO;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class TimeDepFSS_ERF_SimulatorPlot {

	public static void main(String[] args) throws IOException {
		int numTrialsPer = 1000;
//		File dir = new File("/tmp/2013_10_29-erf-audit-cov-0.3");
//		File dir = new File("/tmp/2013_10_31-erf-audit-cov-0.3-dur5");
//		File dir = new File("/tmp/2013_11_12-erf-audit-cov-LOW-dur50");
		File dir = new File("/tmp/2013_11_20-erf-audit-cov-LOW-dur50");
		String prefix = "erf_audit_";
		double magThresh = 7d;
		
		List<File> erfFiles = Lists.newArrayList();
		List<File> simFiles = Lists.newArrayList();
		
		for (File file : dir.listFiles()) {
			String name = file.getName();
			if (!name.startsWith(prefix) || !name.endsWith(".bin"))
				continue;
			if (name.contains("sim_occur"))
				simFiles.add(file);
			else if (name.contains("erf_occur"))
				erfFiles.add(file);
		}
		Preconditions.checkState(erfFiles.size() == simFiles.size());
		
		double[] simulatorOccurances = null;
		double[] forecastOccurances = null;
		
		for (int i=0; i<erfFiles.size(); i++) {
			if (simulatorOccurances == null) {
				simulatorOccurances = MatrixIO.doubleArrayFromFile(simFiles.get(i));
				forecastOccurances = MatrixIO.doubleArrayFromFile(erfFiles.get(i));
			} else {
				double[] newSimulatorOccurances = MatrixIO.doubleArrayFromFile(simFiles.get(i));
				double[] newForecastOccurances = MatrixIO.doubleArrayFromFile(erfFiles.get(i));
				for (int j=0; j<simulatorOccurances.length; j++) {
					simulatorOccurances[j] += newSimulatorOccurances[j];
					forecastOccurances[j] += newForecastOccurances[j];
				}
			}
		}
		
		System.out.println("Loaded in "+simFiles.size()+" files");
		int numTrials = simFiles.size()*numTrialsPer;
		
		double[] simPerWindows = new double[simulatorOccurances.length];
		double[] erfPerWindows = new double[forecastOccurances.length];
		for (int i=0; i<simPerWindows.length; i++) {
			simPerWindows[i] = simulatorOccurances[i]/(double)numTrials;
			erfPerWindows[i] = forecastOccurances[i]/(double)numTrials;
		}
		Arrays.sort(simPerWindows);
		Arrays.sort(erfPerWindows);
		System.out.println("sim per. min="+simPerWindows[0]+"\tmax="+simPerWindows[simPerWindows.length-1]
				+"\tavg="+StatUtils.mean(simPerWindows)+"\tmedian="+DataUtils.median_sorted(simPerWindows));
		System.out.println("erf per. min="+erfPerWindows[0]+"\tmax="+simPerWindows[erfPerWindows.length-1]
				+"\tavg="+StatUtils.mean(erfPerWindows)+"\tmedian="+DataUtils.median_sorted(erfPerWindows));
		
		List<XY_DataSet> elems = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		double max;
		if (magThresh > 0) {
			FaultSystemRupSet rupSet = getRupSet();
			DefaultXY_DataSet dataBelow = getNonZeroData(0d, magThresh, simulatorOccurances, forecastOccurances, rupSet);
			dataBelow.setName("Below M"+(float)magThresh);
			DefaultXY_DataSet dataAbove = getNonZeroData(magThresh, 10d, simulatorOccurances, forecastOccurances, rupSet);
			dataAbove.setName("Above M"+(float)magThresh);
			
			elems.add(dataBelow);
			elems.add(dataAbove);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 1f, Color.GREEN));
			chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 1f, Color.RED));
			
			max = dataBelow.getMaxY();
			max = Math.max(max, dataBelow.getMaxX());
			max = Math.max(max, dataAbove.getMaxY());
			max = Math.max(max, dataAbove.getMaxX());
		} else {
			DefaultXY_DataSet data = getNonZeroData(0d, 10d, simulatorOccurances, forecastOccurances, null);
			
			elems.add(data);
			chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 1f, Color.BLACK));
			
			max = data.getMaxY();
			max = Math.max(max, data.getMaxX());
		}
		
		ArbitrarilyDiscretizedFunc eventRatio = new ArbitrarilyDiscretizedFunc();
		eventRatio.set(0d, 0d);
		eventRatio.set(1e-1, 1e-1);
		eventRatio.set(max, max);
		
		elems.add(eventRatio);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		
		double numOccur = StatUtils.sum(simulatorOccurances);
		double numForecast = StatUtils.sum(forecastOccurances);
		System.out.println("Tot sim occurences: "+numOccur);
		System.out.println("Tot predicted: "+numForecast);
		
		GraphWindow gw = new GraphWindow(elems, "Time Dependent Simulator Audit ("+numTrials+" sims)", chars);
		gw.setX_AxisLabel("Actual RSQSim Occurances (sum="+(float)numOccur+")");
		gw.setY_AxisLabel("FSS ERF Predicted Occurances (sum="+(float)numForecast+")");
		Range range = new Range(0, max);
		gw.setAxisRange(range, range);
		gw.setSize(1000, 800);
		gw.saveAsPNG("/tmp/fss_erf_audit.png");
		gw.setXLog(true);
		gw.setYLog(true);
		range = new Range(1d, max);
		gw.setAxisRange(range, range);
		gw.saveAsPNG("/tmp/fss_erf_audit_log.png");
		
		
	}
	
	private static DefaultXY_DataSet getNonZeroData(double minMag, double maxMag, double[] simulatorOccurances ,
			double[] forecastOccurances, FaultSystemRupSet rupSet) {
		int numNonZeros = 0;
		for (int i=0; i<simulatorOccurances.length; i++) {
			if (rupSet != null) {
				double mag = rupSet.getMagForRup(i);
				if (mag < minMag || mag > maxMag)
					continue;
			}
			if (simulatorOccurances[i] > 0 && forecastOccurances[i] > 0)
				numNonZeros++;
		}
		double[] nonZeroOccur = new double[numNonZeros];
		double[] nonZeroPredict = new double[numNonZeros];
		int cnt = 0;
		for (int i=0; i<simulatorOccurances.length; i++) {
			if (rupSet != null) {
				double mag = rupSet.getMagForRup(i);
				if (mag < minMag || mag > maxMag)
					continue;
			}
			if (simulatorOccurances[i] > 0 && forecastOccurances[i] > 0) {
				nonZeroOccur[cnt] = simulatorOccurances[i];
				nonZeroPredict[cnt++] = forecastOccurances[i];
			}
		}
		Preconditions.checkState(cnt == numNonZeros);
		
		DefaultXY_DataSet ratioData = new DefaultXY_DataSet(nonZeroOccur, nonZeroPredict);
		
		return ratioData;
	}
	
	private static FaultSystemRupSet getRupSet() throws IOException {
		// ugly, copy pasted
		File dataDir = new File("/home/kevin/Simulators");
		File geomFile = new File(dataDir, "ALLCAL2_1-7-11_Geometry.dat");
		System.out.println("Loading geometry from "+geomFile.getAbsolutePath());
		General_EQSIM_Tools tools = new General_EQSIM_Tools(geomFile);
//		File eventFile = new File(dir, "eqs.ALLCAL2_RSQSim_sigma0.5-5_b=0.015.barall");
		File eventFile = new File(dataDir, "eqs.ALLCAL2_RSQSim_sigma0.5-5_b=0.015.long.barall");
		System.out.println("Loading events...");
		boolean supraSeismogenic = true;
		double minMag = 6.5;
		List<RuptureIdentifier> rupIdens = Lists.newArrayList();
		if (supraSeismogenic)
			rupIdens.add(new SupraSeisRupIden(tools));
		if (minMag > 0)
			rupIdens.add(new MagRangeRuptureIdentifier(minMag, 10d));
		LogicalAndRupIden andIden = new LogicalAndRupIden(rupIdens);
		rupIdens = Lists.newArrayList();
		rupIdens.add(andIden);
		
		List<? extends SimulatorEvent> events = EQSIMv06FileReader.readEventsFile(eventFile, tools.getElementsList(), rupIdens);
		
//		Region region = null;
		Region region = new CaliforniaRegions.RELM_SOCAL();
		
		if (region != null) {
			Map<Integer, Boolean> elementsInRegionsCache = Maps.newHashMap();
			
			// just uese centers since they're small enough elements
			for (SimulatorElement elem : tools.getElementsList())
				elementsInRegionsCache.put(elem.getID(), region.contains(elem.getCenterLocation()));
			
			List<SimulatorEvent> eventsInRegion = Lists.newArrayList();
			for (SimulatorEvent e : events) {
				for (int elemID : e.getAllElementIDs()) {
					if (elementsInRegionsCache.get(elemID)) {
						eventsInRegion.add(e);
						break;
					}
				}
			}
			
			System.out.println(eventsInRegion.size()+"/"+events.size()+" events in region: "+region.getName());
			events = eventsInRegion;
		}
		
//		if (supraSeismogenic) {
//			List<EQSIM_Event> supraSeisEvents = Lists.newArrayList();
//			for (EQSIM_Event e : events)
//				if (tools.isEventSupraSeismogenic(e, Double.NaN))
//					supraSeisEvents.add(e);
//			events = supraSeisEvents;
//		}
		
		double durationYears = General_EQSIM_Tools.getSimulationDurationYears(events);
		
		List<SimulatorElement> elements = tools.getElementsList();
		SubSectionBiulder subSectBuilder = new SubSectionBiulder(elements);
		FaultSystemRupSet rupSet = SimulatorFaultSystemSolution.buildRupSet(elements, events, durationYears, subSectBuilder);
		FaultSystemSolution sol = new SimulatorFaultSystemSolution(rupSet, subSectBuilder, events, durationYears);
		System.out.println("Precombine sol has "+rupSet.getNumRuptures()+" rups");
		sol = TimeDepFSS_ERF_Simulator_Test.combineIdenticalRups(sol);
		rupSet = sol.getRupSet();
		System.out.println("Combined sol has "+rupSet.getNumRuptures()+" rups");
		return rupSet;
	}

}
