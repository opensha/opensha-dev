package scratch.kevin.simulators;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.iden.ElementMagRangeDescription;
import org.opensha.sha.simulators.iden.EventsInWindowsMatcher;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.parsers.EQSIMv06FileReader;
import org.opensha.sha.simulators.utils.General_EQSIM_Tools;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public class OmorisLawComparison {
	
	private static void doComparison(List<SimulatorEvent> events, RuptureIdentifier rupIden,
			int maxDays, double magBin, double binWidth) {
		
	}

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		File dir = new File("/home/kevin/Simulators");
		File geomFile = new File(dir, "ALLCAL2_1-7-11_Geometry.dat");
		System.out.println("Loading geometry...");
		General_EQSIM_Tools tools = new General_EQSIM_Tools(geomFile);
//		File eventFile = new File(dir, "eqs.ALLCAL2_RSQSim_sigma0.5-5_b=0.015.barall");
		File eventFile = new File(dir, "eqs.ALLCAL2_RSQSim_sigma0.5-5_b=0.015.long.barall");
		System.out.println("Loading events...");
		List<? extends SimulatorEvent> events = EQSIMv06FileReader.readEventsFile(eventFile, tools.getElementsList());
		System.out.println("Calculating...");
		double minWindowDays = 1d/24d;
		int maxDays = 365;
		double binWidth = 0.1;
		
		Region reg = new CaliforniaRegions.RELM_SOCAL();
		HashSet<Integer> elementsInRegion = MFDCalc.getElementsInsideRegion(tools.getElementsList(), reg);
		
		double totalEventDuration = General_EQSIM_Tools.getSimulationDurationYears(events);
		
		double daysPerYear = BatchPlotGen.DAYS_PER_YEAR;
		
		double minWindowYears = minWindowDays / daysPerYear;
		
		boolean[] randomizes = { false, true };
		
//		List<Double> magBins = Lists.newArrayList(7d, 7.6d);
		List<Double> magBins = Lists.newArrayList(7.6d);
		
		List<RuptureIdentifier> rupIdens = Lists.newArrayList();
		rupIdens.add(new ElementMagRangeDescription("SAF Mojave 7+",
				ElementMagRangeDescription.SAF_MOJAVE_ELEMENT_ID, 7d, 10d));
		rupIdens.add(new ElementMagRangeDescription("SAF Coachella 7+",
				ElementMagRangeDescription.SAF_COACHELLA_ELEMENT_ID, 7d, 10d));
		rupIdens.add(new ElementMagRangeDescription("SAF Carrizo 7+",
				ElementMagRangeDescription.SAF_CARRIZO_ELEMENT_ID, 7d, 10d));
		rupIdens.add(new ElementMagRangeDescription("SAF Cholame 7+",
				ElementMagRangeDescription.SAF_CHOLAME_ELEMENT_ID, 7d, 10d));
		rupIdens.add(new ElementMagRangeDescription("Garlock 7+",
				ElementMagRangeDescription.GARLOCK_WEST_ELEMENT_ID, 7d, 10d));
		rupIdens.add(new ElementMagRangeDescription("San Jacinto 7+",
				ElementMagRangeDescription.SAN_JACINTO__ELEMENT_ID, 7d, 10d));
		
		for (int i=0; i<rupIdens.size(); i++) {
			RuptureIdentifier rupIden = rupIdens.get(i);
			String idenName = rupIden.getName();
			String fsafeName = PeriodicityPlotter.getFileSafeString(idenName);
			
			for (boolean randomize : randomizes) {
				double magBin = 7.6;
				IncrementalMagFreqDist indepMFD = MFDCalc.calcMFD(events, elementsInRegion,
						totalEventDuration*daysPerYear, magBin, 1, binWidth);
				
				ArbitrarilyDiscretizedFunc eventFunc = new ArbitrarilyDiscretizedFunc();
				ArbitrarilyDiscretizedFunc omoriComp = new ArbitrarilyDiscretizedFunc();
				
				Preconditions.checkState(indepMFD.size() == 1);
				Preconditions.checkState(indepMFD.getX(0) == magBin);
				double indepVal = indepMFD.getY(0);
				
				double omoriK = -1;
				
				CSVFile<String> csv = new CSVFile<String>(true);
				double csvMinMag = 6;
				double csvMaxMag = 8;
				double csvDelta = 0.1;
				int csvMagNum = (int)((csvMaxMag-csvMinMag)/csvDelta) + 1;
				
				IncrementalMagFreqDist mfd = new IncrementalMagFreqDist(csvMinMag, csvMaxMag, csvMagNum);
				
				double myCSVMin = csvMinMag-0.5*csvDelta;
				double myCSVMax = mfd.getMaxX()+0.5*csvDelta;
				
				List<String> header = Lists.newArrayList("");
				for (int m=0; m<mfd.size(); m++)
					header.add(""+(float)mfd.getX(m));
				csv.addLine(header);
				
				double secondsPerDay = 60*60*24;
				
				double windowLenYrs = (double)maxDays / daysPerYear;
				EventsInWindowsMatcher match =
						new EventsInWindowsMatcher(events, rupIden, minWindowYears, windowLenYrs, randomize);
				
				List<SimulatorEvent> matches = match.getEventsInWindows();
				List<Double> matchTimes = match.getEventTimesFromWindowStarts();
				
				for (int days=1; days<=maxDays; days++) {
//					if (days > 100) {
//						if (days % 5 != 0)
//							continue;
//					} else if (days > 50) {
//						if (days % 3 != 0)
//							continue;
//					} else if (days > 10) {
//						if (days % 2 != 0)
//							continue;
//					}double windowLenYrs = (double)days / daysPerYear;
					
					double maxDurationSecs = days*secondsPerDay;
					double minDurationSecs = maxDurationSecs-secondsPerDay;
					
					for (int m=0; m<mfd.size(); m++)
						mfd.set(m, 0d);
					
					List<SimulatorEvent> eventsInWindows = Lists.newArrayList();
					
					for (int m=0; m<matches.size(); m++) {
						double timeFromStart = matchTimes.get(m);
						if (timeFromStart >= maxDurationSecs)
							continue;
						SimulatorEvent e = matches.get(m);
						eventsInWindows.add(e);
						if (timeFromStart < minDurationSecs)
							continue;
						double mag = e.getMagnitude();
						if (mag < myCSVMin || mag > myCSVMax)
							continue;
						int ind = mfd.getClosestXIndex(mag);
						mfd.set(ind, mfd.getY(ind)+1);
					}
					
					List<String> line = Lists.newArrayList(days+"");
					for (int m=0; m<mfd.size(); m++)
						line.add((int)mfd.getY(m)+"");
					csv.addLine(line);
					
					IncrementalMagFreqDist depMFD = MFDCalc.calcMFD(eventsInWindows, elementsInRegion,
							0, magBin, 1, binWidth);
					Preconditions.checkState(depMFD.size() == 1);
					Preconditions.checkState(depMFD.getX(0) == magBin);
					
					double t = days;
					double t0 = minWindowDays;
					
					double myEvents = depMFD.getY(0);
//					double myProbGain = depMFD.getY(0)/indepVal;
//					System.out.println("Days:\t"+days+"\tgain:"+myProbGain);
					System.out.println("Days:\t"+days+"\tcml events:"+myEvents);
					eventFunc.set((double)days, myEvents);
					if (omoriK < 0)
						omoriK = myEvents;
//					omoriComp.set((double)days, omoriK/(double)days);
					omoriComp.set((double)days, omoriK + Math.log(t) / t0);
				}
				
				File writeDir = new File("/home/kevin/Simulators/omori_csv");
				
				if (randomize)
					csv.writeToFile(new File(writeDir, "omori_"+fsafeName+"_mfds_randomized.csv"));
				else
					csv.writeToFile(new File(writeDir, "omori_"+fsafeName+"_mfds.csv"));
				
				ArrayList<DiscretizedFunc> funcs = Lists.newArrayList();
				funcs.add(eventFunc);
				funcs.add(omoriComp);
//				GraphWindow gw = new GraphWindow(funcs, "Omori's Law Comparison for M="+magBin);
			}
		}
	}

}
