package scratch.kevin.simulators;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.iden.ElementMagRangeDescription;
import org.opensha.sha.simulators.iden.MagRangeRuptureIdentifier;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.parsers.EQSIMv06FileReader;
import org.opensha.sha.simulators.utils.General_EQSIM_Tools;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class ReasenbergAndJonesComparison {
	
	// simulator times are in seconds
	public static int SECONDS_IN_DAY = 24*60*60;
	
	private List<? extends SimulatorEvent> events;
	private RuptureIdentifier rupIden;
	private double minDays;
	private int maxDays;
	
	private List<? extends SimulatorEvent> matches;
	private EvenlyDiscretizedFunc countFunc;
	
	public ReasenbergAndJonesComparison(
			List<? extends SimulatorEvent> events, RuptureIdentifier rupIden, double minDays, int maxDays) {
		this.events = events;
		this.rupIden = rupIden;
		this.minDays = minDays;
		this.maxDays = maxDays;
		
		update();
	}
	
	private void update() {
		matches = rupIden.getMatches(events);
		
		countFunc = new EvenlyDiscretizedFunc(1d, maxDays, 1d);
		
		double maxDuration = SECONDS_IN_DAY * maxDays;
		for (SimulatorEvent match : matches) {
			double startTime = match.getTime()+match.getDuration();
			double endTime = startTime+maxDuration;
			startTime += minDays*SECONDS_IN_DAY;
			
			double matchMag = match.getMagnitude();
			int eventIndex = Collections.binarySearch(events, match);
			for (int i=eventIndex+1; i<events.size(); i++) {
				SimulatorEvent e = events.get(i);
				
				if (match.getID() == e.getID())
					continue;
				
				double mag = e.getMagnitude();
				if (mag < matchMag)
					continue;
				
				double time = e.getTime();
				
				if (time < startTime)
					continue;
				
				if (time > endTime)
					break;
				
				double delta = time - startTime;
				
//				System.out.println(e.getID()+" delta="+delta);
				for (int index=countFunc.size(); --index>=0;) {
					double secs = countFunc.getX(index) * SECONDS_IN_DAY;
//					System.out.println((index+1)+" days. "+delta+" <= "+secs+" ? "+(delta <= secs));
					if (delta <= secs)
						countFunc.set(index, countFunc.getY(index)+1);
					else
						break;
				}
				break;
			}
		}
		
		String info = null;
		for (int index=0; index<countFunc.size(); index++) {
			double cnt = countFunc.getY(index);
			double rate = cnt / (double)matches.size();
			if (info == null)
				info = "Raw Counts:";
			info += "\n\t";
			int days = (index+1);
			String line = days+" days: "+(int)cnt+"/"+matches.size()+" = "+(float)rate;
			info += line;
//			System.out.println(line);
//			double prob = 1d-Math.exp(-rate);
			countFunc.set(index, rate);
//			countFunc.set(index, prob);
		}
		countFunc.setInfo(info);
	}
	
	public int getNumMatches() {
		return matches.size();
	}

	public EvenlyDiscretizedFunc getCountFunc() {
		return countFunc;
	}
	
	private static double a = -1.67;
//	private static double a_min = -3.3;
//	private static double a_max = -0.7;
	private static double a_min = -2.35;
	private static double a_max = -1.2;
	
	// b is not used because Mm = M
	
	private static double c = 0.05;
	
	private static double p = 1.08;
//	private static double p_min = 0.5;
//	private static double p_max = 1.55;
	private static double p_min = 0.85;
	private static double p_max = 1.3;
	
	public EvenlyDiscretizedFunc getRJRateFunc() {
		EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(1, maxDays, 1d);
		func.setName("Reasenberg & Jones Eqn 3");
		
		for (int day=1; day<=maxDays; day++) {
			double val = Math.pow(10, a)*Math.pow(day+c, -p);
			func.set(day-1, val);
		}
		
		return func;
	}
	
	public EvenlyDiscretizedFunc getRJIntegralFunc() {
		return getRJIntegralFunc(a, c, p, minDays, maxDays);
	}
	
	public static EvenlyDiscretizedFunc getRJIntegralFunc(double a, double c, double p, double minDays, int maxDays) {
		EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(1, maxDays, 1d);
		func.setName("Reasenberg & Jones Eqn 4 (a="+a+", c="+c+", p="+p+", min="+minDays+", max="+maxDays+")");
		
		double integralTAtMinT = getIntegralPart(c, p, minDays);
		
		for (int day=1; day<=maxDays; day++) {
			double val = Math.pow(10, a)*(getIntegralPart(c, p, day)-integralTAtMinT);
			func.set(day-1, val);
		}
		
		return func;
	}
	
	private static ArrayList<DiscretizedFunc> getRJIntegralBounds(double minDays, int maxDays) {
		ArrayList<DiscretizedFunc> funcs = Lists.newArrayList();
		
		funcs.add(getRJIntegralFunc(a_min, c, p_min, minDays, maxDays));
		funcs.add(getRJIntegralFunc(a_min, c, p, minDays, maxDays));
		funcs.add(getRJIntegralFunc(a_min, c, p_max, minDays, maxDays));
		
		funcs.add(getRJIntegralFunc(a, c, p_min, minDays, maxDays));
		funcs.add(getRJIntegralFunc(a, c, p, minDays, maxDays));
		funcs.add(getRJIntegralFunc(a, c, p_max, minDays, maxDays));
		
		funcs.add(getRJIntegralFunc(a_max, c, p_min, minDays, maxDays));
		funcs.add(getRJIntegralFunc(a_max, c, p, minDays, maxDays));
		funcs.add(getRJIntegralFunc(a_max, c, p_max, minDays, maxDays));
		
		return funcs;
	}
	
	private static double getIntegralPart(double c, double p, double t) {
		return Math.pow(t+c, 1d-p) / (1d-p);
	}
	
	private EvenlyDiscretizedFunc findBestFitRJ(double delta) {
		EvenlyDiscretizedFunc best_func = null;
		double best_misfit = Double.MAX_VALUE;
		for (double a=a_min; a<=a_max+delta; a+=delta) {
			for (double p=p_min; p<=p_max+delta; p+=delta) {
//				for (double c=0; c<=1; c+=0.01) {
					EvenlyDiscretizedFunc rj = getRJIntegralFunc(a, c, p, minDays, maxDays);
					
					double misfit = 0;
					for (int i=0; i<rj.size(); i++) {
						double diff = rj.getY(i) - countFunc.getY(i);
						misfit += diff*diff;
					}
					
					if (misfit < best_misfit) {
						best_func = rj;
						best_misfit = misfit;
					}
//				}
			}
		}
		
		return best_func;
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
		double minDays = 0.01;
		int maxDays = 365;
		
		RuptureIdentifier rupIden = new ElementMagRangeDescription("SAF Mojave 7+",
				ElementMagRangeDescription.SAF_MOJAVE_ELEMENT_ID, 7d, 10d);
		ReasenbergAndJonesComparison rj = new ReasenbergAndJonesComparison(
				events, rupIden, minDays, maxDays);
		
		EvenlyDiscretizedFunc rjIntegralComp = rj.getRJIntegralFunc();
		
		ArrayList<DiscretizedFunc> funcs = Lists.newArrayList();
		funcs.add(rj.getCountFunc());
		funcs.add(rjIntegralComp);
		funcs.add(rj.findBestFitRJ(0.05));
		GraphWindow gw = new GraphWindow(funcs,
				"Rate Aftershcoks Larger Than Mainshock (M7+ on SAF Mojave)", null, false);
		gw.setX_AxisLabel("Days");
		gw.setY_AxisLabel("Rate");
		gw.setVisible(true);
		
		rupIden = new ElementMagRangeDescription("SAF Coachella 7+",
				ElementMagRangeDescription.SAF_COACHELLA_ELEMENT_ID, 7d, 10d);
		rj = new ReasenbergAndJonesComparison(
				events, rupIden, minDays, maxDays);
		
		funcs = Lists.newArrayList();
		funcs.add(rj.getCountFunc());
		funcs.add(rjIntegralComp);
		gw = new GraphWindow(funcs,
				"Rate Aftershcoks Larger Than Mainshock (M7+ on SAF Coachella)", null, false);
		gw.setX_AxisLabel("Days");
		gw.setY_AxisLabel("Rate");
		gw.setVisible(true);
		
		rupIden = new ElementMagRangeDescription("SAF Carrizo 7+",
				ElementMagRangeDescription.SAF_CARRIZO_ELEMENT_ID, 7d, 10d);
		rj = new ReasenbergAndJonesComparison(
				events, rupIden, minDays, maxDays);
		
		funcs = Lists.newArrayList();
		funcs.add(rj.getCountFunc());
		funcs.add(rjIntegralComp);
		gw = new GraphWindow(funcs,
				"Rate Aftershcoks Larger Than Mainshock (M7+ on SAF Carrizo)", null, false);
		gw.setX_AxisLabel("Days");
		gw.setY_AxisLabel("Rate");
		gw.setVisible(true);
		
		rupIden = new MagRangeRuptureIdentifier(7d, 10d);
		rj = new ReasenbergAndJonesComparison(
				events, rupIden, minDays, maxDays);
		
		funcs = Lists.newArrayList();
		funcs.add(rj.getCountFunc());
		funcs.add(rjIntegralComp);
		gw = new GraphWindow(funcs,
				"Rate Aftershcoks Larger Than Mainshock (M7+)", null, false);
		gw.setX_AxisLabel("Days");
		gw.setY_AxisLabel("Rate");
		gw.setVisible(true);
		
		rupIden = new MagRangeRuptureIdentifier(6d, 10d);
		rj = new ReasenbergAndJonesComparison(
				events, rupIden, minDays, maxDays);
		
		funcs = Lists.newArrayList();
		funcs.add(rj.getCountFunc());
		funcs.add(rjIntegralComp);
		gw = new GraphWindow(funcs,
				"Rate Aftershcoks Larger Than Mainshock (M6+)", null, false);
		gw.setX_AxisLabel("Days");
		gw.setY_AxisLabel("Rate");
		gw.setVisible(true);
		
		rupIden = new MagRangeRuptureIdentifier(5d, 10d);
		rj = new ReasenbergAndJonesComparison(
				events, rupIden, minDays, maxDays);
		
		funcs = Lists.newArrayList();
		funcs.add(rj.getCountFunc());
		funcs.add(rjIntegralComp);
		gw = new GraphWindow(funcs,
				"Rate Aftershcoks Larger Than Mainshock (M5+)", null, false);
		gw.setX_AxisLabel("Days");
		gw.setY_AxisLabel("Rate");
		gw.setVisible(true);
		
		
		// RJ bounds
		funcs = getRJIntegralBounds(minDays, maxDays);
		gw = new GraphWindow(funcs,
				"Rate Aftershcoks Larger Than Mainshock (R&J Bounds)", null, false);
		gw.setX_AxisLabel("Days");
		gw.setY_AxisLabel("Rate");
		gw.setVisible(true);
	}

}
