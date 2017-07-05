package scratch.kevin.simulators.dists;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;

import org.jfree.chart.plot.XYPlot;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.gui.plot.GraphPanel;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotElement;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotWindow;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.iden.ElementMagRangeDescription;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.parsers.EQSIMv06FileReader;
import org.opensha.sha.simulators.utils.General_EQSIM_Tools;

import scratch.UCERF3.utils.IDPairing;
import scratch.kevin.simulators.PeriodicityPlotter;
import scratch.kevin.simulators.catBuild.RandomCatalogBuilder;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public class FollowerReturnPeriodProvider implements
		ProbabilisticReturnPeriodProvider {
	
	private RuptureIdentifier driver;
	private RuptureIdentifier follower;
//	private HistogramFunction condProbFunc;
	// x is auto-inter event time, y is cross-inter event time
//	private EvenlyDiscrXYZ_DataSet jointProbDataset;
//	private EvenlyDiscrXYZ_DataSet condProbDataset;
	private EvenlyDiscrXYZ_DataSet autoRupProb;
	private EvenlyDiscrXYZ_DataSet autoGivenCrossRupProb;
	private EvenlyDiscrXYZ_DataSet autoHitsData;
	private EvenlyDiscrXYZ_DataSet autoGivenCrossHitsData;
	private EvenlyDiscrXYZ_DataSet totalStateData;
	private double maxVal;
	private HistogramFunction refHist;
	
	private boolean doSmoothing = true;
	private static boolean plotsEnabled = false;
	
	private HistogramFunction followerIndepCumDist;
	
	private ActualDistReturnPeriodProvider driverActualDist;
	private ActualDistReturnPeriodProvider followerActualDist;
	
	// these will be used at the start of simulations where no prior has been found
	private double fakeStartPrevDriverTime;
	private double fakeStartPrevFollowerTime;
	
	private static final double FLAG_REVERT_REG_DIST = -1234.56;
	
	private long dep_count = 0;
	private long fallback_count = 0;
//	private long nan_fallback_count = 0;
	private long tot_count = 0;
	private EvenlyDiscrXYZ_DataSet hitDataset;
	
	public FollowerReturnPeriodProvider(List<? extends SimulatorEvent> events, RuptureIdentifier driver,
			RuptureIdentifier follower, double distDeltaYears, double maxTimeDiff) {
		this(events, driver, driver.getMatches(events), follower, follower.getMatches(events), distDeltaYears, maxTimeDiff);
	}
	
	public FollowerReturnPeriodProvider(List<? extends SimulatorEvent> events, RuptureIdentifier driver, List<? extends SimulatorEvent> driverMatches,
			RuptureIdentifier follower, List<? extends SimulatorEvent> followerMatches, double distDeltaYears, double maxTimeDiff) {
		this.driver = driver;
		this.follower = follower;
		driverActualDist = new ActualDistReturnPeriodProvider(PeriodicityPlotter.getRPs(driverMatches));
		followerActualDist = new ActualDistReturnPeriodProvider(PeriodicityPlotter.getRPs(followerMatches));
		
		// create PDF
		int num = (int)(maxTimeDiff/distDeltaYears - 1);
		double[] discr_vals = null;
		// organized by auto, then each hist contains cross
		HistogramFunction[] hists = new HistogramFunction[num];
		for (int i=0; i<num; i++) {
			HistogramFunction hist = new HistogramFunction(0.5*distDeltaYears, num, distDeltaYears);
			hists[i] = hist;
			if (discr_vals == null) {
				maxVal = hist.getMaxX()+0.5*distDeltaYears;
				discr_vals = new double[num];
				for (int j=0; j<num; j++)
					discr_vals[j] = hist.getX(j);
			}
		}
		refHist = hists[0];
		
		HistogramFunction followerIndepDist = new HistogramFunction(0.5*distDeltaYears, num, distDeltaYears);
		
//		double iVal = PRECISION_SCALE * (x - minX) / delta;
//		int i = (delta == 0) ? 0 : (int) Math.round(iVal);
//		return (i<0) ? 0 : (i>=num) ? num-1 : i;
		
//		List<EQSIM_Event> combEvents = Lists.newArrayList();
//		combEvents.addAll(driverMatches);
//		combEvents.addAll(followerMatches);
//		Collections.sort(combEvents);
		List<? extends SimulatorEvent> combEvents = events;
		
		// TODO lets try this damn one
		double maxTime = events.get(events.size()-1).getTimeInYears();
		double startTime = events.get(0).getTimeInYears();
		int numSteps = (int)((maxTime - startTime)/distDeltaYears);
		
		List<SimulatorEvent> prevEvents = Lists.newArrayList();
		int eventIndex = 0;
		
		autoHitsData = new EvenlyDiscrXYZ_DataSet(num, num, discr_vals[0], discr_vals[0], distDeltaYears);
		autoGivenCrossHitsData = new EvenlyDiscrXYZ_DataSet(num, num, discr_vals[0], discr_vals[0], distDeltaYears);
		totalStateData = new EvenlyDiscrXYZ_DataSet(num, num, discr_vals[0], discr_vals[0], distDeltaYears);
		
		HistogramFunction followerIndepHits = new HistogramFunction(0.5*distDeltaYears, num, distDeltaYears);
		HistogramFunction followerIndepTot = new HistogramFunction(0.5*distDeltaYears, num, distDeltaYears);
		followerIndepCumDist = new HistogramFunction(0.5*distDeltaYears, num, distDeltaYears);
		
		for (int step=0; step<numSteps; step++) {
			double windowStart = startTime + distDeltaYears*step;
			double windowEnd = windowStart + distDeltaYears;
			
			for (int i=eventIndex; i<events.size(); i++) {
				double time = events.get(i).getTimeInYears();
				Preconditions.checkState(time >= windowStart);
				if (time > windowEnd)
					break;
				prevEvents.add(events.get(i));
				eventIndex = i+1;
			}
			
			double prevDriverTime = Double.NaN;
			double prevPrevDriverTime = Double.NaN;
			double prevFollowerTime = Double.NaN;
			double prevPrevFollowerTime = Double.NaN;
			for (int i=prevEvents.size(); --i>=0;) {
				SimulatorEvent e = prevEvents.get(i);
				double eTime = e.getTimeInYears();
				Preconditions.checkState(!Double.isNaN(eTime));
				if (driver.isMatch(e)) {
					if (Double.isNaN(prevDriverTime)) {
						prevDriverTime = eTime;
					} else if (Double.isNaN(prevPrevDriverTime)) {
						prevPrevDriverTime = eTime;
					}
				}
				if (follower.isMatch(e)) {
					if (Double.isNaN(prevFollowerTime)) {
						prevFollowerTime = eTime;
					} else if (Double.isNaN(prevPrevFollowerTime)) {
						prevPrevFollowerTime = eTime;
					}
				}
				if (!(Double.isNaN(prevPrevFollowerTime) || Double.isNaN(prevFollowerTime)
						|| Double.isNaN(prevPrevDriverTime) || Double.isNaN(prevDriverTime)))
					break;
			}
//			System.out.println(prevEvents.size()+", "+eventIndex+"/"+events.size()+", "+prevFollowerTime+", "+prevDriverTime);
			if (Double.isNaN(prevFollowerTime) || Double.isNaN(prevDriverTime))
				continue;
			boolean autoHit = false;
			boolean crossHit = false;
			double origPrevTime = prevFollowerTime;
			if (prevFollowerTime >= windowStart && prevFollowerTime <= windowEnd) {
				autoHit = true;
				prevFollowerTime = prevPrevFollowerTime;
			}
			if (prevDriverTime >= windowStart && prevDriverTime <= windowEnd) {
				crossHit = true;
				prevDriverTime = prevPrevDriverTime;
			}
			double autoDelta = windowEnd - prevFollowerTime;
			double crossDelta = windowEnd - prevDriverTime;
			
//			System.out.println(autoDelta+","+crossDelta);
			
			if (autoDelta < maxVal && crossDelta < maxVal) {
				int autoIndex = refHist.getClosestXIndex(autoDelta);
				int crossIndex = refHist.getClosestXIndex(crossDelta);
				
//				System.out.println(autoDelta+", "+crossDelta+": "+autoIndex+", "+crossIndex
//						+". hit="+hit+", pFolDelta="+(windowEnd-origPrevTime)+", ppFolDelta="+(windowEnd-prevPrevFollowerTime));
				Preconditions.checkState(crossIndex<3 || autoIndex<3 || totalStateData.get(autoIndex-1, crossIndex-1)>0);
//				if (crossIndex>1 && autoIndex>1 && totData.get(autoIndex-1, crossIndex-1)==0) {
//					System.out.println("I'm isolated!!!! WTF?");
//					System.exit(0);
//				}
				
//				System.out.println(autoIndex+", "+crossIndex+" ("+autoDelta+","+crossDelta+")");
				
				totalStateData.set(autoIndex, crossIndex, totalStateData.get(autoIndex, crossIndex)+1d);
				if (autoHit)
					autoHitsData.set(autoIndex, crossIndex, autoHitsData.get(autoIndex, crossIndex)+1d);
				if (autoHit && crossHit)
					autoGivenCrossHitsData.set(autoIndex, crossIndex, autoGivenCrossHitsData.get(autoIndex, crossIndex)+1d);
			}
			if (autoDelta < maxVal) {
				int autoIndex = followerIndepHits.getClosestXIndex(autoDelta);
				followerIndepTot.set(autoIndex, followerIndepTot.getY(autoIndex)+1d);
				if (autoHit)
					followerIndepHits.set(autoIndex, followerIndepHits.getY(autoIndex)+1d);
			}
		}
		
		for (int i=0; i<followerIndepHits.size(); i++) {
			followerIndepCumDist.set(i, followerIndepHits.getY(i)/followerIndepTot.getY(i));
		}
		
		showDist(autoHitsData, (double)Float.MIN_VALUE, autoHitsData.getMaxZ(), follower.getName()+" Hits");
		showDist(autoGivenCrossHitsData, (double)Float.MIN_VALUE, autoHitsData.getMaxZ(), follower.getName()
				+" given "+driver.getName()+" Hits");
		showDist(totalStateData, (double)Float.MIN_VALUE, totalStateData.getMaxZ(), "States");
		
		if (doSmoothing) {
			autoHitsData = getSmoothedDataset(autoHitsData, totalStateData, 50d);
			autoGivenCrossHitsData = getSmoothedDataset(autoGivenCrossHitsData, totalStateData, 50d);
			totalStateData = getSmoothedDataset(totalStateData, totalStateData, 50d);
			
//			showDist(autoHitsData, (double)Float.MIN_VALUE, autoHitsData.getMaxZ(), "Smoothed "+follower.getName()+" Hits");
//			showDist(autoHitsData, (double)Float.MIN_VALUE, autoHitsData.getMaxZ(), "Smoothed "+follower.getName()
//					+" given "+driver.getName()+" Hits");
//			showDist(totalStateData, (double)Float.MIN_VALUE, totalStateData.getMaxZ(), "Smoothed States");
		}

		int numViolations = 0;
		for (int xInd=0; xInd<autoHitsData.getNumX(); xInd++) {
			for (int yInd=0; yInd<autoHitsData.getNumY(); yInd++) {
				if ((float)autoHitsData.get(xInd, yInd) > (float)totalStateData.get(xInd, yInd))
					numViolations++;
			}
		}
		System.out.println("Violations: "+numViolations);
		
//		try {
//			Thread.sleep(1000000000000000l);
//		} catch (InterruptedException e1) {
//			// TODO Auto-generated catch block
//			e1.printStackTrace();
//		}
//		System.exit(0);
		
		autoRupProb = new EvenlyDiscrXYZ_DataSet(num, num, discr_vals[0], discr_vals[0], distDeltaYears);
		autoGivenCrossRupProb = new EvenlyDiscrXYZ_DataSet(num, num, discr_vals[0], discr_vals[0], distDeltaYears);
		for (int xInd=0; xInd<autoRupProb.getNumX(); xInd++) {
			for (int yInd=0; yInd<autoRupProb.getNumY(); yInd++) {
				double hits = autoHitsData.get(xInd, yInd);
				double tot = totalStateData.get(xInd, yInd);
				
				double prob = hits / tot;
				autoRupProb.set(xInd, yInd, prob);
				
				double hitsGivenCross = autoGivenCrossHitsData.get(xInd, yInd);
				double probGivenCross = hitsGivenCross / tot;
				autoGivenCrossRupProb.set(xInd, yInd, probGivenCross);
			}
		}
		
		double medianDriverRP = DataUtils.median(driverActualDist.getRPs());
		double medianFollowerRP = DataUtils.median(followerActualDist.getRPs());
		
		// if the driver RP is longer than the follower RP, then fill "revert" spots
		// in reverse
		boolean flipForFilling = medianDriverRP > medianFollowerRP;
		
		if (flipForFilling) {
			EvenlyDiscrXYZ_DataSet temp =
					new EvenlyDiscrXYZ_DataSet(num, num, discr_vals[0], discr_vals[0], distDeltaYears);
			for (int xInd=0; xInd<num; xInd++)
				for (int yInd=0; yInd<num; yInd++)
					temp.set(xInd, yInd, autoRupProb.get(yInd, xInd));
			autoRupProb = temp;
		}
		
		// fill in top right corner of dist
		// first find x val with val at greatest y
		int highestYInd = 0;
		int startBuffXInd = 0;
		int[] maxYInds = new int[autoRupProb.getNumX()];
		for (int i=0; i<maxYInds.length; i++)
			maxYInds[i] = -1;
		for (int xInd=0; xInd<autoRupProb.getNumX(); xInd++) {
			for (int yInd=autoRupProb.getNumY(); --yInd>0;) {
				if (autoRupProb.get(xInd, yInd) > 0) {
					maxYInds[xInd] = yInd;
					if (yInd > highestYInd) {
						highestYInd = yInd;
						startBuffXInd = xInd;
					}
					break;
				}
			}
		}
		highestYInd++;
		startBuffXInd++;
		for (int i=0; i<maxYInds.length; i++)
			if (maxYInds[i] >= 0)
				maxYInds[i]++;
//		// move up in a diagonal
//		for (int startXInd=startBuffXInd; startXInd<jointProbDataset.getNumX(); startXInd++) {
//			int startYInd = maxYInds[startXInd];
//			for (int xInd=startXInd,yInd=startYInd; ++xInd<jointProbDataset.getNumX() && ++yInd<jointProbDataset.getNumY();)
//				cumCondProbDataset.set(xInd, yInd, 1d);
//		}
		
		// fill vertically
//		for (int x=startBuffXInd; x<jointProbDataset.getNumX(); x++)
//			for (int y=maxYInds[x]; y<jointProbDataset.getNumY(); y++)
//				cumCondProbDataset.set(x, y, 1d);
		
		// move down diagonals
		// first along top
//		int startXProjectedTop = startBuffXInd + (jointProbDataset.getNumX()-1-highestYInd);
//		for (int startXInd=startXProjectedTop; startXInd<jointProbDataset.getNumX(); startXInd++) {
//			int startYInd = jointProbDataset.getNumY()-1;
//			for (int xInd=startXInd,yInd=startYInd; xInd>=startBuffXInd && yInd>=0 && yInd>maxYInds[xInd]; xInd--,yInd--) {
//				// check if we've hit a real value
//				double curVal = cumCondProbDataset.get(xInd, yInd);
//				if (curVal != 0 && !Double.isNaN(curVal))
//					break;
////				// check if we've gone too far
////				int xDiff = xInd - startBuffXInd;
////				int yDiff = highestYInd - yInd;
////				if (xDiff < 0 || yDiff > xDiff)
////					break;
//				cumCondProbDataset.set(xInd, yInd, 1d);
//			}
//		}
		int startXProjectedTop = 0;
		for (int startXInd=startXProjectedTop; startXInd<autoRupProb.getNumX(); startXInd++) {
			int startYInd = autoRupProb.getNumY()-1;
			for (int xInd=startXInd,yInd=startYInd; xInd>=0 && yInd>=0 && yInd>maxYInds[xInd]; xInd--,yInd--) {
				// check if we've hit a real value
				double curVal = autoRupProb.get(xInd, yInd);
				if (curVal != 0 && !Double.isNaN(curVal))
					break;
//				// check if we've gone too far
//				int xDiff = xInd - startBuffXInd;
//				int yDiff = highestYInd - yInd;
//				if (xDiff < 0 || yDiff > xDiff)
//					break;
				if (xInd < startBuffXInd) {
					// filling in top left, see if we should stop
					if (yInd < highestYInd + (startBuffXInd-xInd))
						break;
				}
				autoRupProb.set(xInd, yInd, FLAG_REVERT_REG_DIST);
			}
		}
		// now along right
		int startYProjectedY = autoRupProb.getNumY()-1;
		while (startXProjectedTop >= autoRupProb.getNumX()) {
			startXProjectedTop--;
			startYProjectedY--;
		}
		for (int startYInd=startYProjectedY; startYInd>=0; startYInd--) {
			int startXInd = autoRupProb.getNumX()-1;
			for (int xInd=startXInd,yInd=startYInd; xInd>=startBuffXInd && yInd>=0 && yInd>maxYInds[xInd]; xInd--,yInd--) {
				double curVal = autoRupProb.get(xInd, yInd);
				if (curVal != 0 && !Double.isNaN(curVal))
					break;
				autoRupProb.set(xInd, yInd, FLAG_REVERT_REG_DIST);
			}
		}
//		// finally to the top left
//		for (int startXInd=startBuffXInd; --startXInd>=0;) {
//			int startYInd = highestYInd + (startBuffXInd-startXInd);
//			for (int xInd=startXInd,yInd=startYInd; ++xInd<jointProbDataset.getNumX() && ++yInd<jointProbDataset.getNumY();)
//				cumCondProbDataset.set(xInd, yInd, 1d);
//			for (int yInd=startYInd; yInd<jointProbDataset.getNumY(); yInd++)
//				cumCondProbDataset.set(startXInd, yInd, 1d);
//		}
		
		
//		double sumY = hist.calcSumOfY_Vals();
//		
//		condProbFunc = new HistogramFunction(hist.getMinX(), hist.getNum(), hist.getDelta());
//		runningTots = new double[hist.getNum()];
//		double running = 0;
//		for (int i=0; i<hist.getNum(); i++) {
//			double binVal = hist.getY(i)/sumY;
//			condProbFunc.set(i, binVal);
//			running += binVal;
//			runningTots[i] = running;
//		}
//		showDist(jointProbDataset, jointProbDataset.getMaxZ());
//		showDist(condProbDataset, 1d);
//		showDist(cumCondProbDataset, 1d);
		
//		System.out.println("**********************");
//		System.out.println("EVENTS(0,0): "+hists[0].getY(0));
//		MinMaxAveTracker crossZeroTrack = new MinMaxAveTracker();
//		EvenlyDiscretizedFunc crossZeroHist = new EvenlyDiscretizedFunc(refHist.getMinX(), refHist.getNum(), refHist.getDelta());
//		for (int i=0; i<hists.length; i++) {
//			HistogramFunction hist = hists[i];
//			crossZeroTrack.addValue(hist.getY(0));
//			crossZeroHist.set(i, hist.getY(0));
//		}
//		List<EvenlyDiscretizedFunc> elems = Lists.newArrayList();
//		elems.add(crossZeroHist);
//		new GraphWindow(elems, "Cross Zero Actual Hist");
//		EvenlyDiscretizedFunc crossZeroJointHist = new EvenlyDiscretizedFunc(refHist.getMinX(), refHist.getNum(), refHist.getDelta());
//		EvenlyDiscretizedFunc crossZeroCondHist = new EvenlyDiscretizedFunc(refHist.getMinX(), refHist.getNum(), refHist.getDelta());
//		for (int xInd=0; xInd<jointProbDataset.getNumX(); xInd++) {
//			double jointVal = jointProbDataset.get(xInd, 0);
//			if (jointVal == 0 || Double.isNaN(jointVal))
//				continue;
//			crossZeroJointHist.set(xInd, jointProbDataset.get(xInd, 0));
//			crossZeroCondHist.set(xInd, cumCondProbDataset.get(xInd, 0));
//		}
//		elems = Lists.newArrayList();
//		elems.add(crossZeroJointHist);
//		elems.add(crossZeroCondHist);
//		new GraphWindow(elems, "Cross Zero Actual Hist");
//		System.out.println("Cross 0 stats: "+crossZeroTrack);
//		System.out.println("Joint Prob(0,0): "+jointProbDataset.get(0, 0));
//		System.out.println("Cond Prob(0,0): "+condProbDataset.get(0, 0));
//		System.out.println("Cum Cond Prob(0,0): "+cumCondProbDataset.get(0, 0));
//		System.out.println("**********************");
		
		if (flipForFilling) {
			// flip back
			EvenlyDiscrXYZ_DataSet temp =
					new EvenlyDiscrXYZ_DataSet(num, num, discr_vals[0], discr_vals[0], distDeltaYears);
			for (int xInd=0; xInd<num; xInd++)
				for (int yInd=0; yInd<num; yInd++)
					temp.set(xInd, yInd, autoRupProb.get(yInd, xInd));
			autoRupProb = temp;
		}
		
		showDist(autoRupProb, 0d, 1d, "P("+follower.getName()+" | state)");
		showDist(autoGivenCrossRupProb, 0d, 1d, "P("+follower.getName()+" | state & "+driver.getName()+")");
		
		fakeStartPrevDriverTime = -Math.random()*driverActualDist.getReturnPeriod();
		fakeStartPrevFollowerTime = -Math.random()*followerActualDist.getReturnPeriod();
		hitDataset = new EvenlyDiscrXYZ_DataSet(num, num, discr_vals[0], discr_vals[0], distDeltaYears);
		
//		try {
//			Thread.sleep(1000000000000000l);
//		} catch (InterruptedException e1) {
//			// TODO Auto-generated catch block
//			e1.printStackTrace();
//		}
	}
	
	private static EvenlyDiscrXYZ_DataSet getSmoothedDataset(
			EvenlyDiscrXYZ_DataSet orig, EvenlyDiscrXYZ_DataSet threshBasis, double smoothThreshold) {
		EvenlyDiscrXYZ_DataSet smooth = new EvenlyDiscrXYZ_DataSet(orig.getNumX(), orig.getNumY(),
				orig.getMinX(), orig.getMinY(), orig.getGridSpacingX(), orig.getGridSpacingY());
		
		EvenlyDiscrXYZ_DataSet distances = new EvenlyDiscrXYZ_DataSet(5, 5, -2, -2, 1d);
		for (int i=0; i<distances.size(); i++) {
			Point2D pt = distances.getPoint(i);
			double x = pt.getX();
			double y = pt.getY();
			double dist = Math.sqrt(x*x+y*y);
			distances.set(i, dist);
		}
		
		for (int xInd=0; xInd<orig.getNumX(); xInd++) {
			for (int yInd=0; yInd<orig.getNumY(); yInd++) {
				double val = orig.get(xInd, yInd);
				if (Double.isNaN(val) || val == 0)
					continue;
				double threshVal = threshBasis.get(xInd, yInd);
				// skip x=0 and y=0, those are special cases
				if (threshVal < smoothThreshold && xInd > 0 && yInd > 0) {
//					double keepAmount = val*0.4;
					double keepAmount = val*(threshVal/smoothThreshold);
					smooth.set(xInd, yInd, keepAmount);
					double spreadAmount = val-keepAmount;
					// spread out according to distance,
					EvenlyDiscrXYZ_DataSet spreadData = new EvenlyDiscrXYZ_DataSet(5, 5, -2, -2, 1d);
					double sumZ = 0;
					for (int xAdd=-2; xAdd<=2; xAdd++) {
						for (int yAdd=-2; yAdd<=2; yAdd++) {
							if (Math.abs(xAdd) == 2d && Math.abs(yAdd) == 2d)
								continue;
							if (xAdd == 0 && yAdd == 0)
								continue;
							int myX = xInd+xAdd;
							int myY = yInd+yAdd;
							if (myX < 1 || myX >= orig.getNumX() || myY < 1 || myY >= orig.getNumY())
								continue;
							double scaledVal = distances.get(xAdd+2, yAdd+2)*spreadAmount;
							spreadData.set(xAdd+2, yAdd+2, scaledVal);
							sumZ += scaledVal;
						}
					}
					spreadData.scale(spreadAmount/sumZ);
					double amountSpreaded = 0;
					for (int i=0; i<spreadData.size(); i++) {
						Point2D pt = spreadData.getPoint(i);
						double spreadVal = spreadData.get(i);
						amountSpreaded += spreadVal;
						if (spreadVal > 0) {
							int myX = xInd + (int)pt.getX();
							int myY = yInd + (int)pt.getY();
							smooth.set(myX, myY, spreadVal+smooth.get(myX, myY));
						}
					}
					Preconditions.checkState((float)spreadAmount == (float)amountSpreaded,
							"spread mismatch: "+(float)spreadAmount+" != "+(float)amountSpreaded);
				} else {
					smooth.set(xInd, yInd, val+smooth.get(xInd, yInd));
				}
			}
		}
		
		return smooth;
	}
	
	private void showDist(EvenlyDiscrXYZ_DataSet dataset, double cptMin, double cptMax, String title) {
		if (!plotsEnabled)
			return;
		CPT cpt;
		try {
			cpt = GMT_CPT_Files.MAX_SPECTRUM.instance();
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		cpt = cpt.rescale(cptMin, cptMax);
		cpt.setBelowMinColor(Color.WHITE);
		XYZPlotSpec spec = new XYZPlotSpec(dataset, cpt, title,
				"Auto ("+follower.getName()+") years", "Cross ("+driver.getName()+") years", "z value");
		new XYZPlotWindow(spec);
		if (!"sadf".isEmpty())
			return;
//		List<HistogramFunction> funcs = Lists.newArrayList();
//		List<PlotCurveCh>
		
		double totPixels = 800;
		double funcNum = refHist.size();
		
		int pixelsPerCell = (int)(totPixels / funcNum);
		if (pixelsPerCell < 1)
			pixelsPerCell = 1;
		int sizeX = pixelsPerCell * dataset.getNumX();
		int sizeY = pixelsPerCell * dataset.getNumY();
		BufferedImage img = new BufferedImage(sizeX, sizeY, BufferedImage.TYPE_4BYTE_ABGR);
		
		for (int x=0; x<img.getWidth(); x++) {
			for (int y=0; y<img.getHeight(); y++) {
				// need to flip y
				int myY = img.getHeight()-1-y;
				
				int r = x/pixelsPerCell;
				int c = myY/pixelsPerCell;
				
//				System.out.println("("+x+","+y+") => ("+r+","+c+")");
				
				float val = (float)dataset.get(r, c);
				
//				Color color = cpt.getColor((float)matrix.getEntry(r, c));
//				Color color = cpt.getColor((float)condProbDataset.get(r, c));
				Color color = cpt.getColor(val);
				if (val == 0)
					color = Color.GRAY;
				
				img.setRGB(x, y, color.getRGB());
			}
		}
		
		JFrame frame = new JFrame();
		frame.setTitle("Dist");
		JLabel lblimage = new JLabel(new ImageIcon(img));
		frame.getContentPane().add(lblimage, BorderLayout.CENTER);
		frame.setSize(img.getWidth()+50, img.getHeight()+50);
//		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setVisible(true);
	}
	
	@Override
	public double getReturnPeriod() {
		return followerActualDist.getReturnPeriod();
	}
	
	private enum EventTimeType {
		RELATIVE_CENTERED,
		WINDOW_CENTERED,
		BASIS_RANDOM;
	}

	@Override
	public PossibleRupture getPossibleRupture(List<SimulatorEvent> prevEvents,
			double windowStart, double windowEnd) {
		EventTimeType timeType = EventTimeType.WINDOW_CENTERED;
		tot_count++;
		// find previous driver event
		Preconditions.checkState((float)(windowEnd-windowStart)==(float)getPreferredWindowLength());
		
		boolean crossHappensInWindow = false;
		
		double prevDriverTime = Double.NaN;
		double prevFollowerTime = Double.NaN;
		for (int i=prevEvents.size(); --i>=0;) {
			Preconditions.checkState(Double.isNaN(prevFollowerTime) || Double.isNaN(prevDriverTime));
			SimulatorEvent e = prevEvents.get(i);
			double eTime = e.getTimeInYears();
			Preconditions.checkState(!Double.isNaN(eTime));
			if (Double.isNaN(prevDriverTime) && driver.isMatch(e)) {
				Preconditions.checkState(!follower.isMatch(e),
						"Coruptures still exist in randomized catalog! follower="+follower.getName()+
						", driver="+driver.getName());
				
				if (!crossHappensInWindow && eTime <= windowEnd && eTime >= windowStart) {
					// this means that driver happened during this window. don't store as prev driver
					crossHappensInWindow = true;
				} else {
					prevDriverTime = eTime;
					if (!Double.isNaN(prevFollowerTime))
						// we already have a follower time, done
						break;
				}
			}
			if (Double.isNaN(prevFollowerTime) && follower.isMatch(e)) {
				Preconditions.checkState(!driver.isMatch(e), "Coruptures still exist in randomized catalog!");
				prevFollowerTime = eTime;
				if (!Double.isNaN(prevDriverTime))
					// we already have a driver time, done
					break;
			}
		}
		EvenlyDiscrXYZ_DataSet probDataset;
		if (crossHappensInWindow)
			probDataset = this.autoGivenCrossRupProb;
		else
			probDataset = this.autoRupProb;
		// this is for triggering the first event
//		if (Double.isNaN(prevDriverTime))
//			prevDriverTime = fakeStartPrevDriverTime;
		if (Double.isNaN(prevFollowerTime))
			prevFollowerTime = fakeStartPrevFollowerTime;
//		curTime -= 0.5*getPreferredWindowLength();
		double autoDelta = windowEnd - prevFollowerTime;
		if (Double.isNaN(prevDriverTime) && autoDelta < maxVal) {
			fallback_count++;
			return buildRup(followerIndepCumDist.getY(autoDelta), timeType, windowStart,
					prevFollowerTime+followerIndepCumDist.getX(followerIndepCumDist.getClosestXIndex(autoDelta)), windowStart);
		}
		double crossDelta = windowEnd - prevDriverTime;
//		Preconditions.checkState(autoDelta < maxVal, "Stuck at time="+curTime+", prevFollower="+prevFollowerTime
//				+", autoDelta="+autoDelta+", prevDriver="+prevDriverTime+", crossDelta="+crossDelta);
		double rupWindowBasis;
		if (prevDriverTime < windowStart)
			rupWindowBasis = windowStart;
		else
			rupWindowBasis = prevDriverTime;
		if (autoDelta < maxVal && crossDelta < maxVal) {
			int autoIndex = refHist.getClosestXIndex(autoDelta);
			int crossIndex = refHist.getClosestXIndex(crossDelta);
			double val = probDataset.get(autoIndex, crossIndex);
			
			int numAwayFromNonZero;
			if (val > 0) {
				numAwayFromNonZero = 0;
			} else {
				// now calculate based on distance from last good cell
				int tempAuto = autoIndex;
				int tempCross = crossIndex;
				numAwayFromNonZero = 1;
				while (true) {
					tempAuto--;
					tempCross--;
					if (tempAuto < 0 || tempCross < 0)
						break;
					double testVal = probDataset.get(tempAuto, tempCross);
					if (testVal > 0)
						break;
					numAwayFromNonZero++;
				}
			}
			
			if (val == FLAG_REVERT_REG_DIST || Double.isNaN(val)) {
				fallback_count++;
				double indepVal = followerIndepCumDist.getY(autoDelta);
				
//				// if we're 5 away from any data
//				if (numAwayFromNonZero > 5)
//					numAwayFromNonZero = 5;
//				double distProb = 0.5d+(double)numAwayFromNonZero/10d;
//				distProb = 0d;
//				if (distProb > indepVal)
//					indepVal = distProb;
				return buildRup(indepVal, timeType, windowStart,
						prevFollowerTime+followerIndepCumDist.getX(followerIndepCumDist.getClosestXIndex(autoDelta)), rupWindowBasis);
			}
//			double tots = totData.get(autoIndex, crossIndex);
//			double totThresh = 50d;
//			if (tots < totThresh) {
//				double waterlevel = ((totThresh-tots)/totThresh)*followerIndepCumDist.getY(autoDelta);
////				if (waterlevel < 0.005)
////					waterlevel = 0.005;
//				if (numAwayFromNonZero > 1) {
//					if (numAwayFromNonZero > 10)
//						numAwayFromNonZero = 10;
//					double possibleAuto = (double)numAwayFromNonZero/500d;
//					if (possibleAuto > waterlevel)
//						waterlevel = possibleAuto;
//				}
////				waterlevel = 0d;
////				if (totData.get(autoIndex, crossIndex) == 0)
////					waterlevel = 0.1;
//				if (val < waterlevel)
//					val = waterlevel;
//			}
			dep_count++;
			hitDataset.set(autoIndex, crossIndex, hitDataset.get(autoIndex, crossIndex)+1d);
			return buildRup(val, timeType, windowStart, prevDriverTime+refHist.getX(crossIndex), rupWindowBasis);
		}
		// it's been a long time, wait for a driver event then force corupture
		if (autoDelta < maxVal) {
			fallback_count++;
			return buildRup(followerIndepCumDist.getY(autoDelta),timeType, windowStart,
					prevFollowerTime+followerIndepCumDist.getX(followerIndepCumDist.getClosestXIndex(autoDelta)), rupWindowBasis);
		}
		if (crossDelta < 4*getPreferredWindowLength())
			return buildRup(1d, timeType, windowStart, 0.5d*(windowStart+windowEnd), rupWindowBasis);
		return buildRup(1d, timeType, windowStart, 0.5d*(windowStart+windowEnd), rupWindowBasis);
//		double timeDelta = curTime-prevDriverEvent.getTimeInYears();
//		if (timeDelta<maxVal)
//			return condProbFunc.getY(timeDelta);
//		return 0d;
	}
	
	private PossibleRupture buildRup(double prob, EventTimeType timeType, double windowStart, double relativeCenter, double refBasis) {
		double time;
		switch (timeType) {
		case BASIS_RANDOM:
			time = refBasis + Math.random()*getPreferredWindowLength();
			break;
		case WINDOW_CENTERED:
			time = windowStart + 0.5d*getPreferredWindowLength();
			break;
		case RELATIVE_CENTERED:
			time = relativeCenter;
			break;

		default:
			throw new IllegalStateException();
		}
		return new PossibleRupture(prob, time);
	}

	@Override
	public double getPreferredWindowLength() {
		return autoRupProb.getGridSpacingX();
//		return condProbFunc.getDelta();
	}
	
	public void printStats() {
		System.out.println(follower.getName()+" driven by "+driver.getName());
		System.out.println("Dep: "+dep_count+"/"+tot_count+" ("+(float)(100d*(double)dep_count/(double)tot_count)+" %)");
		System.out.println("Fallback: "+fallback_count+"/"+tot_count+" ("+(float)(100d*(double)fallback_count/(double)tot_count)+" %)");
		long other = tot_count - dep_count - fallback_count;
		System.out.println("Other: "+other+"/"+tot_count+" ("+(float)(100d*(double)other/(double)tot_count)+" %)");
		hitDataset.scale(1d/hitDataset.getMaxZ());
//		showDist(hitDataset, 1d);
	}
	
	private void writeDistPDF(File file, List<SimulatorEvent> randomizedCat) throws IOException {
		String followerName = follower.getName();
		String driverName = driver.getName();
		List<XYZPlotSpec> specs = Lists.newArrayList();
		List<Range> xRanges = Lists.newArrayList();
		List<Range> yRanges = Lists.newArrayList();
		
		xRanges.add(new Range(0d, maxVal));
		
		CPT maxSpect = GMT_CPT_Files.MAX_SPECTRUM.instance();
		maxSpect.setBelowMinColor(Color.WHITE);
		
		CPT probCPT = maxSpect.rescale(0d, 1d);
		CPT hitsCPT = maxSpect.rescale(0d, autoHitsData.getMaxZ());
		CPT totCPT = maxSpect.rescale(0d, totalStateData.getMaxZ());
		
		String xAxisLabel = "Years since prev "+followerName;
		String yAxisLabel = "Years since prev "+driverName;
		String title = followerName+" dists ("+driverName+" driver)";
		
		specs.add(new XYZPlotSpec(autoHitsData, hitsCPT, title, xAxisLabel, yAxisLabel, ""));
		yRanges.add(new Range(0d, 750));
		specs.add(new XYZPlotSpec(totalStateData, totCPT, title, xAxisLabel, yAxisLabel, ""));
		yRanges.add(new Range(0d, 750));
		specs.add(new XYZPlotSpec(autoRupProb, probCPT, title, xAxisLabel, yAxisLabel, ""));
		yRanges.add(new Range(0d, 750));
		specs.add(new XYZPlotSpec(autoGivenCrossHitsData, probCPT, title, xAxisLabel, yAxisLabel, ""));
		yRanges.add(new Range(0d, 750));
		
		HistogramFunction driverHist = new HistogramFunction(refHist.getMinX(), refHist.size(), refHist.getDelta());
		for (double rp : driverActualDist.getRPs())
			if (rp <= maxVal)
				driverHist.add(rp, 1d);
		HistogramFunction followerHist = new HistogramFunction(refHist.getMinX(), refHist.size(), refHist.getDelta());
		for (double rp : followerActualDist.getRPs())
			if (rp <= maxVal)
				followerHist.add(rp, 1d);
		GraphPanel intereventGP = new GraphPanel(PlotPreferences.getDefault());
		List<PlotElement> intereventElems = Lists.newArrayList();
		intereventElems.add(driverHist);
		intereventElems.add(followerHist);
		List<PlotCurveCharacterstics> intereventChars = Lists.newArrayList(
				new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK),
				new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		if (randomizedCat != null) {
			HistogramFunction randFollowHist = new HistogramFunction(refHist.getMinX(), refHist.size(), refHist.getDelta());
			List<SimulatorEvent> randMatches = follower.getMatches(randomizedCat);
			for (double rp : PeriodicityPlotter.getRPs(randMatches))
				if (rp <= maxVal)
					randFollowHist.add(rp, 1d);
			intereventElems.add(randFollowHist);
			intereventChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLUE));
			
			System.out.println(follower.getName()+": rand events: "+randMatches.size()
					+", orig events: "+(followerActualDist.getRPs().length+1));
		}
		PlotSpec intereventSpec = new PlotSpec(intereventElems, intereventChars, title, xAxisLabel, yAxisLabel);
		intereventGP.drawGraphPanel(intereventSpec, false, false);
		List<XYPlot> extraPlots = Lists.newArrayList(intereventGP.getPlot());
		
		XYZGraphPanel xyzGP = new XYZGraphPanel();
		
		xyzGP.drawPlot(specs, false, false, xRanges, yRanges, extraPlots);
		xyzGP.getChartPanel().setSize(1000, 2500);
		xyzGP.saveAsPDF(file.getAbsolutePath());
	}
	
	public static void main(String[] args) throws IOException {
		File dir = new File("/home/kevin/Simulators");
		File geomFile = new File(dir, "ALLCAL2_1-7-11_Geometry.dat");
		System.out.println("Loading geometry...");
		General_EQSIM_Tools tools = new General_EQSIM_Tools(geomFile);
//		File eventFile = new File(dir, "eqs.ALLCAL2_RSQSim_sigma0.5-5_b=0.015.barall");
		File eventFile = new File(dir, "eqs.ALLCAL2_RSQSim_sigma0.5-5_b=0.015.long.barall");
		System.out.println("Loading events...");
		List<? extends SimulatorEvent> events = EQSIMv06FileReader.readEventsFile(eventFile, tools.getElementsList());
		
		boolean include_all = true;
		
		boolean include_sj = false;
		boolean include_coachella = true;
		boolean include_mojave = true;
		
		boolean recover_debug = false;
//		RandomDistType randDistType = RandomDistType.MOJAVE_DRIVER;
//		RandomDistType randDistType = RandomDistType.MOJAVE_COACHELLA_CODRIVER;
//		RandomDistType randDistType = RandomDistType.SAN_JACINTO_COACHELLA_CODRIVER;
		RandomDistType randDistType = RandomDistType.STATE_BASED;
		
		File writeDir = new File(dir, "period_plots");
		if (!writeDir.exists())
			writeDir.mkdir();
		
		List<RuptureIdentifier> rupIdens = Lists.newArrayList();
		List<Color> colors = Lists.newArrayList();
		
		if (include_all) {
			rupIdens.add(new ElementMagRangeDescription("SAF Cholame 7+",
					ElementMagRangeDescription.SAF_CHOLAME_ELEMENT_ID, 7d, 10d));
			colors.add(Color.RED);

			rupIdens.add(new ElementMagRangeDescription("SAF Carrizo 7+",
					ElementMagRangeDescription.SAF_CARRIZO_ELEMENT_ID, 7d, 10d));
			colors.add(Color.BLUE);

			rupIdens.add(new ElementMagRangeDescription("Garlock 7+",
					ElementMagRangeDescription.GARLOCK_WEST_ELEMENT_ID, 7d, 10d));
			colors.add(Color.GREEN);
		}
		
		if (include_all || include_mojave) {
			rupIdens.add(new ElementMagRangeDescription("SAF Mojave 7+",
					ElementMagRangeDescription.SAF_MOJAVE_ELEMENT_ID, 7d, 10d));
			colors.add(Color.BLACK);
		}
		
		if (include_all || include_coachella) {
			rupIdens.add(new ElementMagRangeDescription("SAF Coachella 7+",
					ElementMagRangeDescription.SAF_COACHELLA_ELEMENT_ID, 7d, 10d));
			colors.add(Color.RED);
		}
		
		if (include_all || include_sj) {
			rupIdens.add(new ElementMagRangeDescription("San Jacinto 7+",
					ElementMagRangeDescription.SAN_JACINTO__ELEMENT_ID, 7d, 10d));
			colors.add(Color.CYAN);
		}
		
		writeDir = new File("/tmp");
		
		if (recover_debug) {
			plotsEnabled = false;
			List<? extends SimulatorEvent> rand_events = RandomCatalogBuilder.getRandomResampledCatalog(
					events, rupIdens, randDistType, true);
			
			RandomCatalogBuilder.getRandomResampledCatalog(
					rand_events, rupIdens, randDistType, true);
			plotsEnabled = true;
		} else {
			writeDir = new File(writeDir, randDistType.getFNameAdd()+"_corr_plots");
			if (!writeDir.exists())
				writeDir.mkdir();
			
			Map<IDPairing, HistogramFunction[]> origFuncs =
					PeriodicityPlotter.plotACDF_CCDFs(writeDir, events, rupIdens,
							null, null, 2000d, 10d);
			PeriodicityPlotter.	plotACDF_CCDFs(writeDir, events, rupIdens,
					randDistType, origFuncs, 2000d, 10d);
			
//			File subDir = new File(writeDir, "round2");
//			if (!subDir.exists())
//				subDir.mkdir();
//			PeriodicityPlotter.	plotTimeBetweenAllIdens(subDir, rand_events, rupIdens, rupIdenNames, colors,
//					RandomDistType.MOJAVE_DRIVER, origFuncs, 2000d, 10d);
			
			System.out.println("DONE");
		}
		
		File pdfDir = new File("/tmp/follower_pdfs", randDistType.getFNameAdd());
		if (!pdfDir.exists())
			pdfDir.mkdir();
		
		List<File> pdfs = Lists.newArrayList();
		
		plotsEnabled = false;
		List<SimulatorEvent> randomizedCat = RandomCatalogBuilder.getRandomResampledCatalog(events, rupIdens, randDistType, true);
		plotsEnabled = true;
		
		for (int i=0; i<rupIdens.size(); i++) {
			RuptureIdentifier iden = rupIdens.get(i);
			String name = iden.getName();
			
			RandomReturnPeriodProvider origProv = randDistType.instance(iden, events);
			if (!(origProv instanceof FollowerReturnPeriodProvider))
				continue;
			
			FollowerReturnPeriodProvider prov = (FollowerReturnPeriodProvider)origProv;
			
			String fName = PeriodicityPlotter.getFileSafeString(name);
			File file = new File(pdfDir, fName+".pdf");
			
			prov.writeDistPDF(file, randomizedCat);
			
			pdfs.add(file);
		}
		
		if (!pdfs.isEmpty())
			PeriodicityPlotter.combinePDFs(pdfs, new File(pdfDir, "follower_dists.pdf"));
		
//		System.exit(0);
	}

}
