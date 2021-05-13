package scratch.kevin.ucerf3;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.commons.gui.plot.GraphWindow;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.utils.DeformationModelFetcher;
import scratch.UCERF3.utils.DeformationModelFileParser;
import scratch.UCERF3.utils.DeformationModelFileParser.DeformationSection;
import scratch.UCERF3.utils.UCERF3_DataUtils;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public class ABMSanGorgonioFixer {
	
	private static int[] buildRange(int from, int to) {
		int delta;
		if (from < to)
			delta = 1;
		else
			delta = -1;
		int len = (int)Math.abs(from - to) + 1;
		
		int[] ret = new int[len];
		int i = 0;
		
		System.out.println("from "+from+" to "+to);
		for (int val=from; val-delta != to; val += delta, i++) {
			System.out.println("ret["+i+"] = "+val);
			ret[i] = val;
		}
		return ret;
	}
	
	private static double[] getLengths(int[] range, DeformationSection dmSect) {
		double[] lengths = new double[range.length];
		
		int cnt = 0;
		for (int mini : range) {
			Location loc1 = dmSect.getLocs1().get(mini);
			Location loc2 = dmSect.getLocs2().get(mini);
			
			lengths[cnt++] = LocationUtils.horzDistance(loc1, loc2);
		}
		
		return lengths;
	}
	
	private static double[] getSlips(int[] range, DeformationSection dmSect) {
		double[] slips = new double[range.length];
		
		int cnt = 0;
		for (int mini : range)
			slips[cnt++] = dmSect.getSlips().get(mini);
		
		return slips;
	}
	
	private static double[] getDistsAlong(double[] lengths) {
		double[] ret = new double[lengths.length+1];
		
		double total = 0d;
		for (int i=0; i<lengths.length; i++) {
			total += lengths[i];
			ret[i+1] = total;
		}
		
		return ret;
	}
	
	private static ArbitrarilyDiscretizedFunc buildFunc(double[] distAlongs, double[] slips, double offset) {
		ArbitrarilyDiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
		
		for (int i=0; i<slips.length; i++) {
			double start = offset+distAlongs[i];
			double end = offset+distAlongs[i+1];
			double slip = slips[i];
			
			func.set(start, slip);
			func.set(end-0.01, slip);
		}
		
		return func;
	}
	
	public static EvenlyDiscretizedFunc buildSumFunc(List<DiscretizedFunc> funcs, double max) {
		EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(0d, max, 200);
		
		for (int i=0; i<func.size(); i++) {
			double x = func.getX(i);
			
			double sum = 0;
			for (DiscretizedFunc otherFunc : funcs) {
				if (otherFunc.getMinX() > x || otherFunc.getMaxX() < x)
					continue;
				for (int index=1; index<otherFunc.size(); index++) {
					double myStart = otherFunc.getX(index-1);
					double myEnd = otherFunc.getX(index);
					if (x >= myStart && x <= myEnd) {
						sum += otherFunc.getY(index);
						break;
					}
				}
			}
			
			func.set(i, sum);
		}
		
		return func;
	}

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		int northBranchID = 294;
		double slipToLeave = 2; // leave 2 mm/yr on NBMC
		
		File dataDir = new File(UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR.getParentFile(), "DeformationModels");
		File inputFile = new File(dataDir, "ABM_slip_rake_fm_3_2_2012_08_31.csv");
		File outputFile = new File(dataDir, "ABM_slip_rake_fm_3_2_2012_09_25.csv");
		
		List<int[]> northBranchRanges = Lists.newArrayList();
		List<int[]> newSectRanges = Lists.newArrayList();
		List<Integer> newSectIDs = Lists.newArrayList();
		
		// first mapping is the start of NBMC (0) to San Bernardino N (3)
		newSectIDs.add(282);
		northBranchRanges.add(buildRange(0, 0));
		newSectRanges.add(buildRange(3, 3));
		
		// now map NBMC (1-6) to San Bernardino S (0-4)
		newSectIDs.add(283);
		northBranchRanges.add(buildRange(1, 6));
		newSectRanges.add(buildRange(0, 4));
		
		// now map NBMC (7-13) to San Gorgonio Pass Garnet Hill (7-0)
		newSectIDs.add(284);
		northBranchRanges.add(buildRange(7, 13));
		newSectRanges.add(buildRange(7, 0));
		
		Map<Integer, DeformationSection> dmSects = DeformationModelFileParser.load(inputFile);
		
		DeformationSection northBranchSect = dmSects.get(northBranchID);
		
		Color nbColor = Color.BLUE;
		Color sumColor = Color.BLACK;
		Color[] newSectColors = { Color.RED, Color.GREEN, Color.ORANGE };
		
		int[] nbTotRange = buildRange(0, 13);
		double[] origTotLengths = getLengths(nbTotRange, northBranchSect);
		double[] origTotSlips = getSlips(nbTotRange, northBranchSect);
		double[] origNewTotSlips = new double[origTotSlips.length];
		for (int i=0; i<origNewTotSlips.length; i++)
			origNewTotSlips[i] = slipToLeave;
		double[] origTotDistsAlong = getDistsAlong(origTotLengths);
		
		ArbitrarilyDiscretizedFunc origNBFunc = buildFunc(origTotDistsAlong, origTotSlips, 0d);
		ArbitrarilyDiscretizedFunc newNBFunc = buildFunc(origTotDistsAlong, origNewTotSlips, 0d);
		ArbitrarilyDiscretizedFunc origTot = new ArbitrarilyDiscretizedFunc();
		ArbitrarilyDiscretizedFunc newTot = new ArbitrarilyDiscretizedFunc();
		
		ArrayList<DiscretizedFunc> beforeFuncs = Lists.newArrayList();
		beforeFuncs.add(origNBFunc);
		ArrayList<DiscretizedFunc> afterFuncs = Lists.newArrayList();
		afterFuncs.add(newNBFunc);
		ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList(
				new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, PlotSymbol.FILLED_CIRCLE, 4f, nbColor));
		
		for (int i=0; i<newSectIDs.size(); i++) {
			int newSectID = newSectIDs.get(i);
			int[] northBranchRange = northBranchRanges.get(i);
			int[] newRange = newSectRanges.get(i);
			
			DeformationSection newSect = dmSects.get(newSectID);
			
			double offset = origTotDistsAlong[northBranchRange[0]];
			double[] origLengths = getLengths(northBranchRange, northBranchSect);
			double[] origSlips = getSlips(northBranchRange, northBranchSect);
			double origTotal = StatUtils.sum(origLengths);
			
			double[] newLengths = getLengths(newRange, newSect);
			double[] newSlips = getSlips(newRange, newSect);
			double newTotal = StatUtils.sum(newLengths);
			
			double lengthDelta = Math.abs(origTotal-newTotal);
			System.out.println("Doing mapping "+i+" ("+origTotal+" km => "+newTotal
					+" km, delta = "+lengthDelta+" km)");
			
			// now scale the original lengths to match the new total length
			double lengthScaleFactor = origTotal / newTotal;
			for (int j=0; j<newLengths.length; j++)
				newLengths[j] = newLengths[j] * lengthScaleFactor;
			newTotal = StatUtils.sum(newLengths);
			Preconditions.checkState(Math.abs(newTotal-origTotal) < 1e-10,
					"scaling didn't work! ("+origTotal+" != "+newTotal+")");
			
			double[] origDistsAlong = getDistsAlong(origLengths);
			double[] newDistsAlong = getDistsAlong(newLengths);
			
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, PlotSymbol.FILLED_CIRCLE, 4f, newSectColors[i]));
			beforeFuncs.add(buildFunc(newDistsAlong, newSlips, offset));
			
			for (int j=0; j<newRange.length; j++) {
				double newStart = newDistsAlong[j];
				double newEnd = newDistsAlong[j+1];
				
				List<Double> origSlipParts = Lists.newArrayList();
				List<Double> origLengthParts = Lists.newArrayList();
				
				System.out.println("\tMapping to new sect "+newRange[j]);
				
				// find sections that overlap with the new section
				double runningTotalLenghts = 0;
				for (int origIndex=0; origIndex<origSlips.length; origIndex++) {
					double myStart = origDistsAlong[origIndex];
					double myEnd = origDistsAlong[origIndex+1];
					
					if (myEnd <= newStart)
						// this is completely before our section
						continue;
					if (myStart > newEnd)
						// this is completely after our section
						break;
					
					// this means that there is at least partial overlap
					if (myStart < newStart)
						myStart = newStart;
					if (myEnd > newEnd)
						myEnd = newEnd;
					
					double myNewLength = myEnd - myStart;
					runningTotalLenghts += myNewLength;
					
					float percentage = (float)(100d * (myNewLength/origLengths[origIndex]));
					
					System.out.println("\t\tOrig sect "+northBranchRange[origIndex]+"\t"
							+percentage+"%\tslip: "+origSlips[origIndex]);
					
					origSlipParts.add(origSlips[origIndex]);
					origLengthParts.add(myNewLength);
				}
				
				double slipToMap = DeformationModelFetcher.calcLengthBasedAverage(origLengthParts, origSlipParts);
				
				slipToMap -= slipToLeave;
				Preconditions.checkState(slipToMap > 0, "Mapping negative slip!");
				
				double origSlip = newSlips[j];
				newSlips[j] = origSlip + slipToMap;
				System.out.println("\t\tMapping "+slipToMap+" mm/yr, new slip rate: "+newSlips[j]);
				newSect.getSlips().set(newRange[j], newSlips[j]);
			}
			
			afterFuncs.add(buildFunc(newDistsAlong, newSlips, offset));
		}
		
		// now set NBMC to the minimum
		List<Double> slips = northBranchSect.getSlips();
		for (int i=0; i<slips.size(); i++)
			slips.set(i, slipToLeave);
		
		DeformationModelFileParser.write(dmSects, outputFile);
		
		double xMax = origTotDistsAlong[origTotDistsAlong.length-1];
		beforeFuncs.add(buildSumFunc(beforeFuncs, xMax));
		afterFuncs.add(buildSumFunc(afterFuncs, xMax));
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, sumColor));
		
		ArrayList<DiscretizedFunc> allFuncs = Lists.newArrayList();
		allFuncs.addAll(beforeFuncs);
		allFuncs.addAll(afterFuncs);
		ArrayList<PlotCurveCharacterstics> allChars = Lists.newArrayList();
		for (PlotCurveCharacterstics pchar : chars) {
			pchar = (PlotCurveCharacterstics)pchar.clone();
			pchar.setLineWidth(1f);
			pchar.setSymbolWidth(2f);
			allChars.add(pchar);
		}
		allChars.addAll(chars);
		
		ArrayList<PlotCurveCharacterstics> afterChars = Lists.newArrayList(chars);
		afterFuncs.add(beforeFuncs.get(beforeFuncs.size()-1));
		afterChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
		
		String xAxisLabel = "Distance Along North Branch (km)";
		String yAxisLabel = "Slip Rate (mm/yr)";
		
		GraphWindow gw = new GraphWindow(beforeFuncs, "Slip Rates (Before)", chars);
		gw.setX_AxisLabel(xAxisLabel);
		gw.setY_AxisLabel(yAxisLabel);
		gw.setAxisRange(0, xMax, 0d, 25d);
		gw = new GraphWindow(afterFuncs, "Slip Rates (After)", afterChars);
		gw.setX_AxisLabel(xAxisLabel);
		gw.setY_AxisLabel(yAxisLabel);
		gw.setAxisRange(0, xMax, 0d, 25d);
		gw = new GraphWindow(allFuncs, "Slip Rates (After)", allChars);
		gw.setX_AxisLabel(xAxisLabel);
		gw.setY_AxisLabel(yAxisLabel);
		gw.setAxisRange(0, xMax, 0d, 25d);
	}

}
