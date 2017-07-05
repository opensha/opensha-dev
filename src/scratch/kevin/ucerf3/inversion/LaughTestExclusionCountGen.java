package scratch.kevin.ucerf3.inversion;

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.commons.gui.plot.GraphWindow;

import com.google.common.collect.Lists;

import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.inversion.InversionFaultSystemRupSet;
import scratch.UCERF3.inversion.InversionFaultSystemRupSetFactory;
import scratch.UCERF3.inversion.coulomb.CoulombRatesTester;
import scratch.UCERF3.inversion.laughTest.LaughTestFilter;
import scratch.UCERF3.utils.IDPairing;

public class LaughTestExclusionCountGen {
	
	public static int[] getRupCount(LaughTestFilter filter) {
		InversionFaultSystemRupSet rupSet = InversionFaultSystemRupSetFactory.forBranch(
				filter, 0.1, FaultModels.FM3_1);
		int numPossibleConnections = 0;
		List<FaultSectionPrefData> sects = rupSet.getFaultSectionDataList();
		for (int sect1=0; sect1<rupSet.getNumSections(); sect1++) {
			int parent1 = sects.get(sect1).getParentSectionId();
			for (int sect2 : rupSet.getCloseSectionsList(sect1)) {
				if (sect2 < sect1)
					continue;
				if (parent1 != sects.get(sect2).getParentSectionId())
					numPossibleConnections++;
			}
		}
		HashSet<IDPairing> junctions = new HashSet<IDPairing>();
		for (int rupIndex=0; rupIndex<rupSet.getNumRuptures(); rupIndex++) {
			List<FaultSectionPrefData> rupSects = rupSet.getFaultSectionDataForRupture(rupIndex);
			for (int i=1; i<rupSects.size(); i++) {
				FaultSectionPrefData sect1 = rupSects.get(i-1);
				FaultSectionPrefData sect2 = rupSects.get(i);
				if (sect1.getParentSectionId() != sect2.getParentSectionId()) {
					if (sect1.getSectionId() < sect2.getSectionId())
						junctions.add(new IDPairing(sect1.getSectionId(), sect2.getSectionId()));
					else
						junctions.add(new IDPairing(sect2.getSectionId(), sect1.getSectionId()));
				}
			}
		}
		int numConnections = junctions.size();
		int num = rupSet.getNumRuptures();
		System.out.println(num+" ruptures, "+numConnections+"/"+numPossibleConnections+" connections");
		rupSet = null;
		sects = null;
		junctions = null;
		System.gc();
		int[] ret = { num, numConnections, numPossibleConnections };
		return ret;
	}
	
	

	/**
	 * @param args
	 */
	public static void main(String[] args) {
//		LaughTestFilter filter = LaughTestFilter.getUCERF3p2Filter();
		LaughTestFilter filter = LaughTestFilter.getDefault();
		
		int[] origRups = getRupCount(filter);
		
//		List<Double> jumpDists = Lists.newArrayList(5d, 3d, 1d, 0.1d, 1e-5);
		List<Double> jumpDists = Lists.newArrayList();
		List<int[]> jumpDistCounts = Lists.newArrayList();
		double orig = filter.getMaxJumpDist();
		for (double jumpDist : jumpDists) {
			filter.setMaxJumpDist(jumpDist);
			jumpDistCounts.add(getRupCount(filter));
			filter.clearLaughTests();
		}
		filter.setMaxJumpDist(orig);
		filter.clearLaughTests();
		
		orig = filter.getMaxAzimuthChange();
		filter.setMaxAzimuthChange(Double.POSITIVE_INFINITY);
		int[] azChangeRups = getRupCount(filter);
		filter.setMaxAzimuthChange(orig);
		filter.clearLaughTests();
		
		orig = filter.getMaxTotAzimuthChange();
		filter.setMaxTotAzimuthChange(Double.POSITIVE_INFINITY);
		int[] totAzChangeRups = getRupCount(filter);
		filter.setMaxTotAzimuthChange(orig);
		filter.clearLaughTests();
		
		orig = filter.getMaxCmlAzimuthChange();
		filter.setMaxCmlAzimuthChange(Double.POSITIVE_INFINITY);
		int[] totCmlAzChangeRups = getRupCount(filter);
		filter.setMaxCmlAzimuthChange(orig);
		filter.clearLaughTests();
		
		orig = filter.getMaxCmlRakeChange();
		filter.setMaxCmlRakeChange(Double.POSITIVE_INFINITY);
		int[] totCmlRakeChangeRups = getRupCount(filter);
		filter.setMaxCmlRakeChange(orig);
		filter.clearLaughTests();
		
		CoulombRatesTester origCoulomb = filter.getCoulombFilter();
		filter.setCoulombFilter(null);
		int[] totCoulombRups = getRupCount(filter);
		filter.setCoulombFilter(origCoulomb);
		filter.clearLaughTests();
		
		// now test some coulomb values
		double[] mins = { 0.02, 0.04, 0.06, 0.08, 0.1 };
		double[] excls = { 1.0, 1.25, 1.5, Double.POSITIVE_INFINITY };
		List<String> strs = Lists.newArrayList();
		List<int[]> counts = Lists.newArrayList();
		for (double excl : excls) {
			DiscretizedFunc rupCountsFunc = new ArbitrarilyDiscretizedFunc();
			DiscretizedFunc junctionsFunc = new ArbitrarilyDiscretizedFunc();
			for (double min : mins) {
				origCoulomb.setMinAverageProb(min);
				origCoulomb.setMinIndividualProb(min);
				origCoulomb.setMinimumStressExclusionCeiling(excl);
				
				strs.add("PDCFF: "+(float)min+"; DCFF: "+(float)excl);
				int[] count = getRupCount(filter);
				
				double fractRupsExcluded = (double)(totCoulombRups[0] - count[0])
						/ (double)totCoulombRups[0];
				double fractJunctionsExcluded = (double)(totCoulombRups[1] - count[1])
						/ (double)totCoulombRups[1];
				
				rupCountsFunc.set(min, fractRupsExcluded);
				junctionsFunc.set(min, fractJunctionsExcluded);
				
				counts.add(count);
				filter.clearLaughTests();
			}
			ArrayList<DiscretizedFunc> funcs = Lists.newArrayList();
			funcs.add(rupCountsFunc);
			funcs.add(junctionsFunc);
			ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList();
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
			GraphWindow gw = new GraphWindow(
					funcs, "Fraction of Junctions/Ruptures Excluded", chars);
			gw.setX_AxisLabel("PΔCFF Threshold");
			gw.setY_AxisLabel("Fraction Excluded");
		}
		
		System.out.println("Orig Num Rups: "+origRups[0]+" ("+origRups[1]+"/"+origRups[2]+" junctions)");
		System.out.println("Az filtered: "+getCountsStr(origRups, azChangeRups));
		System.out.println("Tot az filtered: "+getCountsStr(origRups, totAzChangeRups));
		System.out.println("Cumulative az filtered: "+getCountsStr(origRups, totCmlAzChangeRups));
		System.out.println("Cumulative rake filtered: "+getCountsStr(origRups, totCmlRakeChangeRups));
		System.out.println("Coulomb filtered: "+getCountsStr(origRups, totCoulombRups));
		
		System.out.println("\n*** Jump Dists ***");
		for (int i=0; i<jumpDists.size(); i++) {
			int[] jumpCounts = jumpDistCounts.get(i);
			System.out.println("dist="+jumpDists.get(i)+": "+jumpCounts[0]+" rups, "+jumpCounts[1]+"/"+jumpCounts[2]+" connections");
		}
		
		System.out.println("\n*** Coulomb tests ***");
		for (int i=0; i<strs.size(); i++)
			System.out.println("\t"+strs.get(i)+": "+getCountsStr(counts.get(i), totCoulombRups));
	}
	
	private static String getCountsStr(int[] withFilter, int[] withoutFilter) {
		return "ruptures: "+(withoutFilter[0]-withFilter[0])+" excluded ("+withoutFilter[0]+" rups w/o);"
				+"\tconnections: "+withFilter[1]+"/"+withFilter[2]+" with, "+withoutFilter[1]+"/"+withoutFilter[2]
						+" w/o ("+(withoutFilter[1]-withFilter[1])+" excluded)";
	}

}
