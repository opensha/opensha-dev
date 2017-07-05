package scratch.kevin;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.math3.stat.StatUtils;
import org.dom4j.DocumentException;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.FaultTrace;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.utils.FaultSystemIO;

public class NonZeroCalc {
	
	private static boolean isMultiFault(FaultSystemSolution sol, int id) {
		ArrayList<Integer> parents = new ArrayList<Integer>();
		ArrayList<Integer> minIDs = new ArrayList<Integer>();
		ArrayList<Integer> maxIDs = new ArrayList<Integer>();
		for (int secID : sol.getRupSet().getSectionsIndicesForRup(id)) {
			int parent = sol.getRupSet().getFaultSectionData(secID).getParentSectionId();
			int parentIndex = parents.indexOf(parent);
			if (parentIndex < 0) {
				parents.add(parent);
				minIDs.add(secID);
				maxIDs.add(secID);
			} else {
				if (secID < minIDs.get(parentIndex))
					minIDs.set(parentIndex, secID);
				if (secID > maxIDs.get(parentIndex))
					maxIDs.set(parentIndex, secID);
			}
		}
		if (parents.size() == 1)
			return false;
		
		for (int i=0; i<parents.size(); i++) {
			for (int j=0; j<parents.size(); j++) {
				if (i == j)
					continue;
				if (!isNexTo(sol.getRupSet(), minIDs.get(i), maxIDs.get(i), minIDs.get(j), maxIDs.get(j), id))
					return true;
			}
		}
		
		return false;
	}
	
	private static boolean isNexTo(FaultSystemRupSet rupSet, int minID1, int maxID1, int minID2, int maxID2, int rupID) {
		if (calcDist(rupSet.getFaultSectionData(minID1), rupSet.getFaultSectionData(minID2), rupID) < 0.1)
			return true;
		if (calcDist(rupSet.getFaultSectionData(minID1), rupSet.getFaultSectionData(maxID2), rupID) < 0.1)
			return true;
		if (calcDist(rupSet.getFaultSectionData(maxID1), rupSet.getFaultSectionData(minID2), rupID) < 0.1)
			return true;
		if (calcDist(rupSet.getFaultSectionData(maxID1), rupSet.getFaultSectionData(maxID2), rupID) < 0.1)
			return true;
		return false;
	}
	
	private static double calcDist(FaultSectionPrefData data1, FaultSectionPrefData data2, int rupID) {
		FaultTrace trace1 = data1.getFaultTrace();
		FaultTrace trace2 = data2.getFaultTrace();
		
		double[] dists = new double[4];
		
		dists[0] = LocationUtils.horzDistanceFast(trace1.get(0),				trace2.get(0));
		dists[1] = LocationUtils.horzDistanceFast(trace1.get(trace1.size()-1),	trace2.get(0));
		dists[2] = LocationUtils.horzDistanceFast(trace1.get(0),				trace2.get(trace2.size()-1));
		dists[3] = LocationUtils.horzDistanceFast(trace1.get(trace1.size()-1),	trace2.get(trace2.size()-1));
		
		float dipDirDiff = Math.abs(data1.getDipDirection() - data2.getDipDirection());
		if (data1.getAveDip() < 90 && data2.getAveDip() < 90 && dipDirDiff > 45) {
			System.out.println("rupID to look at: "+rupID);
		}
		
		return StatUtils.min(dists);
	}

	/**
	 * @param args
	 * @throws DocumentException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException, DocumentException {
		FaultSystemSolution sol = FaultSystemIO.loadSol(
				new File("dev/scratch/UCERF3/preComputedData/InversionSolutions/2011_10_19-morgan-ALLCAL_Model1.zip"));
//				new File("dev/scratch/UCERF3/preComputedData/InversionSolutions/Model2.xml"));
		int numMultFaults = 0;
		int numMultFaultNonZeros = 0;
		int numNonZeros = 0;
		
		for (int i=0; i<sol.getRupSet().getNumRuptures(); i++) {
			boolean multi = isMultiFault(sol, i);
			boolean nonZero = sol.getRateForRup(i) > 0;
			
			if (multi) {
				numMultFaults++;
				if (nonZero)
					numMultFaultNonZeros++;
			}
			
			if (nonZero)
				numNonZeros++;
		}
		
		System.out.println("Non zero: "+getStr(numNonZeros, sol.getRupSet().getNumRuptures()));
		System.out.println("Multi non zero: "+getStr(numMultFaultNonZeros, numMultFaults));
		System.out.println("Single non zero: "+getStr(numNonZeros-numMultFaultNonZeros,
				sol.getRupSet().getNumRuptures()-numMultFaults));
	}
	
	private static String getStr(int num, int tot) {
		float percent = (float)num / (float)tot * 100f;
		return num+"/"+tot+" ("+percent+" %)";
	}

}
