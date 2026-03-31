package scratch.kevin.nshm27;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Map;

import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionScalingRelationships;
import org.opensha.sha.util.TectonicRegionType;

import gov.usgs.earthquake.nshmp.erf.nshm27.NSHM26_InvConfigFactory;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM26_InterfaceFaultModels;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM26_InterfaceMinSubSects;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM26_LogicTree;

public class InversionMinSubSectMagTests {

	public static void main(String[] args) throws IOException {
		NSHM26_InvConfigFactory factory = new NSHM26_InvConfigFactory();
		Map<NSHM26_InterfaceFaultModels, double[]> countMags = new HashMap<>();
		NSHM26_InterfaceMinSubSects[] minSects = NSHM26_InterfaceMinSubSects.values();
		for (NSHM26_InterfaceFaultModels fm : NSHM26_InterfaceFaultModels.values()) {
			LogicTreeBranch<LogicTreeNode> branch = NSHM26_LogicTree.buildDefault(
					fm.getSeisReg(), TectonicRegionType.SUBDUCTION_INTERFACE, false);
			branch.setValue(PRVI25_SubductionScalingRelationships.AVERAGE);
			FaultSystemRupSet rupSet = factory.buildRuptureSet(branch, FaultSysTools.defaultNumThreads());
			
			int[] counts = new int[minSects.length];
			double[] magSums = new double[minSects.length];
			
			for (int r=0; r<rupSet.getNumRuptures(); r++) {
				double mag = rupSet.getMagForRup(r);
				int numSects = rupSet.getSectionsIndicesForRup(r).size();
				for (int m=0; m<minSects.length; m++) {
					if (numSects == minSects[m].getValue()) {
						counts[m]++;
						magSums[m] += mag;
					}
				}
			}
			double[] avgMags = new double[minSects.length];
			for (int m=0; m<magSums.length; m++)
				avgMags[m] = magSums[m]/(double)counts[m];
			countMags.put(fm, avgMags);
		}
		
		DecimalFormat mDF = new DecimalFormat("0.0");
		
		double[] avgMags = new double[minSects.length];
		for (NSHM26_InterfaceFaultModels fm : countMags.keySet()) {
			double[] mags = countMags.get(fm);
			System.out.println(fm.getName());
			for (int m=0; m<minSects.length; m++) {
				System.out.println("\t"+minSects[m].getName()+":\t"+mDF.format(mags[m]));
				avgMags[m] += mags[m]/(double)countMags.keySet().size();
			}
		}
		System.out.println("Average");
		for (int m=0; m<minSects.length; m++)
			System.out.println("\t"+minSects[m].getName()+":\t"+mDF.format(avgMags[m]));
	}

}
