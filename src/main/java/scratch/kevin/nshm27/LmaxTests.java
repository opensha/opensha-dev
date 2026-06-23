package scratch.kevin.nshm27;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.BitSet;

import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionScalingRelationships;
import org.opensha.sha.util.TectonicRegionType;

import gov.usgs.earthquake.nshmp.erf.nshm27.NSHM27_InvConfigFactory;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_LogicTree;
import gov.usgs.earthquake.nshmp.erf.nshm27.util.NSHM27_RegionLoader.NSHM27_SeismicityRegions;

public class LmaxTests {

	public static void main(String[] args) throws IOException {
		NSHM27_SeismicityRegions reg = NSHM27_SeismicityRegions.AMSAM;
//		NSHM27_SeismicityRegions reg = NSHM27_SeismicityRegions.GNMI;
		
		double[] lMaxs = {150, 175, 200, 250, 300, 400, 500, 600, 700, 800, 900, 950, 1000, 1100, 1200, 1300};
		
		PRVI25_SubductionScalingRelationships scale = PRVI25_SubductionScalingRelationships.LOGA_C4p0;
		
		LogicTreeBranch<LogicTreeNode> branch = NSHM27_LogicTree.buildDefault(reg, TectonicRegionType.SUBDUCTION_INTERFACE, false);
		branch.setValue(scale);
		FaultSystemRupSet rupSet = new NSHM27_InvConfigFactory().buildRuptureSet(branch, 16);
		
		DecimalFormat magDF = new DecimalFormat("0.00");
		
		for (double lMax : lMaxs) {
			BitSet rups = new BitSet(rupSet.getNumRuptures());
			for (int r=0; r<rupSet.getNumRuptures(); r++)
				if (rupSet.getLengthForRup(r)*1e-3 <= lMax)
					rups.set(r);
			MinMaxAveTracker sectMmax = new MinMaxAveTracker();
			for (int s=0; s<rupSet.getNumSections(); s++) {
				double mMax = 0d;
				for (int rupIndex : rupSet.getRupturesForSection(s))
					if (rups.get(rupIndex))
						mMax = Math.max(mMax, rupSet.getMagForRup(rupIndex));
				sectMmax.addValue(mMax);
			}
			System.out.println("\t"+(int)lMax+" km:  \tM"+magDF.format(sectMmax.getAverage()));
//					+"\t["+magDF.format(sectMmax.getMin())+", "+magDF.format(sectMmax.getMax())+"]");
		}
	}

}
