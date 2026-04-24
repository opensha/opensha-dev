package scratch.kevin.nshm23;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.List;

import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.faultSurface.FaultSection;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;

public class FaultSlipRateInvestigations {
	
	public static void main(String[] args) throws IOException {
		NSHM23_FaultModels fm = NSHM23_FaultModels.WUS_FM_v3;
		FaultModels u3FM = FaultModels.FM3_1;
		
//		int nshm23_id = FaultSectionUtils.findSectionID(fm.getFaultSections(), "Silver", "Creek");
//		int u3_id = FaultSectionUtils.findSectionID(u3FM.getFaultSections(), "Silver", "Creek");
		
		int nshm23_id = FaultSectionUtils.findSectionID(fm.getFaultSections(), "Rose", "Canyon");
		int u3_id = FaultSectionUtils.findSectionID(u3FM.getFaultSections(), "Rose", "Canyon");
		
		double u3Avg = avgSlipRate(DeformationModels.MEAN_UCERF3.build(u3FM, DeformationModels.MEAN_UCERF3, null), u3_id);
		NSHM23_DeformationModels[] dms = NSHM23_DeformationModels.values();
		double[] slips = new double[dms.length];
		for (int d=0; d<dms.length; d++)
			slips[d] = avgSlipRate(dms[d].build(fm), nshm23_id);
		
		System.out.println("UCERF3 Avg:\t"+(float)u3Avg);
		DecimalFormat pDF = new DecimalFormat("0.0%");
		for (int d=0; d<dms.length; d++) {
			double fract = (slips[d]-u3Avg)/u3Avg;
			System.out.print("NSHM23 "+dms[d].getShortName()+":\t"+(float)slips[d]+"\t(");
			if (fract >= 0)
				System.out.print("+");
			System.out.println(pDF.format(fract)+" vs U3)");
		}
	}
	
	private static double avgSlipRate(List<? extends FaultSection> subSects, int parentID) {
		return subSects.stream().filter(S->S.getParentSectionId() == parentID)
				.mapToDouble(S->S.getOrigAveSlipRate()).average().getAsDouble();
	}

}
