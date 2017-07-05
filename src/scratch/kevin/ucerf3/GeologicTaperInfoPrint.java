package scratch.kevin.ucerf3;

import java.util.ArrayList;

import org.opensha.commons.geo.Location;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.utils.DeformationModelFetcher;
import scratch.UCERF3.utils.UCERF3_DataUtils;

public class GeologicTaperInfoPrint {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		FaultModels fm = FaultModels.FM3_1;
		DeformationModels dm = DeformationModels.GEOLOGIC;
		
		ArrayList<FaultSectionPrefData> datas =
				new DeformationModelFetcher(fm, dm, UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, 0.1).getSubSectionList();
		
		System.out.println("Name\tOrig Slip Rate\tReduced Slip Rate\tCoupling Coeff\tAseismicity Factor");
		for (FaultSectionPrefData data : datas) {
			if (data.getParentSectionId() == 97 || data.getParentSectionId() == 172) {
				System.out.println(data.getName()+"\t"+data.getOrigAveSlipRate()+"\t"+data.getReducedAveSlipRate()
						+"\t"+data.getCouplingCoeff()+"\t"+data.getAseismicSlipFactor());
			}
		}
		
		for (FaultModels fm1 : FaultModels.values()) {
//			dm = DeformationModels.forFaultModel(fm1).get(0);
//			datas =
//					new DeformationModelFetcher(fm1, dm, UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, 0.1).getSubSectionList();
			datas = fm1.fetchFaultSections();
			System.out.println("\n"+fm1);
			for (FaultSectionPrefData data : datas) {
				if (data.getSectionId() == 13) {
					for (Location loc : data.getFaultTrace())
						System.out.println("\t"+loc.getLatitude()+"\t"+loc.getLongitude());
				}
			}
		}
	}

}
