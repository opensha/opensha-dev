package scratch.kevin.ucerf3;

import java.util.ArrayList;
import java.util.List;

import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import com.google.common.base.Joiner;
import com.google.common.collect.Lists;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.utils.DeformationModelFetcher;
import scratch.UCERF3.utils.UCERF3_DataUtils;

public class ImperialCreepPrint {
	
	public static void main(String[] args) {
		int imperialParentID = 97;
		
		FaultModels fm = FaultModels.FM3_1;
		
		for (DeformationModels dm : DeformationModels.forFaultModel(fm)) {
			if (dm.getRelativeWeight(InversionModels.CHAR_CONSTRAINED) == 0)
				continue;
			ArrayList<FaultSectionPrefData> subSects =
					new DeformationModelFetcher(fm, dm, UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, 0.1).getSubSectionList();
			
			List<Double> aseis = Lists.newArrayList();
			for (FaultSectionPrefData sect : subSects) {
				if (sect.getParentSectionId() == imperialParentID) {
					double myAseis = sect.getAseismicSlipFactor();
					if (!aseis.contains(myAseis))
						aseis.add(myAseis);
				}
			}
			System.out.println(dm+": "+Joiner.on(",").join(aseis));
		}
	}

}
