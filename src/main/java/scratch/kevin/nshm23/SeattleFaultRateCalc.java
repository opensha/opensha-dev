package scratch.kevin.nshm23;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;

import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupMFDsModule;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

public class SeattleFaultRateCalc {

	public static void main(String[] args) throws IOException {
		double[] minMags = {0d, 6d, 6.5d, 7d, 7.5d};
		
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged.zip"));
		FaultSystemRupSet rupSet = sol.getRupSet();
		
//		int[] parentIDs = {
//				FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "seattle", "middle"),
//				FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "seattle", "east"),
//				FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "seattle", "north"),
//				FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "seattle", "south"),
//		};
		
		int[] parentIDs = {
				FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "Turner", "Mill", "Creek"),
				FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "Mount", "Angel"),
				FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "Gales", "Chehalem"),
				FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "Gales", "Parsons"),
				FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "Beaverton"),
				FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "Canby"),
				FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "Bolton"),
				FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "Sylvain", "Oatfield"),
				FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "Portland", "Hills", "(south)"),
				FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "Portland", "Hills", "(north)"),
		};
		
		for (int i=0; i<parentIDs.length; i++)
			Preconditions.checkState(parentIDs[i] >= 0, "No parent ID at position %s", i);
		
		HashSet<Integer> rups = new HashSet<>();
		for (int parentID : parentIDs) {
			String name = null;
			for (FaultSection sect : rupSet.getFaultSectionDataList()) {
				if (sect.getParentSectionId() == parentID) {
					name = sect.getParentSectionName();
					break;
				}
			}
			List<Integer> myRups = rupSet.getRupturesForParentSection(parentID);
			Preconditions.checkNotNull(myRups, "No ruptures found for %s. %s", parentID, name);
			int prevSize = rups.size();
			rups.addAll(myRups);
			int numAdded = rups.size() - prevSize;
			int numDups = myRups.size() - numAdded;
			System.out.println("Found "+myRups.size()+" rups for "+name+" (now have "+rups.size()+" total unique; added "+numAdded+", skipped "+numDups+"  duplicates)");
		}
		
		System.out.println("Found "+rups.size()+" total rups");
		
		RupMFDsModule mfds = sol.requireModule(RupMFDsModule.class);
		for (double minMag : minMags) {
			double rate = 0d;
			for (int rupIndex : rups) {
				DiscretizedFunc mfd = mfds.getRuptureMFD(rupIndex);
				if (mfd == null)
					mfd = new LightFixedXFunc(new double[] {rupSet.getMagForRup(rupIndex)}, new double[] {sol.getRateForRup(rupIndex)});
				for (Point2D pt : mfd)
					if (pt.getX() >= minMag)
						rate += pt.getY();
			}
			if (minMag > 0d)
				System.out.println("M>"+(float)minMag+":\t"+(float)rate+" /year\t(RI: "+(float)(1d/rate)+" years)");
			else
				System.out.println("Supra-seis:\t"+(float)rate+" /year\t(RI: "+(float)(1d/rate)+" years)");
				
		}
	}

}
