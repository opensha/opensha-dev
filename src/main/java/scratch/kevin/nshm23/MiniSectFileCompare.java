package scratch.kevin.nshm23;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;

import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels.MinisectionSlipRecord;

public class MiniSectFileCompare {

	public static void main(String[] args) throws FileNotFoundException, IOException {
//		File refFile = new File("/home/kevin/workspace/opensha/src/main/resources/data/erf/nshm23/def_models/geodetic/"
//				+ "fm_v2/SHEN_BIRD-include_ghost_corr.txt");
//		File compFile = new File("/tmp/nfeg04_02_sig.slp3");
		
//		File refFile = new File("/home/kevin/workspace/opensha/src/main/resources/data/erf/nshm23/def_models/geodetic/"
//				+ "fm_v2/SHEN_BIRD-no_ghost_corr.txt");
//		File compFile = new File("/tmp/nfeg04_02ng_sig.slp3");
		
//		File refFile = new File("/home/kevin/workspace/opensha/src/main/resources/data/erf/nshm23/def_models/geodetic/"
//				+ "fm_v2/EVANS-include_ghost_corr.txt");
//		File compFile = new File("/tmp/Output_withGT_suite_011123.txt");
		
		File refFile = new File("/home/kevin/workspace/opensha/src/main/resources/data/erf/nshm23/def_models/geodetic/"
				+ "fm_v2/EVANS-no_ghost_corr.txt");
		File compFile = new File("/tmp/Output_noGT_suite_011123.txt");
		
		Map<Integer, List<MinisectionSlipRecord>> refMinis = NSHM23_DeformationModels.loadGeodeticModel(
				new InputStreamReader(new FileInputStream(refFile)));
		Map<Integer, List<MinisectionSlipRecord>> compMinis = NSHM23_DeformationModels.loadGeodeticModel(
				new InputStreamReader(new FileInputStream(compFile)));
		
		int compFaultMissing = 0;
		int compSlipChanges = 0;
		int compRakeChanges = 0;
		int compStdDevChanges = 0;
		int compLocMismatches = 0;
		int compMiniCountMismatches = 0;
		int numRef = 0;
		
		for (Integer faultID : refMinis.keySet()) {
			List<MinisectionSlipRecord> refRecs = refMinis.get(faultID);
			List<MinisectionSlipRecord> compRecs = compMinis.get(faultID);
			
			numRef+= refRecs.size();
			
			if (compRecs == null) {
				System.out.println("Comparison is missing fault: "+faultID);
				compFaultMissing++;
				continue;
			}
			
			if (refRecs.size() != compRecs.size()) {
				System.out.println("Comparison has "+compRecs.size()+" minis for "+faultID+", reference has "+refRecs.size());
				int numDiff = refRecs.size() - compRecs.size();
				if (numDiff < 0)
					numDiff = -numDiff;
				compMiniCountMismatches += numDiff;
			}
			
			for (int i=0; i<refRecs.size() && i<compRecs.size(); i++) {
				MinisectionSlipRecord refRec = refRecs.get(i);
				MinisectionSlipRecord compRec = compRecs.get(i);
				if (!LocationUtils.areSimilar(refRec.startLoc, compRec.startLoc) || !LocationUtils.areSimilar(refRec.endLoc, compRec.endLoc)) {
					System.out.println("Location mismatch for "+faultID+" minisection "+refRec.minisectionID);
					System.out.println("\tStart loc: "+refRec.startLoc+" -> "+compRec.startLoc+" (dist="
							+(float)LocationUtils.horzDistanceFast(refRec.startLoc, compRec.startLoc)+")");
					System.out.println("\tEnd loc: "+refRec.endLoc+" -> "+compRec.endLoc+" (dist="
							+(float)LocationUtils.horzDistanceFast(refRec.endLoc, compRec.endLoc)+")");
					compLocMismatches++;
				}
				
				if ((float)refRec.slipRate != (float)compRec.slipRate)
					compSlipChanges++;
				if ((float)refRec.slipRateStdDev != (float)compRec.slipRateStdDev)
					compStdDevChanges++;
				if ((float)refRec.rake != (float)compRec.rake)
					compRakeChanges++;
			}
		}
		
		int compFaultAdded = 0;
		for (Integer faultID : compMinis.keySet()) {
			if (!refMinis.containsKey(faultID)) {
				System.out.println("Comparison has additional fault: "+faultID);
				compFaultAdded++;
			}
		}
		
		System.out.println("Number of reference faults: "+refMinis.size());
		System.out.println("Number of reference minisection: "+numRef);
		System.out.println("Slip rate changes: "+compSlipChanges);
		System.out.println("Std Dev changes: "+compStdDevChanges);
		System.out.println("Rake changes: "+compRakeChanges);
		System.out.println("Comp faults missing: "+compFaultMissing);
		System.out.println("Comp faults added: "+compFaultAdded);
		System.out.println("Comp location mismatches: "+compLocMismatches);
		System.out.println("Comp minisection count mismatches: "+compMiniCountMismatches);
	}

}
