package scratch.kevin.ucerf3;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.estimate.MinMaxPrefEstimate;
import org.opensha.commons.util.FaultUtils;
import org.opensha.refFaultParamDb.dao.db.DB_AccessAPI;
import org.opensha.refFaultParamDb.dao.db.DB_ConnectionPool;
import org.opensha.refFaultParamDb.dao.db.FaultSectionVer2_DB_DAO;
import org.opensha.refFaultParamDb.gui.addEdit.faultSection.EditFaultSection;
import org.opensha.refFaultParamDb.vo.EstimateInstances;
import org.opensha.refFaultParamDb.vo.FaultSectionData;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import com.google.common.collect.Lists;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.utils.DeformationModelFileParser;
import scratch.UCERF3.utils.DeformationModelFileParser.DeformationSection;

public class GeologicRatesInsert {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		try {
//			FileWriter fw = new FileWriter(new File("/tmp/fm2_1_parents.txt"));
//			for (FaultSectionPrefData data : FaultModels.FM2_1.fetchFaultSections())
//				fw.write(data.getSectionId()+"\t"+data.getSectionName()+"\n");
//			fw.close();
//			System.exit(0);
			
			boolean validateOnly = false;
				
			DB_AccessAPI db;
//			DB_AccessAPI db = DB_ConnectionPool.getLatestReadWriteConn();
			if (!validateOnly) {
				db = DB_ConnectionPool.getDirectLatestReadWriteConnection();
				DB_ConnectionPool.authenticateDBConnection(true, false);
			} else {
				db = DB_ConnectionPool.getLatestReadOnlyConn();
			}
			
			FaultSectionVer2_DB_DAO fs2db = new FaultSectionVer2_DB_DAO(db);
			
			List<Map<Integer, DeformationSection>> dm31s = Lists.newArrayList(
					DeformationModelFileParser.load(DeformationModels.GEOLOGIC.getDataFileURL(FaultModels.FM3_1)),
					DeformationModelFileParser.load(DeformationModels.GEOLOGIC_LOWER.getDataFileURL(FaultModels.FM3_1)),
					DeformationModelFileParser.load(DeformationModels.GEOLOGIC_UPPER.getDataFileURL(FaultModels.FM3_1)));
			List<Map<Integer, DeformationSection>> dm32s = Lists.newArrayList(
					DeformationModelFileParser.load(DeformationModels.GEOLOGIC.getDataFileURL(FaultModels.FM3_2)),
					DeformationModelFileParser.load(DeformationModels.GEOLOGIC_LOWER.getDataFileURL(FaultModels.FM3_2)),
					DeformationModelFileParser.load(DeformationModels.GEOLOGIC_UPPER.getDataFileURL(FaultModels.FM3_2)));
			
			System.out.println("Getting all sects!");
			ArrayList<FaultSectionData> sects = fs2db.getAllFaultSections();
			System.out.println("done");
			
			double minProb = 0.025;
			double prefProb = 0.5;
			double maxProb = 0.975;
			for (FaultSectionData sect : sects) {
				
				DeformationSection prefSection = dm31s.get(0).get(sect.getSectionId());
				DeformationSection lowerSection = dm31s.get(1).get(sect.getSectionId());
				DeformationSection upperSection = dm31s.get(2).get(sect.getSectionId());
				
				if (prefSection == null) {
					// try fm 3.2
					prefSection = dm32s.get(0).get(sect.getSectionId());
					lowerSection = dm32s.get(1).get(sect.getSectionId());
					upperSection = dm32s.get(2).get(sect.getSectionId());
				}
				
				if (prefSection == null) {
					System.out.println("No match for "+sect.getName()+", skipping!");
					continue;
				}
				
				double pref = FaultUtils.getLengthBasedAngleAverage(prefSection.getLocsAsTrace(), prefSection.getSlips());
				double min = FaultUtils.getLengthBasedAngleAverage(lowerSection.getLocsAsTrace(), lowerSection.getSlips());
				double max = FaultUtils.getLengthBasedAngleAverage(upperSection.getLocsAsTrace(), upperSection.getSlips());
				
				if (validateOnly) {
					MinMaxPrefEstimate origEst = (MinMaxPrefEstimate)sect.getAveLongTermSlipRateEst().getEstimate();
					if ((float)pref != (float)origEst.getPreferred())
						System.out.println(sect.getName()+" pref discrep: "+pref+" != "+origEst.getPreferred());
					if ((float)min != (float)origEst.getMin())
						System.out.println(sect.getName()+" min discrep: "+min+" != "+origEst.getMin());
					if ((float)max != (float)origEst.getMax())
						System.out.println(sect.getName()+" max discrep: "+max+" != "+origEst.getMax());
					continue;
				}
				
				MinMaxPrefEstimate slipEst = new MinMaxPrefEstimate(min, max, pref, minProb, maxProb, prefProb);
				
				double rake = FaultUtils.getLengthBasedAngleAverage(prefSection.getLocsAsTrace(), prefSection.getRakes());
				
				MinMaxPrefEstimate rakeEst = new MinMaxPrefEstimate(Double.NaN, Double.NaN, rake, Double.NaN, Double.NaN, Double.NaN);
				
				sect.setAveLongTermSlipRateEst(new EstimateInstances(slipEst, EditFaultSection.SLIP_RATE_UNITS));
				
				sect.setAveRakeEst(new EstimateInstances(rakeEst, EditFaultSection.RAKE_UNITS));
				
				System.out.println("Updating: "+sect.getSectionName());
				
				fs2db.update(sect);
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.exit(1);
		}
		System.exit(0);
	}

}
