package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.Location;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.FaultTrace;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.utils.DeformationModelFetcher;
import scratch.UCERF3.utils.UCERF3_DataUtils;

public class WeldonTableWrite {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		List<FaultSectionPrefData> subSects = new DeformationModelFetcher(
				FaultModels.FM3_1, DeformationModels.GEOLOGIC,
				UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, 0.1).getSubSectionList();
		
		CSVFile<String> inCSV = CSVFile.readFile(new File("/tmp/UCERF3_OpenIntervalv2.csv"), true);
		CSVFile<String> csv = new CSVFile<String>(true);
		
		// header
		List<String> header1 = Lists.newArrayList(inCSV.getLine(0));
		header1.add("");
		header1.add("");
		header1.add("");
		header1.add("");
		csv.addLine(header1);
		List<String> header2 = Lists.newArrayList(inCSV.getLine(1));
		header2.add("Lat1");
		header2.add("Lon1");
		header2.add("Lat2");
		header2.add("Lon2");
		csv.addLine(header2);
		
		for (int i=2; i<inCSV.getNumRows(); i++) {
			List<String> line = Lists.newArrayList(inCSV.getLine(i));
			// find subsection
			String parentName = line.get(0).replaceAll("-", " ");
			int subsectNum = Integer.parseInt(line.get(1));
			
			FaultSectionPrefData matchingSect = null;
			
			int cntInParent = -1;
			for (FaultSectionPrefData sect : subSects) {
				if (sect.getParentSectionName().equals(parentName)) {
					cntInParent++;
					if (cntInParent == subsectNum) {
						matchingSect = sect;
						break;
					}
				}
			}
			Preconditions.checkNotNull(matchingSect, "COULDN'T FIND: "+parentName+" (subsect: "+subsectNum+")");
			
			FaultTrace trace = matchingSect.getFaultTrace();
			Location loc1 = trace.get(0);
			Location loc2 = trace.get(trace.size()-1);
			line.add(loc1.getLatitude()+"");
			line.add(loc1.getLongitude()+"");
			line.add(loc2.getLatitude()+"");
			line.add(loc2.getLongitude()+"");
		}
		
		csv.writeToFile(new File("/tmp/UCERF3_OpenIntervalv2_withLocs.csv"));
	}

}
