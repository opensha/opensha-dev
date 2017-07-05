package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import com.google.common.collect.Lists;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.utils.DeformationModelFetcher;
import scratch.UCERF3.utils.UCERF3_DataUtils;

public class SubSectAreaCalc {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		FaultModels fm = FaultModels.FM3_1;
		
		File outFile = new File("D:\\Documents\\temp\\fm3_1_subsect_slips_areas.csv");
		
		DeformationModels[] dms =
				{ DeformationModels.GEOLOGIC, DeformationModels.ABM,
				DeformationModels.NEOKINEMA, DeformationModels.ZENG };
		
		double[][] areas = new double[dms.length][];
		double[][] slips = new double[dms.length][];
		
		for (int d=0; d<dms.length; d++) {
			DeformationModels dm = dms[d];
			DeformationModelFetcher fetch = new DeformationModelFetcher(
					fm, dm, UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, 0.1);
			List<FaultSectionPrefData> sects = fetch.getSubSectionList();
			double[] myAreas = new double[sects.size()];
			double[] mySlips = new double[sects.size()];
			for (int i=0; i<sects.size(); i++) {
				FaultSectionPrefData sect = sects.get(i);
				
				// aseismicity reduces area; km --> m on length & DDW
				double area = sect.getTraceLength()*1e3*sect.getReducedDownDipWidth()*1e3;
				myAreas[i] = area;
				mySlips[i] = sect.getReducedAveSlipRate();
			}
			
			areas[d] = myAreas;
			slips[d] = mySlips;
		}
		
		CSVFile<String> csv = new CSVFile<String>(true);
		
		List<String> header = Lists.newArrayList("SubSect ID");
		for (DeformationModels dm : dms) {
			header.add(dm.getName()+" Area (m^2)");
			header.add(dm.getName()+" Slip Rate (mm/yr)");
		}
		csv.addLine(header);
		
		for (int i=0; i<areas[0].length; i++) {
			List<String> line = Lists.newArrayList(i+"");
			for (int d=0; d<dms.length; d++) {
				line.add(areas[d][i]+"");
				line.add(slips[d][i]+"");
			}
			csv.addLine(line);
		}
		
		csv.writeToFile(outFile);
	}

}
