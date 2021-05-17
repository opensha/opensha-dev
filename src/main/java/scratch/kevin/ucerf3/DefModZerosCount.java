package scratch.kevin.ucerf3;

import java.io.IOException;
import java.util.Map;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.utils.DeformationModelFileParser;
import scratch.UCERF3.utils.DeformationModelFileParser.DeformationSection;

public class DefModZerosCount {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		FaultModels[] fms = { FaultModels.FM3_1, FaultModels.FM3_2 };
		for (FaultModels fm : fms) {
			for (DeformationModels dm : DeformationModels.forFaultModel(fm)) {
				Map<Integer, DeformationSection> def =
					DeformationModelFileParser.load(dm.getDataFileURL(fm));
				int zeros = 0;
				int nans = 0;
				for (DeformationSection sect : def.values()) {
					for (double slip : sect.getSlips()) {
						if (slip == 0) {
							zeros++;
							break;
						}
					}

					for (double slip : sect.getSlips()) {
						if (Double.isNaN(slip)) {
							nans++;
							break;
						}
					}
				}
				System.out.println(fm+": "+dm+": zeros="+zeros+", nans="+nans);
			}
		}
	}

}
