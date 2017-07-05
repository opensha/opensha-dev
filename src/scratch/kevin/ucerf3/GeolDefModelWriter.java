package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.Location;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import com.google.common.base.Joiner;
import com.google.common.collect.Lists;

import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.utils.DeformationModelFileParser;
import scratch.UCERF3.utils.UCERF3_DataUtils;

public class GeolDefModelWriter {
	
	public static void main(String[] args) throws IOException {
		FaultModels fm = FaultModels.FM3_2;
//		File inputFile = new File("D:\\Documents\\temp\\Geologic_Def_model_3_1_2012_08_14.csv");
		File inputFile = new File("/tmp/Geologic_Def_model_3_2_2012_09_06.csv");
		CSVFile<String> inputCSV = CSVFile.readFile(inputFile, true);
		
		Map<Integer, FaultSectionPrefData> fmSects = fm.fetchFaultSectionsMap();
		
		CSVFile<String> dmCSV = new CSVFile<String>(true);
		CSVFile<String> minCSV = new CSVFile<String>(true);
		CSVFile<String> maxCSV = new CSVFile<String>(true);
		
		for (int row=1; row<inputCSV.getNumRows(); row++) {
			List<String> inputLine = inputCSV.getLine(row);
			
//			System.out.println(Joiner.on(",").join(inputLine));
			
			if (inputLine.get(1).isEmpty())
				continue;
			
			String miniStr = inputLine.get(1);
			int[] mini = DeformationModelFileParser.parseMinisectionNumber(miniStr);
			
			if (!fmSects.containsKey(mini[0]))
				continue;
			double rate, min, max;
			
			if (inputLine.get(2).trim().isEmpty()) {
				rate = Double.NaN;
				min = Double.NaN;
				max = Double.NaN;
			} else {
				rate = Double.parseDouble(inputLine.get(2));
				if (inputLine.get(3).trim().isEmpty()) {
					min = Double.NaN;
					max = Double.NaN;
				} else {
					min = Double.parseDouble(inputLine.get(3));
					max = Double.parseDouble(inputLine.get(4));
				}
			}
			double rake = Double.parseDouble(inputLine.get(5));

			FaultSectionPrefData sect = fmSects.get(mini[0]);
			Location startLoc = sect.getFaultTrace().get(mini[1]-1);
			Location endLoc = sect.getFaultTrace().get(mini[1]);

			miniStr = DeformationModelFileParser.getMinisectionString(mini);

			List<String> linePrefix = Lists.newArrayList(miniStr, (float)startLoc.getLongitude()+"",
					(float)startLoc.getLatitude()+"", (float)endLoc.getLongitude()+"",
					(float)endLoc.getLatitude()+"");

			List<String> dmLine = Lists.newArrayList(linePrefix);
			dmLine.add((float)rate+"");
			dmLine.add((float)rake+"");
			List<String> minLine = Lists.newArrayList(linePrefix);
			minLine.add((float)min+"");
			minLine.add((float)rake+"");
			List<String> maxLine = Lists.newArrayList(linePrefix);
			maxLine.add((float)max+"");
			maxLine.add((float)rake+"");

			dmCSV.addLine(dmLine);
			minCSV.addLine(minLine);
			maxCSV.addLine(maxLine);
		}
		
		File outputDir = new File(UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR.getParentFile(), "DeformationModels");
		
		dmCSV.writeToFile(new File(outputDir, "geologic_slip_rake_fm_3_2_2012_09_06.csv"));
		minCSV.writeToFile(new File(outputDir, "geologic_slip_rake_fm_3_2_lowerbound_2012_09_06.csv"));
		maxCSV.writeToFile(new File(outputDir, "geologic_slip_rake_fm_3_2_upperbound_2012_09_06.csv"));
	}

}
