package scratch.kevin.prvi25;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.uncertainty.UncertaintyBoundType;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateFileLoader.RateType;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateModel;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SeismicityRateEpoch;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader.PRVI25_SeismicityRegions;

import com.google.common.base.Preconditions;

public class RateFileScale1973to1900 {

	public static void main(String[] args) throws IOException {
		File ratesDir = new File("/home/kevin/workspace/opensha/src/main/resources/data/erf/prvi25/seismicity/rates");
		File outDir = new File(ratesDir, PRVI25_CrustalSeismicityRate.RATE_DATE+"/1973_scaled_to_1900");
		Preconditions.checkState(outDir.exists() || outDir.mkdir());
		File inDir1900 = new File(ratesDir, PRVI25_CrustalSeismicityRate.RATE_DATE+"/"+PRVI25_SeismicityRateEpoch.FULL.getRateSubDirName());
		File inDir1973 = new File(ratesDir, PRVI25_CrustalSeismicityRate.RATE_DATE+"/"+PRVI25_SeismicityRateEpoch.RECENT.getRateSubDirName());
		
		RateType type = RateType.M1_TO_MMAX;
		UncertaintyBoundType boundType = UncertaintyBoundType.CONF_95;
		PRVI25_SeismicityRegions[] regions = PRVI25_SeismicityRegions.values();
		
		double sumM5_1973 = 0d;
		double sumM5_1900 = 0d;
		List<CSVFile<String>> inCSVs1973 = new ArrayList<>();
		SeismicityRateModel[] inModels1900 = new SeismicityRateModel[regions.length];
		SeismicityRateModel[] inModels1973 = new SeismicityRateModel[regions.length];
		for (int r=0; r<regions.length; r++) {
			String csvName = regions[r].name()+".csv";
			CSVFile<String> csv1973 = CSVFile.readFile(new File(inDir1973, csvName), false);
			inCSVs1973.add(csv1973);
			inModels1973[r] = new SeismicityRateModel(csv1973, type, boundType);
			CSVFile<String> csv1900 = CSVFile.readFile(new File(inDir1900, csvName), false);
			inModels1900[r] = new SeismicityRateModel(csv1900, type, boundType);
			
			sumM5_1973 += inModels1973[r].getMeanRecord().rateAboveM1;
			sumM5_1900 += inModels1900[r].getMeanRecord().rateAboveM1;
		}
		
		double rateScalar = sumM5_1900/sumM5_1973;
		System.out.println("Rate scalar:\t"+(float)sumM5_1900+" / "+(float)sumM5_1973+" = "+(float)rateScalar);
		
		for (int r=0; r<regions.length; r++) {
			PRVI25_SeismicityRegions seisReg = regions[r];
			String csvName = regions[r].name()+".csv";
			
			CSVFile<String> inCSV = inCSVs1973.get(r);
			CSVFile<String> outCSV = new CSVFile<>(false);
			boolean inMBranch = false;
			boolean inExactBranch = false;
			for (int row=0; row<inCSV.getNumRows(); row++) {
				if (inMBranch) {
					if (inCSV.get(row, 0).contains("Branches")) {
						inMBranch = false;
						outCSV.addLine(inCSV.getLine(row));
					} else {
						// scale
						Preconditions.checkState(inCSV.getLine(row).size() == 3);
						List<String> line = new ArrayList<>(3);
						line.add(inCSV.get(row, 0)); // quantile
						line.add((inCSV.getDouble(row, 1)*rateScalar) + ""); // rate
						line.add(inCSV.get(row, 2)); // b
						outCSV.addLine(line);
					}
				} else if (inExactBranch) {
					// scale
					Preconditions.checkState(inCSV.getLine(row).size() == 7);
					List<String> line = new ArrayList<>(7);
					line.add(inCSV.get(row, 0)); // mag
					for (int c=1; c<7; c++)
						line.add((inCSV.getDouble(row, c)*rateScalar) + ""); // rate
					outCSV.addLine(line);
				} else {
					if (inCSV.get(row, 0).equals("quantile") && inCSV.get(row, 1).startsWith("Rate M"))
						inMBranch = true;
					else if (inCSV.getLine(row).size() == 7 && inCSV.get(row, 0).isBlank() && inCSV.get(row, 1).equals("2"))
						inExactBranch = true;
					outCSV.addLine(inCSV.getLine(row));
				}
			}
			outCSV.writeToFile(new File(outDir, csvName));
		}
	}

}
