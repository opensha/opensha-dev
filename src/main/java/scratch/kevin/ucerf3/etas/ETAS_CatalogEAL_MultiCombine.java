package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.UncertainArbDiscDataset;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Table;

public class ETAS_CatalogEAL_MultiCombine {

	public static void main(String[] args) throws IOException {
		File fullTD_csvFile = null;
		File noERT_csvFile = null;
		File outputDir = null;
		String outputPrefix = null;
		if (args.length == 1 && args[0].equals("--hardcoded")) {
			File baseDir = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations");
			File lossCombDir = new File(baseDir, "losses_combined");
			Preconditions.checkState(lossCombDir.exists() || lossCombDir.mkdir());
			
			File fullTD_dir = new File(baseDir,
					"2016_02_19-mojave_m7-10yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14-combined100k/loss_results/combined");
			File noERT_dir = new File(baseDir,
					"2016_02_22-mojave_m7-10yr-no_ert-subSeisSupraNucl-gridSeisCorr-combined100k/loss_results/combined");
			
			outputDir = new File(lossCombDir, "mojave_100k_both_models");
			Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
			outputPrefix = "mojave_combined";
			
			fullTD_csvFile = new File(fullTD_dir, "gmpes_combined_exceed.csv");
			noERT_csvFile = new File(noERT_dir, "gmpes_combined_exceed.csv");
		} else if (args.length != 4) {
			System.err.println("USAGE: <FullTD-exceed-CSV> <NoERT-exceed-CSV> <output-dir> <output-prefix>");
			System.exit(2);
		} else {
			fullTD_csvFile = new File(args[0]);
			Preconditions.checkState(fullTD_csvFile.exists(), "Doesn't exist: %s", fullTD_csvFile.getAbsolutePath());
			noERT_csvFile = new File(args[1]);
			Preconditions.checkState(noERT_csvFile.exists(), "Doesn't exist: %s", noERT_csvFile.getAbsolutePath());
			outputDir = new File(args[2]);
			Preconditions.checkState(outputDir.exists() || outputDir.mkdir(), "Doesn't exist: %s", outputDir.getAbsolutePath());
			outputPrefix = args[3];
		}
		
		boolean isLogX = false;
		boolean triggeredOnly = false;
		double maxX = 200;
		
		CSVFile<String> fullTD_csv = CSVFile.readFile(fullTD_csvFile, true);
		fullTD_csv = checkFixCSV(fullTD_csv, fullTD_csvFile);
		CSVFile<String> noERT_csv = CSVFile.readFile(noERT_csvFile, true);
		noERT_csv = checkFixCSV(noERT_csv, noERT_csvFile);
		
		if (fullTD_csvFile.getAbsolutePath().contains("coachella") || fullTD_csvFile.getAbsolutePath().contains("bernardino"))
			maxX = 50;
		
		Table<String, Double, DiscretizedFunc> exceedFuncs = HashBasedTable.create();
		
		Map<Double, DiscretizedFunc> fullTD = loadCSV(fullTD_csv);
		Map<Double, DiscretizedFunc> noERT = loadCSV(noERT_csv);
		
		for (Double duration : fullTD.keySet())
			exceedFuncs.put("FullTD", duration, fullTD.get(duration));
		for (Double duration : noERT.keySet())
			exceedFuncs.put("NoERT", duration, noERT.get(duration));
		
		String xAxisLabel = fullTD_csv.get(0, 0);
		ETAS_CatalogEALCalculator.writeLossExceed(outputDir, outputPrefix, exceedFuncs, null, isLogX, triggeredOnly, xAxisLabel, maxX, true, false);
	}
	
	private static int realLineSize(List<String> line) {
		int size = 0;
		for (String val : line) {
			if (val == null || val.isEmpty())
				break;
			size++;
		}
		return size;
	}
	
	private static CSVFile<String> checkFixCSV(CSVFile<String> csv, File csvFile) throws IOException {
		int size1 = realLineSize(csv.getLine(0));
		int size2 = realLineSize(csv.getLine(1));
		System.out.println("Sizes: "+size1+" "+size2);
		if (size1 == 1+2*(size2-1)) {
			System.out.println("Fixing offset in: "+csvFile.getAbsolutePath());
			// needs fixing
			CSVFile<String> fixed = new CSVFile<String>(true);
			List<String> header = csv.getLine(0).subList(0, size2);
			fixed.addLine(header);
			for (int row=1; row<csv.getNumRows(); row++) {
				List<String> line = Lists.newArrayList();
				line.add(csv.get(row, 0)); // x value
				// get data from above
				int col;
				if (row == 1) {
					col = size2;
				} else {
					col = 1;
				}
				while (line.size() < size2)
					line.add(csv.get(row-1, col++));
				fixed.addLine(line);
			}
			csv.writeToFile(new File(csvFile.getAbsolutePath()+".bak"));
			fixed.writeToFile(csvFile);
			return fixed;
		}
		return csv;
	}
	
	private static Map<Double, DiscretizedFunc> loadCSV(CSVFile<String> csv) {
		Map<Double, DiscretizedFunc> hists = Maps.newHashMap();
		
		double min = Double.parseDouble(csv.get(1, 0));
		double max = Double.parseDouble(csv.get(csv.getNumRows()-1, 0));
		int num = csv.getNumRows()-1;
		for (double duration : ETAS_CatalogEALCalculator.durations) {
			String durStr = ETAS_CatalogEALCalculator.getDurationLabel(duration);
			int matchCol = -1;
			for (int col = 1; col<csv.getNumCols(); col++) {
				if (csv.get(0, col).equals(durStr)) {
					matchCol = col;
					break;
				}
			}
			if (matchCol < 0) {
				System.out.println("No match for duration: "+durStr);
				continue;
			}
			int lowerCol = -1;
			int upperCol = -1;
			for (int i=matchCol+1; i<csv.getNumCols(); i++) {
				if (lowerCol >= 0 && upperCol >= 0)
					break;
				String name = csv.get(0, i);
				if (name.equals(ETAS_CatalogEALCalculator.conf_lower_str)) {
					Preconditions.checkState(lowerCol < 0, "Duplicate lower col before both found!");
					lowerCol = i;
				} else if (name.equals(ETAS_CatalogEALCalculator.conf_upper_str)) {
					Preconditions.checkState(upperCol < 0, "Duplicate upper col before both found!");
					upperCol = i;
				}
				Preconditions.checkState(!name.contains(" yr") && !name.contains(" mo") && !name.contains(" wk")
						&& !name.contains(" day"), "Encountered new duration before conf bounds!");
			}
			Preconditions.checkState(lowerCol >= 0 && upperCol > lowerCol, "exceed doesn't have conf bounds");
			
			HistogramFunction hist = new HistogramFunction(min, max, num);
			HistogramFunction lower = new HistogramFunction(min, max, num);
			HistogramFunction upper = new HistogramFunction(min, max, num);
			for (int i=0; i<num; i++) {
				int row = i+1;
				double y = Double.parseDouble(csv.get(row, matchCol));
				double yLow = Double.parseDouble(csv.get(row, lowerCol));
				double yUp = Double.parseDouble(csv.get(row, upperCol));
				hist.set(i, getCleaned(y));
				lower.set(i, getCleaned(yLow));
				upper.set(i, getCleaned(yUp));
			}
//			hists.put(duration, hist);
			hists.put(duration, new UncertainArbDiscDataset(hist, lower, upper));
		}
		return hists;
	}
	
	private static double getCleaned(double y) {
		if ((float)y == 0f)
			return 0d;
		if (y < 0 || y < 1e-12) {
			Preconditions.checkState(y > -1e-10);
			return 0d;
		}
		return y;
	}

}
