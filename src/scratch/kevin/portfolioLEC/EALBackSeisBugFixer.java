package scratch.kevin.portfolioLEC;

import java.io.File;
import java.io.IOException;

import org.opensha.commons.data.CSVFile;

import com.google.common.base.Preconditions;

public class EALBackSeisBugFixer {
	
	private static void handleFiles(File origFile, File backSeisFile) throws IOException {
		System.out.println("Fixing: "+backSeisFile.getAbsolutePath());
		
		CSVFile<String> orig = CSVFile.readFile(origFile, true);
		CSVFile<String> backSeis = CSVFile.readFile(backSeisFile, true);
		
		double origERF0_val = Double.parseDouble(orig.get(1, 1)); // this is the value that was erroneously added to each combined value
		
		for (int row=1; row<backSeis.getNumRows(); row++) {
			double combinedVal = Double.parseDouble(backSeis.get(row, 1));
			combinedVal -= origERF0_val;
			Preconditions.checkState(combinedVal > 0);
			backSeis.set(row, 1, combinedVal+"");
		}
		
		String fixedName = backSeisFile.getName();
		if (fixedName.contains(".csv"))
			fixedName = fixedName.substring(0, fixedName.indexOf(".csv"));
		fixedName += "_fixed.csv";
		File fixedFile = new File(backSeisFile.getParentFile(), fixedName);
		
		backSeis.writeToFile(fixedFile);
	}
	
	private static void handleDir(File dir) throws IOException {
		for (File file : dir.listFiles()) {
			if (file.isDirectory()) {
				handleDir(file);
				continue;
			}
			String fileName = file.getName();
			if (!fileName.endsWith(".csv"))
				continue;
			if (fileName.contains("_incl_back_seis"))
				continue;
			// now we have a candidate original file - lets see if there's a corresponding one with back seis added
			String backSeisName = fileName.substring(0, fileName.indexOf(".csv")); // remove extention
			backSeisName += "_incl_back_seis.csv";
			File backSeisFile = new File(dir, backSeisName);
			if (backSeisFile.exists())
				// we have a match!
				handleFiles(file, backSeisFile);
		}
	}

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		File dir = new File("/home/kevin/OpenSHA/portfolio_lec/parallel_eal");
		
		handleDir(dir);
	}

}
