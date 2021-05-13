package scratch.kevin.portfolioLEC;

import java.io.File;
import java.io.IOException;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.util.DataUtils;

public class CSVCompare {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
//		File file1 = new File("/home/kevin/OpenSHA/portfolio_lec/parallel_eal/run4_compare/AS2008_incl_back_seis.csv");
//		File file2 = new File("/home/kevin/OpenSHA/portfolio_lec/parallel_eal/run4/AS2008_incl_back_seis_fixed.csv");
		
		File file1 = new File("/home/kevin/OpenSHA/portfolio_lec/parallel_eal/run4_compare/CB2008_incl_back_seis.csv");
		File file2 = new File("/home/kevin/OpenSHA/portfolio_lec/parallel_eal/run4/CB2008_incl_back_seis_fixed.csv");
		
		CSVFile<String> csv1 = CSVFile.readFile(file1, true);
		CSVFile<String> csv2 = CSVFile.readFile(file2, true);
		
		double maxDiff = 0;
		
		for (int row=1; row<csv1.getNumRows(); row++) {
			double val1 = Double.parseDouble(csv1.get(row, 1));
			double val2 = Double.parseDouble(csv2.get(row, 1));
			
			double pDiff = DataUtils.getPercentDiff(val2, val1);
			if (pDiff > maxDiff)
				maxDiff = pDiff;
		}
		
		System.out.println("Max Diff %: "+maxDiff);
	}

}
