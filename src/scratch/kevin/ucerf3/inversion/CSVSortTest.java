package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;
import java.util.Comparator;

import org.opensha.commons.data.CSVFile;

import com.google.common.base.Joiner;

public class CSVSortTest {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		File file = new File("/tmp/2012_10_01-another-sweep_misfits.csv");
		CSVFile<String> csv = CSVFile.readFile(file, true);
		
		Comparator<String> comparator = new Comparator<String>() {
			
			@Override
			public int compare(String o1, String o2) {
				try {
					double d1 = Double.parseDouble(o1);
					double d2 = Double.parseDouble(o2);
					return Double.compare(d1, d2);
				} catch (NumberFormatException e) {
					return o1.compareTo(o2);
				}
			}
		};
		
		int cols = csv.getNumCols();
		for (int col=cols; --col>=0;)
			csv.sort(col, 1, comparator);
		
		for (int row=0; row<csv.getNumRows(); row++)
			System.out.println(Joiner.on(",").join(csv.getLine(row)));
	}

}
