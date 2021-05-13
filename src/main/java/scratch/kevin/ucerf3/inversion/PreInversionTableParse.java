package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;

import org.opensha.commons.data.CSVFile;

import com.google.common.base.Joiner;

public class PreInversionTableParse {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		CSVFile<String> csv = CSVFile.readFile(new File("/tmp/output.csv"), true);
		
//		System.out.println(csv.get(0, 15));
		
		int cnt = 0;
		int tot = 0;
		
		for (int row=1; row<csv.getNumRows(); row++) {
			if (!csv.get(row, 4).equals("GRConst"))
				continue;
			double impCC = Double.parseDouble(csv.get(row, 15));
			if (impCC > 0.8) {
				cnt++;
				System.out.println(Joiner.on("\t").join(csv.getLine(row).subList(0, 9)));
			}
			tot++;
		}
		
		System.out.println("CNT: "+cnt+"/"+tot);
	}

}
