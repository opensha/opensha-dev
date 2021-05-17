package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;

import org.opensha.commons.data.CSVFile;

public class DefModelMax {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		File dir = new File("/home/kevin/OpenSHA/UCERF3/def_models/2012_02_07-initial");
		
		CSVFile<String> master = null;
		
		for (File file : dir.listFiles()) {
			if (!file.getName().contains("slip"))
				continue;
			CSVFile<String> csv = CSVFile.readFile(file, true);
			if (master == null) {
				master = csv;
			} else {
				for (int row=0; row<master.getNumRows(); row++) {
					double masterSlip = Double.parseDouble(master.get(row, 5));
					double newSlip = Double.parseDouble(csv.get(row, 5));
					if (newSlip > masterSlip)
						master.set(row, 5, newSlip+"");
				}
			}
		}
		
		master.writeToFile(new File(dir, "max_slips.csv"));
	}

}
