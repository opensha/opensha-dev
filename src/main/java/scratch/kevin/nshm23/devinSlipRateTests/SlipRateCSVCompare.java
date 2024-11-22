package scratch.kevin.nshm23.devinSlipRateTests;

import java.io.File;
import java.io.IOException;

import org.opensha.commons.data.CSVFile;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;

public class SlipRateCSVCompare {

	public static void main(String[] args) throws IOException {
		File origFile = new File("/home/kevin/Downloads/slips_for_devin.csv");
		File modFile = new File("/home/kevin/Downloads/slips_modifiedby_devin.csv");
		
		CSVFile<String> origCSV = CSVFile.readFile(origFile, true);
		CSVFile<String> modCSV = CSVFile.readFile(modFile, true);
		
		for (NSHM23_DeformationModels dm : NSHM23_DeformationModels.values()) {
			int col = -1;
			String name = dm.getShortName();
			if (dm == NSHM23_DeformationModels.AVERAGE)
				name = "Branch-Average";
			for (int i=0; i<origCSV.getNumCols(); i++) {
				String heading = origCSV.get(0, i);
				if (heading.startsWith(name) && !heading.contains("Solution")) {
					col = i;
					System.out.println("Mapped "+dm.name()+" to "+heading+" (col="+col+")");
					break;
				}
			}
			if (col < 0) {
				System.out.println("Skipping DM: "+dm+"\n");
				continue;
			}
			
			for (int row=1; row<origCSV.getNumRows(); row++) {
				int sectID = origCSV.getInt(row, 0);
				String sectName = origCSV.get(row, 1);
				
				double origSlip = origCSV.getDouble(row, col);
				double modSlip = modCSV.getDouble(row, col);
				
				boolean same = (float)origSlip == (float)modSlip;
				if (!same && origSlip > 0d)
					same = Math.abs(origSlip - modSlip)/origSlip <= 0.001;
				
				if (!same)
					System.out.println(sectID+". "+sectName+":\t"+(float)origSlip+"\t->\t"+(float)modSlip);
			}
			System.out.println();
		}
	}

}
