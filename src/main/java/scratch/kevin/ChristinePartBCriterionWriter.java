package scratch.kevin;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.sha.imr.attenRelImpl.ngaw2.FaultStyle;

import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig;

public class ChristinePartBCriterionWriter {

	public static void main(String[] args) throws IOException {
		
		double minDist = 1d;
		double maxDist = 20d;
		
		double mag = 7d;
		double vs30 = 760d;
		double z10 = 0.04; // km
		double z25 = 1d; // km
		double rake = 180;
		double dip = 90;
		double width = 20d;
		double zTor = 0d;
		double zHyp = 10d;
		FaultStyle style = FaultStyle.STRIKE_SLIP;
		
		CSVFile<String> csv = new CSVFile<>(true);
		
		for (double rRup=minDist; rRup<=maxDist; rRup++) {
			double rJB = rRup;
			double rX = rJB;
			DiscretizedFunc[] gmmMedians = BBP_PartBValidationConfig.calcNGA2_Medians(mag, rRup, rJB, rX, style, dip, zTor, width, vs30, z10, z25, zHyp);
			UncertainArbDiscFunc criterion = BBP_PartBValidationConfig.calcNGA2_Criterion(gmmMedians);
			
			if (csv.getNumRows() == 0) {
				// add header
				List<String> header = new ArrayList<>();
				header.add("Rrup (km)");
				for (int i=0; i<criterion.size(); i++) {
					header.add((float)criterion.getX(i)+" s Mean");
					header.add((float)criterion.getX(i)+" s Min");
					header.add((float)criterion.getX(i)+" s Max");
				}
				csv.addLine(header);
			}
			
			List<String> line = new ArrayList<>();
			line.add((float)rRup+"");
			for (int i=0; i<criterion.size(); i++) {
				line.add((float)criterion.getY(i)+"");
				line.add((float)criterion.getLowerY(i)+"");
				line.add((float)criterion.getUpperY(i)+"");
			}
			csv.addLine(line);
		}
		
		csv.writeToFile(new File("/tmp/m7_ss_ngaw2_bounds.csv"));
	}

}
