package scratch.kevin.ucerf3.eal;

import java.io.File;
import java.io.IOException;

import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sra.calc.parallel.MPJ_CondLossCalc;

import com.google.common.base.Preconditions;

public class TractCalcVerify {

	public static void main(String[] args) throws IOException {
		if (args.length != 2) {
			System.err.println("Usage: <tract-dir> <combined-file>");
			System.exit(2);
		}
		File tractDir = new File(args[0]);
		File combinedFile = new File(args[1]);
		
		double[][] combined = MPJ_CondLossCalc.loadResults(combinedFile);
		double[][] myResults = null;
		
		int count = 0;
		for (File file : tractDir.listFiles()) {
			if (!file.getName().contains(".bin"))
				continue;
			System.out.println("Processing file "+(count++)+": "+file.getName());
			double[][] tractResults = MPJ_CondLossCalc.loadResults(file);
			if (myResults == null)
				myResults = tractResults;
			else
				MPJ_CondLossCalc.addTo(myResults, tractResults);
		}
		
		MinMaxAveTracker absTracker = new MinMaxAveTracker();
		MinMaxAveTracker percentTracker = new MinMaxAveTracker();
		
		for (int i=0; i<combined.length; i++) {
			if (combined[i] == null)
				continue;
			if (myResults[i] == null) {
				double sum = 0;
				for (double val : combined[i])
					sum += val;
				Preconditions.checkState(sum == 0d, "no tract resutls for "+i+", but combined non-zero: "+sum);
			} else {
				for (int j=0; j<combined[i].length; j++) {
					double v1 = combined[i][j];
					double v2 = myResults[i][j];
					double absDiff = Math.abs(v1 - v2);
					absTracker.addValue(absDiff);
					double pDiff = DataUtils.getPercentDiff(v2, v1);
					percentTracker.addValue(pDiff);
				}
			}
		}
		
		System.out.println("Absolute Differences:\n\t"+absTracker);
		System.out.println("Percent Differences:\n\t"+percentTracker);
	}

}
