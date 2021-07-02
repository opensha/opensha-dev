package scratch.kevin;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.List;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

import scratch.UCERF3.enumTreeBranches.FaultModels;

public class ModPolyTests {

	public static void main(String[] args) throws IOException {
		File baseDir = new File("/tmp/u3_poly");
		
		FaultModels fm = FaultModels.FM3_2;
		File refDir = new File(baseDir, "orig_buffers_"+fm.encodeChoiceString());
//		File refDir = new File(baseDir, "mod_buffers_orig_loc_"+fm.encodeChoiceString());
		File compDir = new File(baseDir, "mod_buffers_new_loc_"+fm.encodeChoiceString());
		
		int numSects = 0;
		int numSizeMismatch = 0;
		int numNotExactlyEqual = 0;
		int numNotApproxEqual = 0;
		
		for (File file : refDir.listFiles()) {
			if (!file.getName().endsWith(".txt"))
				continue;
			List<String> refLines = Files.readLines(file, Charset.defaultCharset());
			List<String> compLines = Files.readLines(new File(compDir, file.getName()), Charset.defaultCharset());
			System.out.println("Comparing "+file.getName());
			numSects++;
			if (refLines.size() != compLines.size()) {
				System.err.println("\tSIZE MISMATCH! "+refLines.size()+" != "+compLines.size());
				numSizeMismatch++;
			} else {
				int exactlyEqual = 0;
				int floatEqual = 0;
				for (int i=0; i<refLines.size(); i++) {
					double[] refPt = load(refLines.get(i));
					double[] compPt = load(compLines.get(i));
					if (refPt[0] == compPt[0] && refPt[1] == compPt[1]) {
						exactlyEqual++;
						floatEqual++;
					} else if ((float)refPt[0] == (float)compPt[0] && (float)refPt[1] == (float)compPt[1]) {
						floatEqual++;
					}
				}
				int num = refLines.size();
				if (exactlyEqual != num) {
					System.err.println("\t"+exactlyEqual+"/"+num+" are exactly equal");
					System.err.println("\t"+floatEqual+"/"+num+" are approximately equal");
					numNotExactlyEqual++;
					if (floatEqual != num)
						numNotApproxEqual++;
				}
			}
		}
		System.out.println(numSizeMismatch+"/"+numSects+" have a different number of points");
		System.out.println("Of those that have the same number of points:");
		System.out.println("\t"+numNotExactlyEqual+"/"+(numSects-numSizeMismatch)+" are not exactly equal (double precision)");
		System.out.println("\t"+numNotApproxEqual+"/"+(numSects-numSizeMismatch)+" are not approximately equal (float precision)");
	}
	
	private static double[] load(String line) {
		String[] split = line.split("\t");
		Preconditions.checkState(split.length == 2);
		double[] ret = new double[2];
		ret[0] = Double.parseDouble(split[0]);
		ret[1] = Double.parseDouble(split[1]);
		return ret;
	}

}
