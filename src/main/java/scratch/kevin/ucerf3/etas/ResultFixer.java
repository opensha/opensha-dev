package scratch.kevin.ucerf3.etas;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.util.ClassUtils;

import com.google.common.base.Preconditions;

import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO.BinarayCatalogsIterable;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;

public class ResultFixer {

	public static void main(String[] args) {
		if (args.length < 3) {
			System.err.println("USAGE: "+ClassUtils.getClassNameWithoutPackage(ResultFixer.class)+
					" <target-num> <output-file> <input 1> [...<input N>]");
			System.exit(2);
		}
		
		int targetNum = Integer.parseInt(args[0]);
		File outputFile = new File(args[1]);
		File[] inputs = new File[args.length-2];
		for (int i=0; i<inputs.length; i++) {
			inputs[i] = new File(args[i+2]);
			Preconditions.checkState(inputs[i].exists(), "Input doesn't exist: %s", inputs[i].getAbsolutePath());
		}
		
		try {
			DataOutputStream out = new DataOutputStream(new BufferedOutputStream(
					new FileOutputStream(outputFile), ETAS_CatalogIO.buffer_len));

			// write number of catalogs as int
			out.writeInt(targetNum);
			
			int curCount = 0;
			
			inputLoop:
			for (File input : inputs) {
				if (input.isDirectory()) {
					// results dir
					for (File subDir : input.listFiles()) {
						String name = subDir.getName();
						if (!subDir.isDirectory() || !name.startsWith("sim"))
							continue;
						File binFile = new File(subDir, "simulatedEvents.bin");
						if (!binFile.exists())
							continue;
						// load it
						List<ETAS_EqkRupture> catalog;
						try {
							catalog = ETAS_CatalogIO.loadCatalogBinary(binFile);
						} catch (Exception e) {
							System.out.println("Bad binary catalog, skipping: "+e.getMessage());
							continue;
						}
						ETAS_CatalogIO.writeCatalogBinary(out, catalog);
						curCount++;
						if (curCount == targetNum)
							break inputLoop;
					}
				} else {
					Preconditions.checkState(input.getName().endsWith(".bin"));
					BinarayCatalogsIterable iterable = ETAS_CatalogIO.getBinaryCatalogsIterable(input, -10);
					for (List<ETAS_EqkRupture> catalog : iterable) {
						ETAS_CatalogIO.writeCatalogBinary(out, catalog);
						curCount++;
						if (curCount == targetNum)
							break inputLoop;
					}
				}
			}
			out.close();
			Preconditions.checkState(curCount == targetNum,
					"Couldn't find enough good catalogs. Wrote "+curCount+"/"+targetNum);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
		System.exit(0);
	}

}
