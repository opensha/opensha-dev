package scratch.kevin.ucerf3.etas;

import java.io.File;

import org.opensha.commons.util.ClassUtils;

import com.google.common.base.Preconditions;

import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;

public class ETAS_BinaryCatalogFilterByMag {

	public static void main(String[] args) {
		if (args.length != 4) {
			System.err.println("USAGE: "+ClassUtils.getClassNameWithoutPackage(ETAS_BinaryCatalogFilterByMag.class)
						+" <input-binary-file> <output-binary-file> <min-mag> <preserve-chain>");
			System.exit(2);
		}
		
		File inputFile = new File(args[0]);
		Preconditions.checkArgument(inputFile.exists(), "Input file doesn't exist: %s", inputFile.getAbsolutePath());
		
		File outputFile = new File(args[1]);
//		Preconditions.checkArgument(outputFile.getParentFile().exists(),
//				"Output file cannot be created: %s", outputFile.getAbsolutePath());
		
		double minMag = Double.parseDouble(args[2]);
		
		boolean preserveChain = Boolean.parseBoolean(args[3]);
		
		try {
			ETAS_CatalogIO.binaryCatalogsFilterByMag(inputFile, outputFile, minMag, preserveChain);
		} catch (Exception e) {
			e.printStackTrace();
			System.err.flush();
			System.exit(1);
		}
		System.exit(0);
	}

}
