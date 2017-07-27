package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.base.Preconditions;

import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;

public class BinaryToASCII {

	public static void main(String[] args) throws IOException {
		if (args.length < 2 || args.length > 3) {
			System.out.println("USAGE: BinaryToASCII <input-file> <output-dir> [<num-to-process>]");
			System.exit(2);
		}
		
		File inputFile = new File(args[0]);
		File outputDir = new File(args[1]);
		int num = -1;
		if (args.length > 2)
			num = Integer.parseInt(args[2]);
		Preconditions.checkArgument(inputFile.exists());
		Preconditions.checkArgument(outputDir.exists() || outputDir.mkdir());
		
		int count = 0;
		for (List<ETAS_EqkRupture> catalog : ETAS_CatalogIO.getBinaryCatalogsIterable(inputFile, -1)) {
			File file = new File(outputDir, "catalog_"+count+".txt");
			ETAS_CatalogIO.writeEventDataToFile(file, catalog);
			
			count++;
			if (num > 0 && count == num)
				break;
		}
	}

}
