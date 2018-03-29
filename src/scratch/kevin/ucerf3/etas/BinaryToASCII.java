package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.util.FileNameComparator;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

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
		
		if (inputFile.isDirectory()) {
			List<File> subDirs = Lists.newArrayList(inputFile.listFiles());
			for (int i=subDirs.size(); --i>=0;) {
				File subDir = subDirs.get(i);
				if (!subDir.isDirectory() || !subDir.getName().startsWith("sim_"))
					subDirs.remove(i);
			}
			subDirs.sort(new FileNameComparator());
			if (num < 0)
				num = subDirs.size();
			else if (subDirs.size() < num)
				num = subDirs.size();
			for (int i=0; i<num; i++) {
				File simInFile = new File(subDirs.get(i), "simulatedEvents.bin");
				Preconditions.checkState(simInFile.exists());
				File outputFile;
				if (inputFile.equals(outputDir))
					outputFile = new File(subDirs.get(i), "simulatedEvents.txt");
				else
					outputFile = new File(outputDir, "catalog_"+i+".txt");
				List<ETAS_EqkRupture> catalog = ETAS_CatalogIO.loadCatalogBinary(simInFile);
				ETAS_CatalogIO.writeEventDataToFile(outputFile, catalog);
			}
		} else {
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

}
