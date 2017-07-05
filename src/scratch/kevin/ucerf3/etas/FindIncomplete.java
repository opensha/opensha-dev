package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import org.opensha.commons.util.FileNameComparator;

import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;

import com.google.common.base.Preconditions;

public class FindIncomplete {

	public static void main(String[] args) throws IOException {
		Preconditions.checkArgument(args.length == 1);
		File resultsDir = new File(args[0]);
		Preconditions.checkArgument(resultsDir.exists() && resultsDir.isDirectory(),
				"Results dir does not exist or is not a directory: %s", resultsDir.getAbsolutePath());
		
		File[] subDirs = resultsDir.listFiles();
		Arrays.sort(subDirs, new FileNameComparator());
		
		int failed = 0;
		int failedButLoads = 0;
		int tot = 0;
		
		for (File subDir : subDirs) {
			if (!subDir.isDirectory())
				continue;
			String name = subDir.getName();
			if (!name.startsWith("sim_"))
				continue;
			
			tot++;
			
			if (!MPJ_ETAS_Simulator.isAlreadyDone(subDir)) {
				failed++;
				
				// try loading
				String loadStr;
				try {
					File asciiFile = new File(subDir, "simulatedEvents.txt");
					File binaryFile = new File(subDir, "simulatedEvents.bin");
					if (binaryFile.exists())
						ETAS_CatalogIO.loadCatalog(binaryFile);
					else
						ETAS_CatalogIO.loadCatalog(asciiFile);
					loadStr = "true";
					failedButLoads++;
				} catch (Exception e) {
					loadStr = "false: "+e.getMessage();
				}
				System.out.println(name+" failed. Loads? "+loadStr);
			}
		}
		
		System.out.println(failed+"/"+tot+" failed ("+failedButLoads+" of which still load)");
	}

}
