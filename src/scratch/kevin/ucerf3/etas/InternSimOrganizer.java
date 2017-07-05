package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Map;

import org.opensha.commons.util.FileNameComparator;

import com.google.common.base.Preconditions;
import com.google.common.collect.Maps;
import com.google.common.io.Files;

public class InternSimOrganizer {
	
	public static void main(String[] args) throws IOException {
		File inputDir = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/2014_interns");
		File outputDir = new File(inputDir, "results");
		if (!outputDir.exists())
			outputDir.mkdir();
		
		int numDigits = 3;
		
//		String[] prefixes = {"2014_07_08-", "2014_07_09-"};
		String[] prefixes = {"2014_07_24-"};
		
		Map<String, Integer> countsMap = Maps.newHashMap();
		
		File[] files = inputDir.listFiles();
		Arrays.sort(files, new FileNameComparator());
		for (String prefix : prefixes) {
			for (File subDir : files) {
				if (!subDir.isDirectory())
					continue;
				String dirName = subDir.getName();
				if (!dirName.startsWith(prefix))
					continue;
				System.out.println("Processing "+dirName);
				String simName = dirName.substring(prefix.length());
				Integer curCount = countsMap.get(simName);
				if (curCount == null)
					curCount = 0;
				File resultsDir = new File(subDir, "results");
				Preconditions.checkState(resultsDir.exists());
				File[] simDirs = resultsDir.listFiles();
				Arrays.sort(simDirs, new FileNameComparator());
				File scenarioDir = new File(outputDir, simName);
				if (!scenarioDir.exists())
					scenarioDir.mkdir();
				for (File simDir : simDirs) {
					if (!simDir.getName().startsWith("sim_"))
						continue;
					try {
						if (!MPJ_ETAS_Simulator.isAlreadyDone(simDir)) {
							System.out.println(simDir.getName()+" for "+simName+" isn't completed, skipping");
							continue;
						}
					} catch (IOException e) {
						System.out.println(simDir.getName()+" for "+simName+" isn't completed, skipping");
						continue;
					}
					File infoFile = new File(simDir, "infoString.txt");
					Preconditions.checkState(infoFile.exists());
					File catFile = new File(simDir, "simulatedEvents.txt");
					Preconditions.checkState(catFile.exists());
					String outputDirName = (curCount++)+"";
					while (outputDirName.length() < numDigits)
						outputDirName = "0"+outputDirName;
					outputDirName = "sim_"+outputDirName;
					File simOutputDir = new File(scenarioDir, outputDirName);
					if (!simOutputDir.exists())
						simOutputDir.mkdir();
					Files.copy(infoFile, new File(simOutputDir, infoFile.getName()));
					Files.copy(catFile, new File(simOutputDir, catFile.getName()));
				}
				
				countsMap.put(simName, curCount);
			}
		}
	}

}
