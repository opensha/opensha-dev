package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.FileUtils;

import scratch.UCERF3.utils.MatrixIO;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class ComoundCopyAddRuns {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		File compoundDir = new File(args[0]); // unzipped compound file
		File newBinsDir = new File(args[1]);
		
		// first copy files over
		for (File file : newBinsDir.listFiles()) {
			if (file.isDirectory())
				continue;
			String name = file.getName();
			if (!name.endsWith(".bin"))
				continue;
			String newName = name.replaceAll("_run", "_rates_");
			FileUtils.copyFile(file, new File(compoundDir, newName));
		}
		
		// now redo the average
		
		// make list per branch
		Map<String, List<File>> branchFilesMap = Maps.newHashMap();
		
		for (File file : compoundDir.listFiles()) {
			if (file.isDirectory())
				continue;
			String name = file.getName();
			if (!name.endsWith(".bin"))
				continue;
			if (!name.contains("_rates_"))
				continue;
			String prefix = name.substring(0, name.indexOf("_rates"));
			List<File> files = branchFilesMap.get(prefix);
			if (files == null) {
				files = Lists.newArrayList();
				branchFilesMap.put(prefix, files);
			}
			files.add(file);
		}
		
		// now build means
		for (String prefix : branchFilesMap.keySet()) {
			List<File> files = branchFilesMap.get(prefix);
			Preconditions.checkState(files.size() == 10);
			double mult = 1d / (double)files.size();
			double[] runningRates = null;
			
			for (File file : files) {
				double[] rates = MatrixIO.doubleArrayFromFile(file);
				if (runningRates == null)
					runningRates = new double[rates.length];
				else
					Preconditions.checkState(runningRates.length == rates.length);
				for (int i=0; i<rates.length; i++)
					runningRates[i] += rates[i]*mult;
			}
			
			File meanFile = new File(compoundDir, prefix+"_rates.bin");
			// this should already exist, overwriting
			Preconditions.checkState(meanFile.exists());
			MatrixIO.doubleArrayToFile(runningRates, meanFile);
		}
	}

}
