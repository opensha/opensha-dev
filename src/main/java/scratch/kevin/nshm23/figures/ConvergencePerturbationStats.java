package scratch.kevin.nshm23.figures;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.AnnealingProgress;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;

public class ConvergencePerturbationStats {

	public static void main(String[] args) throws IOException {
		List<String> modelNames = new ArrayList<>();
		List<File> modelFiles = new ArrayList<>();
		
		File invsDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		
		modelNames.add("UCERF3 Orig Calc Params");
		modelFiles.add(new File(invsDir, "2022_03_25-u3_branches-orig_calc_params-FM3_1/results.zip"));
		
		modelNames.add("UCERF3 New Anneal");
		modelFiles.add(new File(invsDir, "2022_03_24-u3_branches-FM3_1-2000ip/results.zip"));
		
		List<Double> avgIters = new ArrayList<>();
		List<Double> avgPerturbs = new ArrayList<>();
		List<Double> avgEnergies = new ArrayList<>();
		List<Double> avgEnergyAtPrev = new ArrayList<>();
		
		Double prevAvgIters = null;
		
		for (File modelFile : modelFiles) {
			ZipFile zip = new ZipFile(modelFile);
			
			Enumeration<? extends ZipEntry> entries = zip.entries();
			
			MinMaxAveTracker iterTrack = new MinMaxAveTracker();
			MinMaxAveTracker perturbTrack = new MinMaxAveTracker();
			MinMaxAveTracker energyTrack = new MinMaxAveTracker();
			MinMaxAveTracker energyAtPrevItersTrack = prevAvgIters == null ? null : new MinMaxAveTracker();
			
			while (entries.hasMoreElements()) {
				ZipEntry entry = entries.nextElement();
				if (entry.getName().endsWith(AnnealingProgress.PROGRESS_FILE_NAME)) {
					System.out.println("Loading progress from: "+entry.getName());
					CSVFile<String> csv = CSVFile.readStream(zip.getInputStream(entry), true);
					AnnealingProgress progress = new AnnealingProgress(csv);
					
					int last = progress.size()-1;
					long iters = progress.getIterations(last);
					long perturbs = progress.getNumPerturbations(last);
					double energy = progress.getEnergies(last)[0];
					
					iterTrack.addValue(iters);
					perturbTrack.addValue(perturbs);
					energyTrack.addValue(energy);
					
					if (energyAtPrevItersTrack != null) {
						double closestDiff = Double.POSITIVE_INFINITY;
						int closestIndex = -1;
						for (int i=0; i<progress.size(); i++) {
							long myIters = progress.getIterations(i);
							double diff = Math.abs((double)myIters - prevAvgIters);
							if (diff < closestDiff) {
								closestDiff = diff;
								closestIndex = i;
							}
						}
						energyAtPrevItersTrack.addValue(progress.getEnergies(closestIndex)[0]);
					}
				}
			}
			
			System.out.println("\tIterations: "+iterTrack);
			System.out.println("\tPetrubs: "+perturbTrack);
			System.out.println("\tEnergy: "+energyTrack);
			if (prevAvgIters != null) {
				System.out.println("\tEnergy at "+prevAvgIters.floatValue()+" iters: "+energyAtPrevItersTrack);
				avgEnergyAtPrev.add(energyAtPrevItersTrack.getAverage());
			} else {
				avgEnergyAtPrev.add(null);
			}
			
			avgIters.add(iterTrack.getAverage());
			avgPerturbs.add(perturbTrack.getAverage());
			avgEnergies.add(energyTrack.getAverage());
			
			zip.close();
			
			prevAvgIters = iterTrack.getAverage();
		}
		
		for (int m=0; m<modelNames.size(); m++) {
			System.out.println(modelNames.get(m)+":");
			double iters = avgIters.get(m);
			double perturbs = avgPerturbs.get(m);
			System.out.println((float)(iters/perturbs)+" iterations per perturbation");
			double energy = avgEnergies.get(m);
			if (m > 0) {
				double prevEnergy = avgEnergies.get(m-1);
				double change = (energy-prevEnergy)/prevEnergy;
				System.out.println("\tEnergy: "+(float)energy+" ("+new DecimalFormat("0.0%").format(change)+")");
				double energyAtPrevIters = avgEnergyAtPrev.get(m);
				change = (energyAtPrevIters-prevEnergy)/prevEnergy;
				System.out.println("\tEnergy at prev iters: "+(float)energyAtPrevIters+" ("+new DecimalFormat("0.0%").format(change)+")");
			} else {
				System.out.println("\tEnergy: "+(float)energy);
			}
		}
	}

}
