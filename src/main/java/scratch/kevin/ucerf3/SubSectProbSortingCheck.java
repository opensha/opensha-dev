package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import org.dom4j.DocumentException;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.util.DataUtils;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityOptions;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.analysis.FaultSysSolutionERF_Calc;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.logicTree.APrioriBranchWeightProvider;
import scratch.UCERF3.logicTree.BranchWeightProvider;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.MatrixIO;

public class SubSectProbSortingCheck {

	public static void main(String[] args) throws ZipException, IOException, DocumentException {
		// *** these are all inputs ***
		
		// load in rup sets
		Map<FaultModels, FaultSystemRupSet> rupSets = Maps.newHashMap();
		rupSets.put(FaultModels.FM3_1, FaultSystemIO.loadRupSet(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip")));
		rupSets.put(FaultModels.FM3_2, FaultSystemIO.loadRupSet(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_2_MEAN_BRANCH_AVG_SOL.zip")));
		// branch weights
		BranchWeightProvider weightProv = new APrioriBranchWeightProvider();
		
		// probabilities that we are validating
//		File csvFile = new File("/tmp/sub_section_probabilities.csv");
		File csvFile = new File("/tmp/asdf/TimeDependent_AVE_RI_AVE_NORM_TIME_SINCE/all_30yr/sub_section_probabilities.csv");
		CSVFile<String> csv = CSVFile.readFile(csvFile, true);
		
		// branches for processing
		List<LogicTreeBranch> branches = Lists.newArrayList(CompoundFaultSystemSolution.fromZipFile(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip")).getBranches());
		
		// these zip files contain 30 year time dependent probabilities for each branch
		Map<MagDependentAperiodicityOptions, ZipFile> probsZipFiles = Maps.newHashMap();
		File probsDir = new File("/home/kevin/OpenSHA/UCERF3/time_dep_erf_probs");
		probsZipFiles.put(MagDependentAperiodicityOptions.LOW_VALUES,
				new ZipFile(new File(probsDir, "probs_30yr_LOW_VALUES.zip")));
		probsZipFiles.put(MagDependentAperiodicityOptions.MID_VALUES,
				new ZipFile(new File(probsDir, "probs_30yr_MID_VALUES.zip")));
		probsZipFiles.put(MagDependentAperiodicityOptions.HIGH_VALUES,
				new ZipFile(new File(probsDir, "probs_30yr_HIGH_VALUES.zip")));
		probsZipFiles.put(null,
				new ZipFile(new File(probsDir, "probs_30yr_POISSON.zip")));
		
		Map<String, List<Double>> subSectParticProbs = Maps.newHashMap();
		Map<String, List<Double>> subSectParticWeights = Maps.newHashMap();
		
		File outputCSV = new File("/tmp/output.csv");
		
		// assemble probabilities
		int branchCount = 0;
		for (LogicTreeBranch branch : branches) {
			if (branchCount % 10 == 0)
				System.out.println("Branch "+branchCount);
			branchCount++;
			FaultSystemRupSet rupSet = rupSets.get(branch.getValue(FaultModels.class));
			Map<MagDependentAperiodicityOptions, double[]> probsMap = Maps.newHashMap();
			String eName = branch.buildFileName()+".bin";
			for (MagDependentAperiodicityOptions cov : probsZipFiles.keySet()) {
				ZipFile zip = probsZipFiles.get(cov);
				ZipEntry probsEntry = zip.getEntry(eName);
				Preconditions.checkNotNull(probsEntry, "Entry not found in zip: "+eName);
				double[] probs = MatrixIO.doubleArrayFromInputStream(
						zip.getInputStream(probsEntry), probsEntry.getSize());
				Preconditions.checkState(probs.length == rupSet.getNumRuptures(),
						"Prob length mismatch, expected "+rupSet.getNumRuptures()+", got "+probs.length);
				probsMap.put(cov, probs);
			}
			
			double branchWeight = weightProv.getWeight(branch);
			
			for (int sectIndex=0; sectIndex<rupSet.getNumSections(); sectIndex++) {
				String name = rupSet.getFaultSectionData(sectIndex).getName();
				if (!subSectParticProbs.containsKey(name)) {
					subSectParticProbs.put(name, new ArrayList<Double>());
					subSectParticWeights.put(name, new ArrayList<Double>());
				}
				
				for (MagDependentAperiodicityOptions cov : probsMap.keySet()) {
					double weight = branchWeight * FaultSystemSolutionERF.getWeightForCOV(cov);
					double[] probs = probsMap.get(cov);
					
					List<Double> sectRupProbs = Lists.newArrayList();
					for (int rupIndex : rupSet.getRupturesForSection(sectIndex))
						sectRupProbs.add(probs[rupIndex]);
					Preconditions.checkState(!sectRupProbs.isEmpty());
					
					double meanParticProb = FaultSysSolutionERF_Calc.calcSummedProbs(sectRupProbs);
					
					subSectParticProbs.get(name).add(meanParticProb);
					subSectParticWeights.get(name).add(weight);
				}
			}
		}
		
		CSVFile<String> out = new CSVFile<String>(true);
		out.addLine("Name", "File Prob", "Calc Prob", "Abs Diff", "% Diff");
		
		double maxPDiff = 0d;
		String maxPDiffSect = null;
		double maxAbsDiff = 0d;
		String maxAbsDiffSect = null;
		
		for (int row=1; row<csv.getNumRows(); row++) {
			String name = csv.get(row, 0);
			double csvProb = Double.parseDouble(csv.get(row, 3));
			
			List<Double> myProbs =subSectParticProbs.get(name);
			Preconditions.checkNotNull(myProbs);
			List<Double> myWeights =subSectParticWeights.get(name);
			Preconditions.checkState(myProbs.size() == myWeights.size());
			
			double sumWeight = 0d;
			double myProb = 0;
			
			for (int i=0; i<myProbs.size(); i++) {
				myProb += myProbs.get(i)*myWeights.get(i);
				sumWeight += myWeights.get(i);
			}
			
			// normalize
			myProb /= sumWeight;
			
			System.out.println(name+": file="+csvProb+", calc="+myProb);
			
			double absDiff = Math.abs(myProb - csvProb);
			double pDiff = DataUtils.getPercentDiff(csvProb, myProb);

			if (absDiff > maxAbsDiff) {
				maxAbsDiff = absDiff;
				maxAbsDiffSect = name;
			}
			if (pDiff > maxPDiff) {
				maxPDiff = pDiff;
				maxPDiffSect = name;
			}
			
			out.addLine(name, csvProb+"", myProb+"", absDiff+"", pDiff+"");
		}
		
		System.out.println("Max absolute diff of "+maxAbsDiff+" on "+maxAbsDiffSect);
		System.out.println("Max % diff of "+maxPDiff+" % on "+maxPDiffSect);
		
		out.writeToFile(outputCSV);
	}

}
