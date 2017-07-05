package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityOptions;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.logicTree.APrioriBranchWeightProvider;
import scratch.UCERF3.logicTree.BranchWeightProvider;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.utils.MatrixIO;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class SubSectTDProbCICalc {

	public static void main(String[] args) throws ZipException, IOException {
		File compoundFile = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip");
		CompoundFaultSystemSolution cfss = CompoundFaultSystemSolution.fromZipFile(compoundFile);
		
		Region reg = new CaliforniaRegions.RELM_SOCAL();
		
		Map<FaultModels, Collection<Integer>> fmRupsMap = Maps.newHashMap();
		
		List<Double> probs = Lists.newArrayList();
		List<Double> weights = Lists.newArrayList();
		
		BranchWeightProvider weightProv = new APrioriBranchWeightProvider();
		
		MagDependentAperiodicityOptions[] covs = { null, MagDependentAperiodicityOptions.LOW_VALUES,
				MagDependentAperiodicityOptions.MID_VALUES, MagDependentAperiodicityOptions.HIGH_VALUES };
		
		File zipFileDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2014_10_07-ucerf3-erf-probs");
//		ZipFile zipFiles 
		
		for (LogicTreeBranch branch : cfss.getBranches()) {
			FaultModels fm = branch.getValue(FaultModels.class);
			if (!fmRupsMap.containsKey(fm)) {
				FaultSystemRupSet rupSet = cfss.getSolution(branch).getRupSet();
				HashSet<Integer> allRups = new HashSet<Integer>();
				sectLoop:
				for (FaultSectionPrefData sect : rupSet.getFaultSectionDataList()) {
					for (Location loc : sect.getFaultTrace()) {
						if (reg.contains(loc)) {
							allRups.addAll(rupSet.getRupturesForSection(sect.getSectionId()));
							continue sectLoop;
						}
					}
				}
				System.out.println(allRups.size()+" ruptures in region for "+fm.name());
				fmRupsMap.put(fm, allRups);
			}
			Collection<Integer> rups = fmRupsMap.get(fm);
			
			for (MagDependentAperiodicityOptions cov : covs) {
				
			}
			List<Double> rupProbs = Lists.newArrayList();
			
		}
	}
	
	private static double[] getProbs(ZipFile zip, LogicTreeBranch branch) throws IOException {
		// get the rate from the zip file
		String eName = branch.buildFileName()+".bin";
		ZipEntry probsEntry = zip.getEntry(eName);
		Preconditions.checkNotNull(probsEntry, "Entry not found in zip: "+eName);
		double[] probs = MatrixIO.doubleArrayFromInputStream(
				zip.getInputStream(probsEntry), probsEntry.getSize());
		return probs;
	}

}
