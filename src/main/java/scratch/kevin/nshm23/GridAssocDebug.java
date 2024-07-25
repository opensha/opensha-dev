package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.geo.Region;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.modules.ModuleArchive;
import org.opensha.commons.util.modules.OpenSHA_Module;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.mpj.MPJ_GridSeisBranchBuilder;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultGridAssociations;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.MFDGridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.ModelRegion;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.gridded.NSHM23_FaultCubeAssociations;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SeisSmoothingAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.SeismicityRegions;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

public class GridAssocDebug {
	
	public static void main(String[] args) throws IOException {
		FaultSystemSolution inputSol = FaultSystemSolution.load(new File("/data/kevin/nshm23/batch_inversions/"
				+ "2023_04_11-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "node_branch_averaged/DM_GEOLOGIC.zip"));
		FaultSystemRupSet rupSet = inputSol.getRupSet();
		
		rupSet.removeModuleInstances(FaultGridAssociations.class);
		inputSol.removeModuleInstances(GridSourceProvider.class);
		
		File mainDir = new File("/home/kevin/OpenSHA/nshm23/grid_prov_debug");
		
		File beforeDir = new File(mainDir, "before_bugfix");
		File afterDir = new File(mainDir, "after_bugfix");

		// write them out
////		File outDir = new File(mainDir, "before_bugfix");
//		File outDir = new File(mainDir, "after_bugfix");
//		Region region = rupSet.requireModule(ModelRegion.class).getRegion();
//		NSHM23_FaultCubeAssociations assoc = NSHM23_InvConfigFactory.buildFaultCubeAssociations(rupSet, null, region);
//		
//		File assocFile = new File(outDir, MPJ_GridSeisBranchBuilder.GRID_ASSOCIATIONS_ARCHIVE_NAME);
//		ModuleArchive<OpenSHA_Module> archive = new ModuleArchive<>();
//		archive.addModule(assoc);
//		archive.write(assocFile);
//		
//		LogicTreeBranch<?> branch = NSHM23_LogicTreeBranch.DEFAULT_COMBINED;
//		
//		GridSourceProvider gridProv = NSHM23_InvConfigFactory.buildGridSourceProv(
//				inputSol, branch, List.of(SeismicityRegions.CONUS_WEST), assoc);
//		
//		File provFile = new File(outDir, MPJ_GridSeisBranchBuilder.AVG_GRID_SIE_PROV_ARCHIVE_NAME);
//		archive = new ModuleArchive<>();
//		archive.addModule(gridProv);
//		archive.write(provFile);
		
		
		File beforeProvFile = new File(beforeDir, MPJ_GridSeisBranchBuilder.AVG_GRID_SIE_PROV_ARCHIVE_NAME);
		File afterProvFile = new File(afterDir, MPJ_GridSeisBranchBuilder.AVG_GRID_SIE_PROV_ARCHIVE_NAME);
		
		MFDGridSourceProvider provBefore = new ModuleArchive<>(beforeProvFile).requireModule(MFDGridSourceProvider.class);
		MFDGridSourceProvider provAfter = new ModuleArchive<>(afterProvFile).requireModule(MFDGridSourceProvider.class);
		
		MinMaxAveTracker nodeRatePDiffTrack = new MinMaxAveTracker();
		MinMaxAveTracker nodeMoRatePDiffTrack = new MinMaxAveTracker();
		
		MinMaxAveTracker nodeSubSeisRatePDiffTrack = new MinMaxAveTracker();
		MinMaxAveTracker nodeSubSeisMoRatePDiffTrack = new MinMaxAveTracker();
		
		int numDiff = 0;
		for (int n=0; n<provBefore.getNumLocations(); n++) {
			IncrementalMagFreqDist mfdBefore = provBefore.getMFD(n);
			IncrementalMagFreqDist subSeisBefore = provBefore.getMFD_SubSeisOnFault(n);
			IncrementalMagFreqDist mfdAfter = provAfter.getMFD(n);
			IncrementalMagFreqDist subSeisAfter = provAfter.getMFD_SubSeisOnFault(n);
			
			double nodeRatePDiff = ratePDiff(mfdBefore, mfdAfter);
			double nodeMoRatePDiff = moRatePDiff(mfdBefore, mfdAfter);
			if (nodeRatePDiff > 0.01 && nodeMoRatePDiff > 0.01)
				numDiff++;
			
			double subRatePDiff = ratePDiff(subSeisBefore, subSeisAfter);
			double subMoRatePDiff = moRatePDiff(subSeisBefore, subSeisAfter);
			
			nodeRatePDiffTrack.addValue(nodeRatePDiff);
			nodeMoRatePDiffTrack.addValue(nodeMoRatePDiff);
			nodeSubSeisRatePDiffTrack.addValue(subRatePDiff);
			nodeSubSeisMoRatePDiffTrack.addValue(subMoRatePDiff);
		}
		
		System.out.println(numDiff+"/"+provBefore.getNumLocations()+" nodes differ");
		System.out.println("Node total rate % diffs:\t"+nodeRatePDiffTrack);
		System.out.println("Node total moment rate % diffs:\t"+nodeMoRatePDiffTrack);
		System.out.println("Node sub seis rate % diffs:\t"+nodeSubSeisRatePDiffTrack);
		System.out.println("Node sub seis moment rate % diffs:\t"+nodeSubSeisMoRatePDiffTrack);
	}
	
	private static double ratePDiff(IncrementalMagFreqDist mfd1, IncrementalMagFreqDist mfd2) {
		double rate1 = mfd1 == null ? 0d : mfd1.calcSumOfY_Vals();
		double rate2 = mfd2 == null ? 0d : mfd2.calcSumOfY_Vals();
		return pDiff(rate1, rate2);
	}
	
	private static double moRatePDiff(IncrementalMagFreqDist mfd1, IncrementalMagFreqDist mfd2) {
		double rate1 = mfd1 == null ? 0d : mfd1.getTotalMomentRate();
		double rate2 = mfd2 == null ? 0d : mfd2.getTotalMomentRate();
		return pDiff(rate1, rate2);
	}
	
	private static double pDiff(double val1, double val2) {
		if (val1 == 0d && val2 == 0d)
			return 0d;
		return 100d*(val2-val1)/val1;
	}

}
