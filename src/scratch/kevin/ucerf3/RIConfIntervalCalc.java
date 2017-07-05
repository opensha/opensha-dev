package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipException;

import com.google.common.collect.Maps;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.FaultModels;

public class RIConfIntervalCalc {

	public static void main(String[] args) throws ZipException, IOException {
		File compoundFile = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip");
		CompoundFaultSystemSolution cfss = CompoundFaultSystemSolution.fromZipFile(compoundFile);
		
//		int parentID = ;
		
		Map<FaultModels, List<Integer>> fmRupsMap = Maps.newHashMap();
		
		
	}

}
