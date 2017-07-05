package scratch.kevin.cybershake.ucerf3.safWallToWallTests;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.dom4j.DocumentException;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.CompoundSurface;
import org.opensha.sha.faultSurface.RuptureSurface;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.erf.mean.RuptureCombiner;
import scratch.UCERF3.utils.FaultSystemIO;

public class SubsetSolutionGenerator {

	public static void main(String[] args) throws IOException, DocumentException {
		File solFile = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
			+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip");
		File outputDir = new File("/home/kevin/CyberShake/ucerf3/saf_wall_to_wall_tests");
		
		FaultSystemSolution sol = FaultSystemIO.loadSol(solFile);
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		HashSet<Integer> parents = new HashSet(FaultModels.FM3_1.getNamedFaultsMapAlt().get("San Andreas"));
		
		Region soCal = new CaliforniaRegions.RELM_SOCAL();
		
		Map<Integer, Boolean> safSectsInSoCal = Maps.newHashMap();
		for (int sectIndex=0; sectIndex<rupSet.getNumSections(); sectIndex++) {
			FaultSectionPrefData sect = rupSet.getFaultSectionData(sectIndex);
			if (!parents.contains(sect.getParentSectionId()))
				continue;
			boolean inside = false;
			for (Location loc : sect.getFaultTrace()) {
				if (soCal.contains(loc)) {
					inside = true;
					break;
				}
			}
			safSectsInSoCal.put(sectIndex, inside);
		}
		
		List<Integer> rupsWithPartial = Lists.newArrayList();
		
		int numSAF = 0;
		int numPartiallySAFSoCal = 0;
		int numOnlySAFSoCal = 0;
		rupLoop:
		for (int rupIndex = 0; rupIndex<rupSet.getNumRuptures(); rupIndex++) {
			for (int parent : rupSet.getParentSectionsForRup(rupIndex)) {
				if (!parents.contains(parent))
					continue rupLoop;
			}
			numSAF++;
			boolean partial = false;
			boolean only = true;
			for (int sectIndex : rupSet.getSectionsIndicesForRup(rupIndex)) {
				boolean inside = safSectsInSoCal.get(sectIndex);
				partial = partial || inside;
				only = only && inside;
			}
			if (partial) {
				numPartiallySAFSoCal++;
				rupsWithPartial.add(rupIndex);
			}
			if (only)
				numOnlySAFSoCal++;
		}
		
		System.out.println(numSAF+"/"+rupSet.getNumRuptures()+" ruputures are only on SAF");
		System.out.println(numPartiallySAFSoCal+" of those are at least partially in SoCal");
		System.out.println(numOnlySAFSoCal+" of those are entirely in SoCal");
		
		FaultSystemSolution subsetSol = new RuptureCombiner.SubsetSolution(sol, rupsWithPartial);
		System.out.println("Subset sol has "+subsetSol.getRupSet().getNumRuptures()+" ruptures");
		
		int longestSAF = -1;
		double longestSAFLength = 0d;
		int longestSAF_soCal = -1;
		double longestSoCalSAFLength = 0d;
		FaultSystemRupSet subsetRupSet = subsetSol.getRupSet();
		for (int rupIndex=0; rupIndex<subsetRupSet.getNumRuptures(); rupIndex++) {
			boolean partial = false;
			boolean only = true;
			for (int sectIndex : subsetRupSet.getSectionsIndicesForRup(rupIndex)) {
				boolean inside = safSectsInSoCal.get(sectIndex);
				partial = partial || inside;
				only = only && inside;
			}
			double len = subsetRupSet.getLengthForRup(rupIndex)/1000d;
			if (len > longestSAFLength) {
				longestSAFLength = len;
				longestSAF = rupIndex;
			}
			if (only && len > longestSoCalSAFLength) {
				longestSoCalSAFLength = len;
				longestSAF_soCal = rupIndex;
			}
		}
		
		FaultSystemIO.writeSol(subsetSol, new File(outputDir, "saf_subset_sol.zip"));
		
		System.out.println("Longest All SAF: "+longestSAF);
		System.out.println("Longest So Cal SAF: "+longestSAF_soCal);
		
		CompoundSurface surf = (CompoundSurface) rupSet.getSurfaceForRupupture(longestSAF, 1d, false);
		System.out.println();
		System.out.println("Longest Rupture:");
		System.out.println("\tLength: "+longestSAFLength+" (km)");
		System.out.println("\tAve Width: "+surf.getAveWidth()+" (km)");
		double minTopDepth = Double.POSITIVE_INFINITY;
		double maxBotDepth = 0;
		double maxWidth = 0;
		double minWidth = Double.POSITIVE_INFINITY;
		for (RuptureSurface subSurf : surf.getSurfaceList()) {
			minTopDepth = Math.min(minTopDepth, subSurf.getAveRupTopDepth());
			maxBotDepth = Math.max(maxBotDepth, subSurf.getAveRupTopDepth()+subSurf.getAveWidth());
			maxWidth = Math.max(maxWidth, subSurf.getAveWidth());
			minWidth = Math.min(minWidth, subSurf.getAveWidth());
		}
		System.out.println("\tMax Width: "+maxWidth+" (km)");
		System.out.println("\tMin Width: "+minWidth+" (km)");
		System.out.println("\tMin Top Depth: "+minTopDepth+" (km)");
		System.out.println("\tMax Bot Depth: "+maxBotDepth+" (km)");
		System.out.println("\tFirst loc: "+surf.getFirstLocOnUpperEdge());
		System.out.println("\tLast loc: "+surf.getLastLocOnUpperEdge());
		
		System.out.println();
		System.out.println("UCERF3 Surface Point Counts");
		double[] gridSpacings = { 0.1, 0.2 };
		
		for (double gridSpacing : gridSpacings) {
			long count = 0;
			
			for (FaultSectionPrefData sect : rupSet.getFaultSectionDataList()) {
				count += sect.getStirlingGriddedSurface(gridSpacing).size();
			}
			
			System.out.println((int)(gridSpacing*1000d)+"m: "+count);
		}
	}

}
