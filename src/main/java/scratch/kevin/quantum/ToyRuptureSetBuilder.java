package scratch.kevin.quantum;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import org.opensha.commons.data.CSVWriter;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets;
import org.opensha.sha.earthquake.faultSysSolution.modules.ConnectivityClusters;
import org.opensha.sha.earthquake.faultSysSolution.modules.NamedFaults;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.strategies.ExhaustiveBilateralRuptureGrowingStrategy.SecondaryVariations;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.ConnectivityCluster;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.simulators.stiffness.AggregatedStiffnessCalculator;
import org.opensha.sha.simulators.stiffness.SubSectStiffnessCalculator;
import org.opensha.sha.simulators.stiffness.AggregatedStiffnessCalculator.AggregationMethod;
import org.opensha.sha.simulators.stiffness.SubSectStiffnessCalculator.StiffnessType;

import com.google.common.base.Preconditions;

public class ToyRuptureSetBuilder {

	public static void main(String[] args) throws IOException {
		boolean rebuildRupSet = false;
		
		String dirName = new SimpleDateFormat("yyyy_MM_dd").format(new Date());
//		String dirName = "2024_11_06";
		
		dirName += "-quantum_test_rup_set_problem";
		
//		String[] nameSearches = { "Likely" };
//		dirName += "-70_sects";
		
//		String[] nameSearches = { "Cedar", "Mountain", "Mahogany" };
//		dirName += "-124_sects";
		
		String[] nameSearches = { "San", "Felipe", "proxy" };
		dirName += "-213_sects";
		
		System.out.println(dirName);
		
		File outputDir = new File("/home/kevin/markdown/inversions/"+dirName);
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File rsFile = new File(outputDir, "orig_rup_set.zip");
		File exhaustiveFile = new File(outputDir, "exhaustive_rup_set.zip");
		FaultSystemRupSet origRupSet;
		if (rsFile.exists() && !rebuildRupSet) {
			origRupSet = FaultSystemRupSet.load(rsFile);
		} else {
			NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
			factory.setCacheDir(new File("/home/kevin/OpenSHA/nshm23/rup_sets/cache"));
			LogicTreeBranch<LogicTreeNode> branch = NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT.copy();
			branch.setValue(NSHM23_DeformationModels.GEOLOGIC);
			origRupSet = factory.buildRuptureSet(branch, 32);
			System.out.println("Building connectivity clusters");
			ConnectivityClusters clusters = ConnectivityClusters.build(origRupSet);
			
			int targetParentID = FaultSectionUtils.findParentSectionID(origRupSet.getFaultSectionDataList(), nameSearches);
			FaultSection targetSect = null;
			for (FaultSection sect : origRupSet.getFaultSectionDataList()) {
				if (sect.getParentSectionId() == targetParentID) {
					targetSect = sect;
					break;
				}
			}
			Preconditions.checkNotNull(targetSect);
			ConnectivityCluster cluster = null;
			for (ConnectivityCluster testCluster : clusters) {
				if (testCluster.containsSect(targetSect)) {
					cluster = testCluster;
					break;
				}
			}
			
			System.out.println("Cluster for "+targetSect.getSectionName()+" has "+cluster.getNumRuptures()
					+" ruptures on "+cluster.getNumSections()+" sub-sections");
			
			origRupSet = origRupSet.getForSectionSubSet(cluster.getSectIDs());
			origRupSet.removeModuleInstances(NamedFaults.class);
			origRupSet.write(rsFile);
		}
		List<? extends FaultSection> subSects = origRupSet.getFaultSectionDataList();
		
		RuptureSets.CoulombRupSetConfig rsConfig = new RuptureSets.CoulombRupSetConfig(
				subSects, null, NSHM23_ScalingRelationships.AVERAGE);
		// filters out indirect paths
		rsConfig.setNoIndirectPaths(false); // disabled
		// relative slip rate probability
		rsConfig.setSlipRateProb(0f); // disabled
		// fraction of interactions positive
		rsConfig.setCffFractInts(0.75f);
		// number of denominator values for the CFF favorability ratio
		rsConfig.setCffRatioN(2);
		// CFF favorability ratio threshold
		rsConfig.setCffRatioThresh(0.5f);
		// relative CFF probability
		rsConfig.setCffRelativeProb(0.01f);
		// if true, CFF calculations are computed with the most favorable path (up to max jump distance), which may not
		// use the exact jumping point from the connection strategy
		rsConfig.setFavorableJumps(true);
		// cumulative jump probability threshold
		rsConfig.setJumpProbThresh(0f); // disabled
		// cumulative rake change threshold
		rsConfig.setCmlRakeThresh(360f);
		// CONNECTION STRATEGY
		// maximum individual jump distance
		rsConfig.setMaxJumpDist(15d);
		// if true, connections happen at places that actually work and paths are optimized. if false, closest points
		rsConfig.setPlausibleConnections(true);
		// if >0 and <maxDist, connections will only be added above this distance when no other connections exist from
		// a given subsection. e.g., if set to 5, you can jump more than 5 km but only if no <= 5km jumps exist
		rsConfig.setAdaptiveMinDist(0d); // disabled
		// GROWING STRATEGY
		// if nonzero, apply thinning to growing strategy
		rsConfig.setAdaptiveSectFract(0f); // disabled
		// if true, allow bilateral rupture growing (see bilateralMode)
		rsConfig.setBilateral(true); // enabled
		// this controls how the secondary end of bilateral ruptures varies
		rsConfig.setBilateralVariationMode(SecondaryVariations.ALL);
		// if true, allow splays (using default settings)
		rsConfig.setSplays(true); // enabled
		
		FaultSystemRupSet rupSet;
		if (exhaustiveFile.exists() && !rebuildRupSet) {
			rupSet = FaultSystemRupSet.load(exhaustiveFile);
		} else {
			rupSet = rsConfig.build(FaultSysTools.defaultNumThreads());
			
			rupSet.write(exhaustiveFile);
		}
		System.out.println("Original had "+origRupSet.getNumRuptures());
		System.out.println("New has "+rupSet.getNumRuptures());
		
		SubSectStiffnessCalculator stiffnessCalc = rsConfig.getStiffnessCalc();
		
		List<AggregatedStiffnessCalculator> calcs = new ArrayList<>();
		List<String> calcPrefixes = new ArrayList<>();
		
		calcs.add(new AggregatedStiffnessCalculator(StiffnessType.CFF, stiffnessCalc, true,
				AggregationMethod.FLATTEN, AggregationMethod.SUM, AggregationMethod.SUM, AggregationMethod.SUM));
		calcPrefixes.add("cff_sum");
		
		calcs.add(new AggregatedStiffnessCalculator(StiffnessType.CFF, stiffnessCalc, true,
				AggregationMethod.FLATTEN, AggregationMethod.NUM_POSITIVE, AggregationMethod.SUM, AggregationMethod.NORM_BY_COUNT));
		calcPrefixes.add("cff_fraction_positive");
		
		System.out.println("Have "+subSects.size()+" sub-sections");
		
		for (int c=0; c<calcs.size(); c++) {
			AggregatedStiffnessCalculator calc = calcs.get(c);
			String prefix = calcPrefixes.get(c);
			
			System.out.println("Doing "+prefix);

			double[][] full = new double[subSects.size()][subSects.size()];
			double[][] average = new double[subSects.size()][subSects.size()];
			double[][] max = new double[subSects.size()][subSects.size()];
			
			for (int s1=0; s1<subSects.size(); s1++) {
				for (int s2=0; s2<subSects.size(); s2++) {
					if (s1 == s2)
						full[s1][s2] = Double.NaN;
					else
						full[s1][s2] = calc.calc(List.of(subSects.get(s1)), List.of(subSects.get(s2)));
				}
			}
			
			for (int s1=0; s1<subSects.size(); s1++) {
				for (int s2=s1; s2<subSects.size(); s2++) {
					average[s1][s2] = 0.5*(full[s1][s2] + full[s2][s1]);
					average[s2][s1] = average[s1][s2];
					max[s1][s2] = Math.max(full[s1][s2], full[s2][s1]);
					max[s2][s1] = max[s1][s2];
				}
			}

			writeMatrix(full, new File(outputDir, prefix+".csv"));
			writeMatrix(average, new File(outputDir, prefix+"_symmetrical_avg.csv"));
			writeMatrix(max, new File(outputDir, prefix+"_symmetrical_max.csv"));
		}
		
		// now write the actual ruptures
		writeRuptures(rupSet.getSectionIndicesForAllRups(), new File(outputDir, "orig_ruptures.csv"), subSects.size());
		writeRuptures(rupSet.getSectionIndicesForAllRups(), new File(outputDir, "exhaustive_ruptures.csv"), subSects.size());
	}
	
	private static void writeMatrix(double[][] mat, File file) throws IOException {
		BufferedOutputStream bout = new BufferedOutputStream(new FileOutputStream(file));
		CSVWriter writer = new CSVWriter(bout, true);
		
		List<String> header = new ArrayList<>(mat[0].length+1);
		header.add("");
		for (int i=0; i<mat[0].length; i++)
			header.add(i+"");
		writer.write(header);
		
		for (int i=0; i<mat.length; i++) {
			List<String> line = new ArrayList<>(mat[i].length+1);
			line.add(i+"");
			for (int j=0; j<mat[i].length; j++) {
				if (i == j)
					line.add("");
				else
					line.add((float)mat[i][j]+"");
			}
			writer.write(line);
		}
		
		writer.flush();
		
		bout.close();
	}
	
	private static void writeRuptures(List<List<Integer>> rups, File file, int numSects) throws IOException {
		BufferedOutputStream bout = new BufferedOutputStream(new FileOutputStream(file));
		CSVWriter writer = new CSVWriter(bout, true);
		
		List<String> header = new ArrayList<>(numSects+1);
		header.add("Rupture Index");
		for (int i=0; i<numSects; i++)
			header.add("Section "+i);
		writer.write(header);
		
		for (int i=0; i<rups.size(); i++) {
			List<String> line = new ArrayList<>(numSects+1);
			boolean[] rup = new boolean[numSects];
			for (int s : rups.get(i))
				rup[s] = true;
			line.add(i+"");
			for (int j=0; j<rup.length; j++)
				line.add(rup[j] ? "1" : "0");
			writer.write(line);
		}
		
		writer.flush();
		
		bout.close();
	}

}
