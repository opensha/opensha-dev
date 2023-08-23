package scratch.kevin.nshm23.timJonEpistemic;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.Location;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_LogicTreeNode;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.imr.param.SiteParams.Vs30_TypeParam;

import com.google.common.base.Preconditions;

public class TimJonEpistemicCalcSetup {

	public static void main(String[] args) throws IOException {
		File dir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2023_06_23-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
		File inputTreeFile = new File(dir, "logic_tree_full_gridded.json");
		
		int downsample = 70000; // results in approx 50k unique
		File outputDir = new File(dir.getParentFile(), "2023_08_23-tim_jon_nshm23_v3_epistemic_calcs_downsampled");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		/*
		 * write logic tree
		 */
		LogicTree<?> inputLogicTree = LogicTree.read(inputTreeFile);
		
		List<LogicTreeLevel<? extends LogicTreeNode>> origLevels = new ArrayList<>();
		origLevels.addAll(inputLogicTree.getLevels());
		
		boolean includeIdriss = false;
		List<LogicTreeLevel<? extends LogicTreeNode>> addLevels = new ArrayList<>();
		LogicTreeLevel<NGAW2_LogicTreeNode> level = LogicTreeLevel.forEnum(
				NGAW2_LogicTreeNode.class, "NGA-W2 GMM", "NGA-W2");
		addLevels.add(level);
		
//		LogicTreeLevel<NGAW2_EpistemicLogicTreeNode> level = LogicTreeLevel.forEnum(
//				NGAW2_EpistemicLogicTreeNode.class, "NGA-W2 Additional Epistemic Uncertainty", "AddEpi");
//		addLevels.add(level);
		
		LogicTree<?> addTree = LogicTree.buildExhaustive(addLevels, true);
		if (!includeIdriss)
			addTree = addTree.matchingNone(NGAW2_LogicTreeNode.IDRISS_2014);
		
		int numBranches = addTree.size()*inputLogicTree.size();
		System.out.println("Have "+addTree.size()+" new sub-branches, total branches: "+numBranches);
		
		List<LogicTreeLevel<? extends LogicTreeNode>> allLevels = new ArrayList<>();
		allLevels.addAll(origLevels);
		allLevels.addAll(addLevels);
		
		List<LogicTreeBranch<LogicTreeNode>> allBranches = new ArrayList<>(numBranches);
		for (LogicTreeBranch<?> origBranch : inputLogicTree) {
			for (LogicTreeBranch<?> addBranch : addTree) {
				LogicTreeBranch<LogicTreeNode> branch = new LogicTreeBranch<>(allLevels);
				for (LogicTreeNode node : origBranch)
					branch.setValue(node);
				for (LogicTreeNode node : addBranch)
					branch.setValue(node);
				allBranches.add(branch);
			}
		}
		Preconditions.checkState(allBranches.size() == numBranches);
		
		LogicTree<?> fullTree = LogicTree.fromExisting(allBranches.get(0).getLevels(), allBranches);
		
		if (downsample > 0 && downsample < numBranches) {
			fullTree = fullTree.sample(downsample, false, new Random((long)downsample*(long)numBranches));
			System.out.println("Downsampled to "+fullTree.size()+" branches");
		}
		
		fullTree.write(new File(outputDir, "logic_tree.json"));
		
		/*
		 * write sites file
		 */
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("Name", "Latitude", "Longitude", Vs30_Param.NAME, Vs30_TypeParam.NAME);
		
		List<String> siteNames = new ArrayList<>();
		List<Location> siteLocs = new ArrayList<>();
		List<double[]> siteVs30sMeasured = new ArrayList<>();
		List<double[]> siteVs30sInferred = new ArrayList<>();
		
		siteNames.add("Davis");
		siteLocs.add(new Location(38.3210, -121.4530));
		siteVs30sMeasured.add(new double[] {
				216,
				236,
				251,
				266,
				281,
				299,
				327
		});
		siteVs30sInferred.add(new double[] {
				179,
				212,
				238,
				266,
				296,
				333,
				394
		});
		
		siteNames.add("Berkeley");
		siteLocs.add(new Location(37.5216, -122.1527));
		siteVs30sMeasured.add(new double[] {
				514,
				562,
				598,
				633,
				670,
				713,
				779
		}); 
		siteVs30sInferred.add(new double[] {
				288,
				402,
				510,
				633,
				785,
				996,
				1393
		});
		
		for (int s=0; s<siteNames.size(); s++) {
			String sitePrefix = siteNames.get(s);
			Location siteLoc = siteLocs.get(s);
			double[] vs30Measured = siteVs30sMeasured.get(s);
			double[] vs30Inferred = siteVs30sInferred.get(s);
			
			for (boolean inferred : new boolean[] {false,true}) {
				double[] vs30s = inferred ? vs30Inferred : vs30Measured;
				String typeParamValue = inferred ? Vs30_TypeParam.VS30_TYPE_INFERRED : Vs30_TypeParam.VS30_TYPE_MEASURED;
				
				for (double vs30 : vs30s) {
					String siteName = sitePrefix+" "+(int)vs30+" "+typeParamValue;
					csv.addLine(siteName, (float)siteLoc.getLatitude()+"", (float)siteLoc.getLongitude()+"",
							(int)vs30+"", typeParamValue);
				}
			}
		}
		
		csv.writeToFile(new File(outputDir, "hazard_sites.csv"));
	}

}
