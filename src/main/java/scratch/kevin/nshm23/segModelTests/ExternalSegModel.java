package scratch.kevin.nshm23.segModelTests;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.IDPairing;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.Inversions;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Jump;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.PlausibilityConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.JumpProbabilityCalc;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.strategies.ClusterConnectionStrategy;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.RupsThroughCreepingSectBranchNode;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SegmentationModelBranchNode;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

public class ExternalSegModel {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/home/kevin/markdown/nshm23-misc/rsqsim_seg_tests");
		
		File refSolFile = new File(outputDir, "ref_sol.zip");
		
		File inputSolFile = new File("/home/kevin/Simulators/catalogs/bruce/rundir11082/rsqsim_11082_m6_skip10000_sectArea0.5.zip");
		FaultSystemSolution sol = FaultSystemSolution.load(inputSolFile);
		FaultSystemRupSet rupSet = sol.getRupSet();
		int threads = FaultSysTools.defaultNumThreads();
		
		outputDir = new File(outputDir, inputSolFile.getParentFile().getName());
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
		factory.setCacheDir(new File("/home/kevin/OpenSHA/nshm23/rup_sets/cache"));
		LogicTreeBranch<LogicTreeNode> branch = NSHM23_LogicTreeBranch.AVERAGE_COMBINED.copy();
		FaultSystemRupSet applyToRS = factory.buildRuptureSet(branch, threads);
		
		PlausibilityConfiguration applyToConfig = applyToRS.requireModule(PlausibilityConfiguration.class);
		ClusterConnectionStrategy applyToConnStrat = applyToConfig.getConnectionStrategy();
		
		// first figure out average participation rate on each parent
		Map<Integer, List<FaultSection>> sectsByParent = rupSet.getFaultSectionDataList().stream().collect(
				Collectors.groupingBy(S->S.getParentSectionId()));
		Map<Integer, Double> sectAvgRates = new HashMap<>();
		for (Integer parentID : sectsByParent.keySet()) {
			double sumRates = 0d;
			int num = 0;
			for (FaultSection sect : sectsByParent.get(parentID)) {
				sumRates += sol.calcTotParticRateForSect(sect.getSectionId());
				num++;
			}
			Preconditions.checkState(sumRates > 0d);
			sectAvgRates.put(parentID, sumRates/num);
		}
		
		// now figure out corupture rates
		Map<IDPairing, Double> corupRates = new HashMap<>();
		for (int rupIndex=0; rupIndex<rupSet.getNumRuptures(); rupIndex++) {
			List<Integer> parents = rupSet.getParentSectionsForRup(rupIndex);
			double rate = sol.getRateForRup(rupIndex);
			if (rate == 0d)
				continue;
			if (parents.size() > 1) {
				for (int i=0; i<parents.size(); i++) {
					int parent1 = parents.get(i);
					for (int j=i+1; j<parents.size(); j++) {
						int parent2 = parents.get(j);
						Preconditions.checkState(parent1 != parent2);
						IDPairing pair = parent1 > parent2 ? new IDPairing(parent2, parent1) : new IDPairing(parent1, parent2);
						Double prevRate = corupRates.containsKey(pair) ? corupRates.get(pair) : 0d;
						corupRates.put(pair, prevRate+rate);
					}
				}
			}
		}
		
		// convert to probs
		Map<IDPairing, Double> corupProbs = new LinkedHashMap<>();
		List<IDPairing> sortedPairs = new ArrayList<>(corupRates.keySet());
		Collections.sort(sortedPairs);
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("Parent ID 1", "Parent Name 1", "Parent ID 2", "Parent Name 2", "Corupture Rate");
		for (IDPairing pair : sortedPairs) {
			if (!applyToConnStrat.areParentSectsConnected(pair.getID1(), pair.getID2()))
				// not directly connected in the target rupture set, skip
				continue;
			double corupRate = corupRates.get(pair);
			double rate1 = sectAvgRates.get(pair.getID1());
			double rate2 = sectAvgRates.get(pair.getID2());
			String name1 = sectsByParent.get(pair.getID1()).getFirst().getParentSectionName();
			String name2 = sectsByParent.get(pair.getID2()).getFirst().getParentSectionName();
			double passthrough1 = corupRate/rate1;
			double passthrough2 = corupRate/rate2;
			System.out.println(pair+":\tcorup="+(float)corupRate);
			System.out.println("\t"+pair.getID1()+". "+name1+":\tpassthrough="+(float)corupRate+"/"+(float)rate1+"\t=\t"+(float)passthrough1);
			System.out.println("\t"+pair.getID2()+". "+name2+":\tpassthrough="+(float)corupRate+"/"+(float)rate2+"\t=\t"+(float)passthrough2);
			double passthrough = Math.min(1, Math.max(passthrough1, passthrough2)); // TODO min or max here?
			System.out.println("\tPassthrough:\t"+(float)passthrough);
			corupProbs.put(pair, passthrough);
			csv.addLine(pair.getID1()+"", name1, pair.getID2()+"", name2, (float)passthrough+"");
		}
		csv.writeToFile(new File(outputDir, "seg_rates.csv"));
		
		JumpProbabilityCalc.HardcodedJumpProb segModel = new JumpProbabilityCalc.HardcodedJumpProb("RSQSim", corupProbs, true, 0d);
		
		// write it as a map at the section level
		GeographicMapMaker mapMaker = new GeographicMapMaker(NSHM23_RegionLoader.loadFullConterminousWUS());
		mapMaker.setFaultSections(applyToRS.getFaultSectionDataList());
		HashSet<Jump> jumps = new HashSet<>();
		for (Jump jump : applyToConnStrat.getAllPossibleJumps()) {
			if (jump.fromSection.getSectionId() > jump.toSection.getSectionId())
				jump = jump.reverse();
			jumps.add(jump);
		}
		Map<Jump, Double> mapJumps = new HashMap<>();
		for (Jump jump : jumps) {
			double prob = segModel.calcJumpProbability(null, jump, false);
			if (prob > 0d)
				mapJumps.put(jump, prob);
			else
				mapJumps.put(jump, Double.NaN);
		}
		CPT cpt = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(0d, 1d);
		cpt.setNanColor(Color.GRAY);
		mapMaker.plotJumpScalars(mapJumps, cpt, "Passthrough rates");
		mapMaker.plot(outputDir, "seg_model", " ");
		
		String link = GeographicMapMaker.getGeoJSONViewerLink(
				"https://data.opensha.org/ftp/kmilner/markdown/nshm23-misc/rsqsim_seg_tests/"+outputDir.getName()+"/seg_model.geojson");
		System.out.println("GeoJSON link:\n\t"+link);
		
		if (!refSolFile.exists()) {
			System.out.println("Running reference inversion first");
			FaultSystemSolution refSol = Inversions.run(factory, branch, threads);
			refSol.setGridSourceProvider(factory.buildGridSourceProvider(sol, branch));
			refSol.write(refSolFile);
			System.out.println("DONE with ref sol");
		}
		
		CustomSegNode segNode = new CustomSegNode(segModel, false);
		
		// now run our inversion
		List<LogicTreeLevel<? extends LogicTreeNode>> levels = new ArrayList<>();
		List<LogicTreeNode> values = new ArrayList<>();
		for (int l=0; l<branch.size(); l++) {
			LogicTreeLevel<? extends LogicTreeNode> level = branch.getLevel(l);
			LogicTreeNode value = branch.getValue(l);
			if (value instanceof NSHM23_SegmentationModels) {
				levels.add(new CustomSegLevel(segNode));
				values.add(segNode);
			} else {
				// copy over directly
				levels.add(level);
				values.add(value);
			}
		}
		LogicTreeBranch<LogicTreeNode> customBranch = new LogicTreeBranch<>(levels, values);
		
		System.out.println("Running modified inversion");
		FaultSystemSolution invSol = Inversions.run(factory, customBranch, threads);
		invSol.setGridSourceProvider(factory.buildGridSourceProvider(sol, customBranch));
		invSol.write(new File(outputDir, "custom_seg_model_sol.zip"));
		System.out.println("DONE with custom sol");
	}
	
	private static class CustomSegNode implements SegmentationModelBranchNode, RupsThroughCreepingSectBranchNode {
		
		private JumpProbabilityCalc segModel;
		private boolean allowedThroughCreeping;

		public CustomSegNode(JumpProbabilityCalc segModel, boolean allowedThroughCreeping) {
			this.segModel = segModel;
			this.allowedThroughCreeping = allowedThroughCreeping;
		}

		@Override
		public double getNodeWeight(LogicTreeBranch<?> fullBranch) {
			return 1;
		}

		@Override
		public String getFilePrefix() {
			return "RSQSimSeg";
		}

		@Override
		public String getShortName() {
			return "RSQSim-Seg";
		}

		@Override
		public String getName() {
			return "RSQSim seg model";
		}

		@Override
		public boolean isExcludeRupturesThroughCreepingSect() {
			return !allowedThroughCreeping;
		}

		@Override
		public JumpProbabilityCalc getModel(FaultSystemRupSet rupSet, LogicTreeBranch<?> branch) {
			return segModel;
		}
		
	}
	
	private static class CustomSegLevel extends LogicTreeLevel<CustomSegNode> {
		
		private CustomSegNode node;

		public CustomSegLevel(CustomSegNode node) {
			this.node = node;
		}

		@Override
		public String getShortName() {
			return "RSQSim-Seg";
		}

		@Override
		public String getName() {
			return "RSQSim seg model";
		}

		@Override
		public Class<? extends CustomSegNode> getType() {
			return CustomSegNode.class;
		}

		@Override
		public List<? extends CustomSegNode> getNodes() {
			return List.of(node);
		}

		@Override
		public boolean isMember(LogicTreeNode node) {
			return node == this.node;
		}
		
	}

}
