package scratch.kevin.segmentation.normalFault;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.Test;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.InversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.JumpProbabilityConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SegmentationProbabilityConstraint;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.modules.SlipAlongRuptureModel;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Jump;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.SubsectionTransition;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.JumpProbabilityCalc;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.RuptureProbabilityCalc.BinaryRuptureProbabilityCalc;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.TransitionProbabilityCalc;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SectionDistanceAzimuthCalculator;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_ConstraintBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.JumpSegmentationModelBranchNode;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SegmentationMFD_Adjustment;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SegmentationModelBranchNode;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.gson.Gson;

public class NSHM23_NormalFaultSegmentationModelsTest {

	@Test
	public void testReplaceBranchLevelAndRoundTrip() {
		for (NSHM23_NormalFaultSegmentationModels choice : NSHM23_NormalFaultSegmentationModels.values()) {
			LogicTreeBranch<LogicTreeNode> branch = branch(choice);
			assertSame(choice, branch.requireValue(SegmentationModelBranchNode.class));
			assertNull(branch.getValue(NSHM23_SegmentationModels.class));
			assertNull(branch.getValue(JumpSegmentationModelBranchNode.class));
			Gson gson = new Gson();
			LogicTreeBranch<?> restored = gson.fromJson(
					gson.toJson(branch, LogicTreeBranch.class), LogicTreeBranch.class);
			assertSame(choice, restored.requireValue(SegmentationModelBranchNode.class));
			assertEquals(branch.buildFileName(), restored.buildFileName());
		}
	}

	@Test
	public void testFactoryIncludesInternalBendsWithoutJumps() {
		for (NSHM23_NormalFaultSegmentationModels choice : NSHM23_NormalFaultSegmentationModels.values()) {
			FaultSystemRupSet set = rupSet(false);
			assertFalse(NSHM23_InvConfigFactory.hasJumps(set));
			NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
			InversionConfiguration config = factory.buildInversionConfig(set, branch(choice), 1);
			assertNotNull(config);
			int found = 0;
			for (InversionConstraint constraint : config.getConstraints()) {
				assertFalse(constraint instanceof JumpProbabilityConstraint);
				if (constraint instanceof SegmentationProbabilityConstraint) {
					assertEquals(2, constraint.getNumRows());
					found++;
				}
			}
			assertEquals(1, found);
		}
	}

	@Test
	public void testExistingJumpBranchUsesLegacyConstraint() {
		FaultSystemRupSet set = rupSet(true);
		InversionConfiguration config = new NSHM23_InvConfigFactory().buildInversionConfig(
				set, NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT, 1);
		assertNotNull(config);
		assertTrue(config.getConstraints().stream().anyMatch(c -> c instanceof JumpProbabilityConstraint));
		assertFalse(config.getConstraints().stream().anyMatch(c -> c instanceof SegmentationProbabilityConstraint));
		assertSame(NSHM23_SegmentationModels.MID,
				NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT.requireValue(JumpSegmentationModelBranchNode.class));
		assertSame(NSHM23_SegmentationModels.MID,
				NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT.requireValue(SegmentationModelBranchNode.class));
		JumpSegmentationModelBranchNode jumpChoice = NSHM23_SegmentationModels.MID;
		SegmentationModelBranchNode<Jump> genericChoice = jumpChoice;
		JumpProbabilityCalc jumpModel = jumpChoice.getModel(set, null);
		TransitionProbabilityCalc<Jump> genericModel = genericChoice.getModel(set, null);
		ClusterRupture rupture = set.requireModule(ClusterRuptures.class).get(4);
		assertEquals(jumpModel.calcRuptureProb(rupture, false), genericModel.calcRuptureProb(rupture, false), 0d);
	}

	@Test
	public void testMFDAdjustmentIncludesInternalBends() {
		FaultSystemRupSet set = rupSet(false);
		TransitionProbabilityCalc<?> model = new TestNode(false, 0.05).getModel(set, null);
		List<? extends IncrementalMagFreqDist> baseline = targets(set, model, SegmentationMFD_Adjustment.NONE);
		for (SegmentationMFD_Adjustment method : new SegmentationMFD_Adjustment[] {
				SegmentationMFD_Adjustment.JUMP_PROB_THRESHOLD_AVG,
				SegmentationMFD_Adjustment.REL_GR_THRESHOLD_AVG_SINGLE_ITER,
				SegmentationMFD_Adjustment.REL_GR_THRESHOLD_AVG,
				SegmentationMFD_Adjustment.RUP_PROB_THRESHOLD_AVG,
				SegmentationMFD_Adjustment.RUP_MULTIPLY_WORST_JUMP_PROB}) {
			List<? extends IncrementalMagFreqDist> adjusted = targets(set, model, method);
			boolean changed = false;
			for (int s=0; s<baseline.size(); s++) {
				assertEquals(baseline.get(s).getTotalMomentRate(), adjusted.get(s).getTotalMomentRate(),
						baseline.get(s).getTotalMomentRate()*1e-10);
				for (int m=0; m<baseline.get(s).size(); m++) {
					assertTrue(Double.isFinite(adjusted.get(s).getY(m)));
					changed |= Math.abs(baseline.get(s).getY(m)-adjusted.get(s).getY(m)) > 1e-15;
				}
			}
			assertTrue("Internal bends should affect targets for "+method, changed);
		}
	}

	@Test
	public void testGenericTraversalMatchesJumpOnlyMFDSupport() {
		FaultSystemRupSet set = rupSet(true);
		JumpProbabilityCalc jumps = NSHM23_SegmentationModels.MID.getModel(set, null);
		TransitionProbabilityCalc<Jump> generic = new TransitionProbabilityCalc<Jump>() {
			@Override public String getName() { return jumps.getName(); }
			@Override public boolean isDirectional(boolean splayed) { return jumps.isDirectional(splayed); }
			@Override public Iterable<Jump> getTransitions(ClusterRupture rupture) {
				return rupture.getJumpsIterable();
			}
			@Override public double calcTransitionProbability(ClusterRupture rupture, Jump jump, boolean verbose) {
				return jumps.calcJumpProbability(rupture, jump, verbose);
			}
		};
		for (SegmentationMFD_Adjustment method : new SegmentationMFD_Adjustment[] {
				SegmentationMFD_Adjustment.JUMP_PROB_THRESHOLD_AVG,
				SegmentationMFD_Adjustment.REL_GR_THRESHOLD_AVG_SINGLE_ITER,
				SegmentationMFD_Adjustment.REL_GR_THRESHOLD_AVG,
				SegmentationMFD_Adjustment.RUP_PROB_THRESHOLD_AVG,
				SegmentationMFD_Adjustment.RUP_MULTIPLY_WORST_JUMP_PROB}) {
			List<? extends IncrementalMagFreqDist> legacy = targets(set, jumps, method);
			List<? extends IncrementalMagFreqDist> updated = targets(set, generic, method);
			for (int s=0; s<legacy.size(); s++)
				for (int m=0; m<legacy.get(s).size(); m++)
					assertEquals(legacy.get(s).getY(m), updated.get(s).getY(m), 1e-14);
		}
	}

	@Test
	public void testAllAdjustmentModesStillAcceptJumpModels() {
		JumpProbabilityCalc jumps = NSHM23_SegmentationModels.MID.getShawModel();
		for (SegmentationMFD_Adjustment method : SegmentationMFD_Adjustment.values()) {
			if (method == SegmentationMFD_Adjustment.NONE)
				assertNull(method.getAdjustment(jumps));
			else
				assertNotNull(method.getAdjustment(jumps));
		}
	}

	@Test(expected=IllegalArgumentException.class)
	public void testJumpOnlyAdjustmentRejectsInternalTransitionModel() {
		SegmentationMFD_Adjustment.GREEDY.getAdjustment(new TestNode(false, 0.05).getModel(rupSet(false), null));
	}

	@Test
	public void testZeroProbabilityExclusionsUseCombinedBound() {
		FaultSystemRupSet set = rupSet(false);
		ClusterRuptures ruptures = set.requireModule(ClusterRuptures.class);
		for (boolean varying : new boolean[] {false, true}) {
			SegmentationModelBranchNode<?> node = new TestNode(varying);
			BinaryRuptureProbabilityCalc exclusion = node.getExclusionModel(set, null);
			if (varying) {
				// Same boundary has p=0 in the whole-parent rupture but p=1 in the two-subsection rupture.
				assertNull(exclusion);
			} else {
				assertNotNull(exclusion);
				assertTrue(exclusion.isRupAllowed(ruptures.get(3), false)); // stops before bend
				assertFalse(exclusion.isRupAllowed(ruptures.get(4), false)); // crosses bend
				assertFalse(exclusion.isRupAllowed(ruptures.get(5), false));
			}
		}
	}

	@Test
	public void testJumpSpecializationPreservesLegacyZeroExclusions() {
		FaultSystemRupSet set = rupSet(true);
		ClusterRuptures ruptures = set.requireModule(ClusterRuptures.class);
		JumpProbabilityCalc model = new JumpProbabilityCalc() {
			@Override public String getName() { return "Rupture-dependent test jump"; }
			@Override public boolean isDirectional(boolean splayed) { return false; }
			@Override public double calcJumpProbability(ClusterRupture rupture, Jump jump, boolean verbose) {
				return rupture.getTotalNumSects() == 2 ? 1d : 0d;
			}
		};
		BinaryRuptureProbabilityCalc exclusion = new TestJumpNode(model).getExclusionModel(set, null);
		assertNotNull(exclusion);
		assertTrue(exclusion.isRupAllowed(ruptures.get(3), false));
		// Legacy jump rules exclude the pair if any occurrence is zero, rather than averaging first.
		assertFalse(exclusion.isRupAllowed(ruptures.get(4), false));
		assertFalse(exclusion.isRupAllowed(ruptures.get(5), false));
	}

	@Test
	public void testJumpSpecializationRetainsBinaryShortcutAndNullModel() {
		FaultSystemRupSet set = rupSet(true);
		JumpProbabilityCalc.BinaryJumpProbabilityCalc binary = new JumpProbabilityCalc.NoJumps();
		assertSame(binary, new TestJumpNode(binary).getExclusionModel(set, null));
		assertNull(new TestJumpNode(null).getExclusionModel(set, null));
	}

	private static List<? extends IncrementalMagFreqDist> targets(FaultSystemRupSet set,
			TransitionProbabilityCalc<?> model, SegmentationMFD_Adjustment method) {
		NSHM23_ConstraintBuilder builder = new NSHM23_ConstraintBuilder(set, 1d, null, false, false, false);
		builder.adjustForActualRupSlips(false, false).adjustForSegmentationModel(model, method);
		SupraSeisBValInversionTargetMFDs targets = builder.getTargetMFDs();
		return targets.getOnFaultSupraSeisNucleationMFDs();
	}

	private static LogicTreeBranch<LogicTreeNode> branch(NSHM23_NormalFaultSegmentationModels choice) {
		List<LogicTreeLevel<? extends LogicTreeNode>> levels = new ArrayList<>(NSHM23_LogicTreeBranch.levelsOnFault);
		levels.set(levels.indexOf(NSHM23_LogicTreeBranch.SEG), LogicTreeLevel.forEnum(
				NSHM23_NormalFaultSegmentationModels.class, "Normal Fault Segmentation", "NormalSeg"));
		LogicTreeBranch<LogicTreeNode> branch = new LogicTreeBranch<>(levels);
		for (LogicTreeNode node : NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT)
			branch.setValue(node instanceof NSHM23_SegmentationModels ? choice : node);
		branch.setValue(NSHM23_DeformationModels.GEOLOGIC);
		return branch;
	}

	private static FaultSystemRupSet rupSet(boolean jumps) {
		List<FaultSection> sections = Arrays.asList(section(0, 10, 0, 0, 0, 5),
				section(1, 10, 0, 5, 0, 10), section(2, jumps ? 20 : 10, jumps ? 6 : 0, 10, jumps ? 9 : 3, 15));
		List<List<Integer>> ids = Arrays.asList(Arrays.asList(0), Arrays.asList(1), Arrays.asList(2),
				Arrays.asList(0, 1), Arrays.asList(1, 2), Arrays.asList(0, 1, 2));
		List<ClusterRupture> ruptures = new ArrayList<>();
		SectionDistanceAzimuthCalculator dist = new SectionDistanceAzimuthCalculator(sections);
		for (List<Integer> indexes : ids) {
			List<FaultSection> participating = new ArrayList<>();
			for (int index : indexes) participating.add(sections.get(index));
			ruptures.add(ClusterRupture.forOrderedSingleStrandRupture(participating, dist));
		}
		return FaultSystemRupSet.builderForClusterRups(sections, ruptures)
				.forScalingRelationship(NSHM23_ScalingRelationships.WIDTH_LIMITED)
				.rupMags(new double[] {6.05, 6.05, 6.05, 6.35, 6.35, 6.65})
				.slipAlongRupture(new SlipAlongRuptureModel.Uniform()).build();
	}

	private static FaultSection section(int id, int parent, double... xy) {
		FaultSectionPrefData section = new FaultSectionPrefData();
		section.setSectionId(id);
		section.setParentSectionId(parent);
		section.setParentSectionName("Normal fault "+parent);
		section.setSectionName("Normal fault "+parent+", subsection "+id);
		section.setAveDip(50d);
		section.setAveRake(-90d);
		section.setAveUpperDepth(0d);
		section.setAveLowerDepth(15d);
		section.setAveSlipRate(1d);
		section.setSlipRateStdDev(0.1d);
		FaultTrace trace = new FaultTrace(section.getSectionName());
		Location origin = new Location(34, -112);
		for (int i=0; i<xy.length; i+=2)
			trace.add(LocationUtils.location(origin, Math.atan2(xy[i], xy[i+1]), Math.hypot(xy[i], xy[i+1])));
		section.setFaultTrace(trace);
		return section;
	}

	private static class TestJumpNode implements JumpSegmentationModelBranchNode {
		private final JumpProbabilityCalc model;
		TestJumpNode(JumpProbabilityCalc model) { this.model = model; }
		@Override public String getName() { return "Test jump node"; }
		@Override public String getShortName() { return "TestJump"; }
		@Override public String getFilePrefix() { return "TestJump"; }
		@Override public double getNodeWeight(LogicTreeBranch<?> branch) { return 1d; }
		@Override public JumpProbabilityCalc getModel(FaultSystemRupSet set, LogicTreeBranch<?> branch) {
			return model;
		}
	}

	private static class TestNode implements SegmentationModelBranchNode<SubsectionTransition> {
		private final boolean varying;
		private final double probability;
		TestNode(boolean varying) { this(varying, 0d); }
		TestNode(boolean varying, double probability) {
			this.varying = varying;
			this.probability = probability;
		}
		@Override public String getName() { return "Test internal transitions"; }
		@Override public String getShortName() { return "Test"; }
		@Override public String getFilePrefix() { return "Test"; }
		@Override public double getNodeWeight(LogicTreeBranch<?> branch) { return 1d; }
		@Override public TransitionProbabilityCalc<SubsectionTransition> getModel(
				FaultSystemRupSet set, LogicTreeBranch<?> branch) {
			return new TransitionProbabilityCalc<SubsectionTransition>() {
				@Override public String getName() { return "Test internal boundary"; }
				@Override public boolean isDirectional(boolean splayed) { return false; }
				@Override public Iterable<SubsectionTransition> getTransitions(ClusterRupture rupture) {
					List<SubsectionTransition> transitions = new ArrayList<>();
					for (FaultSection from : rupture.buildOrderedSectionList())
						for (FaultSection to : rupture.getTreeNavigator().getDescendants(from))
							if (Math.min(from.getSectionId(), to.getSectionId()) == 1)
								transitions.add(new SubsectionTransition(from, to));
					return transitions;
				}
				@Override public double calcTransitionProbability(
						ClusterRupture rupture, SubsectionTransition t, boolean verbose) {
					return varying && rupture.getTotalNumSects() == 2 ? 1d : probability;
				}
			};
		}
	}
}
