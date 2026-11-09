package scratch.kevin.segmentation.wasatchGeometrical;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.util.IDPairing;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SegmentationProbabilityConstraint;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.RupHistogramPlots.RakeType;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Transition;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Jump;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.NormalFaultGeometricSegmentationModel;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.JumpProbabilityCalc;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.NormalFaultGeometricSegmentationProb.CoefficientSet;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.NormalFaultGeometricSegmentationProb.Probabilities;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SingleStates;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;

import com.google.common.base.Preconditions;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;

/**
 * Usage: OUTPUT_DIR [REFERENCE_RUP_SET.zip]. Always audits the geologic Utah parent
 * sections. With a reference rupture set, uses NSHM23_SingleStates.UT to retain
 * Utah-connected clusters, then checks transition probabilities and matrix encoding.
 * This exercises the new constraint without changing NSHM23 inversion defaults.
 */
public class NormalFaultGeometricSegmentationSmokeTest {

	public static void main(String[] args) throws Exception {
		Preconditions.checkArgument(args.length >= 1 && args.length <= 2, "Usage: OUTPUT_DIR [REFERENCE_RUP_SET.zip]");
		File output = new File(args[0]);
		Preconditions.checkState(output.isDirectory() || output.mkdirs(), "Cannot create output directory");
		writeMechanismAudit(output);
		if (args.length == 2)
			runRegionalSmokeTest(output, new File(args[1]));
	}

	public static void writeMechanismAudit(File output) throws Exception {
		NSHM23_FaultModels faultModel = NSHM23_FaultModels.WUS_FM_v3;
		List<GeoJSONFaultSection> parents = new ArrayList<>();
		for (GeoJSONFaultSection parent : NSHM23_DeformationModels.getGeolFullSects(faultModel).values())
			if (NSHM23_SingleStates.UT.contains(parent)) parents.add(parent);
		parents.sort(Comparator.comparingInt(FaultSection::getSectionId));
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("parent_id", "name", "primary_state", "secondary_state", "geologic_rake_deg", "rake_type",
				"normal_within_10_deg", "normal_within_45_deg", "wasatch");
		Map<String, Integer> counts = new TreeMap<>();
		int normal = 0, broadlyNormal = 0, wasatch = 0, normalNonWasatch = 0;
		for (GeoJSONFaultSection parent : parents) {
			double rake = parent.getAveRake();
			String type = Double.isFinite(rake) ? RakeType.getType(rake).name() : "UNKNOWN";
			counts.merge(type, 1, Integer::sum);
			boolean isNormal = Double.isFinite(rake) && RakeType.NORMAL.isMatch(rake);
			boolean isBroadlyNormal = Double.isFinite(rake) && rake >= -135d && rake <= -45d;
			boolean isWasatch = parent.getSectionName().toLowerCase(java.util.Locale.ROOT).contains("wasatch");
			if (isNormal) normal++;
			if (isBroadlyNormal) broadlyNormal++;
			if (isWasatch) wasatch++;
			if (isNormal && !isWasatch) normalNonWasatch++;
			csv.addLine(""+parent.getSectionId(), parent.getSectionName(), parent.getProperty("PrimState", ""),
					parent.getProperty("SecState", ""), ""+rake, type, ""+isNormal, ""+isBroadlyNormal, ""+isWasatch);
		}
		csv.writeToFile(new File(output, "utah_fault_mechanisms.csv"));
		StringBuilder report = new StringBuilder("# Utah geologic fault mechanisms\n\n")
				.append("NSHM23 FM v3 with GEOLOGIC deformation-model rakes. Counts are parent fault sections, not subsections or named-fault groups.\n\n")
				.append("Selection uses NSHM23_SingleStates.UT.contains(), including secondary-state membership when INCLUDE_SECONDARY=")
				.append(NSHM23_SingleStates.INCLUDE_SECONDARY).append(". This is direct state membership, not the expanded connected-cluster rupture subset.\n\n")
				.append("| Category | Parent sections |\n| --- | ---: |\n")
				.append("| Total | ").append(parents.size()).append(" |\n")
				.append("| Normal: OpenSHA RakeType.NORMAL, rake -100 to -80 degrees inclusive | ").append(normal).append(" |\n")
				.append("| Not normal by that definition | ").append(parents.size()-normal).append(" |\n")
				.append("| Predominantly normal: rake -135 to -45 degrees inclusive | ").append(broadlyNormal).append(" |\n")
				.append("| Wasatch | ").append(wasatch).append(" |\n")
				.append("| Normal outside Wasatch | ").append(normalNonWasatch).append(" |\n\n")
				.append("Mechanism breakdown: ").append(counts).append(".\n");
		Files.writeString(new File(output, "utah_fault_mechanisms.md").toPath(), report, StandardCharsets.UTF_8);
		System.out.println(report);
	}

	private static void runRegionalSmokeTest(File output, File referenceFile) throws Exception {
		FaultSystemRupSet full = FaultSystemRupSet.load(referenceFile);
		FaultSystemRupSet utah = NSHM23_SingleStates.UT.getRuptureSubSet(full);
		ClusterRuptures ruptures = utah.requireModule(ClusterRuptures.class);
		StringBuilder report = new StringBuilder("# Utah normal-fault transition-model smoke test\n\n")
				.append("Reference: ").append(referenceFile.getAbsolutePath()).append("\n\n")
				.append("NSHM23_SingleStates.UT selects connected clusters, so neighboring out-of-state sections can be retained.\n\n")
				.append("Subsections: ").append(utah.getNumSections()).append("; ruptures: ").append(utah.getNumRuptures()).append(".\n\n");
		NSHM23_SegmentationModels backgroundChoice = NSHM23_SegmentationModels.MID;
		JumpProbabilityCalc background = NSHM23_SegmentationModels.buildModel(utah,
				backgroundChoice.getShawR0(), backgroundChoice.getShawShift(), 1d,
				backgroundChoice.getCreepingSectPassthrough(), Double.NaN, false);
		JumpProbabilityCalc original = NSHM23_SegmentationModels.MID.getModel(utah, null);
		int removedWasatchRules = 0;
		for (ClusterRupture rupture : ruptures)
			for (Jump jump : rupture.getJumpsIterable())
				if (original.calcJumpProbability(rupture, jump, false) != background.calcJumpProbability(rupture, jump, false))
					removedWasatchRules++;
		Preconditions.checkState(removedWasatchRules > 0, "Expected active legacy Wasatch rules in this test region");
		report.append("Background: NSHM23 MID without hardcoded Wasatch rules; ")
				.append(removedWasatchRules).append(" jump occurrences changed when removing those rules.\n\n");
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("first_id", "second_id", "first_name", "second_name", "coefficient_set", "bend_probability",
				"gap_probability", "stepover_or_offset_probability", "combined_probability");
		for (CoefficientSet coefficients : CoefficientSet.values()) {
			NormalFaultGeometricSegmentationModel model = new NormalFaultGeometricSegmentationModel(coefficients, 10d, background);
			Map<IDPairing, Transition> unique = new TreeMap<>();
			int outsideRuptures = 0;
			int nonWasatchNormalPairs = 0;
			for (int r=0; r<ruptures.size(); r++) {
				ClusterRupture rupture = ruptures.get(r);
				boolean anyNormal = false;
				for (FaultSection section : rupture.buildOrderedSectionList())
					anyNormal |= NormalFaultGeometricSegmentationModel.isNormal(section);
				List<Transition> transitions = model.getTransitions(rupture);
				if (!anyNormal) outsideRuptures++;
				for (Transition transition : transitions) {
					int first = transition.getFromSection().getSectionId(), second = transition.getToSection().getSectionId();
					IDPairing key = new IDPairing(Math.min(first, second), Math.max(first, second));
					unique.putIfAbsent(key, transition);
				}
			}
			Preconditions.checkState(!unique.isEmpty(), "Smoke test needs transitions");
			for (Transition transition : unique.values()) {
				if (!NormalFaultGeometricSegmentationModel.isNormal(transition.getFromSection())
						|| !NormalFaultGeometricSegmentationModel.isNormal(transition.getToSection())) continue;
				Probabilities p = model.getProbabilities(transition.getFromSection(), transition.getToSection());
				if (!transition.getFromSection().getParentSectionName().toLowerCase(java.util.Locale.ROOT).contains("wasatch")
						|| !transition.getToSection().getParentSectionName().toLowerCase(java.util.Locale.ROOT).contains("wasatch"))
					nonWasatchNormalPairs++;
				csv.addLine(""+transition.getFromSection().getSectionId(), ""+transition.getToSection().getSectionId(),
						transition.getFromSection().getSectionName(), transition.getToSection().getSectionName(), coefficients.name(),
						""+p.bendProbability, ""+p.gapProbability, ""+p.stepoverProbability, ""+p.combinedProbability);
			}
			SegmentationProbabilityConstraint constraint = new SegmentationProbabilityConstraint(1d, utah, model);
			DoubleMatrix2D matrix = new SparseDoubleMatrix2D(constraint.getNumRows(), utah.getNumRuptures());
			long nonzeros = constraint.encode(matrix, new double[matrix.rows()], 0);
			report.append(coefficients).append(": ").append(unique.size()).append(" unique transitions, ")
					.append(matrix.rows()).append(" constraint rows, ").append(nonzeros).append(" nonzero coefficients; ")
					.append(outsideRuptures).append(" ruptures contain no normal subsections; ")
					.append(nonWasatchNormalPairs).append(" normal transition pairs involve a non-Wasatch fault.\n\n");
		}
		csv.writeToFile(new File(output, "utah_normal_fault_transition_probabilities.csv"));
		Files.writeString(new File(output, "utah_normal_fault_smoke_test.md").toPath(), report, StandardCharsets.UTF_8);
		System.out.println(report);
	}
}
