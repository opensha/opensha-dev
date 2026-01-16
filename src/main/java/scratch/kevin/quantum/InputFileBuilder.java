package scratch.kevin.quantum;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import org.apache.commons.math3.util.Precision;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.CSVWriter;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionInputGenerator;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionSolver;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.ConstraintWeightingType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.InversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SlipRateInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.CompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.IterationCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.modules.ConnectivityClusters;
import org.opensha.sha.earthquake.faultSysSolution.modules.NamedFaults;
import org.opensha.sha.earthquake.faultSysSolution.modules.RuptureSubSetMappings;
import org.opensha.sha.earthquake.faultSysSolution.reports.AbstractRupSetPlot;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen.PlotLevel;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.Shaw07JumpDistProb;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.ConnectivityCluster;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.faultSysSolution.util.SubSectionBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.estimators.GRParticRateEstimator;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;
import org.opensha.sha.simulators.stiffness.AggregatedStiffnessCalculator;
import org.opensha.sha.simulators.stiffness.SubSectStiffnessCalculator;
import org.opensha.sha.simulators.stiffness.AggregatedStiffnessCalculator.AggregationMethod;
import org.opensha.sha.simulators.stiffness.SubSectStiffnessCalculator.StiffnessType;
import org.opensha.sha.util.FocalMech;

import com.google.common.base.Preconditions;

public class InputFileBuilder {
	
	private static final LogicTreeBranch<?> REF_BRANCH = NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT;
	
	public static void main(String[] args) throws IOException {
		System.setProperty("java.awt.headless", "true");
		
//		boolean rebuildRupSet = false;
//		File outputDir = new File("/home/kevin/markdown/inversions/2024_10_15-quantum_test_problem");
//		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
//		int targetParentID = FaultSectionUtils.findSectionID(NSHM23_FaultModels.WUS_FM_v3.getFaultSections(), "Likely");
//		writeNSHM23Test(outputDir, targetParentID, rebuildRupSet, "Small NSHM23 Fault Cluster");
		
		File outputDir = new File("/home/kevin/markdown/inversions/2026_01_15-synthetic_quantum_test_problem");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		FocalMech[] mechs = {FocalMech.STRIKE_SLIP, FocalMech.REVERSE};
		double[] dists = {1d, 5d};
		double[] angles = {0d, 30d, 60d, 90d};
		double[] fractAlongs = {1d, 0.5d};
		double length = 50d;
		
		CSVFile<String> csv = new CSVFile<>(true);
		NSHM23_SegmentationModels refSeg = REF_BRANCH.requireValue(NSHM23_SegmentationModels.class);
		Shaw07JumpDistProb refShaw = refSeg.getShawModel();
		csv.addLine("Fault 1 mechanism",
				"Fault 2 mechanism",
				"Distance",
				"Angle",
				"Jumping point location",
				"NSHM23 "+refSeg.getShortName()+" Segmentation max fractional co-rupture",
				"Connected?",
				"Full fault 1 connects?",
				"Full fault 2 connects?",
				"Fraction of 1 corupture rate with 2",
				"Fraction of 2 corupture rate with 1");
		DecimalFormat oDF = new DecimalFormat("0.##");
		int count = 0;
		for (int m1=0; m1<mechs.length; m1++) {
			FocalMech mech1 = mechs[m1];
			for (int m2=m1; m2<mechs.length; m2++) {
				FocalMech mech2 = mechs[m2];
				for (double dist : dists) {
					for (double angle : angles) {
						for (double fractAlong : fractAlongs) {
							if (fractAlong < 1d && angle == 0d)
								continue;
							
							String dirName = mech1.name()+"_to_"+mech2.name();
							String name = mech1.toString()+" -> "+mech2.toString();
							
							dirName += "_"+oDF.format(dist)+"km";
							name += ", "+oDF.format(dist)+" km";
							
							dirName += "_"+oDF.format(angle)+"deg";
							name += ", "+oDF.format(angle)+" degrees";
							
							List<String> line = new ArrayList<>(csv.getNumCols());
							line.add(mech1.toString());
							line.add(mech2.toString());
							line.add((float)dist+"");
							line.add((float)angle+"");
							
							if (fractAlong == 1d) {
								dirName += "_from_end";
								name += ", from end";
								line.add("End");
							} else if (fractAlong == 0d) {
								dirName += "_from_beginning";
								name += ", from beginning";
								line.add("Beginning");
							} else if (fractAlong == 0.5d) {
								dirName += "_from_middle";
								name += ", from middle";
								line.add("Middle");
							} else {
								dirName += "_from_"+oDF.format(fractAlong)+"fract";
								name += ", from "+oDF.format(fractAlong)+" fract";
								line.add(oDF.format(fractAlong)+"x along");
							}
							
							double shawFract = refShaw.calcJumpProbability(dist);
							line.add((float)shawFract+"");
							
							System.out.println("Doing "+name);
							System.out.println("\t"+dirName);
							
							File subDir = new File(outputDir, dirName);
							Preconditions.checkState(subDir.exists() || subDir.mkdir());
							
							List<? extends FaultSection> sects = buildSyntheticFaultPair(mech1, length, mech2, length, dist, angle, fractAlong, angle);
							
							FaultSystemSolution sol = writeSyntheticTest(subDir, sects, name);
							FaultSystemRupSet rupSet = sol.getRupSet();
							int count1 = 0;
							int count2 = 0;
							for (FaultSection sect : rupSet.getFaultSectionDataList()) {
								if (sect.getParentSectionId() == 0)
									count1++;
								else if (sect.getParentSectionId() == 1)
									count2++;
								else
									throw new IllegalStateException("Unexpected parent: "+sect.getParentSectionId());
							}
							
							/**
							 * "Connected?",
				"Full fault 1 connects?",
				"Full fault 2 connects?",
				"Fraction of 1 corupture rate with 2",
				"Fraction of 2 corupture rate with 1");
							 */
							double rate1 = 0d;
							double rate2 = 0d;
							double connRate = 0d;
							boolean connects = false;
							boolean full1connects = false;
							boolean full2connects = false;
							for (int r=0; r<rupSet.getNumRuptures(); r++) {
								int num1 = 0;
								int num2 = 0;
								for (FaultSection sect : rupSet.getFaultSectionDataForRupture(r)) {
									if (sect.getParentSectionId() == 0)
										num1++;
									else if (sect.getParentSectionId() == 1)
										num2++;
									else
										throw new IllegalStateException("Unexpected parent: "+sect.getParentSectionId());
								}Preconditions.checkState(num1 > 0 || num2 > 0);
								double rate = sol.getRateForRup(r);
								if (num1 > 0)
									rate1 += rate;
								if (num2 > 0)
									rate2 += rate;
								if (num1 > 0 && num2 > 0) {
									connects = true;
									full1connects |= num1 == count1;
									full2connects |= num2 == count2;
									connRate += rate;
									
								}
							}
							line.add(connects+"");
							line.add(full1connects+"");
							line.add(full2connects+"");
							line.add((float)(connRate/rate1)+"");
							line.add((float)(connRate/rate2)+"");
							csv.addLine(line);
							
							count++;
						}
					}
				}
			}
		}
		System.out.println("Made "+count+" examples");
		csv.writeToFile(new File(outputDir, "connection_stats.csv"));
	}
	
	private static void writeNSHM23Test(File outputDir, int targetParentID, boolean rebuildRupSet, String name) throws IOException {
		File rsFile = new File(outputDir, "rup_set.zip");
		FaultSystemRupSet rupSet;
		if (rsFile.exists() && !rebuildRupSet) {
			rupSet = FaultSystemRupSet.load(rsFile);
		} else {
			NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
			factory.setCacheDir(new File("/home/kevin/OpenSHA/nshm23/rup_sets/cache"));
			LogicTreeBranch<LogicTreeNode> branch = NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT.copy();
			branch.setValue(NSHM23_DeformationModels.GEOLOGIC);
			rupSet = factory.buildRuptureSet(branch, 32);
			System.out.println("Building connectivity clusters");
			ConnectivityClusters clusters = ConnectivityClusters.build(rupSet);
			
			FaultSection targetSect = null;
			for (FaultSection sect : rupSet.getFaultSectionDataList()) {
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
			
			rupSet = rupSet.getForSectionSubSet(cluster.getSectIDs());
			rupSet.removeModuleInstances(NamedFaults.class);
			rupSet.write(rsFile);
		}
		
		SubSectStiffnessCalculator stiffnessCalc = new RuptureSets.CoulombRupSetConfig(
				rupSet.getFaultSectionDataList(), "test", NSHM23_ScalingRelationships.LOGA_C4p2).getStiffnessCalc();
		
		writeTestFiles(outputDir, rupSet, stiffnessCalc, name);
	}
	
	private static List<? extends FaultSection> buildSyntheticFaultPair(FocalMech mech1, double length1, FocalMech mech2, double length2,
			double distance, double strikeDiff, double jumpLocFractLen, double jumpDirection) {
		List<GeoJSONFaultSection> sects = new ArrayList<>(2);
		
		Location origin = new Location(0d, 0d);
		
		FaultTrace trace1 = new FaultTrace("Fault 1");
		trace1.add(origin);
		Location trace1end = LocationUtils.location(origin, 0d, length1);
		trace1.add(trace1end);
		
		GeoJSONFaultSection.Builder builder = new GeoJSONFaultSection.Builder(0, trace1.getName(), trace1);
		builder.upperDepth(0d).lowerDepth(12d);
		builder.dip(mech1.dip());
		builder.rake(mech1.rake());
		builder.slipRate(1d).slipRateStdDev(0.1d);
		GeoJSONFaultSection sect1 = builder.build();
		sects.add(sect1);
		
		Location jumpingOffPoint;
		if (jumpLocFractLen == 1d)
			jumpingOffPoint = trace1end;
		else if (jumpLocFractLen == 0d)
			jumpingOffPoint = origin;
		else
			jumpingOffPoint = LocationUtils.location(origin, 0d, jumpLocFractLen*length1);
		
		FaultTrace trace2 = new FaultTrace("Fault 2");
		Location trace2start = LocationUtils.location(jumpingOffPoint, Math.toRadians(jumpDirection), distance);
		if (jumpDirection != 0d || strikeDiff != 0d) {
			// probably need to correct the distance
			double testDist = LocationUtils.distanceToLineSegment(origin, trace1end, trace2start);
			double firstTestDist = testDist;
			int iters = 0;
			while (true) {
				if (Precision.equals(testDist, distance, 1e-3))
					break;
				trace2start = LocationUtils.location(trace2start, Math.toRadians(jumpDirection), distance-testDist);
				testDist = LocationUtils.distanceToLineSegment(origin, trace1end, trace2start);
				iters++;
			}
			System.out.println("Corrected off-angle jump distance from "+(float)firstTestDist+" to "
					+(float)testDist+" in "+iters+" iterations");
		}
		trace2.add(trace2start);
		Location trace2end = LocationUtils.location(trace2start, Math.toRadians(strikeDiff), length2);
		trace2.add(trace2end);
		
		builder = new GeoJSONFaultSection.Builder(1, trace2.getName(), trace2);
		builder.upperDepth(0d).lowerDepth(12d);
		builder.dip(mech2.dip());
		builder.rake(mech2.rake());
		builder.slipRate(1d).slipRateStdDev(0.1d);
		GeoJSONFaultSection sect2 = builder.build();
		sects.add(sect2);
		
		System.out.println(sect1.toFeature().toJSON());
		System.out.println(sect2.toFeature().toJSON());
		
		return sects;
	}
	
	private static FaultSystemSolution writeSyntheticTest(File outputDir, List<? extends FaultSection> sects, String name) throws IOException {
		List<FaultSection> subSects = SubSectionBuilder.buildSubSects(sects);
		System.out.println("\tBuilding rupture set");
		RuptureSets.CoulombRupSetConfig rsConfig = new RuptureSets.CoulombRupSetConfig(
				subSects, "test", REF_BRANCH.requireValue(NSHM23_ScalingRelationships.class));
		FaultSystemRupSet rupSet = rsConfig.build(1);
		return writeTestFiles(outputDir, rupSet, rsConfig.getStiffnessCalc(), name);
	}

	public static FaultSystemSolution writeTestFiles(File outputDir, FaultSystemRupSet rupSet,
			SubSectStiffnessCalculator stiffnessCalc, String name) throws IOException {
		// write out a-priori connectivity estimates
		System.out.println("\tBuilding estimated connection coeffs");
		GRParticRateEstimator estimator = new GRParticRateEstimator(rupSet, 0.5, REF_BRANCH.requireValue(NSHM23_SegmentationModels.class).getModel(rupSet, null));
		double[][] estCoeffs = calcConnCoeffs(rupSet, estimator.estimateRuptureRates());
		writeCoeffCSV(estCoeffs, new File(outputDir, "conn_coeffs_a_priori.csv"));
		
		System.out.println("\tConfiguring test inversion");
		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
		factory.setNumItersPerRup(NSHM23_InvConfigFactory.NUM_ITERS_PER_RUP_DEFAULT*10l);
		InversionConfiguration config = factory.buildInversionConfig(rupSet, REF_BRANCH, FaultSysTools.defaultNumThreads());
		
//		List<InversionConstraint> constraints = new ArrayList<>();
//		constraints.add(new SlipRateInversionConstraint(1d, ConstraintWeightingType.NORMALIZED, rupSet));
//		
//		CompletionCriteria completion = new IterationCompletionCriteria(100000000l);
////		CompletionCriteria completion = new IterationCompletionCriteria(1000l);
//		
//		InversionConfiguration config = InversionConfiguration.builder(constraints, completion)
////				.avgThreads(4, new IterationCompletionCriteria(10000l))
////				.threads(16).subCompletion(new IterationCompletionCriteria(1000l))
//				.build();
//		inputGen.generateInputs(true);
//		
//		InversionSolver.Default solver = new InversionSolver.Default();
		InversionSolver solver = factory.getSolver(rupSet, REF_BRANCH);
		
		System.out.println("\tGenerating test inversion inputs");
		InversionInputGenerator inputGen = new InversionInputGenerator(rupSet, config);
		inputGen.generateInputs();
		
		System.out.println("\tRunning test inversion");
		FaultSystemSolution sol = solver.run(rupSet, config, inputGen, null);

		sol.write(new File(outputDir, "solution.zip"));
		
		inputGen.writeArchive(new File(outputDir, "constraints.zip"), sol.getRateForAllRups(), false);
		
		System.out.println("\tWriting connection coeffs");
		double[][] solCoeffs = calcConnCoeffs(rupSet, sol.getRateForAllRups());
		writeCoeffCSV(solCoeffs, new File(outputDir, "conn_coeffs_sol.csv"));
		
		System.out.println("\tBuilding report");
		List<AbstractRupSetPlot> plots = ReportPageGen.getDefaultSolutionPlots(PlotLevel.REVIEW);
		ReportPageGen report = new ReportPageGen(rupSet, sol, name, new File(outputDir, "solution_report"), plots);
		report.setReplot(true);
		
		report.generatePage();
		
		System.out.println("\tBuilding stiffness data");
		writeStiffness(outputDir, rupSet, stiffnessCalc);
		
		System.out.println("\tDONE");
		
		return sol;
	}
	
	public static double[][] calcConnCoeffs(FaultSystemRupSet rupSet, double[] rates) {
		int numSects = rupSet.getNumSections();
		int numRups = rupSet.getNumRuptures();
		Preconditions.checkState(numRups == rates.length);
		double[][] coeffs = new double[numSects][numSects];
		BitSet[] sectRups = new BitSet[numSects];
		for (int s=0; s<numSects; s++) {
			sectRups[s] = new BitSet(numRups);
			for (int r : rupSet.getRupturesForSection(s))
				sectRups[s].set(r);
		}
		
		for (int s1=0; s1<numSects; s1++) {
			for (int s2=0; s2<numSects; s2++) {
				if (s1 == s2) {
					coeffs[s1][s2] = 1d;
				} else {
					double sectParticRate = 0d;
					double coruptureRate = 0d;
					for (int r=0; r<numRups; r++) {
						if (sectRups[s1].get(r)) {
							sectParticRate += rates[r];
							if (sectRups[s2].get(r))
								coruptureRate += rates[r];
						}
					}
					coeffs[s1][s2] = coruptureRate/sectParticRate;
				}
			}
		}
		
		return coeffs;
	}
	
	private static void writeCoeffCSV(double[][] coeffs, File outputFile) throws IOException {
		BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(outputFile));
		CSVWriter csv = new CSVWriter(out, true);
		
		csv.write(List.of("Primary Section Index", "Connected Section Index", "Connection Coefficient"));
		for (int s1=0; s1<coeffs.length; s1++)
			for (int s2=0; s2<coeffs.length; s2++)
				csv.write(List.of(s1+"", s2+"", (float)coeffs[s1][s2]+""));
		
		csv.flush();
		out.close();
	}
	
	public static void writeStiffness(File outputDir, FaultSystemRupSet rupSet, SubSectStiffnessCalculator stiffnessCalc) throws IOException {
		List<? extends FaultSection> subSects = rupSet.getFaultSectionDataList();
		
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

}
