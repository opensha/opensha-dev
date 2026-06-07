package scratch.kevin.mfdInversion;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.jfree.data.Range;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.FileNameUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets.SimpleAzimuthalRupSetConfig;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.inversion.Inversions;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.InversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.JumpProbabilityConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.CompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.IterationCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.params.GenerationFunctionType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.params.NonnegativityConstraintType;
import org.opensha.sha.earthquake.faultSysSolution.modules.AveSlipModule;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.SectSlipRates;
import org.opensha.sha.earthquake.faultSysSolution.modules.SlipAlongRuptureModel;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionSlipRates;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportMetadata;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen.PlotLevel;
import org.opensha.sha.earthquake.faultSysSolution.reports.RupSetMetadata;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Jump;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.JumpProbabilityCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.SubSectionBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_ConstraintBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SegmentationMFD_Adjustment;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SubSectConstraintModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.estimators.GRParticRateEstimator;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;

public class MFD_InversionTests {

	public static void main(String[] args) throws IOException {
		System.setProperty("java.awt.headless", "true");
//		File outputDir = new File("/tmp/mfd_inv_tests");
		File outputDir = new File("C:\\Users\\kmilner\\scratch\\mfd_inv_tests");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir(),
				"Can't create output dir: %s", outputDir.getAbsolutePath());
		// demo Y-shaped fault system
		double upperDepth = 0d;
		double lowerDepth = 12d;
		double dip = 90d;
		double rake = 0d;
		double mainSlipRate = 10d;
		double fractSlipSD = 0.1;
		
//		double splitSlipRate1 = 2d;
//		double splitSlipRate2 = 2d;
		
//		double splitSlipRate1 = 5d;
//		double splitSlipRate2 = 5d;
		
		double splitSlipRate1 = 3d;
		double splitSlipRate2 = 1d;
		
		int mainID = 0;
		int splitID1 = 1;
		int splitID2 = 2;
		
		Location origin = new Location(0d,0d);
		Location splitPoint = LocationUtils.location(origin, 0d, 50d);
		Location splitEnd1 = LocationUtils.location(splitPoint, -Math.PI/4d, 50d);
		Location splitEnd2 = LocationUtils.location(splitPoint, Math.PI/4d, 50d);
		
		Region mapReg = new Region(new Location(-0.1, -0.5), new Location(1, 0.5));
		
		// MFD params
		double b = 0.5d;
		
		double segRate1 = 1d;
		double segRate2 = 1d;
		
//		double segRate1 = 0.1d;
//		double segRate2 = 0.1d;
		
//		double segRate1 = 0.05d;
//		double segRate2 = 0.05d;
		
		SubSectConstraintModels sectConstrModel = SubSectConstraintModels.TOT_NUCL_RATE;
//		SubSectConstraintModels sectConstrModel = SubSectConstraintModels.NUCL_MFD;
		
		// misc
		int threads = 16;
		
		List<FaultSection> sects = new ArrayList<>();
		
		DecimalFormat oDF = new DecimalFormat("0.##");
		
		String mainName = "Main Fault ("+oDF.format(mainSlipRate)+" mm/yr)";
		String splitName1 = "Split Fault 1 ("+oDF.format(splitSlipRate1)+" mm/yr";
		if (segRate1 < 1)
			splitName1 += ", P="+oDF.format(segRate1);
		splitName1 += ")";
		String splitName2 = "Split Fault 2 ("+oDF.format(splitSlipRate2)+" mm/yr";
		if (segRate2 < 1)
			splitName2 += ", P="+oDF.format(segRate2);
		splitName2 += ")";
		
		sects.add(new GeoJSONFaultSection.Builder(mainID, mainName,
				FaultTrace.of(origin, splitPoint))
				.lowerDepth(lowerDepth)
				.upperDepth(upperDepth)
				.dip(dip)
				.rake(rake)
				.slipRate(mainSlipRate)
				.slipRateStdDev(mainSlipRate*fractSlipSD)
				.build());
		sects.add(new GeoJSONFaultSection.Builder(splitID1, splitName1,
				FaultTrace.of(splitPoint, splitEnd1))
				.lowerDepth(lowerDepth)
				.upperDepth(upperDepth)
				.dip(dip)
				.rake(rake)
				.slipRate(splitSlipRate1)
				.slipRateStdDev(splitSlipRate1*fractSlipSD)
				.build());
		sects.add(new GeoJSONFaultSection.Builder(splitID2, splitName2,
				FaultTrace.of(splitPoint, splitEnd2))
				.lowerDepth(lowerDepth)
				.upperDepth(upperDepth)
				.dip(dip)
				.rake(rake)
				.slipRate(splitSlipRate2)
				.slipRateStdDev(splitSlipRate2*fractSlipSD)
				.build());
		
		List<FaultSection> subSects = SubSectionBuilder.buildSubSects(sects);
		
		// colors by distanct to split point (for individual mfd plots)
		double[] sectDistFromSplit = new double[subSects.size()];
		Map<Integer, Double> parentMaxDistFromSplit = new HashMap<>();
		for (FaultSection sect : subSects) {
			FaultTrace trace = sect.getFaultTrace();
			Location middle = new Location(0.5*(trace.first().lat+trace.last().lat),
					0.5*(trace.first().lon+trace.last().lon));
			double dist = LocationUtils.horzDistanceFast(splitPoint, middle);
			sectDistFromSplit[sect.getSectionId()] = dist;
			int parentID = sect.getParentSectionId();
			if (!parentMaxDistFromSplit.containsKey(parentID))
				parentMaxDistFromSplit.put(parentID, dist);
			else
				parentMaxDistFromSplit.put(parentID, Math.max(dist, parentMaxDistFromSplit.get(parentID)));
		}
		CPT distFromSplitCPT = new CPT(0d, 1d, Color.DARK_GRAY, Color.LIGHT_GRAY);
		Color[] sectDistFromSplitColors = new Color[subSects.size()];
		for (int s=0; s<subSects.size(); s++) {
			double parentMax = parentMaxDistFromSplit.get(subSects.get(s).getParentSectionId());
			sectDistFromSplitColors[s] = distFromSplitCPT.getColor(sectDistFromSplit[s]/parentMax);
		}
		
		SimpleAzimuthalRupSetConfig rsConfig = new SimpleAzimuthalRupSetConfig(subSects, NSHM23_ScalingRelationships.AVERAGE);
		
		JumpProbabilityCalc segModel = null;
		if (segRate1 < 1d || segRate2 < 1d) {
			segModel = new JumpProbabilityCalc() {

				@Override
				public boolean isDirectional(boolean splayed) {
					return false;
				}

				@Override
				public String getName() {
					return "Custom Seg model";
				}

				@Override
				public double calcJumpProbability(ClusterRupture fullRupture, Jump jump, boolean verbose) {
					int parent1 = Integer.min(jump.fromCluster.parentSectionID, jump.toCluster.parentSectionID);
					int parent2 = Integer.max(jump.fromCluster.parentSectionID, jump.toCluster.parentSectionID);
					if (parent1 != mainID)
						return 0d;
					if (parent2 == splitID1)
						return segRate1;
					if (parent2 == splitID2)
						return segRate2;
					throw new IllegalStateException("Unexpected jump: "+jump);
				}
				
			};
		}
		FaultSystemSolution origSol = null;
		FaultSystemSolution modSol = null;
		
		for (boolean orig : new boolean[] {true,false}) {
			FaultSystemRupSet rupSet = rsConfig.build(1);
			
			long equivNumVars = Long.max(rupSet.getNumRuptures(), rupSet.getNumSections()*100l);
			// to get some extra stability
			equivNumVars *= 2l;
			NSHM23_ConstraintBuilder constrBuilder = new NSHM23_ConstraintBuilder(rupSet, b);
			
			constrBuilder.magDepRelStdDev(M->NSHM23_InvConfigFactory.MFD_MIN_FRACT_UNCERT*Math.max(1, Math.pow(10, b*0.5*(M-6))));
			if (orig) {
				System.out.println("Running original recipe inversion");
				
				constrBuilder.adjustForActualRupSlips(NSHM23_ConstraintBuilder.ADJ_FOR_ACTUAL_RUP_SLIPS_DEFAULT,
						NSHM23_ConstraintBuilder.ADJ_FOR_SLIP_ALONG_DEFAULT);
				
				if (segModel != null)
					constrBuilder.adjustForSegmentationModel(segModel);
			} else {
				System.out.println("Running updated inversion");
				constrBuilder.applyCustomMFDAdjustment(new MFD_InversionAdjustment(threads, segModel));
			}
			
			double slipWeight = 1d;
			double mfdWeight = sectConstrModel == SubSectConstraintModels.NUCL_MFD ? 1 : 10;
			double nuclWeight = sectConstrModel == SubSectConstraintModels.TOT_NUCL_RATE ? 0.5 : 0d;
			double nuclMFDWeight = sectConstrModel == SubSectConstraintModels.NUCL_MFD ? 0.5 : 0d;
			
			if (slipWeight > 0d)
				constrBuilder.slipRates().weight(slipWeight);
			
			if (mfdWeight > 0d)
				constrBuilder.supraBValMFDs().weight(mfdWeight);
			
			if (nuclWeight > 0d)
				constrBuilder.sectSupraRates().weight(nuclWeight);
			
			if (nuclMFDWeight > 0d)
				constrBuilder.sectSupraNuclMFDs().weight(nuclMFDWeight);
			
			List<InversionConstraint> constraints = constrBuilder.build();
			
			SupraSeisBValInversionTargetMFDs targetMFDs = rupSet.requireModule(SupraSeisBValInversionTargetMFDs.class);
			
			GRParticRateEstimator rateEst = new GRParticRateEstimator(rupSet, targetMFDs);
			
			if (segModel != null) {
				constraints = new ArrayList<>(constraints);
				
//				InitialModelParticipationRateEstimator rateEst = new InitialModelParticipationRateEstimator(
//						rupSet, Inversions.getDefaultVariablePerturbationBasis(rupSet));

//				double weight = 0.5d;
//				boolean ineq = false;
				double weight = 100000d;
				boolean ineq = true;
				
				constraints.add(new JumpProbabilityConstraint.RelativeRate(
						weight, ineq, rupSet, segModel, rateEst));
			}
			
			int avgThreads = Integer.max(1, threads / 4);
			
			CompletionCriteria completion = new IterationCompletionCriteria(equivNumVars*NSHM23_InvConfigFactory.NUM_ITERS_PER_RUP_DEFAULT);
			CompletionCriteria subCompletion = new IterationCompletionCriteria(equivNumVars);
			CompletionCriteria avgCompletion = new IterationCompletionCriteria(equivNumVars*50l);
			
			InversionConfiguration config = InversionConfiguration.builder(constraints, completion)
					.threads(threads)
					.subCompletion(subCompletion)
					.avgThreads(avgThreads, avgCompletion)
					.perturbation(GenerationFunctionType.VARIABLE_EXPONENTIAL_SCALE)
					.nonNegativity(NonnegativityConstraintType.TRY_ZERO_RATES_OFTEN)
//					.sampler(sampler)
					.reweight()
					.variablePertubationBasis(rateEst.estimateRuptureRates())
					.build();
			
			System.out.println("Running the inversion");
			FaultSystemSolution sol = Inversions.run(rupSet, config);
			String slipPrefix;
			String slipTitle;
			if (orig) {
				origSol = sol;
				sol.write(new File(outputDir, "orig_sol.zip"));
				slipPrefix = "slips_orig";
				slipTitle = "Original Recipe";
			} else {
				modSol = sol;
				sol.write(new File(outputDir, "mod_sol.zip"));
				slipPrefix = "slips_mod";
				slipTitle = "MFD Inversion";
			}
			
			SectSlipRates targetSlips = rupSet.getSectSlipRates();
			SolutionSlipRates solSlips = SolutionSlipRates.calc(sol,
					rupSet.requireModule(AveSlipModule.class), rupSet.requireModule(SlipAlongRuptureModel.class));
			
			GeographicMapMaker mapMaker = new GeographicMapMaker(mapReg);
			mapMaker.setFaultSections(subSects);
			
			if (orig) {
				CPT slipCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(0d, mainSlipRate);
				mapMaker.plotSectScalars(S->1e3*targetSlips.getSlipRate(S.getSectionId()),
						slipCPT, "Slip Rates (mm/yr)");
				mapMaker.plot(outputDir, slipPrefix, " ");
			}
			
			CPT diffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-1, 1);
			mapMaker.plotSectScalars(S->1e3*(solSlips.get(S.getSectionId()) - targetSlips.getSlipRate(S.getSectionId())),
					diffCPT, "Slip Rate Misfit (mm/yr)");
			mapMaker.plot(outputDir, slipPrefix+"_diff", slipTitle);
			
			CPT pDiffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-10, 10);
			mapMaker.plotSectScalars(S->100*(solSlips.get(S.getSectionId()) - targetSlips.getSlipRate(S.getSectionId()))
					/targetSlips.getSlipRate(S.getSectionId()),
					pDiffCPT, "Slip Rate Misfit (%)");
			mapMaker.plot(outputDir, slipPrefix+"_pdiff", slipTitle);
		}
		
		Color targetColor = Colors.tab_aqua;
		Color origTargetColor = Colors.tab_olive;
		Color solColor = Color.BLACK;
		Color origSolColor = Colors.tab_brown;
		
		// write MFD plots
		for (int i=-1; i<sects.size(); i++) {
			int parentID;
			String name, prefix;
			if (i < 0) {
				parentID = -1;
				name = "Regional Sum";
				prefix = "mfd_sum";
			} else {
				FaultSection sect = sects.get(i);
				parentID = sect.getSectionId();
				name = sect.getName();
				String nameForPrefix = name;
				if (nameForPrefix.contains("("))
					nameForPrefix = nameForPrefix.substring(0, nameForPrefix.indexOf("(")).trim();
				prefix = "mfd_indv_"+FileNameUtils.simplify(nameForPrefix);
			}
			List<DiscretizedFunc> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			List<? extends IncrementalMagFreqDist> targets = getTargets(modSol, parentID);
			SummedMagFreqDist sum = sum(targets);
			SummedMagFreqDist origSum = sum(getTargets(origSol, parentID));
			
			if (i < 0) {
				// add individual targets
				for (FaultSection sect : sects) {
					SummedMagFreqDist sectSum = sum(getTargets(modSol, sect.getSectionId()));
					sectSum.setName(sect.getName());
					funcs.add(sectSum);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
				}
			} else {
				List<Color> colors = new ArrayList<>(targets.size());
				for (FaultSection sect : subSects)
					if (sect.getParentSectionId() == parentID)
						colors.add(sectDistFromSplitColors[sect.getSectionId()]);
				Preconditions.checkState(colors.size() == targets.size());
				for (int j=0; j<targets.size(); j++) {
					IncrementalMagFreqDist target = targets.get(j);
					target.setName(j == targets.size()/2 ? "Subsections" : null);
					funcs.add(target);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, colors.get(j)));
				}
			}
			
			IncrementalMagFreqDist origSolMFD, solMFD;
			if (parentID < 0) {
				origSolMFD = origSol.calcTotalNucleationMFD(sum.getMinX(), sum.getMaxX(), sum.getDelta());
				solMFD = modSol.calcTotalNucleationMFD(sum.getMinX(), sum.getMaxX(), sum.getDelta());
			} else {
				origSolMFD = origSol.calcNucleationMFD_forParentSect(parentID, sum.getMinX(), sum.getMaxX(), sum.size());
				solMFD = modSol.calcNucleationMFD_forParentSect(parentID, sum.getMinX(), sum.getMaxX(), sum.size());
			}
			
			origSum.setName("Original Target");
			funcs.add(origSum);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, origTargetColor));
			
			origSolMFD.setName("Original Solution");
			funcs.add(origSolMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 2f, origSolColor));
			
			sum.setName("Inverted Target");
			funcs.add(sum);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, targetColor));
			
			solMFD.setName("Solution");
			funcs.add(solMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 4f, solColor));
			
			HeadlessGraphPanel gp = PlotUtils.initScreenHeadless();
			
			PlotSpec plot = new PlotSpec(funcs, chars, name, "Magnitude", "Incremental Nucleation Rate (1/yr)");
			plot.setLegendInset(true);
			
			gp.drawGraphPanel(plot, false, true, new Range(6d, 8d), new Range(1e-7, 1e-2));
			
			PlotUtils.writePlots(outputDir, prefix, gp, 800, 750, true, true, false);
		}
		
		if (segModel != null) {
			HashSet<Jump> jumps = new HashSet<>();
			for (ClusterRupture cRup : origSol.getRupSet().requireModule(ClusterRuptures.class).getAll()) {
				for (Jump jump : cRup.getJumpsIterable()) {
					if (jump.fromSection.getSectionId() > jump.toSection.getSectionId())
						jump = jump.reverse();
					jumps.add(jump);
				}
			}
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			double maxSeg = 0d;
			if (segRate1 < 1d) {
				maxSeg = segRate1;
				funcs.add(new DefaultXY_DataSet(0, segRate1, segRate1, segRate1, segRate1, 0));
				chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.DARK_GRAY));
			}
			if (segRate2 < 1d && segRate2 != segRate1) {
				maxSeg = Math.max(maxSeg, segRate2);
				funcs.add(new DefaultXY_DataSet(0, segRate2, segRate2, segRate2, segRate2, 0));
				chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.DARK_GRAY));
			}
			
			for (boolean orig : new boolean[] {true,false}) {
				FaultSystemSolution sol = orig ? origSol : modSol;
				
				DefaultXY_DataSet xy = new DefaultXY_DataSet();
				for (Jump jump : jumps) {
					double rateSect1 = sol.calcTotParticRateForSect(jump.fromSection.getSectionId());
					double rateSect2 = sol.calcTotParticRateForSect(jump.toSection.getSectionId());
					HashSet<Integer> corups = new HashSet<>(sol.getRupSet().getRupturesForSection(jump.fromSection.getSectionId()));
					corups.retainAll(sol.getRupSet().getRupturesForSection(jump.toSection.getSectionId()));
					double corupRate = 0d;
					for (int rupIndex : corups)
						corupRate += sol.getRateForRup(rupIndex);
					xy.set(corupRate/rateSect1, corupRate/rateSect2);
				}
				
				xy.setName(orig ? "Original Recipe" : "MFD Inversion");
				funcs.add(xy);
				chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 6f, orig ? origSolColor : solColor));
			}
			
			HeadlessGraphPanel gp = PlotUtils.initScreenHeadless();
			
			PlotSpec plot = new PlotSpec(funcs, chars, "Segmentation", "Fract Rate 1", "Fract Rate 2");
			plot.setLegendInset(true);
			
			gp.drawGraphPanel(plot, false, false, new Range(0d, maxSeg*1.1), new Range(0d, maxSeg*1.1));
			
			PlotUtils.writePlots(outputDir, "seg_scatter", gp, 800, false, true, true, false);
		}
		
		ReportMetadata meta = new ReportMetadata(new RupSetMetadata("MFD Inversion Solution", modSol),
				new RupSetMetadata("NSHM23 Recipe Solution", origSol));
		ReportPageGen report = new ReportPageGen(meta, new File(outputDir, "report"),
				ReportPageGen.getDefaultSolutionPlots(PlotLevel.REVIEW));
		report.setReplot(true);
		report.generatePage();
	}
	
	private static List<? extends IncrementalMagFreqDist> getTargets(FaultSystemSolution sol, int parentID) {
		InversionTargetMFDs targets = sol.getRupSet().requireModule(InversionTargetMFDs.class);
		List<? extends IncrementalMagFreqDist> allTargets = targets.getOnFaultSupraSeisNucleationMFDs();
		if (parentID < 0)
			return allTargets;
		List<IncrementalMagFreqDist> ret = new ArrayList<>();
		for (FaultSection sect : sol.getRupSet().getFaultSectionDataList())
			if (sect.getParentSectionId() == parentID)
				ret.add(allTargets.get(sect.getSectionId()));
		return ret;
	}
	
	private static SummedMagFreqDist sum(List<? extends IncrementalMagFreqDist> mfds) {
		SummedMagFreqDist sum = null;
		for (IncrementalMagFreqDist mfd : mfds) {
			if (sum == null) {
				sum = new SummedMagFreqDist(mfd.getMinX(), mfd.size(), mfd.getDelta());
			} else {
				Preconditions.checkState(sum.size() == mfd.size());
				Preconditions.checkState(sum.getMinX() == mfd.getMinX());
			}
			sum.addIncrementalMagFreqDist(mfd);
		}
		return sum;
	}

}
