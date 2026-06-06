package scratch.kevin.mfdInversion;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
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

import com.google.common.base.Preconditions;

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
		double splitSlipRate1 = 2d;
		double splitSlipRate2 = 2d;
		double fractSlipSD = 0.1;
		
		int mainID = 0;
		int splitID1 = 1;
		int splitID2 = 2;
		
		Location origin = new Location(0d,0d);
		Location splitPoint = LocationUtils.location(origin, 0d, 50d);
		Location splitEnd1 = LocationUtils.location(splitPoint, -Math.PI/4d, 50d);
		Location splitEnd2 = LocationUtils.location(splitPoint, Math.PI/4d, 50d);
		
		// MFD params
		double b = 0.5d;
		double segRate1 = 1d;
		double segRate2 = 1d;
		SubSectConstraintModels sectConstrModel = SubSectConstraintModels.TOT_NUCL_RATE;
//		SubSectConstraintModels sectConstrModel = SubSectConstraintModels.NUCL_MFD;
		
		// misc
		int threads = 16;
		
		List<FaultSection> sects = new ArrayList<>();
		
		sects.add(new GeoJSONFaultSection.Builder(mainID, "Main Fault",
				FaultTrace.of(origin, splitPoint))
				.lowerDepth(lowerDepth)
				.upperDepth(upperDepth)
				.dip(dip)
				.rake(rake)
				.slipRate(mainSlipRate)
				.slipRateStdDev(mainSlipRate*fractSlipSD)
				.build());
		sects.add(new GeoJSONFaultSection.Builder(splitID1, "Split Fault 1",
				FaultTrace.of(splitPoint, splitEnd1))
				.lowerDepth(lowerDepth)
				.upperDepth(upperDepth)
				.dip(dip)
				.rake(rake)
				.slipRate(splitSlipRate1)
				.slipRateStdDev(splitSlipRate1*fractSlipSD)
				.build());
		sects.add(new GeoJSONFaultSection.Builder(splitID2, "Split Fault 2",
				FaultTrace.of(splitPoint, splitEnd2))
				.lowerDepth(lowerDepth)
				.upperDepth(upperDepth)
				.dip(dip)
				.rake(rake)
				.slipRate(splitSlipRate2)
				.slipRateStdDev(splitSlipRate2*fractSlipSD)
				.build());
		
		List<FaultSection> subSects = SubSectionBuilder.buildSubSects(sects);
		
		SimpleAzimuthalRupSetConfig rsConfig = new SimpleAzimuthalRupSetConfig(subSects, NSHM23_ScalingRelationships.AVERAGE);
		
		FaultSystemRupSet rupSet = rsConfig.build(1);
		
		long equivNumVars = Long.max(rupSet.getNumRuptures(), rupSet.getNumSections()*100l);
		
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
			if (orig) {
				origSol = sol;
				sol.write(new File(outputDir, "orig_sol.zip"));
			} else {
				modSol = sol;
				sol.write(new File(outputDir, "mod_sol.zip"));
			}
		}
		
		ReportMetadata meta = new ReportMetadata(new RupSetMetadata("MFD Inversion Solution", modSol),
				new RupSetMetadata("NSHM23 Recipe Solution", origSol));
		ReportPageGen report = new ReportPageGen(meta, new File(outputDir, "report"),
				ReportPageGen.getDefaultSolutionPlots(PlotLevel.REVIEW));
		report.setReplot(true);
		report.generatePage();
	}

}
