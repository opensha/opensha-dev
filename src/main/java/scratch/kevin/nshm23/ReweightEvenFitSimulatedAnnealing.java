package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.opensha.commons.util.ExceptionUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration.Builder;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionInputGenerator;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.ConstraintWeightingType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.PaleoRateInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.PaleoSlipInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.ParkfieldInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SectionTotalRateConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.ConstraintRange;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.SerialSimulatedAnnealing;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.SimulatedAnnealing;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.ThreadedSimulatedAnnealing;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.CompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.ProgressTrackingCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.TimeCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfitStats.MisfitStats;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfits;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.estimators.DraftModelConstraintBuilder;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;

import cern.colt.function.tdouble.IntIntDoubleFunction;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;

public class TestReweightInversion extends ThreadedSimulatedAnnealing {
	
	private DoubleMatrix2D origA, origA_ineq, modA, modA_ineq;
	private double[] origD, origD_ineq, modD, modD_ineq;
	private List<ConstraintRange> origRanges;
	
	private double maxAdjustmentFactor = 100d;

	public TestReweightInversion(DoubleMatrix2D A, double[] d, double[] initialState, double relativeSmoothnessWt,
			DoubleMatrix2D A_ineq, double[] d_ineq, int numThreads, CompletionCriteria subCompetionCriteria) {
		super(A, d, initialState, relativeSmoothnessWt, A_ineq, d_ineq, numThreads, subCompetionCriteria);
	}

	public TestReweightInversion(List<? extends SimulatedAnnealing> sas, CompletionCriteria subCompetionCriteria,
			boolean average) {
		super(sas, subCompetionCriteria, average);
	}

	public TestReweightInversion(ThreadedSimulatedAnnealing tsa) {
		super(tsa.getSAs(), tsa.getSubCompetionCriteria(), tsa.isAverage());
		setConstraintRanges(tsa.getConstraintRanges());
	}

//	private static final String targetName = "MAD";
//	private static double getTarget(MisfitStats stats) {
//		return stats.absMean;
//	}
	
	private static final String targetName = "RMSE";
//	private static final double minTargetForPenalty = 1d; // don't penalize anything for being "better" than this value
	
	private static final double avgTargetWeight2 = 2d; // don't mess with anything if the average is above this value
	private static final double avgTargetWeight1 = 1d; // linearly transition to targeted weights up to this average target
	
	private static double getTarget(MisfitStats stats) {
		return stats.rmse;
	}

	@Override
	protected void beforeRound(long curIter, int round) {
		if (round > 0) {
			Stopwatch watch = Stopwatch.createStarted();
			List<ConstraintRange> ranges = getConstraintRanges();
			Preconditions.checkNotNull(ranges);
			Preconditions.checkState(ranges.size() == origRanges.size());
			
			if (modA == null) {
				modA = origA.copy();
				modD = Arrays.copyOf(origD, origD.length);
				if (origA_ineq != null) {
					modA_ineq = origA_ineq.copy();
					modD_ineq = Arrays.copyOf(origD_ineq, origD_ineq.length);
				}
			}
			
			List<ConstraintRange> modRanges = new ArrayList<>();
			
			double[] misfits = getBestMisfit();
			double[] misfits_ineq = getBestInequalityMisfit();
			
			List<MisfitStats> stats = new ArrayList<>();
			double avgTarget = 0d;
			int num = 0;
			for (ConstraintRange range : ranges) {
				if (range.weightingType == ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY) {
					double[] myMisfits = range.inequality ? misfits_ineq : misfits;
					myMisfits = Arrays.copyOfRange(myMisfits, range.startRow, range.endRow);
					for (int i=0; i<myMisfits.length; i++)
						myMisfits[i] /= range.weight;
					MisfitStats myStats = new MisfitStats(myMisfits, range);
					avgTarget += getTarget(myStats);
					stats.add(myStats);
					num++;
				} else {
					stats.add(null);
				}
			}
			System.out.println("Readjusting weights for "+num+" uncertainty-weighted constraints");
			Preconditions.checkState(num > 0,
					"Can't use re-weighted inversion without any uncertainty-weighted constraints!");
			avgTarget /= (double)num;
			Preconditions.checkState(avgTarget > 0d && Double.isFinite(avgTarget), "Bad avg "+targetName+": %s", avgTarget);
			
			System.out.println("\tAverage misfit "+targetName+": "+(float)avgTarget);
//			if (avgTarget > minTargetForPenalty) {
//				System.out.println("\tAverage is above threshold, resetting to: "+(float)minTargetForPenalty);
//				avgTarget = minTargetForPenalty;
//			}
			for (int i=0; i<ranges.size(); i++) {
				MisfitStats myStats = stats.get(i);
				if (myStats == null) {
					modRanges.add(ranges.get(i));
					continue;
				}
				ConstraintRange range = ranges.get(i);
				double prevWeight = range.weight;
				double origWeight = origRanges.get(i).weight;
				double myTarget = getTarget(myStats);
				
				double misfitRatio = myTarget/avgTarget;
//				if (myTarget > minTargetForPenalty)
//					// don't downweight if we're above the threshold
//					misfitRatio = Math.max(misfitRatio, 1d);
				
//				double misfitRatio;
////				if (myTarget < avgTarget && myTarget > minTargetForPenalty)
////					// I'm better than average, but I'm still bad, don't penalize
////					misfitRatio = 1d;
////				else if (myTarget < avgTarget)
////					// I'm better than average and pretty good, compare me to lesser of the average and
////					// the target threshold
////					misfitRatio = myTarget/Math.min(avgTarget, minTargetForPenalty);
//				if (myTarget < avgTarget)
//					// I'm better than average, compare me to lesser of the average and the target threshold
//					misfitRatio = myTarget/Math.min(avgTarget, minTargetForPenalty);
//				else
//					misfitRatio = myTarget/avgTarget;
				
//				if (avgTarget > minTargetForPenalty) {
//					// we're poorly fit on average, don't overly penalize any that are well fit
//					if (myTarget < minTargetForPenalty)
//						// I'm quite well fit, don't encourage me to get worse than the threshold
//						misfitRatio = myTarget/minTargetForPenalty;
//					else
//						// don't penalize
//						misfitRatio = Math.max(1d, misfitRatio);
////					if (myTarget < avgTarget) {
////						if (myTarget < minTargetForPenalty)
////							// I'm quite well fit, don't encourage me to get worse than the threshold
////							misfitRatio = myTarget/minTargetForPenalty;
////						else
////							// I'm in-between well fit and the (poor) average, interpolate between targets
////							misfitRatio = 0.5*(misfitRatio + myTarget/minTargetForPenalty);
////					}
//				}
				
				double calcWeight = misfitRatio * prevWeight;
				double newWeight = Math.max(calcWeight, origWeight/maxAdjustmentFactor);
				newWeight = Math.min(newWeight, origWeight*maxAdjustmentFactor);
//				if (avgTarget > minTargetForPenalty) {
//					if (myTarget > avgTarget)
//						// we're worse than average, don't allow it to expand infinitely, however; bound with the
//						// original weight relative to threshold value
//						newWeight = Math.min(newWeight, origWeight*myTarget/minTargetForPenalty);
//					else
//						// we're better than average
//						newWeight = Math.max(newWeight, origWeight*myTarget/avgTarget);
//				}
//				if (myTarget > minTargetForPenalty)
//					newWeight = origWeight;
				
				System.out.println("\t"+range.shortName+":\t"+targetName+": "+(float)myTarget
						+";\tcalcWeight = "+(float)prevWeight+" x "+(float)misfitRatio+" = "+(float)calcWeight
						+";\tboundedWeight: "+(float)newWeight);
				
				if (avgTarget > avgTargetWeight2) {
					newWeight = origWeight;
					System.out.println("\t\tAbove max avg target, reverting to original weight: "+(float)origWeight);
				} else if (avgTarget > avgTargetWeight1) {
					double fract = (avgTarget - avgTargetWeight1)/(avgTargetWeight2 - avgTargetWeight1);
					Preconditions.checkState(fract >= 0d && fract <= 1d);
					newWeight = origWeight*fract + newWeight*(1-fract);
					System.out.println("\t\tAvg is poorly fit, linearly blending (fract="+(float)fract
							+") calculated weight with orig: "+(float)newWeight);
				}
				
				double scalar = newWeight / prevWeight;
				
				DoubleMatrix2D myA = range.inequality ? modA_ineq : modA;
				double[] myD = range.inequality ? modD_ineq : modD;
				for (int r=range.startRow; r<range.endRow; r++)
					myD[r] *= scalar;
				// could make this faster and only do once for all constraints
				myA.forEachNonZero(new IntIntDoubleFunction() {
					
					@Override
					public double apply(int row, int col, double val) {
						if (row >= range.startRow && row < range.endRow)
							return val*scalar;
						return val;
					}
				});
				
				// update the weight
				modRanges.add(new ConstraintRange(range.name, range.shortName, range.startRow, range.endRow,
						range.inequality, newWeight, range.weightingType));
			}
			
			System.out.println("Re-calculating misfits");
			double[] xbest = getBestSolution();
			double[] misfit = new double[modD.length];
			SerialSimulatedAnnealing.calculateMisfit(modA, modD, null, xbest, -1, Double.NaN, misfit);
			double[] misfit_ineq = null;
			int nRow = modA.rows();
			int nCol = modA.columns();
			int ineqRows = 0;
			if (modA_ineq != null) {
				misfit_ineq = new double[modD_ineq.length];
				ineqRows = misfit_ineq.length;
				SerialSimulatedAnnealing.calculateMisfit(modA_ineq, modD_ineq, null, xbest, -1, Double.NaN, misfit_ineq);
			}

			System.out.println("Re-calculating energies");
			double[] Ebest = SerialSimulatedAnnealing.calculateEnergy(xbest, misfit, misfit_ineq,
					nRow, nCol, ineqRows, ranges, 0d);
			
			setAll(modA, modD, modA_ineq, modD_ineq, Ebest, xbest, misfit, misfit_ineq, getNumNonZero());
			setConstraintRanges(modRanges);
			
			watch.stop();
			
			System.out.println("Took "+timeStr(watch.elapsed(TimeUnit.MILLISECONDS))+" to re-weight");
		}
		super.beforeRound(curIter, round);
	}

	@Override
	public long[] iterate(long startIter, long startPerturbs, CompletionCriteria criteria) {
		this.origA = getA();
		this.origA_ineq = getA_ineq();
		this.origD = getD();
		this.origD_ineq = getD_ineq();
		
		this.origRanges = getConstraintRanges();
		Preconditions.checkNotNull(origRanges, "Re-weigted inversion needs constraint ranges");
		this.origRanges = new ArrayList<>(origRanges);
		// make sure at least one uncert weighted
		boolean found = false;
		for (ConstraintRange range : origRanges) {
			if (range.weightingType == ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY) {
				found = true;
				break;
			}
		}
		Preconditions.checkState(found, "Must supply at least 1 uncertainty-weighted constraint for re-weighted inversion");
		
		long[] ret = super.iterate(startIter, startPerturbs, criteria);
		
		return ret;
	}
	
	public static void main(String[] args) {
		File parentDir = new File("/home/kevin/markdown/inversions");
		
		// run this if I need to attach UCERF3 modules after rebuilding default rupture sets
//		reprocessDefaultRupSets(parentDir);
//		System.exit(0);

		String dirName = new SimpleDateFormat("yyyy_MM_dd").format(new Date());

		dirName += "-coulomb-u3";
		File origRupSetFile = new File(parentDir, "fm3_1_u3ref_uniform_coulomb.zip");
		
		FaultSystemRupSet rupSet;
		try {
			rupSet = FaultSystemRupSet.load(origRupSetFile);
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		
//		rupSet = FaultSystemRupSet.buildFromExisting(rupSet)
//				.u3BranchModules(rupSet.getModule(U3LogicTreeBranch.class)).build();
		
		U3LogicTreeBranch branch = U3LogicTreeBranch.DEFAULT.copy();
		branch.setValue(DeformationModels.ABM);
		branch.setValue(ScalingRelationships.HANKS_BAKUN_08);
		rupSet = FaultSystemRupSet.buildFromExisting(rupSet).forU3Branch(branch).build();
		dirName += "-abm-hb08";
		
		double supraBVal = 0.8;
		dirName += "-nshm23_draft-supra_b_"+oDF.format(supraBVal);
		
		boolean applyDefModelUncertaintiesToNucl = true;
		boolean addSectCountUncertaintiesToMFD = false;
		boolean adjustForIncompatibleData = true;

		DraftModelConstraintBuilder constrBuilder = new DraftModelConstraintBuilder(rupSet, supraBVal,
				applyDefModelUncertaintiesToNucl, addSectCountUncertaintiesToMFD, adjustForIncompatibleData);
		constrBuilder.defaultConstraints();
		
		boolean reweight = true;
		dirName += "-reweight";
		
//		boolean reweight = false;

		CompletionCriteria completion = TimeCompletionCriteria.getInMinutes(30); dirName += "-30m-x1m";
		CompletionCriteria avgCompletion = TimeCompletionCriteria.getInMinutes(1);

//		CompletionCriteria completion = TimeCompletionCriteria.getInMinutes(30); dirName += "-30m-x5m";
//		CompletionCriteria avgCompletion = TimeCompletionCriteria.getInMinutes(5);
		
//		CompletionCriteria completion = TimeCompletionCriteria.getInHours(2); dirName += "-2h-x5m";
//		CompletionCriteria avgCompletion = TimeCompletionCriteria.getInMinutes(5);
		
//		CompletionCriteria completion = TimeCompletionCriteria.getInHours(2); dirName += "-2h-x1m";
//		CompletionCriteria avgCompletion = TimeCompletionCriteria.getInMinutes(1);
		
		Builder builder = InversionConfiguration.builder(constrBuilder.build(), completion);
		builder.sampler(constrBuilder.getSkipBelowMinSampler());
		builder.avgThreads(4, avgCompletion).threads(16);
		
//		builder.except(PaleoSlipInversionConstraint.class).except(PaleoRateInversionConstraint.class)
//			.except(ParkfieldInversionConstraint.class);
//		dirName += "-no_paleo_parkfield";
		
//		builder.except(SectionTotalRateConstraint.class);
//		dirName += "-no_sect";
		
		InversionConfiguration config = builder.build();
		
		InversionInputGenerator inputs = new InversionInputGenerator(rupSet, config);
		inputs.generateInputs(true);
		inputs.columnCompress();
		
		ProgressTrackingCompletionCriteria progress = new ProgressTrackingCompletionCriteria(completion);
		
		SimulatedAnnealing sa = config.buildSA(inputs);
		Preconditions.checkState(sa instanceof ThreadedSimulatedAnnealing);
		if (reweight)
			sa = new TestReweightInversion((ThreadedSimulatedAnnealing)sa);
		
		System.out.println("SA Parameters:");
		System.out.println("\tImplementation: "+sa.getClass().getName());
		System.out.println("\tCompletion Criteria: "+completion);
		System.out.println("\tPerturbation Function: "+sa.getPerturbationFunc());
		System.out.println("\tNon-Negativity Constraint: "+sa.getNonnegativeityConstraintAlgorithm());
		System.out.println("\tCooling Schedule: "+sa.getCoolingFunc());
		if (sa instanceof ThreadedSimulatedAnnealing) {
			ThreadedSimulatedAnnealing tsa = (ThreadedSimulatedAnnealing)sa;
			System.out.println("\tTop-Level Threads: "+tsa.getNumThreads());
			System.out.println("\tSub-Completion Criteria: "+tsa.getSubCompetionCriteria());
			System.out.println("\tAveraging? "+tsa.isAverage());
		}
		
		System.out.println("Annealing!");
		sa.iterate(progress);
		
		System.out.println("DONE. Building solution...");
		double[] rawSol = sa.getBestSolution();
		double[] rates = inputs.adjustSolutionForWaterLevel(rawSol);
		
		FaultSystemSolution sol = new FaultSystemSolution(rupSet, rates);
		// add inversion progress
		sol.addModule(progress.getProgress());
		sol.addModule(config);
		InversionMisfits misfits = new InversionMisfits(sa);
		sol.addModule(misfits);
		sol.addModule(misfits.getMisfitStats());
		
		File outputDir = new File(parentDir, dirName);
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		try {
			sol.write(new File(outputDir, "solution.zip"));
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
	}

	private static final DecimalFormat oDF = new DecimalFormat("0.##");

}
