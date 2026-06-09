package scratch.kevin.mfdInversion;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.IntegerSampler;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionInputGenerator;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.ConstraintWeightingType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.InversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.UncertainDataConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.JumpProbabilityConstraint.SectParticipationRateEstimator;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.ColumnOrganizedAnnealingData;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.ConstraintRange;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.SerialSimulatedAnnealing;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.SimulatedAnnealing;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.ThreadedSimulatedAnnealing;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.CompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.IterationCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.params.CoolingScheduleType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.params.GenerationFunctionType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.params.NonnegativityConstraintType;
import org.opensha.sha.earthquake.faultSysSolution.modules.AveSlipModule;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.modules.SlipAlongRuptureModel;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.FaultSubsectionCluster;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Jump;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Jump.UniqueDistJump;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.JumpProbabilityCalc;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RuptureTreeNavigator;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.estimators.GRParticRateEstimator;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.estimators.SectNucleationMFD_Estimator;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;

public class MFD_ScaledInversionAdjustment extends SectNucleationMFD_Estimator {
	
	private JumpProbabilityCalc segModel;
	
	private RupSetCoruptureMFDStructure structure;
	private double[] solution;
	private int threads;

	public MFD_ScaledInversionAdjustment(int threads, JumpProbabilityCalc segModel) {
		this.threads = threads;
		this.segModel = segModel;
	}

	@Override
	public boolean appliesTo(FaultSection sect) {
		return true;
	}

	@Override
	public void init(FaultSystemRupSet rupSet, List<IncrementalMagFreqDist> origSectSupraSeisMFDs,
			double[] targetSectSupraMoRates, double[] targetSectSupraSlipRates, double[] sectSupraSlipRateStdDevs,
			List<BitSet> sectRupUtilizations, int[] sectMinMagIndexes, int[] sectMaxMagIndexes,
			int[][] sectRupInBinCounts, EvenlyDiscretizedFunc refMFD) {
		boolean verbose = true;
		
		if (verbose) System.out.println("Building corupture structure");
		structure = new RupSetCoruptureMFDStructure(rupSet, origSectSupraSeisMFDs,
				targetSectSupraMoRates, targetSectSupraSlipRates, sectSupraSlipRateStdDevs,
				sectRupUtilizations, sectMinMagIndexes, sectMaxMagIndexes, sectRupInBinCounts, refMFD, segModel);
		
		if (verbose) System.out.println("Building constraints");
		List<InversionConstraint> constraints = new ArrayList<>();
		
		// slip rate constraint: ensure each nucleation MFD matches the original slip rate
		constraints.add(new MFD_SectSlipRateConstraint(structure, 1e2, ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY));
		// slip corupture balancing constraint: ensure slip in large multifault rupture mag bins is available on all faults
		constraints.add(new SectCoruptureBudgetConstraint(structure, 1e1, ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY));
		// scale factor normalization: ensure single fault scale factors are >=1
		constraints.add(new ScaleFactorLimitConstraint(structure, true, 1e4));
		// scale factor normalization: ensure multi fault scale factors are <=1
		constraints.add(new ScaleFactorLimitConstraint(structure, false, 1e4));
		// scale factor normalization: keep-near 1 (weakly)
		constraints.add(new ScaleFactorOneConstraint(structure, 1e0, 5e0));
//		// Minimzation constraint: week constraint to force it to use large magnitude bins if they don't break the other constraints
//		constraints.add(new SectRateMinimizationConstraint(structure, 1e1, ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY)););
		if (segModel != null)
			constraints.add(new MFD_SegmentationConstraint(structure, 5e1));
		
		if (verbose) System.out.println("Building inversion inputs");
		List<ConstraintRange> constraintRowRanges = InversionInputGenerator.buildConstraintRanges(constraints, verbose);
		int numRows = 0;
		int numIneqRows = 0;
		for (ConstraintRange range : constraintRowRanges) {
			if (range.inequality)
				numIneqRows = range.endRow;
			else
				numRows = range.endRow;
		}
		
		int columns = structure.getNumColumns();
		
		DoubleMatrix2D A = null;
		double[] d = null;
		
		if (numRows > 0) {
			A = new SparseDoubleMatrix2D(numRows, columns);
			d = new double[numRows];
		}
		
		DoubleMatrix2D A_ineq = null;
		double[] d_ineq = null;
		
		if (numIneqRows > 0) {
			A_ineq = new SparseDoubleMatrix2D(numIneqRows, columns);
			d_ineq = new double[numIneqRows];
		}
		
		if (verbose) System.out.println("Encoding matrices");
		
		Stopwatch watch = verbose ? Stopwatch.createStarted() : null;
		int numNonZero = 0;
		
		DecimalFormat oneDigit = new DecimalFormat("0.0");
		
		for (int i=0; i<constraints.size(); i++) {
			InversionConstraint constraint = constraints.get(i);
			ConstraintRange rowRange = constraintRowRanges.get(i);
			
			DoubleMatrix2D myA;
			double[] myD;
			if (constraint.isInequality()) {
				myA = A_ineq;
				myD = d_ineq;
			} else {
				myA = A;
				myD = d;
			}
			
			if (verbose)
				System.out.println("\tEncoding "+constraint.getName()
					+", ineq="+constraint.isInequality());
			Stopwatch subWatch = verbose ? Stopwatch.createStarted() : null;
			long myNonZero = constraint.encode(myA, myD, rowRange.startRow);
			if (verbose) {
				long maxNum = (rowRange.endRow - rowRange.startRow)*(long)columns;
				double density = 100d*(double)myNonZero/(double)maxNum;
				System.out.println("\t\tDONE, took "+InversionInputGenerator.getTimeStr(subWatch)+" to encode "
						+myNonZero+" values (density: "+oneDigit.format(density)+" %)");
				subWatch.stop();
			}
			numNonZero += myNonZero;
		}
		
		if (verbose) {
			long maxNum = (numRows+numIneqRows)*(long)columns;
			double density = 100d*(double)numNonZero/(double)maxNum;
			System.out.println("DONE encoding, took "+InversionInputGenerator.getTimeStr(watch)+" to encode "
					+numNonZero+" values (density: "+oneDigit.format(density)+" %)");
			watch.stop();
		}
		
//		for (ConstraintRange range : constraintRowRanges) {
//			DoubleMatrix2D myA = range.inequality ? A_ineq : A;
//			double[] myD = range.inequality ? d_ineq : d;
//			System.out.println("Sol=1 debug for "+range.name);
//			System.out.println("=== A ===");
//			printA(myA, range);
//			System.out.println("=== d ===");
//			printD(myD, range);
//			System.out.println("=== misfits ===");
//			misfitsDebugFor1s(myA, myD, range);
//		}
//		System.exit(0);
		
		// run the inversion
		ColumnOrganizedAnnealingData equalityData = numRows > 0 ? new ColumnOrganizedAnnealingData(A, d) : null;
		ColumnOrganizedAnnealingData inequalityData = numIneqRows > 0 ? new ColumnOrganizedAnnealingData(A_ineq, d_ineq) : null;
		
		double[] initial = new double[columns];
		
		GenerationFunctionType perturb = GenerationFunctionType.EXPONENTIAL_SCALE;
		NonnegativityConstraintType nonneg = NonnegativityConstraintType.TRY_ZERO_RATES_OFTEN;
//		GenerationFunctionType perturb = GenerationFunctionType.UNIFORM_0p0001;
//		NonnegativityConstraintType nonneg = NonnegativityConstraintType.LIMIT_ZERO_RATES;
//		long itersPerVal = 100000l;
//		CoolingScheduleType cool = CoolingScheduleType.FAST_SA;
		long itersPerVal = 100000l;
		CoolingScheduleType cool = CoolingScheduleType.CLASSICAL_SA;
		
		long totalIters = itersPerVal * columns;
		CompletionCriteria completion = new IterationCompletionCriteria(totalIters);
		
//		threads = 1;
		SimulatedAnnealing sa;
		if (threads > 1) {
			
			int avgThreads = Integer.max(threads/4, 2);
			
			int threadsPerAvg = (int)Math.ceil((double)threads/(double)avgThreads);
			Preconditions.checkState(threadsPerAvg <= threads);
			Preconditions.checkState(threadsPerAvg > 0);
			
			CompletionCriteria avgCompletion = new IterationCompletionCriteria(totalIters/10);
			CompletionCriteria subCompletion = new IterationCompletionCriteria(totalIters/100);
			
			int threadsLeft = threads;
			
			// arrange lower-level (actual worker) SAs
			List<SimulatedAnnealing> tsas = new ArrayList<>();
			while (threadsLeft > 0) {
				int myThreads = Integer.min(threadsLeft, threadsPerAvg);
				if (myThreads > 1)
					tsas.add(new ThreadedSimulatedAnnealing(equalityData, inequalityData,
							initial, 0d, myThreads, subCompletion));
				else
					tsas.add(new SerialSimulatedAnnealing(equalityData, inequalityData,
							initial, 0d));
				threadsLeft -= myThreads;
			}
			sa = new ThreadedSimulatedAnnealing(tsas, avgCompletion, true);
		} else {
			sa = new SerialSimulatedAnnealing(equalityData, inequalityData, initial, 0d);
		}
		sa.setConstraintRanges(constraintRowRanges);
		
		sa.setPerturbationFunc(perturb);
		sa.setNonnegativeityConstraintAlgorithm(nonneg);
		sa.setCoolingFunc(cool);
		
		sa.iterate(completion);
		
		solution = sa.getBestSolution();
		
		if (verbose) {
			System.out.println("DONE estimating MFDs");
			DecimalFormat magDF = new DecimalFormat("0.0");
			DecimalFormat fractDF = new DecimalFormat("0.0000");
			for (int s=0; s<rupSet.getNumSections(); s++) {
				StringBuilder line = new StringBuilder();
				line.append(s).append(". ").append(rupSet.getFaultSectionData(s).getSectionName()).append(":");
				List<SectCommonPathwaysMagRange> pathways = structure.getSectCommonPathways(s);
				int colStart = structure.getSectStartColumn(s);
				for (int p=0; p<pathways.size(); p++) {
					SectCommonPathwaysMagRange path = pathways.get(p);
					line.append("    M[").append(magDF.format(refMFD.getX(path.minMagIndex))).append("-")
						.append(magDF.format(refMFD.getX(path.maxMagIndex)));
					line.append("]=").append(fractDF.format(solution[colStart+p]));
				}
				System.out.println(line);
			}
		}
		
//		for (ConstraintRange range : constraintRowRanges) {
//			DoubleMatrix2D myA = range.inequality ? A_ineq : A;
//			double[] myD = range.inequality ? d_ineq : d;
//			System.out.println("Sol=1 debug for "+range.name);
//			System.out.println("=== A ===");
//			printA(myA, range);
//			System.out.println("=== d ===");
//			printD(myD, range);
//			System.out.println("=== misfits ===");
//			misfitsDebug(myA, myD, solution, range);
//		}
//		System.exit(0);
	}

	@Override
	public IncrementalMagFreqDist estimateNuclMFD(FaultSection sect, IncrementalMagFreqDist curSectSupraSeisMFD,
			List<Integer> availableRupIndexes, List<Double> availableRupMags, UncertainDataConstraint sectMomentRate,
			boolean sparseGR) {
		EvenlyDiscretizedFunc refMFD = structure.getRefMFD();
		int sectIndex = sect.getSectionId();
		IncrementalMagFreqDist origMFD = structure.getOrigSectSupraSeisMFDs().get(sectIndex);
		IncrementalMagFreqDist mfd = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		List<SectCommonPathwaysMagRange> pathways = structure.getSectCommonPathways(sectIndex);
		int sectStartCol = structure.getSectStartColumn(sectIndex);
		for (int i=0; i<pathways.size(); i++) {
			SectCommonPathwaysMagRange pathway = pathways.get(i);
			double scalar = solution[sectStartCol + i];
			Preconditions.checkState(Double.isFinite(scalar));
			for (int m=pathway.minMagIndex; m<=pathway.maxMagIndex; m++)
				mfd.set(m, origMFD.getY(m)*scalar);
		}
		return mfd;
	}
	
	private static class RupSetCoruptureMFDStructure {
		
		private Map<Integer, List<FaultSection>> parentSectsMap;
		
		private BitSet includedRups;
		private List<Integer> includedRupIndexes;
		
		// column organization
		private int numCols;
		private int[] sectColIndexes;
		private int[] sectMinMagIndexes;
		private int[] sectMaxMagIndexes;
		
		// rupture pathways
		private List<List<SectMagRupturePathways>> sectRupMagPathways;
		private List<List<SectCommonPathwaysMagRange>> sectCommonPathways;
		private int[] sectMaxSingleFaultMags;

		private FaultSystemRupSet rupSet;

		private List<IncrementalMagFreqDist> origSectSupraSeisMFDs;

		private double[] targetSectSupraMoRates;
		private double[] targetSectSupraSlipRates;
		private double[] sectSupraSlipRateStdDevs;

		private List<BitSet> sectRupUtilizations;

		private int[][] sectRupInBinCounts;

		private EvenlyDiscretizedFunc refMFD;

		private JumpProbabilityCalc segModel;
		
		public RupSetCoruptureMFDStructure(FaultSystemRupSet rupSet, List<IncrementalMagFreqDist> origSectSupraSeisMFDs,
				double[] targetSectSupraMoRates, double[] targetSectSupraSlipRates, double[] sectSupraSlipRateStdDevs,
				List<BitSet> sectRupUtilizations, int[] sectMinMagIndexes, int[] sectMaxMagIndexes,
				int[][] sectRupInBinCounts, EvenlyDiscretizedFunc refMFD, JumpProbabilityCalc segModel) {
			this.rupSet = rupSet;
			
			int numRuptures = rupSet.getNumRuptures();
			int numSections = rupSet.getNumSections();
			
			for (int s=0; s<numSections; s++) {
				IncrementalMagFreqDist mfd = origSectSupraSeisMFDs.get(s);
				Preconditions.checkNotNull(mfd);
				Preconditions.checkState(EvenlyDiscretizedFunc.areXValuesIdentical(mfd, refMFD));
				Preconditions.checkState(sectMinMagIndexes[s] >= 0);
				Preconditions.checkState(sectMaxMagIndexes[s] >= sectMinMagIndexes[s]);
			}
			this.origSectSupraSeisMFDs = origSectSupraSeisMFDs;
			this.targetSectSupraMoRates = targetSectSupraMoRates;
			this.targetSectSupraSlipRates = targetSectSupraSlipRates;
			this.sectSupraSlipRateStdDevs = sectSupraSlipRateStdDevs;
			this.sectRupUtilizations = sectRupUtilizations;
			this.sectRupInBinCounts = sectRupInBinCounts;
			this.refMFD = refMFD;
			this.segModel = segModel;
			List<? extends FaultSection> subSects = rupSet.getFaultSectionDataList();
			
			// bin sections by parent
			parentSectsMap = new HashMap<>();
			for (FaultSection sect : subSects) {
				int parentID = sect.getParentSectionId();
				Preconditions.checkState(parentID >= 0, "ParentIDs are required");
				int subSectIndex = sect.getSubSectionIndex();
				Preconditions.checkState(subSectIndex >= 0, "Sub-section indexes are required");
				if (!parentSectsMap.containsKey(parentID))
					parentSectsMap.put(parentID, new ArrayList<>());
				List<FaultSection> parentSects = parentSectsMap.get(parentID);
				while (parentSects.size() <= subSectIndex)
					parentSects.add(null);
				Preconditions.checkState(parentSects.get(subSectIndex) == null);
				parentSects.set(subSectIndex, sect);
			}
			
			// ruptures we're considering (may not be all if some have been explicitly filtered out beforehand)
			includedRups = new BitSet(numRuptures);
			for (BitSet sectRups : sectRupUtilizations)
				includedRups.or(sectRups);
			
			this.sectMinMagIndexes = sectMinMagIndexes;
			this.sectMaxMagIndexes = sectMaxMagIndexes;
			
			ClusterRuptures cRups = segModel == null ? null : rupSet.requireModule(ClusterRuptures.class);
			
			// set of parent sections involved in each rupture
			List<Set<Integer>> rupParentSets = new ArrayList<>(numRuptures);
			includedRupIndexes = new ArrayList<>(numRuptures);
			for (int rupIndex=0; rupIndex<numRuptures; rupIndex++) {
				if (includedRups.get(rupIndex)) {
					Set<Integer> parents = new HashSet<>();
					for (FaultSection sect : rupSet.getFaultSectionDataForRupture(rupIndex))
						parents.add(sect.getParentSectionId());
					rupParentSets.add(parents);
					includedRupIndexes.add(rupIndex);
				} else {
					rupParentSets.add(null);
				}
			}
			
			sectRupMagPathways = new ArrayList<>(numSections);
			sectCommonPathways = new ArrayList<>(numSections);
			sectMaxSingleFaultMags = new int[numSections];
			for (int s=0; s<numSections; s++) {
				BitSet sectRups = sectRupUtilizations.get(s);
				sectMaxSingleFaultMags[s] = -1;
				
				int myNumMags = 1 + sectMaxMagIndexes[s] - sectMinMagIndexes[s];
				// ruptures on this section, grouped by magnitude bin
				List<List<Integer>> rupsByMagIndex = new ArrayList<>(myNumMags);
				for (int m=0; m<myNumMags; m++)
					rupsByMagIndex.add(null);
				
				int myParentID = subSects.get(s).getParentSectionId();
				
				for (int rupIndex : includedRupIndexes) {
					if (!sectRups.get(rupIndex))
						continue;
					
					int magIndex = refMFD.getClosestXIndex(rupSet.getMagForRup(rupIndex)) - sectMinMagIndexes[s];
					List<Integer> rupsForMag = rupsByMagIndex.get(magIndex);
					if (rupsForMag == null) {
						rupsForMag = new ArrayList<>();
						rupsByMagIndex.set(magIndex, rupsForMag);
					}
					rupsForMag.add(rupIndex);
				}

				List<SectMagRupturePathways> myPathways = new ArrayList<>(myNumMags);
				sectRupMagPathways.add(myPathways);
				
				int curCommonPathwayMmin = sectMinMagIndexes[s];
				int curCommonPathwayMmax = -1;
				Set<Set<Integer>> curCommonPathways = null;
				List<SectCommonPathwaysMagRange> myCommonPathways = new ArrayList<>(Integer.min(4, myNumMags));
				sectCommonPathways.add(myCommonPathways);
				
				for (int m=0; m<myNumMags; m++) {
					List<Integer> rupsForMag = rupsByMagIndex.get(m);
					if (rupsForMag == null) {
						// no ruptures at this magnitude
						myPathways.add(null);
						continue;
					}
					int magIndex = sectMinMagIndexes[m] + m;
					
					// figure out each unique parent combination/pathway
					// for each such pathway, keep a list of weights for each rupture (to determine the overall path weight)
					Map<Set<Integer>, List<Double>> uniqueParentCombWeights = new HashMap<>();
					// also keep track of the ruptures using those pathways (for seg constraints)
					Map<Set<Integer>, List<Integer>> uniqueParentCombRups = new HashMap<>();
					// also keep track of the sections used by those parents for weighting later
					Map<Integer, int[]> parentUsedSectCounts = new HashMap<>();
					boolean anySingleFault = false;
					for (int rupIndex : rupsForMag) {
						Set<Integer> parents = rupParentSets.get(rupIndex);
						if (parents.size() == 1) {
							anySingleFault = true;
							break;
						}
						List<FaultSection> rupSects = rupSet.getFaultSectionDataForRupture(rupIndex);
						// figure out average slip rate for each parent
						Map<Integer, DoubleAverager> parentSlips = new HashMap<>(parents.size());
						for (FaultSection sect : rupSects) {
							int parentID = sect.getParentSectionId();
							if (!parentSlips.containsKey(parentID))
								parentSlips.put(parentID, new DoubleAverager());
							parentSlips.get(parentID).add(targetSectSupraSlipRates[sect.getSectionId()]);
							if (!parentUsedSectCounts.containsKey(parentID))
								parentUsedSectCounts.put(parentID, new int[parentSectsMap.get(parentID).size()]);
							int[] sectsUsed = parentUsedSectCounts.get(parentID);
							sectsUsed[sect.getSubSectionIndex()]++;
						}

						// calculate weight for this rupture, based on lowest segmentation-adjusted slip rate
						double rupWeight = Double.POSITIVE_INFINITY;
						if (segModel == null) {
							// simple, just use smallest slip rate as weight
							for (int parentID : parents) {
								if (parentID == myParentID)
									continue;
								
								rupWeight = Math.min(rupWeight, parentSlips.get(parentID).getAverage());
							}
						} else {
							// also incorporate segmentation rates
							ClusterRupture cRup = cRups.get(rupIndex);
							RuptureTreeNavigator rupNav = cRup.getTreeNavigator();
							
							for (FaultSubsectionCluster cluster : cRup.getClustersIterable()) {
								
								if (cluster.parentSectionID != myParentID) {
									double slipRate  = parentSlips.get(cluster.parentSectionID).getAverage();
									Jump jumpTo = rupNav.getJumpTo(cluster);
									List<Jump> jumpsInvolving = new ArrayList<>(2);
									if (jumpTo != null)
										jumpsInvolving.add(jumpTo);
									for (FaultSubsectionCluster oCluster : rupNav.getDescendants(cluster))
										jumpsInvolving.add(rupNav.getJumpTo(oCluster));
									
									rupWeight = Math.min(rupWeight, slipRate);
									
									for (Jump jump : jumpsInvolving) {
										double segRate = segModel.calcJumpProbability(cRup, jump, false);
										if (segRate < 1d)
											rupWeight = Math.min(rupWeight, segRate*slipRate);
									}
								}
							}
						}
						Preconditions.checkState(Double.isFinite(rupWeight));
						List<Double> pathWeights = uniqueParentCombWeights.get(parents);
						if (pathWeights == null) {
							pathWeights = new ArrayList<>();
							uniqueParentCombWeights.put(parents, pathWeights);
							uniqueParentCombRups.put(parents, new ArrayList<>());
						}
						pathWeights.add(rupWeight);
						uniqueParentCombRups.get(parents).add(rupIndex);
					}
					if (anySingleFault) {
						sectMaxSingleFaultMags[s] = magIndex;
						curCommonPathwayMmax = magIndex;
						// TODO here
						myPathways.add(null);
						// we have single-fault ruptures at this mag, skip the constraint
						continue;
					}
					
					// determine aggregate usage weights for each parent
					Map<Integer, Double> parentWeights = new HashMap<>(parentUsedSectCounts.size()-1);
					double sumPathWeights = 0d;
					List<Set<Integer>> parentPathways = new ArrayList<>(uniqueParentCombWeights.size());
					List<Double> pathwayWeights = new ArrayList<>(uniqueParentCombWeights.size());
					List<int[]> pathwayRups = new ArrayList<>(uniqueParentCombWeights.size());
					List<double[]> pathwayRupWeights = new ArrayList<>(uniqueParentCombWeights.size());
					for (Set<Integer> parents : uniqueParentCombWeights.keySet()) {
						parentPathways.add(parents);
						int[] pathRups = Ints.toArray(uniqueParentCombRups.get(parents));
						double[] pathRupWeights = Doubles.toArray(uniqueParentCombWeights.get(parents));
						pathwayRups.add(pathRups);
						pathwayRupWeights.add(pathRupWeights);
						double pathWeight = StatUtils.sum(pathRupWeights);
						pathwayWeights.add(pathWeight);
						sumPathWeights += pathWeight;
						for (Integer parentID : parents) {
							if (parentID == myParentID)
								continue;
							Double prev = parentWeights.get(parentID);
							if (prev == null)
								prev = 0d;
							parentWeights.put(parentID, prev+pathWeight);
						}
					}
					Preconditions.checkState(sumPathWeights > 0d);
					// normalize
					for (int i=0; i<pathwayWeights.size(); i++)
						pathwayWeights.set(i, pathwayWeights.get(i)/sumPathWeights);
					for (Integer pathParentID : parentWeights.keySet())
						parentWeights.put(pathParentID, parentWeights.get(pathParentID)/sumPathWeights);
					for (double[] weights : pathwayRupWeights) {
						for (int i=0; i<weights.length; i++)
							weights[i] /= sumPathWeights;
					}
					
					myPathways.add(new SectMagRupturePathways(s, sectMinMagIndexes[s]+m, parentPathways, pathwayWeights,
							pathwayRups, pathwayRupWeights, parentWeights, parentUsedSectCounts));
					
					// now see if this magnitude bin used the same set of pathways
					if (curCommonPathways == null || !curCommonPathways.equals(uniqueParentCombWeights.keySet())) {
						// new set of pathways
						
						// store the previous one
						Preconditions.checkState(curCommonPathwayMmax >= curCommonPathwayMmin);
						Preconditions.checkState(curCommonPathwayMmax >= 0);
						myCommonPathways.add(new SectCommonPathwaysMagRange(s, curCommonPathwayMmin, curCommonPathwayMmax, curCommonPathways));
						
						curCommonPathways = uniqueParentCombWeights.keySet();
						curCommonPathwayMmin = magIndex;
						curCommonPathwayMmax = magIndex;
					} else {
						// continuation
						curCommonPathwayMmax = magIndex;
					}
				}
				
				// store the last set of common pathways
				Preconditions.checkState(curCommonPathwayMmax >= curCommonPathwayMmin);
				Preconditions.checkState(curCommonPathwayMmax >= 0);
				myCommonPathways.add(new SectCommonPathwaysMagRange(s, curCommonPathwayMmin, curCommonPathwayMmax, curCommonPathways));
			}
			Preconditions.checkState(sectRupMagPathways.size() == numSections);
			Preconditions.checkState(sectCommonPathways.size() == numSections);
			// determine column indexes of the A matrix, one for each section-pathway combination
			sectColIndexes = new int[numSections];
			numCols = 0;
			for (int s=0; s<numSections; s++) {
				List<SectCommonPathwaysMagRange> pathways = sectCommonPathways.get(s);
				Preconditions.checkNotNull(pathways);
				Preconditions.checkState(!pathways.isEmpty());
				sectColIndexes[s] = numCols;
				numCols += pathways.size();
			}
		}
		
		public int getSectStartColumn(int sectIndex) {
			return sectColIndexes[sectIndex];
		}
		
		public int locateSectMagColumn(int sectIndex, int magIndex) {
			List<SectCommonPathwaysMagRange> pathways = getSectCommonPathways(sectIndex);
			for (int p=0; p<pathways.size(); p++) {
				SectCommonPathwaysMagRange pathway = pathways.get(p);
				if (magIndex >= pathway.minMagIndex && magIndex <= pathway.maxMagIndex)
					return sectColIndexes[sectIndex] + p;
			}
			throw new IllegalStateException("No pathway found for sectIndex="+sectIndex+" containing magIndex="+magIndex);
		}
		
		public double getSectOrigRate(int sectIndex, int magIndex) {
			return origSectSupraSeisMFDs.get(sectIndex).getY(magIndex);
		}
		
		public List<SectCommonPathwaysMagRange> getSectCommonPathways(int sectIndex) {
			return sectCommonPathways.get(sectIndex);
		}
		
		public int getNumColumns() {
			return numCols;
		}
		
		public SectMagRupturePathways getMagBinPathways(int sectIndex, int magIndex) {
			return sectRupMagPathways.get(sectIndex).get(magIndex-sectMinMagIndexes[sectIndex]);
		}

		public FaultSystemRupSet getRupSet() {
			return rupSet;
		}

		public double[] getTargetSectSupraMoRates() {
			return targetSectSupraMoRates;
		}

		public double[] getTargetSectSupraSlipRates() {
			return targetSectSupraSlipRates;
		}

		public double[] getSectSupraSlipRateStdDevs() {
			return sectSupraSlipRateStdDevs;
		}

		public EvenlyDiscretizedFunc getRefMFD() {
			return refMFD;
		}
		
		public BitSet getSectRupUtilizations(int sectIndex) {
			return sectRupUtilizations.get(sectIndex);
		}
		
		public List<Integer> getIncludedRupIndexes() {
			return includedRupIndexes;
		}

		public int[] getSectMinMagIndexes() {
			return sectMinMagIndexes;
		}

		public int[] getSectMaxMagIndexes() {
			return sectMaxMagIndexes;
		}
		
		public List<FaultSection> getSectsForParent(int parentID) {
			return parentSectsMap.get(parentID);
		}

		public List<IncrementalMagFreqDist> getOrigSectSupraSeisMFDs() {
			return origSectSupraSeisMFDs;
		}
		
		public double getMultifaultRupFractContributionToBin(int sectIndex, int magIndex, int rupIndex) {
			SectMagRupturePathways sectPaths = getMagBinPathways(sectIndex, magIndex);
			if (sectPaths == null)
				return 0d;
			for (int i=0; i<sectPaths.parentPathways.size(); i++) {
				int[] rups = sectPaths.pathwayRups.get(i);
				double[] weights = sectPaths.pathwayRupWeights.get(i);
				for (int j=0; j<rups.length; j++)
					if (rups[j] == rupIndex)
						return weights[j];
			}
			return 0d;
		}
		
		public int getSectMaxSingleFaultMagIndex(int sectIndex) {
			return sectMaxSingleFaultMags[sectIndex];
		}
		
		public JumpProbabilityCalc getSegModel() {
			return segModel;
		}
	}
	
	private static record SectMagRupturePathways(int sectIndex, int magIndex, List<Set<Integer>> parentPathways,
			List<Double> pathwayWeights, List<int[]> pathwayRups, List<double[]> pathwayRupWeights,
			Map<Integer, Double> parentWeights, Map<Integer, int[]> parentSectCounts) {};
	
	private static record SectCommonPathwaysMagRange(int sectIndex, int minMagIndex, int maxMagIndex, Set<Set<Integer>> parentPathways) {};
	
	/**
	 * This constraint ensures that each section's individual MFD sums up to match the target slip rate
	 */
	private static class MFD_SectSlipRateConstraint extends InversionConstraint {
		private RupSetCoruptureMFDStructure structure;

		public MFD_SectSlipRateConstraint(RupSetCoruptureMFDStructure structure,
				double weight, ConstraintWeightingType weightingType) {
			super("Section Slip Rate Constraint", "SectSlipRate", weight, false, weightingType);
			this.structure = structure;
		}

		@Override
		public int getNumRows() {
			return structure.getRupSet().getNumSections();
		}

		@Override
		public long encode(DoubleMatrix2D A, double[] d, int startRow) {
			FaultSystemRupSet rupSet = structure.getRupSet();
			
			int numSections = rupSet.getNumSections();
			Preconditions.checkState(A.columns() == structure.getNumColumns());
			
			double[] targetSectSupraSlipRates = structure.getTargetSectSupraSlipRates();
			double[] sectSupraSlipRateStdDevs = structure.getSectSupraSlipRateStdDevs();
			
			double[] weights = new double[numSections];
			for (int s=0; s<numSections; s++)
				weights[s] = this.weight;
			if (weightingType == ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY) {
				for (int s=0; s<numSections; s++)
					if (sectSupraSlipRateStdDevs[s] != 0d)
						weights[s] /= sectSupraSlipRateStdDevs[s];
			}
			
			AveSlipModule aveSlipModule = rupSet.requireModule(AveSlipModule.class);
			SlipAlongRuptureModel slipAlongModule = rupSet.requireModule(SlipAlongRuptureModel.class);
			
			EvenlyDiscretizedFunc refMFD = structure.getRefMFD();
			
			int[] sectMinMagIndexes = structure.getSectMinMagIndexes();
			int[] sectMaxMagIndexes = structure.getSectMaxMagIndexes();
			
			// figure out section slips
			int[][] sectMagRupCounts = new int[numSections][];
			double[][] sectMagSlipConsumptions = new double[numSections][];
			for (int s=0; s<numSections; s++) {
				int numMag = 1 + sectMaxMagIndexes[s] - sectMinMagIndexes[s];
				sectMagRupCounts[s] = new int[numMag];
				sectMagSlipConsumptions[s] = new double[numMag];
			}
			for (int rupIndex : structure.getIncludedRupIndexes()) {
				double[] slips = slipAlongModule.calcSlipOnSectionsForRup(rupSet, aveSlipModule, rupIndex);
				List<Integer> sects = rupSet.getSectionsIndicesForRup(rupIndex);
				Preconditions.checkState(sects.size() == slips.length);
				int magIndex = refMFD.getClosestXIndex(rupSet.getMagForRup(rupIndex));
				double rupArea = rupSet.getAreaForRup(rupIndex);
				for (int s=0; s<slips.length; s++) {
					int sectID = sects.get(s);
					int m = magIndex - sectMinMagIndexes[sectID];
					sectMagRupCounts[sectID][m]++;
					// we're comparing to nucleation MFDs, but slip rate scales with participation
					// this scalar fixes that
					double particToNuclScalar = rupArea/rupSet.getAreaForSection(s);
					sectMagSlipConsumptions[sectID][m] += slips[s]*particToNuclScalar;
				}
			}
			// normalize
			for (int s=0; s<numSections; s++)
				for (int m=0; m<sectMagRupCounts[s].length; m++)
					sectMagSlipConsumptions[s][m] /= (double)sectMagRupCounts[s][m];
		
			long numNonZeroElements = 0l;
			
			// A matrix component of slip-rate constraint
			for (int s=0; s<numSections; s++) {
				int row = startRow + s;
				int colStartIndex = structure.getSectStartColumn(s);
				List<SectCommonPathwaysMagRange> pathways = structure.getSectCommonPathways(s);
//				System.out.println(row+". Debug for slip rate on sect "+s+" ("+rupSet.getFaultSectionData(s).getName()
//						+") w/ slipRate="+(float)targetSectSupraSlipRates[s]);
				for (int p=0; p<pathways.size(); p++) {
					SectCommonPathwaysMagRange pathway = pathways.get(p);
					int col = colStartIndex+p;
//					System.out.println("\tPathway "+p+", Mags=["+(float)refMFD.getX(pathway.minMagIndex)+","
//							+(float)refMFD.getX(pathway.maxMagIndex)+"; parents="+pathway.parentPathways);
					for (int magIndex=pathway.minMagIndex; magIndex<=pathway.maxMagIndex; magIndex++) {
						int m = magIndex - sectMinMagIndexes[s];
						if (sectMagRupCounts[s][m] > 0) {
							double avgSlipConsumption = sectMagSlipConsumptions[s][m];
							double origRate = structure.getSectOrigRate(s, magIndex);
							double val = avgSlipConsumption * origRate;
//							System.out.println("\t\tM"+(float)refMFD.getX(magIndex)+"; avgConsump="
//									+(float)avgSlipConsumption+"\torigRate="+(float)origRate);
							
							if (weightingType == ConstraintWeightingType.NORMALIZED) {
								double target = targetSectSupraSlipRates[s];
								if (target != 0d) {
									// Note that constraints for sections w/ slip rate < 0.1 mm/yr is not normalized by slip rate
									// -- otherwise misfit will be huge (e.g., UCERF3 GEOBOUND model has 10e-13 slip rates that will
									// dominate misfit otherwise)
									if (target < 1e-4 || Double.isNaN(target))
										target = 1e-4;
									val /= target;
								}
							}
							if (!addA(A, row, col, val*weights[s]))
								numNonZeroElements++;
						}
					}
				}
			}
			
			// d vector component of slip-rate constraint
			for (int s=0; s<numSections; s++) {
				double target = targetSectSupraSlipRates[s];
				double val = target;
				if (weightingType == ConstraintWeightingType.NORMALIZED) {
					if (target == 0d)
						// minimize
						val = 0d;
					else if (target < 1E-4 || Double.isNaN(target))
						// For very small slip rates, do not normalize by slip rate
						//  (normalize by 0.0001 instead) so they don't dominate misfit
						val = targetSectSupraSlipRates[s]/1e-4;
					else
						val = 1d;
				}
				int row = startRow+s;
				d[row] = val*weights[s];
				if (Double.isNaN(d[s]) || d[s]<0)
					throw new IllegalStateException("d["+s+"]="+d[s]+" is NaN or 0!  target="+target);
			}
			return numNonZeroElements;
		}
		
	}
	
	private static void printA(DoubleMatrix2D A, ConstraintRange range) {
		for (int row=range.startRow; row<range.endRow; row++) {
			for (int col=0; col<A.columns(); col++) {
				if (col > 0)
					System.out.print(" ");
				System.out.print((float)A.get(row, col));
			}
			System.out.println();
		}
	}
	
	private static void printD(double[] d, ConstraintRange range) {
		for (int row=range.startRow; row<range.endRow; row++)
			System.out.println((float)d[row]);
	}
	
	private static void misfitsDebugFor1s(DoubleMatrix2D A, double[] d, ConstraintRange range) {
		double[] sol1 = new double[A.columns()];
		for (int i=0; i<sol1.length; i++)
			sol1[i] = 1d;
		misfitsDebug(A, d, sol1, range);
	}
	
	private static void misfitsDebug(DoubleMatrix2D A, double[] d, double[] sol, ConstraintRange range) {
		DenseDoubleMatrix1D sol_clone = new DenseDoubleMatrix1D(sol);
		
		DenseDoubleMatrix1D syn = new DenseDoubleMatrix1D(A.rows());
		A.zMult(sol_clone, syn);
		
		System.out.println("\tData\tSnthetic\tMisfit");
		MinMaxAveTracker valTrack = new MinMaxAveTracker();
		MinMaxAveTracker misfitTrack = new MinMaxAveTracker();
		for (int row=range.startRow; row<range.endRow; row++) {
			double val = syn.get(row);
			double misfit = val - d[row];
			if (range.inequality)
				misfit = Math.max(0d, misfit);
			valTrack.addValue(val);
			misfitTrack.addValue(misfit);
			System.out.println("\t"+(float)d[row]+"\t"+(float)val+"\t"+(float)misfit);
		}
		System.out.println("Value stats:\t"+valTrack);
		System.out.println("Misfit stats:\t"+misfitTrack);
	}
	
	/**
	 * This constraint ensures that the co-rupture rate with other sections doesn't break the rate budget in each
	 * multifault rupture MFD bin
	 */
	private static class SectCoruptureBudgetConstraint extends InversionConstraint {

		private RupSetCoruptureMFDStructure structure;
		
		private List<SectParentCoruptureSet> corupSects;

		public SectCoruptureBudgetConstraint(RupSetCoruptureMFDStructure structure,
				double weight, ConstraintWeightingType weightingType) {
			// true means inequality constraint
			super("Section Corupture Budget Constraint", "CorupBudget", weight, true, weightingType);
			this.structure = structure;
			
			FaultSystemRupSet rupSet = structure.getRupSet();
			int numSections = rupSet.getNumSections();
			
			int[] sectMinMagIndexes = structure.getSectMinMagIndexes();
			int[] sectMaxMagIndexes = structure.getSectMaxMagIndexes();
			
			corupSects = new ArrayList<>();
			for (int s=0; s<numSections; s++) {
				for (int magIndex=sectMinMagIndexes[s]; magIndex<=sectMaxMagIndexes[s]; magIndex++) {
					SectMagRupturePathways pathways = structure.getMagBinPathways(s, magIndex);
					if (pathways == null)
						continue;
					
					for (int parentID : pathways.parentWeights().keySet()) {
						double parentWeight = pathways.parentWeights().get(parentID);
						int[] sectCounts = pathways.parentSectCounts().get(parentID);
						int sectCountSum = 0;
						for (int sectCount : sectCounts)
							sectCountSum += sectCount;
						Preconditions.checkState(sectCountSum > 0);
						List<FaultSection> parentSects = structure.getSectsForParent(parentID);
						Preconditions.checkState(parentSects.size() == sectCounts.length);
						double[] sectWeights = new double[sectCounts.length];
						for (int i=0; i<sectCounts.length; i++)
							sectWeights[i] = parentWeight * (double)sectCounts[i]/(double)sectCountSum;
						
						corupSects.add(new SectParentCoruptureSet(s, magIndex, parentID, parentSects, sectWeights));
					}
				}
			}
		}

		@Override
		public int getNumRows() {
			return corupSects.size();
		}

		@Override
		public long encode(DoubleMatrix2D A, double[] d, int startRow) {
			FaultSystemRupSet rupSet = structure.getRupSet();
			
			int numSections = rupSet.getNumSections();
			Preconditions.checkState(A.columns() == structure.getNumColumns());
			List<IncrementalMagFreqDist> origSectSupraSeisMFDs = structure.getOrigSectSupraSeisMFDs();
			double[] targetSectSupraSlipRates = structure.getTargetSectSupraSlipRates();
			double[] sectSupraSlipRateStdDevs = structure.getSectSupraSlipRateStdDevs();
			
			double[] weights = new double[numSections];
			for (int s=0; s<numSections; s++)
				weights[s] = this.weight;
			if (weightingType == ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY) {
				for (int s=0; s<numSections; s++)
					if (sectSupraSlipRateStdDevs[s] != 0d)
						weights[s] /= sectSupraSlipRateStdDevs[s];
			}
			
			long numNonZeroElements = 0;
			
			for (int i=0; i<corupSects.size(); i++) {
				int row = startRow + i;
				SectParentCoruptureSet corups = corupSects.get(i);
				int magIndex = corups.magIndex;
				
				double myWeight = this.weight;
				if (weightingType == ConstraintWeightingType.NORMALIZED) {
					// normalize by the original value for this sect/mag
					myWeight /= origSectSupraSeisMFDs.get(corups.fromSectID).getY(magIndex);
				} else if (weightingType == ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY) {
					if (sectSupraSlipRateStdDevs[corups.fromSectID] > 0d && targetSectSupraSlipRates[corups.fromSectID] > 0d) {
						// apply fractional slip rate uncertainty
						myWeight /= sectSupraSlipRateStdDevs[corups.fromSectID]/targetSectSupraSlipRates[corups.fromSectID];
					}
					// normalize by the original rate for this sect/mag
					myWeight /= origSectSupraSeisMFDs.get(corups.fromSectID).getY(magIndex);
				}
				
				// data vector is always 0. we ensure that the sum across the row is <= 0 when satisfied
				d[row] = 0;
				// the positive side is our section's magnitude (scaled by the total parent weight if <1)
				double sumParentWeight = StatUtils.sum(corups.toSectWeights);
				
				int fromCol = structure.locateSectMagColumn(corups.fromSectID, magIndex);
				double fromOrigRate = structure.getSectOrigRate(corups.fromSectID, magIndex);
				setA(A, row, fromCol, myWeight*sumParentWeight*fromOrigRate);
				numNonZeroElements++;
				// the negative side is the sum across all other sections for this parent (scaled by their weight)
				List<FaultSection> toParentSects = corups.toParentSects();
				for (int j=0; j<toParentSects.size(); j++) {
					if (corups.toSectWeights[j] > 0 ) {
						int toSectID = toParentSects.get(j).getSectionId();
						int toCol = structure.locateSectMagColumn(toSectID, magIndex);
						double toOrigRate = structure.getSectOrigRate(toSectID, magIndex);
						setA(A, row, toCol, -myWeight*corups.toSectWeights[j]*toOrigRate);
						numNonZeroElements++;
					}
				}
			}
			
			return numNonZeroElements;
		}
		
		private static record SectParentCoruptureSet(int fromSectID, int magIndex, int toParentID,
				List<FaultSection> toParentSects, double[] toSectWeights) {};
		
	}
	
	private static class DoubleAverager {
		private double sum = 0d;
		private int count = 0;
		
		public void add(double value) {
			sum += value;
			count++;
		}
		
		public double getAverage() {
			return sum/(double)count;
		}
	}
	
	/**
	 * Segmentation constraint (corup budget only uses seg for path weighting)
	 */
	public static class MFD_SegmentationConstraint extends InversionConstraint {
		
		// never let a weight exceed this value, happens if rupture probability or section rate estimate is exceedingly low 
		private static final double MAX_WEIGHT_SCALAR = 1e5;
		
		private final static boolean D = false;
		
		private transient RupSetCoruptureMFDStructure structure;
		private transient FaultSystemRupSet rupSet;
		private transient Map<UniqueDistJump, List<Integer>> jumpRupsMap;
		
		public MFD_SegmentationConstraint(RupSetCoruptureMFDStructure structure, double weight) {
			super("MFD Segmentation", "MFD-Seg", weight, true, ConstraintWeightingType.NORMALIZED);
			this.structure = structure;
			this.rupSet = structure.getRupSet();		}
		
		private synchronized void checkInitJumpRups() {
			if (jumpRupsMap == null) {
				jumpRupsMap = new HashMap<>();
				
				ClusterRuptures cRups = rupSet.requireModule(ClusterRuptures.class);
				
				for (int r=0; r<cRups.size(); r++) {
					ClusterRupture rup = cRups.get(r);
					for (Jump jump : rup.getJumpsIterable()) {
						UniqueDistJump udJump = new UniqueDistJump(jump);
						List<Integer> jumpRups = jumpRupsMap.get(udJump);
						if (jumpRups == null) {
							jumpRups = new ArrayList<>();
							jumpRupsMap.put(udJump, jumpRups);
						}
						jumpRups.add(r);
						// now add it reversed
						udJump = udJump.reverse();
						jumpRups = jumpRupsMap.get(udJump);
						if (jumpRups == null) {
							jumpRups = new ArrayList<>();
							jumpRupsMap.put(udJump, jumpRups);
						}
						jumpRups.add(r);
					}
				}
			}
		}

		@Override
		public int getNumRows() {
			checkInitJumpRups();
			return jumpRupsMap.size();
		}

		@Override
		public long encode(DoubleMatrix2D A, double[] d, int startRow) {
			long count = 0l;
			
			Preconditions.checkState(A.columns() == structure.getNumColumns());
			
			int row = startRow;
			
			checkInitJumpRups();
			List<UniqueDistJump> allJumps = new ArrayList<>(jumpRupsMap.keySet());
			allJumps.sort(Jump.id_comparator); // sort for consistent row ordering
			
			ClusterRuptures cRups = rupSet.requireModule(ClusterRuptures.class);
			JumpProbabilityCalc segModel = structure.getSegModel();
			
			EvenlyDiscretizedFunc refMFD = structure.getRefMFD();
			
			SectParticipationRateEstimator rateEst = new GRParticRateEstimator(
					rupSet, structure.origSectSupraSeisMFDs, structure.sectRupUtilizations);
			
			int[] sectMinMags = structure.getSectMinMagIndexes();
			int[] sectMaxMags = structure.getSectMaxMagIndexes();
			
			long rawAddCount = 0;
			
			for (UniqueDistJump jump : allJumps) {
				List<Integer> rupsUsingJump = jumpRupsMap.get(jump);
				Preconditions.checkNotNull(rupsUsingJump != null);
				Preconditions.checkState(!rupsUsingJump.isEmpty());
				
				MinMaxAveTracker probTrack = new MinMaxAveTracker();
				for (int r : rupsUsingJump) {
					ClusterRupture rup = cRups.get(r);
					RuptureTreeNavigator nav = rup.getTreeNavigator();
					Jump myJump = nav.getJump(jump.fromSection, jump.toSection);
					Preconditions.checkState(myJump.fromSection.getSectionId() == jump.fromSection.getSectionId());
					
					// see if either end needs to be reversed
					if (!myJump.fromCluster.endSects.contains(myJump.fromSection))
						myJump = new Jump(myJump.fromSection, myJump.fromCluster.reversed(),
								myJump.toSection, myJump.toCluster, myJump.distance);
					Preconditions.checkState(myJump.fromCluster.endSects.contains(myJump.fromSection));
					if (!myJump.toCluster.startSect.equals(myJump.toSection))
						myJump = new Jump(myJump.fromSection, myJump.fromCluster,
								myJump.toSection, myJump.toCluster.reversed(), myJump.distance);
					Preconditions.checkState(myJump.toCluster.startSect.equals(myJump.toSection));
					
					double prob = segModel.calcJumpProbability(rup, myJump, false);
//					if (probTrack.getNum() > 0)
//						Preconditions.checkState((float)prob == (float)probTrack.getAverage(),
//								"%s != %s for jump %s", prob, probTrack.getAverage(), jump);
					probTrack.addValue(prob);
				}
				
				double jumpCondProb = probTrack.getAverage();
				Preconditions.checkState(jumpCondProb >= 0 && jumpCondProb <= 1d);
				if (jumpCondProb == 0) {
					row++;
					continue;
				}
				
				double maxWeight = this.weight*MAX_WEIGHT_SCALAR;
				
				double rateEstWeight = this.weight;
				// scale weight by that estimated total event rate for this section
				double estRate = rateEst.estimateSectParticRate(jump.fromSection.getSectionId());
				
				if (estRate > 0d)
					rateEstWeight /= estRate;
				else
					rateEstWeight = maxWeight;
				
				if (D) System.out.println("Building constraint for jump: "+jump+" with "+rupsUsingJump.size()
					+" rups, prob="+(float)jumpCondProb+" with origWeight="+(float)weight
					+", rateEstWeight="+(float)rateEstWeight);
				
				double effectiveWeight = rateEstWeight/jumpCondProb;
				if (effectiveWeight > maxWeight) {
					if (D) System.err.println("WARNING: capping weight at max="+maxWeight+", would have been "+effectiveWeight);
					rateEstWeight = maxWeight*jumpCondProb;
				}
				
				// sum up the MFD rate on this section, but negative
				double scalarSectBins = -rateEstWeight;
				// then we'll add the MFD bins using this jump on the positive side
				double scalarAllUsing = rateEstWeight/jumpCondProb;
				
				// the inversion (inequality) will ensure that the negative side is net larger
				
				Preconditions.checkState(Double.isFinite(scalarAllUsing),
						"Bad scalarAllUsing=%s for jump %s with jumpCondProb=%s and weight=%s",
						scalarAllUsing, jump, jumpCondProb, rateEstWeight);
				Preconditions.checkState(Double.isFinite(scalarSectBins),
						"Bad scalarSectBins=%s for jump %s with jumpCondProb=%s and weight=%s",
						scalarSectBins, jump, jumpCondProb, rateEstWeight);
				
				// process all mfd bins for this section
				int sectIndex = jump.fromSection.getSubSectionIndex();
				// will need to scale from nuclation to participation
				double sectArea = rupSet.getAreaForSection(sectIndex);
				BitSet sectRups = structure.getSectRupUtilizations(sectIndex);
				int sectNumMags = 1 + sectMaxMags[sectIndex] - sectMinMags[sectIndex];
				double[] magBinnedRupAreaSums = new double[sectNumMags];
				int[] magBinnedRupCounts = new int[sectNumMags];
				for (int r = sectRups.nextSetBit(0); r >= 0; r = sectRups.nextSetBit(r+1)){
					int magIndex = refMFD.getClosestXIndex(rupSet.getMagForRup(r));
					int m = magIndex-sectMinMags[sectIndex];
					magBinnedRupAreaSums[m] += rupSet.getAreaForRup(r);
					magBinnedRupCounts[m]++;
				}
				for (int m=0; m<sectNumMags; m++) {
					if (magBinnedRupCounts[m] > 0) {
						double avgRupArea = magBinnedRupAreaSums[m] / (double)magBinnedRupCounts[m];
						// bin nuclRate = particRate * sectArea / avgRupArea
						// bin particRate = nuclRate * avgRupArea / sectArea
						double particScalar = avgRupArea / sectArea;
						
						int col = structure.locateSectMagColumn(sectIndex, sectMinMags[sectIndex]+m);
						double origRate = structure.getSectOrigRate(sectIndex, sectMinMags[sectIndex]+m);
						
						if (!addA(A, row, col, particScalar*scalarSectBins*origRate))
							count++;
						rawAddCount++;
					}
				}
				
				// now process all the ruptures using this jump
				// we're processing it for each section using the jump, so this will sum to participation rates
				for (int rupIndex : rupsUsingJump) {
					int magIndex = refMFD.getClosestXIndex(rupSet.getMagForRup(rupIndex));
					for (int oSectIndex : rupSet.getSectionsIndicesForRup(rupIndex)) {
						double weight = structure.getMultifaultRupFractContributionToBin(oSectIndex, magIndex, rupIndex);
						if (weight > 0d) {
							int col = structure.locateSectMagColumn(oSectIndex, magIndex);
							double origRate = structure.getSectOrigRate(oSectIndex, magIndex);
							if (!addA(A, row, col, scalarAllUsing*origRate))
								count++;
							rawAddCount++;
						}
					}
				}
				
//				if (D) {
//					System.out.println("Row for jump "+jump+" with P="+(float)jumpCondProb+":\n\t");
//					for (int col=0; col<A.columns(); col++)
//						System.out.print((float)getA(A, row, col)+"\t");
//					System.out.println();
//				}
				
				row++;
			}
			
			int rows = row-startRow;
			long maxPossibleCount = rows * structure.getNumColumns();
			Preconditions.checkState(count <= maxPossibleCount,
					"Count is impossibly-large; have %s, max possible is %s x %s = %s; rawAddCount=%s",
					count, rows, structure.getNumColumns(), maxPossibleCount, rawAddCount);
			
			return count;
		}
	}
	
	public static class ScaleFactorLimitConstraint extends InversionConstraint {

		private RupSetCoruptureMFDStructure structure;
		private boolean singleFault;

		public ScaleFactorLimitConstraint(RupSetCoruptureMFDStructure structure, boolean singleFault, double weight) {
			super(singleFault ? "Single Fault F>=1" : "Multi Fault F<=1", singleFault ? "SingleF>=1" : "MultiF<=1",
					weight, true, ConstraintWeightingType.NORMALIZED);
			this.structure = structure;
			this.singleFault = singleFault;
		}

		@Override
		public int getNumRows() {
			int numSections = structure.getRupSet().getNumSections();
			if (singleFault)
				return numSections;
			int rows = 0;
			for (int s=0; s<numSections; s++)
				rows += structure.getSectCommonPathways(s).size()-1;
			return rows;
		}

		@Override
		public long encode(DoubleMatrix2D A, double[] d, int startRow) {
			Preconditions.checkState(A.columns() == structure.getNumColumns());
			
			int row = startRow;
			long count = 0l;
			
			int numSections = structure.getRupSet().getNumSections();
			
			// if single fault, we're we're ensureing that the value does not go below 1
			double val;
			if (singleFault) {
				// ensure the inverted value does not go below 1
				val = -this.weight;
			} else {
				// ensure that the inverted value does not go above 1
				val = this.weight;
			}
			
			for (int sectIndex=0; sectIndex<numSections; sectIndex++) {
				List<SectCommonPathwaysMagRange> pathways = structure.getSectCommonPathways(sectIndex);
				int numProcessed = 0;
				for (int p=0; p<pathways.size(); p++) {
					SectCommonPathwaysMagRange pathway = pathways.get(p);
					boolean match;
					if (singleFault)
						match = pathway.parentPathways() == null || pathway.parentPathways().isEmpty();
					else
						match = pathway.parentPathways() != null && !pathway.parentPathways().isEmpty();
					
					if (match) {
						numProcessed++;
						
						int col = structure.getSectStartColumn(sectIndex) + p;
						setA(A, row, col, val);
						d[row] = val;
						count++;
						row++;
					}
				}
				if (singleFault)
					Preconditions.checkState(numProcessed == 1);
			}
			
			return count;
		}
	}
	
	public static class ScaleFactorOneConstraint extends InversionConstraint {

		private RupSetCoruptureMFDStructure structure;
		private double weightSingle;
		private double weightMulti;
		
		public ScaleFactorOneConstraint(RupSetCoruptureMFDStructure structure, double weight) {
			this(structure, weight, weight);
		}

		public ScaleFactorOneConstraint(RupSetCoruptureMFDStructure structure, double weightSingle, double weightMulti) {
			super("Scale Factor of 1", "F~1", 0.5*(weightSingle + weightMulti), true, ConstraintWeightingType.NORMALIZED);
			this.structure = structure;
			this.weightSingle = weightSingle;
			this.weightMulti = weightMulti;
		}

		@Override
		public int getNumRows() {
			return structure.getNumColumns();
		}

		@Override
		public long encode(DoubleMatrix2D A, double[] d, int startRow) {
			Preconditions.checkState(A.columns() == structure.getNumColumns());
			
			int row = startRow;
			long count = 0l;
			
			int numSections = structure.getRupSet().getNumSections();
			
			for (int sectIndex=0; sectIndex<numSections; sectIndex++) {
				List<SectCommonPathwaysMagRange> pathways = structure.getSectCommonPathways(sectIndex);
				for (int p=0; p<pathways.size(); p++) {
					double weight = p == 0 ? weightSingle : weightMulti;
					int col = structure.getSectStartColumn(sectIndex) + p;
					setA(A, row, col, weight);
					d[row] = weight;
					count++;
					row++;
				}
			}
			
			return count;
		}
	}
	
	/**
	 * Constraint to try to minimize total rate for each section, to force it to use as much of the large magnitudes as
	 * it can, subject to the inequality and slip rate balancing constraints
	 */
	public static class SectRateMinimizationConstraint extends InversionConstraint {
		
		private transient RupSetCoruptureMFDStructure structure;
		
		public SectRateMinimizationConstraint(RupSetCoruptureMFDStructure structure, double weight, ConstraintWeightingType weightingType) {
			super("Sect Rate Minimization",
					"SectMinimum", weight, false, weightingType);
			this.structure = structure;
		}

		@Override
		public int getNumRows() {
			return structure.getRupSet().getNumSections();
		}

		@Override
		public long encode(DoubleMatrix2D A, double[] d, int startRow) {
			long count = 0l;
			
			Preconditions.checkState(A.columns() == structure.getNumColumns());
			
			int row = startRow;

			int numSections = structure.getRupSet().getNumSections();
			int[] sectMinMagIndexes = structure.getSectMinMagIndexes();
			int[] sectMaxMagIndexes = structure.getSectMaxMagIndexes();
			
			List<IncrementalMagFreqDist> origMFDs = structure.getOrigSectSupraSeisMFDs();
			double[] weights = new double[numSections];
			for (int s=0; s<numSections; s++)
				weights[s] = this.weight;
			if (weightingType == ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY) {
				double[] sectSlipRates = structure.getTargetSectSupraSlipRates();
				double[] sectSlipRateStdDevs = structure.getSectSupraSlipRateStdDevs();
				for (int s=0; s<numSections; s++) {
					double origTotRate = origMFDs.get(s).calcSumOfY_Vals();
					if (sectSlipRates[s] != 0d && sectSlipRateStdDevs[s] != 0d) {
						double fractUncert = sectSlipRateStdDevs[s] / sectSlipRates[s];
						weights[s] /= origTotRate*fractUncert;
					}
				}
			} else if (weightingType == ConstraintWeightingType.NORMALIZED) {
				for (int s=0; s<numSections; s++) {
					double origTotRate = origMFDs.get(s).calcSumOfY_Vals();
					weights[s] /= origTotRate;
				}
			}

			for (int sectIndex=0; sectIndex<numSections; sectIndex++) {
				double weight = weights[sectIndex];
				for (int m=sectMinMagIndexes[sectIndex]; m<=sectMaxMagIndexes[sectIndex]; m++) {
					int col = structure.locateSectMagColumn(sectIndex, m);
					double origRate = structure.getSectOrigRate(sectIndex, m);
					if (!addA(A, row, col, weight*origRate))
						count++;
				}
				d[row] = 0d;
				row++;
			}
			
			return count;
		}
	}

}
