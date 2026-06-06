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
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionInputGenerator;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.ConstraintWeightingType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.InversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.UncertainDataConstraint;
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
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.JumpProbabilityCalc;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RuptureTreeNavigator;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.estimators.SectNucleationMFD_Estimator;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;

public class MFD_InversionAdjustment extends SectNucleationMFD_Estimator {
	
	private JumpProbabilityCalc segModel;
	
	private RupSetCoruptureStructure structure;
	private double[] solution;
	private int threads;
	
	public MFD_InversionAdjustment(int threads, JumpProbabilityCalc segModel) {
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
		structure = new RupSetCoruptureStructure(rupSet, origSectSupraSeisMFDs,
				targetSectSupraMoRates, targetSectSupraSlipRates, sectSupraSlipRateStdDevs,
				sectRupUtilizations, sectMinMagIndexes, sectMaxMagIndexes, sectRupInBinCounts, refMFD, segModel);
		
		if (verbose) System.out.println("Building constraints");
		List<InversionConstraint> constraints = new ArrayList<>();
		
		constraints.add(new MFD_SectSlipRateConstraint(structure, 1d, ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY));
		constraints.add(new SectCoruptureBudgetConstraint(structure, 1d, ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY));
		constraints.add(new RelativeGRInequalityToSingleFaultConstraint(structure, 1e3));
		constraints.add(new RelativeCommonPathGREqualityConstraint(structure, 1e2));
		constraints.add(new SectRateMinimizationConstraint(structure, 1d, ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY));
		if (segModel != null) {
			// TODO: seg constraint
		}
		
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
		
		DoubleMatrix2D A = new SparseDoubleMatrix2D(numRows, columns);
		double[] d = new double[numRows];
		
		DoubleMatrix2D A_ineq = new SparseDoubleMatrix2D(numIneqRows, columns);
		double[] d_ineq = new double[numIneqRows];
		
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
		
		// run the inversion
		ColumnOrganizedAnnealingData equalityData = new ColumnOrganizedAnnealingData(A, d);
		ColumnOrganizedAnnealingData inequalityData = new ColumnOrganizedAnnealingData(A_ineq, d_ineq);
		
		double[] initial = new double[columns];
		long itersPerSectBin = 10000l;
		long totalIters = itersPerSectBin * columns;
		CompletionCriteria completion = new IterationCompletionCriteria(itersPerSectBin);
		
		GenerationFunctionType perturb = GenerationFunctionType.EXPONENTIAL_SCALE;
		NonnegativityConstraintType nonneg = NonnegativityConstraintType.TRY_ZERO_RATES_OFTEN;
		CoolingScheduleType cool = CoolingScheduleType.FAST_SA;
		
		SimulatedAnnealing sa;
		if (threads > 1) {
			
			int avgThreads = Integer.max(threads/4, 2);
			
			int threadsPerAvg = (int)Math.ceil((double)threads/(double)avgThreads);
			Preconditions.checkState(threadsPerAvg <= threads);
			Preconditions.checkState(threadsPerAvg > 0);
			
			CompletionCriteria avgCompletion = new IterationCompletionCriteria(itersPerSectBin/10);
			CompletionCriteria subCompletion = new IterationCompletionCriteria(itersPerSectBin/100);
			
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
		sa.setRuptureSampler(structure.getColSampler());
		
		sa.iterate(completion);
		
		solution = sa.getBestSolution();
		
		if (verbose)
			System.out.println("DONE estimating MFDs");
	}

	@Override
	public IncrementalMagFreqDist estimateNuclMFD(FaultSection sect, IncrementalMagFreqDist curSectSupraSeisMFD,
			List<Integer> availableRupIndexes, List<Double> availableRupMags, UncertainDataConstraint sectMomentRate,
			boolean sparseGR) {
		EvenlyDiscretizedFunc refMFD = structure.getRefMFD();
		int sectIndex = sect.getSectionId();
		int minMag = structure.getSectMinMagIndexes()[sectIndex];
		int maxMag = structure.getSectMaxMagIndexes()[sectIndex];
		IncrementalMagFreqDist mfd = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		for (int magIndex=minMag; magIndex<=maxMag; magIndex++) {
			int col = structure.getColumn(sectIndex, magIndex);
			double val = solution[col];
			Preconditions.checkState(val >= 0d);
			mfd.set(magIndex, val);
		}
		return mfd;
	}
	
	private static class RupSetCoruptureStructure {
		
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
		
		private IntegerSampler colSampler;
		
		public RupSetCoruptureStructure(FaultSystemRupSet rupSet, List<IncrementalMagFreqDist> origSectSupraSeisMFDs,
				double[] targetSectSupraMoRates, double[] targetSectSupraSlipRates, double[] sectSupraSlipRateStdDevs,
				List<BitSet> sectRupUtilizations, int[] sectMinMagIndexes, int[] sectMaxMagIndexes,
				int[][] sectRupInBinCounts, EvenlyDiscretizedFunc refMFD, JumpProbabilityCalc segModel) {
			this.rupSet = rupSet;
			this.origSectSupraSeisMFDs = origSectSupraSeisMFDs;
			this.targetSectSupraMoRates = targetSectSupraMoRates;
			this.targetSectSupraSlipRates = targetSectSupraSlipRates;
			this.sectSupraSlipRateStdDevs = sectSupraSlipRateStdDevs;
			this.sectRupUtilizations = sectRupUtilizations;
			this.sectRupInBinCounts = sectRupInBinCounts;
			this.refMFD = refMFD;
			this.segModel = segModel;
			
			int numRuptures = rupSet.getNumRuptures();
			int numSections = rupSet.getNumSections();
			int numMags = refMFD.size();
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
			
			// determine column indexes of the A matrix, skipping MFD bins before/above sect Mmin/Mmax
			sectColIndexes = new int[numSections];
			this.sectMinMagIndexes = sectMinMagIndexes;
			this.sectMaxMagIndexes = sectMaxMagIndexes;
			numCols = 0;
			HashSet<Integer> excludedColIndexes = new HashSet<>();
			for (int s=0; s<numSections; s++) {
				sectColIndexes[s] = numCols;
				Preconditions.checkState(sectMinMagIndexes[s] >= 0);
				Preconditions.checkState(sectMaxMagIndexes[s] >= sectMinMagIndexes[s]);
				for (int m=sectMinMagIndexes[s]; m<=sectMaxMagIndexes[s]; m++) {
					if (sectRupInBinCounts[s][m] == 0)
						// empty bin, never sample it in the inversion
						excludedColIndexes.add(numCols);
					numCols++;
				}
				
			}
			colSampler = new IntegerSampler.ExclusionIntegerSampler(0, numCols, excludedColIndexes);
			
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
				List<SectCommonPathwaysMagRange> myCommonPathways = new ArrayList<>(myNumMags);
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
					Map<Set<Integer>, List<Double>> uniqueParentCombs = new HashMap<>();
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
						List<Double> pathWeights = uniqueParentCombs.get(parents);
						if (pathWeights == null) {
							pathWeights = new ArrayList<>();
							uniqueParentCombs.put(parents, pathWeights);
						}
						pathWeights.add(rupWeight);
					}
					if (anySingleFault) {
						sectMaxSingleFaultMags[s] = magIndex;
						curCommonPathwayMmax = magIndex;
						myPathways.add(null);
						// we have single-fault ruptures at this mag, skip the constraint
						continue;
					}
					
					// determine aggregate usage weights for each parent
					Map<Integer, Double> parentWeights = new HashMap<>(parentUsedSectCounts.size()-1);
					double sumPathWeights = 0d;
					List<Set<Integer>> parentPathways = new ArrayList<>(uniqueParentCombs.size());
					List<Double> pathwayWeights = new ArrayList<>(uniqueParentCombs.size());
					for (Set<Integer> parents : uniqueParentCombs.keySet()) {
						parentPathways.add(parents);
						double pathWeight = uniqueParentCombs.get(parents).stream().mapToDouble(D->D).average().getAsDouble();
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
					
					myPathways.add(new SectMagRupturePathways(s, sectMinMagIndexes[s]+m, parentPathways, pathwayWeights,
							parentWeights, parentUsedSectCounts));
					
					// now see if this magnitude bin used the same set of pathways
					if (curCommonPathways == null || !curCommonPathways.equals(uniqueParentCombs.keySet())) {
						// new set of pathways
						
						// store the previous one
						Preconditions.checkState(curCommonPathwayMmax >= curCommonPathwayMmin);
						Preconditions.checkState(curCommonPathwayMmax >= 0);
						myCommonPathways.add(new SectCommonPathwaysMagRange(s, curCommonPathwayMmin, curCommonPathwayMmax, curCommonPathways));
						
						curCommonPathways = uniqueParentCombs.keySet();
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
		}
		
		public int getColumn(int sectIndex, int magIndex) {
			Preconditions.checkState(magIndex >= sectMinMagIndexes[sectIndex]);
			Preconditions.checkState(magIndex <= sectMaxMagIndexes[sectIndex]);
			return sectColIndexes[sectIndex] + magIndex - sectMinMagIndexes[sectIndex];
		}
		
		public int getNumColumns() {
			return numCols;
		}
		
		public SectMagRupturePathways getPathways(int sectIndex, int magIndex) {
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
		
		public double getSectCoruptureFraction(int sectIndex, int magIndex, int toSectIndex) {
			SectMagRupturePathways sectPaths = getPathways(sectIndex, magIndex);
			if (sectPaths == null)
				return 0d;
			FaultSection toSect = rupSet.getFaultSectionData(toSectIndex);
			int toParentID = toSect.getParentSectionId();
			Double parentWeight = sectPaths.parentWeights().get(toParentID);
			if (parentWeight == null || parentWeight == 0d)
				return 0d;
			int[] parentSectCounts = sectPaths.parentSectCounts().get(toParentID);
			int toCount = parentSectCounts[toSect.getSubSectionIndex()];
			if (toCount == 0)
				return 0d;
			int sum = 0;
			for (int count : parentSectCounts)
				sum += count;
			return (double)toCount/(double)sum;
		}
		
		public double getParentSectCoruptureFraction(int sectIndex, int magIndex, int toParentID) {
			SectMagRupturePathways sectPaths = getPathways(sectIndex, magIndex);
			if (sectPaths == null)
				return 0d;
			Double parentWeight = sectPaths.parentWeights().get(toParentID);
			if (parentWeight == null)
				return 0d;
			return parentWeight;
		}
		
		public List<SectCommonPathwaysMagRange> getSectCommonPathways(int sectIndex) {
			return sectCommonPathways.get(sectIndex);
		}
		
		public int getSectMaxSingleFaultMagIndex(int sectIndex) {
			return sectMaxSingleFaultMags[sectIndex];
		}
		
		public IntegerSampler getColSampler() {
			return colSampler;
		}
	}
	
	private static record SectMagRupturePathways(int sectIndex, int magIndex, List<Set<Integer>> parentPathways,
			List<Double> pathwayWeights, Map<Integer, Double> parentWeights, Map<Integer, int[]> parentSectCounts) {};
	
	private static record SectCommonPathwaysMagRange(int sectIndex, int minMagIndex, int maxMagIndex, Set<Set<Integer>> parentPathways) {};
	
	/**
	 * This constraint ensures that each section's individual MFD sums up to match the target slip rate
	 */
	private static class MFD_SectSlipRateConstraint extends InversionConstraint {
		private RupSetCoruptureStructure structure;

		public MFD_SectSlipRateConstraint(RupSetCoruptureStructure structure,
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
			double[] targetSectSupraSlipRates = structure.getTargetSectSupraSlipRates();
			double[] sectSupraSlipRateStdDevs = structure.getSectSupraSlipRateStdDevs();
			
			int numSections = rupSet.getNumSections();
			int numCols = structure.getNumColumns();
			Preconditions.checkState(A.columns() == numCols);
			
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
			
			// figure out section slips
			int[] sectMagRupCounts = new int[numCols];
			double[] sectMagSumSlips = new double[numCols];
			for (int rupIndex : structure.getIncludedRupIndexes()) {
				double[] slips = slipAlongModule.calcSlipOnSectionsForRup(rupSet, aveSlipModule, rupIndex);
				List<Integer> sects = rupSet.getSectionsIndicesForRup(rupIndex);
				int magIndex = refMFD.getClosestXIndex(rupSet.getMagForRup(rupIndex));
				for (int s=0; s<slips.length; s++) {
					int sectID = sects.get(s);
					int colIndex = structure.getColumn(sectID, magIndex);
					sectMagRupCounts[colIndex]++;
					sectMagSumSlips[colIndex] += slips[s];
				}
			}
			// normalize
			for (int c=0; c<numCols; c++)
				sectMagSumSlips[c] /= (double)sectMagRupCounts[c];
		
			long numNonZeroElements = 0l;
			
			int[] sectMinMagIndexes = structure.getSectMinMagIndexes();
			int[] sectMaxMagIndexes = structure.getSectMaxMagIndexes();
			
			// A matrix component of slip-rate constraint 
			for (int s=0; s<numSections; s++) {
				int row = startRow + s;
				for (int magIndex=sectMinMagIndexes[s]; magIndex<=sectMaxMagIndexes[s]; magIndex++) {
					int col = structure.getColumn(s, magIndex);
					if (sectMagRupCounts[col] > 0) {
						double avgSlip = sectMagSumSlips[col] / (double)sectMagRupCounts[col];
						double val = avgSlip;
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
						setA(A, row, col, val*weights[s]);
						numNonZeroElements++;
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
	
	/**
	 * This constraint ensures that the co-rupture rate with other sections doesn't break the rate budget in each
	 * multifault rupture MFD bin
	 */
	private static class SectCoruptureBudgetConstraint extends InversionConstraint {

		private RupSetCoruptureStructure structure;
		
		private List<SectParentCoruptureSet> corupSects;

		public SectCoruptureBudgetConstraint(RupSetCoruptureStructure structure,
				double weight, ConstraintWeightingType weightingType) {
			// true means inequality constraint
			super("Section Corupture Budget Constraint", "CorupBudget", weight, true, weightingType);
			this.structure = structure;
			
			FaultSystemRupSet rupSet = structure.getRupSet();
			int numSections = rupSet.getNumSections();
			
			int[] sectMinMagIndexes = structure.getSectMinMagIndexes();
			int[] sectMaxMagIndexes = structure.getSectMaxMagIndexes();
			
//			int numSections = str
			corupSects = new ArrayList<>();
			for (int s=0; s<numSections; s++) {
				for (int magIndex=sectMinMagIndexes[s]; magIndex<=sectMaxMagIndexes[s]; magIndex++) {
					SectMagRupturePathways pathways = structure.getPathways(s, magIndex);
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
				
				double myWeight = this.weight;
				if (weightingType == ConstraintWeightingType.NORMALIZED) {
					// normalize by the original value for this sect/mag
					myWeight /= origSectSupraSeisMFDs.get(corups.fromSectID).getY(corups.magIndex);
				} else if (weightingType == ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY) {
					if (sectSupraSlipRateStdDevs[corups.fromSectID] > 0d && targetSectSupraSlipRates[corups.fromSectID] > 0d) {
						// apply fractional slip rate uncertainty
						myWeight /= sectSupraSlipRateStdDevs[corups.fromSectID]/targetSectSupraSlipRates[corups.fromSectID];
					}
					// normalize by the original rate for this sect/mag
					myWeight /= origSectSupraSeisMFDs.get(corups.fromSectID).getY(corups.magIndex);
				}
				
				// data vector is always 0. we ensure that the sum across the row is <= 0 when satisfied
				d[row] = 0;
				// the positive side is our section's magnitude (scaled by the total parent weight if <1)
				double sumParentWeight = StatUtils.sum(corups.toSectWeights);
				int col = structure.getColumn(corups.fromSectID, corups.magIndex);
				setA(A, row, col, myWeight*sumParentWeight);
				numNonZeroElements++;
				// the negative side is the sum across all other sections for this parent (scaled by their weight)
				List<FaultSection> toParentSects = corups.toParentSects();
				for (int j=0; j<toParentSects.size(); j++) {
					if (corups.toSectWeights[j] > 0 ) {
						setA(A, row, structure.getColumn(toParentSects.get(j).getSectionId(), corups.magIndex),
								-myWeight*corups.toSectWeights[j]);
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
	
	// TODO: segmentation constraint (corup budget only uses seg for path weighting)
	
	/**
	 * Inequality GR constraint for each section's multifault rupture bins, relative to the single-fault rupture bins
	 */
	public class RelativeGRInequalityToSingleFaultConstraint extends InversionConstraint {
		
		private transient RupSetCoruptureStructure structure;
		
		public RelativeGRInequalityToSingleFaultConstraint(RupSetCoruptureStructure structure, double weight) {
			super("GR Inequality Constraint (rel single fault)",
					"GR-Inequality",
					weight, true, ConstraintWeightingType.NORMALIZED);
			this.structure = structure;
		}

		@Override
		public int getNumRows() {
			int rows = 0;
			int numSections = structure.getRupSet().getNumSections();
			for (int sectIndex=0; sectIndex<numSections; sectIndex++) {
				int minSingleFaultMag = structure.getSectMinMagIndexes()[sectIndex];
				int maxSingleFaultMag = structure.getSectMaxSingleFaultMagIndex(sectIndex);
				int numSingleFaultMags = 1 + maxSingleFaultMag - minSingleFaultMag;
				for (SectCommonPathwaysMagRange commonPathways : structure.getSectCommonPathways(sectIndex)) {
					if (commonPathways.parentPathways() == null)
						// this is only for multi-fault rupture mag bins
						continue;
					int mMin = commonPathways.minMagIndex;
					int mMax = commonPathways.maxMagIndex;
					int numPathwaysMags = 1 + mMax - mMin;
					rows += numSingleFaultMags * numPathwaysMags;
				}
			}
			return rows;
		}

		@Override
		public long encode(DoubleMatrix2D A, double[] d, int startRow) {
			long count = 0l;
			
			Preconditions.checkState(A.columns() == structure.getNumColumns());
			
			int row = startRow;

			int numSections = structure.getRupSet().getNumSections();
			for (int sectIndex=0; sectIndex<numSections; sectIndex++) {
				IncrementalMagFreqDist origSectMFD = structure.getOrigSectSupraSeisMFDs().get(sectIndex);
				
				int minSingleFaultMag = structure.getSectMinMagIndexes()[sectIndex];
				int maxSingleFaultMag = structure.getSectMaxSingleFaultMagIndex(sectIndex);
				
				for (SectCommonPathwaysMagRange commonPathways : structure.getSectCommonPathways(sectIndex)) {
					if (commonPathways.parentPathways() == null)
						// this is only for multi-fault rupture mag bins
						continue;
					int mMin = commonPathways.minMagIndex;
					int mMax = commonPathways.maxMagIndex;
					
					for (int singleM=minSingleFaultMag; singleM<=maxSingleFaultMag; singleM++) {
						double val1 = origSectMFD.getY(singleM);
						int col1 = structure.getColumn(sectIndex, singleM);
						
						for (int multiM=mMin; multiM<=mMax; multiM++) {
							double val2 = origSectMFD.getY(multiM);
							int col2 = structure.getColumn(sectIndex, multiM);
							
							if (val1 > 0 && val2 > 0) {
								// ruptures in this row should have a final rate equal to this times the rate of
								// ruptures in the first (comparison) bin (regardless of a-value)
								double ratio = val1/val2;
								
								double aScalar = 1d/val1;
								
								double myWeight = weight*aScalar;
								
								// set the larger magnitude's rate to be negative
								setA(A, row, col1, -myWeight);
								// multiply the smaller magnitude's rate by the expected ratio, and force it to be positive
								setA(A, row, col2, myWeight*ratio);
								count += 2l;
							}
							d[row] = 0d;
							row++;
						}
					}
				}
			}
			
			return count;
		}
	}
	
	/**
	 * Equality GR constraint for each section contiguous magnitude range that uses the same pathways (or is single-fault)
	 */
	public class RelativeCommonPathGREqualityConstraint extends InversionConstraint {
		
		private transient RupSetCoruptureStructure structure;
		
		public RelativeCommonPathGREqualityConstraint(RupSetCoruptureStructure structure, double weight) {
			super("GR Equality Constraint (common pathways)",
					"GR-Equality",
					weight, false, ConstraintWeightingType.NORMALIZED);
			this.structure = structure;
		}

		@Override
		public int getNumRows() {
			int rows = 0;
			int numSections = structure.getRupSet().getNumSections();
			for (int sectIndex=0; sectIndex<numSections; sectIndex++) {
				for (SectCommonPathwaysMagRange commonPathways : structure.getSectCommonPathways(sectIndex)) {
					int mMin = commonPathways.minMagIndex;
					int mMax = commonPathways.maxMagIndex;
					for (int m1=mMin; m1<mMax; m1++)
						for (int m2=m1+1; m2<=mMax; m2++)
							rows++;
				}
			}
			return rows;
		}

		@Override
		public long encode(DoubleMatrix2D A, double[] d, int startRow) {
			long count = 0l;
			
			Preconditions.checkState(A.columns() == structure.getNumColumns());
			
			int row = startRow;

			int numSections = structure.getRupSet().getNumSections();
			for (int sectIndex=0; sectIndex<numSections; sectIndex++) {
				IncrementalMagFreqDist origSectMFD = structure.getOrigSectSupraSeisMFDs().get(sectIndex);
				for (SectCommonPathwaysMagRange commonPathways : structure.getSectCommonPathways(sectIndex)) {
					int mMin = commonPathways.minMagIndex;
					int mMax = commonPathways.maxMagIndex;
					for (int m1=mMin; m1<mMax; m1++) {
						double val1 = origSectMFD.getY(m1);
						int col1 = structure.getColumn(sectIndex, m1);
						for (int m2=m1+1; m2<=mMax; m2++) {
							double val2 = origSectMFD.getY(m2);
							int col2 = structure.getColumn(sectIndex, m2);
							
							if (val1 > 0 && val2 > 0) {
								// ruptures in this row should have a final rate equal to this times the rate of
								// ruptures in the first (comparison) bin (regardless of a-value)
								double ratio = val1/val2;
								
								double aScalar = 1d/val1;
								
								double myWeight = weight*aScalar;
								
								// set the larger magnitude's rate to be negative
								setA(A, row, col1, -myWeight);
								// multiply the smaller magnitude's rate by the expected ratio, and force it to be positive
								setA(A, row, col2, myWeight*ratio);
								count += 2l;
								
								d[row] = 0d;
							}
							
							row++;
						}
					}
				}
			}
			
			return count;
		}
	}
	
	/**
	 * Constraint to try to minimize total rate for each section, to force it to use as much of the large magnitudes as
	 * it can, subject to the inequality and slip rate balancing constraints
	 */
	public class SectRateMinimizationConstraint extends InversionConstraint {
		
		private transient RupSetCoruptureStructure structure;
		
		public SectRateMinimizationConstraint(RupSetCoruptureStructure structure, double weight, ConstraintWeightingType weightingType) {
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
					
					int col = structure.getColumn(sectIndex, m);
					setA(A, row, col, weight);
					count++;
				}
				d[row] = 0d;
				row++;
			}
			
			return count;
		}
	}

}
