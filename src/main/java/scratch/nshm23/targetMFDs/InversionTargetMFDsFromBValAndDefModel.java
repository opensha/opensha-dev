package scratch.nshm23.targetMFDs;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.uncertainty.BoundedUncertainty;
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertainIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.Uncertainty;
import org.opensha.commons.data.uncertainty.UncertaintyBoundType;
import org.opensha.commons.geo.Region;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ClassUtils;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.modules.ArchivableModule;
import org.opensha.commons.util.modules.SubModule;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.MFDInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.UncertainDataConstraint;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.ModSectMinMags;
import org.opensha.sha.earthquake.faultSysSolution.modules.SectSlipRates;
import org.opensha.sha.earthquake.faultSysSolution.modules.SubSeismoOnFaultMFDs;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.SolMFDPlot;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;

import scratch.UCERF3.enumTreeBranches.MaxMagOffFault;
import scratch.UCERF3.enumTreeBranches.TotalMag5Rate;
import scratch.UCERF3.inversion.UCERF3InversionInputGenerator;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;

public class InversionTargetMFDsFromBValAndDefModel extends InversionTargetMFDs implements ArchivableModule {
	
	// inputs
	private double supraSeisBValue;
	
	// things we compute here
	private SectSlipRates targetSlipRates;
	private List<IncrementalMagFreqDist> sectFullMFDs;
	private List<UncertainIncrMagFreqDist> sectSupraSeisMFDs;
	private SubSeismoOnFaultMFDs sectSubSeisMFDs;
	private UncertainIncrMagFreqDist totalOnFaultSupra;
	private IncrementalMagFreqDist totalOnFaultSub;
	private List<IncrementalMagFreqDist> mfdConstraints;
	
	private IncrementalMagFreqDist totalRegional;
	
	private double[] sectFractSupras;
	
	// discretization parameters for MFDs
	public final static double MIN_MAG = 0.05;
	public final static double MAX_MAG = 8.95;
	public final static int NUM_MAG = 90;
	public final static double DELTA_MAG = 0.1;

	private static double TARGET_REL_STD_DEV_DEFAULT = 0.1;
	private static boolean WEIGHT_STD_DEV_BY_PARTICIPATION_DEFAULT = false;
	private static boolean APPLY_DEF_MODEL_UNCERTAINTIES_DEFAULT = false;
	private static boolean USE_EXISTING_TARGET_SLIP_RATES_DEFAULT = false;
	
	private static boolean SPARSE_GR_DEFAULT = true;

	public InversionTargetMFDsFromBValAndDefModel(FaultSystemRupSet rupSet, double supraSeisBValue) {
		this(rupSet, supraSeisBValue, SPARSE_GR_DEFAULT, TARGET_REL_STD_DEV_DEFAULT,
				APPLY_DEF_MODEL_UNCERTAINTIES_DEFAULT, WEIGHT_STD_DEV_BY_PARTICIPATION_DEFAULT, null);
	}
	
	public InversionTargetMFDsFromBValAndDefModel(FaultSystemRupSet rupSet, double supraSeisBValue, boolean sparseGR,
			double defaultRelStdDev, boolean applyDefModelUncertainties, boolean addSectCountUncertainties,
			List<Region> constrainedRegions) {
		this(rupSet, supraSeisBValue, sparseGR, defaultRelStdDev, applyDefModelUncertainties,
				addSectCountUncertainties, USE_EXISTING_TARGET_SLIP_RATES_DEFAULT, null, null, constrainedRegions);
	}
	
	public InversionTargetMFDsFromBValAndDefModel(FaultSystemRupSet rupSet, double supraSeisBValue, boolean sparseGR,
			double defaultRelStdDev, boolean applyDefModelUncertainties, boolean addSectCountUncertainties,
			boolean useExistingTargetSlipRates, List<DataSectNucleationRateEstimator> dataConstraints,
			UncertaintyBoundType expandUncertToDataBound, List<Region> constrainedRegions) {
		super(rupSet);
		this.supraSeisBValue = supraSeisBValue;
		
		ModSectMinMags minMags = rupSet.getModule(ModSectMinMags.class);
		
		int numSects = rupSet.getNumSections();
		
		SectSlipRates inputSlipRates = null;
		if (useExistingTargetSlipRates)
			inputSlipRates = rupSet.requireModule(SectSlipRates.class);
		double[] slipRates = new double[numSects];
		double[] slipRateStdDevs = new double[numSects];
		sectFullMFDs = new ArrayList<>();
		List<IncrementalMagFreqDist> sectSubSeisMFDs = new ArrayList<>();
		
		// supra-seismogenic MFDs for the given b-value
		List<IncrementalMagFreqDist> sectSupraSeisMFDs = new ArrayList<>();
		// standard deviations of those MFDs inferred from deformation model uncertainties
		List<EvenlyDiscretizedFunc> defModMFDStdDevs = applyDefModelUncertainties ? new ArrayList<>() : null;
		// supra-seismogenic MFDs implied by supplied data constraints
		List<Future<ImpliedDataCallable>> dataCalcFutures = null;
		ExecutorService exec = null;
		if (dataConstraints != null && !dataConstraints.isEmpty() && expandUncertToDataBound != null) {
			dataCalcFutures = new ArrayList<>();
			exec = Executors.newFixedThreadPool(FaultSysTools.defaultNumThreads());
		}
		
		SummedMagFreqDist totalOnFaultSub = new SummedMagFreqDist(MIN_MAG, NUM_MAG, DELTA_MAG);
		
		MinMaxAveTracker fractSuprasTrack = new MinMaxAveTracker();
		
		sectFractSupras = new double[numSects];

		for (int s=0; s<numSects; s++) {
			FaultSection sect = rupSet.getFaultSectionData(s);
			double creepReducedSlipRate = sect.getReducedAveSlipRate()*1e-3; // mm/yr -> m/yr
			double creepReducedSlipRateStdDev = sect.getReducedSlipRateStdDev()*1e-3; // mm/yr -> m/yr
			
			double area = rupSet.getAreaForSection(s); // m
			
			// convert it to a moment rate
			// TODO: is this consistent with scaling relationships used?
			double targetMoRate = FaultMomentCalc.getMoment(area, creepReducedSlipRate);
			
			// supra-seismogenic minimum magnitude
			double sectMinMag = rupSet.getMinMagForSection(s);
			
			if (minMags != null) {
				// not below the section minimum magnitude
				sectMinMag = Math.max(sectMinMag, minMags.getMinMagForSection(s));
			}
			
			List<Double> mags = new ArrayList<>();
			List<Integer> rups = new ArrayList<>();
			// make sure we actually have a rupture at that magnitude, otherwise there can be empty bins without
			// any rupture at/above the section minimum magnitude but below the first rupture's bin
			double minAbove = Double.POSITIVE_INFINITY;
			for (int r : rupSet.getRupturesForSection(s)) {
				double mag = rupSet.getMagForRup(r);
				if ((float)mag >= (float)sectMinMag) {
					minAbove = Math.min(mag, minAbove);
					mags.add(mag);
					rups.add(r);
				}
			}
			sectMinMag = minAbove;
			double sectMaxMag = rupSet.getMaxMagForSection(s);
			
			// construct a full G-R including sub-seismogenic ruptures
			int minMagIndex = totalOnFaultSub.getClosestXIndex(sectMinMag);
			int maxMagIndex = totalOnFaultSub.getClosestXIndex(sectMaxMag);
			int targetMFDNum = maxMagIndex+1;
			
			IncrementalMagFreqDist sectFullMFD, subSeisMFD, supraSeisMFD;
			double supraMoRate, subMoRate, fractSupra;
			if (useExistingTargetSlipRates) {
				double supraSlipRate = inputSlipRates.getSlipRate(s); // m/yr
				double supraSlipStdDev = inputSlipRates.getSlipRateStdDev(s); // m/yr
				fractSupra = supraSlipRate/creepReducedSlipRate;
				Preconditions.checkState(fractSupra > 0d && fractSupra <= 1d);
				
				supraMoRate = targetMoRate*fractSupra;
				subMoRate = targetMoRate-supraMoRate;
				
				GutenbergRichterMagFreqDist supraGR = new GutenbergRichterMagFreqDist(
						totalOnFaultSub.getX(minMagIndex), 1+maxMagIndex-minMagIndex, DELTA_MAG, supraMoRate, supraSeisBValue);
				
				supraSeisMFD = new IncrementalMagFreqDist(MIN_MAG, NUM_MAG, DELTA_MAG);
				for (int i=0; i<supraGR.size(); i++)
					supraSeisMFD.set(i+minMagIndex, supraGR.getY(i));
				
				GutenbergRichterMagFreqDist subGR = new GutenbergRichterMagFreqDist(
						totalOnFaultSub.getX(0), minMagIndex, DELTA_MAG, subMoRate, supraSeisBValue);
				subSeisMFD = new IncrementalMagFreqDist(MIN_MAG, NUM_MAG, DELTA_MAG);
				for (int i=0; i<subGR.size(); i++)
					subSeisMFD.set(i, subGR.getY(i));
				
				sectFullMFD = new IncrementalMagFreqDist(MIN_MAG, targetMFDNum, DELTA_MAG);
				for (int i=0; i<targetMFDNum; i++) {
					if (i < minMagIndex)
						sectFullMFD.set(i, subSeisMFD.getY(i));
					else
						sectFullMFD.set(i, supraSeisMFD.getY(i));
				}
				
				// scale target slip rates by the fraction that is supra-seismognic
				slipRates[s] = supraSlipRate;
				slipRateStdDevs[s] = supraSlipStdDev;
			} else {
				sectFullMFD = new GutenbergRichterMagFreqDist(
						MIN_MAG, targetMFDNum, DELTA_MAG, targetMoRate, supraSeisBValue);
				
				// split the target G-R into sub-seismo and supra-seismo parts
				subSeisMFD = new IncrementalMagFreqDist(MIN_MAG, NUM_MAG, DELTA_MAG);
				for (int i=0; i<minMagIndex; i++)
					subSeisMFD.set(i, sectFullMFD.getY(i));
				supraSeisMFD = new IncrementalMagFreqDist(MIN_MAG, NUM_MAG, DELTA_MAG);
				for (int i=minMagIndex; i<sectFullMFD.size(); i++)
					supraSeisMFD.set(i, sectFullMFD.getY(i));
				
				supraMoRate = supraSeisMFD.getTotalMomentRate();
				subMoRate = subSeisMFD.getTotalMomentRate();
				fractSupra = supraMoRate/targetMoRate;
				
				// scale target slip rates by the fraction that is supra-seismognic
				slipRates[s] = creepReducedSlipRate*fractSupra;
				slipRateStdDevs[s] = creepReducedSlipRateStdDev*fractSupra;
			}
			
			if (sparseGR)
				// re-distribute to only bins that actually have ruptures available
				supraSeisMFD = SparseGutenbergRichterSolver.getEquivGR(supraSeisMFD, mags,
						supraSeisMFD.getTotalMomentRate(), supraSeisBValue);
			
			
			double targetMoRateTest = supraMoRate + subMoRate;
			
			Preconditions.checkState((float)targetMoRateTest == (float)targetMoRate,
					"Partitioned moment rate doesn't equal input: %s != %s", (float)targetMoRate, (float)targetMoRateTest);
			
			fractSuprasTrack.addValue(fractSupra);
			sectFractSupras[s] = fractSupra;
			
			sectFullMFDs.add(sectFullMFD);
			sectSubSeisMFDs.add(subSeisMFD);
			sectSupraSeisMFDs.add(supraSeisMFD);
			
			// on fault supra is added later (along with uncertainty propagation)
			
			// add sub seis to total
			totalOnFaultSub.addIncrementalMagFreqDist(subSeisMFD);
			
			double relDefModStdDev = slipRateStdDevs[s]/slipRates[s];
			
			if (applyDefModelUncertainties) {
				EvenlyDiscretizedFunc defModStdDevs = null;
				if (slipRateStdDevs[s]  > 0d) {
					// use the slip rate standard deviation. this simple treatment is confirmed to be the exact same as if
					// we were to construct new GR distributions plus and minus one standard deviation and then calculate
					// a standard deviation from those bounds
					defModStdDevs = new EvenlyDiscretizedFunc(MIN_MAG, NUM_MAG, DELTA_MAG);
					for (int i=0; i<supraSeisMFD.size(); i++)
						defModStdDevs.set(i, supraSeisMFD.getY(i)*relDefModStdDev);
				}
				defModMFDStdDevs.add(defModStdDevs);
			}
			
			System.out.println("Section "+s+". targetMo="+(float)targetMoRate+"\tsupraMo="+(float)supraMoRate
					+"\tsubMo="+(float)subMoRate+"\tfractSupra="+(float)fractSupra
					+"\tsupraRate="+(float)supraSeisMFD.getTotalIncrRate());
			
//			System.out.println("TARGET\n"+sectFullMFD);
//			System.out.println("SUB\n"+subSeisMFD);
//			System.out.println("SUPRA\n"+supraSeisMFD);
			
			if (dataConstraints != null && expandUncertToDataBound != null) {
				// see if we have any data constraints that imply a different MFD for this section,
				// we will later adjust uncertainties accordingly
				List<DataSectNucleationRateEstimator> constraints = new ArrayList<>();
				for (DataSectNucleationRateEstimator constraint : dataConstraints)
					if (constraint.appliesTo(sect))
						constraints.add(constraint);
				
				if (!constraints.isEmpty()) {
					System.out.println("\tCalculating MFDs implied by "+constraints.size()+" data constraints");
					
					dataCalcFutures.add(exec.submit(new ImpliedDataCallable(sect, constraints, mags, rups, sparseGR,
							supraSeisMFD, supraMoRate, relDefModStdDev)));
				}
			}
		}
		this.sectSubSeisMFDs = new SubSeismoOnFaultMFDs(sectSubSeisMFDs);
		
		System.out.println("Fraction supra-seismogenic stats: "+fractSuprasTrack);
		
		List<IncrementalMagFreqDist> dataImpliedSectSupraSeisMFDs = null;
		if (dataCalcFutures != null) {
			dataImpliedSectSupraSeisMFDs = new ArrayList<>();
			for (int s=0; s<numSects; s++)
				dataImpliedSectSupraSeisMFDs.add(null);
			
			System.out.println("Waiting on "+dataCalcFutures.size()+" data constraint calculation futures...");
			
			for (Future<ImpliedDataCallable> future : dataCalcFutures) {
				try {
					ImpliedDataCallable call = future.get();
					dataImpliedSectSupraSeisMFDs.set(call.sect.getSectionId(), call.impliedMFD);
				} catch (InterruptedException | ExecutionException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
			}
			
			System.out.println("Done calculating data constraint MFDs");
			
			exec.shutdown();
		}
		
		if (!useExistingTargetSlipRates) {
			// give the newly computed target slip rates to the rupture set
			this.targetSlipRates = SectSlipRates.precomputed(rupSet, slipRates, slipRateStdDevs);
			rupSet.addModule(targetSlipRates);
		}
		
		this.totalOnFaultSupra = calcRegionalSupraTarget(sectSupraSeisMFDs, rupSet, null, defModMFDStdDevs,
				dataImpliedSectSupraSeisMFDs, expandUncertToDataBound, addSectCountUncertainties, defaultRelStdDev);
		this.totalOnFaultSub = totalOnFaultSub;
		
		if (constrainedRegions != null) {
			Preconditions.checkState(!constrainedRegions.isEmpty(), "Empty region list supplied");
			mfdConstraints = new ArrayList<>();
			for (Region region : constrainedRegions)
				mfdConstraints.add(calcRegionalSupraTarget(sectSupraSeisMFDs, rupSet, region, defModMFDStdDevs,
						dataImpliedSectSupraSeisMFDs, expandUncertToDataBound, addSectCountUncertainties, defaultRelStdDev));
		} else {
			mfdConstraints = List.of(this.totalOnFaultSupra);
		}
		
		// now set individual supra-seis MFD uncertainties
		this.sectSupraSeisMFDs = new ArrayList<>();
		for (int s=0; s<numSects; s++) {
			IncrementalMagFreqDist sectSupraMFD = sectSupraSeisMFDs.get(s);
			
			EvenlyDiscretizedFunc stdDevs = new EvenlyDiscretizedFunc(MIN_MAG, sectSupraMFD.size(), DELTA_MAG);
			for (int i=0; i<stdDevs.size(); i++)
				stdDevs.set(i, sectSupraMFD.getY(i)*defaultRelStdDev);
			
			if (applyDefModelUncertainties) {
				EvenlyDiscretizedFunc defModStdDevs = defModMFDStdDevs.get(s);
				if (defModStdDevs != null) {
					for (int i=0; i<stdDevs.size(); i++)
						stdDevs.set(i, Math.max(stdDevs.getY(i), defModStdDevs.getY(i)));
				}
			}
			
			UncertainIncrMagFreqDist uncertainSectSupraMFD = new UncertainIncrMagFreqDist(sectSupraMFD, stdDevs);
			
			if (dataImpliedSectSupraSeisMFDs != null) {
				IncrementalMagFreqDist impliedMFD = dataImpliedSectSupraSeisMFDs.get(s);
				if (impliedMFD != null) {
					System.out.println("Adjusting MFD for s="+s+" to account for data constraints");
					uncertainSectSupraMFD = adjustForDataImpliedBounds(uncertainSectSupraMFD, impliedMFD, expandUncertToDataBound, false);
				}
			}
			this.sectSupraSeisMFDs.add(uncertainSectSupraMFD);
		}
		
		if (rupSet.hasModule(U3LogicTreeBranch.class)) {
			// create a simple U3 total MFD, but taper it with our supra-seismogenic target
			U3LogicTreeBranch branch = rupSet.requireModule(U3LogicTreeBranch.class);
			double m5Rate = branch.getValue(TotalMag5Rate.class).getRateMag5();
			double maxMagOff = branch.getValue(MaxMagOffFault.class).getMaxMagOffFault();
			
			GutenbergRichterMagFreqDist targetGR = new GutenbergRichterMagFreqDist(MIN_MAG, NUM_MAG, DELTA_MAG);
			targetGR.setAllButTotMoRate(MIN_MAG, targetGR.getMaxX(), 1d, 1d);
			targetGR.scaleToCumRate(targetGR.getClosestXIndex(5.01), m5Rate);
			
			totalRegional = new IncrementalMagFreqDist(MIN_MAG, NUM_MAG, DELTA_MAG);
			int offMaxBin = totalRegional.getClosestXIndex(maxMagOff);
			for (int i=0; i<NUM_MAG; i++) {
				double rate = targetGR.getY(i);
				if (i > offMaxBin)
					// allow it to taper off with our supra-seismogenic target
					rate = Math.min(rate, totalOnFaultSupra.getY(i));
				totalRegional.set(i, rate);
			}
		}
	}

	public static class ImpliedDataCallable implements Callable<ImpliedDataCallable> {
		private FaultSection sect;
		private List<DataSectNucleationRateEstimator> constraints;
		private List<Double> mags;
		private List<Integer> rups;
		private boolean sparseGR;
		private IncrementalMagFreqDist supraSeisMFD;
		private double supraMoRate;
		private double relDefModStdDev;
		
		private IncrementalMagFreqDist impliedMFD;

		public ImpliedDataCallable(FaultSection sect, List<DataSectNucleationRateEstimator> constraints,
				List<Double> mags, List<Integer> rups, boolean sparseGR, IncrementalMagFreqDist supraSeisMFD,
				double supraMoRate, double relDefModStdDev) {
					this.sect = sect;
					this.constraints = constraints;
					this.mags = mags;
					this.rups = rups;
					this.sparseGR = sparseGR;
					this.supraSeisMFD = supraSeisMFD;
					this.supraMoRate = supraMoRate;
					this.relDefModStdDev = relDefModStdDev;
		}

		@Override
		public ImpliedDataCallable call() throws Exception {
			UncertainDataConstraint moRateBounds = new UncertainDataConstraint(null, supraMoRate,
					new Uncertainty(supraMoRate*relDefModStdDev));
			
			List<IncrementalMagFreqDist> impliedMFDs = new ArrayList<>();
			
			for (DataSectNucleationRateEstimator constraint : constraints) {
				// nucleation mfd implied by this constraint
				IncrementalMagFreqDist implied = constraint.estimateNuclMFD(
						sect, supraSeisMFD, rups, mags, moRateBounds, sparseGR);
				Preconditions.checkState(implied.size() == NUM_MAG);
				Preconditions.checkState((float)implied.getMinX() == (float)MIN_MAG);
				Preconditions.checkState((float)implied.getDelta() == (float)DELTA_MAG);
				impliedMFDs.add(implied);
			}
			
			if (impliedMFDs.size() > 1) {
				// need to combine them
				IncrementalMagFreqDist meanMFD = new IncrementalMagFreqDist(MIN_MAG, impliedMFDs.get(0).size(), DELTA_MAG);
				IncrementalMagFreqDist upperMFD = null;
				IncrementalMagFreqDist lowerMFD = null;
				UncertaintyBoundType boundType = null;
				int numBounded = 0;
				for (IncrementalMagFreqDist mfd : impliedMFDs) {
					Preconditions.checkState(mfd.size() == meanMFD.size());
					for (int i=0; i<meanMFD.size(); i++)
						if (mfd.getY(i) > 0)
							meanMFD.add(i, mfd.getY(i));
					if (mfd instanceof UncertainIncrMagFreqDist) {
						UncertainBoundedIncrMagFreqDist bounded;
						if (upperMFD == null) {
							upperMFD = new IncrementalMagFreqDist(MIN_MAG, meanMFD.size(), DELTA_MAG);
							lowerMFD = new IncrementalMagFreqDist(MIN_MAG, meanMFD.size(), DELTA_MAG);
							if (mfd instanceof UncertainBoundedIncrMagFreqDist) {
								bounded = (UncertainBoundedIncrMagFreqDist)mfd;
								boundType = bounded.getBoundType();
							} else {
								boundType = UncertaintyBoundType.ONE_SIGMA;
								bounded = ((UncertainIncrMagFreqDist)mfd).estimateBounds(boundType);
							}
						} else {
							bounded = ((UncertainIncrMagFreqDist)mfd).estimateBounds(boundType);
						}
						numBounded++;
						for (int i=0; i<meanMFD.size(); i++) {
							upperMFD.add(i, bounded.getUpperY(i));
							lowerMFD.add(i, bounded.getLowerY(i));
						}
					}
				}
				meanMFD.scale(1d/(double)impliedMFDs.size());
				if (upperMFD != null) {
					upperMFD.scale(1d/(double)numBounded);
					lowerMFD.scale(1d/(double)numBounded);
					for (int i=0; i<meanMFD.size(); i++) {
						upperMFD.set(i, Math.max(upperMFD.getY(i), meanMFD.getY(i)));
						lowerMFD.set(i, Math.max(0d, Math.min(lowerMFD.getY(i), meanMFD.getY(i))));
					}
					meanMFD = new UncertainBoundedIncrMagFreqDist(meanMFD, lowerMFD, upperMFD, boundType);
				}
				impliedMFD = meanMFD;
			} else {
				impliedMFD = impliedMFDs.get(0);
			}
			return this;
		}
		
	}
	
	/*
	 *  It seems that you are theoretically supposed to sum a true standard deviation in variance space, then take the
	 *  sqrt to get the summed standard deviation,
	 *  e.g. https://en.wikipedia.org/wiki/Variance#Basic_properties
	 *  	https://www.geol.lsu.edu/jlorenzo/geophysics/uncertainties/Uncertaintiespart2.html
	 *  
	 *  But, that leads to incredibly small values such that if every standard deviation is set as a relative
	 *  fraction of the target value, the sum of squared standard deviations divided by the target rate
	 *  will be far less than the input relative value.
	 *  
	 *  Because of this, I think we should just add our "standard deviations" for now
	 */
	private static final boolean SUM_VARIANCES = true;
	
	private UncertainIncrMagFreqDist calcRegionalSupraTarget(List<IncrementalMagFreqDist> sectSupraSeisMFDs,
			FaultSystemRupSet rupSet, Region region, List<EvenlyDiscretizedFunc> defModMFDStdDevs,
			List<IncrementalMagFreqDist> dataImpliedSectSupraSeisMFDs, UncertaintyBoundType expandUncertToDataBound,
			boolean addSectCountUncertainties, double defaultRelStdDev) {
		IncrementalMagFreqDist sumMFD = new IncrementalMagFreqDist(MIN_MAG, NUM_MAG, DELTA_MAG);
		
		double[] binCounts = new double[NUM_MAG];
		sumMFD.setRegion(region);
		double[] fractSectsInRegion = null;
		if (region != null)
			fractSectsInRegion = rupSet.getFractSectsInsideRegion(
				region, MFDInversionConstraint.MFD_FRACT_IN_REGION_TRACE_ONLY);
		
		// first just sum the section supra-seismogenic MFDs
		for (int s=0; s<sectSupraSeisMFDs.size(); s++) {
			double scale = fractSectsInRegion == null ? 1d : fractSectsInRegion[s];
			if (scale > 0d) {
				IncrementalMagFreqDist supraMFD = sectSupraSeisMFDs.get(s);
				// make sure we have the same gridding
				Preconditions.checkState(supraMFD.getMinX() == MIN_MAG);
				Preconditions.checkState(supraMFD.getDelta() == DELTA_MAG);
				for (int i=0; i<supraMFD.size(); i++) {
					double rate = supraMFD.getY(i);
					if (rate > 0d) {
						binCounts[i]++;
						sumMFD.add(i, rate);
					}
				}
			}
		}
		
		EvenlyDiscretizedFunc stdDevs = new EvenlyDiscretizedFunc(MIN_MAG, NUM_MAG, DELTA_MAG);
		if (defModMFDStdDevs != null) {
			// sum up deformation model implied std devs
			Preconditions.checkState(defModMFDStdDevs.size() == sectSupraSeisMFDs.size());
			for (int s=0; s<defModMFDStdDevs.size(); s++) {
				EvenlyDiscretizedFunc dmStdDevs = defModMFDStdDevs.get(s);
				double scale = fractSectsInRegion == null ? 1d : fractSectsInRegion[s];
				if (scale > 0d && dmStdDevs != null) {
					for (int i=0; i<dmStdDevs.size(); i++) {
						if (SUM_VARIANCES)
							// will take sqrt later
							stdDevs.add(i, Math.pow(dmStdDevs.getY(i), 2));
						else
							stdDevs.add(i, dmStdDevs.getY(i));
					}
				}
			}
			if (SUM_VARIANCES)
				for (int i=0; i<NUM_MAG; i++)
					stdDevs.set(i, Math.sqrt(stdDevs.getY(i)));
		}
		
		// now make sure that we're at least at the default uncertainty everywhere
		for (int i=0; i<NUM_MAG; i++)
			stdDevs.set(i, Math.max(stdDevs.getY(i), defaultRelStdDev*sumMFD.getY(i)));
		
		if (addSectCountUncertainties) {
			double refNum = fractSectsInRegion == null ? sectSupraSeisMFDs.size() : StatUtils.sum(fractSectsInRegion);
			System.out.println("Re-weighting target MFD to account for section participation uncertainties.");
			double max = StatUtils.max(binCounts);
			System.out.println("\tMax section participation: "+(float)max);
			System.out.println("\tReference section participation: "+(float)refNum);
			for (int i=0; i<binCounts.length; i++) {
				double rate = sumMFD.getY(i);
				if (binCounts[i] == 0d || rate == 0d)
					continue;
//				double relStdDev = refNum/binCounts[i];
				double relStdDev = Math.sqrt(refNum)/Math.sqrt(binCounts[i]);
				double origStdDev = stdDevs.getY(i);
				double origRel = origStdDev/rate;
				System.out.println("\t\tRelative particpation std dev for M="
						+(float)sumMFD.getX(i)+", binCount="+(float)binCounts[i]+": "+(float)relStdDev
						+"\torigStdDev="+(float)origStdDev+"\torigRelStdDev="+(float)origRel+"\ttotRel="+(float)(origRel*relStdDev));
				// scale existing std dev
				stdDevs.set(i, origStdDev*relStdDev);
			}
		}
		
		UncertainIncrMagFreqDist uncertainMFD = new UncertainIncrMagFreqDist(sumMFD, stdDevs);
		if (dataImpliedSectSupraSeisMFDs != null) {
			// compute data implied std devs
			Preconditions.checkState(dataImpliedSectSupraSeisMFDs.size() == sectSupraSeisMFDs.size());
			IncrementalMagFreqDist sumImpliedMFD = new IncrementalMagFreqDist(MIN_MAG, NUM_MAG, DELTA_MAG);
			IncrementalMagFreqDist sumLowerMFD = new IncrementalMagFreqDist(MIN_MAG, NUM_MAG, DELTA_MAG);
			IncrementalMagFreqDist sumUpperMFD = new IncrementalMagFreqDist(MIN_MAG, NUM_MAG, DELTA_MAG);
			
			for (int s=0; s<dataImpliedSectSupraSeisMFDs.size(); s++) {
				double scale = fractSectsInRegion == null ? 1d : fractSectsInRegion[s];
				if (scale > 0d) {
					IncrementalMagFreqDist impliedMFD = dataImpliedSectSupraSeisMFDs.get(s);
					if (impliedMFD == null)
						// no data constraint for this section, use the regular MFD
						impliedMFD = sectSupraSeisMFDs.get(s);
					// make sure we have the same gridding
					Preconditions.checkState(impliedMFD.getMinX() == MIN_MAG);
					Preconditions.checkState(impliedMFD.getDelta() == DELTA_MAG);
					for (int i=0; i<impliedMFD.size(); i++)
						sumImpliedMFD.add(i, impliedMFD.getY(i));
					if (impliedMFD instanceof UncertainBoundedIncrMagFreqDist) {
						UncertainBoundedIncrMagFreqDist bounded = (UncertainBoundedIncrMagFreqDist)impliedMFD;
						for (int i=0; i<impliedMFD.size(); i++) {
							sumLowerMFD.add(i, bounded.getLowerY(i));
							sumUpperMFD.add(i, bounded.getUpperY(i));
						}
					} else {
						for (int i=0; i<impliedMFD.size(); i++) {
							double rate = impliedMFD.getY(i);
							sumLowerMFD.add(i, rate);
							sumUpperMFD.add(i, rate);
						}
					}
				}
			}
			
			UncertainBoundedIncrMagFreqDist impliedMFD = new UncertainBoundedIncrMagFreqDist(
					sumImpliedMFD, sumLowerMFD, sumUpperMFD, expandUncertToDataBound);
			
			// adjust the standard deviations of the regional mfd to reach the implied MFD
			System.out.println("Adjusting regional MFD to match data bounds");
			uncertainMFD = adjustForDataImpliedBounds(uncertainMFD, impliedMFD, expandUncertToDataBound, true);
		}
		
		System.out.println("Final regional MFD:");
		boolean first = true;
		for (int i=0; i<uncertainMFD.size(); i++) {
			double x = uncertainMFD.getX(i);
			double y = uncertainMFD.getY(i);
			if (first && y == 0d)
				continue;
			first = false;
			double sd = uncertainMFD.getStdDev(i);
			System.out.println("\tM="+(float)x+"\tRate="+(float)y+"\tStdDev="+(float)sd+"\tRelStdDev="+(float)(sd/y));
		}
		return uncertainMFD;
	}
	
	private static UncertainBoundedIncrMagFreqDist adjustForDataImpliedBounds(UncertainIncrMagFreqDist mfd,
			IncrementalMagFreqDist impliedMFD, UncertaintyBoundType expandUncertToDataBound, boolean verbose) {
		Preconditions.checkState(mfd.size() == impliedMFD.size());
		Preconditions.checkState((float)mfd.getMinX() == (float)impliedMFD.getMinX());
		Preconditions.checkState((float)mfd.getDelta() == (float)impliedMFD.getDelta());
		
		UncertainBoundedIncrMagFreqDist boundedImplied = null;
		if (impliedMFD instanceof UncertainBoundedIncrMagFreqDist)
			boundedImplied = (UncertainBoundedIncrMagFreqDist)impliedMFD;
		else if (impliedMFD instanceof UncertainIncrMagFreqDist)
			boundedImplied = ((UncertainIncrMagFreqDist)impliedMFD).estimateBounds(expandUncertToDataBound);

		
		EvenlyDiscretizedFunc stdDevs = new EvenlyDiscretizedFunc(mfd.getMinX(), mfd.size(), mfd.getDelta());

		UncertainBoundedIncrMagFreqDist boundedInput = mfd.estimateBounds(UncertaintyBoundType.ONE_SIGMA);
		IncrementalMagFreqDist lower = boundedInput.getLower();
		IncrementalMagFreqDist upper = boundedInput.getUpper();
		EvenlyDiscretizedFunc inputStdDevs = boundedInput.getStdDevs();
		
		MinMaxAveTracker stdDevTrack = new MinMaxAveTracker();
		MinMaxAveTracker relStdDevTrack = new MinMaxAveTracker();
		for (int i=0; i<NUM_MAG; i++) {
			if (mfd.getY(i) == 0)
				continue;
			double mag = mfd.getX(i);
			double rate = mfd.getY(i);
			double[] dataRates;
			if (boundedImplied == null)
				dataRates = new double[] { impliedMFD.getY(i) };
			else
				dataRates = new double[] { impliedMFD.getY(i), boundedImplied.getLowerY(i), boundedImplied.getUpperY(i) };
			double minImpliedStdDev = Double.POSITIVE_INFINITY;
			double closestData = Double.NaN;
			for (double dataRate : dataRates) {
				double diff = Math.abs(dataRate - rate);
				double impliedStdDev = expandUncertToDataBound.estimateStdDev(rate, rate-diff, rate+diff);
				if (impliedStdDev < minImpliedStdDev) {
					minImpliedStdDev = impliedStdDev;
					closestData = dataRate;
				}
			}
			if (closestData > rate) {
				// adjust the upper bound
				upper.set(i, Math.max(upper.getY(i), rate+minImpliedStdDev));
			} else {
				lower.set(i, Math.max(0, Math.min(lower.getY(i), rate-minImpliedStdDev)));
			}
			double relStdDev = minImpliedStdDev/rate;
			if (verbose) System.out.println("\tM="+(float)mag+"\trate="+(float)rate+"\tdataRate="+(float)dataRates[0]
					+"\tdataBounds=["+(float)dataRates[1]+","+(float)dataRates[2]+"]\timplStdDev="+minImpliedStdDev
					+"\timplRelStdDev="+(float)relStdDev);
			stdDevTrack.addValue(minImpliedStdDev);
			relStdDevTrack.addValue(relStdDev);
			stdDevs.set(i, Math.max(minImpliedStdDev, inputStdDevs.getY(i)));
		}
		System.out.println("\tImplied std dev range:\t["+(float)stdDevTrack.getMin()+","+(float)stdDevTrack.getMax()
				+"], avg="+(float)stdDevTrack.getAverage()+"\t\trel range:\t["+(float)relStdDevTrack.getMin()
				+","+(float)relStdDevTrack.getMax()+"], avg="+(float)relStdDevTrack.getAverage());
		return new UncertainBoundedIncrMagFreqDist(mfd, lower, upper, UncertaintyBoundType.ONE_SIGMA);
	}

	public double getSupraSeisBValue() {
		return supraSeisBValue;
	}

	public List<UncertainIncrMagFreqDist> getSectSupraSeisNuclMFDs() {
		return sectSupraSeisMFDs;
	}

	@Override
	public String getName() {
		return "Supra-Seis B="+(float)supraSeisBValue+" Inversion Target MFDs";
	}

	@Override
	public SubModule<FaultSystemRupSet> copy(FaultSystemRupSet newParent) throws IllegalStateException {
		throw new UnsupportedOperationException();
//		List<Region> regions = null;
//		if (mfdConstraints.size() > 1 || mfdConstraints.get(0).getRegion() != null) {
//			regions = new ArrayList<>();
//			for (IncrementalMagFreqDist constraint : mfdConstraints) {
//				Preconditions.checkNotNull(constraint.getRegion());
//				regions.add(constraint.getRegion());
//			}
//		}
//		return new InversionTargetMFDsFromBValAndDefModel(newParent, supraSeisBValue, regions);
	}

	@Override
	public IncrementalMagFreqDist getTotalRegionalMFD() {
		return totalRegional;
	}

	@Override
	public IncrementalMagFreqDist getTotalOnFaultSupraSeisMFD() {
		return totalOnFaultSupra;
	}

	@Override
	public IncrementalMagFreqDist getTotalOnFaultSubSeisMFD() {
		return totalOnFaultSub;
	}

	@Override
	public IncrementalMagFreqDist getTrulyOffFaultMFD() {
		return null;
	}

	@Override
	public List<? extends IncrementalMagFreqDist> getMFD_Constraints() {
		return mfdConstraints;
	}

	@Override
	public SubSeismoOnFaultMFDs getOnFaultSubSeisMFDs() {
		return sectSubSeisMFDs;
	}

	@Override
	public void writeToArchive(ZipOutputStream zout, String entryPrefix) throws IOException {
		new Precomputed(this).writeToArchive(zout, entryPrefix);
	}

	@Override
	public void initFromArchive(ZipFile zip, String entryPrefix) throws IOException {
		throw new IllegalStateException("Should be deserialized as Precomputed class");
	}

	@Override
	public Class<? extends ArchivableModule> getLoadingClass() {
		return InversionTargetMFDs.Precomputed.class;
	}
	
	public static void main(String[] args) throws IOException {
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(
				new File("/data/kevin/markdown/inversions/fm3_1_u3ref_uniform_reproduce_ucerf3.zip"));
//				new File("/data/kevin/markdown/inversions/fm3_1_u3ref_uniform_coulomb.zip"));
		
		InversionTargetMFDs origTargets = rupSet.getModule(InversionTargetMFDs.class);
		SectSlipRates origSlips = rupSet.getModule(SectSlipRates.class);
		
		RupSetMapMaker mapMaker = new RupSetMapMaker(rupSet, new CaliforniaRegions.RELM_TESTING());
		CPT cpt = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(0d, 1d);
		
		double b = 0.8;
		boolean sparseGR = true;
		double defaultRelStdDev = 0.1;
		boolean applyDefModelUncertainties = true;
		boolean addSectCountUncertainties = false;
		boolean useExistingTargetSlipRates = false;
//		List<DataSectNucleationRateEstimator> dataConstraints = null;
		List<DataSectNucleationRateEstimator> dataConstraints = new ArrayList<>();
		dataConstraints.add(new APrioriSectNuclEstimator(
				rupSet, UCERF3InversionInputGenerator.findParkfieldRups(rupSet), 1d/25d, 0.1d/25d));
		dataConstraints.addAll(PaleoSectNuclEstimator.buildPaleoEstimates(rupSet, true));
		UncertaintyBoundType expandUncertToDataBound = UncertaintyBoundType.ONE_SIGMA;
		InversionTargetMFDsFromBValAndDefModel target = new InversionTargetMFDsFromBValAndDefModel(
				rupSet, b, sparseGR, defaultRelStdDev, applyDefModelUncertainties, addSectCountUncertainties,
				useExistingTargetSlipRates, dataConstraints, expandUncertToDataBound, null);
		rupSet.addModule(target);
		
		mapMaker.plotSectScalars(target.sectFractSupras, cpt, "New Fraction of Moment Supra-Seismogenic");
		
		File outputDir = new File("/tmp");
		mapMaker.plot(outputDir, "supra_seis_fracts", " ");
		
		double[] u3FractSupra = new double[rupSet.getNumSections()];
		MinMaxAveTracker track = new MinMaxAveTracker();
		for (int s=0; s<u3FractSupra.length; s++) {
			double creepReduced = rupSet.getFaultSectionData(s).getReducedAveSlipRate()*1e-3;
			double targetSlip = origSlips.getSlipRate(s);
			double fract = targetSlip/creepReduced;
			track.addValue(fract);
			u3FractSupra[s] = fract;
		}
		System.out.println("U3 fracts: "+track);
		
		mapMaker.plotSectScalars(u3FractSupra, cpt, "U3 Fraction of Moment Supra-Seismogenic");

		mapMaker.plot(outputDir, "supra_seis_fracts_u3", " ");
		
		SolMFDPlot plot = new SolMFDPlot();
		File mfdOutput = new File(outputDir, "test_target_mfds");
		Preconditions.checkState(mfdOutput.exists() || mfdOutput.mkdir());
		plot.writePlot(rupSet, null, "Test Model", mfdOutput);
	}

}
