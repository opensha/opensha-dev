package scratch.nshm23.targetMFDs;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.uncertainty.UncertainIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertaintyBoundType;
import org.opensha.commons.geo.Region;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ClassUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.modules.ArchivableModule;
import org.opensha.commons.util.modules.SubModule;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.MFDInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.ModSectMinMags;
import org.opensha.sha.earthquake.faultSysSolution.modules.SectSlipRates;
import org.opensha.sha.earthquake.faultSysSolution.modules.SubSeismoOnFaultMFDs;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.SolMFDPlot;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
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
	
	private static boolean SPARSE_GR_DEFAULT = true;

	public InversionTargetMFDsFromBValAndDefModel(FaultSystemRupSet rupSet, double supraSeisBValue) {
		this(rupSet, supraSeisBValue, SPARSE_GR_DEFAULT, TARGET_REL_STD_DEV_DEFAULT,
				APPLY_DEF_MODEL_UNCERTAINTIES_DEFAULT, WEIGHT_STD_DEV_BY_PARTICIPATION_DEFAULT, null);
	}
	
	public InversionTargetMFDsFromBValAndDefModel(FaultSystemRupSet rupSet, double supraSeisBValue, boolean sparseGR,
			double defaultRelStdDev, boolean applyDefModelUncertainties, boolean addSectCountUncertainties,
			List<Region> constrainedRegions) {
		this(rupSet, supraSeisBValue, sparseGR, defaultRelStdDev, applyDefModelUncertainties,
				addSectCountUncertainties, null, null, constrainedRegions);
	}
	
	public InversionTargetMFDsFromBValAndDefModel(FaultSystemRupSet rupSet, double supraSeisBValue, boolean sparseGR,
			double defaultRelStdDev, boolean applyDefModelUncertainties, boolean addSectCountUncertainties,
			List<DataSectNucleationRateEstimator> dataConstraints, UncertaintyBoundType expandUncertToDataBound,
			List<Region> constrainedRegions) {
		super(rupSet);
		this.supraSeisBValue = supraSeisBValue;
		
		ModSectMinMags minMags = rupSet.getModule(ModSectMinMags.class);
		
		int numSects = rupSet.getNumSections();
		double[] slipRates = new double[numSects];
		double[] slipRateStdDevs = new double[numSects];
		sectFullMFDs = new ArrayList<>();
		List<IncrementalMagFreqDist> sectSubSeisMFDs = new ArrayList<>();
		sectSupraSeisMFDs = new ArrayList<>();
		
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
			// make sure we actually have a rupture at that magnitude, otherwise there can be empty bins without
			// any rupture at/above the section minimum magnitude but below the first rupture's bin
			double minAbove = Double.POSITIVE_INFINITY;
			for (int r : rupSet.getRupturesForSection(s)) {
				double mag = rupSet.getMagForRup(r);
				if ((float)mag >= (float)sectMinMag) {
					minAbove = Math.min(mag, minAbove);
					mags.add(mag);
				}
			}
			sectMinMag = minAbove;
			double sectMaxMag = rupSet.getMaxMagForSection(s);
			
			// construct a full G-R including sub-seismogenic ruptures
			int minMagIndex = totalOnFaultSub.getClosestXIndex(sectMinMag);
			int maxMagIndex = totalOnFaultSub.getClosestXIndex(sectMaxMag);
			int targetMFDNum = maxMagIndex+1;
			GutenbergRichterMagFreqDist sectFullMFD = new GutenbergRichterMagFreqDist(
					MIN_MAG, targetMFDNum, DELTA_MAG, targetMoRate, supraSeisBValue);
			
			// split the target G-R into sub-seismo and supra-seismo parts
			IncrementalMagFreqDist subSeisMFD = new IncrementalMagFreqDist(MIN_MAG, NUM_MAG, DELTA_MAG);
			for (int i=0; i<minMagIndex; i++)
				subSeisMFD.set(i, sectFullMFD.getY(i));
			IncrementalMagFreqDist supraSeisMFD = new IncrementalMagFreqDist(MIN_MAG, NUM_MAG, DELTA_MAG);
			for (int i=minMagIndex; i<sectFullMFD.size(); i++)
				supraSeisMFD.set(i, sectFullMFD.getY(i));
			
			if (sparseGR)
				// re-distribute to only bins that actually have ruptures available
				supraSeisMFD = SparseGutenbergRichterSolver.getEquivGR(supraSeisMFD, mags,
						supraSeisMFD.getTotalMomentRate(), supraSeisBValue);
			
			double supraMoRate = supraSeisMFD.getTotalMomentRate();
			double subMoRate = subSeisMFD.getTotalMomentRate();
			double targetMoRateTest = supraMoRate + subMoRate;
			
			double fractSupra = supraMoRate/targetMoRate;
			
			fractSuprasTrack.addValue(fractSupra);
			sectFractSupras[s] = fractSupra;
			
			// scale target slip rates by the fraction that is supra-seismognic
			slipRates[s] = creepReducedSlipRate*fractSupra;
			slipRateStdDevs[s] = creepReducedSlipRateStdDev*fractSupra;
			
			UncertainIncrMagFreqDist uncertSupraSeisMFD;
			if (applyDefModelUncertainties && slipRateStdDevs[s]  > 0d) {
				// use the slip rate standard deviation. this simple treatment is confirmed to be the exact same as if
				// we were to construct new GR distributions plus and minus one standard deviation and then calculate
				// a standard deviation from those bounds
				uncertSupraSeisMFD = UncertainIncrMagFreqDist.constantRelStdDev(supraSeisMFD, slipRateStdDevs[s]/slipRates[s]);
			} else {
				// use default relative standard deviation
				uncertSupraSeisMFD = UncertainIncrMagFreqDist.constantRelStdDev(supraSeisMFD, defaultRelStdDev);
			}
			
			System.out.println("Section "+s+". targetMo="+(float)targetMoRate+"\tsupraMo="+(float)supraMoRate
					+"\tsubMo="+(float)subMoRate+"\tfractSupra="+(float)fractSupra+"\tsupraRate="
					+(float)supraSeisMFD.getTotalIncrRate()+"\trelStdDevs: "+getStdDevsStr(uncertSupraSeisMFD, true));
			
//			System.out.println("TARGET\n"+sectFullMFD);
//			System.out.println("SUB\n"+subSeisMFD);
//			System.out.println("SUPRA\n"+supraSeisMFD);
			
			Preconditions.checkState((float)targetMoRateTest == (float)targetMoRate,
					"Partitioned moment rate doesn't equal input: %s != %s", (float)targetMoRate, (float)targetMoRateTest);
			
			if (dataConstraints != null && expandUncertToDataBound != null) {
				List<DataSectNucleationRateEstimator> constraints = new ArrayList<>();
				for (DataSectNucleationRateEstimator constraint : dataConstraints)
					if (constraint.appliesTo(sect))
						constraints.add(constraint);
				
				if (!constraints.isEmpty()) {
					System.out.println("*** Adjusting uncertainties to include "+constraints.size()+" data constraints");
					EvenlyDiscretizedFunc modStdDevs = uncertSupraSeisMFD.getStdDevs().deepClone();
					
					for (DataSectNucleationRateEstimator constraint : constraints) {
						// nucleation rate implied by this constraint
						double impliedNuclRate = constraint.estimateNuclRate(sect, uncertSupraSeisMFD);
						IncrementalMagFreqDist dataTarget;
						if (minMagIndex == maxMagIndex) {
							// simple
							dataTarget = new IncrementalMagFreqDist(minMagIndex, 1, DELTA_MAG);
							dataTarget.set(0, impliedNuclRate);
						} else {
							// create a G-R with that nucleation rate
							GutenbergRichterMagFreqDist impliedGR = new GutenbergRichterMagFreqDist(
									supraSeisMFD.getX(minMagIndex), 1+maxMagIndex-minMagIndex, DELTA_MAG);
							impliedGR.setAllButBvalue(impliedGR.getX(0), impliedGR.getX(impliedGR.size()-1),
									targetMoRate, impliedNuclRate);
							double impliedB = impliedGR.get_bValue();
							System.out.println("\tConstraint of type "
									+ClassUtils.getClassNameWithoutPackage(constraint.getClass())
									+" implies nucl rate of: "+(float)impliedNuclRate+", b-value="+(float)impliedB);
							dataTarget = impliedGR;
							if (sparseGR)
								dataTarget = SparseGutenbergRichterSolver.getEquivGR(impliedGR, mags,
										impliedGR.getTotalMomentRate(), impliedB);
						}
						for (int i=0; i<dataTarget.size(); i++) {
							int index = minMagIndex+i;
							Preconditions.checkState((float)dataTarget.getX(i) == (float)supraSeisMFD.getX(index));
							double primaryEst = supraSeisMFD.getY(index);
							double curStdDev = modStdDevs.getY(index);
							double dataEst = dataTarget.getY(i);
							double diff = Math.abs(dataEst - primaryEst);
							double impliedStdDev = expandUncertToDataBound.estimateStdDev(primaryEst, primaryEst-diff, primaryEst+diff);
							if (impliedStdDev > curStdDev)
								modStdDevs.set(index, impliedStdDev);
						}
					}
					uncertSupraSeisMFD = new UncertainIncrMagFreqDist(supraSeisMFD, modStdDevs);
					System.out.println("\tModified relStdDevs: "+getStdDevsStr(uncertSupraSeisMFD, true));
				}
			}
			
			sectFullMFDs.add(sectFullMFD);
			sectSubSeisMFDs.add(subSeisMFD);
			sectSupraSeisMFDs.add(uncertSupraSeisMFD);
			
			// on fault supra is added later (along with uncertainty propagation)
			
			// add sub seis to total
			totalOnFaultSub.addIncrementalMagFreqDist(subSeisMFD);
		}
		this.sectSubSeisMFDs = new SubSeismoOnFaultMFDs(sectSubSeisMFDs);
		
		System.out.println("Fraction supra-seismogenic stats: "+fractSuprasTrack);
		
		// give the newly computed target slip rates to the rupture set
		this.targetSlipRates = SectSlipRates.precomputed(rupSet, slipRates, slipRateStdDevs);
		rupSet.addModule(targetSlipRates);
		
		this.totalOnFaultSupra = calcRegionalSupraTarget(sectSupraSeisMFDs, rupSet, null, addSectCountUncertainties);
		this.totalOnFaultSub = totalOnFaultSub;
		
		if (constrainedRegions != null) {
			Preconditions.checkState(!constrainedRegions.isEmpty(), "Empty region list supplied");
			mfdConstraints = new ArrayList<>();
			for (Region region : constrainedRegions)
				mfdConstraints.add(calcRegionalSupraTarget(sectSupraSeisMFDs, rupSet, region, addSectCountUncertainties));
		} else {
			mfdConstraints = List.of(this.totalOnFaultSupra);
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
	
	private static String getStdDevsStr(UncertainIncrMagFreqDist mfd, boolean relative) {
		MinMaxAveTracker track = new MinMaxAveTracker();
		for (int i=0; i<mfd.size(); i++) {
			double val = mfd.getY(i);
			if (val == 0d)
				continue;
			double sd = mfd.getStdDev(i);
			if (relative)
				track.addValue(sd/val);
			else
				track.addValue(sd);
		}
		return "["+(float)track.getMin()+","+(float)track.getMax()+"], avg="+(float)track.getAverage();
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
	private static final boolean SUM_VARIANCES = false;
	
	private static UncertainIncrMagFreqDist calcRegionalSupraTarget(List<UncertainIncrMagFreqDist> supraSeisMFDs,
			FaultSystemRupSet rupSet, Region region, boolean addSectCountUncertainties) {
		IncrementalMagFreqDist sumMFD = new IncrementalMagFreqDist(MIN_MAG, NUM_MAG, DELTA_MAG);
		EvenlyDiscretizedFunc stdDevs = new EvenlyDiscretizedFunc(MIN_MAG, NUM_MAG, DELTA_MAG);
		double[] binCounts = new double[NUM_MAG];
		sumMFD.setRegion(region);
		double[] fractSectsInRegion = null;
		if (region != null)
			fractSectsInRegion = rupSet.getFractSectsInsideRegion(
				region, MFDInversionConstraint.MFD_FRACT_IN_REGION_TRACE_ONLY);
		for (int s=0; s<supraSeisMFDs.size(); s++) {
			double scale = fractSectsInRegion == null ? 1d : fractSectsInRegion[s];
			if (scale > 0d) {
				UncertainIncrMagFreqDist supraMFD = supraSeisMFDs.get(s);
				// make sure we have the same gridding
				Preconditions.checkState(supraMFD.getMinX() == MIN_MAG);
				Preconditions.checkState(supraMFD.getDelta() == DELTA_MAG);
				for (int i=0; i<supraMFD.size(); i++) {
					double rate = supraMFD.getY(i);
					if (rate > 0d) {
						binCounts[i]++;
						sumMFD.add(i, rate);
						double stdDev = supraMFD.getStdDev(i);
						Preconditions.checkState(Double.isFinite(stdDev),
								"Bad std dev for m=%s & rate=%s: %s", supraMFD.getX(i), rate, stdDev);
						if (SUM_VARIANCES)
							// will take sqrt later
							stdDevs.add(i, Math.pow(supraMFD.getStdDev(i), 2));
						else
							stdDevs.add(i, supraMFD.getStdDev(i));
					}
				}
			}
		}
		if (SUM_VARIANCES)
			for (int i=0; i<NUM_MAG; i++)
				stdDevs.set(i, Math.sqrt(stdDevs.getY(i)));
		if (addSectCountUncertainties) {
			double refNum = fractSectsInRegion == null ? supraSeisMFDs.size() : StatUtils.sum(fractSectsInRegion);
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
		System.out.println("Final regional MFD:");
		boolean first = true;
		for (int i=0; i<sumMFD.size(); i++) {
			double x = sumMFD.getX(i);
			double y = sumMFD.getY(i);
			if (first && y == 0d)
				continue;
			first = false;
			double sd = stdDevs.getY(i);
			System.out.println("\tM="+(float)x+"\tRate="+(float)y+"\tStdDev="+(float)sd+"\tRelStdDev="+(float)(sd/y));
		}
		return new UncertainIncrMagFreqDist(sumMFD, stdDevs);
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
		boolean addSectCountUncertainties = true;
//		List<DataSectNucleationRateEstimator> dataConstraints = null;
		List<DataSectNucleationRateEstimator> dataConstraints = List.of(new APrioriSecnNuclEstimator(
				rupSet, UCERF3InversionInputGenerator.findParkfieldRups(rupSet), 1d/25d));
		UncertaintyBoundType expandUncertToDataBound = UncertaintyBoundType.TWO_SIGMA;
		InversionTargetMFDsFromBValAndDefModel target = new InversionTargetMFDsFromBValAndDefModel(
				rupSet, b, sparseGR, defaultRelStdDev, applyDefModelUncertainties, addSectCountUncertainties,
				dataConstraints, expandUncertToDataBound, null);
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
