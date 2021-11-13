package scratch.nshm23.targetMFDs;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.uncertainty.UncertainIncrMagFreqDist;
import org.opensha.commons.geo.Region;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
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
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;

import scratch.UCERF3.enumTreeBranches.MaxMagOffFault;
import scratch.UCERF3.enumTreeBranches.TotalMag5Rate;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;

public class InversionTargetMFDsFromBValAndDefModel extends InversionTargetMFDs implements ArchivableModule {
	
	// inputs
	private double supraSeisBValue;
	
	// things we compute here
	private SectSlipRates targetSlipRates;
	private List<IncrementalMagFreqDist> sectFullMFDs;
	private SubSeismoOnFaultMFDs sectSubSeisMFDs;
	private List<IncrementalMagFreqDist> sectSupraSeisMFDs;
	private SummedMagFreqDist totalOnFaultSupra;
	private SummedMagFreqDist totalOnFaultSub;
	private List<IncrementalMagFreqDist> mfdConstraints;
	
	private IncrementalMagFreqDist totalRegional;
	
	private double[] sectFractSupras;
	
	// discretization parameters for MFDs
	public final static double MIN_MAG = 0.05;
	public final static double MAX_MAG = 8.95;
	public final static int NUM_MAG = 90;
	public final static double DELTA_MAG = 0.1;

	public InversionTargetMFDsFromBValAndDefModel(FaultSystemRupSet rupSet, double supraSeisBValue) {
		this(rupSet, supraSeisBValue, null);
	}
	
	public InversionTargetMFDsFromBValAndDefModel(FaultSystemRupSet rupSet, double supraSeisBValue,
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
		
		totalOnFaultSupra = new SummedMagFreqDist(MIN_MAG, NUM_MAG, DELTA_MAG);
		totalOnFaultSub = new SummedMagFreqDist(MIN_MAG, NUM_MAG, DELTA_MAG);
		
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
			double sectMinMag = minMags == null ? rupSet.getMinMagForSection(s) : minMags.getMinMagForSection(s);
			double sectMaxMag = rupSet.getMaxMagForSection(s);
			
			// construct a full G-R including sub-seismogenic ruptures
			int minMagIndex = totalOnFaultSupra.getClosestXIndex(sectMinMag);
			int maxMagIndex = totalOnFaultSupra.getClosestXIndex(sectMaxMag);
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
			
			double supraMoRate = supraSeisMFD.getTotalMomentRate();
			double subMoRate = subSeisMFD.getTotalMomentRate();
			double targetMoRateTest = supraMoRate + subMoRate;
			
			double fractSupra = supraMoRate/targetMoRate;
			
			System.out.println("Section "+s+". target="+(float)targetMoRate+"\tsupra="+(float)supraMoRate
					+"\tsup="+(float)subMoRate+"\tfractSupra="+(float)fractSupra);
			
			fractSuprasTrack.addValue(fractSupra);
			sectFractSupras[s] = fractSupra;
			
//			System.out.println("TARGET\n"+sectFullMFD);
//			System.out.println("SUB\n"+subSeisMFD);
//			System.out.println("SUPRA\n"+supraSeisMFD);
			
			Preconditions.checkState((float)targetMoRateTest == (float)targetMoRate,
					"Partitioned moment rate doesn't equal input: %s != %s", (float)targetMoRate, (float)targetMoRateTest);
			
			sectFullMFDs.add(sectFullMFD);
			sectSubSeisMFDs.add(subSeisMFD);
			sectSupraSeisMFDs.add(supraSeisMFD);
			
			if (supraSeisMFD.size() == 1) {
				Preconditions.checkState(minMagIndex == maxMagIndex);
				totalOnFaultSupra.add(minMagIndex, supraSeisMFD.getY(0));
			} else {
				totalOnFaultSupra.addIncrementalMagFreqDist(supraSeisMFD);
			}
			totalOnFaultSub.addIncrementalMagFreqDist(subSeisMFD);
			
			// scale target slip rates by the fraction that is supra-seismognic
			slipRates[s] = creepReducedSlipRate*fractSupra;
			slipRateStdDevs[s] = creepReducedSlipRateStdDev*fractSupra;
		}
		this.sectSubSeisMFDs = new SubSeismoOnFaultMFDs(sectSubSeisMFDs);
		
		System.out.println("Fraction supra-seismogenic stats: "+fractSuprasTrack);
		
		// give the newly computed target slip rates to the rupture set
		this.targetSlipRates = SectSlipRates.precomputed(rupSet, slipRates, slipRateStdDevs);
		rupSet.addModule(targetSlipRates);
		
		if (constrainedRegions != null) {
			Preconditions.checkState(!constrainedRegions.isEmpty(), "Empty region list supplied");
			mfdConstraints = new ArrayList<>();
			for (Region region : constrainedRegions) {
				SummedMagFreqDist regMFD = new SummedMagFreqDist(MIN_MAG, NUM_MAG, DELTA_MAG);
				regMFD.setRegion(region);
				double[] fractSectsInRegion = rupSet.getFractSectsInsideRegion(
						region, MFDInversionConstraint.MFD_FRACT_IN_REGION_TRACE_ONLY);
				for (int s=0; s<numSects; s++) {
					if (fractSectsInRegion[s] > 0d) {
						IncrementalMagFreqDist supraMFD = sectSubSeisMFDs.get(s);
						if (fractSectsInRegion[s] != 1d) {
							supraMFD = supraMFD.deepClone();
							supraMFD.scale(fractSectsInRegion[s]);
						}
						regMFD.addIncrementalMagFreqDist(supraMFD);
					}
				}
				mfdConstraints.add(regMFD);
			}
		} else {
			mfdConstraints = List.of(totalOnFaultSupra);
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

	public double getSupraSeisBValue() {
		return supraSeisBValue;
	}

	public List<IncrementalMagFreqDist> getSectSupraSeisNuclMFDs() {
		return sectSupraSeisMFDs;
	}

	@Override
	public String getName() {
		return "Supra-Seis B="+(float)supraSeisBValue+" Inversion Target MFDs";
	}

	@Override
	public SubModule<FaultSystemRupSet> copy(FaultSystemRupSet newParent) throws IllegalStateException {
		List<Region> regions = null;
		if (mfdConstraints.size() > 1 || mfdConstraints.get(0).getRegion() != null) {
			regions = new ArrayList<>();
			for (IncrementalMagFreqDist constraint : mfdConstraints) {
				Preconditions.checkNotNull(constraint.getRegion());
				regions.add(constraint.getRegion());
			}
		}
		return new InversionTargetMFDsFromBValAndDefModel(newParent, supraSeisBValue, regions);
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
		InversionTargetMFDsFromBValAndDefModel target = new InversionTargetMFDsFromBValAndDefModel(rupSet, b);
		
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
	}

}
