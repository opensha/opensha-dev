package scratch.kevin.nshm23.dmCovarianceTests;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;

import org.apache.commons.math3.random.Well19937c;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.RupSetDeformationModel;
import org.opensha.sha.earthquake.faultSysSolution.RupSetFaultModel;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SectionDistanceAzimuthCalculator;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

import scratch.kevin.nshm23.dmCovarianceTests.SectionCovarianceSampler.CachedDecomposition;

public abstract class DefModSamplingEnabledInvConfig extends NSHM23_InvConfigFactory {
	
	protected int interpSkip;
	
	private abstract static class PreCachedCorrConfig extends DefModSamplingEnabledInvConfig {
		
		private String cachePrefix;

		protected PreCachedCorrConfig(int interpSkip, String cachePrefix) {
			super(interpSkip);
			this.cachePrefix = cachePrefix;
		}
		
		protected SectionCovarianceSampler buildCovarianceSampler(FaultSystemRupSet rupSet,
				RupSetFaultModel fm, List<? extends FaultSection> subSects)
				throws IOException {
			// need to build one
			Preconditions.checkNotNull(cacheDir, "cache dir should not be null");
			File subDir = getCacheSubDir(fm);
			Preconditions.checkState(subDir.exists(), "cache subdir doesn't exist: %s", subDir.getAbsolutePath());
			Optional<SectionCovarianceSampler> sampler = SectionCovarianceSampler.loadCached(subDir, subSects, cachePrefix, interpSkip);
			Preconditions.checkState(sampler.isPresent(), "No cache found!\n\tDir: %s\n\tPrefix: ", subDir.getAbsolutePath(), cachePrefix);
			return sampler.get();
		}
	}
	
	public static class CachedConnDistCorrInterp1 extends PreCachedCorrConfig {
		
		private static int INTERP_SKIP = 1;
		private static String PREFIX = "covariance_cache_conn_corr_dist100.0km_zeroCoeff0.95_negCorrDist30.0km_5548_sects_19407_trace_locs_683310349448_area_interp1";

		public CachedConnDistCorrInterp1() {
			super(INTERP_SKIP, PREFIX);
		}
		
	}
	
	public static class ConnDistB0p5MidSegCorrInterp1 extends DefModSamplingEnabledInvConfig {
		
		private static int INTERP_SKIP = 1;
		private static double B_VAL = 0.5;
		private static NSHM23_SegmentationModels SEG_MODEL = NSHM23_SegmentationModels.MID;
		private static double MAX_DIST = 200d;
		private static double ZERO_DIST_COEFF = 0.95;
		private static double NEG_CORR_MAX_DIST = 30d;
		private static boolean FROM_MFDS = true;
		
		public ConnDistB0p5MidSegCorrInterp1() {
			super(INTERP_SKIP);
		}

		@Override
		protected SectionCovarianceSampler buildCovarianceSampler(FaultSystemRupSet rupSet, RupSetFaultModel fm,
				List<? extends FaultSection> subSects) throws IOException {
			SectionDistanceAzimuthCalculator distCalc = rupSet.getModule(SectionDistanceAzimuthCalculator.class);
			if (distCalc == null) {
				distCalc = new SectionDistanceAzimuthCalculator(rupSet.getFaultSectionDataList());
				rupSet.addModule(distCalc);
			}
			BvalAndSegConnectivityCorrelationSampler sampler = new BvalAndSegConnectivityCorrelationSampler(
					subSects, rupSet, distCalc, MAX_DIST, ZERO_DIST_COEFF, NEG_CORR_MAX_DIST, B_VAL, SEG_MODEL, FROM_MFDS);
			return sampler;
		}
		
	}
	
	protected DefModSamplingEnabledInvConfig(int interpSkip) {
		this.interpSkip = interpSkip;
	}
	
	private Map<RupSetFaultModel, SectionCovarianceSampler> samplerCache;
	
	public final synchronized SectionCovarianceSampler getCovarianceSampler(
			FaultSystemRupSet rupSet, RupSetFaultModel fm, List<? extends FaultSection> subSects) throws IOException {
		if (samplerCache == null)
			samplerCache = new HashMap<>();
		if (samplerCache.containsKey(fm))
			return samplerCache.get(fm);
		SectionCovarianceSampler sampler = buildCovarianceSampler(rupSet, fm, subSects);
		if (cacheDir != null && cacheDir.exists() && !(sampler instanceof CachedDecomposition)) {
			// see if it's already cached
			File subDir = getCacheSubDir(fm);
			if (subDir.exists()) {
				Optional<SectionCovarianceSampler> cached = sampler.loadCached(subDir, interpSkip);
				String prefix = sampler.getFullCachePrefix(interpSkip);
				System.out.println("Checking for cached version in "+subDir.getAbsolutePath()+" with prefix: "+prefix);
				if (cached.isPresent()) {
					// use the cached version
					sampler = cached.get();
					System.out.println("Found cached version!");
				} else {
					System.out.println("Not found, will build (may be very slow)");
				}
			}
		}
		samplerCache.put(fm, sampler);
		return sampler;
	}
	
	protected File getCacheSubDir(RupSetFaultModel fm) {
		if (cacheDir == null)
			return null;
		File subDir = new File(cacheDir, "rup_sets_"+fm.getFilePrefix()+"_"+fm.getDefaultDeformationModel().getFilePrefix());
		return subDir;
	}
	
	protected abstract SectionCovarianceSampler buildCovarianceSampler(
			FaultSystemRupSet rupSet, RupSetFaultModel fm, List<? extends FaultSection> subSects) throws IOException;

	@Override
	protected List<? extends FaultSection> buildSubSectsForBranch(FaultSystemRupSet rupSet, LogicTreeBranch<?> branch) throws IOException {
		RupSetFaultModel fm = branch.requireValue(RupSetFaultModel.class);
		RupSetDeformationModel dm = branch.requireValue(RupSetDeformationModel.class);
		RandomDefModSampleNode sample = branch.requireValue(RandomDefModSampleNode.class);
		Preconditions.checkState(dm.isApplicableTo(fm),
				"Fault and deformation models are not compatible: %s, %s", fm.getName(), dm.getName());
		
		// build with original std devs
		List<? extends FaultSection> refSects;
		List<? extends FaultSection> origStdDevSects;
		synchronized (NSHM23_DeformationModels.class) {
			double origBound = NSHM23_DeformationModels.HARDCODED_FRACTIONAL_STD_DEV_UPPER_BOUND;
			double origHardcoded = NSHM23_DeformationModels.HARDCODED_FRACTIONAL_STD_DEV;
			NSHM23_DeformationModels.HARDCODED_FRACTIONAL_STD_DEV_UPPER_BOUND = 0d;
			NSHM23_DeformationModels.HARDCODED_FRACTIONAL_STD_DEV = 0d;
			origStdDevSects = super.buildSubSectsForBranch(rupSet, branch);
			NSHM23_DeformationModels.HARDCODED_FRACTIONAL_STD_DEV_UPPER_BOUND = origBound;
			NSHM23_DeformationModels.HARDCODED_FRACTIONAL_STD_DEV = origHardcoded;
			refSects = super.buildSubSectsForBranch(rupSet, branch);
		}
		
		SectionCovarianceSampler sampler = getCovarianceSampler(rupSet, fm, origStdDevSects);
		SlipRateCovarianceSampler slipSampler = new SlipRateCovarianceSampler(sampler);
		List<FaultSection> sampledSects = slipSampler.buildSample(new Well19937c(sample.getSeed()), interpSkip);
		List<FaultSection> ret = new ArrayList<>();
		for (int s=0; s<refSects.size(); s++) {
			FaultSection refSect = refSects.get(s);
			double slipRate = sampledSects.get(s).getOrigAveSlipRate();
			FaultSection modSect = refSect.clone();
			modSect.setAveSlipRate(slipRate);
			ret.add(modSect);
		}
		
		if (cacheDir != null && cacheDir.exists() && !(sampler instanceof CachedDecomposition)) {
			boolean write = true;
			try {
				int mpiRank = mpi.MPI.COMM_WORLD.Rank();
				// if we made it this far, this is an MPI job. make sure we're rank 0
				write = mpiRank == 0;
			} catch (Throwable e) {
				// will throw if MPI not on classpath, ignore
			}
			if (write) {
				File subDir = getCacheSubDir(fm);
				if (subDir.exists() || subDir.mkdir()) {
					sampler.writeCache(subDir);
				}
			}
		}
		return ret;
	}

}
