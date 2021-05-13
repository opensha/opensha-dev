package scratch.kevin.ucerf3.eal;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;

import org.dom4j.DocumentException;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sra.calc.parallel.MPJ_CondLossCalc;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.erf.mean.TrueMeanBuilder;
import scratch.UCERF3.griddedSeismicity.AbstractGridSourceProvider;
import scratch.UCERF3.griddedSeismicity.GridSourceProvider;
import scratch.UCERF3.logicTree.APrioriBranchWeightProvider;
import scratch.UCERF3.logicTree.BranchWeightProvider;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.utils.FaultSystemIO;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;

public class UCERF3_BranchAvgLossFetcher {
	
	private File[] ealDataDirs;
	
	private FaultSystemSolution trueMeanSol;
	private Map<LogicTreeBranch, List<Integer>> branchMappings;
	
	private CompoundFaultSystemSolution cfss;
	
	/**
	 * @param trueMeanSolFile True mean fault system solution, this must be the exact same file used with
	 * MPJ_CondLossCalc to calculate losses, and must include rupture mappings in the zip file.
	 * @param cfss Compound fault system solution used to retrieve branch specific mags/rates
	 * @param ealDataDirs Directory/directories that contain(s) EAL result binary files.
	 * @throws IOException
	 * @throws DocumentException
	 */
	public UCERF3_BranchAvgLossFetcher(File trueMeanSolFile, CompoundFaultSystemSolution cfss, File... ealDataDirs)
			throws IOException, DocumentException {
		this(FaultSystemIO.loadSol(trueMeanSolFile), TrueMeanBuilder.loadRuptureMappings(trueMeanSolFile), cfss, ealDataDirs);
	}
	
	public UCERF3_BranchAvgLossFetcher(FaultSystemSolution trueMeanSol, Map<LogicTreeBranch, List<Integer>> branchMappings,
			CompoundFaultSystemSolution cfss, File... ealDataDirs) throws IOException, DocumentException {
		Preconditions.checkState(ealDataDirs.length > 0, "Must supply at least one data dir");
		this.ealDataDirs = ealDataDirs;
		this.trueMeanSol = trueMeanSol;
		this.branchMappings = branchMappings;
		this.cfss = cfss;
	}
	
	/**
	 * This returns a loss distribution for each rupture in the given fault model for the specified IMR. It will first
	 * search through the previously supplied data directories looking for a file containing losses for that IMR.<br>
	 * 
	 * A distribution is returned for each rupture rather than a single value as different logic tree branches have
	 * different magnitudes/rakes/areas for each rupture which affect the conditional loss. The distributions can
	 * optionally be normalized (thus making them conditional losses).
	 * @param attenRelRef desired attenuation relationship - note that this must already be calculated and stored in
	 * one of the data directories.
	 * @param fm fault model for which to calculate results
	 * @param normalize if true losses will be normalized to conditional losses, otherwise the y values will sum to the total
	 * time independent annual rate (adjusted for branches with ruptures below section minimum mag).
	 * @return distribution of loss (y values are weights, optionally normalized) for each rupture in the given FM.
	 * @throws IOException if there is an error loading the data files
	 * @throws IllegalStateException if no file is found for the given IMR.
	 */
	public DiscretizedFunc[] getFaultLosses(AttenRelRef attenRelRef, FaultModels fm, boolean normalize) throws IOException {
		File rupLossesFile = getFaultLossFile(attenRelRef, false);
		
		System.out.println("Loading losses from: "+rupLossesFile.getAbsolutePath());
		double[][] rupLosses = MPJ_CondLossCalc.loadResults(rupLossesFile);
		DiscretizedFunc[] totRupMFDs = trueMeanSol.getRupMagDists();
		
		// calc true mean eal as a check
		double trueMeanEAL = 0d;
		for (int r=0; r<rupLosses.length; r++)
			for (int m=0; m<totRupMFDs[r].size(); m++)
				trueMeanEAL += rupLosses[r][m] * totRupMFDs[r].getY(m);
		System.out.println("True mean fault based EAL (both FMs): "+trueMeanEAL);
		
		BranchWeightProvider weightProv = new APrioriBranchWeightProvider();
		
		double totWeight = 0d;
		
		DiscretizedFunc[] rupLossDists = null;
		for (LogicTreeBranch branch : branchMappings.keySet()) {
			FaultModels branchFM = branch.getValue(FaultModels.class);
			if (branchFM != fm)
				continue;
			
			List<Integer> branchMapping = branchMappings.get(branch);
			Preconditions.checkNotNull(branchMapping);
			
			double[] rates = cfss.getRates(branch);
			double[] mags = cfss.getMags(branch);
			
			double weight = weightProv.getWeight(branch);
			totWeight += weight;
			
			if (rupLossDists == null) {
				rupLossDists = new DiscretizedFunc[rates.length];
				for (int i=0; i<rupLossDists.length; i++)
					rupLossDists[i] = new ArbitrarilyDiscretizedFunc();
			} else {
				Preconditions.checkState(rupLossDists.length == rates.length);
			}
			
			for (int r=0; r<rates.length; r++) {
				double rate = rates[r];
				double mag = mags[r];
				
				int trueMeanIndex = branchMapping.get(r);
				if (rate == 0 || trueMeanIndex < 0)
					// zero rate or below sub seismo
					continue;
				
				// MFD in true mean sol
				DiscretizedFunc mfd = totRupMFDs[trueMeanIndex];
				Preconditions.checkState(mfd.size() == rupLosses[trueMeanIndex].length);
				
				int lossIndex = getMatchingXIndexFloatPrecision(mag, mfd);
				Preconditions.checkState(lossIndex >= 0);
				double loss = rupLosses[trueMeanIndex][lossIndex];
				if (loss == 0d)
					continue;
				
				DiscretizedFunc rupDist = rupLossDists[r];
				int indexInLossDist = getMatchingXIndexFloatPrecision(loss, rupDist);
				if (indexInLossDist < 0)
					// first time we have this loss
					rupDist.set(loss, weight*rate);
				else
					rupDist.set(loss, rupDist.getY(indexInLossDist)+weight*rate);
			}
		}
		
		if (normalize) {
			for (DiscretizedFunc rupDist : rupLossDists)
				rupDist.scale(1d/((ArbitrarilyDiscretizedFunc)rupDist).calcSumOfY_Vals());
		} else if (totWeight != 1d) {
			for (DiscretizedFunc rupDist : rupLossDists)
				rupDist.scale(1d/totWeight);
		}
		
		return rupLossDists;
	}
	
	private File getFaultLossFile(AttenRelRef attenRelRef, boolean gridded) {
		String fName;
		if (gridded)
			fName = attenRelRef.name()+"_fss_gridded.bin";
		else
			fName = attenRelRef.name()+"_fss_index.bin";
		for (File dir : ealDataDirs) {
			File file = new File(dir, fName);
			if (file.exists())
				return file;
		}
		throw new IllegalStateException("File '"+fName+"' not found in data directories");
	}
	
	/**
	 * This returns a magnitude loss distribution for gridded seismicity at every node in the given region. It will first
	 * search through the previously supplied data directories looking for a file containing losses for that IMR.<br>
	 * 
	 * Note that magnitudes may differ in double precision with magnitudes from your grid source provider and should be
	 * mapped in floating precision.
	 * 
	 * @param attenRelRef desired attenuation relationship - note that this must already be calculated and stored in
	 * one of the data directories.
	 * @param region - can be null, only used for indexing verification
	 * @return
	 * @throws IOException
	 */
	public DiscretizedFunc[] getGriddedMagLossDists(AttenRelRef attenRelRef, GriddedRegion region) throws IOException {
		File griddedLossesFile = getFaultLossFile(attenRelRef, true);
		
		System.out.println("Loading losses from: "+griddedLossesFile.getAbsolutePath());
		return MPJ_CondLossCalc.loadGridSourcesFile(griddedLossesFile, region);
	}
	
	public static int getMatchingXIndexFloatPrecision(double x, DiscretizedFunc func) {
		for (int i=0; i<func.size(); i++)
			if ((float)func.getX(i) == (float)x)
				return i;
		return -1;
	}
	
	public static void main(String[] args) throws IOException, DocumentException {
		Stopwatch watch = Stopwatch.createStarted();
		// true mean FSS which includes rupture mapping information. this must be the exact file used to calulate EALs
		File trueMeanSolFile = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_TRUE_HAZARD_MEAN_SOL_WITH_MAPPING.zip");
		
		// directory which contains EAL data
		File dataDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2014_05_05-ucerf3-eal-calc-wald-vs30");
		
		// IMR for which EAL data has already been computed
		AttenRelRef attenRelRef = AttenRelRef.ASK_2014;
		
		// Fault model of interest
		FaultModels fm = FaultModels.FM3_1;
		
		// Branch averaged FSS
		FaultSystemSolution baSol = FaultSystemIO.loadSol(new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		
		// Compound fault system solution
		CompoundFaultSystemSolution cfss = CompoundFaultSystemSolution.fromZipFile(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip"));
		
		UCERF3_BranchAvgLossFetcher fetch = new UCERF3_BranchAvgLossFetcher(trueMeanSolFile, cfss, dataDir);
		
		watch.stop();
		System.out.println("Initial setup took "+watch.elapsed(TimeUnit.SECONDS)+" s");
		watch.reset();
		watch.start();
		DiscretizedFunc[] faultLosses = fetch.getFaultLosses(attenRelRef, fm, false);
		watch.stop();
		System.out.println("Fault losses took "+watch.elapsed(TimeUnit.SECONDS)+" s");
		watch.reset();
		watch.start();
		// calculate total fault based EAL
		double totFaultLosses = 0;
		int numRateDiscrep = 0;
		for (int rupID=0; rupID<faultLosses.length; rupID++) {
			DiscretizedFunc faultLossDist = faultLosses[rupID];
			double totRate = 0d;
			for (Point2D pt : faultLossDist) {
				totFaultLosses += pt.getX()*pt.getY();
				totRate += pt.getY();
			}
			if ((float)totRate != (float)baSol.getRateForRup(rupID))
				numRateDiscrep++;
		}
		System.out.println("Total Fault Losses: "+totFaultLosses);
		System.out.println(numRateDiscrep+"/"+faultLosses.length+" ("
				+(float)(100d*numRateDiscrep/(double)faultLosses.length)+" %) rate discrepancies. "
						+ "Should be low, only for ruptures below section minimum magnitude which isn't accounted for in BA sol");
		
		// grid sources
		GridSourceProvider prov = baSol.getGridSourceProvider();
		// this conditional loss as a function of magnitude at each grid node 
		DiscretizedFunc[] griddedDists = fetch.getGriddedMagLossDists(attenRelRef, prov.getGriddedRegion());
		watch.stop();
		System.out.println("Gridded losses took "+watch.elapsed(TimeUnit.SECONDS)+" s");
		watch.reset();
		double totGriddedLosses = 0;
		for (int i=0; i<griddedDists.length; i++) {
			IncrementalMagFreqDist mfd = prov.getNodeMFD(i, AbstractGridSourceProvider.SOURCE_MIN_MAG_CUTOFF);
			for (Point2D pt : mfd) {
				double mag = pt.getX();
				double rate = pt.getY();
				if (rate == 0d)
					continue;
				// must get loss with this method as mags in loss distribution may not match our mag to double precision
				int lossIndex = getMatchingXIndexFloatPrecision(mag, griddedDists[i]);
				Preconditions.checkState(lossIndex >= 0, "Loss not found for node "+i+", mag="+mag);
				double loss = griddedDists[i].getY(lossIndex);
				
				totGriddedLosses += loss*rate;
			}
		}
		System.out.println("Total Gridded Losses: "+totGriddedLosses);
		System.out.println("Total Losses: "+(totGriddedLosses+totFaultLosses));
	}

}
