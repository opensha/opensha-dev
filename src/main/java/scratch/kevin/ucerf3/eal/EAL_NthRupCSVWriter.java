package scratch.kevin.ucerf3.eal;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupMFDsModule;
import org.opensha.sha.earthquake.param.BackgroundRupParam;
import org.opensha.sha.earthquake.param.BackgroundRupType;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sra.calc.parallel.MPJ_CondLossCalc;

import com.google.common.base.Preconditions;

import scratch.UCERF3.U3CompoundFaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.mean.TrueMeanBuilder;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;

public class EAL_NthRupCSVWriter {

	public static void main(String[] args) throws IOException {
		
		FaultSystemSolution sol = FaultSystemSolution.load(
				new File("/home/kevin/git/ucerf3-etas-launcher/inputs/"
						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_SpatSeisU3_MEAN_BRANCH_AVG_SOL.zip"));
		sol.removeModuleInstances(RupMFDsModule.class);
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.INCLUDE);
		erf.getTimeSpan().setDuration(1d);
		
		File ealBaseDir = new File("/home/kevin/OpenSHA/UCERF3/eal");
		
//		File ealSubDir = new File(ealBaseDir, "2024_05_24-ucerf3-fm31_ba-pt_gridded-ngaw2-cea-proxy-100pct-wills2015");
//		erf.setParameter(BackgroundRupParam.NAME, BackgroundRupType.POINT);
		File ealSubDir = new File(ealBaseDir, "2024_05_24-ucerf3-fm31_ba-ngaw2-cea-proxy-100pct-wills2015");
		erf.setParameter(BackgroundRupParam.NAME, BackgroundRupType.CROSSHAIR);
		
		String binFileName = "NGAWest_2014_AVG_NOIDRISS.bin";
		
		erf.updateForecast();
		
		double[] nthLosses = new double[erf.getTotNumRups()];
		int[] nthFSSIndexes = new int[nthLosses.length];
		int[] nthGridIndexes = new int[nthLosses.length];
		double[] covs = new double[nthLosses.length];
		
		for (int n=0; n<nthLosses.length; n++) {
			if (n < erf.getTotNumRupsFromFaultSystem()) {
				nthFSSIndexes[n] = erf.getFltSysRupIndexForNthRup(n);
				nthGridIndexes[n] = -1;
			} else {
				nthFSSIndexes[n] = -1;
				nthGridIndexes[n] = erf.getSrcIndexForNthRup(n) - erf.getNumFaultSystemSources();
			}
		}
		
		File binFile = new File(ealSubDir, binFileName);
		Preconditions.checkState(binFile.exists());
		
		double[][] results = MPJ_CondLossCalc.loadResults(binFile);
		
		Preconditions.checkState(results.length == erf.getNumSources());
		
		double baEAL = 0d;
		for (int sourceID=0; sourceID<erf.getNumSources(); sourceID++) {
			ProbEqkSource source = erf.getSource(sourceID);
			for (int rupID=0; rupID<source.getNumRuptures(); rupID++) {
				double rate = source.getRupture(rupID).getMeanAnnualRate(1d);
				baEAL += rate*results[sourceID][rupID];
			}
		}
		
		System.out.println("\tBranch Avg EAL:\t"+(float)baEAL);
		
		Preconditions.checkState(results.length == erf.getNumSources(),
				"Source count mismatch: file=%s, erf=%s", results.length, erf.getNumSources());
		for (int sourceID=0; sourceID<erf.getNumSources(); sourceID++) {
			int numRups = erf.getSource(sourceID).getNumRuptures();
			Preconditions.checkState(results[sourceID].length == numRups,
					"Rupture count mismatch for source %s: file=%s, erf=%s", sourceID, results.length, erf.getNumSources());
			int fssIndex, gridIndex;
			if (sourceID < erf.getNumFaultSystemSources()) {
				fssIndex = erf.getFltSysRupIndexForSource(sourceID);
				gridIndex = -1;
			} else {
				fssIndex = -1;
				gridIndex = sourceID - erf.getNumFaultSystemSources();
			}
			for (int rupID=0; rupID<results[sourceID].length; rupID++) {
				int n = erf.get_nthRupIndicesForSource(sourceID)[rupID];
				Preconditions.checkState(nthFSSIndexes[n] == fssIndex);
				Preconditions.checkState(nthGridIndexes[n] == gridIndex);
				nthLosses[n] += results[sourceID][rupID];
				covs[n] = LossCOV_Model.PORTER_POWER_LAW_2020_09_01.getCOV(nthLosses[n]);
			}
		}
		
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("Nth ERF Index", "FSS Index", "Grid Node Index", "Mag", "Conditional Average Loss", "Porter (2020) Loss COV");
		for (int n=0; n<nthLosses.length; n++)
			csv.addLine(n+"", nthFSSIndexes[n]+"", nthGridIndexes[n]+"", (float)erf.getNthRupture(n).getMag()+"", nthLosses[n]+"", covs[n]+"");
		
		csv.writeToFile(new File(ealSubDir, "average_nth_rup_losses.csv"));
	}
	
	private static void calcMeanToBAMappings(FaultSystemSolutionERF trueMeanERF, Map<U3LogicTreeBranch, List<Integer>> mappings,
			U3CompoundFaultSystemSolution cfss, FaultSystemSolutionERF baERF, FaultModels fm, int[] mappedNthIDs, double[] fractMappings) {
		Preconditions.checkState(mappedNthIDs.length == trueMeanERF.getTotNumRups());
		Preconditions.checkState(fractMappings.length == trueMeanERF.getTotNumRups());
		FaultSystemSolution baSol = baERF.getSolution();
		
		FaultSystemSolution trueMeanSol = trueMeanERF.getSolution();
		DiscretizedFunc[] rupMFDs = trueMeanSol.requireModule(RupMFDsModule.class).getRuptureMFDs();
		
		for (int i=0; i<mappedNthIDs.length; i++)
			mappedNthIDs[i] = -1;
		
		double totWeight = 0d;
		for (U3LogicTreeBranch branch : mappings.keySet()) {
			if (!branch.hasValue(fm))
				continue;
			System.out.println("Processing mappings for branch "+branch);
			List<Integer> branchIDs = mappings.get(branch);
			Preconditions.checkState(branchIDs.size() == baSol.getRupSet().getNumRuptures());
			double[] mags = cfss.getMags(branch);
			double weight = branch.getBranchWeight();
			totWeight += weight;
			
			for (int r=0; r<branchIDs.size(); r++) {
				Integer meanRupIndex = branchIDs.get(r);
				if (meanRupIndex == null || meanRupIndex < 0)
					continue;
				DiscretizedFunc mfd = rupMFDs[meanRupIndex];
				int trueMeanSourceID = trueMeanERF.getSrcIndexForFltSysRup(meanRupIndex);
				ProbEqkSource trueMeanSource = trueMeanERF.getSource(trueMeanSourceID);
				Preconditions.checkState(mfd.size() == trueMeanSource.getNumRuptures(),
						"True mean %s MFD has %s mags but source has %s rups",
						meanRupIndex, mfd.size(), trueMeanSource.getNumRuptures());
				int mfdIndex = mfd.getXIndex(mags[r]);
				if (mfdIndex < 0) {
					// try again with some extra tolerance
					for (int i=0; i<mfd.size(); i++) {
						if ((float)mags[r] == (float)mfd.getX(i) || (float)Math.abs(mags[r] - mfd.getX(i)) <= 0.01f) {
							mfdIndex = i;
							break;
						}
					}
				}
				Preconditions.checkState(mfdIndex >= 0);
				int trueMeanN = trueMeanERF.get_nthRupIndicesForSource(trueMeanSourceID)[mfdIndex];
				
				int baSourceID = baERF.getSrcIndexForFltSysRup(r);
				if (baSourceID < 0)
					// skipped rupture
					continue;
				ProbEqkSource baSource = baERF.getSource(baSourceID);
				Preconditions.checkState(baSource.getNumRuptures() == 1,
						"BA Source %s has %s ruptures?", baSourceID, baSource.getNumRuptures());
				int baN = baERF.get_nthRupIndicesForSource(baSourceID)[0];
				
				if (mappedNthIDs[trueMeanN] >= 0)
					Preconditions.checkState(mappedNthIDs[trueMeanN] == baN);
				else
					mappedNthIDs[trueMeanN] = baN;
				fractMappings[trueMeanN] += weight;
			}
		}
		
		if (totWeight != 1d)
			for (int i=0; i<fractMappings.length; i++)
				fractMappings[i] /= totWeight;
		
		// now do gridded
		int firstTrueMeanGridSourceID = trueMeanERF.getNumFaultSystemSources();
		int firstBAGridSoruceID = baERF.getNumFaultSystemSources();
		int numGridSources = trueMeanERF.getGridSourceProvider().size();
		Preconditions.checkState(numGridSources == baERF.getGridSourceProvider().size());
		for (int g=0; g<numGridSources; g++) {
			int trueMeanSourceID = firstTrueMeanGridSourceID + g;
			int baSourceID = firstBAGridSoruceID + g;
			
			ProbEqkSource trueMeanSource = trueMeanERF.getSource(trueMeanSourceID);
			ProbEqkSource baSource = baERF.getSource(baSourceID);
			Location trueMeanLoc = ((PointSurface)trueMeanSource.getSourceSurface()).getLocation();
			trueMeanLoc = new Location(trueMeanLoc.lat, trueMeanLoc.lon);
			Location baLoc = ((PointSurface)baSource.getSourceSurface()).getLocation();
			baLoc = new Location(baLoc.lat, baLoc.lon);
			Preconditions.checkState(LocationUtils.areSimilar(trueMeanLoc, baLoc),
					"Grid location mismatch: %s != %s", trueMeanLoc, baLoc);
			
			// loop over BA ruptures, which may map to multiple true mean (because of finite surfaces in true mean)
			for (int baRupID=0; baRupID<baSource.getNumRuptures(); baRupID++) {
				ProbEqkRupture rup = baSource.getRupture(baRupID);
				double mag = rup.getMag();
				double rake = rup.getAveRake();
				
				List<Integer> matches = new ArrayList<>(2);
				for (int trueMeanRupID=0; trueMeanRupID<trueMeanSource.getNumRuptures(); trueMeanRupID++) {
					ProbEqkRupture oRup = trueMeanSource.getRupture(trueMeanRupID);
					double oMag = oRup.getMag();
					double oRake = oRup.getAveRake();
					if ((float)mag == (float)oMag && (float)rake == (float)oRake)
						matches.add(trueMeanRupID);
				}
				
				Preconditions.checkState(!matches.isEmpty());
				double weightEach = 1d/matches.size();
				
				int baN = baERF.get_nthRupIndicesForSource(baSourceID)[baRupID];
				
				for (int match : matches) {
					int trueMeanN = trueMeanERF.get_nthRupIndicesForSource(trueMeanSourceID)[match];
					mappedNthIDs[trueMeanN] = baN;
					fractMappings[trueMeanN] = weightEach;
				}
			}
		}
	}

}
