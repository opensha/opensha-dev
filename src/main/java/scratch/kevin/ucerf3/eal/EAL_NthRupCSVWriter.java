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
		// these loss results use a true mean solution
		File trueMeanSolFile = new File("/home/kevin/OpenSHA/UCERF3/rup_sets/orig/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_TRUE_HAZARD_MEAN_SOL_WITH_MAPPING.zip");
		File compoundSolFile = new File("/home/kevin/OpenSHA/UCERF3/rup_sets/orig/full_ucerf3_compound_sol.zip");
		FaultSystemSolution trueMeanSol = FaultSystemSolution.load(trueMeanSolFile);
		Map<U3LogicTreeBranch, List<Integer>> mappings = TrueMeanBuilder.loadRuptureMappings(trueMeanSolFile);
		U3CompoundFaultSystemSolution cfss = U3CompoundFaultSystemSolution.fromZipFile(compoundSolFile);
		FaultSystemSolutionERF trueMeanERF = new FaultSystemSolutionERF(trueMeanSol);
		trueMeanERF.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.INCLUDE);
		trueMeanERF.setParameter(BackgroundRupParam.NAME, BackgroundRupType.CROSSHAIR);
		trueMeanERF.getTimeSpan().setDuration(1d);
		trueMeanERF.updateForecast();
		
		// we'll use this to calculate the FM3.1 only true mean EAL for validation
		BitSet uniqueToFM31 = new BitSet(trueMeanSol.getRupSet().getNumRuptures());
		// set all  FM3.1
		for (U3LogicTreeBranch branch : mappings.keySet()) {
			if (!branch.hasValue(FaultModels.FM3_1))
				continue;
			for (Integer rup : mappings.get(branch)) {
				if (rup != null && rup >= 0)
					uniqueToFM31.set(rup);
			}
		}
		// unset all FM3.2
		for (U3LogicTreeBranch branch : mappings.keySet()) {
			if (branch.hasValue(FaultModels.FM3_1))
				continue;
			for (Integer rup : mappings.get(branch)) {
				if (rup != null && rup >= 0)
					uniqueToFM31.clear(rup);
			}
		}
		
		FaultSystemSolution baSol = FaultSystemSolution.load(
				new File("/home/kevin/OpenSHA/UCERF3/rup_sets/modular/FM3_1_SpatSeisU3_branch_averaged.zip"));
		baSol.removeModuleInstances(RupMFDsModule.class);
		FaultSystemSolutionERF baERF = new FaultSystemSolutionERF(baSol);
		baERF.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.INCLUDE);
		baERF.setParameter(BackgroundRupParam.NAME, BackgroundRupType.POINT);
		baERF.getTimeSpan().setDuration(1d);
		baERF.updateForecast();
		
		int[] mappedNthIDs = new int[trueMeanERF.getTotNumRups()];
		double[] fractMappings = new double[trueMeanERF.getTotNumRups()];
		calcMeanToBAMappings(trueMeanERF, mappings, cfss, baERF, FaultModels.FM3_1, mappedNthIDs, fractMappings);
		
		AttenRelRef[] gmpeRefs = {
				AttenRelRef.CB_2014,
				AttenRelRef.ASK_2014,
				AttenRelRef.BSSA_2014,
				AttenRelRef.CY_2014,
				AttenRelRef.IDRISS_2014
		};
		double[] gmpeWeights = {
				0.22,
				0.22,
				0.22,
				0.22,
				0.12
		};
		
		Preconditions.checkState(gmpeWeights.length == gmpeRefs.length);
		
		File ealBaseDir = new File("/home/kevin/OpenSHA/UCERF3/eal");
		
		File[] ealSubDirs = {
				new File(ealBaseDir, "2020_03_17-ucerf3-ngaw2-cea-proxy-100pct-wald2007"),
				new File(ealBaseDir, "2020_03_17-ucerf3-ngaw2-cea-proxy-100pct-wills2015")
		};
		double[] ealWeights = {
				0.5d,
				0.5d
		};
		Preconditions.checkState(ealWeights.length == ealSubDirs.length);
		
		File outputDir = new File(ealBaseDir, "2020_03_27-ucerf3-ngaw2-cea-100pct-consolidate");
		
		double[] nthLosses = new double[baERF.getTotNumRups()];
		int[] nthFSSIndexes = new int[nthLosses.length];
		int[] nthGridIndexes = new int[nthLosses.length];
		
		for (int n=0; n<nthLosses.length; n++) {
			if (n < baERF.getTotNumRupsFromFaultSystem()) {
				nthFSSIndexes[n] = baERF.getFltSysRupIndexForNthRup(n);
				nthGridIndexes[n] = -1;
			} else {
				nthFSSIndexes[n] = -1;
				nthGridIndexes[n] = baERF.getSrcIndexForNthRup(n) - baERF.getNumFaultSystemSources();
			}
		}
		
		double totWeight = 0d;
		
		for (int g=0; g<gmpeRefs.length; g++) {
			String binFileName = gmpeRefs[g].name()+".bin";
			for (int e=0; e<ealSubDirs.length; e++) {
				File binFile = new File(ealSubDirs[e], binFileName);
				Preconditions.checkState(binFile.exists());
				double weight = gmpeWeights[g] * ealWeights[e];
				System.out.println("Processing "+binFile.getAbsolutePath()+" with weight="+(float)weight);
				totWeight += weight;
				
				double[][] trueMeanResults = MPJ_CondLossCalc.loadResults(binFile);
				
				System.out.println("Mapping True Mean results to BA");
				// now map from true mean results to ba results
				double[][] results = new double[baERF.getNumSources()][];
				for (int sourceID=0; sourceID<baERF.getNumSources(); sourceID++)
					results[sourceID] = new double[baERF.getSource(sourceID).getNumRuptures()];
				
				double trueMeanEAL = 0d;
				double trueMeanEstFM31_EAL = 0d;
				for (int trueMeanN=0; trueMeanN<mappedNthIDs.length; trueMeanN++) {
					int trueMeanSourceID = trueMeanERF.getSrcIndexForNthRup(trueMeanN);
					int trueMeanRupID = trueMeanERF.getRupIndexInSourceForNthRup(trueMeanN);

					double trueMeanRate = trueMeanERF.getRupture(trueMeanSourceID, trueMeanRupID).getMeanAnnualRate(1d);
					trueMeanEAL += trueMeanResults[trueMeanSourceID][trueMeanRupID] * trueMeanRate;
					
					if (mappedNthIDs[trueMeanN] < 0)
						// not mapped, probably from another fault model
						continue;

					if (trueMeanSourceID < trueMeanERF.getNumFaultSystemSources() && uniqueToFM31.get(trueMeanERF.getFltSysRupIndexForNthRup(trueMeanN)))
						// this ruptures was only on FM3.1, so it's rate was cut in half; double it
						trueMeanRate *= 2;
					trueMeanEstFM31_EAL += trueMeanResults[trueMeanSourceID][trueMeanRupID] * trueMeanRate;
					
					int baSourceID = baERF.getSrcIndexForNthRup(mappedNthIDs[trueMeanN]);
					int baRupID = baERF.getRupIndexInSourceForNthRup(mappedNthIDs[trueMeanN]);
					
					results[baSourceID][baRupID] += fractMappings[trueMeanN] * trueMeanResults[trueMeanSourceID][trueMeanRupID];
				}
				
				double baEAL = 0d;
				for (int sourceID=0; sourceID<baERF.getNumSources(); sourceID++) {
					ProbEqkSource source = baERF.getSource(sourceID);
					for (int rupID=0; rupID<source.getNumRuptures(); rupID++) {
						double rate = source.getRupture(rupID).getMeanAnnualRate(1d);
						baEAL += rate*results[sourceID][rupID];
					}
				}

				System.out.println("\tTrue mean EAL:\t"+(float)trueMeanEAL);
				System.out.println("\tTrue mean est FM3.1 EAL:\t"+(float)trueMeanEstFM31_EAL);
				System.out.println("\tBranch Avg EAL:\t"+(float)baEAL);
				
				Preconditions.checkState(results.length == baERF.getNumSources(),
						"Source count mismatch: file=%s, erf=%s", results.length, baERF.getNumSources());
				for (int sourceID=0; sourceID<baERF.getNumSources(); sourceID++) {
					int numRups = baERF.getSource(sourceID).getNumRuptures();
					Preconditions.checkState(results[sourceID].length == numRups,
							"Rupture count mismatch for source %s: file=%s, erf=%s", sourceID, results.length, baERF.getNumSources());
					int fssIndex, gridIndex;
					if (sourceID < baERF.getNumFaultSystemSources()) {
						fssIndex = baERF.getFltSysRupIndexForSource(sourceID);
						gridIndex = -1;
					} else {
						fssIndex = -1;
						gridIndex = sourceID - baERF.getNumFaultSystemSources();
					}
					for (int rupID=0; rupID<results[sourceID].length; rupID++) {
						int n = baERF.get_nthRupIndicesForSource(sourceID)[rupID];
						Preconditions.checkState(nthFSSIndexes[n] == fssIndex);
						Preconditions.checkState(nthGridIndexes[n] == gridIndex);
						nthLosses[n] += results[sourceID][rupID]*weight;
					}
				}
			}
		}
		
		if ((float)totWeight != 1f)
			for (int n=0; n<nthLosses.length; n++)
				nthLosses[n] /= totWeight;
		
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("Nth ERF Index", "FSS Index", "Grid Node Index", "Mag", "Conditional Average Loss");
		for (int n=0; n<nthLosses.length; n++)
			csv.addLine(n+"", nthFSSIndexes[n]+"", nthGridIndexes[n]+"", (float)baERF.getNthRupture(n).getMag()+"", nthLosses[n]+"");
		
		csv.writeToFile(new File(outputDir, "average_nth_rup_losses.csv"));
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
