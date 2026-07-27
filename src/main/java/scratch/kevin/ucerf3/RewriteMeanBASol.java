package scratch.kevin.ucerf3;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.util.modules.OpenSHA_Module;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchParentSectParticMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchRegionalMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchSectBVals;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchSectNuclMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchSectParticMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.ModSectMinMags;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupMFDsModule;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.UniqueRupture;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

import scratch.UCERF3.erf.mean.MeanUCERF3;
import scratch.UCERF3.erf.mean.MeanUCERF3.Presets;

public class RewriteMeanBASol {

	public static void main(String[] args) throws IOException {
		// takes the mean UCERF3 solution in the old format and with all branch-specific mod-sect-min-mag imprinted,
		// and maps it to the modular re-built solution
		// so we'll have a modern solution file that matches the original branch-averaged hazard (including TD RI calcs)
		
		File baseDir = new File("/home/kevin/OpenSHA/fss_inversions/2021_11_30-u3_branches-orig_calcs-5h");
		
//		File modInSolFile = new File(baseDir, "FM3_1_branch_averaged_full_modules.zip");
//		Presets preset = Presets.FM3_1_BRANCH_AVG;
//		String outPrefix = "FM3_1_branch_averaged_modular_from_erf";
		File modInSolFile = new File(baseDir, "FM3_2_branch_averaged_full_modules.zip");
		Presets preset = Presets.FM3_2_BRANCH_AVG;
		String outPrefix = "FM3_2_branch_averaged_modular_from_erf";
		
		File outSolFile = new File(baseDir, outPrefix+".zip");
		
		FaultSystemSolution modInSol = FaultSystemSolution.load(modInSolFile);
		FaultSystemRupSet modInRS = modInSol.getRupSet();
		
		// remove logic tree modules
		modInSol.removeModuleInstances(BranchSectParticMFDs.class);
		modInSol.removeModuleInstances(BranchSectNuclMFDs.class);
		modInSol.removeModuleInstances(BranchParentSectParticMFDs.class);
		modInSol.removeModuleInstances(BranchSectBVals.class);
		modInSol.removeModuleInstances(BranchRegionalMFDs.class);
		
		// get any prints out the way
		modInRS.getModules(true);
		modInSol.getModules(true);
		
		MeanUCERF3 erf = new MeanUCERF3();
		erf.setPreset(preset);
		erf.updateForecast();
		FaultSystemSolution erfSol = erf.getSolution();
		FaultSystemRupSet erfRS = erfSol.getRupSet();
		
		// get any prints out the way
		erfRS.getModules(true);
		erfSol.getModules(true);
		
		Map<String, List<FaultSection>> inSectsByParentName = modInRS.getFaultSectionDataList().stream().collect(
				Collectors.groupingBy(S->S.getParentSectionName()));
		
		Preconditions.checkState(erfRS.getNumSections() == modInRS.getNumSections());
		int[] sectsERFtoMod = new int[erfRS.getNumSections()];
		int[] sectsModToERF = new int[erfRS.getNumSections()];
		for (FaultSection erfSect : erfRS.getFaultSectionDataList()) {
			String erfParentName = erfSect.getParentSectionName();
			Preconditions.checkState(inSectsByParentName.containsKey(erfParentName),
					"ERF parent not found in mod: %s", erfParentName);
			List<FaultSection> candidateSects = inSectsByParentName.get(erfParentName);
			String erfSectName = erfSect.getName();
			FaultSection match = null;
			for (FaultSection candidate : candidateSects) {
				String candidateName = candidate.getName();
				if (candidateName.equals(erfSectName)) {
					Preconditions.checkState(match == null,
							"Multiple sections match name: %s\n\tFirst: %s\n\tSubsequent: %s",
							erfSectName, match, candidate);
					match = candidate;
				}
			}
			Preconditions.checkNotNull(match, "No match found for sect: %s", erfSectName);
			sectsERFtoMod[erfSect.getSectionId()] = match.getSectionId();
			sectsModToERF[match.getSectionId()] = erfSect.getSectionId();
		}
		
		double[] updatedMags = new double[modInRS.getNumRuptures()];
		double[] updatedRakes = new double[modInRS.getNumRuptures()];
		// initialize mags and rakes to original
		// and build sections-to-rup map
		Map<UniqueRupture, Integer> modUniques = new HashMap<>();
		for (int r=0; r<modInRS.getNumRuptures(); r++) {
			updatedMags[r] = modInRS.getMagForRup(r);
			updatedRakes[r] = modInRS.getAveRakeForRup(r);
			modUniques.put(UniqueRupture.forIDs(modInRS.getSectionsIndicesForRup(r)), r);
		}
		double[] updatedRates = new double[modInRS.getNumRuptures()];
		boolean[] rupMapped = new boolean[modInRS.getNumRuptures()];
		int numRupsMapped = 0;
		DiscretizedFunc[] updatedMFDs = new DiscretizedFunc[modInRS.getNumRuptures()];
		RupMFDsModule erfMFDs = erfSol.requireModule(RupMFDsModule.class);
		System.out.println("***************************************************");
		System.out.println("Doing rate mapping");
		for (int erfRupIndex=0; erfRupIndex<erfRS.getNumRuptures(); erfRupIndex++) {
			List<Integer> erfIndexes = erfRS.getSectionsIndicesForRup(erfRupIndex);
			List<Integer> modIndexes = new ArrayList<>(erfIndexes.size());
			for (int erfIndex : erfIndexes)
				modIndexes.add(sectsERFtoMod[erfIndex]);
			UniqueRupture unique = UniqueRupture.forIDs(modIndexes);
			int modRupIndex = modUniques.get(unique);
			Preconditions.checkState(!rupMapped[modRupIndex]);
			rupMapped[modRupIndex] = true;
			numRupsMapped++;
			updatedRates[modRupIndex] = erfSol.getRateForRup(erfRupIndex);
			updatedMags[modRupIndex] = erfRS.getMagForRup(erfRupIndex);
			updatedRakes[modRupIndex] = erfRS.getAveRakeForRup(erfRupIndex);
			updatedMFDs[modRupIndex] = erfMFDs.getRuptureMFD(erfRupIndex);
			if (updatedMFDs[modRupIndex] != null)
				Preconditions.checkState((float)updatedRates[modRupIndex] == (float)updatedMFDs[modRupIndex].calcSumOfY_Vals());
		}
		System.out.println("***************************************************");
		
		FaultSystemRupSet updatedRS = FaultSystemRupSet.buildFromExisting(modInRS, false)
				.rupMags(updatedMags)
				.rupRakes(updatedRakes)
				.build();
		
		// update aseismicity
		int numAseisUpdated = 0;
		for (int s=0; s<updatedRS.getNumSections(); s++) {
			int erfSectID = sectsModToERF[s];
			FaultSection modSect = updatedRS.getFaultSectionData(s);
			FaultSection erfSect = erfRS.getFaultSectionData(erfSectID);
			double modAseis = modSect.getAseismicSlipFactor();
			double erfAseis = erfSect.getAseismicSlipFactor();
			if ((float)modAseis != (float)erfAseis) {
				numAseisUpdated++;
				modSect.setAseismicSlipFactor(erfAseis);
			}
		}
		
		List<String> copiedRSModules = new ArrayList<>();
		for (OpenSHA_Module module : modInRS.getModules(true)) {
			if (!(module instanceof ModSectMinMags)) {
				updatedRS.addModule(module);
				copiedRSModules.add(module.getName());
			}
		}
		
		FaultSystemSolution updatedSol = new FaultSystemSolution(updatedRS, updatedRates);
		List<String> copiedSolModules = new ArrayList<>();
		for (OpenSHA_Module module : modInSol.getModules(true)) {
			if (!(module instanceof RupMFDsModule)) {
				updatedSol.addModule(module);
				copiedSolModules.add(module.getName());
			}
		}
		updatedSol.addModule(new RupMFDsModule(updatedSol, updatedMFDs));
		
		CSVFile<String> sectsCSV = new CSVFile<>(true);
		sectsCSV.addLine("Subsection Index", "Subsection Name",
				"ERF", "Partic Rate", "Partic RI", "Min BA Mag (rate>0)", "Min MFD Mag",
				"BA Sol", "Partic Rate", "Partic RI", "Min BA Mag (rate>0)", "Min MFD Mag",
				"BA Sol w/ ModMinMags applied", "Partic Rate", "Partic RI", "Min BA Mag (rate>0)", "Min MFD Mag");
		
		double[] modInParticRates = modInSol.calcTotParticRateForAllSects();
		double[] modInParticRatesModMin = calcModSectMinMagSectParticRates(modInSol);
		double[] updatedParticRates = updatedSol.calcTotParticRateForAllSects();
		ModSectMinMags inMinMags = modInRS.requireModule(ModSectMinMags.class);
		RupMFDsModule inMFDs = modInSol.requireModule(RupMFDsModule.class);
		
		System.out.println("***************************************************");
		int numSectRatesChanged = 0;
		int numIncreased = 0;
		int numDecreased = 0;
		for (int s=0; s<modInRS.getNumSections(); s++) {
			double inRate = modInParticRates[s];
			double inRateModMin = modInParticRatesModMin[s];
			double updatedRate = updatedParticRates[s];
			if ((float)inRate != (float)updatedRate || (float)inRateModMin != (float)updatedRate) {
				System.out.println(s+". "+modInRS.getFaultSectionData(s).getName()+":");
				System.out.println("\tFull rate:\t"+(float)inRate+" -> "+(float)updatedRate
						+"\t("+(float)(1d/inRate)+" -> "+(float)(1/updatedRate)+")");
				System.out.println("\tModMin rate:\t"+(float)inRateModMin+" -> "+(float)updatedRate
						+"\t("+(float)(1d/inRateModMin)+" -> "+(float)(1/updatedRate)+")");
				numSectRatesChanged++;
				if (updatedRate > inRate)
					numIncreased++;
				else
					numDecreased++;
				List<String> line = new ArrayList<>(sectsCSV.getNumCols());
				line.add(s+"");
				line.add(modInRS.getFaultSectionData(s).getSectionName());
				for (int i=0; i<3; i++) {
					line.add("");
					double particRate = switch (i) {
					case 0:
						yield updatedRate;
					case 1:
						yield inRate;
					case 2:
						yield inRateModMin;
					default:
						throw new IllegalArgumentException("Unexpected value: " + i);
					};
					line.add((float)particRate+"");
					line.add((float)(1d/particRate)+"");
					double minMag = Double.POSITIVE_INFINITY;
					double minMFDMag = Double.POSITIVE_INFINITY;
					for (int r : modInRS.getRupturesForSection(s)) {
						if (i == 3 && inMinMags.isRupBelowSectMinMag(r))
							continue;
						double rate;
						DiscretizedFunc mfd;
						double mag;
						if (i == 0) {
							rate = updatedSol.getRateForRup(r);
							mfd = updatedMFDs[r];
							mag = updatedRS.getMagForRup(r);
						} else {
							rate = modInSol.getRateForRup(r);
							mfd = inMFDs.getRuptureMFD(r);
							mag = modInRS.getMagForRup(r);
						}
						if (rate > 0d) {
							minMag = Math.min(minMag, mag);
							if (mfd != null) {
								for (Point2D pt : mfd) {
									if (pt.getY() > 0d) {
										minMFDMag = Math.min(minMFDMag, pt.getX());
										break;
									}
								}
							}
						}
					}
					line.add((float)minMag+"");
					line.add((float)minMFDMag+"");
				}
				sectsCSV.addLine(line);
			}
		}
		System.out.println("Modified participation rates on "+numSectRatesChanged+"/"+updatedRS.getNumSections()+" sections");
		System.out.println("\t"+numIncreased+" increased, "+numDecreased+" decreased");
		sectsCSV.writeToFile(new File(baseDir, outPrefix+"_sects.csv"));
		
		System.out.println("***************************************************");
		System.out.println("Mapped "+numRupsMapped+" ruptures");
		System.out.println("Updated aseismicity on "+numAseisUpdated+"/"+updatedRS.getNumSections()+" sections");
		System.out.println("Copied these rupture set modules:");
		for (String name : copiedRSModules)
			System.out.println("\t"+name);
		System.out.println("Copied these solution modules:");
		for (String name : copiedSolModules)
			System.out.println("\t"+name);
		
		System.out.println("Mod sol has:");
		printSolStats(modInSol);
		
		System.out.println("ERF sol has:");
		printSolStats(erfSol);
		
		System.out.println("Updated sol has:");
		printSolStats(updatedSol);
		System.out.println("***************************************************");
		
		updatedSol.write(outSolFile);
	}
	
	private static void printSolStats(FaultSystemSolution sol) {
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		System.out.println("\tRS class: "+rupSet.getClass());
		System.out.println("\tSol class: "+sol.getClass());
		System.out.println("\t"+rupSet.getNumSections()+" sections");
		System.out.println("\t"+rupSet.getNumRuptures()+" ruptures");
		System.out.println("\t"+sol.getTotalRateForAllFaultSystemRups()+" total rate");
		RupMFDsModule mfds = sol.requireModule(RupMFDsModule.class);
		double sumMFDrate = 0d;
		for (DiscretizedFunc mfd : mfds.getRuptureMFDs())
			if (mfd != null)
				sumMFDrate += mfd.calcSumOfY_Vals();
		System.out.println("\t"+sumMFDrate+" total RupMFD rate");
		ModSectMinMags modMinMags = rupSet.getModule(ModSectMinMags.class);
		System.out.println("\tModSectMinMags:\t"+modMinMags);
		if (modMinMags != null) {
			double modRate = 0d;
			for (int r=0; r<rupSet.getNumRuptures(); r++)
				if (!modMinMags.isRupBelowSectMinMag(r))
					modRate += sol.getRateForRup(r);
			System.out.println("\t"+(float)modRate+" total rate after appling ModSectMinMags");
		}
	}
	
	private static double[] calcModSectMinMagSectParticRates(FaultSystemSolution sol) {
		FaultSystemRupSet rupSet = sol.getRupSet();
		ModSectMinMags minMags = rupSet.getModule(ModSectMinMags.class);
		if (minMags == null)
			return sol.calcTotParticRateForAllSects();
		double[] ret = new double[rupSet.getNumSections()];
		for (int r=0; r<rupSet.getNumRuptures(); r++) {
			if (!minMags.isRupBelowSectMinMag(r)) {
				double rate = sol.getRateForRup(r);
				for (int s : rupSet.getSectionsIndicesForRup(r))
					ret[s] += rate;
			}
		}
		return ret;
	}

}
