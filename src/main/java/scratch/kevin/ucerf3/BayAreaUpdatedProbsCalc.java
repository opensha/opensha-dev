package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.param.AleatoryMagAreaStdDevParam;
import org.opensha.sha.earthquake.param.HistoricOpenIntervalParam;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityOptions;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

import cern.colt.Arrays;
import scratch.UCERF3.analysis.FaultSysSolutionERF_Calc;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.mean.MeanUCERF3;
import scratch.UCERF3.erf.mean.MeanUCERF3.Presets;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;

public class BayAreaUpdatedProbsCalc {
	
	static Map<String, List<Integer>> loadFaultParentSectMappings() {
		Map<String, List<Integer>> map = new HashMap<>();
		
		map.put("San Andreas", List.of(654, 655, 657, 658));
		map.put("Hayward-Rodgers Creek", List.of(637, 638, 639, 651));
		map.put("Calaveras", List.of(601, 602, 603, 621));
		map.put("Concord", List.of(635, 636, 711, 713, 2, 622, 623, 640, 677));
		map.put("San Gregorio", List.of(660, 661));
		map.put("Maacama", List.of(644));
		// these aren't bay area, but from Ruth's 8/15/24 e-mail
//		map.put("SAF Creeping", List.of(658));
//		map.put("SAF Parkfield", List.of(32));
//		map.put("Imperial", List.of(97));
//		map.put("Superstition Hills", List.of(98));
//		map.put("S. SAF", List.of(285, 300, 287, 286, 301, 282, 283, 284, 295));
		
		return map;
	}

	public static void main(String[] args) throws IOException {
		int year = 2026;
//		int year = 2014;
		int duration = 30;
		
//		File outputDir = new File("/home/kevin/OpenSHA/UCERF3/bay_area_updated_probs");
		File outputDir = new File("/tmp/ruth_probs");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		double[] minMags = {6.7, 7, 7.5};
		
		Presets[] presets = {Presets.FM3_1_MAG_VAR, Presets.FM3_2_MAG_VAR};
//		Presets[] presets = {Presets.FM3_1_BRANCH_AVG, Presets.FM3_2_BRANCH_AVG};
		double[] presetWeights = {0.5d, 0.5d};
		
		MeanUCERF3 erf = new MeanUCERF3();
		erf.setPreset(presets[0]);
		
		String regName = "Bay Area Region";
		Region reg = new CaliforniaRegions.SF_BOX();
		boolean clipFaultsToReg = true;
//		String regName = null;
//		Region reg = null;
//		boolean clipFaultsToReg = false;
		Map<String, List<Integer>> faultIDs = loadFaultParentSectMappings();
		
		List<String> faultsSorted = new ArrayList<>();
		faultsSorted.addAll(faultIDs.keySet());
		Collections.sort(faultsSorted);
		
		Map<String, double[]> faultProbs = new HashMap<>();
		for (String faultName : faultsSorted)
			faultProbs.put(faultName, new double[minMags.length]);
		double[] regProbs = reg == null ? null : new double[minMags.length];
		
		DecimalFormat pDF = new DecimalFormat("0.00%");
		
		MagDependentAperiodicityOptions[] covs = {MagDependentAperiodicityOptions.HIGH_VALUES,
				MagDependentAperiodicityOptions.MID_VALUES, MagDependentAperiodicityOptions.LOW_VALUES, null};
		
		Map<String, double[]> faultRIs = new HashMap<>();
		Map<String, Map<Long, Integer>> faultDateLasts = new HashMap<>();
		Map<String, int[]> faultYears = new HashMap<>();
		for (String faultName : faultsSorted)
			faultRIs.put(faultName, new double[presets.length]);
		
		for (int p=0; p<presets.length; p++) {
			erf.setPreset(presets[p]);
			
			double presetWeight = presetWeights[p];
			
			for (int c=0; c<covs.length; c++) {
				if (covs[c] == null) {
					erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
				} else {
					erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_BPT);
					erf.setParameter(MagDependentAperiodicityParam.NAME, covs[c]);
					erf.setParameter(AleatoryMagAreaStdDevParam.NAME, 0.0);
					erf.setParameter(HistoricOpenIntervalParam.NAME, year-1875d);
					erf.getTimeSpan().setStartTime(year);
				}
				erf.getTimeSpan().setDuration(duration);
				
				double weight = presetWeight*FaultSystemSolutionERF.getWeightForCOV(covs[c]);
				
				erf.updateForecast();
				
				FaultSystemSolution sol = erf.getSolution();
				FaultSystemRupSet rupSet = sol.getRupSet();
				
				if (c == 0) {
					for (String faultName : faultsSorted) {
						HashSet<Integer> parents = new HashSet<>(faultIDs.get(faultName));
						double sumRIs = 0d;
						int numSects = 0;
						Map<Long, Integer> dateLasts = new HashMap<>();
						for (int s=0; s<rupSet.getNumSections(); s++) {
							FaultSection sect = rupSet.getFaultSectionData(s);
							if (parents.contains(sect.getParentSectionId())) {
								double ri = 1d / sol.calcTotParticRateForSect(s);
								sumRIs += ri;
								numSects++;
								long dateLast = sect.getDateOfLastEvent();
								if (dateLast > Long.MIN_VALUE) {
									if (dateLasts.containsKey(dateLast))
										dateLasts.put(dateLast, dateLasts.get(dateLast)+1);
									else
										dateLasts.put(dateLast, 1);
								}
							}
						}
						double avgRI = sumRIs/numSects;
						faultRIs.get(faultName)[p] = avgRI;
						if (!dateLasts.isEmpty())
							faultDateLasts.put(faultName, dateLasts);
					}
				}
				
				// make sure all faults exist
				HashSet<Integer> parentIDs = new HashSet<>();
				Map<Integer, List<Integer>> parentSubSectIDs = new HashMap<>();
				// Match the legacy fact sheet: keep ruptures that touch a subsection with a trace point in the region.
				for (FaultSection sect : rupSet.getFaultSectionDataList()) {
					if (!clipFaultsToReg || reg == null || sect.getFaultTrace().stream().anyMatch(reg::contains)) {
						parentIDs.add(sect.getParentSectionId());
						parentSubSectIDs.computeIfAbsent(sect.getParentSectionId(), k -> new ArrayList<>())
								.add(sect.getSectionId());
					}
				}
				
				List<BitSet> faultRupMappings = new ArrayList<>();
				
				List<List<Double>> regionalRupProbs = null;
				if (reg != null) {
					regionalRupProbs = new ArrayList<>();
					for (int m=0; m<minMags.length; m++)
						regionalRupProbs.add(new ArrayList<>());
				}
				List<List<List<Double>>> faultRupProbs = new ArrayList<>();
				
				for (String faultName : faultsSorted) {
					BitSet rupMappings = new BitSet(rupSet.getNumRuptures());
					faultRupMappings.add(rupMappings);
					List<List<Double>> myFaultProbs = new ArrayList<>();
					for (int m=0; m<minMags.length; m++)
						myFaultProbs.add(new ArrayList<>());
					faultRupProbs.add(myFaultProbs);
					for (int parentID : faultIDs.get(faultName)) {
						if (parentIDs.contains(parentID)) {
							for (int subSectID : parentSubSectIDs.get(parentID))
								for (int rupIndex : rupSet.getRupturesForSection(subSectID))
								rupMappings.set(rupIndex);
						} else {
							System.err.println("WARNING: "+faultName+" parent "+parentID
									+" doesn't exist or has no subsections in the region for "+presets[p]);
						}
					}
				}
				
				double[] rupsInRegion = reg == null ? null : rupSet.getFractRupsInsideRegion(reg, false);
				
				double minMinMag = StatUtils.min(minMags);
				
				int numFSS = erf.getNumFaultSystemSources();
				for (int sourceID=0; sourceID<erf.getNumSources(); sourceID++) {
					ProbEqkSource source = erf.getSource(sourceID);
					if (sourceID < numFSS) {
						// fault-based
						int rupIndex = erf.getFltSysRupIndexForSource(sourceID);
						// see if it matches any of our faults
						for (int i=0; i<faultRupMappings.size(); i++) {
							if (faultRupMappings.get(i).get(rupIndex)) {
								// match
								for (ProbEqkRupture rup : source) {
									double mag = rup.getMag();
									if (mag < minMinMag)
										continue;
									for (int m=0; m<minMags.length; m++)
										if ((float)mag >= (float)minMags[m])
											faultRupProbs.get(i).get(m).add(rup.getProbability());
								}
							}
						}
						if (rupsInRegion != null && rupsInRegion[rupIndex] == 0d)
							// not in region
							continue;
					} else {
						Preconditions.checkState(source.getSourceSurface().isPointSurface());
						Location loc = source.getSourceSurface().getLocationsIterator().next();
						if (reg == null || !reg.contains(loc))
							continue;
					}
					if (reg != null) {
						// if we're here, it's contained
						for (ProbEqkRupture rup : source) {
							double mag = rup.getMag();
							if (mag < minMinMag)
								continue;
							for (int m=0; m<minMags.length; m++)
								if ((float)mag >= (float)minMags[m])
									regionalRupProbs.get(m).add(rup.getProbability());
						}
					}
				}
				
				// sum probs and average in
				System.out.println("Probabilities for "+presets[p]+", "+covs[c]);
				for (int m=0; m<minMags.length; m++) {
					System.out.println("M>="+(float)minMags[m]);
					if (reg != null) {
						double regProb = FaultSysSolutionERF_Calc.calcSummedProbs(regionalRupProbs.get(m));
						System.out.println("Region:\t"+pDF.format(regProb));
						regProbs[m] += weight*regProb;
					}
					for (int i=0; i<faultsSorted.size(); i++) {
						double faultProb = FaultSysSolutionERF_Calc.calcSummedProbs(faultRupProbs.get(i).get(m));
						System.out.println(faultsSorted.get(i)+":\t"+pDF.format(faultProb));
						faultProbs.get(faultsSorted.get(i))[m] += weight*faultProb;
					}
				}
			}
		}
		
		// write CSV file
		CSVFile<String> csv = new CSVFile<>(true);
		
		List<String> header = new ArrayList<>();
		header.add("");
		for (double minMag : minMags)
			header.add("M>="+(float)minMag);
		csv.addLine(header);
		
		if (reg != null) {
			List<String> regLine = new ArrayList<>();
			regLine.add(regName);
			for (int m=0; m<minMags.length; m++)
				regLine.add(pDF.format(regProbs[m]));
			csv.addLine(regLine);
		}
		
		for (String faultName : faultsSorted) {
			List<String> line = new ArrayList<>();
			line.add(faultName);
			double[] probs = faultProbs.get(faultName);
			for (int m=0; m<minMags.length; m++)
				line.add(pDF.format(probs[m]));
			csv.addLine(line);
		}
		
		csv.writeToFile(new File(outputDir, "probs_"+year+"_"+duration+"yr.csv"));
		
		csv = new CSVFile<>(false);
		
		header = new ArrayList<>();
		header.add("");
		header.add("Avg. RI (yrs)");
		header.add("Avg. Time-Since-Last (yrs)");
		header.add("Avg. Date of Last (yr)");
		header.add("DOLE 1");
		header.add("DOLE 1 Count");
		header.add("...");
		header.add("...");
		header.add("DOLE N");
		header.add("DOLE N Count");
		csv.addLine(header);
		
		DecimalFormat yearDF = new DecimalFormat("0.0");
		
		// fault RIs and NTs
		GregorianCalendar cal = new GregorianCalendar();
		cal.clear();
		cal.set(GregorianCalendar.YEAR, year);
		long startMillis = cal.getTimeInMillis();
		System.out.println("Start millis: "+startMillis);
		
		for (String faultName : faultsSorted) {
			List<String> line = new ArrayList<>();
			line.add(faultName);
			double[] ris = faultRIs.get(faultName);
			double avgRI = StatUtils.mean(ris);
			System.out.println(faultName+" RIs:\t"+Arrays.toString(ris)+" = "+(float)avgRI);
			line.add(yearDF.format(avgRI));
			
			Map<Long, Integer> dateLasts = faultDateLasts.get(faultName);
			if (dateLasts == null) {
				line.add("N/A");
				line.add("N/A");
			} else {
				double sumTS = 0d;
				double sumYRs = 0d;
				int count = 0;
				List<Integer> years = new ArrayList<>();
				List<Integer> yearCounts = new ArrayList<>();
				for (long dateLast : dateLasts.keySet()) {
					long millisSince = startMillis - dateLast;
					double timeSince = (double)millisSince / ProbabilityModelsCalc.MILLISEC_PER_YEAR;
					int myCount = dateLasts.get(dateLast);
					sumTS += timeSince*myCount;
					count += myCount;
					GregorianCalendar cal2 = new GregorianCalendar();
					cal2.clear();
					cal2.setTimeInMillis(dateLast);
					int myYear = cal2.get(GregorianCalendar.YEAR)+1;
					years.add(myYear);
					yearCounts.add(myCount);
					sumYRs += myYear*myCount;
				}
				line.add(yearDF.format(sumTS/count));
				line.add(yearDF.format(sumYRs/count));
				for (int y=0; y<years.size(); y++) {
					line.add(years.get(y)+"");
					line.add(yearCounts.get(y)+"");
				}
			}
			csv.addLine(line);
		}
		
		csv.writeToFile(new File(outputDir, "ris_"+year+".csv"));
	}

}
