package scratch.kevin.simulators.bayesian;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.dom4j.DocumentException;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;
import com.google.common.collect.Range;
import com.google.common.io.Files;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.kevin.simulators.RSQSimMarkovChainBuilder;

public class U3CompareEventCalc extends FaultStateEventCalc {
	
	private List<String> faultNames;
	
	public U3CompareEventCalc(File u3SimFile, FaultSystemSolution sol, List<String[]> parentSectNames,
			Range<Double> magRange, double timeDiscretizationYears, double minAreaFract, boolean middleSubSect)
					throws IOException {
		Preconditions.checkState(!middleSubSect, "middle sub sect not yet implemented");
		List<HashSet<Integer>> matchingRupsList = new ArrayList<>();
		
		int numFaults = parentSectNames.size();
		faultNames = new ArrayList<>();

		FaultSystemRupSet rupSet = sol.getRupSet();
		Map<String, List<FaultSection>> sectsForParents = new HashMap<>();
		for (FaultSection subSect : rupSet.getFaultSectionDataList()) {
			List<FaultSection> sectsForParent = sectsForParents.get(subSect.getParentSectionName());
			if (sectsForParent == null) {
				sectsForParent = new ArrayList<>();
				sectsForParents.put(subSect.getParentSectionName(), sectsForParent);
			}
			sectsForParent.add(subSect);
		}
		
//		System.out.println("Rupture Counts");
		for (String[] parentNames : parentSectNames) {
			List<FaultSection> subSectsForParent = new ArrayList<>();
			for (String parentName : parentNames)
				subSectsForParent.addAll(sectsForParents.get(parentName));
			double totArea = 0d;
			HashSet<Integer> subSectIDs = new HashSet<>();
			for (FaultSection sect : subSectsForParent) {
				totArea += rupSet.getAreaForSection(sect.getSectionId());
				subSectIDs.add(sect.getSectionId());
			}
			
			HashSet<Integer> matchingRups = new HashSet<>();
			matchingRupsList.add(matchingRups);
			
			for (int r=0; r<rupSet.getNumRuptures(); r++) {
				if (!magRange.contains(rupSet.getMagForRup(r)))
					continue;
				double areaOnParentsRuptured = 0d;
				for (int s : rupSet.getSectionsIndicesForRup(r)) {
					if (subSectIDs.contains(s))
						areaOnParentsRuptured += rupSet.getAreaForSection(s);
				}
				if (areaOnParentsRuptured == 0d)
					continue;
				double fractRuptured = areaOnParentsRuptured/totArea;
				Preconditions.checkState((float)fractRuptured > 0f && (float)fractRuptured <= 1f, "bad fractRuptured=%s", fractRuptured);
				if (fractRuptured >= minAreaFract)
					matchingRups.add(r);
			}
			String faultName = RSQSimMarkovChainBuilder.getCombinedParentsName(rupSet.getFaultSectionDataList(), parentNames);
			faultNames.add(faultName);
//			System.out.println("\t"+faultName+":\t"+matchingRups.size()+" ruptures");
		}
		
		List<U3Rup> rups = new ArrayList<>();
		
		for (String line : Files.readLines(u3SimFile, Charset.defaultCharset())) {
			line = line.trim();
			if (line.isEmpty() || line.startsWith("nthRup"))
				continue; // header
			String[] split = line.split("\t");
			if (split.length < 4)
				// partial line
				break;
			int fssIndex = Integer.parseInt(split[1]);
			double timeYears = Double.parseDouble(split[2]);
			
			rups.add(new U3Rup(fssIndex, timeYears));
		}
		
		Map<String, Integer> countsMap = new HashMap<>();
		int[] marginalCounts = new int[numFaults];
		
		int eventStartIndex = 0;
		double startTimeYears = rups.get(0).timeYears;
		double endTimeYears = rups.get(rups.size()-1).timeYears;
		
		double durationYears = endTimeYears - startTimeYears;
		System.out.println("Loaded "+(float)durationYears+" years");
		
		int stateCount = 0;
		
		for (double windowStart=startTimeYears; windowStart+timeDiscretizationYears <= endTimeYears; windowStart += timeDiscretizationYears) {
			double windowEnd = windowStart+timeDiscretizationYears;
			
			boolean[] eventVector = new boolean[numFaults];
			
			for (int i=eventStartIndex; i<rups.size(); i++) {
				U3Rup event = rups.get(i);
				double eventTime = event.timeYears;
				if (eventTime < windowStart) {
					eventStartIndex = i;
					continue;
				} else if (eventTime >= windowEnd) {
					break;
				}
				// this means that the event is within the time window
				for (int n=0; n<numFaults; n++)
					if (matchingRupsList.get(n).contains(event.fssIndex))
						eventVector[n] = true;
			}
			
			for (int n=0; n<numFaults; n++)
				if (eventVector[n])
					marginalCounts[n]++;
			
			String binaryRep = getStringRep(eventVector);
			
			Integer prevCount = countsMap.get(binaryRep);
			if (prevCount == null)
				prevCount = 0;
			countsMap.put(binaryRep, prevCount+1);
			stateCount++;
		}
		init(numFaults, countsMap, marginalCounts, stateCount);
	}
	
	private class U3Rup {
		private int fssIndex;
		private double timeYears;
		
		public U3Rup(int fssIndex, double timeYears) {
			this.fssIndex = fssIndex;
			this.timeYears = timeYears;
		}
	}

	@Override
	public String getFaultName(int index) {
		return faultNames.get(index);
	}

	public static void main(String[] args) throws IOException, DocumentException {
		File u3BaseSimDir = new File("/home/kevin/OpenSHA/UCERF3/time_dep_catalogs/"
				+ "2018_07_17-MID_VALUES-fm3_1-geol/batch0");
		File[] u3SimDirs = {
				new File(u3BaseSimDir, "1000000yr_run0"),
				new File(u3BaseSimDir, "1000000yr_run1"),
				new File(u3BaseSimDir, "1000000yr_run2"),
				new File(u3BaseSimDir, "1000000yr_run3"),
				new File(u3BaseSimDir, "1000000yr_run4"),
				new File(u3BaseSimDir, "1000000yr_run5"),
				new File(u3BaseSimDir, "1000000yr_run6"),
				new File(u3BaseSimDir, "1000000yr_run7"),
				new File(u3BaseSimDir, "1000000yr_run8"),
				new File(u3BaseSimDir, "1000000yr_run9")
		};
		File solFile = new File("/home/kevin/OpenSHA/UCERF3/FM3_1_GEOL_MEAN_BRANCH_AVG_SOL.zip");
		FaultSystemSolution sol = FaultSystemIO.loadSol(solFile);
		
		List<Range<Double>> magRanges = new ArrayList<>();
		magRanges.add(Range.closed(7d, Double.POSITIVE_INFINITY));
		double[] timeDiscretizationsYears =  { 10d, 5d, 1d };
		double minAreaFract = 0.2;
		boolean middleSubSect = false; // else any
		List<String[]> parentSectNames = CatalogEventCalc.getParentSectsCajonPass();
		
//		double[] minMags = { 6d, 7d };
//		double[] timeDiscretizationsYears =  { 10d };
//		double minAreaFract = 0.2;
//		boolean middleSubSect = false; // else any
//		List<String[]> parentSectNames = CatalogEventCalc.getParentSectsSetOf4_SAF();
		
		for (File u3SimDir : u3SimDirs) {
			File u3SimFile = new File(u3SimDir, "sampledEventsData.txt");
			
			System.out.println("Loading: "+u3SimFile.getAbsolutePath());
			
			for (Range<Double> magRange : magRanges) {
				for (double timeDiscretizationYears : timeDiscretizationsYears) {
					U3CompareEventCalc calc = new U3CompareEventCalc(u3SimFile, sol, parentSectNames, magRange,
							timeDiscretizationYears, minAreaFract, middleSubSect);
					
					String csvPrefix = "u3_event_probs_"+parentSectNames.size()+"faults_m"
							+optionalDigitDF.format(magRange.lowerEndpoint());
					if (Double.isFinite(magRange.upperEndpoint()))
						csvPrefix += "_"+optionalDigitDF.format(magRange.upperEndpoint());
					csvPrefix += "_"+optionalDigitDF.format(timeDiscretizationYears)+"yr";
					
					System.out.println("Writing "+csvPrefix+".csv");
					calc.writeStatesCSV(new File(u3SimDir, csvPrefix+".csv"));
				}
			}
		}
	}

}
