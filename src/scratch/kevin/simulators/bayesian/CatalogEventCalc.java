package scratch.kevin.simulators.bayesian;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.CSVFile;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.iden.RuptureIdentifier;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimMarkovChainBuilder;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class CatalogEventCalc extends FaultStateEventCalc {
	
	private List<RuptureIdentifier> faultIdentifiers;
	
	public CatalogEventCalc(List<? extends SimulatorEvent> events, double timeDiscretizationYears, List<RuptureIdentifier> faultIdentifiers) {
		this.faultIdentifiers = faultIdentifiers;
		
		int numFaults = faultIdentifiers.size();
		Map<String, Integer> countsMap = new HashMap<>();
		int[] marginalCounts = new int[numFaults];
		
		int eventStartIndex = 0;
		double startTimeYears = events.get(0).getTimeInYears();
		double endTimeYears = events.get(events.size()-1).getTimeInYears();
		
		int stateCount = 0;
		
		for (double windowStart=startTimeYears; windowStart+timeDiscretizationYears <= endTimeYears; windowStart += timeDiscretizationYears) {
			double windowEnd = windowStart+timeDiscretizationYears;
			
			boolean[] eventVector = new boolean[numFaults];
			
			for (int i=eventStartIndex; i<events.size(); i++) {
				SimulatorEvent event = events.get(i);
				double eventTime = event.getTimeInYears();
				if (eventTime < windowStart) {
					eventStartIndex = i;
					continue;
				} else if (eventTime >= windowEnd) {
					break;
				}
				// this means that the event is within the time window
				for (int n=0; n<numFaults; n++)
					if (faultIdentifiers.get(n).isMatch(event))
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
	
	public int getMarginalCount(RuptureIdentifier faultIdentifier) {
		for (int i=0; i<faultIdentifiers.size(); i++)
			if (faultIdentifiers.get(i) == faultIdentifier)
				return getMarginalCount(i);
		throw new IllegalStateException("Identifier not found: "+faultIdentifier.getName());
	}
	
	public double getMarginalProb(RuptureIdentifier faultIdentifier) {
		return (double)getMarginalCount(faultIdentifier)/(double)getTotalStateCount();
	}
	
	public String getFaultName(int index) {
		return faultIdentifiers.get(index).getName();
	}
	
//	private Iterable<boolean[]> getVisitedStatesIterableForFaults(Map<RuptureIdentifier, Boolean> faultStates) {
//		Boolean[] fixedStates = new Boolean[numFaults];
//		for (int i=0; i<numFaults; i++)
//			fixedStates[i] = faultStates.get(faultIdentifiers.get(i)); // will be null if absent from map, which is ok
//		return getVisitedStatesIterable(fixedStates);
//	}
//	
//	private Iterable<boolean[]> getVisitedStatesIterable() {
//		return getVisitedStatesIterable(null);
//	}
//	
//	private Iterable<boolean[]> getVisitedStatesIterable(Boolean[] fixedStates) {
//		final Iterator<boolean[]> it;
//		if (fixedStates == null) {
//			it = countsMap.keySet().iterator();
//		} else {
//			Preconditions.checkState(fixedStates.length == numFaults);
//			it = new VisitedStatesIterator(fixedStates);
//		}
//		
//		return new Iterable<boolean[]>() {
//			
//			@Override
//			public Iterator<boolean[]> iterator() {
//				return it;
//			}
//		};
//	}
//	
//	private class VisitedStatesIterator implements Iterator<boolean[]> {
//		
//		private final Boolean[] fixedStates;
//		
//		private Iterator<boolean[]> allStatesIterator;
//		private boolean[] next = null;
//		
//		public VisitedStatesIterator(Boolean[] fixedStates) {
//			if (fixedStates == null)
//				fixedStates = new Boolean[numFaults];
//			else
//				Preconditions.checkState(fixedStates.length == numFaults);
//			this.fixedStates = fixedStates;
//			
//			allStatesIterator = countsMap.keySet().iterator();
//		}
//
//		@Override
//		public boolean hasNext() {
//			if (next == null)
//				next = buildNext();
//			return next != null;
//		}
//
//		@Override
//		public boolean[] next() {
//			if (next != null) {
//				boolean[] ret = next;
//				next = null;
//				return ret;
//			}
//			return buildNext();
//		}
//		
//		private boolean[] buildNext() {
//			allStatesLoop:
//			while (allStatesIterator.hasNext()) {
//				boolean[] state = allStatesIterator.next();
//				for (int i=0; i<fixedStates.length; i++)
//					if (fixedStates[i] != null && fixedStates[i].booleanValue() != state[i])
//						continue allStatesLoop;
//				return state;
//			}
//			return null;
//		}
//		
//	}
	
	static List<String[]> getParentSectsSetOf9() {
		List<String[]> parentSectBundles = new ArrayList<>();
		
		parentSectBundles.add(new String[] {"San Andreas (Carrizo) rev", "San Andreas (Cholame) rev"});
		parentSectBundles.add(new String[] {"San Andreas (Big Bend)"});
		parentSectBundles.add(new String[] {"San Andreas (Mojave N)", "San Andreas (Mojave S)"});
		parentSectBundles.add(new String[] {"San Andreas (San Bernardino N)", "San Andreas (San Bernardino S)"});
		parentSectBundles.add(new String[] {"San Andreas (San Gorgonio Pass-Garnet HIll)"});
		parentSectBundles.add(new String[] {"San Andreas (Coachella) rev"});
		parentSectBundles.add(new String[] {"Garlock (West)"});
		parentSectBundles.add(new String[] {"San Jacinto (San Bernardino)", "San Jacinto (San Jacinto Valley) rev"});
		parentSectBundles.add(new String[] {"San Jacinto (Stepovers Combined)", "San Jacinto (Anza) rev"});
		
		return parentSectBundles;
	}
	
	static List<String[]> getParentSectsSetOf6() {
		List<String[]> parentSectBundles = new ArrayList<>();
		
		parentSectBundles.add(new String[] {"San Andreas (Carrizo) rev"});
		parentSectBundles.add(new String[] {"San Andreas (Cholame) rev"});
		parentSectBundles.add(new String[] {"San Andreas (Mojave S)"});
		parentSectBundles.add(new String[] {"San Andreas (Coachella) rev"});
		parentSectBundles.add(new String[] {"Garlock (West)"});
		parentSectBundles.add(new String[] {"San Jacinto (Anza) rev"});
		
		return parentSectBundles;
	}
	
	public static void main(String[] args) throws IOException {
		File catalogsDir = new File("/data/kevin/simulators/catalogs");
		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance(catalogsDir);
		
		double minMag = 7d;
		double distSpacing = 1d; // years
		double skipYears = 5000;
		boolean middleSubSect = false; // else any
		double minAreaFract = 0.2;
		
		List<String[]> parentSectBundles = getParentSectsSetOf9();
		
		List<RuptureIdentifier> faultIdens = new ArrayList<>();
		
		for (String[] parentSectNames : parentSectBundles)
			faultIdens.add(RSQSimMarkovChainBuilder.getU3_SectionIdentifier(
					catalog, minAreaFract, middleSubSect, parentSectNames));
		
		System.out.println();
		System.out.println("Faults Names");
		for (RuptureIdentifier iden : faultIdens)
			System.out.println("\t"+iden.getName());
		System.out.println();
		
		System.out.println("Loading catalog...");
		List<RSQSimEvent> events = catalog.loader().minMag(minMag).skipYears(skipYears).load();
		System.out.println("Loaded "+events.size()+" events");
		
		CatalogEventCalc calc = new CatalogEventCalc(events, distSpacing, faultIdens);
		
		String csvPrefix = "catalog_event_probs_"+faultIdens.size()+"faults_m"
				+optionalDigitDF.format(minMag)+"_"+optionalDigitDF.format(distSpacing)+"yr";
		
		calc.writeStatesCSV(new File(catalog.getCatalogDir(), csvPrefix+".csv"));
		
		// now for each half
		double origStart = events.get(0).getTime();
		double origEnd = events.get(events.size()-1).getTime();
		double midpoint = origStart + 0.5*(origEnd-origStart);
		List<RSQSimEvent> firstHalf = new ArrayList<>(events);
		firstHalf.removeIf((RSQSimEvent e) -> e.getTime() >= midpoint);
		List<RSQSimEvent> secondHalf = new ArrayList<>(events);
		secondHalf.removeIf((RSQSimEvent e) -> e.getTime() < midpoint);
		
		calc = new CatalogEventCalc(firstHalf, distSpacing, faultIdens);
		calc.writeStatesCSV(new File(catalog.getCatalogDir(), csvPrefix+"_first_half.csv"));
		calc = new CatalogEventCalc(secondHalf, distSpacing, faultIdens);
		calc.writeStatesCSV(new File(catalog.getCatalogDir(), csvPrefix+"_second_half.csv"));
	}

}
