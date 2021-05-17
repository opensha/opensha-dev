package scratch.kevin.simulators;

import java.util.Collections;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.iden.RuptureIdentifier;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import scratch.kevin.markov.EmpiricalMarkovChain;
import scratch.kevin.markov.SparseNDimensionalHashDataset;

public class MarkovChainBuilder {
	
	public static EmpiricalMarkovChain build(double distSpacing, List<? extends SimulatorEvent> events,
			List<RuptureIdentifier> rupIdens) {
		return build(distSpacing, events, rupIdens, 0d);
	}
	
	public static EmpiricalMarkovChain build(double distSpacing, List<? extends SimulatorEvent> events,
			List<RuptureIdentifier> rupIdens, double startTimeOffset) {
		return new EmpiricalMarkovChain(getStatesPath(distSpacing, events, rupIdens, startTimeOffset), distSpacing);
	}
	
	public static EmpiricalMarkovChain build(double distSpacing, List<List<? extends SimulatorEvent>> matchesLists) {
		return build(distSpacing, matchesLists, 0d);
	}
	
	public static EmpiricalMarkovChain build(double distSpacing, List<List<? extends SimulatorEvent>> matchesLists,
			double startTimeOffset) {
		return new EmpiricalMarkovChain(getStatesPath(distSpacing, matchesLists, startTimeOffset), distSpacing);
	}
	
	public static List<List<? extends SimulatorEvent>> getMatchesLists(List<? extends SimulatorEvent> events,
			List<RuptureIdentifier> rupIdens) {
		List<List<? extends SimulatorEvent>> matchesLists = Lists.newArrayList();
		for (int i=0; i<rupIdens.size(); i++) {
			List<SimulatorEvent> matches = Lists.newArrayList(rupIdens.get(i).getMatches(events));
			matchesLists.add(matches);
			double[] matchRIs = new double[matches.size()-1];
			for (int j=1; j<matches.size(); j++)
				matchRIs[j-1] = matches.get(j).getTimeInYears()-matches.get(j-1).getTimeInYears();
			System.out.println(rupIdens.get(i).getName()+" mean RI: "+StatUtils.mean(matchRIs));
		}
		return matchesLists;
	}
	
	public static List<int[]> getStatesPath(double distSpacing, List<? extends SimulatorEvent> events,
			List<RuptureIdentifier> rupIdens, double startTimeOffset) {
		return getStatesPath(distSpacing, getMatchesLists(events, rupIdens), startTimeOffset);
	}
	
	public static List<int[]> getStatesPath(double distSpacing, List<List<? extends SimulatorEvent>> matchesLists,
			double startTimeOffset) {
		return getStatesPath(distSpacing, matchesLists, startTimeOffset, null);
	}
	
	public static List<int[]> getStatesPath(double distSpacing, List<List<? extends SimulatorEvent>> matchesLists,
			double startTimeOffset, List<List<? extends SimulatorEvent>> eventsForStatesList) {
		int nDims = matchesLists.size();
		
		// only used for utility methods in binning
		SparseNDimensionalHashDataset<Double> totalStatesDataset = new SparseNDimensionalHashDataset<Double>(nDims,
				distSpacing*0.5, distSpacing);
		
		Preconditions.checkState(nDims > 0);
		
		Preconditions.checkState(eventsForStatesList == null || eventsForStatesList.isEmpty(),
				"eventsForStatesList, if supplied, must be empty. will be filled in for each returned state");
		
		double maxTime = 0d;
		double startTime = Double.POSITIVE_INFINITY;
		for (List<? extends SimulatorEvent> matches : matchesLists) {
			maxTime = Math.max(maxTime, matches.get(matches.size()-1).getTimeInYears());
			startTime = Math.min(startTime, matches.get(0).getTimeInYears() + startTimeOffset);
		}
		
		Preconditions.checkState(startTimeOffset >= 0);
		if (startTimeOffset > 0) {
			List<List<? extends SimulatorEvent>> myMatches = Lists.newArrayList();
			
			for (List<? extends SimulatorEvent> matches : matchesLists) {
				List<SimulatorEvent> newList = Lists.newArrayList(matches);
				while (newList.size() > 0) {
					if (newList.get(0).getTimeInYears() < startTime)
						newList.remove(0);
					else
						break;
				}
				myMatches.add(newList);
				matchesLists = myMatches;
			}
		}
		int numSteps = (int)((maxTime - startTime)/distSpacing);
		
		List<int[]> fullPath = Lists.newArrayList();
		
		int[] lastMatchIndexBeforeWindowEnd = new int[nDims];
		for (int i=0; i<nDims; i++)
			lastMatchIndexBeforeWindowEnd[i] = -1;
		
		int[] prevState = null;
		
		System.out.println("Binning catalog into state vectors");
		
		int skippedSteps = 0;
		
		int startStep = 0;
		
		double startWindowStart = startTime + distSpacing*startStep;
		for (int n=0; n<nDims && startStep>0; n++) {
			List<? extends SimulatorEvent> myMatches = matchesLists.get(n);
			for (int i=lastMatchIndexBeforeWindowEnd[n]+1; i<myMatches.size(); i++) {
				double time = myMatches.get(i).getTimeInYears();
				if (time > startWindowStart)
					break;
				lastMatchIndexBeforeWindowEnd[n] = i;
			}
		}
		
		stepLoop:
		for (int step=startStep; step<numSteps; step++) {
//			if (step % 1000 == 0)
//				System.out.println("Markov Step "+step);
			double windowStart = startTime + distSpacing*step;
			double windowEnd = windowStart + distSpacing;
			
			List<SimulatorEvent> eventsOccurringInWindow = null;
			if (eventsForStatesList != null)
				eventsOccurringInWindow = Lists.newArrayList();
			
			for (int n=0; n<nDims; n++) {
				List<? extends SimulatorEvent> myMatches = matchesLists.get(n);
				for (int i=lastMatchIndexBeforeWindowEnd[n]+1; i<myMatches.size(); i++) {
					SimulatorEvent event = myMatches.get(i);
					double time = event.getTimeInYears();
					Preconditions.checkState(time >= windowStart);
					if (time > windowEnd)
						break;
					if (eventsOccurringInWindow != null && time >= windowStart)
						eventsOccurringInWindow.add(event);
					lastMatchIndexBeforeWindowEnd[n] = i;
				}
			}
			
			int[] curState = new int[nDims];
			
			for (int n=0; n<nDims; n++) {
				List<? extends SimulatorEvent> myMatches = matchesLists.get(n);
				
				double prevEvent;
				if (lastMatchIndexBeforeWindowEnd[n] >= 0) {
					prevEvent = myMatches.get(lastMatchIndexBeforeWindowEnd[n]).getTimeInYears();
				} else {
					// skip places at start where state not defined
					skippedSteps++;
					Preconditions.checkState(prevState == null);
					continue stepLoop;
				}
				
				double myDelta = windowEnd - prevEvent;
				curState[n] = totalStatesDataset.indexForDimVal(n, myDelta);
			}
			
			fullPath.add(curState);
			if (eventsOccurringInWindow != null) {
				Collections.sort(eventsOccurringInWindow);
				eventsForStatesList.add(eventsOccurringInWindow);
			}
			
			prevState = curState;
		}
		
		System.out.println("DONE binning catalog into state vectors (skipped "+skippedSteps+" steps)");
		
		return fullPath;
	}

}
