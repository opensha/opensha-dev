package scratch.kevin.markov;

import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

/**
 * This creates an empirical Markov chain for the given catalog
 * @author kevin
 *
 */
public class EmpiricalMarkovChain extends MarkovChain {
	
	private double firstBinCenter;
	
	private SparseNDimensionalHashDataset<Double> totalStatesDataset;
	private SparseNDimensionalHashDataset<PossibleStates> stateTransitionDataset;
	
	// this stores the actual path
	private List<int[]> fullPath;
	private Map<IndicesKey, List<Integer>> stateIndexesMap;
	
	private transient Map<IndicesKey, Collection<int[]>> parentStatesMap;
	
	public EmpiricalMarkovChain(List<int[]> fullPath, double distSpacing) {
		init(fullPath, distSpacing);
	}
	
	private void init(List<int[]> path, double distSpacing) {
		Preconditions.checkArgument(!path.isEmpty(), "path cannot be empty");
		int nDims = path.get(0).length;
		for (int[] state : path)
			Preconditions.checkState(state.length == nDims);
		this.firstBinCenter = distSpacing*0.5;
		this.fullPath = Lists.newArrayList();
		
		totalStatesDataset = new SparseNDimensionalHashDataset<Double>(nDims, firstBinCenter, distSpacing);
		stateTransitionDataset = new SparseNDimensionalHashDataset<PossibleStates>(nDims, firstBinCenter, distSpacing);
		
		stateIndexesMap = Maps.newHashMap();
		
		PossibleStates possibleInitialStates = new PossibleStates(null);
		
		super.init(nDims, distSpacing, possibleInitialStates);
		
		System.out.println("Assembling state transition probabilities");
		
		for (int index=0; index<path.size(); index++) {
			int[] curState = path.get(index);
			
			addState(curState);
		}
		
		System.out.println("DONE assembling state transition probabilities");
	}
	
	@Override
	public void addState(int[] fromState, int[] toState) {
		// register current state
		Double stateCount = totalStatesDataset.get(toState);
		if (stateCount == null)
			stateCount = 0d;
		stateCount += 1d;
		totalStatesDataset.set(toState, stateCount);
		
		PossibleStates possibleInitialStates = getOccupancy();

		// register this state as a transition from the previous state
		if (fromState != null) {
			PossibleStates possibilities = stateTransitionDataset.get(fromState);
			if (possibilities == null) {
				possibilities = new PossibleStates(fromState);
				stateTransitionDataset.set(fromState, possibilities);
			}
			possibilities.add(toState, 1d);
			possibleInitialStates.add(toState, 1d);

			IndicesKey key = new IndicesKey(toState);
			List<Integer> indexesForState = stateIndexesMap.get(key);
			if (indexesForState == null) {
				indexesForState = Lists.newArrayList();
				stateIndexesMap.put(key, indexesForState);
			}
			indexesForState.add(fullPath.size());
		}
		
		fullPath.add(toState);
	}
	
	public void addState(int[] curState) {
		int[] prevState;
		if (fullPath.isEmpty())
			prevState = null;
		else
			prevState = fullPath.get(fullPath.size()-1);
		
		addState(prevState, curState);
	}

	/**
	 * Just for use interally generating collapsed chains
	 * 
	 * @param nDims
	 * @param distSpacing
	 * @param firstBinCenter
	 * @param totalStatesDataset
	 * @param stateTransitionDataset
	 * @param possibleInitialStates
	 */
	private EmpiricalMarkovChain(
			int nDims,
			double distSpacing,
			double firstBinCenter,
			SparseNDimensionalHashDataset<Double> totalStatesDataset,
			SparseNDimensionalHashDataset<PossibleStates> stateTransitionDataset,
			PossibleStates possibleInitialStates,
			List<int[]> fullPath,
			Map<IndicesKey, List<Integer>> stateIndexesMap) {
		this.firstBinCenter = firstBinCenter;
		this.totalStatesDataset = totalStatesDataset;
		this.stateTransitionDataset = stateTransitionDataset;
		this.fullPath = fullPath;
		this.stateIndexesMap = stateIndexesMap;
		
		super.init(nDims, distSpacing, possibleInitialStates);
	}

	public double getFirstBinCenter() {
		return firstBinCenter;
	}

	public SparseNDimensionalHashDataset<Double> getTotalStatesDataset() {
		return totalStatesDataset;
	}

	public SparseNDimensionalHashDataset<PossibleStates> getStateTransitionDataset() {
		return stateTransitionDataset;
	}
	
	@Override
	public EmpiricalMarkovChain getCollapsedChain(int... indices) {
		int nDims = indices.length;
		
		SparseNDimensionalHashDataset<Double> collapsedTotalStatesDataset =
				new SparseNDimensionalHashDataset<Double>(nDims, firstBinCenter, getDistSpacing());
		SparseNDimensionalHashDataset<PossibleStates> collapsedStateTransitionDataset =
				new SparseNDimensionalHashDataset<PossibleStates>(nDims, firstBinCenter, getDistSpacing());
		
		PossibleStates initials = new PossibleStates(null);
		
		for (int[] state : this.stateTransitionDataset.getPopulatedIndices()) {
			int[] collapsedState = getCollapsedState(state, indices);
			double tot = this.totalStatesDataset.get(state);
			PossibleStates states = this.stateTransitionDataset.get(state);
			
			Double collapsedTot = collapsedTotalStatesDataset.get(collapsedState);
			PossibleStates collapsedPossible = collapsedStateTransitionDataset.get(collapsedState);
			if (collapsedTot == null) {
				// first time in state
				collapsedTot = 0d;
				Preconditions.checkState(collapsedPossible == null);
				collapsedPossible = new PossibleStates(collapsedState);
			}
			collapsedTot += tot;
			for (int i=0; i<states.getNumStates(); i++) {
				int[] pState = states.getStates().get(i);
				double pFreq = states.getFrequency(pState);
				int[] collapsedPState = getCollapsedState(pState, indices);
				collapsedPossible.add(collapsedPState, pFreq);
				initials.add(collapsedPState, pFreq);
			}
			
			collapsedTotalStatesDataset.set(collapsedState, collapsedTot);
			collapsedStateTransitionDataset.set(collapsedState, collapsedPossible);
		}
		
		List<int[]> collapsedFullPath = Lists.newArrayList();
		Map<IndicesKey, List<Integer>> collapsedStateIndexesMap = Maps.newHashMap();
		for (int[] state : fullPath) {
			int[] collapsedState = getCollapsedState(state, indices);
			collapsedFullPath.add(collapsedState);
			IndicesKey key = new IndicesKey(collapsedState);
			List<Integer> indexesForState = collapsedStateIndexesMap.get(key);
			if (indexesForState == null) {
				indexesForState = Lists.newArrayList();
				collapsedStateIndexesMap.put(key, indexesForState);
			}
			indexesForState.add(collapsedFullPath.size()-1);
		}
		
		return new EmpiricalMarkovChain(nDims, getDistSpacing(), this.firstBinCenter,
				collapsedTotalStatesDataset, collapsedStateTransitionDataset, initials,
				collapsedFullPath, collapsedStateIndexesMap);
	}
	
	/**
	 * Returns a new Markov chain where the states have been shifted by the given amount.
	 * @param shifts
	 * @return
	 */
	public EmpiricalMarkovChain getShiftedChain(int... shifts) {
		int nDims = getNDims();
		Preconditions.checkArgument(shifts.length == nDims,
				"must supply shift for each dimension (0 means no shift in that dimension)");
		List<int[]> newPath = Lists.newArrayList();
		
		Map<IndicesKey, List<Integer>> newStateIndexesMap = Maps.newHashMap();
		
		stateLoop:
		for (int i=0; i<fullPath.size(); i++) {
			int[] newState = new int[nDims];
			for (int n=0; n<nDims; n++) {
				int index = i + shifts[n];
				if (index < 0)
					continue stateLoop;
				if (index >= fullPath.size())
					break stateLoop;
				newState[n] = fullPath.get(index)[n];
			}
			newPath.add(newState);
			
			IndicesKey key = new IndicesKey(newState);
			List<Integer> indexesForState = newStateIndexesMap.get(key);
			if (indexesForState == null) {
				indexesForState = Lists.newArrayList();
				newStateIndexesMap.put(key, indexesForState);
			}
			indexesForState.add(newPath.size()-1);
		}
		
//		System.out.println("Shifted chain has "+newPath.size()+" states (orig had "+fullPath.size()+")");
		Preconditions.checkState(newPath.size() > 1);
		
		double distSpacing = getDistSpacing();
		
		// now buid the chain
		SparseNDimensionalHashDataset<Double> newTotalStatesDataset =
				new SparseNDimensionalHashDataset<Double>(nDims, firstBinCenter, distSpacing);
		SparseNDimensionalHashDataset<PossibleStates> newStateTransitionDataset =
				new SparseNDimensionalHashDataset<PossibleStates>(nDims, firstBinCenter, distSpacing);
		
		int[] lastMatchIndexBeforeWindowEnd = new int[nDims];
		for (int i=0; i<nDims; i++)
			lastMatchIndexBeforeWindowEnd[i] = -1;
		
		int[] prevState = newPath.get(0);
		
		PossibleStates newPossibleInitialStates = new PossibleStates(null);
		
//		System.out.println("Assembling state transition probabilities");
		
		for (int i=1; i<newPath.size(); i++) {
			int[] curState = newPath.get(i);
			
			// register current state
			Double stateCount = newTotalStatesDataset.get(curState);
			if (stateCount == null)
				stateCount = 0d;
			stateCount += 1d;
			newTotalStatesDataset.set(curState, stateCount);

			// register this state as a transition from the previous state
			if (prevState != null) {
				PossibleStates possibilities = newStateTransitionDataset.get(prevState);
				if (possibilities == null) {
					possibilities = new PossibleStates(prevState);
					newStateTransitionDataset.set(prevState, possibilities);
				}
				possibilities.add(curState, 1d);
				newPossibleInitialStates.add(curState, 1d);
			}
			
			prevState = curState;
		}
		
		// debug
//		List<int[]> populatedIndices = newStateTransitionDataset.getPopulatedIndices();
//		Collections.shuffle(populatedIndices);
//		for (int i=0; i<3; i++) {
//			int[] state = populatedIndices.get(i);
//			System.out.print("Debug for ["+state[0]+","+state[1]+"]. occup="+newTotalStatesDataset.get(state)+", dest states: ");
//			PossibleStates possible = newStateTransitionDataset.get(state);
//			for (int[] dest : possible.getStates())
//				System.out.print(" "+possible.getFrequency(dest)+"["+dest[0]+","+dest[1]+"]");
//			System.out.println();
//		}
		
		return new EmpiricalMarkovChain(nDims, distSpacing, firstBinCenter, newTotalStatesDataset,
				newStateTransitionDataset, newPossibleInitialStates, newPath, newStateIndexesMap);
	}
	
	// TODO must implement fullPath to be re-enabled
//	/**
//	 * This gets a reversed Markov chian with origination frequencies, instead of transition frequencies
//	 * @return
//	 */
//	public synchronized MarkovChainBuilder getReversedChain() {
//		if (reversed == null) {
//			Map<IndicesKey, Collection<int[]>> parentStatesMap = getParentStatesMap();
//			SparseNDimensionalHashDataset<PossibleStates> reversedStateTransitionDataset =
//					new SparseNDimensionalHashDataset<PossibleStates>(nDims, firstBinCenter, distSpacing);
//			
//			for (int[] state : this.stateTransitionDataset.getPopulatedIndices()) {
//				Collection<int[]> parentStates = parentStatesMap.get(new IndicesKey(state));
//				PossibleStates possible = new PossibleStates(state);
//				if (parentStates != null) {
//					for (int[] parentState : parentStates) {
//						// this is the number of times we went from the parent state to the current state
//						double transFreq = stateTransitionDataset.get(parentState).getFrequency(state);
//						possible.add(parentState, transFreq);
//					}
//				}
//				reversedStateTransitionDataset.set(state, possible);
//			}
//			
//			reversed = new MarkovChainBuilder(nDims, this.distSpacing, this.firstBinCenter,
//					totalStatesDataset, reversedStateTransitionDataset, possibleInitialStates);
//		}
//		return reversed;
//	}
	
	/**
	 * This calculates a map from child states to all possible parent states
	 * @return
	 */
	public synchronized Map<IndicesKey, Collection<int[]>> getParentStatesMap() {
		if (parentStatesMap == null) {
			parentStatesMap = Maps.newHashMap();
			
			for (int[] state : stateTransitionDataset.getPopulatedIndices()) {
//				System.out.println("Parent state: ["+state[0]+","+state[1]+"]");
//				if (state[0] == 2 && state[1] == 2)
//					System.out.println("I'm at [2,2]");
				PossibleStates poss = stateTransitionDataset.get(state);
				for (int[] toState : poss.getStates()) {
//					System.out.println("\tChild state: ["+toState[0]+","+t/oState[1]+"]");
//					if (toState[0] == 3 && toState[1] == 3)
//						System.out.println("Found a parent for [3,3]: ["+state[0]+","+state[1]+"]");
//					else if (state[0] == 2 && state[1] == 2)
//						System.out.println("Found a different child for [2,2]: ["+toState[0]+","+toState[1]+"]");
					IndicesKey key = new IndicesKey(toState);
					Collection<int[]> fromStates = parentStatesMap.get(key);
					if (fromStates == null) {
//						fromStates = new HashSet<int[]>();
						fromStates = Lists.newArrayList();
						parentStatesMap.put(key, fromStates);
					}
					fromStates.add(state);
				}
			}
		}
		return parentStatesMap;
	}

	@Override
	public PossibleStates getDestinationStates(int[] fromState) {
		return stateTransitionDataset.get(fromState);
	}
	
	public double getActualTransPathsProbBetweenStates(int[] fromState, int[] toState, int numSteps,
			int[]... requiredSubsequentStates) {
		Preconditions.checkState(requiredSubsequentStates.length <= (int)Math.abs(numSteps));
		List<Integer> indexes = stateIndexesMap.get(new IndicesKey(fromState));
		int count = 0;
		int numChecked = 0;
		stateLoop:
		for (int fromIndex : indexes) {
			// numSteps can be negative
			int toIndex = fromIndex + numSteps;
			if (toIndex < 0 || toIndex > fullPath.size())
				continue;
			// now check subsequent states if applicable
			if (requiredSubsequentStates.length > 0) {
				int cnt = 0;
				if (numSteps < 0) {
					// backwards
					for (int i=fromIndex; --i>=toIndex && cnt < requiredSubsequentStates.length;)
						if (!Arrays.equals(fullPath.get(i), requiredSubsequentStates[cnt++]))
							continue stateLoop;
				} else {
					// forwards
					for (int i=fromIndex; ++i<=toIndex && cnt < requiredSubsequentStates.length;)
						if (!Arrays.equals(fullPath.get(i), requiredSubsequentStates[cnt++]))
							continue stateLoop;
				}
			}
			numChecked++;
			int[] testToState = fullPath.get(toIndex);
			for (int i=0; i<toState.length; i++) {
				if (toState[i] >= 0 && testToState[i] != toState[i]) {
					// not a match
					continue stateLoop;
				}
			}
			count++;
		}
		Preconditions.checkState(count <= indexes.size());
		double prob = (double)count/(double)numChecked;
		if (Double.isNaN(prob))
			return 0d;
		return prob;
	}
	
	public List<int[]> getFullPath() {
		return Collections.unmodifiableList(fullPath);
	}

}
