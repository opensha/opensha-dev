package scratch.kevin.markov;

import java.util.List;
import java.util.Map;

import scratch.kevin.markov.MarkovPath.LoopCounter;

import com.google.common.collect.Lists;

public abstract class MarkovChain {
	
	private int nDims;
	
	private double distSpacing;
	
	private PossibleStates occupancy;

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
	protected void init(
			int nDims,
			double distSpacing,
			PossibleStates possibleInitialStates) {
		this.nDims = nDims;
		this.distSpacing = distSpacing;
		this.occupancy = possibleInitialStates;
	}

	public int getNDims() {
		return nDims;
	}

	public double getDistSpacing() {
		return distSpacing;
	}

	public PossibleStates getOccupancy() {
		return occupancy;
	}
	
	/**
	 * Returns a collapsed chain that only considers the given indices
	 * 
	 * @param indices
	 * @return
	 */
	public abstract MarkovChain getCollapsedChain(int... indices);
	
	/**
	 * Adds a state transition with unit weight from the given state to the given state.
	 * 
	 * @param fromState
	 * @param toState
	 */
	public abstract void addState(int[] fromState, int[] toState);
	
	protected static int[] getCollapsedState(int[] state, int... indices) {
		int[] collapsedState = new int[indices.length];
		for (int i=0; i<indices.length; i++)
			collapsedState[i] = state[indices[i]];
		return collapsedState;
	}
	
	/**
	 * Gets all of the possible chains of the given length between the given states.
	 * Returned lists will not include the fromState, but will include the toState.
	 * Note that toState can be reached multiple times along the path.
	 * 
	 * @param fromState
	 * @param toState -1 to indicate any state at that index
	 * @param numSteps
	 * @return
	 */
	public List<MarkovPath> getTheoreticalPathsBetweenStates(int[] fromState, int[] toState,
			int numSteps, int maxLoops, double minProb) {
		List<MarkovPath> res = getTheoreticalPathsBetweenStates(
				new LoopCounter(), fromState, fromState, toState,numSteps, maxLoops, minProb);
		pathsFinalized += res.size();
		return res;
	}
	
	/**
	 * Calculates the transition probability from fromState to toState
	 * @param fromState
	 * @param toState
	 * @return transition probability. Special cases: NaN if fromState has not been occupied, 0 if toState
	 * is not a destination state from fromState.
	 */
	public double getTransitionProb(int[] fromState, int[] toState) {
		PossibleStates states = getDestinationStates(fromState);
		if (states == null)
			return Double.NaN;
		return states.getFrequency(toState) / states.tot;
	}
	
	/**
	 * Returns all possible destination states from the given states
	 * 
	 * @param fromState
	 * @return
	 */
	public abstract PossibleStates getDestinationStates(int[] fromState);
	
	private static long pathsFinalized = 0;
	private static long pathsRejected = 0;
	
	private List<MarkovPath> getTheoreticalPathsBetweenStates(LoopCounter counter, int[] origFromState,
			int[] prevState, int[] toState, int numSteps, int maxLoops, double minProb) {
//		if (numSteps == 0) {
//			// we're checking if this state is a match
//			for (int i=0; i<fromState.length; i++) {
//				if (fromState[i] != toState[i] && toState[i] >= 0)
//					// fails
//					return null;
//			}
//			// this means it passes
//			List<List<int[]>> ret = Lists.newArrayList();
//			ret.add(Lists.newArrayList(fromState));
//			return ret;
//		}
		PossibleStates states = getDestinationStates(prevState);
		
		List<MarkovPath> paths = Lists.newArrayList();
		for (int[] possibleState : states.getStates()) {
			counter.cloneResgister(possibleState);
			if (counter.getMaxLoops() > maxLoops)
				continue;
			if (numSteps == 1) {
				// we're at the end here
				boolean match = true;
				for (int i=0; i<prevState.length; i++) {
					if (possibleState[i] != toState[i] && toState[i] >= 0) {
						// fails
						match = false;
						break;
					}
				}
				if (match) {
					MarkovPath path = new MarkovPath(origFromState);
					path.addToStart(possibleState, getTransitionProb(prevState, possibleState));
					paths.add(path);
				}
			} else {
				// in the middle
				// is going to this state too many loops?
				List<MarkovPath> subPaths = getTheoreticalPathsBetweenStates(counter, origFromState, possibleState,
						toState, numSteps-1, maxLoops, minProb);
				for (MarkovPath subPath : subPaths) {
					subPath = subPath.cloneAddToStart(possibleState, getTransitionProb(prevState, possibleState));
					if (subPath.getMaxLoops() <= maxLoops && subPath.getProbability() >= minProb) {
//						if (Math.random() < 0.00001) {
//							System.out.println("fin="+pathsFinalized+", rej="+pathsRejected+". rand: "+subPath.getPathStr());
//						}
						paths.add(subPath);
					} else {
						pathsRejected++;
					}
				}
//				int loops = countLoops(counts, possibleState);
//				if (loops <= maxLoops) {
//					Map<IndicesKey, Integer> newCounts = Maps.newHashMap();
//					newCounts.putAll(counts);
//					newCounts.put(new IndicesKey(possibleState), loops+1);
//					List<MarkovPath> subPaths = getPathsBetweenStates(origFromState, possibleState,
//							toState, numSteps-1, maxLoops, newCounts);
//					for (List<int[]> subPath : subPaths) {
//						subPath = Lists.newArrayList(subPath);
//						subPath.add(0, possibleState);
//						paths.add(subPath);
//					}
//				} else {
//					System.out.println("Bailing after too many loops!");
//				}
			}
		}
		
		return paths;
	}
	
	private static final int countLoops(Map<IndicesKey, Integer> counts, int[] newState) {
		Integer count = counts.get(new IndicesKey(newState));
		if (count == null)
			// first time this state has been encountered
			return 0;
		return count;
	}

}
