package scratch.kevin.markov;

import java.util.List;
import java.util.Map;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class PossibleStates {
	List<int[]> states = Lists.newArrayList();
	private Map<IndicesKey, Integer> stateIndexMap = Maps.newHashMap();
	List<Double> frequencies = Lists.newArrayList();
	double tot = 0d;
	int[] fromState;
	
	public PossibleStates(int[] fromState) {
		this.fromState = fromState;
	}
	
	public void add(int[] state, double frequency) {
		Preconditions.checkState(frequency >= 0, "Frequency cannot be negative! freq="+frequency);
		IndicesKey key = new IndicesKey(state);
		Integer index = stateIndexMap.get(key);
		if (index == null) {
			stateIndexMap.put(key, states.size());
			states.add(state);
			frequencies.add(frequency);
		} else {
			frequencies.set(index, frequencies.get(index)+frequency);
		}
		tot += frequency;
	}
	
	public double getFrequency(int[] indices) {
		Integer index = stateIndexMap.get(new IndicesKey(indices));
		if (index == null)
			return 0d;
		return frequencies.get(index);
	}
	
	public int[] drawState() {
		double rand = Math.random()*tot;
		double runningTot = 0d;
		
		for (int i=0; i<states.size(); i++) {
			runningTot += frequencies.get(i);
			if (rand <= runningTot)
				return states.get(i);
		}
		throw new IllegalStateException("Frequencies don't add up...");
	}
	
	public double getTot() {
		return tot;
	}
	
	public List<int[]> getStates() {
		return states;
	}
	
	public List<int[]> getStatesMatching(int[] input) {
		List<int[]> myStates = Lists.newArrayList();
		stateLoop:
		for (int[] state : states) {
			for (int i=0; i<state.length; i++)
				if (input[i] >= 0 && input[i] != state[i])
					continue stateLoop;
			myStates.add(state);
		}
		return myStates;
	}
	
	public int getNumStates() {
		return states.size();
	}
	
	public int[] getFromState() {
		return fromState;
	}
	
	public PossibleStates getMarginal(int index) {
		int[] newFromState;
		if (fromState == null)
			newFromState = null;
		else
			newFromState = new int[] {fromState[index]};
		PossibleStates marginal = new PossibleStates(newFromState);
		
		for (int[] state : states) {
			int[] margState = {state[index]};
			marginal.add(margState, getFrequency(state));
		}
		
		return marginal;
	}
}