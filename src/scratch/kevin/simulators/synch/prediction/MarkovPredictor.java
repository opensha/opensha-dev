package scratch.kevin.simulators.synch.prediction;

import java.util.List;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.primitives.Doubles;

import scratch.kevin.markov.EmpiricalMarkovChain;
import scratch.kevin.markov.MarkovChain;
import scratch.kevin.markov.PossibleStates;

public class MarkovPredictor implements Predictor {
	
	MarkovChain chain;
	private int[] prevState;
	private int numMisses = 0;
	private int totPredictions = 0;
	
	private Predictor backupPredictor;
	
	public MarkovPredictor() {
		this(null);
	}
	
	public MarkovPredictor(Predictor backupPredictor) {
		this.backupPredictor = backupPredictor;
	}

	@Override
	public String getName() {
		if (chain == null)
			return "Markov";
		if (backupPredictor == null)
			return chain.getNDims()+"D Markov";
		return chain.getNDims()+"D Markov (back="+backupPredictor.getShortName()+")";
	}

	@Override
	public String getShortName() {
		if (chain == null)
			return "Markov";
		if (backupPredictor == null)
			return chain.getNDims()+"DMarkov";
		return chain.getNDims()+"DMarkov_back"+backupPredictor.getShortName();
	}

	@Override
	public void init(List<int[]> path, double distSpacing) {
		chain = new EmpiricalMarkovChain(path, distSpacing);
		if (backupPredictor != null)
			backupPredictor.init(path, distSpacing);
		prevState = path.get(path.size()-1);
	}

	@Override
	public void addState(int[] state) {
		chain.addState(prevState, state);
		if (backupPredictor != null)
			backupPredictor.addState(state);
		prevState = state;
	}

	@Override
	public double[] getRuptureProbabilities() {
		return getRuptureProbabilities(prevState);
	}
	
	@Override
	public double[] getRuptureProbabilities(int[] prevState) {
		if (chain.getDestinationStates(prevState) == null)
			numMisses++;
		totPredictions++;
		return getRuptureProbabilities(chain, backupPredictor, prevState);
	}
	
	public static double[] getRuptureProbabilities(MarkovChain chain, Predictor backupPredictor, int[] prevState) {
		double[] ret = new double[prevState.length];
		
		PossibleStates possible = chain.getDestinationStates(prevState);
		if (possible != null) {
			// first fill ret with frequencies
			List<int[]> states = possible.getStates();
			double totFreq = 0;
			for (int i=0; i<states.size(); i++) {
				int[] state = states.get(i);
				double freq = possible.getFrequency(state);
				Preconditions.checkState(Doubles.isFinite(freq) && freq >= 0d);
				
				for (int j=0; j<state.length; j++)
					if (state[j] == 0)
						ret[j] += freq;
				
				totFreq += freq;
			}
			// now normalize
			if (totFreq > 0d)
				for (int i=0; i<ret.length; i++)
					ret[i] = ret[i]/totFreq;
		} else {
			if (backupPredictor != null)
				ret = backupPredictor.getRuptureProbabilities();
		}
		
		return ret;
	}

	@Override
	public void printDiagnostics() {
		String str = "Markov Misses: "+numMisses+"/"+totPredictions+" ("
				+100f*((float)numMisses/(float)totPredictions)+" %)";
		if (backupPredictor != null)
			str += ", replaced with probs from "+backupPredictor.getName();
		System.out.println(str);
	}

	@Override
	public Predictor getCollapsed(int... indexes) {
		Predictor b = null;
		if (this.backupPredictor != null)
			b = this.backupPredictor.getCollapsed(indexes);
		MarkovPredictor p = new MarkovPredictor(b);
		p.chain = chain.getCollapsedChain(indexes);
		return p;
	}

}
