package scratch.kevin.simulators.synch.prediction;

import java.util.List;

public class PoissonPredictor implements Predictor {
	
	private int stateCount;
	private int[] rupCounts;
	private int nDims;
	private double[] rupProbs;

	@Override
	public String getName() {
		return "Poisson";
	}

	@Override
	public String getShortName() {
		return getName();
	}

	@Override
	public void init(List<int[]> path, double distSpacing) {
		nDims = path.get(0).length;
		rupCounts = new int[nDims];
		stateCount = 0;
		
		for (int[] state : path)
			doAddState(state);
		
		updateProbs();
	}

	public void doAddState(int[] state) {
		stateCount++;
		for (int i=0; i<nDims; i++)
			if (state[i] == 0)
				rupCounts[i]++;
	}

	@Override
	public void addState(int[] state) {
		doAddState(state);
		updateProbs();
	}
	
	private synchronized void updateProbs() {
		rupProbs = new double[nDims];
		
		for (int i=0; i<nDims; i++)
			rupProbs[i] = (double)rupCounts[i]/(double)stateCount;
	}

	@Override
	public double[] getRuptureProbabilities() {
		return rupProbs;
	}
	
	@Override
	public double[] getRuptureProbabilities(int[] prevState) {
		return getRuptureProbabilities();
	}

	@Override
	public void printDiagnostics() {
		// do nothing
	}

	@Override
	public Predictor getCollapsed(int... indexes) {
		PoissonPredictor p = new PoissonPredictor();
		p.stateCount = stateCount;
		p.rupCounts = new int[indexes.length];
		for (int i=0; i<indexes.length; i++)
			p.rupCounts[i] = rupCounts[indexes[i]];
		p.nDims = indexes.length;
		p.updateProbs();
		return p;
	}

}
