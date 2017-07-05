package scratch.kevin.simulators.synch.prediction;

import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.function.EvenlyDiscretizedFunc;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public class RecurrIntervalPredictor implements Predictor {
	
	private List<List<Double>> rupFreqs;
	private List<List<Double>> occFreqs;
	private int[] prevState;
	private int nDims;
	
	private double distSpacing;
	
	private int stateCount = 0;

	@Override
	public String getName() {
		return "Single Fault Recurrence Interval";
	}

	@Override
	public String getShortName() {
		return "RI";
	}

	@Override
	public void init(List<int[]> initialPath, double distSpacing) {
		rupFreqs = Lists.newArrayList();
		occFreqs = Lists.newArrayList();
		nDims = initialPath.get(0).length;
		this.distSpacing = distSpacing;
		
		for (int i=0; i<nDims; i++) {
			rupFreqs.add(new ArrayList<Double>());
			occFreqs.add(new ArrayList<Double>());
		}
		
		for (int[] state : initialPath)
			addState(state);
	}

	@Override
	public void addState(int[] state) {
		Preconditions.checkArgument(state != null);
		if (prevState != null) {
			for (int i=0; i<nDims; i++) {
				boolean rupture = state[i] == 0;
				int faultState = prevState[i];
				ensureListsCapacity(i, faultState+1);
				increment(occFreqs.get(i), faultState);
				if (rupture)
					increment(rupFreqs.get(i), faultState);
			}
		}
		
		prevState = state;
		
		stateCount++;
	}
	
	private void increment(List<Double> list, int index) {
		double prev = list.get(index);
		list.set(index, prev+1d);
	}
	
	private void ensureListsCapacity(int index, int size) {
		ensureListCapacity(size, rupFreqs.get(index));
		ensureListCapacity(size, occFreqs.get(index));
	}
	
	private void ensureListCapacity(int size, List<Double> list) {
		while (list.size() < size)
			list.add(0d);
	}

	@Override
	public double[] getRuptureProbabilities() {
		return getRuptureProbabilities(prevState);
	}
	
	@Override
	public double[] getRuptureProbabilities(int[] curState) {
		Preconditions.checkNotNull(curState);
		Preconditions.checkNotNull(occFreqs);
		double[] ret = new double[nDims];
		
		for (int i=0; i<nDims; i++) {
			if (curState[i] >= occFreqs.get(i).size())
				// never been here
				continue;
			double occFreq = occFreqs.get(i).get(curState[i]);
			double rupFreq = rupFreqs.get(i).get(curState[i]);
			ret[i] = rupFreq/occFreq;
		}
		
		return ret;
	}
	
	public EvenlyDiscretizedFunc getOccupancyDist(int index, int size) {
		EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(0.5*distSpacing, size, distSpacing);
		
		List<Double> occs = occFreqs.get(index);
		
		for (int i=0; i<size; i++) {
			double occ;
			if (i >= occs.size())
				occ = 0d;
			else
				occ = occs.get(i);
			
			func.set(i, occ);
		}
		
		func.scale(1d/func.calcSumOfY_Vals());
		
		return func;
	}

	@Override
	public void printDiagnostics() {
		// do nothing
	}

	@Override
	public Predictor getCollapsed(int... indexes) {
		RecurrIntervalPredictor p = new RecurrIntervalPredictor();
		p.rupFreqs = Lists.newArrayList();
		p.occFreqs = Lists.newArrayList();
		p.prevState = new int[indexes.length];
		p.nDims = indexes.length;
		for (int i=0; i<indexes.length; i++) {
			int index = indexes[i];
			p.rupFreqs.add(rupFreqs.get(index));
			p.occFreqs.add(occFreqs.get(index));
			p.prevState[i] = prevState[index];
		}
		return p;
	}
	
	public int getStateCount() {
		return stateCount;
	}

}
