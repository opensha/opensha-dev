package scratch.kevin.simulators.synch.prediction;

import java.util.HashSet;
import java.util.List;

import com.google.common.base.Joiner;
import com.google.common.collect.Lists;

public class SplitPredictor implements Predictor {
	
	private Predictor[] predictors;
	
	int[] prevState;
	
	public SplitPredictor(Predictor... predictors) {
		this.predictors = predictors;
	}

	@Override
	public String getShortName() {
		HashSet<String> names = new HashSet<String>();
		for (Predictor p : predictors)
			names.add(p.getShortName());
		return "SPLIT_"+Joiner.on("_").join(names);
	}

	@Override
	public String getName() {
		List<String> names = Lists.newArrayList();
		for (Predictor p : predictors)
			names.add(p.getShortName());
		return "SPLIT: "+Joiner.on(", ").join(names);
	}

	@Override
	public void init(List<int[]> path, double distSpacing) {
		for (Predictor p : predictors)
			p.init(path, distSpacing);
		
		prevState = path.get(path.size()-1);
	}

	@Override
	public void addState(int[] state) {
		for (Predictor p : predictors)
			p.addState(state);
		
		prevState = state;
	}

	@Override
	public double[] getRuptureProbabilities() {
		return getRuptureProbabilities(prevState);
	}

	@Override
	public double[] getRuptureProbabilities(int[] state) {
		double[] ret = new double[state.length];
		for (int i=0; i<state.length; i++) {
			ret[i] = predictors[i].getRuptureProbabilities(state)[i];
		}
		return ret;
	}

	@Override
	public void printDiagnostics() {
		// TODO Auto-generated method stub

	}

	@Override
	public Predictor getCollapsed(int... indexes) {
		Predictor[] predictors = new Predictor[indexes.length];
		for (int i=0; i<indexes.length; i++) {
			predictors[i] = this.predictors[indexes[i]].getCollapsed(indexes);
		}
		return new SplitPredictor(predictors);
	}

}
