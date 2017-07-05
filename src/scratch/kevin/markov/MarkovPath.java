package scratch.kevin.markov;

import java.util.Iterator;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class MarkovPath implements Iterable<int[]> {
	
	private int[] fromState;
	private List<int[]> path;
	private double probability;
	private LoopCounter counter;
	
	public MarkovPath(int[] fromState) {
		this.fromState = fromState;
		path = Lists.newArrayList();
		probability = 1;
		counter = new LoopCounter();
		counter.register(fromState);
	}
	
	public void addToStart(int[] newState, double transProb) {
		add(0, newState, transProb);
	}
	
	public void add(int index, int[] newState, double transProb) {
		path.add(index, newState);
		probability *= transProb;
		counter.register(newState);
	}
	
	public Object clone() {
		MarkovPath other = new MarkovPath(fromState);
		other.path.addAll(path);
		other.probability = probability;
		other.counter.counts.putAll(counter.counts);
		other.counter.maxLoops = counter.maxLoops;
		return other;
	}
	
	public MarkovPath cloneAddToStart(int[] newState, double transProb) {
		MarkovPath other = (MarkovPath)clone();
		other.addToStart(newState, transProb);
		return other;
	}
	
	public int getMaxLoops() {
		return counter.maxLoops;
	}
	
	public double getProbability() {
		return probability;
	}
	
	public int size() {
		return path.size();
	}

	@Override
	public Iterator<int[]> iterator() {
		return path.iterator();
	}
	
	public int[] get(int index) {
		return path.get(index);
	}
	
	public String getPathStr() {
		String str = "PATH: ["+fromState[0]+","+fromState[1]+"]: ";
		for (int[] elem : path)
			str += " ["+elem[0]+","+elem[1]+"]";
		return str;
	}
	
	public static class LoopCounter {
		private Map<IndicesKey, Integer> counts;
		private int maxLoops = 0;
		public LoopCounter() {
			counts = Maps.newHashMap();
		}
		
		public void register(int[] state) {
			IndicesKey key = new IndicesKey(state);
			Integer count = counts.get(key);
			if (count == null)
				count = 1;
			else
				count += 1;
			counts.put(key, count);
			if (count-1 > maxLoops)
				maxLoops = count-1;
		}
		
		public int getMaxLoops() {
			return maxLoops;
		}
		
		public LoopCounter cloneResgister(int[] state) {
			LoopCounter o = new LoopCounter();
			o.counts.putAll(counts);
			o.maxLoops = maxLoops;
			o.register(state);
			return o;
		}
	}
}
