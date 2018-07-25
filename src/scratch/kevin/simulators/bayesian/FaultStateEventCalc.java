package scratch.kevin.simulators.bayesian;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.CSVFile;
import org.opensha.sha.simulators.iden.RuptureIdentifier;

import com.google.common.base.Preconditions;

public abstract class FaultStateEventCalc {
	
	static final DecimalFormat optionalDigitDF = new DecimalFormat("0.##");

	private int numFaults;
	private Map<String, Integer> countsMap;
	private int[] marginalCounts;
	private int stateCount;

	protected void init(int numFaults, Map<String, Integer> countsMap, int[] marginalCounts, int stateCount) {
		this.numFaults = numFaults;
		this.countsMap = countsMap;
		this.marginalCounts = marginalCounts;
		this.stateCount = stateCount;
	}

	protected static String getStringRep(boolean[] state) {
		String binaryRep = "";
		for (boolean val : state) {
			if (val)
				binaryRep += "1";
			else
				binaryRep += "0";
		}
		return binaryRep;
	}
	
	protected static boolean[] parseStringRep(String str) {
		boolean[] ret = new boolean[str.length()];
		for (int i=0; i<ret.length; i++) {
			char c = str.charAt(i);
			if (c == '1')
				ret[i] = true;
			else
				Preconditions.checkState(c == '0', "unexpected character parsing binary: %s", str);
		}
		return ret;
	}

	public final int getCount(boolean[] state) {
		Integer count = countsMap.get(getStringRep(state));
		if (count == null)
			return 0;
		return count;
	}
	
	public final int getTotalStateCount() {
		return stateCount;
	}
	
	public final int getNumFaults() {
		return numFaults;
	}

	public final double getCatalogProb(boolean[] state) {
		return (double)getCount(state)/(double)stateCount;
	}

	/**
	 * Calculates a probability using additive (Laplace) smoothing, with indepentent probabilities as the prior. As alpha increases,
	 * the resultant probability approaches the prior.
	 * @param state
	 * @param alpha alpha parameter, positive (or zero for catalog probability)
	 * @return
	 */
	public final double getLaplaceProb(boolean[] state, double alpha) {
		Preconditions.checkState(alpha >= 0);
		int count = getCount(state);
		return (count + alpha*getIndependentProb(state))/(stateCount + alpha);
	}
	
	public final int getMarginalCount(int index) {
		return marginalCounts[index];
	}
	
	public final double getMarginalProb(int index) {
		return (double)getMarginalCount(index)/(double)getTotalStateCount();
	}

	public final double getIndependentProb(boolean[] state) {
		double prob = 1;
		for (int n=0; n<numFaults; n++) {
			int faultCount = marginalCounts[n];
			if (!state[n])
				faultCount = stateCount - faultCount;
			prob *= (double)faultCount/(double)stateCount;
		}
		return prob;
	}
	
	public abstract String getFaultName(int index);
	
	public void writeStatesCSV(File outputFile) throws IOException {
		CSVFile<String> csv = new CSVFile<>(true);
		
		List<String> header = new ArrayList<>();
		header.add("Index");
		for (int i=0; i<getNumFaults(); i++)
			header.add(getFaultName(i));
		header.add("Binary Representation");
		header.add("Catalog Count");
		header.add("Total State Count");
		header.add("Catalog Prob");
		header.add("Independent Prob");
		csv.addLine(header);
		
		addLines(csv, new boolean[numFaults], 0);
		
		csv.sort(getNumFaults()+1, 1, csvBinaryComparator);
		
		csv.writeToFile(outputFile);
	}
	
	public void addLines(CSVFile<String> csv, boolean[] state, int pos) {
		if (pos == numFaults) {
			// actually do it
			List<String> line = new ArrayList<>();
			String binary = getStringRep(state);
			line.add(Integer.parseInt(binary, 2)+"");
			for (boolean val : state)
				line.add(val+"");
			line.add(binary);
			line.add(getCount(state)+"");
			line.add(stateCount+"");
			line.add(getCatalogProb(state)+"");
			line.add(getIndependentProb(state)+"");
			csv.addLine(line);
		} else {
			boolean[] stateFalse = Arrays.copyOf(state, numFaults);
			stateFalse[pos] = false;
			addLines(csv, stateFalse, pos+1);
			boolean[] stateTrue = Arrays.copyOf(state, numFaults);
			stateTrue[pos] = true;
			addLines(csv, stateTrue, pos+1);
		}
	}
	
	private static Comparator<String> csvBinaryComparator = new Comparator<String>() {

		@Override
		public int compare(String o1, String o2) {
			boolean[] state1 = parseStringRep(o1);
			boolean[] state2 = parseStringRep(o2);
			return Integer.compare(countTrues(state1), countTrues(state2));
		}
	};
	
	private static int countTrues(boolean[] state) {
		int count = 0;
		for (boolean val : state)
			if (val)
				count++;
		return count;
	}
	
	public static FaultStateEventCalc loadStatesCSV(File csvFile) throws IOException {
		return new LoadedFaultStateEventCalc(CSVFile.readFile(csvFile, true));
	}
	
	private static class LoadedFaultStateEventCalc extends FaultStateEventCalc {
		
		private List<String> faultNames;
		
		public LoadedFaultStateEventCalc(CSVFile<String> csv) {
			int numFaults = csv.getNumCols()-6;
			Map<String, Integer> countsMap = new HashMap<>();
			int[] marginalCounts = new int[numFaults];
			int stateCount = 0;
			
			int faultColOffset = 1; // first column is index
			int countCol = faultColOffset + numFaults + 1;
			int totalStateCol = countCol + 1;
			
			// init names
			faultNames = new ArrayList<>();
			for (int i=0; i<numFaults; i++)
				faultNames.add(csv.get(0, i+1));
			
			int checkTotalState = -1;
			
			for (int row=1; row<csv.getNumRows(); row++) {
				int count = Integer.parseInt(csv.get(row, countCol));
				boolean[] state = new boolean[numFaults];
				for (int i=0; i<state.length; i++) {
					state[i] = Boolean.parseBoolean(csv.get(row, i+faultColOffset));
					if (state[i])
						marginalCounts[i] += count;
				}
				stateCount += count;
				if (row == 1)
					checkTotalState = Integer.parseInt(csv.get(row, totalStateCol));
				String key = getStringRep(state);
				Preconditions.checkState(!countsMap.containsKey(key));
				countsMap.put(key, count);
			}
			
			Preconditions.checkState(stateCount == checkTotalState,
					"State count mismatch. File is %s, reconstructed is %s", checkTotalState, stateCount);
			
			init(numFaults, countsMap, marginalCounts, stateCount);
		}

		@Override
		public String getFaultName(int index) {
			return faultNames.get(index);
		}
		
	}
	
//	public static void main(String[] args) throws IOException {
//		File intputFile = new File("/data/kevin/simulators/catalogs/rundir2585_1myr/catalog_event_probs.csv");
//		File outputFile = new File("/tmp/parsed_catalog_event_probs.csv");
//		LoadedFaultStateEventCalc loaded = new LoadedFaultStateEventCalc(CSVFile.readFile(intputFile, true));
//		loaded.writeStatesCSV(outputFile);
//	}

}