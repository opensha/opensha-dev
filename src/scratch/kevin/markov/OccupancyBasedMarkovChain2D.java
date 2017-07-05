package scratch.kevin.markov;

import java.util.Map;

import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;

import com.google.common.base.Preconditions;
import com.google.common.collect.Maps;
import com.google.common.primitives.Doubles;

/**
 * This Markov chain uses state occupancy to determine 2D state transition probabiltiies
 * @author kevin
 *
 */
public class OccupancyBasedMarkovChain2D extends MarkovChain {
	
	private Map<IndicesKey, PossibleStates> stateTransitionDataset;
	
	private PossibleStates[] marginals;
	
	public OccupancyBasedMarkovChain2D(double distSpacing, EvenlyDiscrXYZ_DataSet occ) {
		this(distSpacing, getPossibleStatesFromDiscr(occ));
	}
	
	public OccupancyBasedMarkovChain2D(double distSpacing, PossibleStates possibleStates) {
		Preconditions.checkState(possibleStates.getStates().get(0).length == 2, "Must be 2D!");
		init(2, distSpacing, possibleStates);
		
		stateTransitionDataset = Maps.newHashMap();
	}
	
	private static PossibleStates getPossibleStatesFromDiscr(EvenlyDiscrXYZ_DataSet occ) {
		PossibleStates poss = new PossibleStates(null);
		
		for (int xInd=0; xInd<occ.getNumX(); xInd++)
			for (int yInd=0; yInd<occ.getNumY(); yInd++)
				poss.add(new int[] {xInd, yInd}, occ.get(xInd, yInd));
		
		return poss;
	}

	@Override
	public MarkovChain getCollapsedChain(int... indices) {
		Preconditions.checkState(indices.length == 2 && indices[0] == 0 && indices[1] == 1,
				"Two D chain, so collapsing not relevant");
		return this;
	}

	@Override
	public synchronized void addState(int[] fromState, int[] toState) {
		PossibleStates possibleStates = getOccupancy();
		possibleStates.add(toState, 1d);
		// must clear transition map because all states that can transition to this state are now invalid
		stateTransitionDataset.clear();
		marginals = null;
	}

	@Override
	public synchronized PossibleStates getDestinationStates(int[] fromState) {
		IndicesKey key = new IndicesKey(fromState);
		PossibleStates dests = stateTransitionDataset.get(key);
		if (dests == null) {
			PossibleStates occupancy = getOccupancy();
			if (marginals == null) {
				marginals = new PossibleStates[getNDims()];
				for (int i=0; i<getNDims(); i++)
					marginals[i] = occupancy.getMarginal(i);
			}
			dests = calcDestStates(fromState, occupancy, marginals);
			stateTransitionDataset.put(key, dests);
		}
		return dests;
	}
	
	private static PossibleStates calcDestStates(int[] fromState, PossibleStates occupancy, PossibleStates[] marginals) {
		double occ = occupancy.getTot();
		
		// calculate marginal probabilities
		double pE1givenS1 = calcMarginalProb(marginals[0], fromState[0]);
		double pE2givenS2 = calcMarginalProb(marginals[1], fromState[1]);
		double pNE1givenS1 = 1d-pE1givenS1;
		double pNE2givenS2 = 1d-pE2givenS2;
		
		checkValidProb(pE1givenS1);
		checkValidProb(pE2givenS2);
		checkValidProb(pNE1givenS1);
		checkValidProb(pNE2givenS2);

		// calculate P(s1|s2) and P(s2|s1)
		double nS1 = 0;
		double nS2 = 0;
		double nS1S2 = 0;
		for (int[] state : occupancy.getStates()) {
			double freq = occupancy.getFrequency(state);
			if (state[0] == fromState[0])
				nS1 += freq;
			if (state[1] == fromState[1])
				nS2 += freq;
			if (state[0] == fromState[0] && state[1] == fromState[1])
				nS1S2 += freq;
		}
		double pS1givenS2 = nS1S2/nS2;
		checkValidProb(pS1givenS2);
		double pS2givenS1 = nS1S2/nS1;
		checkValidProb(pS2givenS1);
		double pS1 = nS1/occ;
		checkValidProb(pS1);
		double pS2 = nS2/occ;
		checkValidProb(pS2);
		// N(S1,S2)
//		double nS1S2 = occupancy.getFrequency(fromState);
//		if (nS1S2 == 0d)
//			return null;
//		double pS1S2 = nS1S2/occ;
//		checkValidProb(pS1S2);
//		// N(S1)
//		double nS1 = marginals[0].getFrequency(new int[] {fromState[0]});
//		double pS1 = nS1/occ;
//		checkValidProb(pS1);
//		// N(S2)
//		double nS2 = marginals[1].getFrequency(new int[] {fromState[1]});
//		double pS2 = nS2/occ;
//		checkValidProb(pS2);
		System.out.println("State: ["+fromState[0]+","+fromState[1]+"]");
		System.out.println("pS1givenS2="+(float)pS1givenS2+", pS2givenS1="+(float)pS2givenS1+", pS1="+(float)pS1+", pS2="+(float)pS2);
		System.out.println("pS1S2="+(float)(occupancy.getFrequency(fromState)/occ));
		System.out.println("calc 1: "+(float)((occupancy.getFrequency(fromState)/occ)/(pS1*pS2)));
		System.out.println("calc 2: "+(float)(pS2givenS1/pS2));
		System.out.println("calc 3: "+(float)(pS1givenS2/pS1));
		// N(S1+1,S2+1)
		double nS1p1S2p1 = occupancy.getFrequency(new int[] {fromState[0]+1, fromState[1]+1});
		
		// prob neither: P(!E1,!E2|S1,S2)
		double pNE1_NE2givenS1S2 = nS1p1S2p1/nS1S2;
		checkValidProb(pNE1_NE2givenS1S2);
//		Preconditions.checkState(Doubles.isFinite(pNE1_NE2givenS1S2),
//				"Not finite: "+nS1p1S2p1+"/"+nS1S2+" = "+pNE1_NE2givenS1S2);
		
//		System.out.println("pNE2givenS2="+pNE2givenS2+", pS1S2="+pS1S2+", pS1="+pS1+", pS2="+pS2);
		// P(!E2|S1,S2)		
//		double pNE2givenS1S2 = pNE2givenS2 * pS1S2 / (pS1*pS2);
		double pNE2givenS1S2 = pNE2givenS2 * pS1givenS2 / pS1;
		System.out.println((float)pNE2givenS1S2+" = "+(float)pNE2givenS2+" * "+(float)pS1givenS2+" / "+(float)pS1);
		checkValidProb(pNE2givenS1S2);
		
		// P(E1,!E2|S1,S2)
		Preconditions.checkState(pNE2givenS1S2 >= pNE1_NE2givenS1S2,
				"pNE2givenS1S2="+pNE2givenS1S2+" < pNE1_NE2givenS1S2="+pNE1_NE2givenS1S2);
		double pE1_NE2givenS1S2 = pNE2givenS1S2 - pNE1_NE2givenS1S2;
		checkValidProb(pE1_NE2givenS1S2);
		
		// P(!E1|S1,S2)
//		double pNE1givenS1S2 = pNE1givenS1 * pS1S2 / (pS1*pS2);
		double pNE1givenS1S2 = pNE1givenS1 * pS2givenS1 / pS2;
		checkValidProb(pNE1givenS1S2);
		
		// P(!E1,E2|S1,S2)
		double pNE1_E2givenS1S2 = pNE1givenS1S2 - pNE1_NE2givenS1S2;
		checkValidProb(pNE1_E2givenS1S2);
		
		// P(E1,E2|S1,S2)
		double pE1_E2givenS1S2 = 1d - (pNE1_NE2givenS1S2+pE1_NE2givenS1S2+pNE1_E2givenS1S2);
		checkValidProb(pE1_E2givenS1S2);
		
		PossibleStates dests = new PossibleStates(fromState);
		dests.add(new int[] {fromState[0]+1, fromState[1]+1}, occ*pNE1_NE2givenS1S2);
		dests.add(new int[] {0, fromState[1]+1}, occ*pE1_NE2givenS1S2);
		dests.add(new int[] {fromState[0]+1, 0}, occ*pNE1_E2givenS1S2);
		dests.add(new int[] {0, 0}, occ*pE1_E2givenS1S2);
		
		return dests;
	}
	
	private static void checkValidProb(double prob) {
		Preconditions.checkState(Doubles.isFinite(prob) && prob >= 0d && prob <= 1d, "invalid prob: "+prob);
	}
	
	public static double calcMarginalProb(PossibleStates marginal, int startIndex) {
		double freqAt = 0d;
		double freqAtOrAbove = 0d;
		
		for (int[] state : marginal.getStates()) {
			Preconditions.checkState(state.length == 1); // make sure it's actually a marginal
			double freq = marginal.getFrequency(state);
			if (state[0] == startIndex)
				freqAt += freq;
			if (state[0] >= startIndex)
				freqAtOrAbove += freq;
		}
		
		Preconditions.checkState(freqAtOrAbove >= freqAt);
		if (freqAtOrAbove == 0d)
			return 1d;
		
		return freqAt / freqAtOrAbove;
	}

}
