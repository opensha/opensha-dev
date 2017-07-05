package scratch.kevin.markov;

import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

public class XYZBasedMarkovChain extends MarkovChain {

	private EvenlyDiscrXYZ_DataSet occupancy;
	private EvenlyDiscrXYZ_DataSet prob_E1_E2;
	private EvenlyDiscrXYZ_DataSet prob_E1_nE2;
	private EvenlyDiscrXYZ_DataSet prob_nE1_E2;
	private EvenlyDiscrXYZ_DataSet prob_nE1_nE2;
	
	public XYZBasedMarkovChain(EvenlyDiscrXYZ_DataSet occupancy,
			EvenlyDiscrXYZ_DataSet prob_E1_E2, EvenlyDiscrXYZ_DataSet prob_E1_nE2,
			EvenlyDiscrXYZ_DataSet prob_nE1_E2, EvenlyDiscrXYZ_DataSet prob_nE1_nE2) {
		this.occupancy = occupancy;
		this.prob_E1_E2 = prob_E1_E2;
		this.prob_E1_nE2 = prob_E1_nE2;
		this.prob_nE1_E2 = prob_nE1_E2;
		this.prob_nE1_nE2 = prob_nE1_nE2;
		
		PossibleStates possible = new PossibleStates(null);
		for (int xInd=0; xInd<prob_E1_E2.getNumX(); xInd++) {
			for (int yInd=0; yInd<prob_E1_E2.getNumY(); yInd++) {
				possible.add(new int[] {xInd, yInd}, occupancy.get(xInd, yInd));
			}
		}
		init(2, occupancy.getGridSpacingX(), possible);
	}

	@Override
	public MarkovChain getCollapsedChain(int... indices) {
		if (indices.length != 2)
			throw new UnsupportedOperationException("XYZ Based Chains Must Be 2 Dimensionsal");
		Preconditions.checkState(indices[0] == 0 && indices[1] == 1);
		return this;
	}

	@Override
	public void addState(int[] fromState, int[] toState) {
		Preconditions.checkState((float)occupancy.getSumZ() > 1f, "Occupancy must be actual count not normalized to add a state");
		int xIndStart = fromState[0];
		int yIndStart = fromState[1];
		int xIndDest = toState[0];
		int yIndDest = toState[1];
		
		if (xIndStart >= occupancy.getNumX() || xIndDest >= occupancy.getNumX()
				|| yIndStart >= occupancy.getNumY() || yIndDest >= occupancy.getNumY())
			// do nothing, outside of bounds
			return;
		// register this state as a transition from the previous state
		if (fromState != null) {
			// do minus 1 here because in a previous addState(...) call the occupancy has already been incremented
			// to reflect this occupation, but the transition probability not yet updated
			double prevFromOcc = occupancy.get(xIndStart, yIndStart)-1;
			if (prevFromOcc < 0 || !Doubles.isFinite(prevFromOcc))
				prevFromOcc = 0;
			double p_E1_E2 = prob_E1_E2.get(xIndStart, yIndStart)*prevFromOcc;
			double p_E1_nE2 = prob_E1_nE2.get(xIndStart, yIndStart)*prevFromOcc;
			double p_nE1_E2 = prob_nE1_E2.get(xIndStart, yIndStart)*prevFromOcc;
			double p_nE1_nE2 = prob_nE1_nE2.get(xIndStart, yIndStart)*prevFromOcc;
			if (!Doubles.isFinite(p_E1_E2)) p_E1_E2 = 0d;
			if (!Doubles.isFinite(p_E1_nE2)) p_E1_E2 = 0d;
			if (!Doubles.isFinite(p_nE1_E2)) p_E1_E2 = 0d;
			if (!Doubles.isFinite(p_nE1_nE2)) p_E1_E2 = 0d;
			if (xIndDest == 0 && yIndDest == 0)
				p_E1_E2 += 1d;
			else if (xIndDest == 0 && yIndDest == yIndStart+1)
				p_E1_nE2 += 1d;
			else if (xIndDest == xIndStart+1 && yIndDest == 0)
				p_nE1_E2 += 1d;
			else if (xIndDest == xIndStart+1 && yIndDest == yIndStart+1)
				p_nE1_nE2 += 1d;
			else
				throw new IllegalStateException("Impossible to transition from ["+xIndStart+","+yIndStart
						+"] to ["+xIndDest+","+yIndDest+"]");
			// re normalize
			double sum = p_E1_E2 + p_E1_nE2 + p_nE1_E2 + p_nE1_nE2;
			prob_E1_E2.set(xIndStart, yIndStart, p_E1_E2/sum);
			prob_E1_nE2.set(xIndStart, yIndStart, p_E1_nE2/sum);
			prob_nE1_E2.set(xIndStart, yIndStart, p_nE1_E2/sum);
			prob_nE1_nE2.set(xIndStart, yIndStart, p_nE1_nE2/sum);
		}
		// now increment toState
		double prevOcc = occupancy.get(xIndDest, yIndDest);
		if (!Doubles.isFinite(prevOcc))
			prevOcc = 0d;
		occupancy.set(xIndDest, yIndDest, prevOcc+1);
	}

	@Override
	public PossibleStates getDestinationStates(int[] fromState) {
		PossibleStates possible = new PossibleStates(fromState);
		int xInd = fromState[0];
		int yInd = fromState[1];
		if (xInd >= occupancy.getNumX() || yInd >= occupancy.getNumY())
			return null;
		if (occupancy.get(xInd, yInd) == 0d || !Doubles.isFinite(occupancy.get(xInd, yInd)))
			return null;
		try {
			possible.add(new int[] {xInd+1, yInd+1},  nanAsZero(prob_nE1_nE2.get(xInd, yInd)));
			possible.add(new int[] {0, yInd+1},  nanAsZero(prob_E1_nE2.get(xInd, yInd)));
			possible.add(new int[] {xInd+1, 0},  nanAsZero(prob_nE1_E2.get(xInd, yInd)));
			possible.add(new int[] {0, 0},  nanAsZero(prob_E1_E2.get(xInd, yInd)));
		} catch (RuntimeException e) {
			System.out.println("DEBUG: "+xInd+","+yInd+" occ="+occupancy.get(xInd, yInd));
			System.out.println("\t"+prob_nE1_nE2.get(xInd, yInd));
			System.out.println("\t"+prob_E1_nE2.get(xInd, yInd));
			System.out.println("\t"+prob_nE1_E2.get(xInd, yInd));
			System.out.println("\t"+prob_E1_E2.get(xInd, yInd));
			System.out.flush();
			throw e;
		}
		
		return possible;
	}
	
	private static double nanAsZero(double value) {
		if (Double.isNaN(value))
			return 0d;
		return value;
	}

	public EvenlyDiscrXYZ_DataSet getOccupancyXYZ() {
		return occupancy;
	}

	public EvenlyDiscrXYZ_DataSet getProb_E1_E2() {
		return prob_E1_E2;
	}

	public EvenlyDiscrXYZ_DataSet getProb_E1_nE2() {
		return prob_E1_nE2;
	}

	public EvenlyDiscrXYZ_DataSet getProb_nE1_E2() {
		return prob_nE1_E2;
	}

	public EvenlyDiscrXYZ_DataSet getProb_nE1_nE2() {
		return prob_nE1_nE2;
	}

}
