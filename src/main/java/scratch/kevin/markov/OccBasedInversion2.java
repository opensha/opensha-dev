package scratch.kevin.markov;

import java.io.IOException;
import java.util.List;
import java.util.Random;

import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.iden.RuptureIdentifier;

import scratch.kevin.simulators.MarkovChainBuilder;
import scratch.kevin.simulators.SimAnalysisCatLoader;
import scratch.kevin.simulators.SynchIdens;
import scratch.kevin.simulators.SynchIdens.SynchFaults;
import scratch.kevin.simulators.synch.StateSpacePlotter;

import com.google.common.primitives.Doubles;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public class OccBasedInversion2 {
	
	public static MarkovChain invert(EvenlyDiscrXYZ_DataSet occupancy, int numIters) {
		occupancy = fixOccupanyForMonoIncrease(occupancy);
		
		CPT cpt;
		try {
			cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, occupancy.getMaxZ());
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		
		// starting model is time independent
		EvenlyDiscrXYZ_DataSet prob_E1 = newSameSize(occupancy);
		EvenlyDiscrXYZ_DataSet prob_E2 = newSameSize(occupancy);
		
		EvenlyDiscretizedFunc occMarginalX = occupancy.calcMarginalXDist();
		EvenlyDiscretizedFunc occMarginalY = occupancy.calcMarginalYDist();
		
		for (int xInd=0; xInd<occupancy.getNumX(); xInd++) {
			double probX = probFromMarginal(occMarginalX, xInd);
			for (int yInd=0; yInd<occupancy.getNumY(); yInd++) {
				double probY = probFromMarginal(occMarginalY, yInd);
				
				prob_E1.set(xInd, yInd, probX);
				prob_E2.set(xInd, yInd, probY);
			}
		}
		
		int occSamples = 100000;
		
		EvenlyDiscrXYZ_DataSet curOcc = calcOcc(prob_E1, prob_E2, occupancy,
				occSamples, new int[] {0,0});
		try {
			StateSpacePlotter.plot2D(occupancy, cpt, "X", "Y", "Orig Occupancy", "", true, null);
			
			StateSpacePlotter.plot2D(curOcc, cpt, "X", "Y", "Start Invert Occupancy", "", true, null);
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		double energy = energy(occupancy, curOcc);
		
		int numKept = 0;
		
		for (int i=0; i<numIters; i++) {
			EvenlyDiscrXYZ_DataSet new_prob_E1 = prob_E1.copy();
			EvenlyDiscrXYZ_DataSet new_prob_E2 = prob_E2.copy();
			perturb(new_prob_E1, new_prob_E2);
			
			EvenlyDiscrXYZ_DataSet newOcc = calcOcc(new_prob_E1, new_prob_E2, occupancy,
					occSamples, new int[] {0,0});
			
			double newEnergy = energy(occupancy, newOcc);
			double T = 1 / (i+1);
			double P;
			if (newEnergy < energy)
				P = 1d;
			else
				P = Math.exp((energy - newEnergy) /  T); 
			
			if (P > r.nextDouble()) {
				// keep it
				prob_E1 = new_prob_E1;
				prob_E2 = new_prob_E2;
				energy = newEnergy;
				curOcc = newOcc;
				numKept++;
//				try {
//					StateSpacePlotter.plot2D(curOcc, cpt, "X", "Y", "Iter "+i+" Occupancy", "", true, null);
//				} catch (IOException e) {
//					throw ExceptionUtils.asRuntimeException(e);
//				}
			}
			if (i % 100 == 0)
				System.out.println("Iteration "+i+". E="+energy+", "+numKept+" kept");
		}
		
		try {
			StateSpacePlotter.plot2D(occupancy, cpt, "X", "Y", "Orig Occupancy", "", true, null);
			
			StateSpacePlotter.plot2D(curOcc, cpt, "X", "Y", "New Occupancy", "", true, null);
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		
		// TODO build chain
		return null;
	}
	
	private static EvenlyDiscrXYZ_DataSet fixOccupanyForMonoIncrease(EvenlyDiscrXYZ_DataSet occ) {
		occ = occ.copy();
		
		// list of all cells along the top/right side
		List<int[]> startCells = Lists.newArrayList();
		for (int xInd=0; xInd<occ.getNumX(); xInd++)
			startCells.add(new int[] {xInd, occ.getNumY()-1});
		for (int yInd=0; yInd<occ.getNumY()-1; yInd++)
			startCells.add(new int[] {occ.getNumX()-1, yInd});
		
		for (int[] startCell : startCells) {
			int xInd = startCell[0];
			int yInd = startCell[1];
			
			
			while (--xInd >= 0 && --yInd >=0) {
				double prevVal = occ.get(xInd+1, yInd+1);
				double curVal = occ.get(xInd, yInd);
				if (curVal < prevVal)
					occ.set(xInd, yInd, prevVal);
			}
		}
		
		return occ;
	}
	
	private static void assertValidProb(double prob) {
		Preconditions.checkState(Doubles.isFinite(prob));
		Preconditions.checkState(prob >= 0d && prob <= 1d, "bad prob: "+prob);
	}
	
	private static double probFromMarginal(EvenlyDiscretizedFunc marginal, int index) {
		if (index == marginal.size()-1)
			return 1d;
		double occAt = marginal.getY(index);
		double occNext = marginal.getY(index+1);
		double prob = (occAt - occNext)/occAt;
		if (prob < 0)
			prob = 0; // marginals can  actually be barely bad
		Preconditions.checkState(prob >= 0 && prob <= 1,
				"Bad prob: "+prob+" = ("+occAt+" - "+occNext+") / "+occAt);
		return prob;
	}
	
	private static EvenlyDiscrXYZ_DataSet calcOcc(
			EvenlyDiscrXYZ_DataSet prob_E1, EvenlyDiscrXYZ_DataSet prob_E2,
			EvenlyDiscrXYZ_DataSet origOcc, int samples, int[] initialState) {
		EvenlyDiscrXYZ_DataSet occ = newSameSize(prob_E1);
		
		int[] curState = initialState;
		
		for (int i=0; i<samples; i++) {
			if (curState[0] >= prob_E1.getNumX()-1)
				curState[0] = 0;
			if (curState[1] >= prob_E1.getNumY()-1)
				curState[1] = 0;
			occ.set(curState[0], curState[1], 1d+occ.get(curState[0], curState[1]));
			double rand = r.nextDouble();
			
			double prob_nE1_nE2 =
					origOcc.get(curState[0]+1, curState[1]+1)/origOcc.get(curState[0], curState[1]);
			assertValidProb(prob_nE1_nE2);
			double cmlProb = prob_nE1_nE2;
			if (rand < cmlProb) {
				curState = new int[] { curState[0]+1, curState[1]+1};
				continue;
			}
			
			double pE1_NE2givenS1S2 = (1d-prob_E2.get(curState[0], curState[1])) - prob_nE1_nE2;
			assertValidProb(pE1_NE2givenS1S2);
			cmlProb += pE1_NE2givenS1S2;
			if (rand < cmlProb) {
				curState = new int[] { 0, curState[1]+1};
				continue;
			}
			
			double pNE1_E2givenS1S2 = (1d-prob_E1.get(curState[0], curState[1])) - prob_nE1_nE2;
			assertValidProb(pNE1_E2givenS1S2);
			cmlProb += pNE1_E2givenS1S2;
			if (rand < cmlProb) {
				curState = new int[] { curState[0]+1, 0};
				continue;
			}
			
			curState = new int[] { 0, 0 };
		}
		occ.scale(1d/occ.getSumZ());
		return occ;
	}
	
	private static double energy(EvenlyDiscrXYZ_DataSet target, EvenlyDiscrXYZ_DataSet synthetic) {
		double e = 0;
		for (int index=0; index<target.size(); index++)
			e += Math.pow(100d*target.get(index) - 100d*synthetic.get(index), 2d);
		return e;
	}
	
	private static void perturb(EvenlyDiscrXYZ_DataSet prob_E1, EvenlyDiscrXYZ_DataSet prob_E2) {
		
		int xInd = r.nextInt(prob_E1.getNumX());
		int yInd = r.nextInt(prob_E1.getNumY());
		
		double factor;
		if (r.nextBoolean())
			factor = 1.1;
		else
			factor = 0.9;
		
		if (r.nextBoolean())
			prob_E1.set(xInd, yInd, prob_E1.get(xInd, yInd)*factor);
		else
			prob_E2.set(xInd, yInd, prob_E2.get(xInd, yInd)*factor);
	}
	
	private static EvenlyDiscrXYZ_DataSet newSameSize(EvenlyDiscrXYZ_DataSet o) {
		return new EvenlyDiscrXYZ_DataSet(o.getNumX(), o.getNumY(), o.getMinX(), o.getMinY(),
				o.getGridSpacingX(), o.getGridSpacingY());
	}
	
	private static Random r = new Random();

	public static void main(String[] args) throws IOException {
		double distSpacing = 10d;
		List<RuptureIdentifier> rupIdens = SynchIdens.getIndividualFaults(7d, 10d,
//				SynchFaults.SAF_MOJAVE, SynchFaults.SAF_COACHELLA, SynchFaults.SAN_JACINTO, SynchFaults.SAF_CARRIZO);
//				SynchFaults.SAF_MOJAVE, SynchFaults.SAF_COACHELLA, SynchFaults.SAN_JACINTO);
				SynchFaults.SAF_MOJAVE, SynchFaults.SAF_COACHELLA);
//				SynchFaults.SAF_CARRIZO, SynchFaults.SAF_CHOLAME);
//				SynchFaults.SAF_CARRIZO, SynchFaults.SAF_COACHELLA);
//				SynchFaults.SAF_COACHELLA, SynchFaults.SAN_JACINTO);
//				SynchFaults.SAF_COACHELLA, SynchFaults.SAF_CHOLAME);
//				SynchFaults.SAF_MOJAVE, SynchFaults.GARLOCK_WEST);
		
		List<? extends SimulatorEvent> events = new SimAnalysisCatLoader(true, rupIdens, true).getEvents();
		
		List<int[]> fullPath = MarkovChainBuilder.getStatesPath(distSpacing, events, rupIdens, 0d);
		
		EmpiricalMarkovChain origChain = new EmpiricalMarkovChain(fullPath, distSpacing);
		
		EvenlyDiscrXYZ_DataSet occupancy =
				new StateSpacePlotter(origChain, rupIdens, null).getOccupancy(0, 1);
		
		invert(occupancy, 10000);
	}

}
