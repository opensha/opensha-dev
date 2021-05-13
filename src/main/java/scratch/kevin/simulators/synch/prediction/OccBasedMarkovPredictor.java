package scratch.kevin.simulators.synch.prediction;

import java.util.List;

import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;

import com.google.common.base.Preconditions;

import scratch.kevin.markov.OccBasedIterativeSolver;
import scratch.kevin.markov.PossibleStates;
import scratch.kevin.simulators.synch.StateSpacePlotter;
import scratch.kevin.simulators.synch.TestOccFromSynch;

public class OccBasedMarkovPredictor extends MarkovPredictor {
	
	private boolean inferOcc;
	private double fractOfOccToKeep;
	
	public OccBasedMarkovPredictor(boolean inferOcc, double fractOfOccToKeep) {
		this(null, inferOcc, fractOfOccToKeep);
	}
	
	public OccBasedMarkovPredictor(Predictor backupPredictor, boolean inferOcc, double fractOfOccToKeep) {
		super(backupPredictor);
		this.inferOcc = inferOcc;
		this.fractOfOccToKeep = fractOfOccToKeep;
	}

	@Override
	public void init(List<int[]> path, double distSpacing) {
		super.init(path, distSpacing);
		// validate dimensions
		Preconditions.checkState(chain.getNDims() == 2, "Must be 2-dimensional for occupancy based");
		
		// calculate size of XYZ grids
		int numX = StateSpacePlotter.getNBinsForFract(chain, 0, fractOfOccToKeep);
		int numY = StateSpacePlotter.getNBinsForFract(chain, 1, fractOfOccToKeep);
		
		// calculate actual occupancy
		double binStart = 0.5*distSpacing;
		EvenlyDiscrXYZ_DataSet occXYZ = new EvenlyDiscrXYZ_DataSet(numX, numY, binStart, binStart, distSpacing);
		
		PossibleStates occupancy = chain.getOccupancy();
		double occScale = (double)path.size()/occupancy.getTot();
		for (int[] state : occupancy.getStates()) {
			int xInd = state[0];
			int yInd = state[1];
			
			if (xInd >= occXYZ.getNumX() || yInd >= occXYZ.getNumY())
				continue;
			
			double freq = occupancy.getFrequency(state)*occScale;
			
			occXYZ.set(xInd, yInd, occXYZ.get(xInd, yInd)+freq);
		}
		
		// get inferred occupancy, if applicable
		if (inferOcc) {
			EvenlyDiscrXYZ_DataSet newOccXYZ = TestOccFromSynch.buildOcc(occXYZ.getRow(0), occXYZ.getCol(0),
					occXYZ.calcMarginalXDist(), occXYZ.calcMarginalYDist(), 100);
			newOccXYZ.scale(occXYZ.getSumZ()/newOccXYZ.getSumZ());
			occXYZ = newOccXYZ;
		}
		
		chain = OccBasedIterativeSolver.calc(occXYZ, 0);
	}

	@Override
	public String getName() {
		String name = "";
		if (inferOcc)
			name += "Inferred ";
		name += "Occ Based "+super.getName();
		return name;
	}

}
