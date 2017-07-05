package scratch.kevin.simulators.synch;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotWindow;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.iden.RuptureIdentifier;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

import scratch.kevin.markov.EmpiricalMarkovChain;
import scratch.kevin.markov.OccupancyBasedMarkovChain2D;
import scratch.kevin.markov.PossibleStates;
import scratch.kevin.simulators.MarkovChainBuilder;
import scratch.kevin.simulators.PeriodicityPlotter;
import scratch.kevin.simulators.SimAnalysisCatLoader;
import scratch.kevin.simulators.SynchIdens;
import scratch.kevin.simulators.SynchIdens.SynchFaults;
import scratch.kevin.simulators.synch.StateSpacePlotter.MarkovProb;

public class TestOccFromSynch {
	
	public static void main(String[] args) throws IOException {
		double minMag = 7d; 
		double maxMag = 10d;
		double distSpacing = 10d; // years
		
		boolean do2DPlots = true;
		File predictDir = new File("/home/kevin/Simulators/predict");
		File plot2DOutputDir = new File(predictDir, "plots_2d");
		File synchPlotOutputDir = new File(predictDir, "synch_plots");
		
		List<RuptureIdentifier> rupIdens = SynchIdens.getIndividualFaults(minMag, maxMag,
//				SynchFaults.SAF_MOJAVE, SynchFaults.SAF_COACHELLA, SynchFaults.SAN_JACINTO, SynchFaults.SAF_CARRIZO);
//				SynchFaults.SAF_MOJAVE, SynchFaults.SAF_COACHELLA, SynchFaults.SAN_JACINTO);
//				SynchFaults.SAF_MOJAVE, SynchFaults.SAF_COACHELLA);
//				SynchFaults.SAF_CARRIZO, SynchFaults.SAF_CHOLAME);
//				SynchFaults.SAF_CARRIZO, SynchFaults.SAF_COACHELLA);
				SynchFaults.SAF_COACHELLA, SynchFaults.SAN_JACINTO);
//				SynchFaults.SAF_COACHELLA, SynchFaults.SAF_CHOLAME);
//				SynchFaults.SAF_MOJAVE, SynchFaults.GARLOCK_WEST);
		String name1 = rupIdens.get(0).getName();
		String name2 = rupIdens.get(1).getName();
		
		List<? extends SimulatorEvent> events = new SimAnalysisCatLoader(true, rupIdens, true).getEvents();
		
		List<int[]> fullPath = MarkovChainBuilder.getStatesPath(distSpacing, events, rupIdens, 0d);
		
		EmpiricalMarkovChain origChain = new EmpiricalMarkovChain(fullPath, distSpacing);
		PossibleStates origOccupancy = origChain.getOccupancy();
		
		EvenlyDiscrXYZ_DataSet origOccXYZ = new EvenlyDiscrXYZ_DataSet(
				50, 50, 0.5*distSpacing, 0.5*distSpacing, distSpacing);
//		EvenlyDiscrXYZ_DataSet newOccXYZ = new EvenlyDiscrXYZ_DataSet(
//				50, 50, 0.5*distSpacing, 0.5*distSpacing, distSpacing);
		
		for (int[] state : origOccupancy.getStates()) {
			int xInd = state[0];
			int yInd = state[1];
			
			if (xInd >= origOccXYZ.getNumX() || yInd >= origOccXYZ.getNumY())
				continue;
			
			double freq = origOccupancy.getFrequency(state);
			
			origOccXYZ.set(xInd, yInd, origOccXYZ.get(xInd, yInd)+freq);
		}
		origOccXYZ.scale(1d/origOccXYZ.getSumZ());
		
		EvenlyDiscrXYZ_DataSet indepOccXYZ = new EvenlyDiscrXYZ_DataSet(
				50, 50, 0.5*distSpacing, 0.5*distSpacing, distSpacing);
		
		
		EvenlyDiscretizedFunc marginalX = origOccXYZ.calcMarginalXDist();
		EvenlyDiscretizedFunc marginalY = origOccXYZ.calcMarginalYDist();
		
		// multiply marginals for starting
		for (int xInd=0; xInd<indepOccXYZ.getNumX(); xInd++)
			for (int yInd=0; yInd<indepOccXYZ.getNumY(); yInd++)
				indepOccXYZ.set(xInd, yInd, marginalX.getY(xInd)*marginalY.getY(yInd));
		
//		// fill in edges
//		for (int xInd=0; xInd<newOccXYZ.getNumX(); xInd++)
//			newOccXYZ.set(xInd, 0, origOccXYZ.get(xInd, 0));
//		for (int yInd=1; yInd<newOccXYZ.getNumY(); yInd++)
//			newOccXYZ.set(0, yInd, origOccXYZ.get(0, yInd));
		
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, origOccXYZ.getMaxZ());
		File outputDir = new File("/home/kevin/Simulators/inferred_occupancy");
		
//		newOccXYZ = getScaledToDiags(newOccXYZ, origOccXYZ);
		
		EvenlyDiscrXYZ_DataSet newOccXYZ = buildOcc(origOccXYZ.getRow(0), origOccXYZ.getCol(0),
				origOccXYZ.calcMarginalXDist(), origOccXYZ.calcMarginalYDist(), 100);
		
		// now try to fix marginals
//		Preconditions.checkState(validateMonotonical(origOccXYZ));
//		Preconditions.checkState(validateMonotonical(newOccXYZ));
//		for (int i=0; i<100 ; i++) {
////			StateSpacePlotter.plot2D(newOccXYZ, cpt, name1, name2, "Alg. Occupancy",
////					"occ_alg_it_"+i, true, outputDir);
//			newOccXYZ = getScaledToDiags(newOccXYZ, origOccXYZ);
//			EvenlyDiscrXYZ_DataSet scaleMargeX = getScaledToMarginals(newOccXYZ, marginalX, true);
//			EvenlyDiscrXYZ_DataSet scaleMargeY = getScaledToMarginals(newOccXYZ, marginalY, false);
//			newOccXYZ = blend(new double[] {0.5, 0.25, 0.25},
//					new EvenlyDiscrXYZ_DataSet[] {newOccXYZ, scaleMargeX, scaleMargeY});
//			newOccXYZ = getMonotonicallyDecreaseDiags(newOccXYZ);
//			newOccXYZ.scale(1d/newOccXYZ.getSumZ());
////			replaceAxis(newOccXYZ, origOccXYZ);
//		}
		System.out.println("Sum Z: orig="+origOccXYZ.getSumZ()+", new="+newOccXYZ.getSumZ());
		
//		double maxZ = Math.max(newOccXYZ.getMaxZ(), origOccXYZ.getMaxZ());
		boolean marginals = false;
		StateSpacePlotter.plot2D(origOccXYZ, cpt, name1, name2, "Orig. Occupancy",
				"occ_orig", marginals, outputDir);
		StateSpacePlotter.plot2D(indepOccXYZ, cpt, name1, name2, "Indep. Occupancy",
				"occ_indep", marginals, outputDir);
		StateSpacePlotter.plot2D(newOccXYZ, cpt, name1, name2, "Alg. Occupancy",
				"occ_alg", marginals, outputDir);
		
//		File origOutput = new File(outputDir, "markov_probs_orig");
//		Preconditions.checkState((origOutput.exists() && origOutput.isDirectory()) || origOutput.mkdir());
//		File newOutput = new File(outputDir, "markov_probs_alg");
//		Preconditions.checkState((newOutput.exists() && newOutput.isDirectory()) || newOutput.mkdir());
//		StateSpacePlotter origPlotter = new StateSpacePlotter(origChain,
//				rupIdens, origOutput);
//		newOccXYZ.scale(fullPath.size());
////		StateSpacePlotter newPlotter = new StateSpacePlotter(new OccupancyBasedMarkovChain2D(distSpacing, newOccXYZ),
////				rupIdens, newOutput);
//		StateSpacePlotter newPlotter = new StateSpacePlotter(new OccupancyBasedMarkovChain2D(
//				distSpacing, origChain.getOccupancy()), rupIdens, newOutput);
//		
//		origPlotter.plotOccupancies();
//		newPlotter.plotOccupancies();
//		
//		for (MarkovProb prob : MarkovProb.values()) {
//			origPlotter.plotProb(prob);
//			newPlotter.plotProb(prob);
//		}
	}
	
	public static EvenlyDiscrXYZ_DataSet buildOcc(
			EvenlyDiscretizedFunc firstRow, EvenlyDiscretizedFunc firstCol,
			EvenlyDiscretizedFunc marginalX, EvenlyDiscretizedFunc marginalY, int numIters) {
		// validate inputs
		for (int i=0; i<firstRow.size(); i++)
			Preconditions.checkState(Doubles.isFinite(firstRow.getY(i)));
		for (int i=0; i<firstCol.size(); i++)
			Preconditions.checkState(Doubles.isFinite(firstCol.getY(i)));
		for (int i=0; i<marginalX.size(); i++)
			Preconditions.checkState(Doubles.isFinite(marginalX.getY(i)));
		for (int i=0; i<marginalY.size(); i++)
			Preconditions.checkState(Doubles.isFinite(marginalY.getY(i)));
		Preconditions.checkArgument(firstRow.size() == marginalX.size());
		Preconditions.checkArgument(firstCol.size() == marginalY.size());
		
		// first calculate assuming independence
		EvenlyDiscrXYZ_DataSet newOccXYZ = new EvenlyDiscrXYZ_DataSet(
				firstRow.size(), firstCol.size(), firstRow.getMinX(), firstCol.getMinX(),
				firstRow.getDelta(), firstCol.getDelta());
		for (int xInd=0; xInd<newOccXYZ.getNumX(); xInd++)
			for (int yInd=0; yInd<newOccXYZ.getNumY(); yInd++)
				newOccXYZ.set(xInd, yInd, marginalX.getY(xInd)*marginalY.getY(yInd));
		
//		checkFinite(newOccXYZ);
		
		for (int i=0; i<numIters ; i++) {
//			StateSpacePlotter.plot2D(newOccXYZ, cpt, name1, name2, "Alg. Occupancy",
//					"occ_alg_it_"+i, true, outputDir);
//			newOccXYZ = getScaledToDiags(newOccXYZ, origOccXYZ);
			newOccXYZ = getScaledToDiags(newOccXYZ, firstRow, firstCol);
//			checkFinite(newOccXYZ);
			EvenlyDiscrXYZ_DataSet scaleMargeX = getScaledToMarginals(newOccXYZ, marginalX, true);
//			checkFinite(scaleMargeX);
			EvenlyDiscrXYZ_DataSet scaleMargeY = getScaledToMarginals(newOccXYZ, marginalY, false);
//			checkFinite(scaleMargeY);
			newOccXYZ = blend(new double[] {0.5, 0.25, 0.25},
					new EvenlyDiscrXYZ_DataSet[] {newOccXYZ, scaleMargeX, scaleMargeY});
//			checkFinite(newOccXYZ);
			newOccXYZ = getMonotonicallyDecreaseDiags(newOccXYZ);
//			checkFinite(newOccXYZ);
			newOccXYZ.scale(1d/newOccXYZ.getSumZ());
//			checkFinite(newOccXYZ);
//			replaceAxis(newOccXYZ, origOccXYZ);
		}
		
		checkFinite(newOccXYZ);
		
		return newOccXYZ;
	}
	
	private static void checkFinite(EvenlyDiscrXYZ_DataSet xyz) {
		for (int xInd=0; xInd<xyz.getNumX(); xInd++)
			for (int yInd=0; yInd<xyz.getNumY(); yInd++)
				Preconditions.checkState(Doubles.isFinite(xyz.get(xInd, yInd)));
	}
	
	private static void replaceAxis(EvenlyDiscrXYZ_DataSet target, EvenlyDiscrXYZ_DataSet source) {
		for (int xInd=0; xInd<target.getNumX(); xInd++)
			target.set(xInd, 0, source.get(xInd, 0));
		for (int yInd=0; yInd<target.getNumY(); yInd++)
			target.set(0, yInd, source.get(0, yInd));
	}
	
//	private static void scaleDiags1(EvenlyDiscrXYZ_DataSet target, EvenlyDiscrXYZ_DataSet source) {
//		for (int xInd=0; xInd<target.getNumX(); xInd++) {
//			int myYInd = 0;
//			int myXInd = xInd;
//			double startVal = source.get(myXInd, myYInd);
//			
//			while (myXInd < source.getNumX() && myYInd < source.getNumY()) {
//				target.set(myXInd, myYInd, target.get(myXInd, myYInd)*startVal);
//				myXInd++;
//				myYInd++;
//			}
//		}
//		
//		for (int yInd=1; yInd<target.getNumY(); yInd++) {
//			int myXInd = 0;
//			int myYInd = yInd;
//			double startVal = source.get(myXInd, myYInd);
//			
//			while (myXInd < source.getNumX() && myYInd < source.getNumY()) {
//				target.set(myXInd, myYInd, target.get(myXInd, myYInd)*startVal);
//				myXInd++;
//				myYInd++;
//			}
//		}
//		
//		double newSumZ = 0d;
//		for (int index=0; index<target.size(); index++)
//			newSumZ += target.get(index);
//		target.scale(1d/newSumZ);
//	}
	
	private static EvenlyDiscrXYZ_DataSet getScaledToDiags(
			EvenlyDiscrXYZ_DataSet current, EvenlyDiscrXYZ_DataSet source) {
		return getScaledToDiags(current, source.getRow(0), source.getCol(0));
	}
	
	private static EvenlyDiscrXYZ_DataSet getScaledToDiags(EvenlyDiscrXYZ_DataSet current,
			EvenlyDiscretizedFunc firstRow, EvenlyDiscretizedFunc firstCol) {
		current = current.copy();
		for (int xInd=0; xInd<current.getNumX(); xInd++) {
			int myYInd = 0;
			int myXInd = xInd;
			double targetStartVal = firstRow.getY(myXInd);
			double curStartVal = current.get(myXInd, myYInd);
			double scalar = targetStartVal/curStartVal;
			if (!Doubles.isFinite(scalar))
				scalar = 0;
			
			while (myXInd < current.getNumX() && myYInd < current.getNumY()) {
				current.set(myXInd, myYInd, current.get(myXInd, myYInd)*scalar);
				myXInd++;
				myYInd++;
			}
		}
		
		for (int yInd=1; yInd<current.getNumY(); yInd++) {
			int myXInd = 0;
			int myYInd = yInd;
			double targetStartVal = firstCol.getY(myYInd);
			double curStartVal = current.get(myXInd, myYInd);
			double scalar = targetStartVal/curStartVal;
			if (!Doubles.isFinite(scalar))
				scalar = 0;
			
			while (myXInd < current.getNumX() && myYInd < current.getNumY()) {
				current.set(myXInd, myYInd, current.get(myXInd, myYInd)*scalar);
				myXInd++;
				myYInd++;
			}
		}
		return current;
	}
	
	private static EvenlyDiscrXYZ_DataSet getScaledToMarginals(
			EvenlyDiscrXYZ_DataSet xyz, EvenlyDiscretizedFunc marginal, boolean x) {
		xyz = xyz.copy();
		int maxI, maxJ;
		if (x) {
			maxI = xyz.getNumX();
			maxJ = xyz.getNumY();
		} else {
			maxI = xyz.getNumY();
			maxJ = xyz.getNumX();
		}
		
		EvenlyDiscretizedFunc startMarginal;
		if (x)
			startMarginal = xyz.calcMarginalXDist();
		else
			startMarginal = xyz.calcMarginalYDist();
		
		for (int i=1; i<maxI; i++) {
			double targetSum = marginal.getY(i);
			double startSum = startMarginal.getY(i);
			double scalar = targetSum/startSum;
			if (startSum == 0d)
				scalar = 0d;
			for (int j=1; j<maxJ; j++) {
				int xInd, yInd;
				if (x) {
					xInd = i;
					yInd = j;
				} else {
					xInd = j;
					yInd = i;
				}
				xyz.set(xInd, yInd, xyz.get(xInd, yInd)*scalar);
			}
		}
		return xyz;
	}
	
	private static EvenlyDiscrXYZ_DataSet blend(double[] fracts, EvenlyDiscrXYZ_DataSet[] xyzs) {
		Preconditions.checkState(fracts.length == xyzs.length);
		EvenlyDiscrXYZ_DataSet sum = xyzs[0].copy();
		for (int xInd=0; xInd<sum.getNumX(); xInd++) {
			for (int yInd=0; yInd<sum.getNumY(); yInd++) {
				double val = 0d;
				for (int i=0; i<xyzs.length; i++)
					val += fracts[i]*xyzs[i].get(xInd, yInd);
				sum.set(xInd, yInd, val);
			}
		}
		return sum;
	}
	
	private static EvenlyDiscrXYZ_DataSet getMonotonicallyDecreaseDiags(EvenlyDiscrXYZ_DataSet xyz) {
		xyz = xyz.copy();
		
		for (int x=1; x<xyz.getNumX(); x++) {
			for (int y=1; y<xyz.getNumY(); y++) {
				double myZ = xyz.get(x, y);
				double prevZ = xyz.get(x-1, y-1);
				if (myZ > prevZ) {
					double scalar = prevZ/myZ;
					int xInd = x;
					int yInd = y;
					while (xInd < xyz.getNumX() && yInd < xyz.getNumY()) {
						xyz.set(xInd, yInd, xyz.get(xInd, yInd)*scalar);
						xInd++;
						yInd++;
					}
				}
			}
		}
		
		return xyz;
	}
	
//	private static boolean validateMonotonical(EvenlyDiscrXYZ_DataSet xyz) {
//		if (!isMonotonicallyDecreasing(xyz.calcMarginalXDist()))
//			return false;
//		System.out.println("Passed X marg");
//		if (!isMonotonicallyDecreasing(xyz.calcMarginalYDist()))
//			return false;
//		System.out.println("Passed Y marg");
//		for (int xInd=0; xInd<xyz.getNumX(); xInd++)
//			if (!isMonotonicallyDecreasing(xyz.getDiag(xInd, 0)))
//				return false;
//		System.out.println("Passed X diag");
//		for (int yInd=0; yInd<xyz.getNumY(); yInd++)
//			if (!isMonotonicallyDecreasing(xyz.getDiag(0, yInd)))
//				return false;
//		System.out.println("Passed Y diag");
//		return true;
//	}
//	
//	private static boolean isMonotonicallyDecreasing(DiscretizedFunc func) {
//		for (int i=1; i<func.getNum(); i++) {
//			if ((float)func.getY(i) > (float)func.getY(i-1)) {
//				System.out.println("i="+i+", y("+i+")="+func.getY(i)+", y("+(i-1)+")="+func.getY(i-1));
//				return false;
//			}
//		}
//		return true;
//	}

}
