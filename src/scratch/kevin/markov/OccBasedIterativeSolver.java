package scratch.kevin.markov;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.Named;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.iden.RuptureIdentifier;

import scratch.kevin.simulators.MarkovChainBuilder;
import scratch.kevin.simulators.PeriodicityPlotter;
import scratch.kevin.simulators.SimAnalysisCatLoader;
import scratch.kevin.simulators.SynchIdens;
import scratch.kevin.simulators.SynchIdens.SynchFaults;
import scratch.kevin.simulators.synch.StateSpacePlotter;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.primitives.Doubles;

public class OccBasedIterativeSolver {
	
	public static MarkovChain calc(EvenlyDiscrXYZ_DataSet occupancy, int iterations) {
		CPT occCPT;
		try {
			occCPT = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, occupancy.getMaxZ());
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		CPT probCPT = occCPT.rescale(0d, 1d);
		
		// starting model is independent
		EvenlyDiscrXYZ_DataSet prob_E1_E2 = newSameSize(occupancy);
		EvenlyDiscrXYZ_DataSet prob_E1_nE2 = newSameSize(occupancy);
		EvenlyDiscrXYZ_DataSet prob_nE1_E2 = newSameSize(occupancy);
		EvenlyDiscrXYZ_DataSet prob_nE1_nE2 = calcProbNone(occupancy);
		
		calcAssumingIndependence(occupancy, prob_E1_E2, prob_E1_nE2, prob_nE1_E2, prob_nE1_nE2);
		
		int occSamples = 100000;
		
		EvenlyDiscrXYZ_DataSet curOcc = calcOcc(occupancy, prob_E1_E2,
				prob_E1_nE2, prob_nE1_E2, prob_nE1_nE2, occSamples, new int[] {0,0});
		
		System.out.println("Doing "+iterations+" iterations");
		for (int i=0; i<iterations; i++) {
			EvenlyDiscrXYZ_DataSet prob_E1_E2_bothFix = prob_E1_E2.copy();
			EvenlyDiscrXYZ_DataSet prob_E1_nE2_bothFix = prob_E1_nE2.copy();
			EvenlyDiscrXYZ_DataSet prob_nE1_E2_bothFix = prob_nE1_E2.copy();
			EvenlyDiscrXYZ_DataSet prob_nE1_nE2_bothFix = prob_nE1_nE2.copy();
			// now adjust such that the occupancy probability of [0,0] will be correct
			adjustForTotProbBoth(occupancy, curOcc, prob_E1_E2_bothFix, prob_E1_nE2_bothFix,
					prob_nE1_E2_bothFix, prob_nE1_nE2_bothFix);
			
			// adjust to fit occ x axis
			EvenlyDiscrXYZ_DataSet prob_E1_E2_xFix = prob_E1_E2.copy();
			EvenlyDiscrXYZ_DataSet prob_E1_nE2_xFix = prob_E1_nE2.copy();
			EvenlyDiscrXYZ_DataSet prob_nE1_E2_xFix = prob_nE1_E2.copy();
			EvenlyDiscrXYZ_DataSet prob_nE1_nE2_xFix = prob_nE1_nE2.copy();
			adjustForEdge(true, occupancy, curOcc, prob_E1_E2_xFix, prob_E1_nE2_xFix,
					prob_nE1_E2_xFix, prob_nE1_nE2_xFix);
			
			// adjust to fit occ y axis
			EvenlyDiscrXYZ_DataSet prob_E1_E2_yFix = prob_E1_E2.copy();
			EvenlyDiscrXYZ_DataSet prob_E1_nE2_yFix = prob_E1_nE2.copy();
			EvenlyDiscrXYZ_DataSet prob_nE1_E2_yFix = prob_nE1_E2.copy();
			EvenlyDiscrXYZ_DataSet prob_nE1_nE2_yFix = prob_nE1_nE2.copy();
			adjustForEdge(false, occupancy, curOcc, prob_E1_E2_yFix, prob_E1_nE2_yFix,
					prob_nE1_E2_yFix, prob_nE1_nE2_yFix);
			
			double[] blendFracts = {0.3, 0.5, 0.1, 0.1};
//			double[] blendFracts = {0.25, 0.45, 0.15, 0.15};
//			double[] blendFracts = {0.1, 0.5, 0.2, 0.2};
//			double[] blendFracts = {0.1, 0.4, 0.25, 0.25};
			
			// now blend everything
			prob_E1_E2 = blend(blendFracts, new EvenlyDiscrXYZ_DataSet[] {
					prob_E1_E2, prob_E1_E2_bothFix, prob_E1_E2_xFix, prob_E1_E2_yFix});
			prob_E1_nE2 = blend(blendFracts, new EvenlyDiscrXYZ_DataSet[] {
					prob_E1_nE2, prob_E1_nE2_bothFix, prob_E1_nE2_xFix, prob_E1_nE2_yFix});
			prob_nE1_E2 = blend(blendFracts, new EvenlyDiscrXYZ_DataSet[] {
					prob_nE1_E2, prob_nE1_E2_bothFix, prob_nE1_E2_xFix, prob_nE1_E2_yFix});
			prob_nE1_nE2 = blend(blendFracts, new EvenlyDiscrXYZ_DataSet[] {
					prob_nE1_nE2, prob_nE1_nE2_bothFix, prob_nE1_nE2_xFix, prob_nE1_nE2_yFix});
			
			curOcc = calcOcc(occupancy,
					prob_E1_E2, prob_E1_nE2, prob_nE1_E2, prob_nE1_nE2, occSamples, new int[] {0,0});
			curOcc.scale(occupancy.getSumZ()/curOcc.getSumZ());
			
			assertValidProbs(occupancy, prob_E1_E2, prob_E1_nE2, prob_nE1_E2, prob_nE1_nE2);
		}
		if (iterations <= 0) {
			// just adjust such that the occupancy probability of [0,0] will be correct
			adjustForTotProbBoth(occupancy, curOcc, prob_E1_E2, prob_E1_nE2, prob_nE1_E2, prob_nE1_nE2);
			
			curOcc = calcOcc(occupancy,
					prob_E1_E2, prob_E1_nE2, prob_nE1_E2, prob_nE1_nE2, occSamples, new int[] {0,0});
			curOcc.scale(occupancy.getSumZ()/curOcc.getSumZ());
		}
		
		try {
			StateSpacePlotter.plot2D(occupancy, occCPT, "X", "Y", "Orig Occupancy", "", true, null);
			
			StateSpacePlotter.plot2D(curOcc, occCPT, "X", "Y", "New Occupancy", "", true, null);
			
//			StateSpacePlotter.plot2D(prob_E1_E2, probCPT, "X", "Y", "P(E1 E2)", "", false, null);
//			StateSpacePlotter.plot2D(prob_E1_nE2, probCPT, "X", "Y", "P(E1 !E2)", "", false, null);
//			StateSpacePlotter.plot2D(prob_nE1_E2, probCPT, "X", "Y", "P(!E1 E2)", "", false, null);
//			StateSpacePlotter.plot2D(prob_nE1_nE2, probCPT, "X", "Y", "P(!E1 !E2)", "", false, null);
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		
		XYZBasedMarkovChain xyzChain = new XYZBasedMarkovChain(curOcc, prob_E1_E2, prob_E1_nE2, prob_nE1_E2, prob_nE1_nE2);
		
		return xyzChain;
	}
	
	private static EvenlyDiscrXYZ_DataSet newSameSize(EvenlyDiscrXYZ_DataSet o) {
		return new EvenlyDiscrXYZ_DataSet(o.getNumX(), o.getNumY(), o.getMinX(), o.getMinY(),
				o.getGridSpacingX(), o.getGridSpacingY());
	}
	
	private static EvenlyDiscrXYZ_DataSet calcProbNone(EvenlyDiscrXYZ_DataSet occupancy) {
		EvenlyDiscrXYZ_DataSet pNone = newSameSize(occupancy);
		
		for (int xInd=0; xInd<occupancy.getNumX(); xInd++) {
			for (int yInd=0; yInd<occupancy.getNumY(); yInd++) {
				double occ = occupancy.get(xInd, yInd);
				if (occ == 0d) {
					pNone.set(xInd, yInd, Double.NaN);
					continue;
				}
				double occNext;
				if (xInd < occupancy.getNumX()-1 && yInd < occupancy.getNumY()-1)
					occNext = occupancy.get(xInd+1, yInd+1);
				else
					occNext = 0;
				double prob = occNext/occ;
				if (prob < 0d)
					prob = 0d;
				assertValidProb(prob);
				pNone.set(xInd, yInd, prob);
			}
		}
		
		return pNone;
	}
	
	private static void calcAssumingIndependence(EvenlyDiscrXYZ_DataSet occupancy,
			EvenlyDiscrXYZ_DataSet prob_E1_E2, EvenlyDiscrXYZ_DataSet prob_E1_nE2,
			EvenlyDiscrXYZ_DataSet prob_nE1_E2, EvenlyDiscrXYZ_DataSet prob_nE1_nE2) {
		EvenlyDiscretizedFunc occMarginalX = occupancy.calcMarginalXDist();
		EvenlyDiscretizedFunc occMarginalY = occupancy.calcMarginalYDist();
		
		System.out.println("Calculating assuming independence");
		
		// calculate assuming independence
		for (int xInd=0; xInd<occupancy.getNumX(); xInd++) {
			for (int yInd=0; yInd<occupancy.getNumY(); yInd++) {
				if (occupancy.get(xInd, yInd) == 0d) {
					prob_E1_E2.set(xInd, yInd, Double.NaN);
					prob_E1_nE2.set(xInd, yInd, Double.NaN);
					prob_nE1_E2.set(xInd, yInd, Double.NaN);
					continue;
				}
				double pNone = prob_nE1_nE2.get(xInd, yInd);
				double probEither = 1d-pNone;
				if ((float)probEither == 0f) {
					prob_E1_E2.set(xInd, yInd, 0d);
					prob_E1_nE2.set(xInd, yInd, 0d);
					prob_nE1_E2.set(xInd, yInd, 0d);
				} else {
					double margProbE1 = calcMarginalProb(occMarginalX, xInd);
					double margProbE2 = calcMarginalProb(occMarginalY, yInd);
					
					if ((float)margProbE1 == 0f && (float)margProbE2 == 0f) {
						// we can get here with bad occupancy distributions
						// just calculate assuming marginals eaqual, so that
						// 			margProbE1 + margProbE1 + margProbE1*margProbE2 = probEither
						margProbE1 = Math.sqrt(probEither+1d) - 1d;
						margProbE2 = margProbE1;
					}
					
					double indep_E1_E2 = margProbE1*margProbE2;
					assertValidProb(indep_E1_E2);
//					double indep_E1_nE2 = margProbE1 - indep_E1_E2;
					double indep_E1_nE2 = margProbE1 * (1d - margProbE2);
					assertValidProb(indep_E1_nE2);
//					double indep_nE1_E2 = margProbE2 - indep_E1_E2;
					double indep_nE1_E2 = margProbE2 * (1d - margProbE1);
					assertValidProb(indep_nE1_E2);
					
					double indep_total_either = indep_E1_E2 + indep_E1_nE2 + indep_nE1_E2;
					Preconditions.checkState(indep_total_either > 0d,
							"Bad indep tot: "+indep_E1_E2+" + "+indep_E1_nE2+" + "+indep_nE1_E2+" = "+indep_total_either);
					double scalar = probEither/indep_total_either;
					Preconditions.checkState(Doubles.isFinite(scalar),
							"Bad scalar: "+probEither+"/"+indep_total_either+" = "+scalar);
					
					prob_E1_E2.set(xInd, yInd, indep_E1_E2*scalar);
					prob_E1_nE2.set(xInd, yInd, indep_E1_nE2*scalar);
					prob_nE1_E2.set(xInd, yInd, indep_nE1_E2*scalar);
				}
			}
		}
		
		assertValidProbs(occupancy, prob_E1_E2, prob_E1_nE2, prob_nE1_E2, prob_nE1_nE2);
	}
	
	private static double calcMarginalProb(EvenlyDiscretizedFunc occMarginal, int index) {
		double occ = occMarginal.getY(index);
		double occNext;
		if (index < occMarginal.size()-1)
			occNext = occMarginal.getY(index+1);
		else
			occNext = 0d;
		if (occNext > occ)
			occNext = occ;
		double prob = (occ - occNext)/occ;
		if (prob < 0d)
			prob = 0d;
		assertValidProb(prob);
		return prob;
	}
	
	private static void adjustForTotProbBoth(
			EvenlyDiscrXYZ_DataSet occupancy, EvenlyDiscrXYZ_DataSet curOcc,
			EvenlyDiscrXYZ_DataSet prob_E1_E2, EvenlyDiscrXYZ_DataSet prob_E1_nE2,
			EvenlyDiscrXYZ_DataSet prob_nE1_E2, EvenlyDiscrXYZ_DataSet prob_nE1_nE2) {
		double real_00 = occupancy.get(0, 0)/occupancy.getSumZ();
		double cur_00 = curOcc.get(0, 0)/curOcc.getSumZ();
		double both_scalar = real_00 / cur_00;
		for (int xInd=0; xInd<occupancy.getNumX(); xInd++) {
			for (int yInd=0; yInd<occupancy.getNumY(); yInd++) {
				if (occupancy.get(xInd, yInd) == 0d)
					continue;
				double p_E1_E2 = prob_E1_E2.get(xInd, yInd);
				if (p_E1_E2 == 0d)
					continue;
				double p_E1_nE2 = prob_E1_nE2.get(xInd, yInd);
				double p_nE1_E2 = prob_nE1_E2.get(xInd, yInd);
				double p_nE1_nE2 = prob_nE1_nE2.get(xInd, yInd);
				
				p_E1_E2 *= both_scalar;
				if (p_E1_E2 + p_nE1_nE2 > 1d)
					p_E1_E2 = 1d - p_nE1_nE2;
				assertValidProb(p_E1_E2);
				assertValidProb(p_E1_E2+p_nE1_nE2);
				
				double origSingles = p_E1_nE2 + p_nE1_E2;
				double leftForSingles = 1d - (p_E1_E2 + p_nE1_nE2);
				
				if (origSingles == 0d) {
//					System.out.println("Error case! ["+xInd+","+yInd+"]");
//					System.out.println("\tOcc: "+occupancy.get(xInd, yInd));
//					System.out.println("P(E1 E2) = "+p_E1_E2);
//					System.out.println("P(E1 !E2) = "+p_E1_nE2);
//					System.out.println("P(!E1 E2) = "+p_nE1_E2);
//					System.out.println("P(!E1 !E2) = "+p_nE1_nE2);
					
					p_E1_nE2 = 0d;
					p_nE1_E2 = 0d;
				} else {
					p_E1_nE2 = p_E1_nE2*leftForSingles/origSingles;
					p_nE1_E2 = p_nE1_E2*leftForSingles/origSingles;
				}
				
				prob_E1_E2.set(xInd, yInd, p_E1_E2);
				prob_E1_nE2.set(xInd, yInd, p_E1_nE2);
				prob_nE1_E2.set(xInd, yInd, p_nE1_E2);
			}
		}
		assertValidProbs(occupancy, prob_E1_E2, prob_E1_nE2, prob_nE1_E2, prob_nE1_nE2);
	}
	
	private static void adjustForEdge(boolean xAxis,
			EvenlyDiscrXYZ_DataSet occupancy, EvenlyDiscrXYZ_DataSet curOcc,
			EvenlyDiscrXYZ_DataSet prob_E1_E2, EvenlyDiscrXYZ_DataSet prob_E1_nE2,
			EvenlyDiscrXYZ_DataSet prob_nE1_E2, EvenlyDiscrXYZ_DataSet prob_nE1_nE2) {
		int num1,num2;
		if (xAxis) {
			num1 = occupancy.getNumX();
			num2 = occupancy.getNumY();
		} else {
			num1 = occupancy.getNumY();
			num2 = occupancy.getNumX();
		}
		for (int i=0; i<num1-1; i++) {
			double target, current;
			if (xAxis) {
				target = occupancy.get(i+1, 0);
				current = curOcc.get(i+1, 0);
			} else {
				target = occupancy.get(0, i+1);
				current = curOcc.get(0, i+1);
			}
			if (current == 0d || target == 0d)
				continue;
			double scalar = target / current;
			for (int j=0; j<num2; j++) {
				int xInd, yInd;
				if (xAxis) {
					xInd = i;
					yInd = j;
				} else {
					xInd = j;
					yInd = i;
				}
				if (occupancy.get(xInd, yInd) == 0d)
					continue;
				
				double p_E1_E2 = prob_E1_E2.get(xInd, yInd);
				double p_E1_nE2 = prob_E1_nE2.get(xInd, yInd);
				double p_nE1_E2 = prob_nE1_E2.get(xInd, yInd);
				double p_nE1_nE2 = prob_nE1_nE2.get(xInd, yInd);
				
				double curSum = p_E1_E2+p_E1_nE2+p_nE1_E2+p_nE1_nE2;
				Preconditions.checkState((float)curSum == 1f,
						"Bad input sum to edge: "+p_E1_E2+"+"+p_E1_nE2+"+"+p_nE1_E2+"+"+p_nE1_nE2+" = "+curSum);
				
				if (xAxis) {
					if (p_nE1_E2 == 0d)
						continue;
					p_nE1_E2 *= scalar;
					if (p_nE1_E2 + p_nE1_nE2 > 1d)
						p_nE1_E2 = 1d - p_nE1_nE2;
					
					double origOthers = p_E1_nE2 + p_E1_E2;
					double leftOthers = 1d - (p_nE1_E2 + p_nE1_nE2);
					
					if (origOthers == 0d) {
						p_E1_nE2 = 0d;
						p_E1_E2 = 0d;
						p_nE1_E2 = 1d - p_nE1_nE2;
					} else {
						p_E1_nE2 = p_E1_nE2*leftOthers/origOthers;
						p_E1_E2 = p_E1_E2*leftOthers/origOthers;
					}
					curSum = p_E1_E2+p_E1_nE2+p_nE1_E2+p_nE1_nE2;
					Preconditions.checkState((float)curSum == 1f,
							"Bad corrected sum to X edge: "+p_E1_E2+"+"+p_E1_nE2+"+"+p_nE1_E2+"+"+p_nE1_nE2+" = "+curSum
							+"\norigOthers: "+origOthers+", leftOther: "+leftOthers);
				} else {
					if (p_E1_nE2 == 0d)
						continue;
					p_E1_nE2 *= scalar;
					if (p_E1_nE2 + p_nE1_nE2 > 1d)
						p_E1_nE2 = 1d - p_nE1_nE2;
					
					double origOthers = p_nE1_E2 + p_E1_E2;
					double leftOthers = 1d - (p_E1_nE2 + p_nE1_nE2);
					
					if (origOthers == 0d) {
						p_nE1_E2 = 0d;
						p_E1_E2 = 0d;
						p_E1_nE2 = 1d - p_nE1_nE2;
					} else {
						p_nE1_E2 = p_nE1_E2*leftOthers/origOthers;
						p_E1_E2 = p_E1_E2*leftOthers/origOthers;
					}
					curSum = p_E1_E2+p_E1_nE2+p_nE1_E2+p_nE1_nE2;
					Preconditions.checkState((float)curSum == 1f,
							"Bad corrected sum to Y edge: "+p_E1_E2+"+"+p_E1_nE2+"+"+p_nE1_E2+"+"+p_nE1_nE2+" = "+curSum
							+"\norigOthers: "+origOthers+", leftOther: "+leftOthers);
				}
				
				prob_E1_E2.set(xInd, yInd, p_E1_E2);
				prob_E1_nE2.set(xInd, yInd, p_E1_nE2);
				prob_nE1_E2.set(xInd, yInd, p_nE1_E2);
			}
		}
		assertValidProbs(occupancy, prob_E1_E2, prob_E1_nE2, prob_nE1_E2, prob_nE1_nE2);
	}
	
	/**
	 * checks all probabilities are valid (0 <= P <= 1) and all for state transition probs sum to 1
	 */
	private static void assertValidProbs(EvenlyDiscrXYZ_DataSet occupancy,
			EvenlyDiscrXYZ_DataSet prob_E1_E2, EvenlyDiscrXYZ_DataSet prob_E1_nE2,
			EvenlyDiscrXYZ_DataSet prob_nE1_E2, EvenlyDiscrXYZ_DataSet prob_nE1_nE2) {
		for (int xInd=0; xInd<prob_E1_E2.getNumX(); xInd++) {
			for (int yInd=0; yInd<prob_E1_E2.getNumY(); yInd++) {
				if (occupancy.get(xInd, yInd) == 0d)
					continue;
				double p_E1_E2 = prob_E1_E2.get(xInd, yInd);
				double p_E1_nE2 = prob_E1_nE2.get(xInd, yInd);
				double p_nE1_E2 = prob_nE1_E2.get(xInd, yInd);
				double p_nE1_nE2 = prob_nE1_nE2.get(xInd, yInd);
				
				assertValidProb(p_E1_E2);
				assertValidProb(p_E1_nE2);
				assertValidProb(p_nE1_E2);
				assertValidProb(p_nE1_nE2);
				
				double sum = p_E1_E2 + p_E1_nE2 + p_nE1_E2 + p_nE1_nE2;
				Preconditions.checkState((float)sum == 1f, "Bad sum: "+sum);
			}
		}
	}
	
	private static Random r = new Random();
	
	private static void assertValidProb(double prob) {
		Preconditions.checkState(Doubles.isFinite(prob) && prob >= 0d && prob <= 1d, "Bad prob: "+prob);
	}
	
	private static EvenlyDiscrXYZ_DataSet calcOcc(EvenlyDiscrXYZ_DataSet origOcc,
			EvenlyDiscrXYZ_DataSet prob_E1_E2, EvenlyDiscrXYZ_DataSet prob_E1_nE2,
			EvenlyDiscrXYZ_DataSet prob_nE1_E2, EvenlyDiscrXYZ_DataSet prob_nE1_nE2,
			int samples, int[] initialState) {
		EvenlyDiscrXYZ_DataSet occ = newSameSize(prob_E1_E2);
		
		int[] curState = initialState;
		
		for (int i=0; i<samples; i++) {
			if (curState[0] >= occ.getNumX()-1)
				curState[0] = 0;
			if (curState[1] >= occ.getNumY()-1)
				curState[1] = 0;
			occ.set(curState[0], curState[1], 1d+occ.get(curState[0], curState[1]));
			double rand = r.nextDouble();
			
			double cmlProb = prob_nE1_nE2.get(curState[0], curState[1]);
			if (rand < cmlProb) {
				curState = new int[] { curState[0]+1, curState[1]+1};
				continue;
			}
			
			cmlProb += prob_E1_nE2.get(curState[0], curState[1]);
			if (rand < cmlProb) {
				curState = new int[] { 0, curState[1]+1};
				continue;
			}
			
			cmlProb += prob_nE1_E2.get(curState[0], curState[1]);
			if (rand < cmlProb) {
				curState = new int[] { curState[0]+1, 0};
				continue;
			}
			
			curState = new int[] { 0, 0 };
		}
		// now mask out any area that are zero in the input occupancy
		for (int xInd=0; xInd<occ.getNumX(); xInd++) {
			for (int yInd=0; yInd<occ.getNumY(); yInd++) {
				if (origOcc.get(xInd, yInd) == 0d)
					occ.set(xInd, yInd, 0d);
			}
		}
		occ.scale(1d/occ.getSumZ());
		return occ;
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
	
	private static void generateScatterComparisonPlots(MarkovChain origChain, MarkovChain calcChain, File outputDir,
			List<? extends Named> names) throws IOException {
		int numOcc = 10;
		DefaultXY_DataSet[] scatter_E1_E2 = new DefaultXY_DataSet[numOcc];
		DefaultXY_DataSet[] scatter_E1_nE2 = new DefaultXY_DataSet[numOcc];
		DefaultXY_DataSet[] scatter_nE1_E2 = new DefaultXY_DataSet[numOcc];
		DefaultXY_DataSet[] scatter_nE1_nE2 = new DefaultXY_DataSet[numOcc];
		for (int i=0; i<numOcc; i++) {
			scatter_E1_E2[i] = new DefaultXY_DataSet();
			scatter_E1_nE2[i] = new DefaultXY_DataSet();
			scatter_nE1_E2[i] = new DefaultXY_DataSet();
			scatter_nE1_nE2[i] = new DefaultXY_DataSet();
		}
		double maxOcc = 0d;
		for (int[] state : origChain.getOccupancy().getStates())
			maxOcc = Math.max(maxOcc, origChain.getOccupancy().getFrequency(state));
		EvenlyDiscretizedFunc occFunc = new EvenlyDiscretizedFunc(0d, maxOcc, numOcc);
		for (int[] state : calcChain.getOccupancy().getStates()) {
			double occ = origChain.getOccupancy().getFrequency(state);
			double orig_E1_E2 = origChain.getTransitionProb(state, new int[] {0,0});
			double calc_E1_E2 = calcChain.getTransitionProb(state, new int[] {0,0});
			double orig_E1_nE2 = origChain.getTransitionProb(state, new int[] {0,state[1]+1});
			double calc_E1_nE2 = calcChain.getTransitionProb(state, new int[] {0,state[1]+1});
			double orig_nE1_E2 = origChain.getTransitionProb(state, new int[] {state[0]+1,0});
			double calc_nE1_E2 = calcChain.getTransitionProb(state, new int[] {state[0]+1,0});
			double orig_nE1_nE2 = origChain.getTransitionProb(state, new int[] {state[0]+1,state[1]+1});
			double calc_nE1_nE2 = calcChain.getTransitionProb(state, new int[] {state[0]+1,state[1]+1});
			
			int index = occFunc.getClosestXIndex(occ);
			
			scatter_E1_E2[index].set(orig_E1_E2, calc_E1_E2);
			scatter_E1_nE2[index].set(orig_E1_nE2, calc_E1_nE2);
			scatter_nE1_E2[index].set(orig_nE1_E2, calc_nE1_E2);
			scatter_nE1_nE2[index].set(orig_nE1_nE2, calc_nE1_nE2);
		}
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, maxOcc);
		List<PlotCurveCharacterstics> scatterChars = Lists.newArrayList();
		DefaultXY_DataSet ref_line = new DefaultXY_DataSet();
		ref_line.set(0d, 0d);
		ref_line.set(1d, 1d);
		scatterChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.GRAY));
		for (int i=0; i<numOcc; i++)
			scatterChars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, cpt.getColor((float)occFunc.getX(i))));
		List<PlotSpec> specs = Lists.newArrayList();
		
		String title = names.get(0).getName()+" and "+names.get(1).getName();
		
		List<XY_DataSet> elems_E1_E2 = Lists.newArrayList();
		elems_E1_E2.add(ref_line);
		elems_E1_E2.addAll(Lists.newArrayList(scatter_E1_E2));
		PlotSpec spec_E1_E2 = new PlotSpec(elems_E1_E2, scatterChars, title, "Original", "Calculated");
		XYTextAnnotation ann = new XYTextAnnotation("P(E1 E2 | s1 s2)", 0.1, 0.9);
		ann.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 18));
		ann.setTextAnchor(TextAnchor.BASELINE_LEFT);
		spec_E1_E2.setPlotAnnotations(Lists.newArrayList(ann));
		specs.add(spec_E1_E2);
		
		List<XY_DataSet> elems_E1_nE2 = Lists.newArrayList();
		elems_E1_nE2.add(ref_line);
		elems_E1_nE2.addAll(Lists.newArrayList(scatter_E1_nE2));
		PlotSpec spec_E1_nE2 = new PlotSpec(elems_E1_nE2, scatterChars, title, "Original", "Calculated");
		ann = new XYTextAnnotation("P(E1 !E2 | s1 s2)", 0.1, 0.9);
		ann.setTextAnchor(TextAnchor.BASELINE_LEFT);
		ann.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 18));
		spec_E1_nE2.setPlotAnnotations(Lists.newArrayList(ann));
		specs.add(spec_E1_nE2);
		
		List<XY_DataSet> elems_nE1_E2 = Lists.newArrayList();
		elems_nE1_E2.add(ref_line);
		elems_nE1_E2.addAll(Lists.newArrayList(scatter_nE1_E2));
		PlotSpec spec_nE1_E2 = new PlotSpec(elems_nE1_E2, scatterChars, title, "Original", "Calculated");
		ann = new XYTextAnnotation("P(!E1 E2 | s1 s2)", 0.1, 0.9);
		ann.setTextAnchor(TextAnchor.BASELINE_LEFT);
		ann.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 18));
		spec_nE1_E2.setPlotAnnotations(Lists.newArrayList(ann));
		specs.add(spec_nE1_E2);
		
		List<XY_DataSet> elems_nE1_nE2 = Lists.newArrayList();
		elems_nE1_nE2.add(ref_line);
		elems_nE1_nE2.addAll(Lists.newArrayList(scatter_nE1_nE2));
		PlotSpec spec_nE1_nE2 = new PlotSpec(elems_nE1_nE2, scatterChars, title, "Original", "Calculated");
		ann = new XYTextAnnotation("P(!E1 !E2 | s1 s2)", 0.1, 0.9);
		ann.setTextAnchor(TextAnchor.BASELINE_LEFT);
		ann.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 18));
		spec_nE1_nE2.setPlotAnnotations(Lists.newArrayList(ann));
		specs.add(spec_nE1_nE2);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		gp.setBackgroundColor(Color.WHITE);
		
		Range xRange = new Range(0d, 1d);
		Range yRange = new Range(0d, 1d);
		List<Range> xRanges = Lists.newArrayList(xRange);
		List<Range> yRanges = Lists.newArrayList();
		for (int i=0; i<specs.size(); i++)
			yRanges.add(yRange);
		
//		gp.setUserBounds(getXRange(), getYRange());
		gp.drawGraphPanel(specs, false, false, xRanges, yRanges);
		gp.getChartPanel().setSize(800, 1800);
		
		String prefix = "markov_scatter_"+PeriodicityPlotter.getFileSafeString(names.get(0).getName())
				+"_"+PeriodicityPlotter.getFileSafeString(names.get(1).getName());
		
		gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
		gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
	}

	public static void main(String[] args) throws IOException {
		double distSpacing = 10d;
		List<RuptureIdentifier> rupIdens = SynchIdens.getIndividualFaults(7d, 10d,
//				SynchFaults.SAF_MOJAVE, SynchFaults.SAF_COACHELLA, SynchFaults.SAN_JACINTO, SynchFaults.SAF_CARRIZO);
//				SynchFaults.SAF_MOJAVE, SynchFaults.SAF_COACHELLA, SynchFaults.SAN_JACINTO);
//				SynchFaults.SAF_MOJAVE, SynchFaults.SAF_COACHELLA);
//				SynchFaults.SAF_CARRIZO, SynchFaults.SAF_CHOLAME);
//				SynchFaults.SAF_CARRIZO, SynchFaults.SAF_COACHELLA);
				SynchFaults.SAF_COACHELLA, SynchFaults.SAN_JACINTO);
//				SynchFaults.SAF_COACHELLA, SynchFaults.SAF_CHOLAME);
//				SynchFaults.SAF_MOJAVE, SynchFaults.GARLOCK_WEST);
		
		List<? extends SimulatorEvent> events = new SimAnalysisCatLoader(true, rupIdens, true).getEvents();
		
		List<int[]> fullPath = MarkovChainBuilder.getStatesPath(distSpacing, events, rupIdens, 0d);
		
		EmpiricalMarkovChain origChain = new EmpiricalMarkovChain(fullPath, distSpacing);
		
		EvenlyDiscrXYZ_DataSet occupancy =
				new StateSpacePlotter(origChain, rupIdens, null).getOccupancy(0, 1);
		occupancy.scale((double)fullPath.size()/occupancy.getSumZ());
		
		MarkovChain calcChain = calc(occupancy, 0);
		
		File outputDir = new File("/tmp");
		StateSpacePlotter origPlot = new StateSpacePlotter(origChain, rupIdens, new File(outputDir, "occ_based_orig"));
		StateSpacePlotter calcPlot = new StateSpacePlotter(calcChain, rupIdens, new File(outputDir, "occ_based_calc"));
		
		origPlot.plotStateTrans(true, true);
		calcPlot.plotStateTrans(true, true);
		
		origPlot.plotOccupancies(false);
		calcPlot.plotOccupancies(false);
		
		generateScatterComparisonPlots(origChain, calcChain, outputDir, rupIdens);
	}

}
