package scratch.ned.ETAS_ERF.sandbox;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.ListIterator;

import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.param.AleatoryMagAreaStdDevParam;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.UCERF2;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.magdist.ArbIncrementalMagFreqDist;

import scratch.UCERF3.erf.FaultSystemSolutionPoissonERF;
import scratch.UCERF3.erf.ETAS.IntegerPDF_FunctionSampler;

/**
 * This class store information about all ruptures that nucleate inside this geographic block.
 * @author field
 *
 */
public class EqksAtPoint {
	
	ArrayList<Integer> rupIndexN_List;	// this stores the Nth index of ERF ruptures that nucleate at this point 
	ArrayList<Double> rupRateList;	// this holds the nucleation rate for each rupture at this point
	ArrayList<Double> rupFractList;	// this holds the fraction of the rupture surface at this point
	ArrayList<Integer> srcIndexList;	// this stores the index of each point source that nucleate at this point
	ArrayList<Double> srcRateList;	// this holds the nucleation rate for each rupture inside the block
	ArrayList<Double> srcFractList;	// this holds the fraction of the source that necleates at this point

	int[] rupIndexN_Array;
	double[] rupRateArray;
	double[] rupFractArray;
	int[] srcIndexArray;
	double[] srcRateArray;
	double[] srcFractArray;

	double totalRateInside = -1;
	
	IntegerPDF_FunctionSampler randomSampler;	// this is for random sampling of ruptures


	/**
	 * Constructor
	 */
	public EqksAtPoint() {
		rupIndexN_List = new ArrayList<Integer>();
		rupRateList = new ArrayList<Double>();
		rupFractList = new ArrayList<Double>();
		srcIndexList = new ArrayList<Integer>();
		srcRateList = new ArrayList<Double>();
		srcFractList = new ArrayList<Double>();

	}	
	
	
	public EqksAtPoint(	int[] rupIndexN_Array, double[] rupRateInsideArray, double[] rupFractInsideArray,
			int[] srcIndexN_Array, double[] srcRateInsideArray, double[] srcFractInsideArray) {

		this.rupIndexN_Array = rupIndexN_Array;
		this.rupRateArray = rupRateInsideArray;
		this.rupFractArray = rupFractInsideArray;
		this.srcIndexArray = srcIndexN_Array;
		this.srcRateArray = srcRateInsideArray;
		this.srcFractArray = srcFractInsideArray;

	}
	
	/**
	 * This gives the total rate at which ruptures nucleate inside the block
	 * @return
	 */
	public double getTotalRateInside() {
		// check to see whether it's already been calculated
		this.getRandomSampler();
		if(randomSampler != null)
			return randomSampler.calcSumOfY_Vals();
		else
			return 0;
//		if(totalRateInside == -1) {
//			totalRateInside=0;
//			for(double rate:rupRateInsideArray) 
//				totalRateInside += rate;
//			for(double rate:srcRateInsideArray) 
//				totalRateInside += rate;
//		}
//		return totalRateInside;
	}
	
		
	
	/**
	 * This adds the rate for the given nthRupIndex
	 * @param rate
	 * @param nthRupIndex
	 * @param fract - the fraction of the rupture this given rate represents
	 */
	public void addRupRate(double rate, int nthRupIndex, double fract) {
		int localIndex = rupRateList.indexOf(nthRupIndex);		// faster with hashmap?
		if(localIndex<0) {	// index does not exist
			int size = rupIndexN_List.size();
			if(size>0) 
				if(nthRupIndex < rupIndexN_List.get(size-1))
					throw new RuntimeException("Ruptures must be entered in order");
			rupRateList.add(rate);
			rupIndexN_List.add(nthRupIndex);
			rupFractList.add(fract);
		}
		else {	// index exists; add the rate
			double newRate = rupRateList.get(localIndex)+rate;
			rupRateList.set(localIndex, newRate);
			double newFract = rupFractList.get(localIndex)+fract;
			rupFractList.set(localIndex, newFract);
		}
	}

	
	/**
	 * This adds the rate for the given point src
	 * @param rate
	 * @param nthRupIndex
	 */
	public void addSrcRate(double rate, int srcIndex, double fract) {
		int localIndex = srcRateList.indexOf(srcIndex);		// faster with hashmap?
		if(localIndex<0) {	// index does not exist
			int size = srcIndexList.size();
			if(size>0) 
				if(srcIndex < srcIndexList.get(size-1))
					throw new RuntimeException("Sources must be entered in order");
			if(rate>0) {
				srcRateList.add(rate);
				srcIndexList.add(srcIndex);
				srcFractList.add(fract);
			}			
		}
		else {	// index exists; add the rate
			double newRate = srcRateList.get(localIndex)+rate;
			srcRateList.set(localIndex, newRate);
			double newFract = srcFractList.get(localIndex)+fract;
			srcFractList.set(localIndex, newFract);
		}
	}

	/**
	 * This changes the rate for the specified rupture
	 * @param totRupRate - total rate, which will get reduced by the faction inside value
	 * @param nthRupIndex - the index of the nth rupture in the ERF
	 */
	public void changeRupRate(double totRupRate, int nthRupIndex) {
		int localIndex = Arrays.binarySearch(rupIndexN_Array, nthRupIndex);;
		if(localIndex < 0)	{
			throw new RuntimeException("nthRupIndex="+nthRupIndex+" not found (was rate zero when this object was created?)");
		}
		double oldRate = rupRateArray[localIndex];
		double newRate = totRupRate*rupFractArray[localIndex];
		// update totalRate
		if(totalRateInside != -1)
			totalRateInside += newRate-oldRate;
		// update sampler
		if(randomSampler != null)
				randomSampler.set(localIndex,newRate);
		rupRateArray[localIndex] = newRate;
	}
	
	
	/**
	 * This returns a 2-element array of ints, where the first value indicates
	 * the index type (0 for rup index, or 1 for src index), and the second element
	 * is the index (nth rup or ith src).
	 * @return
	 */
	public int[] getRandomRupOrSrc() {
		// make random sampler if it doesn't already exist
		getRandomSampler();
		if(randomSampler == null)
			throw new RuntimeException("No ruptures at this point");
		int[] toReturn = new int[2];
		int localIndex = randomSampler.getRandomInt();
		if (localIndex < rupIndexN_Array.length) {
			toReturn[0] = 0;
			toReturn[1] = rupIndexN_Array[localIndex];
		}
		else {
			toReturn[0] = 1;
			toReturn[1] = srcIndexArray[localIndex-rupIndexN_Array.length];
		}
		return toReturn;
	}
	
	/**
	 * This creates (if not already existent) and returns the randomEqkRupSampler
	 * @return
	 */
	public IntegerPDF_FunctionSampler getRandomSampler() {
		if(randomSampler == null) {
			int numPts = rupIndexN_Array.length+srcIndexArray.length;
			if(numPts>0) {
				randomSampler = new IntegerPDF_FunctionSampler(numPts);
				for(int i=0;i<rupIndexN_Array.length;i++) 
					randomSampler.set(i,rupRateArray[i]);
				for(int i=0;i<srcIndexArray.length;i++) 
					randomSampler.set(i+rupIndexN_Array.length,srcRateArray[i]);				
			}
			else
				randomSampler = null;
		}
		return randomSampler;
	}
	

	/**
	 * This converts the array lists to arrays (to reduce memory usage)
	 */
	public void finishAndShrinkSize() {
		
		int num = rupIndexN_List.size();
		rupIndexN_Array = new int[num];
		rupRateArray = new double[num];
		rupFractArray = new double[num];
		
		for(int i=0;i<num;i++) {
			rupIndexN_Array[i] = rupIndexN_List.get(i);
			rupRateArray[i] = rupRateList.get(i);
			rupFractArray[i] = rupFractList.get(i);
		}
				
		num = srcIndexList.size();
		srcIndexArray = new int[num];
		srcRateArray = new double[num];
		srcFractArray = new double[num];
		
		for(int i=0;i<num;i++) {
			srcIndexArray[i] = srcIndexList.get(i);
			srcRateArray[i] = srcRateList.get(i);
			srcFractArray[i] = srcFractList.get(i);
		}
				
		rupIndexN_List = null;
		rupRateList = null;
		srcIndexList = null;
		srcRateList = null;

	}
	
	public int[] getRupIndexN_Array() {
		return rupIndexN_Array;
	}
	
	public double[] getRupRateInsideArray() {
		return rupRateArray;
	}
	
	public double[] getRupFractInsideArray() {
		return rupFractArray;
	}
	
	public int[] getSrcIndexN_Array() {
		return srcIndexArray;
	}
	
	public double[] getSrcRateInsideArray() {
		return srcRateArray;
	}
	
	public double[] getSrcFractInsideArray() {
		return srcFractArray;
	}


	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		

	}
	
	
	/** This computes the mag-freq dist for this point
	 * 
	 * @return
	 */
	public ArbIncrementalMagFreqDist getMagFreqDist(FaultSystemSolutionPoissonERF erf) {
		getRandomSampler(); // do this to make sure it exists
		ArbIncrementalMagFreqDist magDist = new ArbIncrementalMagFreqDist(2.05, 8.95, 70);
		for(int j=0; j<rupIndexN_Array.length; j++) {
			double mag = erf.getNthRupture(rupIndexN_Array[j]).getMag();
			magDist.addResampledMagRate(mag, randomSampler.getY(j), true);
		}
		double duration = erf.getTimeSpan().getDuration();
		for(int j=rupIndexN_Array.length; j<rupIndexN_Array.length+srcIndexArray.length; j++) {
			ProbEqkSource src = erf.getSource(srcIndexArray[j-rupIndexN_Array.length]);
			double totSrcRate = src.computeTotalEquivMeanAnnualRate(duration);
			double srcRateAtPt = randomSampler.getY(j);
			for(ProbEqkRupture rup : src) {
				magDist.addResampledMagRate(rup.getMag(), rup.getMeanAnnualRate(duration)*srcRateAtPt/totSrcRate, true);
			}
		}
		return magDist;
	}

	
	/** This computes the expected, normalized mag-freq dist for the block (total rate is 1.0)
	 * 
	 * @return
	 */
	public ArbIncrementalMagFreqDist getMagProbDist(FaultSystemSolutionPoissonERF erf) {
		throw new RuntimeException("Method needs work");
//		ArbIncrementalMagFreqDist magDist = new ArbIncrementalMagFreqDist(2.05, 8.95, 70);
//		for(int j=0; j<rupIndexN_Array.length; j++) {
//			double mag = erf.getNthRupture(rupIndexN_Array[j]).getMag();
//			magDist.addResampledMagRate(mag, rupRateInsideArray[j], true);
//		}
//		magDist.scaleToCumRate(2.05, 1);
//		return magDist;
	}
}
