package scratch.ned.ETAS_ERF;

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
public class EqksInGeoBlock {
	
	double minLat, maxLat, minLon,maxLon, minDepth, maxDepth; // the dimensions of the block
	FaultSystemSolutionPoissonERF erf;					// reference to the erf this block is used with
	
	// these lists hold information about each rupture that nucleates in this block
	ArrayList<Integer> rupIndexN_List;	// this stores the Nth index of the rupture inside the ERF
	ArrayList<Double> rateInsideList;	// this holds the nucleation rate for each rupture inside the block
	ArrayList<Double> fractInsideList;	// this holds the fraction of the rupture that's inside the block
	ArrayList<Double> magList;			// this holds the magnitude for the ruptures that nucleate inside
	
	Location blockCenterLoc;
	
	double totalRateInside = -1;
	double blockVolume = -1;
	
	IntegerPDF_FunctionSampler randomEqkRupSampler;	// this is for random sampling of ruptures


	/**
	 * Constructor that takes limits to define the blocks
	 * @param minLat
	 * @param maxLat
	 * @param minLon
	 * @param maxLon
	 * @param minDepth
	 * @param maxDepth
	 */
	public EqksInGeoBlock(double minLat, double maxLat, double minLon, double maxLon, 
			double minDepth, double maxDepth) {
		this.minLat=minLat;
		this.maxLat=maxLat;
		this.minLon=minLon;
		this.maxLon=maxLon;
		this.minDepth=minDepth;
		this.maxDepth=maxDepth;
		rupIndexN_List = new ArrayList<Integer>();
		rateInsideList = new ArrayList<Double>();
		fractInsideList = new ArrayList<Double>();
		magList = new ArrayList<Double>();
		blockCenterLoc = new Location((minLat+maxLat)/2,(minLon+maxLon)/2,(minDepth+maxDepth)/2);
		
	}
	
	
	/**
	 * Constructor that takes bin center and width info to define the lat/lon bounds of the block
	 * @param blockCenerLoc
	 * @param latLonBlockWidth
	 * @param minDepth
	 * @param maxDepth
	 */
	public EqksInGeoBlock(Location blockCenerLoc, double latLonBlockWidth, double minDepth, double maxDepth) {
		this(blockCenerLoc.getLatitude()-latLonBlockWidth/2, 
				blockCenerLoc.getLatitude()+latLonBlockWidth/2, 
				blockCenerLoc.getLongitude()-latLonBlockWidth/2,
				blockCenerLoc.getLongitude()+latLonBlockWidth/2, 
				minDepth, maxDepth);
	}

	
	/**
	 * This constructor takes and processes an ERF too
	 * @param minLat
	 * @param maxLat
	 * @param minLon
	 * @param maxLon
	 * @param minDepth
	 * @param maxDepth
	 * @param erf
	 */
	public EqksInGeoBlock(double minLat, double maxLat, double minLon, double maxLon, 
			double minDepth, double maxDepth, FaultSystemSolutionPoissonERF erf) {
		
		this(minLat, maxLat, minLon, maxLon, minDepth, maxDepth);
		this.erf=erf;
		processERF();
		
	}
	
	/**
	 * This returns the center location of the block
	 * @return
	 */
	public Location getBlockCenterLoc() {return blockCenterLoc;}
	
	
	/**
	 * This writes to system the following for each rupture that nucleates inside the block:
	 * rupIndexN_List, rateInside, fractInside, randomEqkRupSampler.getY(i), and srcName (if erf given)
	 */
	public void writeResults() {
		// do the following to make sure randomEqkRupSampler has been created
		getRandomSampler();
		System.out.println("TotalRateInside="+getTotalRateInside());
		for(int i=0; i<rupIndexN_List.size();i++) {
			System.out.print("\t"+i+"\t"+rupIndexN_List.get(i)+"\t"+
					rateInsideList.get(i)+"\t"+fractInsideList.get(i)+"\t"+
					magList.get(i)+"\t"+randomEqkRupSampler.getY(i));
			if(erf != null)
				System.out.print("\t"+erf.getSource(erf.getSrcIndexForNthRup(rupIndexN_List.get(i))).getName()+"\n");
			else
				System.out.print("\n");
		}
	}
	
	/**
	 * This returns the approximate volume of the block
	 * (useful for normalizing rates into rate densities)
	 * @return volume in sq km
	 */
	public double getBlockVolume() {
		if(blockVolume == -1) {
			Location loc1 = new Location(minLat,minLon);
			Location loc2 = new Location(minLat,maxLon);
			Location loc3 = new Location(maxLat,minLon);
			double distLat = LocationUtils.horzDistanceFast(loc1, loc3);
			double distLon = LocationUtils.horzDistanceFast(loc1, loc2);
			blockVolume = distLat*distLon*(maxDepth-minDepth);
		}
		return blockVolume;
	}
	
	
	/**
	 * This returns the cube-root of the block volume
	 * (ave block size)
	 * @return
	 */
	public double getAveBlockSize() {
		return Math.pow(getBlockVolume(),0.3333333);
	}
	
	
	
	/**
	 * This gives the total rate at which ruptures nucleate inside the block
	 * @return
	 */
	public double getTotalRateInside() {
		// check to see whether it's already been calculated
		if(totalRateInside == -1) {
			totalRateInside=0;
			for(Double rate:this.rateInsideList) totalRateInside += rate;
		}
		return totalRateInside;
	}
	
	
	/**
	 * This processes the ERF to generate the information for  ruptures that nucleate inside the block
	 */
	private void processERF() {
		int numSrc = erf.getNumSources();
		double forecastDuration = erf.getTimeSpan().getDuration();
		for(int s=0;s<numSrc;s++) {
			ProbEqkSource src = erf.getSource(s);
			for(int r=0; r<src.getNumRuptures();r++) {
				ProbEqkRupture rup = src.getRupture(r);
				int nthRup = erf.getIndexN_ForSrcAndRupIndices(s, r);
				processRupture(rup, nthRup, forecastDuration);
			}
		}
	}
	
	
	/**
	 * This processes a given rupture (determines whether it 
	 * nucleates inside and saves info if it does)
	 * @param rup
	 * @param srcIndex
	 * @param rupIndex
	 * @param forecastDuration
	 */
	public double processRupture(ProbEqkRupture rup, int nthRup, double forecastDuration) {
		double rate = rup.getMeanAnnualRate(forecastDuration);
		double fractInside = 0;
		LocationList locList = rup.getRuptureSurface().getEvenlyDiscritizedListOfLocsOnSurface();
		int numLoc = locList.size();
		int numLocInside = 0;
		for(Location loc : locList)
			if(isLocInside(loc)) numLocInside+=1;
		fractInside = (double)numLocInside/(double)numLoc;	// return this regardless of rate
		if(numLocInside>0 && rate > 0) {	// only fill these in if it's a non-zero rate
			fractInsideList.add(fractInside);
			rateInsideList.add(rate*fractInside);
			rupIndexN_List.add(nthRup);
			magList.add(rup.getMag());
		}			
		return fractInside;
	}
	
	
	public void testThisBlock(FaultSystemSolutionPoissonERF erf) {
		
		EqksInGeoBlock testBlock = new EqksInGeoBlock(minLat, maxLat, minLon, maxLon, minDepth, maxDepth);
		double duration = erf.getTimeSpan().getDuration();
		for(int n=0; n<erf.getTotNumRups();n++) {
			testBlock.processRupture(erf.getNthRupture(n), n, duration);
		}
		
		ArrayList<Double> thisRates = this.getRateInsideList();
		ArrayList<Double> testRates = testBlock.getRateInsideList();
		ArrayList<Double> thisFracts = this.getFractInsideList();
		ArrayList<Double> testFracts = testBlock.getFractInsideList();
		if(thisRates.size() != testRates.size())
			throw new RuntimeException("Error\tthisRates.size()="+thisRates.size()+"\ttestRates.size()="+testRates.size());
		for(int i=0;i<thisRates.size();i++) {
			double ratioRates = testRates.get(i)/thisRates.get(i);
			double ratioFracts = testFracts.get(i)/thisFracts.get(i);
			if(ratioRates <0.999 || ratioRates > 1.001) {
				System.out.println("PROBLEM: "+rupIndexN_List.get(i)+"\t"+ratioRates+"\t"+testRates.get(i)+"\t"+thisRates.get(i)+"\t"+ratioFracts+"\t"+testFracts.get(i)+"\t"+thisFracts.get(i)+"\t");

				
				//				testBlock.writeResults();
//				this.writeResults();
//				System.exit(0);
			}

		}
		
	}
	
	
	/**
	 * This adds to the given info to the lists (e.g., if determined externally)
	 * @param rate
	 * @param fracInside
	 * @param nthRupIndex
	 * @param magnitude
	 */
	public void processRate(double rate, double fracInside, int nthRupIndex, double magnitude) {
		if(rate > 0) {
			fractInsideList.add(fracInside);
			rateInsideList.add(rate);
			rupIndexN_List.add(nthRupIndex);
			magList.add(magnitude);
		}			
	}
	
	/**
	 * This changes the rate for the specified rupture
	 * @param totRupRate - total rate, which will get reduced by the faction inside value
	 * @param nthRupIndex - the index of the nth rupture in the ERF
	 */
	public void changeRate(double totRupRate, int nthRupIndex) {
		int localIndex = findLocalIndex(nthRupIndex);
		if(localIndex < 0)	// return if index not found (finite rupture not hitting sub-block)
			return;
		double oldRate = rateInsideList.get(localIndex);
		double newRate = totRupRate*fractInsideList.get(localIndex);
		// update totalRate
		if(totalRateInside != -1)
			totalRateInside += newRate-oldRate;
		// update sampler
		if(randomEqkRupSampler != null)
				randomEqkRupSampler.set(localIndex,newRate);
		rateInsideList.set(localIndex, newRate);
	}
	
	
	public double tempGetRandomEqkRupSamplerY_Val(int nthRupIndex) {
		getRandomSampler();
		int localIndex = findLocalIndex(nthRupIndex);
		if(localIndex >= 0)
			return randomEqkRupSampler.getY(localIndex);
		else
			return -1;
	}
	
	/**
	 * This returns the local index of the given nthRupIndex.
	 * 
	 * A negative value is returned if it's not found
	 * 
	 * Use a hashmap instead?  Presume ordering to speed up.
	 * @param nthRupIndex
	 * @return
	 */
	private int findLocalIndex(int nthRupIndex) {
		// duplicate rupIndexN_List to an array
		int[] rupIndexN_Array = new int[rupIndexN_List.size()];
		int last_n = -1;
		for(int i=0;i<rupIndexN_List.size();i++) {
			int n = rupIndexN_List.get(i);
			if(n>last_n) {
				rupIndexN_Array[i] = n;
				last_n = n;
			}
			else throw new RuntimeException("rupIndexN_List must be orderd");
			
		}
		return Arrays.binarySearch(rupIndexN_Array, nthRupIndex);
	}
	
	
	/**
	 * This divides the block into equal-sized sub-blocks, where the original block is sliced in all three
	 * dimensions, making the final number of blocks = numAlongLatLon*numAlongLatLon*numAlongDepth 
	 * 
	 * Important: this assumes that any point sources within the block should be equally divided among
	 * the sub blocks.
	 * @param numAlongLatLon
	 * @param numAlongDepth
	 * @return
	 */
	public ArrayList<EqksInGeoBlock> getSubBlocks(int numAlongLatLon, int numAlongDepth, FaultSystemSolutionPoissonERF erf) {
		ArrayList<EqksInGeoBlock> subBlocks = new ArrayList<EqksInGeoBlock>();
		double forecastDuration = erf.getTimeSpan().getDuration();
		int numSubBlocks = numAlongLatLon*numAlongLatLon*numAlongDepth;
		
		// make the sub-block boundaries
		double[] latBoundaries = new double[numAlongLatLon+1];
		double[] lonBoundaries = new double[numAlongLatLon+1];
		double[] depthBoundaries = new double[numAlongDepth+1];
		for(int i=0; i< numAlongLatLon+1; i++) {
			latBoundaries[i] = minLat+(double)i*(maxLat-minLat)/(double)numAlongLatLon;
			lonBoundaries[i] = minLon+(double)i*(maxLon-minLon)/(double)numAlongLatLon;
		}
		for(int i=0; i< numAlongDepth+1; i++) {
			depthBoundaries[i] = minDepth+(double)i*(maxDepth-minDepth)/(double)numAlongDepth;
		}
		
		// Make the sub blocks
		for(int latSlice=0; latSlice<numAlongLatLon; latSlice++) {
			for(int lonSlice=0; lonSlice<numAlongLatLon; lonSlice++) {
				for(int depSlice=0; depSlice<numAlongDepth; depSlice++) {
					double subMinLat = latBoundaries[latSlice];
					double subMaxLat = latBoundaries[latSlice+1];
					double subMinLon = lonBoundaries[lonSlice];
					double subMaxLon = lonBoundaries[lonSlice+1];
					double subMinDepth = depthBoundaries[depSlice];
					double subMaxDepth = depthBoundaries[depSlice+1];
					EqksInGeoBlock subBlock = new EqksInGeoBlock(subMinLat, subMaxLat, subMinLon, subMaxLon, subMinDepth, subMaxDepth);
					subBlocks.add(subBlock);
				}
			}
		}
		
		// put the rupture rates here into the sub-blocks
		for(int r=0; r<getNumRupsInside();r++) {
			int nthRup = rupIndexN_List.get(r);
			ProbEqkRupture rup = erf.getNthRupture(nthRup);
			if(!rup.getRuptureSurface().isPointSurface()) {	// if not point surface
				double testFraction=0;
				for(EqksInGeoBlock subBlock:subBlocks) {
					testFraction += subBlock.processRupture(rup, nthRup, forecastDuration);	
				}
				
				// test
				double ratio = testFraction/fractInsideList.get(r);
				if(ratio<0.999 || ratio > 1.001) {
					System.out.println("\tfracton diff:\ttestFraction="+testFraction+"\tfractInsideList.get(r)="+fractInsideList.get(r));
					System.out.println("\trup rate="+rup.getMeanAnnualRate(forecastDuration));
					LocationList list=rup.getRuptureSurface().getEvenlyDiscritizedListOfLocsOnSurface();
					for(Location loc: list) {
						if(this.isLocInside(loc)) {	// first make sure the point is in this block
							boolean gotOne=false;
							for(EqksInGeoBlock subBlock:subBlocks) {
								if(subBlock.isLocInside(loc))
									gotOne=true;
							}
							if(gotOne==false) {	// this loc was not found
								System.out.println("\tLoc not found: "+loc);
								System.out.println("\t"+nthRup+"\t"+isLocInside(loc)+"\t"+minLat+"\t"+maxLat+"\t"+minLon+"\t"+maxLon+"\t"+minDepth+"\t"+maxDepth);
								for(EqksInGeoBlock subBlock:subBlocks) {
									System.out.println("\t\t"+"\t"+subBlock.minLat+"\t"+subBlock.maxLat+"\t"+
											subBlock.minLon+"\t"+subBlock.maxLon+"\t"+subBlock.minDepth+"\t"+subBlock.maxDepth+"\t"+subBlock.isLocInside(loc));
								}
								System.out.println("\n");
								testThisBlock(erf);
								System.exit(0);
							}							
						}
					}
					System.exit(0);
				}
			}
			else { // assume point sources equally divided
				for(EqksInGeoBlock subBlock:subBlocks) {
					double rate = rup.getMeanAnnualRate(forecastDuration)/(double)numSubBlocks;
					double fracInside = 1.0/(double)numSubBlocks;
					subBlock.processRate(rate, fracInside, nthRup, rup.getMag());
				}
			}
		}

		// check total rates
		double totRate=0;
		for(EqksInGeoBlock block : subBlocks)
			totRate += block.getTotalRateInside();
		double testRate2 = getTotalRateInside();
		double ratio = totRate/testRate2;
		if(Math.abs(totRate) < 1e-15 && Math.abs(testRate2) < 1e-15)
			ratio = 1;
		if(ratio<0.999 || ratio>1.001) {
			System.out.println("PROBLEM: ratio="+ratio+";\ttotRate="+totRate+"\ttestRate2="+testRate2+"\n");
			// TEST
			for(int r=0; r<rupIndexN_List.size();r++) {
				double targetRate = rateInsideList.get(r);
				double summedRate=0;
				for(int b=0; b<subBlocks.size();b++) {
					EqksInGeoBlock blk = subBlocks.get(b);
					if(blk.getRateInsideList().size()>0) {
						double rate = blk.tempGetRandomEqkRupSamplerY_Val(rupIndexN_List.get(r));
						if(rate >= 0)
							summedRate+=rate;					
					}
				}
				String srcName= erf.getSource(erf.getSrcIndexForNthRup(rupIndexN_List.get(r))).getName();
				System.out.println("\t"+r+"\t"+rupIndexN_List.get(r)+"\t"+(summedRate/targetRate)+
						"\t"+summedRate+"\t"+targetRate+"\t"+srcName);
			}
			
			System.out.println("\ntestThisBlock:\n");
			testThisBlock(erf);
			System.exit(0);

		}

//		System.out.println("/nRate Check: "+totRate+" vs "+this.getTotalRateInside());
		
		return subBlocks;
	}


	
	/**
	 * This divides the block into equal-sized sub-blocks, where the original block is sliced in all three
	 * dimensions, making the final number of blocks = numAlongLatLon*numAlongLatLon*numAlongDepth 
	 * 
	 * Important: this assumes that any point sources within the block should be equally divided among
	 * the sub blocks.
	 * @param numAlongLatLon
	 * @param numAlongDepth
	 * @return
	 */
	public ArrayList<EqksInGeoBlock> getSubBlocksOld(int numAlongLatLon, int numAlongDepth, FaultSystemSolutionPoissonERF erf) {
		ArrayList<EqksInGeoBlock> subBlocks = new ArrayList<EqksInGeoBlock>();
		double forecastDuration = erf.getTimeSpan().getDuration();
		int numSubBlocks = numAlongLatLon*numAlongLatLon*numAlongDepth;
		for(int latSlice=0; latSlice<numAlongLatLon; latSlice++)
			for(int lonSlice=0; lonSlice<numAlongLatLon; lonSlice++)
				for(int depSlice=0; depSlice<numAlongDepth; depSlice++) {
					double subMinLat = minLat+latSlice*(maxLat-minLat)/(double)numAlongLatLon;
					double subMaxLat = subMinLat+(maxLat-minLat)/(double)numAlongLatLon;
					double subMinLon = minLon+lonSlice*(maxLon-minLon)/(double)numAlongLatLon;
					double subMaxLon = subMinLon+(maxLon-minLon)/(double)numAlongLatLon;
					double subMinDepth = minDepth+depSlice*(maxDepth-minDepth)/(double)numAlongDepth;
					double subMaxDepth = subMinDepth+(maxDepth-minDepth)/(double)numAlongDepth;
					EqksInGeoBlock subBlock = new EqksInGeoBlock(subMinLat, subMaxLat, subMinLon, subMaxLon, subMinDepth, subMaxDepth);
					for(int r=0; r<getNumRupsInside();r++) {
						int nthRup = rupIndexN_List.get(r);
						ProbEqkRupture rup = erf.getNthRupture(nthRup);
						if(!rup.getRuptureSurface().isPointSurface())
							subBlock.processRupture(rup, nthRup, forecastDuration);	
						else { // assume point sources equally divided
							double rate = rup.getMeanAnnualRate(forecastDuration)/(double)numSubBlocks;
							double fracInside = 1.0/(double)numSubBlocks;
							subBlock.processRate(rate, fracInside, nthRup, rup.getMag());
						}
							
					}
					subBlocks.add(subBlock);
				}
		
		// check total rates
		double totRate=0;
		for(int b=0; b<subBlocks.size();b++) {
			EqksInGeoBlock block = subBlocks.get(b);
			double rate=block.getTotalRateInside();
			totRate += rate;
//			System.out.println("/nBlock "+b+" rate = "+rate);
		}
		double testRate2 = getTotalRateInside();
		double ratio = totRate/testRate2;
		if(ratio<0.999 || ratio>1.001) {
			System.out.println("PROBLEM: ratio="+ratio+";\ttotRate="+totRate+"\ttestRate2="+testRate2+"\n");
			// TEST
			for(int r=0; r<rupIndexN_List.size();r++) {
				double targetRate = rateInsideList.get(r);
				double summedRate=0;
				for(int b=0; b<subBlocks.size();b++) {
					EqksInGeoBlock blk = subBlocks.get(b);
					if(blk.getRateInsideList().size()>0) {
						double rate = blk.tempGetRandomEqkRupSamplerY_Val(rupIndexN_List.get(r));
						if(rate >= 0)
							summedRate+=rate;					
					}
				}
				System.out.println("\t"+r+"\t"+rupIndexN_List.get(r)+"\t"+(summedRate/targetRate));
			}

		}

//		System.out.println("/nRate Check: "+totRate+" vs "+this.getTotalRateInside());
		
		return subBlocks;
	}
	
	/**
	 * This returns the number of ruptures that nucleate inside the block
	 * @return
	 */
	public int getNumRupsInside() {
		return this.rupIndexN_List.size();
	}
	
	/**
	 * Still need to modify location if point source or set hypocenter if finite source?
	 * (or do that in what calls this?).  Should also clone the rupture and set it as an
	 * obs and/or aftershock?
	 * @return
	 */
	public ProbEqkRupture getRandomRupture() {
		int n = getRandomRuptureIndexN();
		int srcIndex = erf.getSrcIndexForNthRup(n);
		int rupIndex = erf.getRupIndexInSourceForNthRup(n);
		return erf.getRupture(srcIndex, rupIndex);
	}
	
	
	/**
	 * This returns the index N of a randomly sampled rupture.
	 * @return
	 */
	public int getRandomRuptureIndexN() {
		// make random sampler if it doesn't already exist
		getRandomSampler();
		
		int localRupIndex = randomEqkRupSampler.getRandomInt();
		return rupIndexN_List.get(localRupIndex);
	}
	
	/**
	 * This creates (if not already existent) and returns the randomEqkRupSampler
	 * @return
	 */
	public IntegerPDF_FunctionSampler getRandomSampler() {
		if(randomEqkRupSampler == null) {
			randomEqkRupSampler = new IntegerPDF_FunctionSampler(rupIndexN_List.size());
			for(int i=0;i<rupIndexN_List.size();i++) 
				randomEqkRupSampler.set(i,rateInsideList.get(i));
		}
		return randomEqkRupSampler;
	}
	
	
	/**
	 * This assigns a random hypocentral location for the passed-in rupture.  
	 * If the rupture is a point source the location is chosen at random from
	 * within the block (assuming a uniform distribution).  If the rupture has a
	 * finite surface, the hypocentral location is chosen randomly among those surface
	 * locations that are inside the block (again, a uniform distribution).
	 * @param rupture
	 * @return
	 */
	public void setRandomHypocenterLoc(EqkRupture rupture) {
		RuptureSurface rupSurface = rupture.getRuptureSurface();
		if(rupSurface.isPointSurface()) { // randomly assign a point inside the block assuming uniform distribution
			double lat = minLat + Math.random()*(maxLat-minLat);
			double lon = minLon + Math.random()*(maxLon-minLon);
			double depth = minDepth + Math.random()*(maxDepth-minDepth);
			rupture.setHypocenterLocation(new Location(lat,lon,depth));
		}
		else {
			ArrayList<Location> locsInside = new ArrayList<Location>();
			for(Location loc:rupSurface.getEvenlyDiscritizedListOfLocsOnSurface()) {
				if(isLocInside(loc)) 
					locsInside.add(loc);
			}
			// now choose a random index from the locsInside list
			double randIndex = Math.round(locsInside.size()*Math.random()-0.5);
			rupture.setHypocenterLocation(locsInside.get((int)randIndex));
		}
	}
	
	
	/**
	 * This tells whether the given location is inside the block
	 * @param loc
	 * @return
	 */
	public boolean isLocInside(Location loc) {
		if(     loc.getLatitude()>=minLat && loc.getLatitude()< maxLat && 
				loc.getLongitude()>=minLon && loc.getLongitude()< maxLon &&
				loc.getDepth()>=minDepth && loc.getDepth()<maxDepth)
			return true;
		else
			return false;
	}




	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		long startTime=System.currentTimeMillis();

		FaultSystemSolutionPoissonERF erf = new FaultSystemSolutionPoissonERF("/Users/field/ALLCAL_UCERF2.zip");
		erf.getAdjustableParameterList().getParameter(AleatoryMagAreaStdDevParam.NAME).setValue(0.12);
		erf.updateForecast();

		
		double runtime = (System.currentTimeMillis()-startTime)/1000;
		System.out.println("ERF instantiation took "+runtime+" seconds");


		startTime=System.currentTimeMillis();
		EqksInGeoBlock eqksInGeoBlock = new EqksInGeoBlock(33.7-0.05, 33.7+0.05, -116.1-0.05, 
				-116.1+0.05, 0, 16, erf);
		
		runtime = (System.currentTimeMillis()-startTime)/1000;
		System.out.println("processing ERF  took "+runtime+" seconds");
		
		eqksInGeoBlock.writeResults();	
		System.out.println("Total Rate Inside = "+eqksInGeoBlock.getTotalRateInside());
		
		
		// Do 2 sub blocks
		startTime=System.currentTimeMillis();
		eqksInGeoBlock.getSubBlocks(2,2,erf);
		runtime = (System.currentTimeMillis()-startTime)/1000;
		System.out.println("2 slices sub blocks took "+runtime+" seconds");

		// Do 4 sub blocks
		startTime=System.currentTimeMillis();
		eqksInGeoBlock.getSubBlocks(2,2,erf);
		runtime = (System.currentTimeMillis()-startTime)/1000;
		System.out.println("4 slices sub blocks took "+runtime+" seconds");

	}
	
	/** This computes the expected, normalized mag-freq dist for the block (total rate is 1.0)
	 * 
	 * @return
	 */
	public ArbIncrementalMagFreqDist getMagProbDist() {
		ArbIncrementalMagFreqDist magDist = new ArbIncrementalMagFreqDist(2.05, 8.95, 70);
		for(int j=0; j<magList.size(); j++)
			magDist.addResampledMagRate(magList.get(j), rateInsideList.get(j), true);
		magDist.scaleToCumRate(2.05, 1);
		return magDist;
	}

	
	
	public ArrayList<Double> getRateInsideList() {return rateInsideList; }
	
	public ArrayList<Double> getFractInsideList() {return fractInsideList; }
	
	public ArrayList<Integer> getRupIndexN_List() {return rupIndexN_List; }
	
	public ArrayList<Double> getMagList() {return magList; }

}
