package scratch.ned.ETAS_Tests;

import java.util.ArrayList;
import java.util.ListIterator;

import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.UCERF2;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;
import org.opensha.sha.faultSurface.AbstractEvenlyGriddedSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.sha.magdist.ArbIncrementalMagFreqDist;

import scratch.UCERF3.erf.ETAS.IntegerPDF_FunctionSampler;

/**
 * This class store information about all ruptures that nucleate inside this geographic block.
 * @author field
 *
 */
public class EqksInGeoBlock {
	
	double minLat, maxLat, minLon,maxLon, minDepth, maxDepth; // the dimensions of the block
	AbstractERF erf;					// reference to the erf this block is used with
	
	// these lists hold information about each rupture that nucleates in this block
	ArrayList<Integer> srcIndexList;	// this stores the source index for each rupture that nucleates inside the block
	ArrayList<Integer> rupIndexList;	// this stores the ruptures index (inside the source)
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
		srcIndexList = new ArrayList<Integer>();
		rupIndexList = new ArrayList<Integer>();
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
			double minDepth, double maxDepth, AbstractERF erf) {
		
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
	 * srcIndex, rupIndex, rateInside, fractInside, randomEqkRupSampler.getY(i), and srcName (if erf given)
	 */
	public void writeResults() {
		// do the following to make sure randomEqkRupSampler has been created
		getRandomRuptureIndices();
		for(int i=0; i<srcIndexList.size();i++) {
			System.out.print(i+"\t"+srcIndexList.get(i)+"\t"+rupIndexList.get(i)+"\t"+
					rateInsideList.get(i)+"\t"+fractInsideList.get(i)+"\t"+
					magList.get(i)+"\t"+randomEqkRupSampler.getY(i));
			if(erf != null)
				System.out.print("\t"+erf.getSource(srcIndexList.get(i)).getName()+"\n");
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
		return Math.pow(getBlockVolume(),1/3);
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
			int numRups = erf.getNumRuptures(s);
			for(int r=0; r<numRups;r++) {
				ProbEqkRupture rup = erf.getRupture(s, r);
				processRupture(rup, s, r, forecastDuration);
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
	public void processRupture(ProbEqkRupture rup, int srcIndex, int rupIndex, double forecastDuration) {
		double rate = rup.getMeanAnnualRate(forecastDuration);
		if(rate > 0) {
			LocationList locList = rup.getRuptureSurface().getEvenlyDiscritizedListOfLocsOnSurface();
			int numLoc = locList.size();
			int numLocInside = 0;
			for(Location loc : locList)
				if(isLocInside(loc)) numLocInside+=1;
			if(numLocInside>0) {
				fractInsideList.add((double)numLocInside/(double)numLoc);
				rateInsideList.add((double)rate*(double)numLocInside/(double)numLoc);
				srcIndexList.add(srcIndex);
				rupIndexList.add(rupIndex);
				magList.add(rup.getMag());
			}			
		}
	}
	
	
	/**
	 * This adds to the given info to the lists (e.g., if determined externally)
	 * @param rate
	 * @param fracInside
	 * @param srcIndex
	 * @param rupIndex
	 * @param magnitude
	 */
	public void processRate(double rate, double fracInside, int srcIndex, int rupIndex, double magnitude) {
		if(rate > 0) {
			fractInsideList.add(fracInside);
			rateInsideList.add(rate);
			srcIndexList.add(srcIndex);
			rupIndexList.add(rupIndex);
			magList.add(magnitude);
		}			
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
	public ArrayList<EqksInGeoBlock> getSubBlocks(int numAlongLatLon, int numAlongDepth, AbstractERF erf) {
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
						int iSrc = srcIndexList.get(r);
						int iRup = rupIndexList.get(r);
						ProbEqkRupture rup = erf.getRupture(iSrc, iRup);
						if(!rup.getRuptureSurface().isPointSurface())
							subBlock.processRupture(rup, iSrc, iRup, forecastDuration);	
						else { // assume point sources equally divided
							double rate = rup.getMeanAnnualRate(forecastDuration)/numSubBlocks;
							double fracInside = 1/numSubBlocks;
							subBlock.processRate(rate, fracInside, iSrc, iRup, rup.getMag());
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
//		System.out.println("/nRate Check: "+totRate+" vs "+this.getTotalRateInside());
		
		return subBlocks;
	}
	
	/**
	 * This returns the number of ruptures that nucleate inside the block
	 * @return
	 */
	public int getNumRupsInside() {
		return this.srcIndexList.size();
	}
	
	
	/**
	 * Still need to modify location if point source or set hypocenter if finite source?
	 * (or do that in what calls this?).  Should also clone the rupture and set it as an
	 * obs and/or aftershock?
	 * @return
	 */
	public ProbEqkRupture getRandomRupture() {
		// make random sampler if it doesn't already exist
		getRandomSampler();
		
		int localRupIndex = randomEqkRupSampler.getRandomInt();
		int iSrc=srcIndexList.get(localRupIndex);
		int iRup = rupIndexList.get(localRupIndex);
		ProbEqkRupture rup = erf.getRupture(iSrc, iRup);
//		rup.setRuptureIndexAndSourceInfo(iSrc, "ETAS Event", iRup);
		return rup;
	}
	
	
	/**
	 * This returns the source index (in the 0th array element) and rupture index 
	 * (in the 1st array element) for a randomly sampled rupture.
	 * @return
	 */
	public int[] getRandomRuptureIndices() {
		// make random sampler if it doesn't already exist
		if(randomEqkRupSampler == null) {
			randomEqkRupSampler = new IntegerPDF_FunctionSampler(srcIndexList.size());
			for(int i=0;i<srcIndexList.size();i++) 
				randomEqkRupSampler.set(i,rateInsideList.get(i));
		}
		
		int[] indices = new int[2];
		int localRupIndex = randomEqkRupSampler.getRandomInt();
		indices[0] = srcIndexList.get(localRupIndex);
		indices[1] = rupIndexList.get(localRupIndex);
		return indices;
	}
	
	
	/**
	 * This creates (if not already existent) and returns the randomEqkRupSampler
	 * @return
	 */
	public IntegerPDF_FunctionSampler getRandomSampler() {
		if(randomEqkRupSampler == null) {
			randomEqkRupSampler = new IntegerPDF_FunctionSampler(srcIndexList.size());
			for(int i=0;i<srcIndexList.size();i++) 
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
			ListIterator<Location> locs = rupSurface.getLocationsIterator();
			while(locs.hasNext()) {
				Location loc = locs.next();
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
	private boolean isLocInside(Location loc) {
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

		// Create UCERF2 instance
		int duration = 1;
		MeanUCERF2 meanUCERF2 = new MeanUCERF2();
		meanUCERF2.setParameter(UCERF2.RUP_OFFSET_PARAM_NAME, new Double(10.0));
		meanUCERF2.getParameter(UCERF2.PROB_MODEL_PARAM_NAME).setValue(UCERF2.PROB_MODEL_POISSON);
		meanUCERF2.setParameter(UCERF2.BACK_SEIS_NAME, UCERF2.BACK_SEIS_INCLUDE);
		meanUCERF2.setParameter(UCERF2.BACK_SEIS_RUP_NAME, UCERF2.BACK_SEIS_RUP_POINT);
		meanUCERF2.getTimeSpan().setDuration(duration);
		meanUCERF2.updateForecast();
		
		double runtime = (System.currentTimeMillis()-startTime)/1000;
		System.out.println("ERF instantiation took "+runtime+" seconds");


		startTime=System.currentTimeMillis();
		EqksInGeoBlock eqksInGeoBlock = new EqksInGeoBlock(33.7-0.05, 33.7+0.05, -116.1-0.05, 
				-116.1+0.05, 0, 16, meanUCERF2);
		
		runtime = (System.currentTimeMillis()-startTime)/1000;
		System.out.println("processing ERF  took "+runtime+" seconds");
		
		eqksInGeoBlock.writeResults();	
		System.out.println("Total Rate Inside = "+eqksInGeoBlock.getTotalRateInside());
		
		
		// Do 2 sub blocks
		startTime=System.currentTimeMillis();
		eqksInGeoBlock.getSubBlocks(2,2,meanUCERF2);
		runtime = (System.currentTimeMillis()-startTime)/1000;
		System.out.println("2 slices sub blocks took "+runtime+" seconds");

		// Do 4 sub blocks
		startTime=System.currentTimeMillis();
		eqksInGeoBlock.getSubBlocks(2,2,meanUCERF2);
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
	
	public ArrayList<Integer> getSrcIndexList() {return srcIndexList; }
	
	public ArrayList<Double> getMagList() {return magList; }

}
