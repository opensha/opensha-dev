package scratch.ned.nshm23;

import java.util.ArrayList;

import org.opensha.commons.calc.magScalingRelations.MagScalingRelationship;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.sha.earthquake.FocalMechanism;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.griddedForecast.HypoMagFreqDistAtLoc;
import org.opensha.sha.earthquake.rupForecastImpl.PointToFiniteSource;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

/**
 * <p>Title: PoissonAreaSourceNSHM23 </p>
 * <p>Description: This implements the area sources used in the 2023 NSHM. </p>
 * 
 * @author Edward Field
 * @version 1.0
 */

public class AreaSourceNSHM23 extends ProbEqkSource implements java.io.Serializable{

	private static final long serialVersionUID = 1L;
	// for Debug purposes
	private static String C = new String("PoissonAreaSource");
	private static String NAME = "Poisson Area Source";
	private boolean D = false;

	int id;
	String name;
	Region region;
	double strike;
	double gridResolution;
	GriddedRegion griddedRegion;
	double[] nodeWeights;
	double maxRupLength;
	ArrayList<ProbEqkSource> gridSrcListArray;
	ArrayList<ProbEqkRupture> allRupturesList;


	/**
	 * 	
	 * 
	 * STILL NEED TO HANDLE THE strike=NaN CASE PROPERLY (PRESENTLY GETTING A SINGLE RANDOM STRIKE.  EMAIL
	 * REQUEST SENT TO PETER ON HOW TO DO THIS.
	 * @param id
	 * @param name
	 * @param region
	 * @param gridResolution
	 * @param magFreqDist
	 * @param upperSesDepth
	 * @param lowerSeisDepth
	 * @param strike
	 * @param duration
	 * @param lineSource
	 */
	public AreaSourceNSHM23(int id, String name, Region region, double gridResolution, IncrementalMagFreqDist 
			magFreqDist, double upperSeisDepth, double lowerSeisDepth, FocalMechanism focalMech,
			double duration, boolean lineSource) {
		
		this.id = id;
		this.name = name;
		NAME = name;
		this.region = region;
		griddedRegion = new GriddedRegion(region, gridResolution, new Location(0.0,0.0)); 
		computeNodeWeights(); // this accounts for cell area as a function of latitude
		
		double minMag = 0;

		this.isPoissonian = true;
				
		MagScalingRelationship magScalingRel = new WC1994_MagLengthRelationship();
		
		allRupturesList = new ArrayList<ProbEqkRupture>();
		
		gridSrcListArray = new ArrayList<ProbEqkSource>(); 
		for(int i=0;i<griddedRegion.getNumLocations();i++) {
			IncrementalMagFreqDist mfd = magFreqDist.deepClone();
			mfd.scale(nodeWeights[i]);
			
			PointToFiniteSource ptSrc;
			if(focalMech.getDip() == 90d) {
				ptSrc = new PointToFiniteSource(mfd, griddedRegion.getLocation(i), focalMech,
						upperSeisDepth, magScalingRel,lowerSeisDepth, duration, minMag, lineSource);			}
			else { // do both dip directions
				double strike2 = focalMech.getStrike()+180;
				if(strike2>360) strike2 -=360; // not sure this is needed
				FocalMechanism focalMech2 = new FocalMechanism(strike2,focalMech.getDip(),focalMech.getRake());
				mfd.scale(0.5);  // half rate for each strike
				IncrementalMagFreqDist[] magDistArray = {mfd, mfd};
                FocalMechanism[] focalMechArray = {focalMech, focalMech2};
				HypoMagFreqDistAtLoc hypoMagFreqDistAtLoc = new HypoMagFreqDistAtLoc(magDistArray,griddedRegion.getLocation(i),focalMechArray);
				ArbitrarilyDiscretizedFunc aveRupTopVersusMag = new ArbitrarilyDiscretizedFunc();
				aveRupTopVersusMag.set(0.0,upperSeisDepth); // no mag dependence
				aveRupTopVersusMag.set(10.0,upperSeisDepth);
				ptSrc = new PointToFiniteSource(hypoMagFreqDistAtLoc, aveRupTopVersusMag, magScalingRel,lowerSeisDepth, 
						duration, minMag, lineSource);
			}
			
			gridSrcListArray.add(ptSrc);
			allRupturesList.addAll(ptSrc.getRuptureList());
		}
		
		maxRupLength = PointToFiniteSource.getRupLength(magFreqDist.getMaxMagWithNonZeroRate(), upperSeisDepth, 
				lowerSeisDepth, focalMech.getDip(), magScalingRel);
	}
	
	
	/**
	 * This computes the weight for each node as one over the total number of nodes
	 * multiplied by the area of the node (which changes with lat), and then renormalized
	 */
	private void computeNodeWeights() {
		int numPts = griddedRegion.getNodeCount();
		nodeWeights = new double[numPts];
		double tot=0;
		for(int i=0;i<numPts;i++) {
			double latitude = griddedRegion.locationForIndex(i).getLatitude();
			nodeWeights[i] = Math.cos(latitude*Math.PI/180);
			tot += nodeWeights[i];
		}
		for(int i=0;i<numPts;i++) nodeWeights[i] /= (tot);
	}
		
	
	/**
	 * This makes and returns the nth probEqkRupture for this source.
	 */
	public ProbEqkRupture getRupture(int nthRupture){
		return allRupturesList.get(nthRupture);
	}

		
	public Region getRegion() { return region;}

	public GriddedRegion getGriddedRegion() { return griddedRegion;}

	/**
	 * This returns the shortest horizontal dist to the point source (minus half the length of the 
	 * longest rupture).
	 * 
	 * @param site
	 * @return minimum distance
	 */
	public double getMinDistance(Site site) {
		double dist = region.distanceToLocation(site.getLocation()) - maxRupLength/2.0;
		if(dist < 0) dist=0;
		return dist;
	}


	@Override
	public LocationList getAllSourceLocs() {
		// TODO Auto-generated method stub
		return null;
	}


	@Override
	public RuptureSurface getSourceSurface() {
		// TODO Auto-generated method stub
		return null;
	}


	@Override
	public int getNumRuptures() {
		return allRupturesList.size();
	}
	
	@Override
	public String getName() {
		return name;
	}


}

