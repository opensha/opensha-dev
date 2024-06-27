package scratch.kevin.pointSources;

import java.util.Random;

import org.opensha.commons.calc.magScalingRelations.MagAreaRelationship;
import org.opensha.commons.calc.magScalingRelations.MagLengthRelationship;
import org.opensha.commons.calc.magScalingRelations.MagScalingRelationship;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.util.FaultUtils;
import org.opensha.sha.faultSurface.EvenlyGriddedSurface;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.FrankelGriddedSurface;
import org.opensha.sha.faultSurface.QuadSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.utils.PtSrcDistCorr;

import com.google.common.base.Preconditions;
import com.google.common.collect.Range;

public class PointSurfaceBuilder {
	
	// required/set on input
	private Location loc;
	private double zTop;
	private double zBot;
	
	// optional
	private Region sampleFromCell = null;
	private Random rand;
	private double mag = Double.NaN; 
	private double strike = Double.NaN;
	private double dip = 90d;
	private double length = Double.NaN;
	private boolean footwall = true;
	private boolean traceCentered = true;
	
	private MagScalingRelationship scale = WC94;
	private double gridSpacing = 1d;
	
	// calculated on the fly
	private double width = Double.NaN;
	
	private static final WC1994_MagLengthRelationship WC94 = new WC1994_MagLengthRelationship();

	public PointSurfaceBuilder(Location loc) {
		this.loc = loc;
		zTop = loc.getDepth();
		zBot = loc.getDepth();
	}
	
	public PointSurfaceBuilder sampleFromCell(Region cell) {
		this.sampleFromCell = cell;
		return this;
	}
	
	public PointSurfaceBuilder random(Random rand) {
		this.rand = rand;
		return this;
	}
	
	private Random getRand() {
		if (rand == null)
			rand = new Random(Double.doubleToLongBits(loc.lat) + Double.doubleToLongBits(loc.lon));
		return rand;
	}
	
	private Location getLoc() {
		if (sampleFromCell != null) {
			double minLat = sampleFromCell.getMinLat();
			double maxLat = sampleFromCell.getMaxLat();
			double latSpan = maxLat - minLat;
			Preconditions.checkState(latSpan > 0d);
			double minLon = sampleFromCell.getMinLon();
			double maxLon = sampleFromCell.getMaxLon();
			double lonSpan = maxLon - minLon;
			Preconditions.checkState(lonSpan > 0d);
			Random rand = getRand();
			boolean rectangular = sampleFromCell.isRectangular();
			int maxNumTries = 100;
			int tries = 0;
			while (true) {
				double lat = minLat + rand.nextDouble(latSpan);
				double lon = minLon + rand.nextDouble(lonSpan);
				Location randLoc = new Location(lat, lon, loc.depth);
				if (rectangular || sampleFromCell.contains(randLoc))
					return randLoc;
				// outside
				tries++;
				Preconditions.checkState(tries <= maxNumTries,
						"Couldn't randomly sample a location in the grid cell after %s tries", tries);
			}
		}
		return loc;
	}
	
	public PointSurfaceBuilder dip(double dip) {
		FaultUtils.assertValidDip(dip);
		this.dip = dip;
		this.width = Double.NaN;
		return this;
	}
	
	/**
	 * Set the depth (both upper and lower to the same) in km
	 * @param depth
	 * @return
	 */
	public PointSurfaceBuilder singleDepth(double depth) {
		FaultUtils.assertValidDepth(depth);
		this.zTop = depth;
		this.zBot = depth;
		this.width = Double.NaN;
		return this;
	}
	
	/**
	 * Set the upper depth in km
	 * @param zTop
	 * @return
	 */
	public PointSurfaceBuilder upperDepth(double zTop) {
		FaultUtils.assertValidDepth(zTop);
		this.zTop = zTop;
		this.width = Double.NaN;
		return this;
	}
	
	/**
	 * Set the lower depth in km
	 * @param zBot
	 * @return
	 */
	public PointSurfaceBuilder lowerDepth(double zBot) {
		FaultUtils.assertValidDepth(zBot);
		this.zBot = zBot;
		this.width = Double.NaN;
		return this;
	}
	
	/**
	 * Sets the magnitude, used to infer the length if the length is not explicitly set
	 * @param mag
	 * @return
	 */
	public PointSurfaceBuilder magnitude(double mag) {
		this.mag = mag;
		return this;
	}
	
	/**
	 * Sets the magnitude scaling relationship used to infer lengths (if length is not explicitly set)
	 * @param scale
	 * @return
	 */
	public PointSurfaceBuilder scaling(MagScalingRelationship scale) {
		Preconditions.checkNotNull(scale);
		this.scale = scale;
		return this;
	}
	
	/**
	 * Set the length in km, or NaN to infer from magnitude and scaling
	 * @param length
	 * @return
	 */
	public PointSurfaceBuilder length(double length) {
		this.length = length;
		return this;
	}
	
	/**
	 * Sets the strike direction in decimal degrees, or NaN for no direction
	 * 
	 * @param strike
	 * @return
	 */
	public PointSurfaceBuilder strike(double strike) {
		this.strike = strike;
		return this;
	}
	
	/**
	 * Sets the footwall parameter, used with point representations
	 * @param footwall
	 * @return
	 */
	public PointSurfaceBuilder footwall(boolean footwall) {
		this.footwall = footwall;
		return this;
	}
	
	/**
	 * Sets the grid spacing to be used when a gridded surface is built
	 * @param gridSpacing
	 * @return
	 */
	public PointSurfaceBuilder gridSpacing(double gridSpacing) {
		this.gridSpacing = gridSpacing;
		return this;
	}
	
	/**
	 * Sets if the trace should be centered on the grid node (true), or if the middle of the surface should be (false).
	 * Only affects dipping ruptures.
	 * @param traceCentered
	 * @return
	 */
	public PointSurfaceBuilder traceCentered(boolean traceCentered) {
		this.traceCentered = traceCentered;
		return this;
	}
	
	private double getCalcLength() {
		if (Double.isFinite(length))
			return length;
		if (Double.isFinite(mag)) {
			// calculate from scaling relationship
			if (scale instanceof MagLengthRelationship) {
				return ((MagLengthRelationship)scale).getMedianLength(mag);
			} else {
				Preconditions.checkState(scale instanceof MagAreaRelationship);
				double area = ((MagAreaRelationship)scale).getMedianArea(mag);
				double width = getCalcWidth();
				if (width > 0)
					return area/width;
				else
					// zero width, return zero
					return 0d;
			}
		} else {
			// can't calculate, set to zero
			return 0d;
		}
	}
	
	private double getCalcWidth() {
		Preconditions.checkState(zBot >= zTop, "zBOT must be >= zTOR");
		if (Double.isNaN(width)) {
			if (dip == 90d)
				width = zBot-zTop;
			else
				width = (zBot-zTop)/Math.sin(Math.toRadians(dip));
		}
		return width;
	}
	
	/**
	 * Builds a point surface representation where rJB is calculated according to the chosen {@link PtSrcDistCorr},
	 * and other distances are calculated using the (possibly corrected) rJB, the footwall setting, and zTop/zBot/dip. 
	 * @return
	 */
	public FiniteApproxPointSurface buildPointSurface() {
		Preconditions.checkState(zBot >= zTop, "zBOT must be >= zTOR"); 
		
		double length = getCalcLength();
		
		return new FiniteApproxPointSurface(getLoc(), dip, zTop, zBot, footwall, length);
	}
	
	private FaultTrace buildTrace(double strike) {
		Preconditions.checkState(Double.isFinite(strike), "Can't build finite surface because strike=%s", strike);
		double length = getCalcLength();
		Preconditions.checkState(length > 0, "Can't build finite surface because length=%s; "
				+ "set magnitude to infer length from scaling relationship", length);
		double halfLen = 0.5*length;
		double strikeRad = Math.toRadians(strike);
		Location loc = getLoc();
		Location l0 = LocationUtils.location(loc, strikeRad-Math.PI, halfLen);
		Location l1 = LocationUtils.location(loc, strikeRad, halfLen);
		if (!traceCentered && zBot > zTop && dip < 90) {
			// translate it so that the surface is centered rather than the trace
			double horzWidth = (zBot-zTop)/Math.tan(Math.toRadians(dip));
			// move to the left (so that it dips to the right)
			double transAz = strikeRad - 0.5*Math.PI;
			l0 = LocationUtils.location(l0, transAz, 0.5*horzWidth);
			l1 = LocationUtils.location(l1, transAz, 0.5*horzWidth);
		}
		l0 = new Location(l0.lat, l0.lon, zTop);
		l1 = new Location(l1.lat, l1.lon, zTop);
		FaultTrace trace = new FaultTrace(null);
		trace.add(l0);
		trace.add(l1);
		return trace;
	}
	
	private double[] getRandStrikes(int num, Range<Double> strikeRange) {
		double[] strikes = new double[num];
		Random rand = getRand();
		if (strikeRange == null) {
			// pick a random strike as the initial orientation, then evenly space relatively to that
			double origStrike = rand.nextDouble(360d);
			double delta = 360d/(double)num;
			for (int i=0; i<num; i++)
				strikes[i] = origStrike + i*delta;
		} else {
			// randomly sample within the given range
			double lower = strikeRange.lowerEndpoint();
			double upper = strikeRange.upperEndpoint();
			double span = upper - lower;
			Preconditions.checkState(span > 0d);
			for (int i=0; i<num; i++)
				strikes[i] = lower + rand.nextDouble(span);
		}
		return strikes;
	}
	
	/**
	 * Builds a {@link QuadSurface} representation of this point surface. The strike direction must be set. This
	 * representation is very efficient with distance calculations, regardless of fault size. Even for very small
	 * surfaces (e.g., M5), it still performs slightly better than a 1km gridded surface (and it is much faster for larger
	 * surfaces).
	 * @return
	 */
	public QuadSurface buildQuadSurface()  {
		return buildQuadSurface(strike);
	}
	
	/**
	 * Builds a {@link QuadSurface} representation of this point surface using the passed in strike direction. This
	 * representation is very efficient with distance calculations, regardless of fault size. Even for very small
	 * surfaces (e.g., M5), it still performs slightly better than a 1km gridded surface (and it is much faster for larger
	 * surfaces).
	 * @param strike
	 * @return
	 */
	public QuadSurface buildQuadSurface(double strike)  {
		FaultTrace trace = buildTrace(strike);
		
		return new QuadSurface(trace, dip, getCalcWidth());
	}
	
	/**
	 * Builds the given number of random strike quad surfaces. The initial orientation will be randomly sampled, then
	 * if num>1, additional strikes will be evenly distributed.
	 * @param num
	 * @return
	 */
	public QuadSurface[] buildRandQuadSurfaces(int num) {
		return buildRandQuadSurfaces(num, null);
	}
	
	/**
	 * Builds the given number of random strike quad surfaces. If strikeRange is non null, orientations will be randomly
	 * sampled from the given range.
	 * @param num
	 * @param strikeRange
	 * @return
	 */
	public QuadSurface[] buildRandQuadSurfaces(int num, Range<Double> strikeRange) {
		QuadSurface[] ret = new QuadSurface[num];
		double[] strikes = getRandStrikes(num, strikeRange);
		for (int i=0; i<num; i++)
			ret[i] = buildQuadSurface(strikes[i]);
		return ret;
	}
	
	/**
	 * Builds a gridded surface representation. Distance calculations will always performs worse than
	 * {@link #buildQuadSurface()}, so use this only if you actually need a gridded surface.
	 * @return
	 */
	public EvenlyGriddedSurface buildGriddedSurface() {
		return buildGriddedSurface(strike);
	}
	
	/**
	 * Builds a gridded surface representation. Distance calculations will always performs worse than
	 * {@link #buildQuadSurface()}, so use this only if you actually need a gridded surface.
	 * @return
	 */
	public EvenlyGriddedSurface buildGriddedSurface(double strike) {
		Preconditions.checkState(zBot >= zTop, "zBOT must be >= zTOR"); 
		FaultTrace trace = buildTrace(strike);
		
		return new FrankelGriddedSurface(trace, dip, zTop, zBot, gridSpacing);
	}
	
	/**
	 * Builds the given number of random strike gridded surfaces. The initial orientation will be randomly sampled, then
	 * if num>1, additional strikes will be evenly distributed.
	 * @param num
	 * @return
	 */
	public EvenlyGriddedSurface[] buildRandGriddedSurfaces(int num) {
		return buildRandGriddedSurfaces(num, null);
	}
	
	/**
	 * Builds the given number of random strike gridded surfaces. If strikeRange is non null, orientations will be randomly
	 * sampled from the given range.
	 * @param num
	 * @param strikeRange
	 * @return
	 */
	public EvenlyGriddedSurface[] buildRandGriddedSurfaces(int num, Range<Double> strikeRange) {
		EvenlyGriddedSurface[] ret = new EvenlyGriddedSurface[num];
		double[] strikes = getRandStrikes(num, strikeRange);
		for (int i=0; i<num; i++)
			ret[i] = buildGriddedSurface(strikes[i]);
		return ret;
	}
	
	/**
	 * Builds a surface for the given inputs. This returns {@link #buildQuadSurface()} if the strike direction has
	 * been set, and {@link #buildPointSurface()} otherwise.
	 * @return
	 */
	public RuptureSurface build() {
		if (Double.isFinite(strike))
			return buildQuadSurface();
		return buildPointSurface();
	}
	
	public static void main(String[] args) {
		Location center = new Location(0d, 0d);
		PointSurfaceBuilder builder = new PointSurfaceBuilder(center);
		builder.magnitude(7.05d);
		builder.upperDepth(1d);
		builder.lowerDepth(14d);
		builder.dip(90d);
		builder.strike(0d);
		QuadSurface surf = builder.buildQuadSurface();
		
		System.out.println("Quad rJB at colocated point: "+surf.getDistanceJB(center));
		System.out.println("Trace:\t"+surf.getUpperEdge());
	}
}