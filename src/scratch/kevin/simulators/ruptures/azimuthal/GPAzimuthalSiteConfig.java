package scratch.kevin.simulators.ruptures.azimuthal;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.Somerville_2006_MagAreaRel;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.earthquake.FocalMechanism;

import scratch.kevin.bbp.BBP_SourceFile;
import scratch.kevin.bbp.BBP_SourceFile.BBP_PlanarSurface;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;

public class GPAzimuthalSiteConfig extends AzimuthalSiteConfig<Integer> {
	
	private double length;
	private double width;
	private double mag;
	
	private List<BBP_SourceFile> sources;
	private List<Double> hypoDASs;
	private List<Double> hypoDDWs;
	
	private List<Location> locs;
	
	public GPAzimuthalSiteConfig(Scenario scenario, int numRuptures, double elemArea,
			double buffer, double spacing, boolean positiveX) {
		super(scenario, getIDs(numRuptures));
		
		mag = scenario.getMagnitude();
		double rake;
		switch (scenario.getFaultStyle()) {
		case NORMAL:
			rake = -90;
			break;
		case REVERSE:
			rake = 90;
			break;
		case STRIKE_SLIP:
			rake = 180;
			break;

		default:
			throw new IllegalStateException();
		}
		
		WC1994_MagLengthRelationship wc = new WC1994_MagLengthRelationship();
		length = wc.getMedianLength(mag, rake);
		System.out.println("WC 1994 Length for M="+(float)mag+": "+(float)length+" km");
		Somerville_2006_MagAreaRel magArea = new Somerville_2006_MagAreaRel();
		double area = magArea.getMedianArea(mag);
		System.out.println("Somerville 2006 Area for M="+(float)mag+": "+(float)area+" km");
		width = area/length;
		System.out.println("Calculated wdith for M="+(float)mag+": "+(float)width+" km");
		
		double sqrtElemArea = Math.sqrt(elemArea);
		
		FocalMechanism mech = new FocalMechanism(0d, scenario.getDip(), rake);
		
		Random r = new Random((long)area*(1l + numRuptures));

		sources = new ArrayList<>();
		hypoDASs = new ArrayList<>();
		hypoDDWs = new ArrayList<>();
		
		Location topCenter = new Location(34d, -118d);
		
		for (int i=0; i<numRuptures; i++) {
			BBP_PlanarSurface surface = new BBP_PlanarSurface(topCenter, length, width, mech);
			double hypoAlongStrike = length*(r.nextDouble()-0.5); // 0 is center
			double hypoDownDip = width*r.nextDouble(); // TODO
			BBP_SourceFile src = new BBP_SourceFile(surface, mag, hypoAlongStrike, hypoDownDip,
					sqrtElemArea, sqrtElemArea, 0.15, 1000+i);
			sources.add(src);
			hypoDASs.add(hypoAlongStrike + 0.5*length);
			hypoDDWs.add(hypoDownDip);
		}
		
		init(spacing, buffer, length, positiveX);
		
		locs = new ArrayList<>();
		
		Location[] rect = sources.get(0).getSurface().getRectangle();
		Location p1 = rect[0];
		Location p2 = rect[1];
		
		EvenlyDiscrXYZ_DataSet gc2xyz = getGC2XYZ();
		for (int i=0; i<gc2xyz.size(); i++) {
			Point2D pt = gc2xyz.getPoint(i);
			locs.add(SimpleGC2SiteLocationCalc.gc2ToLoc(p1, p2, pt.getX(), pt.getY()));
		}
	}
	
	private static List<Integer> getIDs(int numRuptures) {
		List<Integer> ids = new ArrayList<>();
		for (int i=0; i<numRuptures; i++)
			ids.add(i);
		return ids;
	}

	@Override
	public List<Location> getRuptureSiteLocs(Integer rupture) {
		return locs;
	}

	@Override
	public double getHypocenterDAS(Integer rupture) {
		return hypoDASs.get(rupture);
	}

	@Override
	public double getHypocenterDDW(Integer rupture) {
		return hypoDDWs.get(rupture);
	}

	@Override
	public double getLength(Integer rupture) {
		return length;
	}

	@Override
	public double getWidth(Integer rupture) {
		return width;
	}

	@Override
	public double getDip(Integer rupture) {
		return getScenario().getDip();
	}

	@Override
	public double getMag(Integer rupture) {
		return mag;
	}

	@Override
	public int getID(Integer rupture) {
		return rupture;
	}
	
	public BBP_SourceFile getBBP_Source(Integer rupture) {
		return sources.get(rupture);
	}

}
