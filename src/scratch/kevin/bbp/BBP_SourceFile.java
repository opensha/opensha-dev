package scratch.kevin.bbp;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.text.DecimalFormat;
import java.util.HashMap;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.sha.earthquake.FocalMechanism;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.QuadSurface;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

public class BBP_SourceFile {
	
	private BBP_PlanarSurface surface;
	private double mag;
	private double hypoAlongStrike;
	private double hypoDownDip;
	private double dWid;
	private double dLen;
	private double cornerFreq;
	private int seed;
	
	public static class BBP_PlanarSurface {
		private Location topCenter;
		private double length;
		private double width;
		private FocalMechanism mech;
		
		public BBP_PlanarSurface(Location topCenter, double length, double width, FocalMechanism mech) {
			this.topCenter = topCenter;
			this.length = length;
			this.width = width;
			this.mech = mech;
		}
		
		public Location getTopCenter() {
			return topCenter;
		}

		public double getLength() {
			return length;
		}

		public double getWidth() {
			return width;
		}

		public FocalMechanism getFocalMechanism() {
			return mech;
		}
	}

	public BBP_SourceFile(BBP_PlanarSurface surface, double mag, double hypoAlongStrike, double hypoDownDip,
			double dWid, double dLen, double cornerFreq, int seed) {
		this.surface = surface;
		this.mag = mag;
		this.hypoAlongStrike = hypoAlongStrike;
		this.hypoDownDip = hypoDownDip;
		this.dWid = dWid;
		this.dLen = dLen;
		this.cornerFreq = cornerFreq;
		Preconditions.checkState(seed >= 0);
		this.seed = seed;
	}
	
	public void writeToFile(File file) throws IOException {
		FileWriter fw = new FileWriter(file);
		
		FocalMechanism mech = getFocalMechanism();
		
		fw.write("MAGNITUDE = "+(float)mag+"\n");
		fw.write("FAULT_LENGTH = "+(float)surface.length+"\n");
		fw.write("FAULT_WIDTH = "+(float)surface.width+"\n");
		fw.write("DEPTH_TO_TOP = "+(float)surface.topCenter.getDepth()+"\n");
		fw.write("STRIKE = "+(float)mech.getStrike()+"\n");
		fw.write("RAKE = "+(float)mech.getRake()+"\n");
		fw.write("DIP = "+(float)mech.getDip()+"\n");
		fw.write("LAT_TOP_CENTER = "+(float)surface.topCenter.getLatitude()+"\n");
		fw.write("LON_TOP_CENTER = "+(float)surface.topCenter.getLongitude()+"\n");
		fw.write("HYPO_ALONG_STK = "+(float)hypoAlongStrike+"\n");
		fw.write("HYPO_DOWN_DIP = "+(float)hypoDownDip+"\n");
		fw.write("DWID = "+lenDF.format(dWid)+"\n");
		fw.write("DLEN = "+lenDF.format(dLen)+"\n");
		fw.write("CORNER_FREQ = "+(float)cornerFreq+"\n");
		fw.write("SEED = "+seed+"\n");
		
		fw.close();
	}
	
	private static final DecimalFormat lenDF = new DecimalFormat("0.0#");
	
	public static BBP_SourceFile readFile(File file) throws IOException {
		HashMap<String, String> vals = new HashMap<>();
		for (String line : Files.readLines(file, Charset.defaultCharset())) {
			line = line.trim();
			if (line.isEmpty() || line.startsWith("#"))
				continue;
			Preconditions.checkState(line.contains("="));
			int eqIndex = line.indexOf("=");
			String key = line.substring(0, eqIndex).trim();
			String val = line.substring(eqIndex+1).trim();
			vals.put(key, val);
		}
		
		double mag = Double.parseDouble(vals.get("MAGNITUDE"));
		double lat = Double.parseDouble(vals.get("LAT_TOP_CENTER"));
		double lon = Double.parseDouble(vals.get("LON_TOP_CENTER"));
		double length = Double.parseDouble(vals.get("FAULT_LENGTH"));
		double width = Double.parseDouble(vals.get("FAULT_WIDTH"));
		double depthTop = Double.parseDouble(vals.get("DEPTH_TO_TOP"));
		Location topCenter = new Location(lat, lon, depthTop);
		double strike = Double.parseDouble(vals.get("STRIKE"));
		double rake = Double.parseDouble(vals.get("RAKE"));
		double dip = Double.parseDouble(vals.get("DIP"));
		FocalMechanism mech = new FocalMechanism(strike, dip, rake);
		double hypoAlongStrike = Double.parseDouble(vals.get("HYPO_ALONG_STK"));
		double hypoDownDip = Double.parseDouble(vals.get("HYPO_DOWN_DIP"));
		double dWid = Double.parseDouble(vals.get("DWID"));
		double dLen = Double.parseDouble(vals.get("DLEN"));
		double cornerFreq = Double.parseDouble(vals.get("CORNER_FREQ"));
		int seed = Integer.parseInt(vals.get("SEED"));
		
		BBP_PlanarSurface surface = new BBP_PlanarSurface(topCenter, length, width, mech);
		
		return new BBP_SourceFile(surface, mag, hypoAlongStrike, hypoDownDip, dWid, dLen, cornerFreq, seed);
	}

	public int getSeed() {
		return seed;
	}

	public void setSeed(int seed) {
		this.seed = seed;
	}

	public double getMag() {
		return mag;
	}

	public FocalMechanism getFocalMechanism() {
		return surface.getFocalMechanism();
	}

	public double getHypoAlongStrike() {
		return hypoAlongStrike;
	}

	public double getHypoDownDip() {
		return hypoDownDip;
	}

	public double getdWid() {
		return dWid;
	}

	public double getdLen() {
		return dLen;
	}

	public double getCornerFreq() {
		return cornerFreq;
	}
	
	public Location[] buildRectangle() {
		Location[] rect = new Location[4];
		
		double strikeRad = Math.toRadians(surface.mech.getStrike());
		Location topRight = LocationUtils.location(surface.topCenter, strikeRad, surface.length*0.5);
		Location topLeft = LocationUtils.location(surface.topCenter, strikeRad+Math.PI, surface.length*0.5);
		
		LocationVector downDipVect = getDownDipVector(surface.width);
		Location botRight = LocationUtils.location(topRight, downDipVect);
		Location botLeft = LocationUtils.location(topLeft, downDipVect);
		
		rect[0] = topLeft;
		rect[1] = topRight;
		rect[2] = botRight;
		rect[3] = botLeft;
		
		return rect;
	}
	
	private LocationVector getDownDipVector(double ddw) {
		double downDipRad = Math.toRadians(surface.mech.getStrike()) + Math.PI/2d;
		double dipRad = Math.toRadians(surface.mech.getDip());
		double downDipHorz = Math.cos(dipRad)*ddw;
		double downDipVert = Math.sin(dipRad)*ddw;
		// az in degrees here
		return new LocationVector(Math.toDegrees(downDipRad), downDipHorz, downDipVert);
	}
	
	public Location getHypoLoc() {
		double strikeRad = Math.toRadians(surface.mech.getStrike());
		Location hypoAlongTop = LocationUtils.location(surface.topCenter, strikeRad, hypoAlongStrike);
		LocationVector downDipVect = getDownDipVector(hypoDownDip);
		return LocationUtils.location(hypoAlongTop, downDipVect);
	}
	
	public BBP_PlanarSurface getSurface() {
		return surface;
	}
	
	public double getHorizontalDistance(Location loc) {
		QuadSurface surf = getQuadSurface();
		return surf.getDistanceJB(loc);
	}
	
	public double getLinearDistance(Location loc) {
		QuadSurface surf = getQuadSurface();
		return surf.getDistanceRup(loc);
	}
	
	private QuadSurface quad = null;
	
	public synchronized QuadSurface getQuadSurface() {
		if (quad == null) {
			Location[] rect = buildRectangle();
			FaultTrace trace = new FaultTrace("");
			trace.add(rect[0]);
			trace.add(rect[1]);
			quad = new QuadSurface(trace, getFocalMechanism().getDip(), getSurface().getWidth());
		}
		return quad;
	}

}
