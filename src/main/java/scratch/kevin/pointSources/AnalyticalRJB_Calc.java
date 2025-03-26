package scratch.kevin.pointSources;

import java.text.DecimalFormat;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.QuadSurface;
import org.opensha.sha.faultSurface.utils.PointSurfaceBuilder;

import com.google.common.base.Stopwatch;

public class AnalyticalRJB_Calc {

	public static void main(String[] args) {
//		double rEpi = 50d;
//		double rupLen = 100d;
		double rEpi = 100d;
		double rupLen = 100d;
		System.out.println("rEpi="+(float)rEpi+", rupLen="+(float)rupLen);
		System.out.println("Vertical:");
		DecimalFormat distDF = new DecimalFormat("0.00");
		for (double alphaDeg=0; alphaDeg<=180.01; alphaDeg +=10d) {
			double alphaRad = Math.toRadians(alphaDeg);
			double rJB = calcRJB_VerticalCentered(rEpi, rupLen, alphaRad);
			double rJBquad = calcRJB_VerticalCenteredFromQuadSurface(rEpi, rupLen, alphaRad);
			System.out.println("\t"+(int)alphaDeg+":\trJB_an="+distDF.format(rJB)+"\trJB_quad="
					+distDF.format(rJBquad)+"\t(diff="+distDF.format(rJB-rJBquad)+")");
		}
		double dip = 50d;
		double dipRad = Math.toRadians(dip);
		double width = 10d;
		System.out.println("dip="+(int)dip+", width="+(int)width);
		for (double alphaDeg=0; alphaDeg<=180.01; alphaDeg +=10d) {
			double alphaRad = Math.toRadians(alphaDeg);
			double rJB = calcRJB_DippingCentered(rEpi, rupLen, width, dipRad, alphaRad);
			double rJBquad = calcRJB_DippingCenteredFromQuadSurface(rEpi, rupLen, width, dipRad, alphaRad);
			System.out.println("\t"+(int)alphaDeg+":\trJB_an="+distDF.format(rJB)+"\trJB_quad="
					+distDF.format(rJBquad)+"\t(diff="+distDF.format(rJB-rJBquad)+")");
		}
		
		dip = 90d;
		dipRad = Math.toRadians(dip);
		System.out.println("Vertical using dipping form");
		for (double alphaDeg=0; alphaDeg<=180.01; alphaDeg +=10d) {
			double alphaRad = Math.toRadians(alphaDeg);
			double rJB = calcRJB_DippingCentered(rEpi, rupLen, width, dipRad, alphaRad);
			double rJBquad = calcRJB_DippingCenteredFromQuadSurface(rEpi, rupLen, width, dipRad, alphaRad);
			System.out.println("\t"+(int)alphaDeg+":\trJB_an="+distDF.format(rJB)+"\trJB_quad="
					+distDF.format(rJBquad)+"\t(diff="+distDF.format(rJB-rJBquad)+")");
		}
		
		System.out.println("Vertical with fractDAS=0");
		double fractDAS = 0d;
		double fractDD = 0.5;
		for (double alphaDeg=0; alphaDeg<=360.01; alphaDeg +=10d) {
			double alphaRad = Math.toRadians(alphaDeg);
			double rJB = calcRJB(rEpi, rupLen, width, dipRad, fractDAS, fractDD, alphaRad);
			double rJBquad = calcRJB_quad(rEpi, rupLen, width, dipRad, fractDAS, fractDD, alphaRad);
			System.out.println("\t"+(int)alphaDeg+":\trJB_an="+distDF.format(rJB)+"\trJB_quad="
					+distDF.format(rJBquad)+"\t(diff="+distDF.format(rJB-rJBquad)+")");
		}
		
		System.out.println("Vertical with fractDAS=0.25");
		fractDAS = 0.25;
		for (double alphaDeg=0; alphaDeg<=360.01; alphaDeg +=10d) {
			double alphaRad = Math.toRadians(alphaDeg);
			double rJB = calcRJB(rEpi, rupLen, width, dipRad, fractDAS, fractDD, alphaRad);
			double rJBquad = calcRJB_quad(rEpi, rupLen, width, dipRad, fractDAS, fractDD, alphaRad);
			System.out.println("\t"+(int)alphaDeg+":\trJB_an="+distDF.format(rJB)+"\trJB_quad="
					+distDF.format(rJBquad)+"\t(diff="+distDF.format(rJB-rJBquad)+")");
		}
		
		dip = 50d;
		dipRad = Math.toRadians(dip);
		System.out.println("Dipping with fractDAS=0.25 & fractDD=0.5");
		for (double alphaDeg=0; alphaDeg<=360.01; alphaDeg +=10d) {
			double alphaRad = Math.toRadians(alphaDeg);
			double rJB = calcRJB(rEpi, rupLen, width, dipRad, fractDAS, fractDD, alphaRad);
			double rJBquad = calcRJB_quad(rEpi, rupLen, width, dipRad, fractDAS, fractDD, alphaRad);
			System.out.println("\t"+(int)alphaDeg+":\trJB_an="+distDF.format(rJB)+"\trJB_quad="
					+distDF.format(rJBquad)+"\t(diff="+distDF.format(rJB-rJBquad)+")");
		}
		
		fractDAS = 0.5;
		fractDD = 0d;
		System.out.println("Dipping with fractDAS=0.5 & fractDD=0");
		for (double alphaDeg=0; alphaDeg<=360.01; alphaDeg +=10d) {
			double alphaRad = Math.toRadians(alphaDeg);
			double rJB = calcRJB(rEpi, rupLen, width, dipRad, fractDAS, fractDD, alphaRad);
			double rJBquad = calcRJB_quad(rEpi, rupLen, width, dipRad, fractDAS, fractDD, alphaRad);
			System.out.println("\t"+(int)alphaDeg+":\trJB_an="+distDF.format(rJB)+"\trJB_quad="
					+distDF.format(rJBquad)+"\t(diff="+distDF.format(rJB-rJBquad)+")");
		}
		
		fractDAS = 0.25;
		fractDD = 0d;
		System.out.println("Dipping with fractDAS=0.25 & fractDD=0");
		for (double alphaDeg=0; alphaDeg<=360.01; alphaDeg +=10d) {
			double alphaRad = Math.toRadians(alphaDeg);
			double rJB = calcRJB(rEpi, rupLen, width, dipRad, fractDAS, fractDD, alphaRad);
			double rJBquad = calcRJB_quad(rEpi, rupLen, width, dipRad, fractDAS, fractDD, alphaRad);
			System.out.println("\t"+(int)alphaDeg+":\trJB_an="+distDF.format(rJB)+"\trJB_quad="
					+distDF.format(rJBquad)+"\t(diff="+distDF.format(rJB-rJBquad)+")");
		}
		
		// now benchmark them both
		int numCases = 1000;
		double[] rEips = new double[numCases];
		double[] lengths = new double[numCases];
		double[] widths = new double[numCases];
		double[] dipRads = new double[numCases];
		double[] fractDASs = new double[numCases];
		double[] fractDDs = new double[numCases];
		double[] alphas = new double[numCases];
		for (int i = 0; i < numCases; i++) {
			rEips[i] = Math.random() * 400d;
			lengths[i] = Math.random() * 150d;
			widths[i] = 1 + Math.random() * 9d;
			dipRads[i] = Math.toRadians(10+Math.random()*80d);
			fractDASs[i] = Math.random();
			fractDDs[i] = Math.random();
			alphas[i] = Math.random() * 2d*Math.PI;
		}
		int numBenchmark = 100000000;
		System.out.println("Benchmarking "+numBenchmark+" cases each in "+numCases+" sets");
		Stopwatch watch = Stopwatch.createStarted();
		for (int n=0; n<numBenchmark; n++) {
			int index = n%numCases;
			calcRJB(rEips[index], lengths[index], widths[index], dipRads[index], fractDASs[index], fractDDs[index], alphas[index]);
		}
		watch.stop();
		double secsAna = watch.elapsed().toMillis()/1000d;
		System.out.println("Analytical:\t"+secsAna+" secs ("+(float)(numBenchmark/secsAna)+" /s)");
		watch = Stopwatch.createStarted();
		for (int n=0; n<numBenchmark; n++) {
			int index = n%numCases;
			calcRJB_quad(rEips[index], lengths[index], widths[index], dipRads[index], fractDASs[index], fractDDs[index], alphas[index]);
		}
		watch.stop();
		double secsQuad = watch.elapsed().toMillis()/1000d;
		System.out.println("QuadSurface:\t"+secsQuad+" secs ("+(float)(numBenchmark/secsQuad)+" /s)");
		System.out.println("Analytical is "+(float)(secsQuad/secsAna)+" x faster");
	}

	public static double calcRJB_VerticalCentered(double rEpi, double rupLength, double alphaRad) {
		// Half the total length of the fault
		double halfLen = rupLength / 2.0;

		// Trig values
		double cosA = Math.cos(alphaRad);
		double sinA = Math.sin(alphaRad);

		// Check if perpendicular lands on the segment
		// i.e., |rEpi * cos(alpha)| <= halfLen
		boolean usePerp = (Math.abs(rEpi * cosA) <= halfLen);

		if (usePerp) {
			// The perpendicular distance = rEpi * |sinA|
			// (If alpha is strictly in [0, pi], sinA >= 0, so we could skip the abs)
			return rEpi * Math.abs(sinA);
		} else {
			// Need distance to whichever endpoint is closer.
			// Endpoint E+ = (halfLen*cosA, halfLen*sinA)
			// Vector from site to E+: (rEpi - halfLen*cosA, 0 - halfLen*sinA)
			double dxPlus = rEpi - halfLen * cosA;
			double dyPlus = - halfLen * sinA;
			double distPlus = Math.sqrt(dxPlus * dxPlus + dyPlus * dyPlus);

			// Endpoint E- = (-halfLen*cosA, -halfLen*sinA)
			// Vector from site to E-: (rEpi + halfLen*cosA, 0 + halfLen*sinA)
			double dxMinus = rEpi + halfLen * cosA;
			double dyMinus = halfLen * sinA;
			double distMinus = Math.sqrt(dxMinus * dxMinus + dyMinus * dyMinus);

			return Math.min(distPlus, distMinus);
		}
	}
	
	private static double calcRJB_VerticalCenteredFromQuadSurface(double rEpi, double rupLength, double alphaRad) {
//		Location gridNode = new Location(0d, 0d);
//		Location siteLoc = LocationUtils.location(gridNode, 0d, rEpi);
//		Location traceEnd = LocationUtils.location(gridNode, alphaRad, 0.5*rupLength);
//		Location traceStart = LocationUtils.location(gridNode, alphaRad+Math.PI, 0.5*rupLength);
//		FaultTrace trace = new FaultTrace(null);
//		trace.add(traceStart);
//		trace.add(traceEnd);
//		QuadSurface quad = new QuadSurface(trace, 90d, 1d);
//		return quad.getDistanceJB(siteLoc);
		Location gridNode = new Location(0d, 0d);
		Location siteLoc = LocationUtils.location(gridNode, 0d, rEpi);
		
		PointSurfaceBuilder builder = new PointSurfaceBuilder(gridNode);
		builder.length(rupLength);
		builder.upperDepthWidthAndDip(0d, 1d, 90d);
		builder.fractionalDAS(0.5).fractionalHypocentralDepth(0.5d);
		builder.strike(Math.toDegrees(alphaRad));
		QuadSurface quad = builder.buildQuadSurface();
		return quad.getDistanceJB(siteLoc);
	}

	public static double calcRJB_DippingCentered(double rEpi, double rupLength, double rupWidth, double dipRad, double alphaRad) {
		double widthProj = Math.cos(dipRad)*rupWidth;
	    // half-dimensions of fault footprint
	    double halfL = rupLength / 2.0;
	    double halfW = widthProj / 2.0;

	    // local coordinates of the site
	    // local x-axis = (cos alpha, sin alpha)
	    // local y-axis = (-sin alpha, cos alpha)
	    double xLoc = rEpi * Math.cos(alphaRad);
	    double yLoc = -rEpi * Math.sin(alphaRad);

	    // how far outside the fault rectangle are we, along each axis?
	    double dx = Math.max(0.0, Math.abs(xLoc) - halfL);
	    double dy = Math.max(0.0, Math.abs(yLoc) - halfW);

	    // Euclidean distance outside that rectangle
	    return Math.sqrt(dx * dx + dy * dy);
	}
	
	private static double calcRJB_DippingCenteredFromQuadSurface(double rEpi, double rupLength, double rupWidth, double dipRad, double alphaRad) {
		Location gridNode = new Location(0d, 0d);
		Location siteLoc = LocationUtils.location(gridNode, 0d, rEpi);
		
		PointSurfaceBuilder builder = new PointSurfaceBuilder(gridNode);
		builder.length(rupLength);
		builder.upperDepthWidthAndDip(0d, rupWidth, Math.toDegrees(dipRad));
		builder.fractionalDAS(0.5).fractionalHypocentralDepth(0.5d);
		builder.strike(Math.toDegrees(alphaRad));
		QuadSurface quad = builder.buildQuadSurface();
		return quad.getDistanceJB(siteLoc);
	}
	
	/**
	 * Calculates an analytical Joyner-Boore distance (rJB) for a rectangular fault with the given parameters
	 * 
	 * @param rEpi epicentral distance from the site to the grid node
	 * @param rupLength length of the fault that intersects that grid node
	 * @param rupWidth down-dip width (3D) of the fault that intersects that grid node
	 * @param dipRad dip of the fault (radians)
	 * @param gridNodeFractDAS fractional distance along-strike of the rupture where the grid node lies. A value of 0.5
	 * indicates that the rupture is centered along-strike and extends 0.5*rupLen in either direction. A value of 0
	 * indicates that the rupture begins at the grid node and extends rupLength in the alpha direction  
	 * @param gridNodeFractDepth fractional depth of the fault below the grid node. A value of 0.5 indicates that the
	 * rupture is centered down-dip about the grid node. A value of 0 indicates that the upper edge of the rupture is
	 * directly below the grid node.
	 * @param alphaRad strike angle of the fault relative to the site (radians). A value of 0 indicates that the strike
	 * direction is from the grid node directly toward the site (parallel); a value of PI/2 indicates that the strike
	 * direction is perpendicular to the site.
	 * @return
	 */
	public static double calcRJB(double rEpi, double rupLength, double rupWidth, double dipRad,
			double gridNodeFractDAS, double gridNodeFractDepth, double alphaRad) {
		// 1) Horizontal dimension of the down-dip direction
		//    (the fault extends rupWidth in 3D, so horizontally it's rupWidth*cos(dip))
		double wHorz = rupWidth * Math.cos(dipRad);

		// 2) Fault rectangle in local (strike,dip) coords is:
		//       X in [Xmin, Xmax], with total length = rupLength
		//       Y in [Ymin, Ymax], with total width = wHorz
		//    where (0,0) is the grid node in local coordinates.
		double xMin = -gridNodeFractDAS * rupLength;
		double xMax = xMin + rupLength; // = (1 - gridNodeFractDAS)*rupLength
		double yMin = -gridNodeFractDepth * wHorz;
		double yMax = yMin + wHorz;     // = (1 - gridNodeFractDepth)*wHorz

		// 3) Convert the site's global coords (rEpi, 0) -> local (xLoc, yLoc)
		//    local X-axis = strike = (cos(alpha), sin(alpha))
		//    local Y-axis = dip in map = (-sin(alpha), cos(alpha))
		//    node is at (0,0), site is at (rEpi, 0).
		double cosA = Math.cos(alphaRad);
		double sinA = Math.sin(alphaRad);
		double xLoc = rEpi * cosA;   // = x*cosA + y*sinA, but site y=0
		double yLoc = -rEpi * sinA;  // = -x*sinA + y*cosA

		// 4) Distance from (xLoc, yLoc) to that axis-aligned bounding box
		double dx = 0.0;
		if (xLoc < xMin) {
			dx = xMin - xLoc;
		} else if (xLoc > xMax) {
			dx = xLoc - xMax;
		}

		double dy = 0.0;
		if (yLoc < yMin) {
			dy = yMin - yLoc;
		} else if (yLoc > yMax) {
			dy = yLoc - yMax;
		}

		return Math.sqrt(dx * dx + dy * dy);
	}
	
	public static double calcRJB_quad(double rEpi, double rupLength, double rupWidth, double dipRad,
			double gridNodeFractDAS, double gridNodeFractDepth, double alphaRad) {
		Location gridNode = new Location(0d, 0d);
		Location siteLoc = LocationUtils.location(gridNode, 0d, rEpi);
		
		PointSurfaceBuilder builder = new PointSurfaceBuilder(gridNode);
		builder.length(rupLength);
		builder.upperDepthWidthAndDip(0d, rupWidth, Math.toDegrees(dipRad));
		builder.fractionalDAS(gridNodeFractDAS).fractionalHypocentralDepth(gridNodeFractDepth);
		builder.strike(Math.toDegrees(alphaRad));
		QuadSurface quad = builder.buildQuadSurface();
		return quad.getDistanceJB(siteLoc);
	}

}
