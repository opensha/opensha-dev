package scratch.kevin.pointSources;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.Precision;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.faultSurface.FiniteApproxPointSurface;
import org.opensha.sha.faultSurface.QuadSurface;
import org.opensha.sha.faultSurface.utils.PointSurfaceBuilder;

import com.google.common.base.Preconditions;

public class AnalyticalRRupCalc {

	public static void main(String[] args) {
		Location gridLoc = new Location(0d, 0d);
		double zTop = 1d;
		double ddw = 10d;
		double dip = 15d;
		double dipRad = Math.toRadians(dip);
		double length = 10d;
		PointSurfaceBuilder builder = new PointSurfaceBuilder(gridLoc);
		builder.upperDepthWidthAndDip(zTop, ddw, dip);
		builder.length(length);
		QuadSurface quad = builder.buildQuadSurface(0d);
		double zBot = quad.getEvenlyDiscritizedLowerEdge().first().depth;
		double horzWidth = Math.cos(dipRad)*ddw;
//		FiniteApproxPointSurface hwSurf = builder.footwall(false).buildFiniteApproxPointSurface();
//		FiniteApproxPointSurface fwSurf = builder.footwall(true).buildFiniteApproxPointSurface();
		
		EvenlyDiscretizedFunc rJBs = new EvenlyDiscretizedFunc(0d, 11, 10d);
		
		EvenlyDiscretizedFunc alphaDiscr = HistogramFunction.getEncompassingHistogram(0d, 89.9d, 1d);
//		System.out.println("Alphas:\n"+alphaDiscr);
		
		for (int r=0; r<rJBs.size(); r++) {
			double rJB = rJBs.getX(r);
			System.out.println("rJB="+(int)rJB+" km");
			
			double[] rRups = new double[alphaDiscr.size()];
			for (int a=0; a<alphaDiscr.size(); a++) {
				double alphaRad = Math.toRadians(alphaDiscr.getX(a));
				
				Location siteLoc;
				if (rJB < 0.0001) {
					// special case, pick a random location on the 2D projection of the surface
					double randAlong = (Math.random()-0.5)*length;
					double randDown = (Math.random()-0.5)*horzWidth;
					siteLoc = LocationUtils.location(gridLoc, 0d, randAlong);
					siteLoc = LocationUtils.location(siteLoc, Math.PI*0.5, randDown);
					double calcRJB = quad.getDistanceJB(siteLoc);
					Preconditions.checkState(calcRJB<0.001, "rJB calc failed, calcRJB=%s, rJB=%s, randAlong=%s, randDown=%s",
							calcRJB, rJB, randAlong, randDown);
				} else {
					// figure out the distance to go
					double rEpi = rJB;
					siteLoc = LocationUtils.location(gridLoc, alphaRad, rEpi);
					double calcRJB = Double.NaN;
					for (int i=0; i<10; i++) {
						calcRJB = quad.getDistanceJB(siteLoc);
						double delta = calcRJB - rJB;
//						System.out.println("Iter "+i+" with rEpi="+(float)rEpi+": calcRJB="+(float)calcRJB+", delta="+(float)delta);
						rEpi -= delta;
						siteLoc = LocationUtils.location(gridLoc, alphaRad, rEpi);
					}
					Preconditions.checkState(Precision.equals(rJB, calcRJB, 0.001));
				}
				rRups[a] = quad.getDistanceRup(siteLoc);
			}
			double minRrup = StatUtils.min(rRups);
			double maxRrup = StatUtils.max(rRups);
			double avgRrup = StatUtils.mean(rRups);
			System.out.println("\tQuad Rrup mean="+(float)avgRrup+", range=["+(float)minRrup+", "+(float)maxRrup+")");
			double fwRrup = FiniteApproxPointSurface.getCorrDistRup(rJB, zTop, zBot, dipRad, length, horzWidth, true);
			double hwRrup = FiniteApproxPointSurface.getCorrDistRup(rJB, zTop, zBot, dipRad, length, horzWidth, false);
			System.out.println("\tfw="+(float)fwRrup+"\thw="+(float)hwRrup);
		}
		
//		for (double rEpi=0d; rEpi<101d; rEpi += 10d) {
//			System.out.println((int)rEpi+" km");
//			// figure out rJB range
//			double[] rJBs = new double[alphaDiscr.size()];
//			Location[] siteLocs = new Location[alphas.length];
//			for (int i=0; i<alphas.length; i++) {
//				siteLocs[i] = LocationUtils.location(gridLoc, Math.toRadians(alphas[i]), rEpi);
//				rJBs[i] = quad.getDistanceJB(siteLocs[i]);
//			}
//			
//			for (int i=0; i<alphas.length; i++) {
//				System.out.println("\talpha="+(int)alpha+"\trJB="+(float)rJB);
//				double rRup = quad.getDistanceRup(siteLoc);
//				double fwRrup = FiniteApproxPointSurface.getCorrDistRup(rJB, zTop, zBot, dipRad, length, horzWidth, true);
//				double hwRrup = FiniteApproxPointSurface.getCorrDistRup(rJB, zTop, zBot, dipRad, length, horzWidth, false);
//				System.out.println("\t\trRup_quad="+(float)rRup+"\tfw="+(float)fwRrup+"\thw="+(float)hwRrup);
//			}
//			
//			for (double alpha=0d; alpha<90.01; alpha += 10d) {
//				Location siteLoc = LocationUtils.location(gridLoc, Math.toRadians(alpha), rEpi);
//				double rJB = quad.getDistanceJB(siteLoc);
//				System.out.println("\talpha="+(int)alpha+"\trJB="+(float)rJB);
//				double rRup = quad.getDistanceRup(siteLoc);
//				double fwRrup = FiniteApproxPointSurface.getCorrDistRup(rJB, zTop, zBot, dipRad, length, horzWidth, true);
//				double hwRrup = FiniteApproxPointSurface.getCorrDistRup(rJB, zTop, zBot, dipRad, length, horzWidth, false);
//				System.out.println("\t\trRup_quad="+(float)rRup+"\tfw="+(float)fwRrup+"\thw="+(float)hwRrup);
//			}
//		}
	}

	public static double calcRrup(double rJB, double rupLength, double rupWidth, double dipRad,
			double zTop, double zBot, double gridNodeFractDAS, double gridNodeFractDepth, double alphaRad) {
		// clamp alpha to [0, pi]
		alphaRad = alphaRad % Math.PI;
		
		double horzWidth = rupWidth * Math.cos(dipRad);

		// special case: if rJB=0 => site is above x=0
		double xSite;
		if (rJB == 0.0) {
			xSite = 0.0;
		} else {
			// Decide footwall/hanging-wall from alpha
			// For example:
			if (alphaRad <= Math.PI/2) {
				// hanging wall
				xSite = horzWidth + rJB;
			} else {
				// footwall
				xSite = -rJB;
			}
		}

		// line segment from (0,zTop)->(horzWidth,zBot)
		double dx = horzWidth;
		double dz = (zBot - zTop);

		// vector from top edge -> site
		double vx = xSite;     // site.x - 0
		double vz = 0.0 - zTop;

		// direction of segment
		double segLenSq = dx*dx + dz*dz;
		double dot       = vx*dx + vz*dz;
		double tStar     = dot / segLenSq;

		// clamp t to [0,1]
		if (tStar < 0.0) tStar = 0.0;
		if (tStar > 1.0) tStar = 1.0;

		// nearest point
		double xFault = dx * tStar;
		double zFault = zTop + dz * tStar;

		double dxSite = xSite - xFault;
		double dzSite = 0.0   - zFault;

		return Math.sqrt(dxSite*dxSite + dzSite*dzSite);
	}

}
