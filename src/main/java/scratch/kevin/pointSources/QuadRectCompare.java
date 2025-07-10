package scratch.kevin.pointSources;

import org.apache.commons.math3.util.Precision;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.faultSurface.EvenlyGriddedSurface;
import org.opensha.sha.faultSurface.QuadSurface;
import org.opensha.sha.faultSurface.RectangularSurface;
import org.opensha.sha.faultSurface.utils.PointSurfaceBuilder;

import com.google.common.base.Preconditions;

public class QuadRectCompare {

	public static void main(String[] args) {
		double[] zTops = {0d, 1d, 3d, 5d };
		double[] zBots = {6d, 10d};
		double[] dips = {90d, 60d, 30d};
		double[] rEpis = {0d, 10d, 50d, 100d, 200d, 300d};
		double[] lengths = {5d, 10d, 50d, 100d};
		double[] alphas = { 0d, 30d, 60d, 90d, 120d, 150d, 180d, 210d, 240d, 270d, 300d, 330d };
		double[] betas = { 0d, 30d, 60d, 90d, 120d, 150d, 180d, 210d, 240d, 270d, 300d, 330d };
		Location[] gridLocs = {
				new Location(0d, 0d),
				new Location(20d, 10d),
				new Location(40d, 20d)
		};
		
		double quadRectPrecision = 0.1;
		double quadGridPrecision = 1d;
		
		boolean doJB = true;
		boolean doRup = true;
		boolean doX = true;

		for (Location gridNode : gridLocs) {
			PointSurfaceBuilder builder = new PointSurfaceBuilder(gridNode);
			for (double length : lengths) {
				builder.length(length);
				for (double zTop : zTops) {
					for (double zBot : zBots) {
						builder.upperDepth(zTop).lowerDepth(zBot);
						for (double dip : dips) {
							builder.dip(dip);
							
							for (double alpha : alphas) {
								builder.strike(alpha);
								
								QuadSurface quad = builder.buildQuadSurface();
								RectangularSurface rect = builder.buildRectSurface();
								EvenlyGriddedSurface grid = builder.gridSpacing(0.1).buildGriddedSurface();
								
								String surfDescr = "len=" + (float) length + ", zTop=" + (float) zTop + ", zBot="
										+ (float) zBot + ", dip=" + (float) dip + ", alpha=" + (float) alpha;
								System.out.println("Testing "+surfDescr);
								
								for (double beta : betas) {
									for (double rEpi : rEpis) {
										String str = "beta="+(float)beta+", rEpi="+(float)rEpi;
										Location siteLoc = LocationUtils.location(gridNode, Math.toRadians(beta), rEpi);
										
										if (doJB) {
											double rJBquad = quad.getDistanceJB(siteLoc);
											double rJBrect = rect.getDistanceJB(siteLoc);
											double rJBgrid = grid.getDistanceJB(siteLoc);
											Preconditions.checkState(Precision.equals(rJBquad, rJBgrid, quadGridPrecision),
													"%s; JBquad=%s, JBgrid=%s", str, rJBquad, rJBgrid);
											Preconditions.checkState(Precision.equals(rJBquad, rJBrect, quadRectPrecision),
													"%s; JBquad=%s, JBrect=%s", str, rJBquad, rJBgrid);
										}
										if (doRup) {
											double rRupquad = quad.getDistanceRup(siteLoc);
											double rRuprect = rect.getDistanceRup(siteLoc);
											double rRupgrid = grid.getDistanceRup(siteLoc);
											Preconditions.checkState(
													Precision.equals(rRupquad, rRupgrid, quadGridPrecision),
													"%s; RUPquad=%s, RUPgrid=%s", str, rRupquad, rRupgrid);
											Preconditions.checkState(
													Precision.equals(rRupquad, rRuprect, quadRectPrecision),
													"%s; RUPquad=%s, RUPrect=%s", str, rRupquad, rRuprect);
										}
										if (doX) {
											double xQuad = quad.getDistanceX(siteLoc);
											double xRect = rect.getDistanceX(siteLoc);
											double xGrid = grid.getDistanceX(siteLoc);
											double refX = Math.max(Math.abs(xQuad), Math.abs(xRect));
											double refDist = rEpi-0.5*length;
											double refPrecision = Math.max(1d, refDist*2e-2);
											boolean approxZero = refX < refPrecision;
											Preconditions.checkState(Precision.equals(xQuad, xGrid, refPrecision)
													&& (approxZero || (xQuad >= 0) == (xGrid >= 0)),
//													"%s; Xquad=%s, Xgrid=%s", str, xQuad, xGrid);
													"%s; Xquad=%s, Xgrid=%s, refX=%s, refDist=%s, refPrecision=%s", str, xQuad, xGrid, refX, refDist, refPrecision);
											Preconditions.checkState(Precision.equals(xQuad, xRect, refPrecision)
													&& (approxZero || (xQuad >= 0) == (xRect >= 0)),
//													"%s; Xquad=%s, Xrect=%s", str, xQuad, xRect);
													"%s; Xquad=%s, Xrect=%s, refX=%s, refDist=%s, refPrecision=%s", str, xQuad, xRect, refX, refDist, refPrecision);
										}
									}
								}
							}
						}
					}
				}
			}
		}
		System.out.println("All pass!");
	}

}
