package scratch.kevin.pointSources.paperFigs2026;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CompletableFuture;

import org.jfree.data.Range;
import org.opensha.commons.calc.magScalingRelations.MagLengthRelationship;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.Leonard2010_MagLengthRelationship;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.faultSurface.RectangularSurface;
import org.opensha.sha.faultSurface.utils.PointSurfaceBuilder;
import org.opensha.sha.util.FocalMech;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;

public class HangingWallFractFigure {

	public static void main(String[] args) throws IOException {
		File outputDir = new File(ConstantsAndSettings.FIGURES_DIR, "hw_fract");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		EvenlyDiscretizedFunc distances = new EvenlyDiscretizedFunc(0d, 100d, 101);
		
		double mag = 7.05d;
		FocalMech mech = FocalMech.REVERSE;
		MagLengthRelationship ml = Leonard2010_MagLengthRelationship.DIP_SLIP;
		String prefix = "m7_rev_hw_fract";
		String label = "M7.05, Reverse";

		double rake = mech.rake();
		double dip = mech.dip();
		double upperDepth = 1d;
		double lowerDepth = 14d;
		double length = ml.getMedianLength(mag);
		
		int numCenteredCalcSurfs = 3600;
		int numUncenteredCalcSurfs = numCenteredCalcSurfs*10;
		
		Location gridLoc = new Location(0d, 0d);
		
		PointSurfaceBuilder surfBuilder = new PointSurfaceBuilder(gridLoc)
				.magnitude(mag).dip(dip).upperDepth(upperDepth).lowerDepth(lowerDepth).length(length);
		
		surfBuilder.strike(0d);
		surfBuilder.fractionalDAS(0.5d);
		surfBuilder.fractionalHypocentralDepth(0.5d);
		PointSurface ptSurf = surfBuilder.buildPointSurface();
		RectangularSurface[] centeredCalcSurfs = surfBuilder.buildRandRectSurfaces(numCenteredCalcSurfs);
		
		surfBuilder.sampleDASs().sampleHypocentralDepths();
		RectangularSurface[] uncenteredCalcSurfs = surfBuilder.buildRandRectSurfaces(numUncenteredCalcSurfs);
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		EvenlyDiscretizedFunc centeredFunc = calc(gridLoc, distances, centeredCalcSurfs);
		centeredFunc.setName("Centered");
		funcs.add(centeredFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Colors.tab_blue));
		
		EvenlyDiscretizedFunc uncenteredFunc = calc(gridLoc, distances, uncenteredCalcSurfs);
		uncenteredFunc.setName("Uncentered");
		funcs.add(uncenteredFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
		
		PlotSpec plot = new PlotSpec(funcs, chars, label, "Grid Node Distance (km)", "Fraction on Hanging Wall");
		plot.setLegendInset(true);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(plot, false, false, new Range(0d, distances.getMaxX()), new Range(0.5, 1d));
		
		PlotUtils.writePlots(outputDir, prefix, gp, 800, 600, true, true, false);
	}
	
	private static EvenlyDiscretizedFunc calc(Location gridLoc, EvenlyDiscretizedFunc distances, RectangularSurface[] surfs) {
		List<CompletableFuture<Double>> futures = new ArrayList<>(distances.size());
		for (int d=0; d<distances.size(); d++) {
			double distance = distances.getX(d);
			
			futures.add(CompletableFuture.supplyAsync(()->{
				Location siteLoc = LocationUtils.location(gridLoc, 0d, distance);
				
				int numHW = 0;
				for (RectangularSurface surf : surfs) {
					if (surf.getDistanceX(siteLoc) >= 0d)
						numHW++;
				}
				return (double)numHW/(double)surfs.length;
			}));
		}
		
		EvenlyDiscretizedFunc ret = new EvenlyDiscretizedFunc(distances.getMinX(), distances.size(), distances.getDelta());
		
		for (int i=0; i<distances.size(); i++)
			ret.set(i, futures.get(i).join());
		
		return ret;
	}

}
