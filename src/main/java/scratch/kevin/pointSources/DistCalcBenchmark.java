package scratch.kevin.pointSources;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.jfree.data.Range;
import org.opensha.commons.calc.magScalingRelations.MagLengthRelationship;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.utils.PointSurfaceBuilder;

import com.google.common.base.Stopwatch;

class DistCalcBenchmark {
	
	private enum SurfType {
		PT_SURF_FOOTWALL("Point Surface (footwall)"),
		PT_SURF_HANGING_WALL("Point Surface (hanging wall)"),
		FINITE_QUAD("Finite Quad Surface"),
		FINITE_GRIDDED("Finite Gridded Surface");
		
		private String label;

		private SurfType(String label) {
			this.label = label;
		}
		
		public RuptureSurface buildSurface(PointSurfaceBuilder builder) {
			switch (this) {
			case PT_SURF_FOOTWALL:
				builder.footwall(true);
				return builder.buildPointSurface();
			case PT_SURF_HANGING_WALL:
				builder.footwall(false);
				return builder.buildPointSurface();
			case FINITE_GRIDDED:
				return builder.buildGriddedSurface();
			case FINITE_QUAD:
				return builder.buildQuadSurface();

			default:
				throw new IllegalStateException();
			}
		}
	}

	public static void main(String[] args) throws IOException {
		PointSurfaceBuilder builder = new PointSurfaceBuilder(new Location(0d, 0d));

		EvenlyDiscretizedFunc refFunc = new EvenlyDiscretizedFunc(5d, 8.5d, 30);
		
		SurfType[] types = SurfType.values();
		EvenlyDiscretizedFunc[] timeFuncs = new EvenlyDiscretizedFunc[types.length];
		for (int i=0; i<types.length; i++) {
			timeFuncs[i] = refFunc.deepClone();
			timeFuncs[i].setName(types[i].label);
		}
		
		MagLengthRelationship WC94 = new WC1994_MagLengthRelationship();
		
		builder.dip(90d);
		builder.strike(0d);
		
		Location[] testLocs = new Location[100000];
		for (int i=0; i<testLocs.length; i++)
			testLocs[i] = new Location(2*Math.random() - 1d, 2*Math.random() - 1d);

		// start at -1 to do a burn in run first that we discard
		for (int i=-1; i<refFunc.size(); i++) {
			double mag = i >= 0 ? refFunc.getX(i) : refFunc.getMaxX();
			System.out.println("Doing M"+(mag));

			double length = WC94.getMedianLength(mag);
			double aspectWidth = length / 1.5;
			
			builder.length(length);
			builder.upperDepth(0);
			builder.lowerDepth(Math.min(aspectWidth, 14d));
			
			for (int t=0; t<types.length; t++) {
				RuptureSurface surf = types[t].buildSurface(builder);
				
				Stopwatch watch = Stopwatch.createStarted();
				for (Location loc : testLocs) {
					surf.getDistanceJB(loc);
					surf.getDistanceRup(loc);
					surf.getDistanceX(loc);
				}
				watch.stop();
				
				double secs = watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
				System.out.println(types[t].label+":\t"+(float)secs+" s");
				if (i >= 0)
					timeFuncs[t].set(i, secs);
			}
		}
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();

		CPT cpt = GMT_CPT_Files.CATEGORICAL_BATLOW_UNIFORM.instance();
		double maxY = 0d;
		for (int i=0; i<timeFuncs.length; i++) {
			maxY = Math.max(maxY, timeFuncs[i].getMaxY());
			funcs.add(timeFuncs[i]);
			Color color = cpt.get(i).minColor;
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, color));
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, " ", "Magnitude", "Time For "+testLocs.length+" Distance Calcs (s)");
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		Range xRange = new Range(refFunc.getMinX(), refFunc.getMaxX());
		Range yRange = new Range(0d, maxY*1.2);
		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		
		PlotUtils.writePlots(new File("/tmp"), "pt_src_dist_calcs", gp, 800, 800, true, false, false);
	}

}
