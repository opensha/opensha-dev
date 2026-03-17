package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.util.Precision;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationUtils.LocationAverager;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.faultSurface.CompoundSurface;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.cache.SurfaceCachingPolicy;
import org.opensha.sha.faultSurface.cache.SurfaceDistances;
import org.opensha.sha.faultSurface.cache.SurfaceCachingPolicy.CacheTypes;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;

import net.mahdilamb.colormap.Colors;

public class NewCompoundSurfaceComparisons {

	public static void main(String[] args) throws IOException {
		SurfaceCachingPolicy.force(CacheTypes.THREAD_LOCAL);
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged.zip"));
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		surfBuildTime(rupSet, 50);
		System.exit(0);
		
		GriddedRegion sitesGrid = new GriddedRegion(NSHM23_RegionLoader.loadFullConterminousWUS(), 0.25, GriddedRegion.ANCHOR_0_0);
		
		System.out.println("Have "+sitesGrid.getNodeCount()+" sites");
		
		CompoundSurface[] origSurfs = new CompoundSurface[rupSet.getNumRuptures()];
		CompoundSurface[] newSurfs = new CompoundSurface[rupSet.getNumRuptures()];
		
		File outputDir = new File("/tmp/compound_surf_test");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		int numDebug = 0;
		int maxNumDebug = 100;
		
		for (int r=0; r<rupSet.getNumRuptures(); r++) {
			if (r < 10000 && r % 1000 == 0 || r % 10000 == 0)
				System.out.println("Processing rupture "+r);
			RuptureSurface surf = rupSet.getSurfaceForRupture(r, 1d);
			Preconditions.checkState(surf instanceof CompoundSurface);
			
			newSurfs[r] = (CompoundSurface)surf;
			origSurfs[r] = new CompoundSurface.Simple(newSurfs[r].getSurfaceList()); // retains original ordering
			
			if (numDebug < maxNumDebug) {
				if (!LocationUtils.areSimilar(origSurfs[r].getFirstLocOnUpperEdge(), newSurfs[r].getFirstLocOnUpperEdge())
						|| !LocationUtils.areSimilar(origSurfs[r].getLastLocOnUpperEdge(), newSurfs[r].getLastLocOnUpperEdge())
						|| !LocationUtils.areSimilar(origSurfs[r].getFirstLocOnLowerEdge(), newSurfs[r].getFirstLocOnLowerEdge())
						|| !LocationUtils.areSimilar(origSurfs[r].getLastLocOnLowerEdge(), newSurfs[r].getLastLocOnLowerEdge())) {
					System.out.println("Ordering mismatch for "+r);
					debugSurf(outputDir, rupSet, r, origSurfs[r], newSurfs[r]);
					numDebug++;
					if (numDebug == maxNumDebug) {
						System.out.println("Bailing after "+numDebug+" debugs");
					}
				}
			}
			
			origSurfs[r] = origSurfs[r];
			newSurfs[r] = newSurfs[r];
		}
		
		System.out.println("DONE building");
		
		System.out.println("Testing distances for "+sitesGrid.getNodeCount()+" sites");
		final double tol = 1e-4;
		List<CompletableFuture<Void>> futures = new ArrayList<>();
		for (Location loc : sitesGrid.getNodeList()) {
			futures.add(CompletableFuture.runAsync(() -> {
				for (int r=0; r<origSurfs.length; r++) {
					SurfaceDistances origDists = origSurfs[r].getDistances(loc);
					SurfaceDistances newDists = newSurfs[r].getDistances(loc);
					Preconditions.checkState(Precision.equals(origSurfs[r].getQuickDistance(loc), newSurfs[r].getQuickDistance(loc), tol),
							"Rquick mismatch for %s: %s != %s", r, origSurfs[r].getQuickDistance(loc), newSurfs[r].getQuickDistance(loc));
					Preconditions.checkState(Precision.equals(origDists.getDistanceRup(), newDists.getDistanceRup(), tol),
							"Rrup mismatch for %s: %s != %s", r, origDists.getDistanceRup(), newDists.getDistanceRup());
					Preconditions.checkState(Precision.equals(origDists.getDistanceJB(), newDists.getDistanceJB(), tol),
							"Rjb mismatch for %s: %s != %s", r, origDists.getDistanceJB(), newDists.getDistanceJB());
					Preconditions.checkState(Precision.equals(origDists.getDistanceX(), newDists.getDistanceX(), tol),
							"Rx mismatch for %s: %s != %s", r, origDists.getDistanceX(), newDists.getDistanceX());
				}
			}));
		}
		for (int i=0; i<futures.size(); i++) {
			CompletableFuture<Void> future = futures.get(i);
			future.join();
			System.out.print(".");
			if (i > 0 && i % 100 == 0)
				System.out.println(" "+i+"/"+futures.size());
		}
		System.out.println(" DONE "+futures.size());
	}
	
	private static void debugSurf(File outputDir, FaultSystemRupSet rupSet, int rupIndex,
			CompoundSurface origSurf, CompoundSurface newSurf) throws IOException {
		List<FaultSection> sects = rupSet.getFaultSectionDataForRupture(rupIndex);
		GeographicMapMaker mapMaker = new GeographicMapMaker(sects);
		mapMaker.setWriteGeoJSON(false);
		mapMaker.setWritePDFs(false);
		
		

//		System.out.println("\tOrig dip: "+origSurf.getAveDip());
//		System.out.println("\tNew dip: "+newSurf.getAveDip());
//		System.out.println("\tOrig strike: "+origSurf.getAveStrike());
//		System.out.println("\tNew strike: "+newSurf.getAveStrike());
		
		LocationList scatters = new LocationList();
		List<PlotCurveCharacterstics> scatterChars = new ArrayList<>();
		
		scatters.add(origSurf.getFirstLocOnUpperEdge());
		scatterChars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 12f, Colors.tab_green));
		scatters.add(origSurf.getLastLocOnUpperEdge());
		scatterChars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 12f, Colors.tab_red));
		scatters.add(origSurf.getFirstLocOnLowerEdge());
		scatterChars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 10f, Colors.tab_green));
		scatters.add(origSurf.getLastLocOnLowerEdge());
		scatterChars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 10f, Colors.tab_red));
		scatters.add(newSurf.getFirstLocOnUpperEdge());
		scatterChars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_SQUARE, 7f, Colors.tab_green.darker().darker()));
		scatters.add(newSurf.getLastLocOnUpperEdge());
		scatterChars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_SQUARE, 7f, Colors.tab_red.darker().darker()));
		scatters.add(newSurf.getFirstLocOnLowerEdge());
		scatterChars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_SQUARE, 5f, Colors.tab_green.darker().darker()));
		scatters.add(newSurf.getLastLocOnLowerEdge());
		scatterChars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_SQUARE, 5f, Colors.tab_red.darker().darker()));
		
		LocationList origUpper = origSurf.getUpperEdge();
		LocationList newUpper = newSurf.getUpperEdge();
		
//		double strike = newSurf.getAveStrike();
//		CPT orderCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance();
//		orderCPT = orderCPT.rescale(0d, origUpper.size()-1d);
//		LocationVector offset = new LocationVector(strike+90d, 2d, 0d);
//		for (int i=0; i<origUpper.size(); i++) {
//			scatters.add(LocationUtils.location(origUpper.get(i), offset));
//			scatterChars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 1.5f, orderCPT.getColor((double)i)));
//		}
//		orderCPT = orderCPT.rescale(0d, newUpper.size()-1d);
//		offset = new LocationVector(strike-90d, 2d, 0d);
//		for (int i=0; i<newUpper.size(); i++) {
//			scatters.add(LocationUtils.location(newUpper.get(i), offset));
//			scatterChars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_SQUARE, 1.5f, orderCPT.getColor((double)i)));
//		}
		
		List<LocationList> lines = new ArrayList<>();
		List<PlotCurveCharacterstics> lineChars = new ArrayList<>();
		lines.add(origUpper);
		lineChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Colors.tab_orange));
		lines.add(newUpper);
		lineChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Colors.tab_blue));
		lines.add(origUpper);
		lineChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Colors.tab_orange));
		mapMaker.plotLines(lines, lineChars);
		
		List<LocationList> arrows = new ArrayList<>();
		for (int s=1; s<sects.size(); s++) {
			FaultSection prev = sects.get(s-1);
			FaultSection sect = sects.get(s);
			if (sect.getParentSectionId() != prev.getParentSectionId()) {
				// draw jump
				arrows.add(LocationList.of(middle(origSurf.getSurfaceList().get(s-1)), middle(origSurf.getSurfaceList().get(s))));
			}
		}
		mapMaker.plotArrows(arrows, 2d, Colors.tab_brown, 1f);
		
		mapMaker.plotScatters(scatters, scatterChars);
		
		mapMaker.plot(outputDir, "rupture_"+rupIndex, " ");
	}
	
	private static void surfBuildTime(FaultSystemRupSet rupSet, int trials) {
		Stopwatch watch = Stopwatch.createUnstarted();
		for (int i=0; i<trials; i++) {
			System.out.println("Surf build trial "+i);
			rupSet.clearCache();
			watch.start();
			for (int r=0; r<rupSet.getNumRuptures(); r++) {
				rupSet.getSurfaceForRupture(r, 1d);
			}
			watch.stop();
		}
		double timeSecs = watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
		System.out.println("Built "+rupSet.getNumRuptures()+" rutures "+trials+" time(s) in "+(float)timeSecs+" s");
	}
	
	private static Location middle(RuptureSurface surf) {
		LocationAverager avg = new LocationAverager();
		for (Location loc : surf.getEvenlyDiscritizedListOfLocsOnSurface())
			avg.add(loc, 1d);
		return avg.getAverage();
	}

}
