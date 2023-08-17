package scratch.kevin.miscFigures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.function.Supplier;

import org.opensha.commons.data.comcat.ComcatAccessor;
import org.opensha.commons.data.comcat.ComcatRegionAdapter;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.PoliticalBoundariesData;
import org.opensha.commons.util.FaultUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.modules.NamedFaults;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;

public class SCECStatewide {

	public static void main(String[] args) throws IOException {
		Region relm = new CaliforniaRegions.RELM_TESTING();
		Region ca = new Region(new Location(relm.getMinLat(), relm.getMinLon()),
				new Location(relm.getMaxLat(), relm.getMaxLon()));
		XY_DataSet[] caOutlines = PoliticalBoundariesData.loadCAOutlines();
		GeographicMapMaker mapMaker = new GeographicMapMaker(ca, caOutlines);
		
		File outputDir = new File("/home/kevin/SCEC/2023_annual_meeting_program"); 
		
		NSHM23_FaultModels fm = NSHM23_FaultModels.NSHM23_v2;
		boolean subSects = false;
		List<? extends FaultSection> sects;
		if (subSects)
			sects = fm.getDefaultDeformationModel().build(fm);
		else
			sects = fm.getFaultSections();
		NamedFaults namedFaults = fm.getNamedFaults();
		
		boolean includeEQs = true;
		boolean distFadeEQs = false;
		
		HashSet<Integer> namedSects = new HashSet<>();
		namedSects.addAll(namedFaults.getParentIDsForFault("San Andreas"));
		namedSects.addAll(namedFaults.getParentIDsForFault("San Jacinto"));
		namedSects.addAll(namedFaults.getParentIDsForFault("Hayward-Rodgers Creek"));
		for (FaultSection sect : sects) {
			if (sect.getName().contains("Mendocino") || sect.getName().contains("Calaveras")) {
				if (subSects)
					namedSects.add(sect.getParentSectionId());
				else
					namedSects.add(sect.getSectionId());
			}
		}
		
		Color scecRed = new Color(153, 0, 0);
		
		HashSet<FaultSection> highlightSects = new HashSet<>();
		List<Location> highlightLocs = new ArrayList<>();
		for (FaultSection sect : sects) {
			int id = subSects ? sect.getParentSectionId() : sect.getSectionId();
			if (namedSects.contains(id)) {
				highlightSects.add(sect);
				FaultTrace trace = sect.getFaultTrace();
				if (!subSects && trace.getTraceLength() > 10d)
					trace = FaultUtils.resampleTrace(trace, (int)(trace.getTraceLength()/5d + 0.5));
				highlightLocs.addAll(trace);
			}
		}
		
		double minSaturateDist = 25d;
		double maxPlotDist = 400d;
		CPT cpt = new CPT(minSaturateDist, maxPlotDist, new Color(127, 127, 127, 120), new Color(127, 127, 127, 0));
		cpt.setBelowMinColor(cpt.getMinColor());
		cpt.setAboveMaxColor(cpt.getMaxColor());
		CPT eqCPT = new CPT(minSaturateDist, maxPlotDist, new Color(127, 127, 127, 255), new Color(127, 127, 127, 0));
		eqCPT.setBelowMinColor(eqCPT.getMinColor());
		eqCPT.setAboveMaxColor(eqCPT.getMaxColor());
		
		System.out.println("Have "+highlightLocs.size()+" highlight trace distance locs");
		
		List<CompletableFuture<Double>> distFutures = new ArrayList<>(sects.size());
		for (FaultSection sect : sects) {
			if (!highlightSects.contains(sect)) {
				distFutures.add(CompletableFuture.supplyAsync(new Supplier<Double>() {

					@Override
					public Double get() {
						double minDist = Double.POSITIVE_INFINITY;
						for (Location loc1 : sect.getFaultTrace())
							for (Location loc2 : highlightLocs)
								minDist = Math.min(minDist, LocationUtils.horzDistanceFast(loc1, loc2));
						if (minDist > maxPlotDist)
							return Double.NaN;
						return minDist;
					}
				}));
				
			} else {
				distFutures.add(null);
			}
		}
		System.out.println("Waiting on "+distFutures.size()+" dist futures");
		double[] scalars = new double[sects.size()];
		for (int s=0; s<sects.size(); s++) {
			CompletableFuture<Double> future = distFutures.get(s);
			if (future != null)
				scalars[s] = future.join();
		}
		System.out.println("Done, now plotting");
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		mapMaker.setPoliticalBoundaryChar(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		writePlot(outputDir, "ca_only", gp, mapMaker);
		
		mapMaker.setFaultSections(sects);
		
		mapMaker.setSectOutlineChar(null);
		mapMaker.setSectHighlights(highlightSects, new PlotCurveCharacterstics(PlotLineType.SOLID, 7f, scecRed));
		mapMaker.plotSectScalars(scalars, cpt, null);
		mapMaker.setScalarThickness(3f);
		mapMaker.setSkipNaNs(true);
		
		mapMaker.setPoliticalBoundaries(null);
		writePlot(outputDir, "faults_only", gp, mapMaker);
		
		if (includeEQs) {
			ComcatAccessor comcat = new ComcatAccessor();
//			long startTime = -3786778739000l; // 1850
			long startTime = 665413261000l; // 2/2/91 -> SCEC founded
			long endTime = System.currentTimeMillis();
			double minMag = 5d;
			ObsEqkRupList events = comcat.fetchEventList(null, startTime, endTime, -10, 50d,
					new ComcatRegionAdapter(ca), false, false, minMag, 1000, 100);
			events.sortByMag();
			System.out.println("Fetched "+events.size()+" events");
			List<Location> eqLocs = new ArrayList<>();
			List<PlotCurveCharacterstics> eqChars = new ArrayList<>();
			for (ObsEqkRupture event : events) {
				Location hypo = event.getHypocenterLocation();
//				float thickness = (float)(2d + 6*(Math.min(8d, event.getMag()) - minMag));
				float thickness = (float)(3d + 6*(Math.min(8d, event.getMag()) - minMag));
				Color fillColor, outlineColor;
				if (distFadeEQs) {
					double minDist = Double.POSITIVE_INFINITY;
					for (Location loc : highlightLocs)
						minDist = Math.min(minDist, LocationUtils.horzDistance(loc, hypo));
					if (minDist > maxPlotDist)
						continue;
					fillColor = eqCPT.getColor((float)minDist);
					outlineColor = new Color(0, 0, 0, fillColor.getAlpha());
				} else {
					fillColor = Color.GRAY;
					outlineColor = Color.BLACK;
				}
				eqLocs.add(hypo);
				eqChars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, thickness, fillColor));
//				// add outline on top
//				eqLocs.add(hypo);
//				eqChars.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE, thickness, outlineColor));
			}
			mapMaker.plotScatters(eqLocs, eqChars, null);
			
			mapMaker.clearSectScalars();
			mapMaker.clearFaultSections();
			writePlot(outputDir, "eqs_only", gp, mapMaker);
		}
		
		mapMaker.setPoliticalBoundaries(caOutlines);
		mapMaker.setFaultSections(sects);
		mapMaker.setSectHighlights(highlightSects, new PlotCurveCharacterstics(PlotLineType.SOLID, 7f, scecRed));
		mapMaker.plotSectScalars(scalars, cpt, null);
		
		writePlot(outputDir, "full", gp, mapMaker);
	}
	
	private static void writePlot(File outputDir, String prefix, HeadlessGraphPanel gp, GeographicMapMaker mapMaker) throws IOException {
		System.out.println("Writing plot: "+prefix);
		
		PlotSpec spec = mapMaker.buildPlot(" ");
		
		gp.drawGraphPanel(spec, false, false, mapMaker.getXRange(), mapMaker.getYRange());
		
		PlotUtils.setAxisVisible(gp, false, false);
		PlotUtils.setXTick(gp, 1000);
		PlotUtils.setYTick(gp, 1000);
		PlotUtils.writePlots(outputDir, prefix, gp, 1200, true, true, true, false);
	}

}
