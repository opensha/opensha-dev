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
import org.opensha.commons.util.Interpolate;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.NamedFaults;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.earthquake.param.HistoricOpenIntervalParam;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;

import com.google.common.base.Preconditions;

import scratch.UCERF3.analysis.FaultSysSolutionERF_Calc;
import scratch.UCERF3.erf.FaultSystemSolutionERF;

public class SCECStatewide {

	public static void main(String[] args) throws IOException {
		Region relm = new CaliforniaRegions.RELM_TESTING();
		Region ca = new Region(new Location(relm.getMinLat(), relm.getMinLon()),
				new Location(relm.getMaxLat(), relm.getMaxLon()));
		XY_DataSet[] caOutlines = PoliticalBoundariesData.loadCAOutlines();
		GeographicMapMaker mapMaker = new GeographicMapMaker(ca, caOutlines);
		
		boolean faultProbs = true;
		
//		File outputDir = new File("/home/kevin/SCEC/2023_annual_meeting_program");
		File outputDir = new File("/home/kevin/SCEC/2023_statewide_announcement");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		NSHM23_FaultModels fm = NSHM23_FaultModels.WUS_FM_v2;
		boolean subSects = faultProbs;
		List<? extends FaultSection> sects;
		FaultSystemSolution sol = null;
		if (subSects) {
			if (faultProbs) {
				sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
						+ "2023_06_23-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
						+ "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip"));
				sects = sol.getRupSet().getFaultSectionDataList();
			} else {
				sects = fm.getDefaultDeformationModel().build(fm);
			}
		} else {
			sects = fm.getFaultSections();
		}
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
		
		double minSaturateDist = faultProbs ? 150 : 25d;
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
		
		Runnable faultCall;
		
		if (faultProbs) {
			FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
			erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
			erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
			erf.getTimeSpan().setDuration(30d);
//			erf.setParameter(HistoricOpenIntervalParam.NAME, 2022d-1875d);
			erf.updateForecast();
			
			double minMag = 6.5;
			
			boolean logProbs = true;
//			CPT probCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-2d, 0d);
			CPT probCPT = new CPT(Math.log10(0.01), Math.log10(0.5),
					new Color(254,224,144),
					new Color(253,174,97),
					new Color(244,109,67),
					new Color(215,48,39),
					new Color(165,0,38),
					new Color(82,0,19));
			double minThickness = 1.5;
			double maxThickness = 8;
			probCPT.setBelowMinColor(probCPT.getMinColor());
			probCPT.setAboveMaxColor(probCPT.getMaxColor());
			probCPT.setNanColor(new Color(255, 255, 255, 0));
			
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			List<Double> comps = new ArrayList<>();
			for (int s=0; s<sects.size(); s++) {
				FaultSection sect = sects.get(s);
				double prob = FaultSysSolutionERF_Calc.calcParticipationProbForSect(
						erf, minMag, sect.getSectionId());
				comps.add(prob);
				if (logProbs)
					prob = Math.log10(prob);
				if (!Double.isFinite(prob) || prob == 0d) {
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.WHITE));
				} else {
					double thickness;
					if (highlightSects.contains(sect))
						thickness = 5d;
					else if (prob <= probCPT.getMinValue())
//					if (prob <= probCPT.getMinValue())
						thickness = minThickness;
					else if (prob >= probCPT.getMaxValue())
						thickness = maxThickness;
					else
						thickness = Interpolate.findY(probCPT.getMinValue(), minThickness,
								probCPT.getMaxValue(), maxThickness, prob);
					thickness = Math.min(thickness, 5d);
					Color color = probCPT.getColor((float)prob);
					if (!highlightSects.contains(sect)) {
						double dist = scalars[s];
						if (Double.isNaN(dist))
							dist = Double.POSITIVE_INFINITY;
						if (dist > minSaturateDist) {
							double weight;
							if (dist >= maxPlotDist) {
								weight = 1d;
								color = Color.WHITE;
							} else {
								weight = (dist-minSaturateDist)/(maxPlotDist-minSaturateDist);
								color = blend(Color.WHITE, color, weight);
							}
//							Preconditions.checkState(weight >= 0d && weight <= 1d, "Bad weight: %s", weight);
//							color = new Color(color.getRed(), color.getGreen(), color.getBlue(), (int)(255d*weight));
						}
//						int alpha = cpt.getColor((float)distScalar).getAlpha();
////						color = new Color(color.getRed(), color.getGreen(), color.getBlue(), alpha);
//						double weight = (double)alpha
//						color = blend(scecRed, color, alpha)
					}
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, (float)thickness, color));
				}
			}
			faultCall = new Runnable() {
				
				@Override
				public void run() {
					mapMaker.setSectOutlineChar(null);
					mapMaker.plotSectChars(chars, probCPT, "30 Year Probability", comps);
				}
			};
			
		} else {
			faultCall = new Runnable() {
				
				@Override
				public void run() {
					mapMaker.setSectOutlineChar(null);
					mapMaker.setSectHighlights(highlightSects, new PlotCurveCharacterstics(PlotLineType.SOLID, 7f, scecRed));
					mapMaker.plotSectScalars(scalars, cpt, null);
					mapMaker.setScalarThickness(3f);
				}
			};
		}
		faultCall.run();
		
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
			
			mapMaker.clearSectChars();
			mapMaker.clearSectScalars();
			mapMaker.clearFaultSections();
			writePlot(outputDir, "eqs_only", gp, mapMaker);
		}
		
		mapMaker.setPoliticalBoundaries(caOutlines);
		mapMaker.setFaultSections(sects);
		faultCall.run();
		
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
	
	private static Color blend(Color c1, Color c2, double weight) {
		float r = (float)((weight*c1.getRed() + (1d-weight)*c2.getRed())/255d);
		float g = (float)((weight*c1.getGreen() + (1d-weight)*c2.getGreen())/255d);
		float b = (float)((weight*c1.getBlue() + (1d-weight)*c2.getBlue())/255d);
		return new Color(r, g, b);
	}

}
