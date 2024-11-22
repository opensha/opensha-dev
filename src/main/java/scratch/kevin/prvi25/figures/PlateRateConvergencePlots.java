package scratch.kevin.prvi25.figures;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.geo.json.FeatureCollection;
import org.opensha.commons.geo.json.Geometry.DepthSerializationType;
import org.opensha.commons.geo.json.Geometry.MultiLineString;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalDeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalFaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionDeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionFaultModels;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;

import com.google.common.base.Preconditions;
import com.google.gson.Gson;

import net.mahdilamb.colormap.Colors;
import scratch.kevin.prvi25.FaultSystemLineIntegralCalculator;
import scratch.kevin.prvi25.FaultSystemLineIntegralCalculator.LineIntegralResult;
import scratch.kevin.prvi25.FaultSystemLineIntegralCalculator.VectorComponent;

public class PlateRateConvergencePlots {

	public static void main(String[] args) throws IOException {
		PRVI25_CrustalFaultModels crustalFM = PRVI25_CrustalFaultModels.PRVI_CRUSTAL_FM_V1p1;
		PRVI25_CrustalDeformationModels crustalDM = PRVI25_CrustalDeformationModels.GEOLOGIC_DIST_AVG;
		
		List<? extends FaultSection> crustalSects = crustalDM.build(crustalFM);
		
		String plateConvergencePath = "/data/erf/prvi25/def_models/plate_convergence/2024_11_20.csv";
		InputStream convergeIS = PRVI25_CrustalDeformationModels.class.getResourceAsStream(plateConvergencePath);
		Preconditions.checkNotNull(convergeIS);
		CSVFile<String> convergeCSV = CSVFile.readStream(convergeIS, true);
		
		String refDMGeoJsonPath = "/data/erf/prvi25/fault_models/subduction/inputs/PRVI_subduction_FullRate_LargePolys_drape_removeFields.geojson";
		InputStream subDmIS = PRVI25_CrustalDeformationModels.class.getResourceAsStream(refDMGeoJsonPath);
		Preconditions.checkNotNull(subDmIS);
		Gson inGSON = FeatureCollection.buildGson(DepthSerializationType.DEPTH_M);
		BufferedReader read = new BufferedReader(new InputStreamReader(subDmIS));
		List<Feature> inFeatures = inGSON.fromJson(read, FeatureCollection.class).features;
		read.close();
		
		File outputDir = new File("/tmp");
		
		Map<Integer, FaultTrace> subRefTraces = new HashMap<>();
		for (Feature feature : inFeatures) {
			int id = feature.properties.getInt("id", -1);
			Preconditions.checkState(id > 0);
			System.out.println("Loading trace for "+id+". "+feature.properties.getString("FaultDesc"));
			MultiLineString mls = (MultiLineString)feature.geometry;
			System.out.println("\tHave "+mls.lines.size()+" line strings");
			FaultTrace upper = new FaultTrace(null);
			upper.addAll(mls.lines.get(0));
			System.out.println("\tUpper strike: "+upper.getAveStrike());
			subRefTraces.put(id, upper);
		}
		
		boolean includeFullInAll = false;
		VectorComponent[] components = VectorComponent.values();
		
//		boolean includeFullInAll = true;
//		VectorComponent[] components = {
//				VectorComponent.PARALLEL,
//				VectorComponent.PERPENDICULAR
//		};
		
		double[] coupledExtraFracts = { 0.2, 0.4, 0.6, 0.8 };
		
		Map<VectorComponent, DiscretizedFunc> plateRateFuncs = new EnumMap<>(VectorComponent.class);
		Map<VectorComponent, DiscretizedFunc> plateRateCoupledFuncs = new EnumMap<>(VectorComponent.class);
		Map<VectorComponent, DiscretizedFunc[]> plateRateCoupledExtraFuncs = new EnumMap<>(VectorComponent.class);
		
		for (VectorComponent component : VectorComponent.values()) {
			plateRateFuncs.put(component, new ArbitrarilyDiscretizedFunc());
			plateRateCoupledFuncs.put(component, new ArbitrarilyDiscretizedFunc());
			DiscretizedFunc[] coupledExtra = new DiscretizedFunc[coupledExtraFracts.length];
			for (int i=0; i<coupledExtraFracts.length; i++)
				coupledExtra[i] = new ArbitrarilyDiscretizedFunc();
			plateRateCoupledExtraFuncs.put(component, coupledExtra);
		}
		
		for (int row=1; row<convergeCSV.getNumRows(); row++) {
			int id = convergeCSV.getInt(row, 0);
			if (id > 7550) {
				System.out.println("Skipping "+id);
				continue;
			}
			System.out.println("Building convergence vector for "+id);
			FaultTrace trace = subRefTraces.get(id);
			Preconditions.checkNotNull(trace);
			
			double east = convergeCSV.getDouble(row, 16);
			double north = convergeCSV.getDouble(row, 17);
			double full = Math.sqrt(east*east + north*north);
			double coupling = convergeCSV.getDouble(row, 18);
			System.out.println("\tVector: north="+(float)north+", east="+(float)east+", coupling="+(float)coupling);
			
			double len = trace.getTraceLength();
			double cmlLen = 0d;
			double halfLen = 0.5*len;
			Location prev = trace.first();
			Location middleLoc = null;
			for (int i=1; i<trace.size(); i++) {
				Location loc = trace.get(i);
				double newCmlLen = cmlLen + LocationUtils.horzDistance(prev, loc);
				if (newCmlLen >= halfLen) {
					double delta = halfLen - cmlLen;
					Preconditions.checkState(delta > 0,
							"bad delta=%s for i=%s/%s, halfLen=%s, prevCmlLen=%s, newCmlLen=%s",
							(float)delta, i, trace.size(), (float)halfLen, (float)cmlLen, (float)newCmlLen);
					middleLoc = LocationUtils.location(prev, LocationUtils.azimuthRad(prev, loc), delta);
					break;
				}
				cmlLen = newCmlLen;
				prev = loc;
			}
			Preconditions.checkNotNull(middleLoc);
			System.out.println("\tMiddle loc: "+middleLoc);
			
			for (VectorComponent comp : VectorComponent.values()) {
				double plateRate;
				switch (comp) {
				case FULL_HORIZONTAL:
					plateRate = full;
					break;
				case PARALLEL:
					plateRate = north;
					break;
				case PERPENDICULAR:
					plateRate = east;
					break;

				default:
					throw new IllegalStateException();
				}
				plateRateFuncs.get(comp).set(middleLoc.lon, plateRate);
				plateRateCoupledFuncs.get(comp).set(middleLoc.lon, plateRate*coupling);
				for (int i=0; i<coupledExtraFracts.length; i++)
					plateRateCoupledExtraFuncs.get(comp)[i].set(middleLoc.lon, plateRate*coupledExtraFracts[i]);
			}
		}

		double minLat = 16;
		double maxLat = 21;
		double minLon = -71d;
		double maxLon = -60.5d;
		double delta = 0.1d;
		List<Location[]> integralLocs = new ArrayList<>();
		for (double lon=minLon; (float)lon<=(float)maxLon; lon += delta)
			integralLocs.add(new Location[] {new Location(minLat, lon), new Location(maxLat, lon)});
		List<Location[]> mapLineLocs = new ArrayList<>();
		for (double lon=minLon+1; (float)lon<=(float)maxLon; lon += 1)
			mapLineLocs.add(new Location[] {new Location(minLat, lon), new Location(maxLat, lon)});
		
		FaultSystemLineIntegralCalculator crustalCalc = new FaultSystemLineIntegralCalculator(crustalSects, true);
		
		double mapDelta = 1d;
		List<LineIntegralResult> mapFakeIntegrals = new ArrayList<>();
		for (double lon=minLon; (float)lon<=(float)maxLon; lon += mapDelta)
			mapFakeIntegrals.add(crustalCalc.calcLineIntegral(new Location(minLat, lon), new Location(maxLat, lon)));
		
		List<LineIntegralResult> crustalIntegrals = new ArrayList<>(integralLocs.size());
		
		for (Location[] loc : integralLocs)
			crustalIntegrals.add(crustalCalc.calcLineIntegral(loc[0], loc[1]));
		
		Map<VectorComponent, DiscretizedFunc> crustalDataFuncs = new EnumMap<>(VectorComponent.class);
		for (VectorComponent comp : components)
			crustalDataFuncs.put(comp, crustalCalc.buildIntegralFunction(false, crustalIntegrals, comp));
		
		for (PRVI25_SubductionFaultModels subFM : PRVI25_SubductionFaultModels.values()) {
			if (subFM.getNodeWeight(null) == 0d)
				continue;
			for (PRVI25_SubductionDeformationModels subDM : PRVI25_SubductionDeformationModels.values()) {
				if (subDM.getNodeWeight(null) == 0d)
					continue;
				System.out.println("Plotting "+subFM.getShortName()+", "+subDM.getShortName());
				List<? extends FaultSection> subductionSects = subDM.build(subFM);
				
				List<FaultSection> combinedSects = FaultSystemLineIntegralCalculator.unionSubSectLists(
						List.of(crustalSects, subductionSects));
				
				FaultSystemLineIntegralCalculator calc = new FaultSystemLineIntegralCalculator(combinedSects, true);
				
				List<LineIntegralResult> integrals = new ArrayList<>(integralLocs.size());
				
				for (Location[] loc : integralLocs)
					integrals.add(calc.calcLineIntegral(loc[0], loc[1]));
				Map<VectorComponent, DiscretizedFunc> combDataFuncs = new EnumMap<>(VectorComponent.class);
				for (VectorComponent comp : components)
					combDataFuncs.put(comp, calc.buildIntegralFunction(false, integrals, comp));
				
				FaultSystemLineIntegralCalculator subductionCalc = new FaultSystemLineIntegralCalculator(subductionSects, true);
				
				List<LineIntegralResult> subductionIntegrals = new ArrayList<>(integralLocs.size());
				
				for (Location[] loc : integralLocs)
					subductionIntegrals.add(subductionCalc.calcLineIntegral(loc[0], loc[1]));
				
				Map<VectorComponent, DiscretizedFunc> subductionDataFuncs = new EnumMap<>(VectorComponent.class);
				for (VectorComponent comp : components)
					subductionDataFuncs.put(comp, subductionCalc.buildIntegralFunction(false, subductionIntegrals, comp));
				
				List<PlotSpec> plots = new ArrayList<>();
				List<Range> xRanges = List.of(new Range(minLon, maxLon));
				List<Range> yRanges = new ArrayList<>();
				
				for (int c=0; c<components.length; c++) {
					VectorComponent component = components[c];
					List<DiscretizedFunc> funcs = new ArrayList<>();
					List<PlotCurveCharacterstics> chars = new ArrayList<>();
					
					DiscretizedFunc refFull = plateRateFuncs.get(component);
					DiscretizedFunc refCoupled = plateRateCoupledFuncs.get(component);
					DiscretizedFunc combData = combDataFuncs.get(component);
					DiscretizedFunc crustalData = crustalDataFuncs.get(component);
					DiscretizedFunc subductionData = subductionDataFuncs.get(component);
					
					PlotCurveCharacterstics extraChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GRAY);
					
					for (DiscretizedFunc extraFunc : plateRateCoupledExtraFuncs.get(component)) {
						funcs.add(extraFunc);
						chars.add(extraChar);
					}
					
					crustalData.setName("Crustal only");
					funcs.add(crustalData);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Colors.tab_green));
					
					subductionData.setName("Subduction ("+subFM.getShortName()+"-"+subDM.getShortName()+") only");
					funcs.add(subductionData);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Colors.tab_blue));

					combData.setName("Crustal + Subduction");
					funcs.add(combData);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Colors.tab_orange));
					
					refFull.setName("Block model plate rate");
					funcs.add(refFull);
					if (!includeFullInAll || component == VectorComponent.FULL_HORIZONTAL)
						chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, PlotSymbol.FILLED_CIRCLE, 6f, Color.BLACK));
					else
						chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
					
					refCoupled.setName("Interface-coupled");
					funcs.add(refCoupled);
					chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 4f, Colors.BLACK));
					
					if (coupledExtraFracts.length > 0) {
						DecimalFormat pDF = new DecimalFormat("0.#%");
						DiscretizedFunc fakeFunc = new EvenlyDiscretizedFunc(-1000d, -100d, 10);
						String name = null;
						for (double coupling : coupledExtraFracts) {
							if (name == null)
								name = "";
							else
								name += ", ";
							name += pDF.format(coupling);
						}
						fakeFunc.setName(name+" Coupled");
						funcs.add(fakeFunc);
						chars.add(extraChar);
					}
					
					if (includeFullInAll && component != VectorComponent.FULL_HORIZONTAL) {
						DiscretizedFunc fullPlate = plateRateFuncs.get(VectorComponent.FULL_HORIZONTAL);
						fullPlate.setName("|Full horizontal plate rate|");
						funcs.add(fullPlate);
						chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Colors.tab_purple));
						DiscretizedFunc fullCoupledPlate = plateRateCoupledFuncs.get(VectorComponent.FULL_HORIZONTAL);
						fullCoupledPlate.setName("Coupled");
						funcs.add(fullCoupledPlate);
						chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 4f, Colors.tab_purple));
					}
					
					Range yRange = null;
					String yAxisLabel;
					switch (component) {
					case FULL_HORIZONTAL:
						yAxisLabel = "|Horizontal|";
						yRange = new Range(0, 20);
						break;
					case PARALLEL:
						yAxisLabel = "N-S";
						yRange = new Range(-10, 10);
						break;
					case PERPENDICULAR:
						yAxisLabel = "E-W";
						yRange = new Range(-20, 0);
						break;

					default:
						throw new IllegalStateException();
					}
					yAxisLabel += " (mm/yr)";
					
					PlotSpec plot = new PlotSpec(funcs, chars,
							"Crustal and Subduction ("+subFM.getShortName()+", "+subDM.getShortName()+")", "Longitude", yAxisLabel);
					plot.setLegendVisible(c == components.length-1);
					
					plots.add(plot);
//					yRanges.add(FaultSystemLineIntegralCalculator.getPlotYRange(plot));
					yRanges.add(yRange);
				}
				
				HeadlessGraphPanel gp = PlotUtils.initHeadless();
				
				gp.drawGraphPanel(plots, false, false, xRanges, yRanges);
				
				PlotUtils.writePlots(outputDir, "prvi_plate_convergence_"+subFM.getFilePrefix()+"_"+subDM.getFilePrefix(),
						gp, 900, 1200, true, true, false);
				
				// now add the map
//				List<LineIntegralResult> mapLines = List.of();
				List<LineIntegralResult> mapLines = new ArrayList<>(mapLineLocs.size());
//				for (Location[] mapLoc : mapLineLocs.subList(0, 1))
				for (Location[] mapLoc : mapLineLocs)
					mapLines.add(calc.calcLineIntegral(mapLoc[0], mapLoc[1]));
				GeographicMapMaker mapMaker = calc.buildMapPlot(false, false, 10d, mapLines);
				PlotSpec mapPlot = mapMaker.buildPlot(plots.get(0).getTitle());
				mapPlot.setSubtitles(null);
				plots.add(0, mapPlot);
				yRanges.add(0, mapMaker.getYRange());
				
				gp.drawGraphPanel(plots, false, false, xRanges, yRanges);
				
//				double aspectRatio = PlotUtils.calcAspectRatio(xRanges.get(0), yRanges.get(0), true);
//				Preconditions.checkState(aspectRatio > 0d);
//				int mapWeight = 20;
//				int chartWeight = (int)(mapWeight / aspectRatio + 0.5);
//				int[] weights = new int[plots.size()];
//				weights[0] = mapWeight;
//				for (int i=1; i<weights.length; i++)
//					weights[i] = chartWeight;
//				PlotUtils.setSubPlotWeights(gp, weights);
				
				int mapWeight = 20;
				int chartWeight = 13;
				int[] weights = new int[plots.size()];
				weights[0] = mapWeight;
				for (int i=1; i<weights.length; i++)
					weights[i] = chartWeight;
				PlotUtils.setSubPlotWeights(gp, weights);
				
				PlotUtils.writePlots(outputDir, "prvi_plate_convergence_map_"+subFM.getFilePrefix()+"_"+subDM.getFilePrefix(),
						gp, 900, true, true, true, false);
			}
		}
	}

}
