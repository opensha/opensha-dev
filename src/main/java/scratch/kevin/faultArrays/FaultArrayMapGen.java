package scratch.kevin.faultArrays;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.siteData.impl.SRTM30PlusTopography;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.elements.CoastAttributes;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.mapping.gmt.elements.PSText;
import org.opensha.commons.mapping.gmt.elements.PSText.Justify;
import org.opensha.commons.mapping.gmt.elements.PSXYPolygon;
import org.opensha.commons.mapping.gmt.elements.PSXYSymbol;
import org.opensha.commons.mapping.gmt.elements.TopographicSlopeFile;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Ints;

import scratch.UCERF3.U3FaultSystemSolution;
import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.UCERF3.analysis.FaultSysSolutionERF_Calc;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.erf.mean.MeanUCERF3;
import scratch.UCERF3.erf.mean.MeanUCERF3.Presets;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;
import scratch.kevin.miscFigures.U3vsPopulationMap;

public class FaultArrayMapGen {
	
	private static enum ScalarType {
		SUBSECT_PARTIC_PROB("Subsection M>="+(float)magForProbs+" "+(int)durationForProbs+"yr Participation Prob"),
		PARENT_PARTIC_PROB("Section M>="+(float)magForProbs+" "+(int)durationForProbs+"yr Participation Prob"),
		TIME_SINCE_LAST("Time Since Last Event (years)"),
		SUBSECT_NORM_RI("Subsection Normalized Recurrence Interval"),
		SOLID("Solid"),
		ALL_GRAY("All Gray");
		
		private String label;

		private ScalarType(String label) {
			this.label = label;
		}
	}
	
//	private static Region region = new Region(new Location(35.1, -114.5), new Location(32, -120));
//	private static Region region = new Region(new Location(35.5, -114), new Location(31.5, -121));
//	private static Region region = new Region(new Location(36.5, -114), new Location(31.5, -121));
//	private static Region region = new Region(new Location(36, -114), new Location(31.5, -120.5)); // max upper left of Yehuda's GRD file
	private static Region region = new Region(new Location(31.75, -125), new Location(42.5, -114));
	
	private static boolean srtm = true;
	private static GriddedGeoDataSet topoXYZ = null;
	
	private static synchronized GriddedGeoDataSet fetchTopo(TopographicSlopeFile topoRes) throws IOException {
		if (topoXYZ != null)
			return topoXYZ;
		SRTM30PlusTopography topo = new SRTM30PlusTopography();
		double discr = Math.max((double)topoRes.resolution()/3600d, topo.getResolution());
		System.out.println("Topo discretization: "+discr);
		GriddedRegion gridReg = new GriddedRegion(region, discr, null);
		System.out.println("Fetching "+gridReg.getNodeCount()+" topo values");
		ArrayList<Double> vals = new SRTM30PlusTopography().getValues(gridReg.getNodeList());
		topoXYZ = new GriddedGeoDataSet(gridReg, false);
		for (int i=0; i<topoXYZ.size(); i++)
			topoXYZ.set(i, vals.get(i));
		return topoXYZ;
	}
	
	private static CPT topoCPT = null;
	
	private static synchronized CPT getTopoCPT() throws FileNotFoundException, IOException {
		if (topoCPT == null) {
			topoCPT = CPT.loadFromFile(new File("/home/kevin/SCEC/2019_fault_array_proposal/origGMT/test.cpt"));
		}
		return topoCPT;
	}
	
	private static final double magForProbs = 6.5;
	private static final double durationForProbs = 30d;
	
	private static GriddedGeoDataSet cachedPopData = null;
	
	private static void plotMap(File outputDir, String prefix, double arraySpacingSolid, double arraySpacingDashed,
			FaultModels fm, MeanUCERF3 meanU3, ScalarType scalarType, boolean plotLinearTrace,
			boolean topo, boolean population) throws IOException, GMT_MapException {
		CPT topoCPT = getTopoCPT();
		TopographicSlopeFile topoRes = topo ? TopographicSlopeFile.SRTM_30_PLUS : null;
//		TopographicSlopeFile topoRes = TopographicSlopeFile.US_SIX;
//		GriddedGeoDataSet gridData = fetchTopo(topoRes);
		GriddedGeoDataSet gridData = null;
		String label = null;
		CPT mapCPT = topoCPT;
		double mapSpacing = 3d/3600d;
		if (population) {
			if (cachedPopData == null) {
//				mapSpacing = U3vsPopulationMap.cellsize*5d;
				mapSpacing = U3vsPopulationMap.cellsize;
				gridData = U3vsPopulationMap.fetchPopData(region, mapSpacing);
				cachedPopData = gridData.copy();
			} else {
				gridData = cachedPopData.copy();
			}
			mapCPT = U3vsPopulationMap.getLogPopCPT();
			
			gridData.log10();
			for (int i=0; i<gridData.size(); i++)
				if (!Double.isFinite(gridData.get(i)))
					gridData.set(i, -1d);
			
			label = "Log10 Population Per Cell";
		} else if (srtm && topo) {
			gridData = fetchTopo(topoRes);
		}
//		TopographicSlopeFile topoRes = null;
		GMT_Map map = new GMT_Map(region, gridData, mapSpacing, mapCPT);
		
		map.setBlackBackground(false);
		map.setRescaleCPT(false);
		map.setCustomScaleMin((double)mapCPT.getMinValue());
		map.setCustomScaleMax((double)mapCPT.getMaxValue());
//		map.setCoast(new CoastAttributes(Color.BLACK, 0.6d, new Color(160, 200, 200)));
		map.setCoast(new CoastAttributes(Color.BLACK, 0.6d, new Color(112, 131, 147)));
		map.setCustomLabel(null);
		map.setUseGMTSmoothing(true);
		if (topo) {
			if (gridData != null) {
				map.setTopoResolution(topoRes);
			} else {
				map.setCustomIntenPath("/home/kevin/SCEC/2019_fault_array_proposal/origGMT/socal_mex.grad");
				map.setCustomGRDPath("/home/kevin/SCEC/2019_fault_array_proposal/origGMT/socal_mex.grd");
			}
		} else {
			map.setTopoResolution(null);
		}
//		map.setLabelSize(MAP_LABEL_SIZE);
//		map.setLabelTickSize(MAP_LABEL_TICK_SIZE);
		if (label != null)
			map.setCustomLabel(label);
		map.setDpi(300);
//		map.setImageWidth(imageWidth);
		map.setDrawScaleKM(false);
		map.setGMT_Param("MAP_FRAME_TYPE", "FANCY");
		map.setGMT_Param("MAP_FRAME_WIDTH", "0.04i");
		map.setGMT_Param("MAP_TICK_LENGTH_PRIMARY", "0.04i");
		
		for (FaultSection sect : meanU3.getSolution().getRupSet().getFaultSectionDataList())
			for (PSXYPolygon poly : FaultBasedMapGen.getPolygons(sect.getFaultTrace(), Color.BLACK, 0.3))
				map.addPolys(poly);
		
		Map<Integer, FaultSection> allParents = fm.getFaultSectionIDMap();
		
		CPT cpt = null;
		
		if (scalarType != null) {
			List<FaultSection> sects = new ArrayList<>();
			HashSet<Integer> parentsToPlot = new HashSet<>();
			parentsToPlot.addAll(Ints.asList(FaultArrayCalc.S_SAF_PARENTS));
			parentsToPlot.addAll(Ints.asList(FaultArrayCalc.SJC_PARENTS));
//			parentsToPlot.addAll(Ints.asList(FaultArrayCalc.ELSINORE_PARENTS));
			if (region.getMaxLat() > 38) {
				parentsToPlot.addAll(Ints.asList(FaultArrayCalc.N_SAF_PARENTS));
				parentsToPlot.addAll(Ints.asList(FaultArrayCalc.CALAVERAS_HAYWARD_PARENTS));
			}
			for (FaultSection sect : meanU3.getSolution().getRupSet().getFaultSectionDataList())
				if (parentsToPlot.contains(sect.getParentSectionId()))
					sects.add(sect);
			
//			CPT cpt = GMT_CPT_Files.GMT_DRYWET.instance().reverse();
			cpt = GMT_CPT_Files.MAX_SPECTRUM.instance();
			cpt.setBelowMinColor(cpt.getMinColor());
			cpt.setAboveMaxColor(cpt.getMaxColor());
			Map<Integer, Double> parentProbs = new HashMap<>();
			Map<Integer, Double> parentColorVals = new HashMap<>();
			double[] sectPartRates = null;
			switch (scalarType) {
			case SUBSECT_PARTIC_PROB:
				cpt = cpt.rescale(-4, 0);
				break;
			case PARENT_PARTIC_PROB:
				cpt = cpt.rescale(-4, 0);
				for (Integer parentID : parentsToPlot)
					parentProbs.put(parentID, FaultSysSolutionERF_Calc.calcParticipationProbForParentSects(meanU3, magForProbs, parentID));
				break;
			case TIME_SINCE_LAST:
				cpt = cpt.rescale(0, 2019-1600);
				break;
			case SUBSECT_NORM_RI:
				cpt = cpt.rescale(0, 2);
				sectPartRates = meanU3.getSolution().calcTotParticRateForAllSects();
				break;
			case SOLID:
				cpt = cpt.rescale(0, 5);
				for (Integer parentID : FaultArrayCalc.ELSINORE_PARENTS)
					parentColorVals.put(parentID, 2d);
				for (Integer parentID : FaultArrayCalc.SJC_PARENTS)
					parentColorVals.put(parentID, 3d);
				for (Integer parentID : FaultArrayCalc.S_SAF_PARENTS)
					parentColorVals.put(parentID, 4d);
				break;
			case ALL_GRAY:
				cpt = new CPT(0, 1, Color.GRAY, Color.GRAY);
				break;

			default:
				throw new IllegalStateException("Not Yet Implemented: "+scalarType);
			}
			
			for (FaultSection sect : sects) {
				Color color;
				if (scalarType == ScalarType.ALL_GRAY) {
					color = Color.DARK_GRAY;
				} else {
					double val;
					switch (scalarType) {
					case SUBSECT_PARTIC_PROB:
						val = Math.log10(FaultSysSolutionERF_Calc.calcParticipationProbForSect(meanU3, magForProbs, sect.getSectionId()));
						break;
					case PARENT_PARTIC_PROB:
						val = Math.log10(parentProbs.get(sect.getParentSectionId()));
						break;
					case TIME_SINCE_LAST:
						Long timeLast = sect.getDateOfLastEvent();
						if (timeLast == Long.MIN_VALUE)
							val = 1000d;
						else
							val = (System.currentTimeMillis() - timeLast)/ProbabilityModelsCalc.MILLISEC_PER_YEAR;
						break;
					case SUBSECT_NORM_RI:
						if (sect.getDateOfLastEvent() == Long.MIN_VALUE) {
							val = 1d;
						} else {
							double timeSince = (System.currentTimeMillis() - sect.getDateOfLastEvent())/ProbabilityModelsCalc.MILLISEC_PER_YEAR;
							double ri = 1d/sectPartRates[sect.getSectionId()];
							val = timeSince / ri;
						}
						break;
					case SOLID:
						val = parentColorVals.get(sect.getParentSectionId());
						break;

					default:
						throw new IllegalStateException("Not Yet Implemented: "+scalarType);
					}
					color = cpt.getColor((float)val);
				}
				for (PSXYPolygon poly : FaultBasedMapGen.getPolygons(sect.getFaultTrace(), color, 2.5))
					map.addPolys(poly);
			}
		}
		
		if (plotLinearTrace) {
			for (PSXYPolygon poly : FaultBasedMapGen.getPolygons(FaultArrayCalc.getS_SAF_LinearTrace(allParents), Color.BLACK, 1))
				map.addPolys(poly);
			for (PSXYPolygon poly : FaultBasedMapGen.getPolygons(FaultArrayCalc.getSJC_LinearTrace(allParents), Color.BLACK, 1))
				map.addPolys(poly);
			for (PSXYPolygon poly : FaultBasedMapGen.getPolygons(FaultArrayCalc.getElsinoreLinearTrace(allParents), Color.BLACK, 1))
				map.addPolys(poly);
		}
		
		List<LocationList> arrayLocLists = new ArrayList<>();
		arrayLocLists.add(FaultArrayCalc.getS_SAF_LinearTrace(allParents));
		arrayLocLists.add(FaultArrayCalc.getSJC_LinearTrace(allParents));
		if (region.getMaxLat() > 38) {
			arrayLocLists.add(FaultArrayCalc.getN_SAF_LinearTrace(allParents));
			arrayLocLists.add(FaultArrayCalc.getHaywardCalaverasLinearTrace(allParents));
		}
		
		if (arraySpacingDashed > 0) {
			for (LocationList trace : arrayLocLists) {
				addArrayPolys(map, trace, arraySpacingDashed, Color.BLACK, PlotLineType.DASHED, dashed_array_plot_len_km);
			}
//			addArrayPolys(map, FaultArrayCalc.getS_SAF_LinearTrace(allParents), arraySpacingDashed, Color.GREEN.darker(), PlotLineType.DASHED);
//			addArrayPolys(map, FaultArrayCalc.getSJC_LinearTrace(allParents), arraySpacingDashed, Color.MAGENTA.darker(), PlotLineType.DASHED);
//			addArrayPolys(map, FaultArrayCalc.getElsinoreLinearTrace(allParents), arraySpacing, Color.BLUE);
		}
		if (arraySpacingSolid > 0d) {
			for (LocationList trace : arrayLocLists) {
				addArrayPolys(map, trace, arraySpacingSolid, Color.GREEN.darker(), PlotLineType.SOLID, array_plot_len_km);
			}
		}
		
		if (region.getMaxLat() < 38) {
			// so cal zoomed
			for (String city : cities.keySet()) {
				Location loc = cities.get(city);
				Point2D.Double pt = new Point2D.Double(loc.getLongitude(), loc.getLatitude());
				Point2D.Double textPT = new Point2D.Double(loc.getLongitude()+0.08, loc.getLatitude());
				
				map.addSymbol(new PSXYSymbol(pt, PSXYSymbol.Symbol.CIRCLE, 0.05, 0.01, Color.BLACK, Color.BLACK));
				map.addText(new PSText(textPT, Color.BLACK, 10, city, Justify.LEFT_BOTTOM));
			}
		}
		
		FaultBasedMapGen.plotMap(outputDir, prefix, false, map);
		convertHiRes(outputDir, prefix);
		if (scalarType != null) {
			// also write version with just CPT files
			map.setGriddedData(null);
			map.setCpt(cpt);
			map.setCustomGRDPath(null);
			map.setCustomIntenPath(null);
			map.setCustomScaleMin((double)cpt.getMinValue());
			map.setCustomScaleMax((double)cpt.getMaxValue());
			map.setCustomLabel(scalarType.label);
			FaultBasedMapGen.plotMap(outputDir, prefix+"_cpt", false, map);
			convertHiRes(outputDir, prefix+"_cpt");
		}
	}
	
	private static void convertHiRes(File outputDir, String prefix) throws IOException {
		File pdfFile = new File(outputDir, prefix+".pdf");
		File hiresFile = new File(outputDir, prefix+"_hires.png");
		// convert hi res
		String[] command = { "/bin/bash", "-c", "convert -density 300 "
				+pdfFile.getAbsolutePath()+" "+hiresFile.getAbsolutePath() };
		
		Process p = Runtime.getRuntime().exec(command);
		try {
			int exit = p.waitFor();
			System.out.println("Converted, exit code: "+exit);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
	}
	
	static Map<String, Location> cities;
	static {
		cities = new HashMap<>();
		// San Diego, Los Angeles, Santa Barbara, Victorville, Palm Springs
		cities.put("San Diego", new Location(32.70, -117.15));
		cities.put("Los Angeles", new Location(34.05, -118.25));
		cities.put("Santa Barbara", new Location(34.45, -119.70));
		cities.put("Victorville", new Location(34.536376, -117.292734));
		cities.put("Palm Springs", new Location(33.830368, -116.545593));
		cities.put("San Luis Obispo", new Location(35.282638, -120.659932));
		
	}
	
	private static final double dashed_array_plot_len_km = 15;
	private static final double array_plot_len_km = 30;
//	private static final Color array_color = Color.RED.darker();
//	private static final Color array_color = new Color(0, 200, 70);
	
	private static void addArrayPolys(GMT_Map map, LocationList linearTrace, double arraySpacing, Color color,
			PlotLineType lineType, double arrayLen) {
		List<Location> arrayLocs = FaultArrayCalc.calcArrayLocations(linearTrace, arraySpacing);
		
		for (Location arrayLoc : arrayLocs) {
			double normalAz = FaultArrayCalc.calcFaultNormalAzimuth(linearTrace, arrayLoc, 5d);
			
			LocationList line = new LocationList();
			line.add(LocationUtils.location(arrayLoc, Math.toRadians(normalAz), 0.5*arrayLen));
			line.add(LocationUtils.location(arrayLoc, Math.toRadians(normalAz), -0.5*arrayLen));
			for (PSXYPolygon poly : FaultBasedMapGen.getPolygons(line, color, 1, lineType))
				map.addPolys(poly);
		}
	}

	public static void main(String[] args) throws IOException, GMT_MapException {
//		File outputDir = new File("/home/kevin/SCEC/2021_fault_array_proposal/plots");
		File outputDir = new File("/home/kevin/SCEC/2023_fault_array_proposal/plots");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdirs());
		FaultBasedMapGen.LOCAL_MAPGEN = true;
		
		FaultModels fm = FaultModels.FM3_1;
		MeanUCERF3 meanU3 = new MeanUCERF3();
		meanU3.setPreset(Presets.FM3_1_MAG_VAR);
		meanU3.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_PREF_BLEND);
		meanU3.getTimeSpan().setStartTime(2019);
		meanU3.getTimeSpan().setDuration(durationForProbs);
		meanU3.updateForecast();
		
		boolean topo = false;
		boolean popData = true;

//		plotMap(outputDir, "linearized_trace_test", 0d, fm, null, null, true);
//		plotMap(outputDir, "linearized_trace_20km", 20d, fm, null, null, true);
//		plotMap(outputDir, "linearized_trace_40km", 40d, fm, null, null, true);
		
		ScalarType[] types = {
				ScalarType.PARENT_PARTIC_PROB,
//				ScalarType.SUBSECT_PARTIC_PROB,
//				ScalarType.SOLID,
//				ScalarType.TIME_SINCE_LAST,
//				ScalarType.ALL_GRAY,
				};
		for (ScalarType type : types) {
			String prefix = type.name().toLowerCase();
//			plotMap(outputDir, prefix+"_20km", 20d, fm, meanU3, type, false, srtm, popData);
//			plotMap(outputDir, prefix+"_30km", 30d, fm, meanU3, type, false, srtm, popData);
//			plotMap(outputDir, prefix+"_40km", 40d, fm, meanU3, type, false, srtm, popData);
//			plotMap(outputDir, prefix+"_60km", 60d, fm, meanU3, type, false, srtm, popData);
			
			plotMap(outputDir, prefix+"_80_40km", 80d, 40d, fm, meanU3, type, false, topo, popData);
		}
	}

}
