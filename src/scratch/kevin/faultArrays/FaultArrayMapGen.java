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
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.elements.CoastAttributes;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.mapping.gmt.elements.PSText;
import org.opensha.commons.mapping.gmt.elements.PSText.Justify;
import org.opensha.commons.mapping.gmt.elements.PSXYPolygon;
import org.opensha.commons.mapping.gmt.elements.PSXYSymbol;
import org.opensha.commons.mapping.gmt.elements.TopographicSlopeFile;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Ints;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.UCERF3.analysis.FaultSysSolutionERF_Calc;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.erf.mean.MeanUCERF3;
import scratch.UCERF3.erf.mean.MeanUCERF3.Presets;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;

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
	private static Region region = new Region(new Location(35.5, -114), new Location(31.5, -121));
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
		if (topoCPT == null)
			topoCPT = CPT.loadFromFile(new File("/home/kevin/SCEC/2019_fault_array_proposal/origGMT/test.cpt"));
		return topoCPT;
	}
	
	private static final double magForProbs = 6.5;
	private static final double durationForProbs = 30d;
	
	private static void plotMap(File outputDir, String prefix, double arraySpacing, FaultModels fm,
			MeanUCERF3 meanU3, ScalarType scalarType, boolean plotLinearTrace) throws IOException, GMT_MapException {
//		TopographicSlopeFile topoRes = TopographicSlopeFile.SRTM_30_PLUS;
//		GriddedGeoDataSet gridData = fetchTopo(topoRes);
		CPT topoCPT = getTopoCPT();
		GMT_Map map = new GMT_Map(region, null, 3d/3600d, topoCPT);
		
		map.setBlackBackground(false);
		map.setRescaleCPT(false);
		map.setCustomScaleMin((double)topoCPT.getMinValue());
		map.setCustomScaleMax((double)topoCPT.getMaxValue());
//		map.setCoast(new CoastAttributes(Color.BLACK, 0.6d, new Color(160, 200, 200)));
		map.setCoast(new CoastAttributes(Color.BLACK, 0.6d, new Color(112, 131, 147)));
		map.setCustomLabel(null);
		map.setUseGMTSmoothing(true);
//		map.setTopoResolution(topoRes);
		map.setCustomIntenPath("/home/kevin/SCEC/2019_fault_array_proposal/origGMT/socal_mex.grad");
		map.setCustomGRDPath("/home/kevin/SCEC/2019_fault_array_proposal/origGMT/socal_mex.grd");
//		map.setLabelSize(MAP_LABEL_SIZE);
//		map.setLabelTickSize(MAP_LABEL_TICK_SIZE);
//		map.setlabel
		map.setDpi(150);
//		map.setImageWidth(imageWidth);
		map.setDrawScaleKM(false);
		map.setGMT_Param("MAP_FRAME_TYPE", "FANCY");
		map.setGMT_Param("MAP_FRAME_WIDTH", "0.04i");
		map.setGMT_Param("MAP_TICK_LENGTH_PRIMARY", "0.04i");
		
		for (FaultSectionPrefData sect : meanU3.getSolution().getRupSet().getFaultSectionDataList())
			for (PSXYPolygon poly : FaultBasedMapGen.getPolygons(sect.getFaultTrace(), Color.BLACK, 0.3))
				map.addPolys(poly);
		
		Map<Integer, FaultSectionPrefData> allParents = fm.fetchFaultSectionsMap();
		
		CPT cpt = null;
		
		if (scalarType != null) {
			List<FaultSectionPrefData> sects = new ArrayList<>();
			HashSet<Integer> parentsToPlot = new HashSet<>();
			parentsToPlot.addAll(Ints.asList(FaultArrayCalc.SAF_PARENTS));
			parentsToPlot.addAll(Ints.asList(FaultArrayCalc.SJC_PARENTS));
			parentsToPlot.addAll(Ints.asList(FaultArrayCalc.ELSINORE_PARENTS));
			for (FaultSectionPrefData sect : meanU3.getSolution().getRupSet().getFaultSectionDataList())
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
				for (Integer parentID : FaultArrayCalc.SAF_PARENTS)
					parentColorVals.put(parentID, 4d);
				break;
			case ALL_GRAY:
				cpt = new CPT(0, 1, Color.GRAY, Color.GRAY);
				break;

			default:
				throw new IllegalStateException("Not Yet Implemented: "+scalarType);
			}
			
			for (FaultSectionPrefData sect : sects) {
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
			for (PSXYPolygon poly : FaultBasedMapGen.getPolygons(FaultArrayCalc.getSAF_LinearTrace(allParents), Color.BLACK, 1))
				map.addPolys(poly);
			for (PSXYPolygon poly : FaultBasedMapGen.getPolygons(FaultArrayCalc.getSJC_LinearTrace(allParents), Color.BLACK, 1))
				map.addPolys(poly);
			for (PSXYPolygon poly : FaultBasedMapGen.getPolygons(FaultArrayCalc.getElsinoreLinearTrace(allParents), Color.BLACK, 1))
				map.addPolys(poly);
		}
		
		if (arraySpacing > 0) {
			addArrayPolys(map, FaultArrayCalc.getSAF_LinearTrace(allParents), arraySpacing, Color.GREEN.darker());
			addArrayPolys(map, FaultArrayCalc.getSJC_LinearTrace(allParents), arraySpacing, Color.MAGENTA.darker());
			addArrayPolys(map, FaultArrayCalc.getElsinoreLinearTrace(allParents), arraySpacing, Color.BLUE);
		}
		
		for (String city : cities.keySet()) {
			Location loc = cities.get(city);
			Point2D.Double pt = new Point2D.Double(loc.getLongitude(), loc.getLatitude());
			Point2D.Double textPT = new Point2D.Double(loc.getLongitude()+0.08, loc.getLatitude());
			
			map.addSymbol(new PSXYSymbol(pt, PSXYSymbol.Symbol.CIRCLE, 0.05, 0.01, Color.BLACK, Color.BLACK));
			map.addText(new PSText(textPT, Color.BLACK, 10, city, Justify.LEFT_BOTTOM));
		}
		
		FaultBasedMapGen.plotMap(outputDir, prefix, false, map);
		if (scalarType != null) {
			map.setGriddedData(null);
			map.setCpt(cpt);
			map.setCustomGRDPath(null);
			map.setCustomIntenPath(null);
			map.setCustomScaleMin((double)cpt.getMinValue());
			map.setCustomScaleMax((double)cpt.getMaxValue());
			map.setCustomLabel(scalarType.label);
			FaultBasedMapGen.plotMap(outputDir, prefix+"_cpt", false, map);
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
	
	private static final double array_plot_len_km = 15;
//	private static final Color array_color = Color.RED.darker();
//	private static final Color array_color = new Color(0, 200, 70);
	
	private static void addArrayPolys(GMT_Map map, LocationList linearTrace, double arraySpacing, Color color) {
		List<Location> arrayLocs = FaultArrayCalc.calcArrayLocations(linearTrace, arraySpacing);
		
		for (Location arrayLoc : arrayLocs) {
			double normalAz = FaultArrayCalc.calcFaultNormalAzimuth(linearTrace, arrayLoc, 5d);
			
			LocationList line = new LocationList();
			line.add(LocationUtils.location(arrayLoc, Math.toRadians(normalAz), 0.5*array_plot_len_km));
			line.add(LocationUtils.location(arrayLoc, Math.toRadians(normalAz), -0.5*array_plot_len_km));
			for (PSXYPolygon poly : FaultBasedMapGen.getPolygons(line, color, 1))
				map.addPolys(poly);
		}
	}

	public static void main(String[] args) throws IOException, GMT_MapException {
		File outputDir = new File("/home/kevin/SCEC/2019_fault_array_proposal/plots");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		FaultBasedMapGen.LOCAL_MAPGEN = true;
		
		FaultModels fm = FaultModels.FM3_1;
		MeanUCERF3 meanU3 = new MeanUCERF3();
		meanU3.setPreset(Presets.FM3_1_MAG_VAR);
		meanU3.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_PREF_BLEND);
		meanU3.getTimeSpan().setStartTime(2019);
		meanU3.getTimeSpan().setDuration(durationForProbs);
		meanU3.updateForecast();

//		plotMap(outputDir, "linearized_trace_test", 0d, fm, null, null, true);
//		plotMap(outputDir, "linearized_trace_20km", 20d, fm, null, null, true);
//		plotMap(outputDir, "linearized_trace_40km", 40d, fm, null, null, true);
		
		ScalarType[] types = {
				ScalarType.PARENT_PARTIC_PROB,
				ScalarType.SUBSECT_PARTIC_PROB,
				ScalarType.SOLID,
				ScalarType.TIME_SINCE_LAST,
				ScalarType.ALL_GRAY,
				};
		for (ScalarType type : types) {
			String prefix = type.name().toLowerCase();
			plotMap(outputDir, prefix+"_20km", 20d, fm, meanU3, type, false);
			plotMap(outputDir, prefix+"_40km", 40d, fm, meanU3, type, false);
		}
	}

}
