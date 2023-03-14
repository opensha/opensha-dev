package scratch.kevin.miscFigures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.elements.CoastAttributes;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.mapping.gmt.elements.PSXYPolygon;
import org.opensha.commons.mapping.gmt.elements.TopographicSlopeFile;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.param.HistoricOpenIntervalParam;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.UCERF3.analysis.FaultSysSolutionERF_Calc;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.utils.LastEventData;

public class U3vsPopulationMap {
	
	private enum Backgrounds {
		POPULATION,
		POPULATION_PLUS_TOPO,
		NONE,
	}
	
	private enum Faults {
		NONE,
		BLACK,
		ALL_COLORED,
		HIGH_PROB_COLORED,
		MAJOR
	}
	
	private static final File popDataFile = new File("/home/kevin/Downloads/gpw-v4/gpw_v4_population_density_rev11_2020_30_sec_1.asc");
	private static final int ncols = 10800;
	private static final int nrows = 10800;
	private static final double xll = -180;
	private static final double yll = -4.2632564145606e-14;
	public static final double cellsize = 0.0083333333333333;
	
	public static GriddedGeoDataSet fetchPopData(Region mapRegion, double mapSpacing) throws NumberFormatException, IOException {
		double[][] data = new double[nrows][];
		int row = nrows-1;
		
		int skippedRows = 0;
		
		MinMaxAveTracker track = new MinMaxAveTracker();
		int numFirstDataVals = 0;
		for (String line : Files.readLines(popDataFile, Charset.defaultCharset())) {
			if (row == 0 && !line.startsWith("-9999"))
				continue;
			
			double lat = yll + cellsize*(double)row;
			
			if (!line.startsWith("-9999"))
				numFirstDataVals++;
			
			boolean skip = lat > mapRegion.getMaxLat()+5*cellsize || lat < mapRegion.getMinLat()-5*cellsize;
			if (skip) {
				skippedRows++;
			} else {
				data[row] = new double[ncols];
				String[] strs = line.trim().split(" ");
				Preconditions.checkState(strs.length == ncols);
				for (int col=0; col<ncols; col++) {
					double val = Double.parseDouble(strs[col]);
					if (val < 0)
						val = 0;
					track.addValue(val);
					data[row][col] = val;
				}
			}
			row--;
		}
		System.out.println("done reading data, skipped "+skippedRows+"/"+nrows+" rows");
		System.out.println("\t"+numFirstDataVals+" lines start with a real data val");
		System.out.println("\n\t"+track);
		
		GriddedRegion gridReg = new GriddedRegion(mapRegion, mapSpacing, GriddedRegion.ANCHOR_0_0);
		System.out.println("Region has "+gridReg.getNodeCount()+" grid vals");
		GriddedGeoDataSet popData = new GriddedGeoDataSet(gridReg, false);
		System.out.println("Building population gridded");
		for (int i=0; i<popData.size(); i++) {
			Location loc = gridReg.getLocation(i);
			
			int col = ((int)((loc.getLongitude() - xll) / cellsize + 0.5));
			row = ((int)((loc.getLatitude() - yll) / cellsize + 0.5));
			
			popData.set(i, data[row][col]);
		}
		
		return popData;
	}
	
	public static CPT getLogPopCPT() {
		CPT popCPT = new CPT(0d, 4d,
				new Color(224,243,248),
				new Color(171,217,233),
				new Color(116,173,209),
				new Color(69,117,180),
				new Color(49,54,149),
				new Color(24,27,74));
		Color zeroColor = new Color(230, 230, 230);
//		Color zeroColor = new Color(210, 210, 210);
//		Color zeroColor = Color.LIGHT_GRAY;
		popCPT.setNanColor(zeroColor);
		popCPT.setBelowMinColor(zeroColor);
		popCPT.setAboveMaxColor(popCPT.getMaxColor());
		return popCPT;
	}

	public static void main(String[] args) throws IOException, GMT_MapException {
		File outputDir = new File("/home/kevin/SCEC/2022_scec_proposal/population_faults_fig");
		
		Region mapRegion = new Region(new Location(31.75, -125), new Location(42.5, -114));
//		Region mapRegion = new Region(new Location(32.5, -121), new Location(36, -115));
		
		double mapSpacing = cellsize;
//		mapSpacing *= 5;
		
		GriddedGeoDataSet popData = fetchPopData(mapRegion, mapSpacing);
		GriddedRegion gridReg = popData.getRegion();
		
//		CPT popCPT = GMT_CPT_Files.BLACK_RED_YELLOW_UNIFORM.instance().reverse();
//		popCPT = popCPT.rescale(0, 4);
//		for (int i=0; i<popCPT.size(); i++) {
//			double weight =  (double)i/(double)(popCPT.size()-1);
//			CPTVal val = popCPT.get(i);
//			val.minColor = blend(val.minColor, saturate(val.minColor), weight);
//			val.maxColor = blend(val.maxColor, saturate(val.maxColor), weight);
//		}
		CPT popCPT = getLogPopCPT();
		
		popData.log10();
		for (int i=0; i<popData.size(); i++)
			if (!Double.isFinite(popData.get(i)))
				popData.set(i, -1d);
		
//		CPT popCPT = new CPT(1e-10, 20000d, Color.GRAY, Color.YELLOW, Color.RED, Color.RED.darker().darker().darker());
//		popCPT.add(0, new CPTVal(0f, Color.GRAY, (float)popCPT.getMinValue(), Color.GRAY));
//		popCPT.setAboveMaxColor(popCPT.getMaxColor());
		
		FaultSystemSolution sol = FaultSystemSolution.load(
				new File("/home/kevin/OpenSHA/UCERF3/rup_sets/modular/FM3_1_branch_averaged.zip"));
		LastEventData.populateSubSects(sol.getRupSet().getFaultSectionDataList(), LastEventData.load());
		
		boolean logProbs = true;
//		CPT probCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-2d, 0d);
		CPT probCPT = new CPT(Math.log10(0.05), Math.log10(0.5),
				new Color(254,224,144),
				new Color(253,174,97),
				new Color(244,109,67),
				new Color(215,48,39),
				new Color(165,0,38),
				new Color(82,0,19));
		probCPT.setBelowMinColor(probCPT.getMinColor());
		probCPT.setAboveMaxColor(probCPT.getMaxColor());
		probCPT.setNanColor(Color.WHITE);
		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_PREF_BLEND);
		erf.getTimeSpan().setDuration(30d);
		erf.setParameter(HistoricOpenIntervalParam.NAME, 2022d-1875d);
		erf.updateForecast();
		
		double minMag = 6.5;
		
		double[] faultProbs = new double[sol.getRupSet().getNumSections()];
		for (FaultSection sect : sol.getRupSet().getFaultSectionDataList())
			faultProbs[sect.getSectionId()] = FaultSysSolutionERF_Calc.calcParticipationProbForSect(
					erf, minMag, sect.getSectionId());
		
		boolean includeFullParents = true;
		HashMap<Integer, Double> maxProbForParents = new HashMap<>();
		for (FaultSection sect : sol.getRupSet().getFaultSectionDataList()) {
			int parent = sect.getParentSectionId();
			if (maxProbForParents.containsKey(parent))
				maxProbForParents.put(parent, Math.max(maxProbForParents.get(parent), faultProbs[sect.getSectionId()]));
			else
				maxProbForParents.put(parent, faultProbs[sect.getSectionId()]);
		}
		
		List<FaultSection> sectsSorted = new ArrayList<>(sol.getRupSet().getFaultSectionDataList());
		sectsSorted.sort(new Comparator<FaultSection>() {

			@Override
			public int compare(FaultSection o1, FaultSection o2) {
				return Double.compare(faultProbs[o1.getSectionId()], faultProbs[o2.getSectionId()]);
			}
		});
		
		HashSet<Integer> majorParents = new HashSet<>();
		Map<String, List<Integer>> namedMap = FaultModels.FM3_1.getNamedFaultsMapAlt();
		for (String name : namedMap.keySet()) {
			if (name.equals("Compton"))
				continue;
			if (name.equals("Puente Hills"))
				continue;
			if (name.equals("Little Salmon"))
				continue;
			if (name.equals("Green Valley"))
				continue;
			if (name.equals("San Gregorio"))
				continue;
			majorParents.addAll(namedMap.get(name));
		}
		for (FaultSection sect : sectsSorted) {
			String parentName = sect.getParentSectionName().toLowerCase();
			if (parentName.contains("andreas") || parentName.contains("imperial") || parentName.contains("brawley")
					|| parentName.contains("mendocino"))
				majorParents.add(sect.getParentSectionId());
		}
		
		GriddedGeoDataSet fakeData = new GriddedGeoDataSet(new GriddedRegion(mapRegion, mapSpacing, GriddedRegion.ANCHOR_0_0), false);
		for (int i=0; i<fakeData.size(); i++)
			fakeData.set(i, Double.NaN);
		
//		double maxLat = Math.rint(((maxTempLat-minLat)/mapSpacing))*mapSpacing +mapRegion.getMaxLat();
//		double maxLon = Math.rint(((maxTempLon-minLon)/mapSpacing))*mapSpacing +mapRegion.getMinLat();
		
		for (Backgrounds bkg : Backgrounds.values()) {
			for (Faults flts : Faults.values()) {
				if (flts == Faults.HIGH_PROB_COLORED || flts == Faults.MAJOR)
					probCPT = probCPT.rescale(Math.log10(0.01d), Math.log10(0.5d));
				else
					probCPT = probCPT.rescale(-3d, Math.log10(0.5d));
				GMT_Map map;
				CPT cpt;
				if (bkg == Backgrounds.POPULATION_PLUS_TOPO) {
					map = new GMT_Map(mapRegion, popData, gridReg.getSpacing(), popCPT);
					map.setTopoResolution(TopographicSlopeFile.SRTM_30_PLUS);
					cpt = popCPT;
					map.setCPTCustomInterval(1d);
				} else if (bkg == Backgrounds.POPULATION) {
					map = new GMT_Map(mapRegion, popData, gridReg.getSpacing(), popCPT);
					map.setTopoResolution(null);
//					map.setUseGMTSmoothing(false);
					cpt = popCPT;
					map.setCPTCustomInterval(1d);
				} else {
					map = new GMT_Map(mapRegion, fakeData, fakeData.getRegion().getSpacing(), probCPT);
					map.setTopoResolution(null);
					map.setUseGMTSmoothing(false);
					cpt = probCPT;
					map.setCPTCustomInterval(0.5d);
				}
				
				map.setDrawScaleKM(false);
				map.setBlackBackground(false);
				map.setRescaleCPT(false);
				map.setCustomScaleMin((double)cpt.getMinValue());
				map.setCustomScaleMax((double)cpt.getMaxValue());
//				map.setCoast(new CoastAttributes(Color.BLACK, 0.6d, new Color(160, 200, 200)));
				map.setCoast(new CoastAttributes(Color.BLACK, 0.6d, new Color(112, 131, 147)));
				map.setCustomLabel(null);
				map.setDpi(600);
				map.setLogPlot(false);
				if (flts != Faults.NONE) {
					for (FaultSection sect : sectsSorted) {
						double prob = FaultSysSolutionERF_Calc.calcParticipationProbForSect(erf, minMag, sect.getSectionId());
						if (logProbs)
							prob = Math.log10(prob);
						
						double thickness;
						Color color;
						
						boolean black = flts == Faults.BLACK || !Double.isFinite(prob);
						
						if (!black && flts == Faults.HIGH_PROB_COLORED) {
							if (includeFullParents) {
								double parentProb = maxProbForParents.get(sect.getParentSectionId());
								if (logProbs)
									parentProb = Math.log10(parentProb);
								black = parentProb < probCPT.getMinValue();
							} else {
								black = (float)prob < probCPT.getMinValue();
							}
						} else if (flts == Faults.MAJOR) {
							black = !majorParents.contains(sect.getParentSectionId());
						}
						
						if (black) {
							// simple line
//							thickness = 0.3;
							thickness = 0.8d;
							color = Color.BLACK;
						} else {
							thickness = 2.5;
							color = probCPT.getColor((float)prob);
						}
						
						for (PSXYPolygon poly : FaultBasedMapGen.getPolygons(sect.getFaultTrace(), color, thickness))
							map.addPolys(poly);
					}
				}
				if (flts == Faults.NONE && bkg == Backgrounds.NONE)
					continue;
				String label;
				if (bkg == Backgrounds.NONE)
					label = "Log10 M"+(float)minMag+" 30 year Probability";
				else
					label = "Log10 Population Per Cell";
				map.setCustomLabel(label);
				String prefix = "map_background_"+bkg.name()+"_faults_"+flts.name();
				System.out.println("Plotting map: "+prefix);
				FaultBasedMapGen.plotMap(outputDir, prefix, false, map);
				
				File pdfFile = new File(outputDir, prefix+".pdf");
				File hiresFile = new File(outputDir, prefix+"_hires.png");
				
				// convert hi res
				String[] command = { "/bin/bash", "-c", "convert -density 600 "
						+pdfFile.getAbsolutePath()+" "+hiresFile.getAbsolutePath() };
				
				Process p = Runtime.getRuntime().exec(command);
				try {
					int exit = p.waitFor();
					System.out.println("Converted, exit code: "+exit);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			
		}
		
		System.out.println("DONE");
	}
	
	private static final int saturation_steps = 1;
	
	private static Color saturate(Color c) {
		int r = c.getRed();
		int g = c.getGreen();
		int b = c.getBlue();
		
		for (int i=0; i<saturation_steps; i++) {
			r = (int)(0.5d*(r + 255d)+0.5);
			g = (int)(0.5d*(g + 255d)+0.5);
			b = (int)(0.5d*(b + 255d)+0.5);
		}
		
		return new Color(r, g, b, c.getAlpha());
	}
	
	private static Color blend(Color c1, Color c2, double weight) {
		float r = (float)((weight*c1.getRed() + (1d-weight)*c2.getRed())/255d);
		float g = (float)((weight*c1.getGreen() + (1d-weight)*c2.getGreen())/255d);
		float b = (float)((weight*c1.getBlue() + (1d-weight)*c2.getBlue())/255d);
		return new Color(r, g, b);
	}

}
