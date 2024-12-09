package scratch.kevin.prvi25.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.geo.json.FeatureProperties;
import org.opensha.commons.geo.json.Geometry;
import org.opensha.commons.geo.json.Geometry.LineString;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.modules.ProxyFaultSectionInstances;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.UniqueRupture;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.PRVI25_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTreeBranch;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;
import org.opensha.sha.faultSurface.RuptureSurface;

import com.google.common.base.Preconditions;

public class ProxyFaultRepresentationFigures {
	
	private enum RupType {
		PROXY,
		SPLIT_PROXIES,
		GRIDDED,
		FINITE,
		RAND_STRIKE_FINITE
	}

	public static void main(String[] args) throws IOException {
//		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/"
//				+ "batch_inversions/2024_05_21-prvi25_crustal_branches-proxyGriddedTests/"
//				+ "results_PRVI_FM_INITIAL_branch_averaged.zip"));
//		FaultSystemRupSet rupSet = sol.getRupSet();
		
		FaultSystemRupSet rupSet = new PRVI25_InvConfigFactory().buildRuptureSet(
				PRVI25_LogicTreeBranch.DEFAULT_CRUSTAL_ON_FAULT, FaultSysTools.defaultNumThreads());
		ProxyFaultSectionInstances proxyModule = ProxyFaultSectionInstances.build(rupSet, 5, 5d);
		rupSet.addModule(proxyModule);
		
		File outputDir = new File("/home/kevin/Documents/papers/2024_PRVI_ERF/prvi25-erf-paper/Figures/proxy_faults");
		
		boolean titles = false;
		
//		int singleParent = FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "Anegada", "SW");
		int singleParent = -1;
		
		double lengthLimit = 75d;
		
//		RupType[] types = RupType.values();
		RupType[] types = {
				RupType.PROXY,
				RupType.SPLIT_PROXIES
		};
		
		CPT colors = GMT_CPT_Files.CATEGORICAL_BATLOW_UNIFORM.instance();

		Map<Integer, MinMaxAveTracker[]> parentProxyRegTracks = new HashMap<>();
		Map<Integer, String> parentProxyNames = new HashMap<>();
		for (FaultSection sect : rupSet.getFaultSectionDataList()) {
			Region poly = sect.getZonePolygon();
			if (!sect.isProxyFault() || poly == null)
				continue;
			int parentID = sect.getParentSectionId();
			if (singleParent >= 0 && parentID != singleParent)
				continue;
			
			MinMaxAveTracker[] parentTracks = parentProxyRegTracks.get(parentID);
			if (parentTracks == null) {
				parentTracks = new MinMaxAveTracker[2];
				parentTracks[0] = new MinMaxAveTracker();
				parentTracks[1] = new MinMaxAveTracker();
				parentProxyRegTracks.put(parentID, parentTracks);
				parentProxyNames.put(parentID, sect.getParentSectionName());
			}
			parentTracks[0].addValue(poly.getMinLat());
			parentTracks[0].addValue(poly.getMaxLat());
			parentTracks[1].addValue(poly.getMinLon());
			parentTracks[1].addValue(poly.getMaxLon());
//			for (FaultSection proxySect : proxyModule.getProxySects()) {
//				if (proxySect.getParentSectionId() == sect.getSectionId()) {
//					for (Location loc : proxySect.getFaultSurface(1d).getPerimeter()) {
//						parentTracks[0].addValue(loc.lat);
//						parentTracks[1].addValue(loc.lon);
//					}
//				}
//			}
		}
		
		for (int parentID : parentProxyNames.keySet()) {
			String parentName = parentProxyNames.get(parentID);
			MinMaxAveTracker[] regTracks = parentProxyRegTracks.get(parentID);
			Location lowerLeft = new Location(regTracks[0].getMin(), regTracks[1].getMin());
			Location upperRight = new Location(regTracks[0].getMax(), regTracks[1].getMax());
			
			String prefix = parentName.replaceAll("\\W+", "_");
			while (parentName.startsWith("_"))
				parentName = parentName.substring(1);
			while (parentName.endsWith("_"))
				parentName = parentName.substring(0, parentName.length()-1);
			while (parentName.contains("__"))
				parentName = parentName.replace("__", "_");
			
			List<FaultSection> proxySects = new ArrayList<>();
			for (FaultSection sect : rupSet.getFaultSectionDataList())
				if (sect.getParentSectionId() == parentID)
					proxySects.add(sect);
			
			int rupID = -1;
			double maxRupMag = 0d;
			List<Integer> allRups = new ArrayList<>(rupSet.getRupturesForParentSection(parentID));
			Collections.reverse(allRups);
			for (int rupIndex : allRups) {
				double mag = rupSet.getMagForRup(rupIndex);
				if ((float)mag < (float)maxRupMag)
					continue;
				if (rupSet.getLengthForRup(rupIndex) > lengthLimit*1e3) // km to m
					continue;
				List<FaultSection> rupSects = rupSet.getFaultSectionDataForRupture(rupIndex);
//				if (rupSects.size() != proxySects.size())
//					continue;
				boolean allMatch = true;
				for (FaultSection sect : rupSects) {
					if (sect.getParentSectionId() != parentID) {
						allMatch = false;
						break;
					}
				}
//				System.out.println("Rup "+rupIndex+" is M"+(float)mag+" has "+rupSects.size()+" sects; allMatch? "+allMatch);
				if (allMatch) {
					rupID = rupIndex;
					maxRupMag = mag;
				}
			}
			Preconditions.checkState(rupID >= 0, "No matching rup found for %s. %s; tried %s on parent.",
					parentID, parentName, allRups.size());
			
			System.out.println(parentName);
			System.out.println("Using an M"+(float)maxRupMag+" rupture, len="+(float)rupSet.getLengthForRup(rupID)*1e-3);
			
			double rupLen = rupSet.getLengthForRup(rupID)*1e-3;
			double regBuffer = Math.min(60d, Double.max(30d, rupLen+5d));
			lowerLeft = LocationUtils.location(lowerLeft, 5d*Math.PI/4d, regBuffer);
			upperRight = LocationUtils.location(upperRight, Math.PI/4d, regBuffer);
			
			Region region = new Region(lowerLeft, upperRight);
			
			UniqueRupture rupUnique = UniqueRupture.forIDs(rupSet.getSectionsIndicesForRup(rupID));
			
			LocationList gridLocs = new LocationList();
			for (FaultSection sect : proxySects) {
				if (rupUnique.contains(sect.getSectionId())) {
					GriddedRegion grid = new GriddedRegion(sect.getZonePolygon(), 0.05d, GriddedRegion.ANCHOR_0_0);
					gridLocs.addAll(grid.getNodeList());
				}
			}
			
			RuptureSurface rup = rupSet.getSurfaceForRupture(rupID, 1d);
			
			for (RupType type : types) {
				GeographicMapMaker mapMaker = new GeographicMapMaker(region);
				mapMaker.setWriteGeoJSON(false);
				mapMaker.setScalarThickness(3f);
				
				if (type == RupType.PROXY) {
					mapMaker.setFaultSections(proxySects);
					Color color = color(colors, 0);
					List<Color> sectColors = new ArrayList<>();
					for (int i=0; i<proxySects.size(); i++) {
						if (rupUnique.contains(proxySects.get(i).getSectionId()))
							sectColors.add(color);
						else
							sectColors.add(null);
					}
					mapMaker.plotSectColors(sectColors);
					mapMaker.setScalarThickness(5f);
					mapMaker.plot(outputDir, prefix+"_proxy_fault", titles ? "Original Proxy Fault Rupture" : " ");
				} else if (type == RupType.GRIDDED) {
					mapMaker.setFaultSections(proxySects);
					List<PlotCurveCharacterstics> chars = new ArrayList<>();
					for (int i=0; i<gridLocs.size(); i++)
						chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 3f, color(colors, i)));
					mapMaker.plotScatters(gridLocs, chars, null);
					mapMaker.plot(outputDir, prefix+"_gridded", titles ? "Gridded Ruptures" : " ");
				} else if (type == RupType.SPLIT_PROXIES) {
					List<FaultSection> combFaultSects = new ArrayList<>();
					List<Color> sectColors = new ArrayList<>();
					Map<Integer, Integer> subProxyIDs = new HashMap<>();
					for (FaultSection sect : proxyModule.getProxySects()) {
						if (sect.getParentSectionName().startsWith(parentName)) {
							subProxyIDs.put(sect.getSectionId(), combFaultSects.size());
							combFaultSects.add(sect);
						}
					}
					combFaultSects.addAll(proxySects);
					mapMaker.setFaultSections(combFaultSects);
					for (int i=0; i<combFaultSects.size(); i++)
						sectColors.add(null);
					
					List<List<Integer>> rupProxyIndexes = proxyModule.getRupProxySectIndexes(rupID);
					for (int i=0; i<rupProxyIndexes.size(); i++) {
						Color color = color(colors, i);
						for (int index : rupProxyIndexes.get(i))
							sectColors.set(subProxyIDs.get(index), color);
					}
					mapMaker.setSkipNaNs(true);
					mapMaker.plotSectColors(sectColors);
					mapMaker.plot(outputDir, prefix+"_split_proxy_faults", titles ? "Distributed Proxy Fault Rupture Instances" : " ");
				} else if (type == RupType.FINITE) {
					List<FaultSection> combFaultSects = new ArrayList<>();
					List<Color> sectColors = new ArrayList<>();
					mapMaker.setSkipNaNs(true);
					
					double len = rup.getAveLength();
					double dip = rup.getAveDip();
					double strike = rup.getAveStrike();
					for (int i=0; i<gridLocs.size(); i++) {
						Location loc = gridLocs.get(i);
						Location first = LocationUtils.location(loc, new LocationVector(strike, -len/2d, 0d));
						Location last = LocationUtils.location(loc, new LocationVector(strike, len/2d, 0d));
						Geometry geometry = new LineString(first, last);
						FeatureProperties props = new FeatureProperties();
						props.set(GeoJSONFaultSection.FAULT_ID, i);
						props.set(GeoJSONFaultSection.UPPER_DEPTH, rup.getAveRupTopDepth());
						props.set(GeoJSONFaultSection.LOW_DEPTH, rup.getEvenlyDiscritizedLowerEdge().get(0).depth);
						props.set(GeoJSONFaultSection.DIP, dip);
						Feature feature = new Feature(i, geometry, props);
						GeoJSONFaultSection sect = GeoJSONFaultSection.fromFeature(feature);
						combFaultSects.add(sect);
						sectColors.add(color(colors, i));
					}
					combFaultSects.addAll(proxySects);
					while (sectColors.size() < combFaultSects.size())
						sectColors.add(null);
					mapMaker.setSkipNaNs(true);
					mapMaker.setFaultSections(combFaultSects);
					mapMaker.plotSectColors(sectColors);
					mapMaker.plot(outputDir, prefix+"_finite", titles ? "Finite Fault Ruptures" : " ");
				} else if (type == RupType.RAND_STRIKE_FINITE) {
					List<FaultSection> combFaultSects = new ArrayList<>();
					List<Color> sectColors = new ArrayList<>();
					mapMaker.setSkipNaNs(true);
					
					double len = rup.getAveLength();
					double dip = rup.getAveDip();
					double strike = rup.getAveStrike();
					
					Random rand = new Random(proxySects.size());
					for (int i=0; i<gridLocs.size(); i++) {
						Location loc = gridLocs.get(i);
						strike = rand.nextDouble()*360d;
						for (int j=0; j<2; j++) {
							strike += 90d;
							int id = combFaultSects.size();
							Location first = LocationUtils.location(loc, new LocationVector(strike, -len/2d, 0d));
							Location last = LocationUtils.location(loc, new LocationVector(strike, len/2d, 0d));
							Geometry geometry = new LineString(first, last);
							FeatureProperties props = new FeatureProperties();
							props.set(GeoJSONFaultSection.FAULT_ID, id);
							props.set(GeoJSONFaultSection.UPPER_DEPTH, rup.getAveRupTopDepth());
							props.set(GeoJSONFaultSection.LOW_DEPTH, rup.getEvenlyDiscritizedLowerEdge().get(0).depth);
							props.set(GeoJSONFaultSection.DIP, dip);
							Feature feature = new Feature(id, geometry, props);
							GeoJSONFaultSection sect = GeoJSONFaultSection.fromFeature(feature);
							sectColors.add(color(colors, id));
							combFaultSects.add(sect);
						}
					}
					combFaultSects.addAll(proxySects);
					while (sectColors.size() < combFaultSects.size())
						sectColors.add(null);
					mapMaker.setSkipNaNs(true);
					mapMaker.setFaultSections(combFaultSects);
					mapMaker.plotSectColors(sectColors);
					mapMaker.plot(outputDir, prefix+"_rand_strike_finite", titles ? "Finite Random Strike Fault Ruptures" : " ");
				}
			}
		}
	}
	
	private static Color color(CPT colors, int index) {
		index = index % colors.size();
		return colors.get(index).minColor;
	}

}
