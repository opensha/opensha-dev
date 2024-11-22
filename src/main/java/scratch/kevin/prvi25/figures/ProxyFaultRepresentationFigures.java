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
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
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
		rupSet.addModule(ProxyFaultSectionInstances.build(rupSet, 5, 5d));
		
//		int parentID = FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "Anegada", "SW");
//		Region region = new Region(new Location(17.25, -66.25), new Location(18.5, -64));
//		int parentID = FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "SW", "Puerto", "Rico");
//		Region region = new Region(new Location(17.7, -67.05), new Location(18.1, -66.6));
		int parentID = FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "Cerro", "Goden");
		Region region = new Region(new Location(18.08, -67.4), new Location(18.4, -66.95));
		
		String parentName = null;
		
		CPT colors = GMT_CPT_Files.CATEGORICAL_BATLOW_UNIFORM.instance();
		
		List<FaultSection> proxySects = new ArrayList<>();
		for (FaultSection sect : rupSet.getFaultSectionDataList()) {
			if (sect.getParentSectionId() == parentID) {
				proxySects.add(sect);
				if (parentName == null)
					parentName = sect.getParentSectionName();
			}
		}
		
		double lengthLimit = PRVI25_InvConfigFactory.MAX_PROXY_FAULT_RUP_LEN_DEFAULT*1e3; // km -> m
		
		int rupID = -1;
		double maxRupMag = 0d;
		List<Integer> allRups = new ArrayList<>(rupSet.getRupturesForParentSection(parentID));
		Collections.reverse(allRups);
		for (int rupIndex : allRups) {
			double mag = rupSet.getMagForRup(rupIndex);
			if ((float)mag < (float)maxRupMag)
				continue;
			if (rupSet.getLengthForRup(rupIndex) > lengthLimit)
				continue;
			List<FaultSection> rupSects = rupSet.getFaultSectionDataForRupture(rupIndex);
//			if (rupSects.size() != proxySects.size())
//				continue;
			boolean allMatch = true;
			for (FaultSection sect : rupSects) {
				if (sect.getParentSectionId() != parentID) {
					allMatch = false;
					break;
				}
			}
			if (allMatch) {
				rupID = rupIndex;
				maxRupMag = mag;
			}
		}
		Preconditions.checkState(rupID >= 0);
		
		System.out.println(parentName);
		System.out.println("Using an M"+(float)maxRupMag+" rupture, len="+(float)rupSet.getLengthForRup(rupID)*1e-3);
		
		File outputDiur = new File("/tmp/");
		
		UniqueRupture rupUnique = UniqueRupture.forIDs(rupSet.getSectionsIndicesForRup(rupID));
		
		LocationList gridLocs = new LocationList();
		for (FaultSection sect : proxySects) {
			if (rupUnique.contains(sect.getSectionId())) {
				GriddedRegion grid = new GriddedRegion(sect.getZonePolygon(), 0.05d, GriddedRegion.ANCHOR_0_0);
				gridLocs.addAll(grid.getNodeList());
			}
		}
		
		ProxyFaultSectionInstances proxyModule = rupSet.requireModule(ProxyFaultSectionInstances.class);
		
		RuptureSurface rup = rupSet.getSurfaceForRupture(rupID, 1d);
		
		for (RupType type : RupType.values()) {
			GeographicMapMaker mapMaker = new GeographicMapMaker(region);
			
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
				mapMaker.plot(outputDiur, "zone_rups_proxy_fault", "Proxy Fault Rupture");
			} else if (type == RupType.GRIDDED) {
				mapMaker.setFaultSections(proxySects);
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				for (int i=0; i<gridLocs.size(); i++)
					chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 3f, color(colors, i)));
				mapMaker.plotScatters(gridLocs, chars, null);
				mapMaker.plot(outputDiur, "zone_rups_gridded", "Gridded Ruptures");
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
				mapMaker.plot(outputDiur, "zone_rups_split_proxy_faults", "Distributed Proxy Fault Rupture Instances");
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
				mapMaker.plot(outputDiur, "zone_rups_finite", "Finite Fault Ruptures");
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
				mapMaker.plot(outputDiur, "zone_rups_rand_strike_finite", "Finite Random Strike Fault Ruptures");
			}
		}
		System.out.println(colors.size());
	}
	
	private static Color color(CPT colors, int index) {
		index = index % colors.size();
		return colors.get(index).minColor;
	}

}
