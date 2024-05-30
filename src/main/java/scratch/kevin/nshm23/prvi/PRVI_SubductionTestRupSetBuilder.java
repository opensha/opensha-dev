package scratch.kevin.nshm23.prvi;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.geo.json.FeatureCollection;
import org.opensha.commons.geo.json.FeatureProperties;
import org.opensha.commons.geo.json.Geometry.MultiPoint;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets.*;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionDeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionFaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionScalingRelationships;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;

public class PRVI_SubductionTestRupSetBuilder {

	public static void main(String[] args) throws IOException {
		PRVI25_SubductionFaultModels fm = PRVI25_SubductionFaultModels.PRVI_SUB_FM_LARGE;
		PRVI25_SubductionDeformationModels dm = PRVI25_SubductionDeformationModels.FULL;
		// first just write out the fault model
		List<? extends FaultSection> sects = fm.getFaultSections();
		for (FaultSection sect : sects) {
			FaultTrace upper = sect.getFaultTrace();
			FaultTrace lower = sect.getLowerFaultTrace();
			System.out.println(sect.getSectionName());
			System.out.println("\tUpper trace: depth="+depthStr(upper)+"; strike="+oneDigit.format(upper.getAveStrike()));
			System.out.println("\tLower trace: depth="+depthStr(lower)+"; strike="+oneDigit.format(lower.getAveStrike()));
			System.out.println("\tdip="+(float)sect.getAveDip()+"; dipDir="+sect.getDipDirection());
		}
		GeographicMapMaker mapMaker = new GeographicMapMaker(sects);
		mapMaker.plot(new File("/tmp"), "subduction_sections", "PRVI Subduction Fault Model");
		PRVI25_SubductionScalingRelationships scale = PRVI25_SubductionScalingRelationships.LOGA_C4p0;
		List<? extends FaultSection> subSects = dm.build(fm);
//		RupSetConfig config = new CoulombRupSetConfig(subSects, scale);
//		RupSetConfig config = new SimpleAzimuthalRupSetConfig(subSects, scale);
		RupSetConfig config = new SimpleSubductionRupSetConfig(subSects, scale);
		
		FaultSystemRupSet rupSet = config.build(FaultSysTools.defaultNumThreads());
		
		for (FaultSection sect : subSects) {
			FaultTrace upper = sect.getFaultTrace();
			FaultTrace lower = sect.getLowerFaultTrace();
			System.out.println(sect.getSectionName());
			System.out.println("\tUpper trace: depth="+depthStr(upper)+"; strike="+oneDigit.format(upper.getAveStrike()));
			System.out.println("\tLower trace: depth="+depthStr(lower)+"; strike="+oneDigit.format(lower.getAveStrike()));
		}
		mapMaker = new GeographicMapMaker(subSects);
		mapMaker.plot(new File("/tmp"), "subduction_subsections", "PRVI Subduction Fault Model");
		
		rupSet.write(new File("/tmp/subduction_rup_set_"+dm.getFilePrefix()+".zip"));
		
		double minMag = Double.POSITIVE_INFINITY;
		double maxMag = Double.NEGATIVE_INFINITY;
		for (int rupIndex=0; rupIndex<rupSet.getNumRuptures(); rupIndex++) {
			double mag = rupSet.getMagForRup(rupIndex);
			minMag = Math.min(minMag, mag);
			maxMag = Math.max(maxMag, mag);
		}
		System.out.println("Rupture count: "+rupSet.getNumRuptures());
		System.out.println("Sect minimums:");
		for (int s=0; s<subSects.size(); s++) {
			int minNumSects = Integer.MAX_VALUE;
			double sectMinMag = Double.POSITIVE_INFINITY;
			for (int rupIndex : rupSet.getRupturesForSection(s)) {
				minNumSects = Integer.min(minNumSects, rupSet.getSectionsIndicesForRup(rupIndex).size());
				sectMinMag = Math.min(sectMinMag, rupSet.getMagForRup(rupIndex));
			}
			System.out.println("\t"+s+". "+rupSet.getFaultSectionData(s).getSectionName()+":\tM"+twoDigits.format(sectMinMag)+", "+minNumSects+" sects");
		}
		System.out.println("Mag range: ["+twoDigits.format(minMag)+", "+twoDigits.format(maxMag)+"]");
		
		// now write surface debug GeoJSON
		List<Feature> features = new ArrayList<>();
		features.addAll(FeatureCollection.read(new File("/tmp/subduction_subsections.geojson")).features);
		for (FaultSection sect : rupSet.getFaultSectionDataList()) {
			MultiPoint geom = new MultiPoint(sect.getFaultSurface(15d, false, false).getEvenlyDiscritizedListOfLocsOnSurface());
			Feature points = new Feature(geom, new FeatureProperties());
			features.add(points);
		}
		
		FeatureCollection.write(new FeatureCollection(features), new File("/tmp/subduction_subsection_points.geojson"));
	}
	
	private static String depthStr(FaultTrace trace) {
		MinMaxAveTracker track = new MinMaxAveTracker();
		for (Location loc : trace)
			track.addValue(loc.depth);
		return "avg="+oneDigit.format(track.getAverage())+" ["+oneDigit.format(track.getMin())+", "+oneDigit.format(track.getMax())+"]";
	}
	private static DecimalFormat oneDigit = new DecimalFormat("0.0");
	private static DecimalFormat twoDigits = new DecimalFormat("0.00");

}
