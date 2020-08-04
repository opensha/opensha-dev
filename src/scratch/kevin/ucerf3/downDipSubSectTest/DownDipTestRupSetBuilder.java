package scratch.kevin.ucerf3.downDipSubSectTest;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.geo.Location;
import org.opensha.commons.util.FaultUtils;
import org.opensha.commons.util.IDPairing;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRuptureBuilder;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.PlausibilityFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.strategies.ClusterConnectionStrategy;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.strategies.ClusterPermutationStrategy;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.strategies.DistCutoffClosestSectClusterConnectionStrategy;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SectionDistanceAzimuthCalculator;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.SimpleFaultData;

import com.google.common.base.Preconditions;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.inversion.SectionCluster;
import scratch.UCERF3.inversion.SectionClusterList;
import scratch.UCERF3.inversion.SectionConnectionStrategy;
import scratch.UCERF3.utils.DeformationModelFetcher;
import scratch.UCERF3.utils.FaultSystemIO;

public class DownDipTestRupSetBuilder {

	public static void main(String[] args) throws IOException {
		File outputFile = new File("/tmp/down_dip_sub_sect_rup_set.zip");
		
		// we're going to manually build faults here
		// first, build a big fault with down-dip subsections
		String sectName = "Test SubSect Down-Dip Fault";
		int sectID = 0;
		int startID = 0;
		double upperDepth = 0d;
		double lowerDepth = 30d;
		double dip = 20d;
//		int numDownDip = 5;
//		int numAlongStrike = 6;
		int numDownDip = 20;
		int numAlongStrike = 30;
		FaultTrace trace = new FaultTrace(sectName);
		trace.add(new Location(34, -119, upperDepth));
		trace.add(new Location(34.1, -118.75, upperDepth));
		trace.add(new Location(34.15, -118.5, upperDepth));
		trace.add(new Location(34.1, -118.25, upperDepth));
		trace.add(new Location(34, -118, upperDepth));
		
		SimpleFaultData faultData = new SimpleFaultData(dip, lowerDepth, upperDepth, trace);
		double aveRake = 90d;
		
		DownDipSubSectBuilder downDipBuilder = new DownDipSubSectBuilder(sectName, sectID, startID,
				faultData, aveRake, numAlongStrike, numDownDip);
		
		List<FaultSection> subSections = new ArrayList<>();
		subSections.addAll(downDipBuilder.getSubSectsList());
		System.out.println("Have "+subSections.size()+" sub-sections for "+sectName);
		startID = subSections.size();
		
		// now lets add a nearby crustal fault that it can interact with
		sectName = "Test Crustal Fault";
		sectID = 1;
		upperDepth = 0d;
		lowerDepth = 14d;
		dip = 90d;
		trace = new FaultTrace(sectName);
		trace.add(new Location(33.5, -117.5, upperDepth));
		trace.add(new Location(34.1, -118.25, upperDepth));
		FaultSectionPrefData crustalSect = new FaultSectionPrefData();
		crustalSect.setSectionId(sectID);
		crustalSect.setSectionName(sectName);
		crustalSect.setFaultTrace(trace);
		crustalSect.setAveUpperDepth(upperDepth);
		crustalSect.setAveLowerDepth(lowerDepth);
		crustalSect.setAseismicSlipFactor(0d);
		crustalSect.setAveDip(dip);
		double maxSectLength = 0.5*crustalSect.getOrigDownDipWidth();
		
		subSections.addAll(crustalSect.getSubSectionsList(maxSectLength, startID, 2));
		System.out.println("Have "+subSections.size()+" sub-sections in total");
		
		for (int s=0; s<subSections.size(); s++)
			Preconditions.checkState(subSections.get(s).getSectionId() == s,
				"section at index %s has ID %s", s, subSections.get(s).getSectionId());
		
		// instantiate plausibility filters
		List<PlausibilityFilter> filters = new ArrayList<>();
		filters.add(new RectangularityFilter(downDipBuilder, 1));
		
		SectionDistanceAzimuthCalculator distAzCalc = new SectionDistanceAzimuthCalculator(subSections);
		
		// this creates rectangular permutations only for our down-dip fault to speed up rupture building
		ClusterPermutationStrategy permutationStrategy = new DownDipTestPermutationStrategy(downDipBuilder);
		// connection strategy: parent faults connect at closest point, and only when dist <=5 km
		ClusterConnectionStrategy connectionStrategy = new DistCutoffClosestSectClusterConnectionStrategy(5d);
		int maxNumSplays = 0; // don't allow any splays
		
		ClusterRuptureBuilder builder = new ClusterRuptureBuilder(subSections, connectionStrategy,
				distAzCalc, filters, maxNumSplays);
		
		List<ClusterRupture> ruptures = builder.build(permutationStrategy);
		
		System.out.println("Built "+ruptures.size()+" total ruptures");
		
		// build a rupture set (doing this manually instead of creating an inversion fault system rup set,
		// mostly as a demonstration)
		double[] sectSlipRates = new double[subSections.size()];
		double[] sectAreasReduced = new double[subSections.size()];
		double[] sectAreasOrig = new double[subSections.size()];
		for (int s=0; s<sectSlipRates.length; s++) {
			FaultSection sect = subSections.get(s);
			sectAreasReduced[s] = sect.getArea(true);
			sectAreasOrig[s] = sect.getArea(false);
			sectSlipRates[s] = sect.getReducedAveSlipRate()*1e-3; // mm/yr => m/yr
		}
		
		double[] rupMags = new double[ruptures.size()];
		double[] rupRakes = new double[ruptures.size()];
		double[] rupAreas = new double[ruptures.size()];
		double[] rupLengths = new double[ruptures.size()];
		ScalingRelationships scale = ScalingRelationships.SHAW_2009_MOD;
		List<List<Integer>> rupsIDsList = new ArrayList<>();
		for (int r=0; r<ruptures.size(); r++) {
			ClusterRupture rup = ruptures.get(r);
			List<FaultSection> rupSects = rup.buildOrderedSectionList();
			List<Integer> sectIDs = new ArrayList<>();
			double totLength = 0d;
			double totArea = 0d;
			double totOrigArea = 0d; // not reduced for aseismicity
			List<Double> sectAreas = new ArrayList<>();
			List<Double> sectRakes = new ArrayList<>();
			for (FaultSection sect : rupSects) {
				sectIDs.add(sect.getSectionId());
				double length = sect.getTraceLength()*1e3;	// km --> m
				totLength += length;
				double area = sectAreasReduced[sect.getSectionId()];	// sq-m
				totArea += area;
				totOrigArea += sectAreasOrig[sect.getSectionId()];	// sq-m
				sectAreas.add(area);
				sectRakes.add(sect.getAveRake());
			}
			rupAreas[r] = totArea;
			rupLengths[r] = totLength;
			rupRakes[r] = FaultUtils.getInRakeRange(FaultUtils.getScaledAngleAverage(sectAreas, sectRakes));
			double origDDW = totOrigArea/totLength;
			rupMags[r] = scale.getMag(totArea, origDDW);
			rupsIDsList.add(sectIDs);
		}
		
		String info = "Test down-dip subsectioning rup set";
		
		FaultSystemRupSet rupSet = new FaultSystemRupSet(subSections, sectSlipRates, null, sectAreasReduced,
				rupsIDsList, rupMags, rupRakes, rupAreas, rupLengths, info);
		FaultSystemIO.writeRupSet(rupSet, outputFile);
	}

}
