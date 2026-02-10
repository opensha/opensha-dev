package scratch.kevin.nshm23.bbpScaling;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.util.Precision;
import org.netlib.util.doubleW;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.FaultUtils;
import org.opensha.sha.earthquake.FocalMechanism;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.RectangularSurface;
import org.opensha.sha.faultSurface.RuptureSurface;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.BBP_SourceFile;
import scratch.kevin.bbp.BBP_SourceFile.BBP_PlanarSurface;

public class ExampleRupSourceFileWriter {

	public static void main(String[] args) throws IOException {
		List<int[]> parentCombs = new ArrayList<>();
		
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged.zip"));
		
		File mainOutputDir = new File("/home/kevin/CyberShake/nshm23_scaling");
		
		boolean creepReduced = true;
		
		List<? extends FaultSection> allSubSects = rupSet.getFaultSectionDataList();
		
		Map<Integer, List<FaultSection>> groupedByParent = allSubSects.stream()
				.collect(java.util.stream.Collectors.groupingBy(FaultSection::getParentSectionId));
		
		List<String> sectAbbrevs = new ArrayList<>();
		List<Integer> sectIDs = new ArrayList<>();
		
		double[] dists = {10d, 20d, 30d, 40d, 50d, 75d, 100d, 125, 150d, 175d, 200d};
		boolean[] rJBs = {true,false};
		
		FaultSection singleSect = null;
		
		File outputDir = new File(mainOutputDir, "southern_saf");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		sectAbbrevs.add("BB");
		sectIDs.add(FaultSectionUtils.findParentSectionID(allSubSects, "Andreas", "Big Bend"));
		singleSect = groupedByParent.get(sectIDs.get(0)).get(0);
		
		sectAbbrevs.add("MN");
		sectIDs.add(FaultSectionUtils.findParentSectionID(allSubSects, "Andreas", "Mojave", "north"));
		
		sectAbbrevs.add("MS");
		sectIDs.add(FaultSectionUtils.findParentSectionID(allSubSects, "Andreas", "Mojave", "south"));
		
		sectAbbrevs.add("SBN");
		sectIDs.add(FaultSectionUtils.findParentSectionID(allSubSects, "Andreas", "Bernardino", "north"));
		
		sectAbbrevs.add("NB");
		sectIDs.add(FaultSectionUtils.findParentSectionID(allSubSects, "Andreas", "North", "Branch"));
		
//		sectAbbrevs.add("SBS");
//		sectIDs.add(FaultSectionUtils.findParentSectionID(allSubSects, "Andreas", "Bernardino", "south"));
		
//		sectAbbrevs.add("SGP");
//		sectIDs.add(FaultSectionUtils.findParentSectionID(allSubSects, "Andreas", "Gorgonio", "Garnet"));
		
		sectAbbrevs.add("CC");
		sectIDs.add(FaultSectionUtils.findParentSectionID(allSubSects, "Andreas", "Coachella"));
		
//		File outputDir = new File(mainOutputDir, "northridge");
//		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
//		
//		sectAbbrevs.add("NR");
//		sectIDs.add(189);
//		singleSect = groupedByParent.get(sectIDs.get(0)).get(0);
		
//		File outputDir = new File(mainOutputDir, "compton");
//		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
//		
//		sectAbbrevs.add("CMPT");
//		sectIDs.add(FaultSectionUtils.findParentSectionID(allSubSects, "Compton"));
//		singleSect = groupedByParent.get(sectIDs.get(0)).get(0);
		
		
		
		
		for (int i=0; i<sectAbbrevs.size(); i++) {
			String sectName = sectAbbrevs.get(i);
			int sectID = sectIDs.get(i);
			String fullName = groupedByParent.get(sectID).get(0).getParentSectionName();
			System.out.println(sectName+" is "+sectID+": "+fullName);
		}
		
		List<Integer> rupsToSearch = rupSet.getRupturesForParentSection(sectIDs.get(0));
		
		List<String> rupNames = new ArrayList<>();
		List<Integer> rupIndexes = new ArrayList<>();
		
		Set<Integer> curParentIDs = new HashSet<>();
		String curName = null;
		for (int i=0; i<sectAbbrevs.size(); i++) {
			String sectName = sectAbbrevs.get(i);
			int sectID = sectIDs.get(i);
			Preconditions.checkState(sectID >= 0);
			
			if (curName == null)
				curName = "";
			else
				curName += "_";
			
			curParentIDs.add(sectID);
			curName += sectName;
			
			double maxLen = 0d;
			int idForMaxLen = -1;
			
			for (int rupIndex : rupsToSearch) {
				List<FaultSection> sectsForRup = rupSet.getFaultSectionDataForRupture(rupIndex);
				boolean allMatch = true;
				Set<Integer> rupParentIDs = new HashSet<>(curParentIDs.size());
				for (FaultSection sect : sectsForRup) {
					if (!curParentIDs.contains(sect.getParentSectionId())) {
						allMatch = false;
						break;
					}
					rupParentIDs.add(sect.getParentSectionId());
				}
				if (allMatch && rupParentIDs.size() == curParentIDs.size()) {
					double length = rupSet.getLengthForRup(rupIndex) * 1e-3;
					if (length > maxLen) {
						maxLen = length;
						idForMaxLen = rupIndex;
					}
				}
			}
			
			Preconditions.checkState(idForMaxLen >= 0,
					"Didn't find any ruptures that use each of (and nothing else): %s\n\t%s", curName, curParentIDs);
			
			double avgMag = rupSet.getMagForRup(idForMaxLen);
			System.out.println("Longest rup found for "+curName+": "+(float)maxLen+" km, M"+(float)avgMag);
			
			if (singleSect != null && singleSect.getParentSectionId() == sectID) {
				// also add single fault ruptures
				for (int rupIndex : rupSet.getRupturesForSection(singleSect.getSectionId())) {
					if (rupIndex == idForMaxLen)
						continue;
					List<FaultSection> sectsForRup = rupSet.getFaultSectionDataForRupture(rupIndex);
					boolean allMatch = true;
					for (FaultSection sect : sectsForRup) {
						if (sect.getParentSectionId() != sectID) {
							allMatch = false;
							break;
						}
					}
					if (allMatch) {
						String rupName = curName+"_partial_"+sectsForRup.size()+"sects";
						double length = rupSet.getLengthForRup(rupIndex) * 1e-3;
						System.out.println("Also adding single-fault "+rupName+": "
								+(float)length+" km, M"+(float)rupSet.getMagForRup(rupIndex));
						rupNames.add(rupName);
						rupIndexes.add(rupIndex);
					}
				}
			}
			
			rupNames.add(curName);
			rupIndexes.add(idForMaxLen);
		}
		
		int fixedSeed = 123456789;
		double dWid = 0.1;
		double dLen = 0.1;
		double cornerFreq = 0.15;
		
		NSHM23_ScalingRelationships[] scales = {
				NSHM23_ScalingRelationships.AVERAGE,
				NSHM23_ScalingRelationships.LOGA_C4p1,
				NSHM23_ScalingRelationships.LOGA_C4p2,
				NSHM23_ScalingRelationships.LOGA_C4p3,
				NSHM23_ScalingRelationships.WIDTH_LIMITED,
				NSHM23_ScalingRelationships.WIDTH_LIMITED_CSD,
				NSHM23_ScalingRelationships.LOGA_C4p2_SQRT_LEN
		};
		
		VelocityModel vm = VelocityModel.LA_BASIN_500;
		
		for (int r=0; r<rupIndexes.size(); r++) {
			int rupIndex = rupIndexes.get(r);
			String name = rupNames.get(r);
			
			BBP_PlanarSurface bbpSurf = planarEquivalentSurface(rupSet, rupIndex, creepReduced);
			
			double hypoAlong = 0.5*bbpSurf.getLength();
			double hypoDown = 0.5*bbpSurf.getWidth();
			
			// these are in SI units, not km
			double area = rupSet.getAreaForRup(rupIndex);
			double length = rupSet.getLengthForRup(rupIndex);
			double width = area/length;
			double totOrigArea = 0;
			for (FaultSection sect : rupSet.getFaultSectionDataForRupture(rupIndex)) {
				totOrigArea += sect.getArea(false);	// sq-m
			}
			double origDDW = totOrigArea / length;
			System.out.println(name+" ("+rupIndex+"):\tA="+(float)area+"; L="+(float)length+"; W="+(float)width+"; origW="+(float)origDDW);
			double aveRake = rupSet.getAveRakeForRup(rupIndex);
//			mags[r] = scale.getMag(rupAreas[r], rupLengths[r], rupAreas[r]/rupLengths[r], origDDW, rakes[r]);
			
			CSVFile<String> csv = new CSVFile<>(true);
			csv.addLine("Scaling Relationship", "NSHM23 Magnitude", "NSHM23 Average Slip (m)");
			
			for (NSHM23_ScalingRelationships scale : scales) {
				double mag = scale.getMag(area, length, width, origDDW, aveRake);
				double slip = scale.getAveSlip(area, length, width, origDDW, aveRake);
				csv.addLine(scale.getShortName(), (float)mag+"", (float)slip+"");
			}
			
			csv.writeToFile(new File(outputDir, name+"_scaling.csv"));
			
			BBP_SourceFile src = new BBP_SourceFile(bbpSurf, rupSet.getMagForRup(rupIndex), hypoAlong, hypoDown, dWid, dLen, cornerFreq, fixedSeed);
			src.writeToFile(new File(outputDir, name+".src"));

			RectangularSurface rectSurf = bbpSurf.getRectSurface();
			for (boolean rJB : rJBs) {
				List<LocationList> racetracks = new ArrayList<>();
				MinMaxAveTracker latTrack = new MinMaxAveTracker();
				MinMaxAveTracker lonTrack = new MinMaxAveTracker();
				String trackPrefix = name+"_racetrack";
				if (rJB)
					trackPrefix += "_rJB";
				else
					trackPrefix += "_rRup";
				for (double dist : dists) {
					LocationList racetrack = buildDiscrRacetrack(rectSurf, dist, 50, rJB);
					CSVFile<String> trackCSV = new CSVFile<>(true);
					trackCSV.addLine("Site Index", "Latitude", "Longitude", "Rjb", "Rrup", "RX");
					List<BBP_Site> sites = new ArrayList<>();
					for (int i=0; i<racetrack.size(); i++) {
						Location loc = racetrack.get(i);
						latTrack.addValue(loc.lat);
						lonTrack.addValue(loc.lon);
						
						trackCSV.addLine(i+"", (float)loc.lat+"", (float)loc.lon+"",
								(float)rectSurf.getDistanceJB(loc)+"",
								(float)rectSurf.getDistanceRup(loc)+"",
								(float)rectSurf.getDistanceX(loc)+"");
						sites.add(new BBP_Site("s"+i, loc, vm.getVs30(), 0.15, 100d));
					}
					String prefix = trackPrefix+"_"+(int)dist+"km";
					BBP_Site.writeToFile(new File(outputDir, prefix+".stl"), sites);
					trackCSV.writeToFile(new File(outputDir, prefix+".csv"));
					racetracks.add(racetrack);
				}
				GeographicMapMaker mapMaker = new GeographicMapMaker(new Region(
						new Location(latTrack.getMin()-0.5, lonTrack.getMin()-0.5),
						new Location(latTrack.getMax()+0.5, lonTrack.getMax()+0.5)));
				mapMaker.setWritePDFs(false);
				mapMaker.setWriteGeoJSON(false);
				List<LocationList> lines = new ArrayList<>();
				List<PlotCurveCharacterstics> lineChars = new ArrayList<>();
				
				lines.add(rectSurf.getPerimeter());
				lineChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.DARK_GRAY));
				lines.add(rectSurf.getUpperEdge());
				lineChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
				
				LocationList scatters = new LocationList();
				for (LocationList racetrack : racetracks) {
					LocationList closed = new LocationList(racetracks.size()+1);
					closed.addAll(racetrack);
					closed.add(racetrack.first());
					lines.add(closed);
					lineChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Colors.tab_blue));
					
					scatters.addAll(racetrack);
				}
				
				mapMaker.plotLines(lines, lineChars);
				
				mapMaker.plotScatters(scatters, Colors.tab_orange);
				mapMaker.setScatterSymbol(PlotSymbol.FILLED_CIRCLE, 1f);
				
				mapMaker.plot(outputDir, trackPrefix, " ");
			}
		}
	}
	
	public static BBP_PlanarSurface planarEquivalentSurface(FaultSystemRupSet rupSet, int rupIndex, boolean creepReduced) {
		RuptureSurface surf = rupSet.getSurfaceForRupture(rupIndex, 1d, creepReduced);
		
		Location firstLoc = surf.getFirstLocOnUpperEdge();
		Location lastLoc = surf.getLastLocOnUpperEdge();
		
		LocationVector vector = LocationUtils.vector(firstLoc, lastLoc);
		double strike = vector.getAzimuth();
		double rake = rupSet.getAveRakeForRup(rupIndex);
		double dip = surf.getAveDip();
		
		FocalMechanism mech = new FocalMechanism(strike, dip, rake);
		
		double vectorLen = vector.getHorzDistance();
		
		LocationVector halfVector = new LocationVector(vector.getAzimuth(), 0.5*vectorLen, 0d);
		Location topCenter = LocationUtils.location(firstLoc, halfVector);
		topCenter = new Location(topCenter.getLatitude(), topCenter.getLongitude(), surf.getAveRupTopDepth());
		
		double length = surf.getAveLength();
		double width = surf.getAveWidth();
		
		return new BBP_PlanarSurface(topCenter, length, width, mech);
	}
	
	private static final double[] QUARTER_ARC_AZIMUTHS;
	
	static {
		QUARTER_ARC_AZIMUTHS = new double[100];
		double deltaEach = 0.5*Math.PI/QUARTER_ARC_AZIMUTHS.length;
		for (int i=0; i<QUARTER_ARC_AZIMUTHS.length; i++)
			QUARTER_ARC_AZIMUTHS[i] = (0.5+i)*deltaEach;
	}
	
	private static LocationList buildRacetrack(RectangularSurface surf, double distance, boolean rJB) {
		// will evenly discretize later
		LocationList arbDiscrRacetrack = new LocationList();
		
		// these might have Z>0
		Location firstUpper3D = surf.getFirstLocOnUpperEdge();
		Location lastUpper3D = surf.getLastLocOnUpperEdge();
		Preconditions.checkState(Precision.equals(firstUpper3D.depth, lastUpper3D.depth, 1e-4));
		Location firstUpper = new Location(firstUpper3D.lat, firstUpper3D.lon);
		Location lastUpper = new Location(lastUpper3D.lat, lastUpper3D.lon);
		// degrees
		double strike = surf.getAveStrike();
		double strikeRad = Math.toRadians(strike);
		double leftAz = strikeRad - 0.5*Math.PI;
		double rightAz = strikeRad + 0.5*Math.PI;
		double backAz = strikeRad + Math.PI;
		
		double surfDistance;
		if (rJB || Precision.equals(firstUpper3D.depth, 0d, 1e-4))
			surfDistance = distance;
		else
			surfDistance = Math.sqrt(distance*distance - firstUpper3D.depth*firstUpper3D.depth);
		
		// always a straight line to the left
		arbDiscrRacetrack.add(LocationUtils.location(firstUpper, leftAz, surfDistance));
		arbDiscrRacetrack.add(LocationUtils.location(lastUpper, leftAz, surfDistance));
		
		// now a clockwise quarter circle (this leaves off the first and last point of the semi circle
		for (double deltaAz : QUARTER_ARC_AZIMUTHS)
			arbDiscrRacetrack.add(LocationUtils.location(lastUpper, leftAz+deltaAz, surfDistance));
		
		// now straight ahead
		arbDiscrRacetrack.add(LocationUtils.location(lastUpper, strikeRad, surfDistance));
		if (surf.getAveDip() <89.999) {
			// dipping, complicated
			// first -> last still in the along strike direction
			// these have Z>0
			Location firstLower3D = surf.getFirstLocOnLowerEdge();
			Location lastLower3D = surf.getLastLocOnLowerEdge();
			Location firstLower = new Location(firstLower3D.lat, firstLower3D.lon);
			Location lastLower = new Location(lastLower3D.lat, lastLower3D.lon);
			if (rJB) {
				// we need to move about the surface projection of the box
				// now along strike from bottom
				arbDiscrRacetrack.add(LocationUtils.location(lastLower, strikeRad, surfDistance));
				
				// now arc
				for (double deltaAz : QUARTER_ARC_AZIMUTHS)
					arbDiscrRacetrack.add(LocationUtils.location(lastLower, strikeRad+deltaAz, surfDistance));
				
				// now right side
				arbDiscrRacetrack.add(LocationUtils.location(lastLower, rightAz, surfDistance));
				arbDiscrRacetrack.add(LocationUtils.location(firstLower, rightAz, surfDistance));
				
				// now arc
				for (double deltaAz : QUARTER_ARC_AZIMUTHS)
					arbDiscrRacetrack.add(LocationUtils.location(firstLower, rightAz+deltaAz, surfDistance));
				
				// now bottom side
				arbDiscrRacetrack.add(LocationUtils.location(firstLower, backAz, surfDistance));
			} else {
				// rRup, complicated
				double rRupTol = 1e-7;
				
				// arc to the right (this will be a wide arc)
				for (double deltaAz : QUARTER_ARC_AZIMUTHS)
					arbDiscrRacetrack.add(findRrupLoc(surf, distance, lastUpper, strikeRad+deltaAz, surfDistance, rRupTol));
				// right side
				arbDiscrRacetrack.add(findRrupLoc(surf, distance, lastUpper, rightAz, surfDistance, rRupTol));
				arbDiscrRacetrack.add(findRrupLoc(surf, distance, firstUpper, rightAz, surfDistance, rRupTol));
				// arc to the bottom (wide)
				for (double deltaAz : QUARTER_ARC_AZIMUTHS)
					arbDiscrRacetrack.add(findRrupLoc(surf, distance, firstUpper, rightAz+deltaAz, surfDistance, rRupTol));
			}
		} else {
			// vertical, simple
			
			// another simple quarter circle
			for (double deltaAz : QUARTER_ARC_AZIMUTHS)
				arbDiscrRacetrack.add(LocationUtils.location(lastUpper, strikeRad+deltaAz, surfDistance));
			
			// to the right
			arbDiscrRacetrack.add(LocationUtils.location(lastUpper, rightAz, surfDistance));
			arbDiscrRacetrack.add(LocationUtils.location(firstUpper, rightAz, surfDistance));
			
			// half curcle below
			for (double deltaAz : QUARTER_ARC_AZIMUTHS)
				arbDiscrRacetrack.add(LocationUtils.location(firstUpper, rightAz+deltaAz, surfDistance));
			
		}
		// all cases end with an arc from the back azimuth
		arbDiscrRacetrack.add(LocationUtils.location(firstUpper, backAz, surfDistance));
		for (double deltaAz : QUARTER_ARC_AZIMUTHS)
			arbDiscrRacetrack.add(LocationUtils.location(firstUpper, backAz+deltaAz, surfDistance));
		
		// and close it
		arbDiscrRacetrack.add(arbDiscrRacetrack.first());
		
		return arbDiscrRacetrack;
	}
	
	private static LocationList buildDiscrRacetrack(RectangularSurface surf, double distance, int num, boolean rJB) {
		// will evenly discretize later
		LocationList arbDiscrRacetrack = buildRacetrack(surf, distance, rJB);
		
		// now discretize it
		// this method takes a fault trace, convert it
		FaultTrace trace = new FaultTrace(null, arbDiscrRacetrack.size());
		trace.addAll(arbDiscrRacetrack);
		// this returns a trace of size num+1, which is what we want; we'll throw out the last (duplicate) location
		FaultTrace resampled = FaultUtils.resampleTrace(trace, num);
		Preconditions.checkState(resampled.size() == num+1);
		Preconditions.checkState(LocationUtils.areSimilar(resampled.first(), resampled.last()));
		LocationList ret = new LocationList(num);
		for (int i=0; i<num; i++)
			ret.add(resampled.get(i));
		return ret;
	}
	
	private static Location findRrupLoc(RectangularSurface surf, double rRup, Location startLoc, double azimuthRad, double guessDistance, double tol) {
		int iters = 0;
		double distance = guessDistance;
		LocationVector vect = new LocationVector(Math.toDegrees(azimuthRad), distance, 0d);
		Location loc = LocationUtils.location(startLoc, vect);
		double calcRrup = surf.getDistanceRup(loc);
		while (iters < 100 || !Precision.equals(rRup, calcRrup, tol)) {
//			System.out.println("Iter "+iters+"; calcRrup="+calcRrup+", distance="+distance);
			double delta = calcRrup-rRup;
			distance -= delta;
			vect.setHorzDistance(distance);
			loc = LocationUtils.location(startLoc, vect);
			calcRrup = surf.getDistanceRup(loc);
			iters++;
			if ((float)rRup == (float)calcRrup && Precision.equals(rRup, calcRrup, 1e-10))
				// drop iter requirement, we've nailed it
				break;
		}
//		System.exit(0);
		return loc;
	}

}
