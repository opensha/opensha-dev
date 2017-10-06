package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.commons.util.FaultUtils;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.FocalMechanism;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.srf.RSQSimEventSlipTimeFunc;
import org.opensha.sha.simulators.srf.RSQSimSRFGenerator;
import org.opensha.sha.simulators.srf.RSQSimSRFGenerator.SRFInterpolationMode;
import org.opensha.sha.simulators.srf.SRF_PointData;
import org.opensha.sha.simulators.utils.RSQSimUtils;
import org.opensha.sha.simulators.utils.RupturePlotGenerator;
import org.opensha.sha.simulators.utils.SimulatorUtils;

import com.google.common.base.Preconditions;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.kevin.bbp.BBP_Module.Method;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_SourceFile;
import scratch.kevin.bbp.BBP_SourceFile.BBP_PlanarSurface;
import scratch.kevin.bbp.BBP_Wrapper;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class RSQSimBBP_Config {
	
	private static BBP_PlanarSurface estimateBBP_PlanarSurface(SimulatorEvent event) {
		double length = SimulatorUtils.estimateRuptureLength(event);
		double zTOR = Double.POSITIVE_INFINITY;
		double zBOR = Double.NEGATIVE_INFINITY;
		ArrayList<SimulatorElement> elems = event.getAllElements();
		for (SimulatorElement e : elems) {
			for (Location loc : e.getVertices()) {
				double z = loc.getDepth();
				zTOR = Math.min(zTOR, z);
				zBOR = Math.max(zBOR, z);
			}
		}
		if (zTOR < 0)
			zTOR = 0;
		System.out.println("zTOR = "+zTOR);
		System.out.println("zBOR = "+zBOR);
		List<Double> strikes = new ArrayList<>();
		List<Double> rakes = new ArrayList<>();
		List<Double> dips = new ArrayList<>();
		for (SimulatorElement e : elems) {
			FocalMechanism mech = e.getFocalMechanism();
			double strike = mech.getStrike();
			double rake = mech.getRake();
			double dip = mech.getDip();
			if (dip == 90 && (rake == 180 || rake == -180)) {
				
			}
			strikes.add(strike);
			rakes.add(rake);
			dips.add(dip);
		}
		double aveRake = FaultUtils.getInRakeRange(FaultUtils.getAngleAverage(rakes));
		double aveStrike = FaultUtils.getAngleAverage(strikes);
		double aveDip = 0d;
		for (double dip : dips)
			aveDip += dip;
		aveDip /= dips.size();
		FocalMechanism mech = new FocalMechanism(aveStrike, aveDip, aveRake);
		
		double width = (zBOR - zTOR)/Math.sin(Math.toRadians(aveDip));
		
		double aveLat = 0d;
		double aveLon = 0d;
		for (SimulatorElement elem : elems) {
			Location loc = elem.getCenterLocation();
			aveLat += loc.getLatitude();
			aveLon += loc.getLongitude();
		}
		aveLat /= elems.size();
		aveLon /= elems.size();
		
		Location topCenter = new Location(aveLat, aveLon, zTOR);
		
		return new BBP_PlanarSurface(topCenter, length, width, mech);
	}
	
	private static BBP_PlanarSurface planarEquivalentU3Surface(RSQSimCatalog catalog, RSQSimEvent event,
			double minFractForInclusion, boolean adjWidthMatchArea) {
		EqkRupture rup = catalog.getGMPE_Rupture(event, minFractForInclusion);
		RuptureSurface u3Surf = rup.getRuptureSurface();
		
		Location firstLoc = u3Surf.getFirstLocOnUpperEdge();
		Location lastLoc = u3Surf.getLastLocOnUpperEdge();
		
		LocationVector vector = LocationUtils.vector(firstLoc, lastLoc);
		double strike = vector.getAzimuth();
		double rake = rup.getAveRake();
		double dip = u3Surf.getAveDip();
		
		FocalMechanism mech = new FocalMechanism(strike, dip, rake);
		
		double length = vector.getHorzDistance();
		
		LocationVector halfVector = new LocationVector(vector.getAzimuth(), 0.5*length, 0d);
		Location topCenter = LocationUtils.location(firstLoc, halfVector);
		topCenter = new Location(topCenter.getLatitude(), topCenter.getLongitude(), u3Surf.getAveRupTopDepth());
		
		double width = u3Surf.getAveWidth();
		
		if (adjWidthMatchArea) {
			double eventArea = event.getArea()*1e-6; // to km^2
			width = eventArea/length;
		}
		
		return new BBP_PlanarSurface(topCenter, length, width, mech);
	}
	
	private static double absAngleDiff(double angle1, double angle2) {
		double angleDiff = Math.abs(angle1 - angle2);
		while (angleDiff > 270)
			angleDiff -= 360;
		angleDiff = Math.abs(angleDiff);
		return angleDiff;
	}
	
	private static BBP_SourceFile buildBBP_Source(SimulatorEvent event, RSQSimEventSlipTimeFunc func,
			BBP_PlanarSurface surface, int seed) {
		ArrayList<SimulatorElement> elems = event.getAllElements();
		
		Location hypo = null;
		double earliestTime = Double.POSITIVE_INFINITY;
		for (SimulatorElement elem : elems) {
			double time = func.getTimeOfFirstSlip(elem.getID());
			if (time < earliestTime) {
				earliestTime = time;
				hypo = elem.getCenterLocation();
			}
		}
		LocationVector hypoVector = LocationUtils.vector(surface.getTopCenter(), hypo);
		double hypoAlongStrike = hypoVector.getHorzDistance();
		if (hypoAlongStrike > surface.getLength())
			hypoAlongStrike = surface.getLength();
		double angleDiff = absAngleDiff(hypoVector.getAzimuth(), surface.getFocalMechanism().getStrike());
		if (angleDiff > 90)
			// it's reverse
			hypoAlongStrike = -hypoAlongStrike;
		double hypoDownDip = Math.abs(hypoVector.getVertDistance());
		if (hypoDownDip > surface.getWidth())
			hypoDownDip = surface.getWidth();
		
		double meanArea = 0d;
		for (SimulatorElement e : elems)
			meanArea += e.getArea()/1000000d; // to km^2
		meanArea /= elems.size();
		double sqSideLen = Math.sqrt(meanArea);
		
		System.out.println("Mean area: "+meanArea+", dx: "+sqSideLen);
		
		return new BBP_SourceFile(surface, event.getMagnitude(), hypoAlongStrike, hypoDownDip, sqSideLen, sqSideLen, 0.15, seed);
//		return new BBP_SourceFile(event.getMagnitude(), topCenter, length, width, mech, hypoAlongStrike, hypoDownDip,
//				sqSideLen, sqSideLen, 0.15, seed);
	}

	public static void main(String[] args) throws IOException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		
		RSQSimCatalog catalog = Catalogs.BRUCE_2194_LONG.instance(baseDir);
//		int[] eventIDs = { 399681 };
		int[] eventIDs = { 136704, 145982 };
		
//		RSQSimCatalog catalog = Catalogs.JG_UCERF3_millionElement.instance(baseDir);
//		int[] eventIDs = { 4099020 };
		
		File bbpSiteFile = new File("/home/kevin/bbp/bbp_data/run/stations_cs_sites.stl");
		
		boolean runBBP = false;
		boolean u3Surfaces = true;
		
		VelocityModel vm = VelocityModel.LA_BASIN;
		Method method = Method.GP;
		
		File outputDir = new File(catalog.getCatalogDir(), "event_srfs");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		System.out.println("Loading geometry...");
		List<SimulatorElement> elements = catalog.getElements();
		double meanArea = 0d;
		for (SimulatorElement e : elements)
			meanArea += e.getArea()/1000000d; // to km^2
		meanArea /= elements.size();
		System.out.println("Loaded "+elements.size()+" elements. Mean area: "+(float)meanArea+" km^2");
		List<FaultSectionPrefData> subSects = catalog.getU3SubSects();
		RSQSimUtils.cleanVertFocalMechs(elements, subSects);
		System.out.println("Loading events...");
		List<RSQSimEvent> events = catalog.loadEventsByID(eventIDs);
		System.out.println("Loaded "+events.size()+" events");
		Preconditions.checkState(events.size() == eventIDs.length);
		
		SRFInterpolationMode mode = SRFInterpolationMode.ADJ_VEL;
		double dt = 0.05;
		double srfVersion = 1.0;
		
		int seed = 12345;
		
		for (RSQSimEvent event : events) {
			System.out.println("Event: "+event.getID()+", M"+(float)event.getMagnitude());
			int eventID = event.getID();
			RSQSimEventSlipTimeFunc func = catalog.getSlipTimeFunc(event);
			String eventStr = "event_"+eventID;
			
			BBP_PlanarSurface surface;
			if (u3Surfaces)
				surface = planarEquivalentU3Surface(catalog, event, 0.2, true);
			else
				surface = estimateBBP_PlanarSurface(event);
			BBP_SourceFile bbpSource = buildBBP_Source(event, func, surface, seed);
			File srcFile = new File(outputDir, eventStr+".src");
			bbpSource.writeToFile(srcFile);
			
			SimulatorUtils.estimateVertexDAS(event);
			Location[] sourceRect = bbpSource.buildRectangle();
			Location sourceHypo = bbpSource.getHypoLoc();
			RupturePlotGenerator.writeSlipPlot(event, func, outputDir, eventStr, sourceRect, sourceHypo, null);
			RupturePlotGenerator.writeMapPlot(elements, event, func, outputDir, eventStr+"_map", sourceRect, sourceHypo, null);
			
			File eventOutputDir = new File(outputDir, eventStr+"_"+(float)dt+"s");
			System.out.println("dt="+dt+" => "+outputDir.getAbsolutePath());
			Preconditions.checkState(eventOutputDir.exists() || eventOutputDir.mkdir());
			
			ArrayList<SimulatorElement> patches = event.getAllElements();
			
			File srfFile = new File(eventOutputDir.getAbsolutePath()+"_"+mode.name()+".srf");
			System.out.println("Generating SRF for dt="+(float)dt+", "+mode);
			List<SRF_PointData> srf = RSQSimSRFGenerator.buildSRF(func, patches, dt, mode);
			SRF_PointData.writeSRF(srfFile, srf, srfVersion);
			
			if (runBBP) {
				File bbpOutputDir = new File(eventOutputDir.getAbsolutePath()+"_"+mode.name()+"_bbp");
				BBP_Wrapper bbpWrap = new BBP_Wrapper(vm, method, srcFile, null, srfFile, bbpSiteFile, bbpOutputDir);
				bbpWrap.setDoHF(false);
				bbpWrap.run();
			}
		}
	}

}
