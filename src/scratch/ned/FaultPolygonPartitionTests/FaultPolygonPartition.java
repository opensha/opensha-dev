package scratch.ned.FaultPolygonPartitionTests;

import java.awt.Color;


import java.util.ArrayList;

import org.jfree.data.Range;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import scratch.ned.FSS_Inversion2019.PlottingUtils;

public class FaultPolygonPartition {
	
	
	public static void doit(double polygonHalfWidthKm, boolean normalizeRatesByArea) {
		
		double totMgt5_Rate = 10;
		double faultMinSupraMag = 6.55;
		double maxOffFaultMag = 7.55;
		
		double minLat = 33;
		double maxLat = 34;
		double midLat = (minLat+maxLat)/2.0;

		double minLon = -125.0;
		double maxLon = -124.0;
		double midLon = (minLon+maxLon)/2.0;
		
		Location midLoc = new Location(midLat,midLon);
		
		Location tempLoc = LocationUtils.location(midLoc, (90.0*Math.PI/180.), polygonHalfWidthKm);
		double testDistance = LocationUtils.horzDistance(midLoc,tempLoc);
//		System.out.println("testHalfWidth="+(float)testDistance);
//		System.out.println("midLoc="+midLoc);
//		System.out.println("tempLoc="+tempLoc);

		double deltaLon = tempLoc.getLongitude()-midLon;
//		System.out.println("deltaLon="+(float)deltaLon);

		double minLonPolygon = midLon-deltaLon;
		double maxLonPolygon = midLon+deltaLon;
		
		Region region = new Region(new Location(minLat,minLon),new Location(maxLat,maxLon));
		
		Region faultPolygonRegion = new Region(new Location(minLat,minLonPolygon),new Location(maxLat,maxLonPolygon));
		
		
		double regionWidth = LocationUtils.horzDistance(new Location(midLat,minLon),new Location(midLat,maxLon));
		
		double polyRegionWidth = LocationUtils.horzDistance(new Location(midLat,minLonPolygon),new Location(midLat,maxLonPolygon));
		
		double polyFractArea = polyRegionWidth/regionWidth;

		
		System.out.println("polygonHalfWidthKm="+(float)polygonHalfWidthKm);
		System.out.println("regionHalfWidth="+(float)(regionWidth/2.0));
		System.out.println("polygonHalfWidthTest="+(float)(polyRegionWidth/2.0));

		
		if(polygonHalfWidthKm>regionWidth/2.0)
			throw new RuntimeException("polygonHalfWidthKm too large");
		
		FaultTrace faultTrace = new FaultTrace("parentFaultTrace");
		faultTrace.add(new Location(minLat,midLon));
		faultTrace.add(new Location(maxLat,midLon));
		
		double traceLength = faultTrace.getTraceLength();
		
		System.out.println("traceLength="+(float)traceLength);
		
		double regionArea = regionWidth*traceLength;
		double polyArea = polyRegionWidth*traceLength;
		double offFaultArea = regionArea-polyArea;


		
		GutenbergRichterMagFreqDist totMFD = new GutenbergRichterMagFreqDist(5.05,30,0.1);
		totMFD.setAllButTotMoRate(5.05, 7.95, totMgt5_Rate, 1.0);
		totMFD.setName("totMFD");
		
		double totOnFltRate = totMgt5_Rate*polygonHalfWidthKm/(regionWidth/2.0);
		
		GutenbergRichterMagFreqDist subSeisOnFaultMFD = new GutenbergRichterMagFreqDist(5.05,30,0.1);
		subSeisOnFaultMFD.setAllButTotMoRate(5.05, faultMinSupraMag-0.1, totOnFltRate, 1.0);
		subSeisOnFaultMFD.setName("subSeisOnFaultMFD");

		
		IncrementalMagFreqDist totOnFaultMFD = new IncrementalMagFreqDist(5.05,30,0.1);
		for(double mag = 5.05;mag<faultMinSupraMag-0.05; mag+=0.1)
			totOnFaultMFD.set(mag,subSeisOnFaultMFD.getY(mag));
		for(double mag = maxOffFaultMag+0.1;mag<8.0; mag+=0.1)
			totOnFaultMFD.set(mag,totMFD.getY(mag));
		
		double mag1 = faultMinSupraMag-0.1;
		double logRate1 = Math.log10(totOnFaultMFD.getY(mag1));
		double mag2 = maxOffFaultMag+0.1;
		double logRate2 = Math.log10(totOnFaultMFD.getY(mag2));
		for(double mag = faultMinSupraMag;mag<maxOffFaultMag+0.05; mag+=0.1) {
			double logRate = logRate1 - ((mag-mag1)/(mag2-mag1))*(logRate1-logRate2);
			totOnFaultMFD.set(mag,Math.pow(10, logRate));
		}
		totOnFaultMFD.setName("totOnFaultMFD");
		
		
		IncrementalMagFreqDist totOffFaultMFD = new IncrementalMagFreqDist(5.05,30,0.1);
		for(int i=0;i<totOffFaultMFD.size();i++)
			totOffFaultMFD.set(i,totMFD.getY(i)-totOnFaultMFD.getY(i));
		totOffFaultMFD.setName("totOffFaultMFD");
		
		// Normalize MFDs by region area
		String namePrefix="";
		if(normalizeRatesByArea) {
			totMFD.scale(1.0/regionArea);
			subSeisOnFaultMFD.scale(1.0/polyArea);
			totOnFaultMFD.scale(1.0/polyArea); // don't normalize because it's all put on fault
			totOffFaultMFD.scale(1.0/offFaultArea);		
			namePrefix = "Area Normalized; ";
		}
		
		System.out.println("regionArea="+regionArea);
		System.out.println("polyArea="+polyArea);
		System.out.println("offFaultArea="+offFaultArea);
		System.out.println("offFaultArea+polyArea="+(offFaultArea+polyArea));

		
		ArrayList<XY_DataSet> mfdList = new ArrayList<XY_DataSet>();
		mfdList.add(totMFD);
		mfdList.add(subSeisOnFaultMFD);
		mfdList.add(totOnFaultMFD);
		mfdList.add(totOffFaultMFD);
	
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();	
		
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GREEN));
		
		Range xAxisRange = null;
		Range yAxisRange = new Range(1e-9,10);
		boolean logX = false;
		boolean logY = true;

		PlottingUtils.writeAndOrPlotFuncs(mfdList, plotChars, namePrefix+" Fract Area: "+(float)polyFractArea, "Mag", "Rate", 
				xAxisRange, yAxisRange, logX, logY, null, true);




		
	}

	public static void main(String[] args) {
//		doit(12, true);
//		doit(1, true);
//		doit(25, true);
		doit(46, true);
		
//		doit(0.001, false);
//		doit(12, false);
//		doit(1, false);
//		doit(25, false);
		doit(46, false);

	}

}
