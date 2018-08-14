package scratch.aftershockStatisticsETAS;

import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.faultSurface.FaultTrace;

import com.google.common.base.Stopwatch;

import wContour.Contour;
import wContour.Global.Border;
import wContour.Global.PolyLine;

public class ETAS_RateModel2D {

	private ETAS_AftershockModel forecastModel;
	private GriddedGeoDataSet rateModel;
	private Boolean D = false;
	
	public ETAS_RateModel2D(ETAS_AftershockModel forecastModel){
		this.forecastModel = forecastModel;
	}
	
	public ETAS_RateModel2D(ETAS_AftershockModel forecastModel, double plotDuration, double spacing, double stressDrop, double mainshockFitDuration, String fitType, FaultTrace faultTrace){
		this.forecastModel = forecastModel;
		this.rateModel = calculateRateModel( plotDuration, spacing, stressDrop,  mainshockFitDuration,  fitType,  faultTrace);
	}
	
	/**
	 * Returns a 2d grid of earthquake rates based on epicenters. 
	 * 
	 */
	public GriddedGeoDataSet calculateRateModel(double plotDuration, double spacing, double stressDrop, double mainshockFitDuration, String fitType, FaultTrace faultTrace){

		// TODO: assign these parameters somewhere, rather than hardcoded.
		double seismogenicDepth = 10;	//in km
		double mapSizeInRuptureLengths = 10;
		double minMapDistance = 20;	
		double forecastMaxDays = forecastModel.getForecastMinDays() + plotDuration;
		double maxDeltaMag = 3.0;	//minimum aftershock magnitude that contributes to rate estmate

		ObsEqkRupList aftershockPlotList = new ObsEqkRupList();
		ObsEqkRupList aftershockFitList = new ObsEqkRupList();
		ETASEqkRupture equivalentMainshock = new ETASEqkRupture(forecastModel.mainShock, stressDrop);
		
		aftershockPlotList = forecastModel.aftershockList.getRupsBefore((long) (forecastModel.mainShock.getOriginTime() + forecastModel.getForecastMinDays()* 24*60*60*1000));
		aftershockFitList = forecastModel.aftershockList.getRupsBefore((long) (forecastModel.mainShock.getOriginTime() + mainshockFitDuration*24*60*60*1000));

		// set up grid centered on mainshock
		
		Location centerLocation = ETAS_StatsCalc.getCentroid(forecastModel.mainShock, aftershockPlotList);
		double lat0 = centerLocation.getLatitude();
		double lon0 = centerLocation.getLongitude();
		double geomFactor = Math.cos(Math.toRadians(lat0)); //we use a single central estimate of degLon-->km mapping

		// make the readius set to the largest earthquake in the catalog
		double maxMag = forecastModel.mainShock.getMag();
		ObsEqkRupList largeAftershocks = forecastModel.aftershockList.getRupsAboveMag(maxMag);
		if (!largeAftershocks.isEmpty()) {
			for (ObsEqkRupture rup : aftershockPlotList) {
				if (rup.getMag() > maxMag)
					maxMag = rup.getMag();
			}
		}
		
		double mainshockRadius = ETAS_StatsCalc.magnitude2radius(maxMag, stressDrop);	//in km

		// make grid extend to n times source radius or at least 1 degree 
		double latmin, latmax, lonmin, lonmax;
		if (mainshockRadius*mapSizeInRuptureLengths < minMapDistance){
			latmin = lat0 - minMapDistance/111.111d;
			latmax = lat0 + minMapDistance/111.111d;
			lonmin = lon0 - minMapDistance/111.111/geomFactor;
			lonmax = lon0 + minMapDistance/111.111/geomFactor;
		}else{
			latmin = lat0 - mainshockRadius*mapSizeInRuptureLengths/111.111;
			latmax = lat0 + mainshockRadius*mapSizeInRuptureLengths/111.111;
			lonmin = lon0 - mainshockRadius*mapSizeInRuptureLengths/geomFactor/111.111;
			lonmax = lon0 + mainshockRadius*mapSizeInRuptureLengths/geomFactor/111.111;
		}
		
		GriddedRegion griddedRegion = new GriddedRegion(new Location(latmin, lonmin),
				new Location(latmax, lonmax), spacing, null);
		GriddedGeoDataSet gridData = new GriddedGeoDataSet(griddedRegion, false);

		if(D) System.out.println("Region:" + latmin + " " + latmax + " " + lonmin + " " + lonmax + " " + griddedRegion.getNodeCount() +" "+ geomFactor);

		if (fitType.equals("aftershocks") && aftershockFitList.size() >= 3){
			// fit finite mainshock source to early aftershocks
			if(D) System.out.println("Fitting " + aftershockFitList.size() + " early aftershocks, out of " + forecastModel.aftershockList.size() + " total aftershocks.");
			equivalentMainshock = ETAS_StatsCalc.fitMainshockLineSource(forecastModel.mainShock, aftershockFitList, stressDrop);
			
		} else if (fitType.equals("shakemap") && faultTrace != null && faultTrace.size() > 1){
			// fit finite mainshock source to shakemap source (fits line to rupture geometry...)
			equivalentMainshock = ETAS_StatsCalc.fitMainshockLineSource(forecastModel.mainShock, faultTrace, stressDrop);
		} else if (fitType.equals("custom") && faultTrace != null && faultTrace.size() > 1) {
			equivalentMainshock = ETAS_StatsCalc.fitMainshockLineSource(forecastModel.mainShock, faultTrace, stressDrop);
		}
		//else {equivalentMainshock was already constructed with a faultTrace made up of the mainshock hypocenter}

		if(D){
			FaultTrace trace = equivalentMainshock.getFaultTrace();
			for (int i = 0; i < trace.size(); i++){
				Location loc = trace.get(i);
				System.out.println(loc);
			}
		}

		// compute rates at each point in the rate map for the mainshock source
		if(D) System.out.println("computing MS rate integral from day " + forecastModel.forecastMinDays + " to " + forecastMaxDays);

		double x, y, x0, y0, newVal;
		double t0 = 0;
		double mag0 = equivalentMainshock.getMag();

		for (int j=0; j<gridData.size(); j++) {
			Location gridLoc = gridData.getLocation(j);
			x = gridLoc.getLongitude();
			y = gridLoc.getLatitude();

			newVal = rateXY(x,y,t0,equivalentMainshock, forecastModel.forecastMinDays, forecastMaxDays, seismogenicDepth);
			gridData.set(j, gridData.get(j) + newVal);
		}

		// compute rates at each point in the rate map for the aftershock sources
		if(D) System.out.println("computing AS rate integral for time " + forecastModel.forecastMinDays + " to " + forecastMaxDays);

		// trim the fluff (don't compute rate for very small aftershocks
		Iterator<ObsEqkRupture> iter = aftershockFitList.listIterator();
		while (iter.hasNext()){
			ObsEqkRupture as = iter.next();
			if (as.getMag() < forecastModel.mainShock.getMag() - maxDeltaMag)
				iter.remove();
		}

		for (ObsEqkRupture rup : aftershockFitList){
			x0 = rup.getHypocenterLocation().getLongitude();
			y0 = rup.getHypocenterLocation().getLatitude();
			mag0 = rup.getMag();
			t0 = (rup.getOriginTime() - forecastModel.mainShock.getOriginTime()) / ETAS_StatsCalc.MILLISEC_PER_DAY;

			for (int i=0; i<gridData.size(); i++) {
				Location gridLoc = gridData.getLocation(i);
				x = gridLoc.getLongitude();
				y = gridLoc.getLatitude();

				newVal = rateXY(x,y,t0,mag0,x0,y0, stressDrop, forecastModel.forecastMinDays, forecastMaxDays, seismogenicDepth);
				gridData.set(i, gridData.get(i) + newVal);
			}
		}

		//normalize the gridData rate map to give the correct total forecast number
		if(D) System.out.println("Mc: " + forecastModel.magComplete);
		double gridSum = 0; 

		//		double rateTotal = this.getModalNumEvents(this.magComplete, forecastMinDays, forecastMaxDays);
		double rateTotal = forecastModel.getCumNumFractileWithAleatory(new double[]{0.5}, forecastModel.magComplete, forecastModel.forecastMinDays, forecastMaxDays)[0];
		gridSum = gridData.getSumZ();
		gridData.scale(rateTotal/gridSum);

		rateModel = gridData;
		return gridData;
	}

	/** smooth the rate map to give the rate/probability of earthquakes of a given size and distance
	 *  "spacing" must be specified and must match the spacing in gridData 
	 */
	public GriddedGeoDataSet calculateSmoothRateModel(double magPlot, double distance){
		if (rateModel==null){
			System.err.println("Gridded rate forecastModel must be calculated first.");
			return null;
		}
		
		GriddedGeoDataSet smoothGridData = new GriddedGeoDataSet(rateModel.getRegion(), false);

		Location pt1;
		Location pt2;
		double dx, dy, r2;
		double d2 = distance*distance/111.111/111.111; //km --> deg
		double lon, lat, lon0, lat0, latFactor;
		double rateSum;

		// set up timer/time estimator
		double toc, timeEstimate;
		Stopwatch watch = Stopwatch.createStarted();
		int warnTime = 3;
		String initialMessageString = "Computing smoothed rate map. ";
		
		for (int i=0 ; i < smoothGridData.size(); i++){
			pt1 = smoothGridData.getLocation(i);
			lon0 = pt1.getLongitude();
			lat0 = pt1.getLatitude();
			latFactor = Math.cos(Math.toRadians(lat0));

			//find the grid points near enough to include
			rateSum = 0;
			for (int j = 0; j < smoothGridData.size(); j++){
				//start with the quick estimate
				pt2 = smoothGridData.getLocation(j);
				lon = pt2.getLongitude();
				lat = pt2.getLatitude();
				dx = (lon-lon0)*latFactor;
				dy = (lat-lat0);

				r2 = dx*dx + dy*dy;
				if (r2 < d2){
					rateSum += rateModel.get(j);
				}
			}

			// run the timer to see how long this is going to take
			toc = watch.elapsed(TimeUnit.SECONDS);
			if(toc > warnTime){
				warnTime += 10;

				long count = i;
				long total = smoothGridData.size();
				timeEstimate = toc * (double)total/ (double)count;
				System.out.format(initialMessageString + "Approximately %d seconds remaining...\n", (int) ((timeEstimate - toc)));
				initialMessageString = "...";
				
				if (forecastModel.progress != null){
					//									progress.updateProgress(count, total, String.format("%d%% complete. %d seconds remaining", (int) (((double) count)/((double) total) * 100), (int) ((timeEstimate - toc)/1000)));
					forecastModel.progress.setProgressMessage(String.format("%d%% complete. %d seconds remaining", (int) (((double) count)/((double) total) * 100),(int) ((timeEstimate - toc))));
					forecastModel.progress.pack();
				}
			}
			// scale to new magnitude reference
			rateSum *= Math.pow(10, -forecastModel.get_b()*(magPlot - forecastModel.magComplete));
			double probSum = 1d - Math.exp(-rateSum);

			smoothGridData.set(i, probSum);
		}
		
		return smoothGridData;
	}

	/** smooth the rate map to give the rate/probability of earthquakes of a given size and distance
	 *  "spacing" must be specified and must match the spacing in gridData 
	 */
	public GriddedGeoDataSet calculateMMIRateModel(double mmiRef){
		
		GriddedGeoDataSet smoothGridData = new GriddedGeoDataSet(rateModel.getRegion(), false);

		Location pt1;
		Location pt2;
		double dx, dy, r, mag;
		double lon, lat, lon0, lat0, latFactor;
		double rateSum;
		double magRef = forecastModel.magComplete;

		// set up timer/time estimator
		double toc, timeEstimate;
		Stopwatch watch = Stopwatch.createStarted();
		int warnTime = 3;
		String initialMessageString = "Computing smoothed rate map. ";
		
		mmiModel = null; //reset the mmi-magnitude interpolation forecastModel
		
		for (int i=0 ; i < smoothGridData.size(); i++){
			pt1 = smoothGridData.getLocation(i);
			lon0 = pt1.getLongitude();
			lat0 = pt1.getLatitude();
			latFactor = Math.cos(Math.toRadians(lat0));

			//find the grid points near enough to include
			rateSum = 0;
			for (int j = 0; j < smoothGridData.size(); j++){
				//start with the quick estimate
				pt2 = smoothGridData.getLocation(j);
				lon = pt2.getLongitude();
				lat = pt2.getLatitude();
				dx = (lon-lon0)*latFactor;
				dy = (lat-lat0);

				r = Math.sqrt(dx*dx + dy*dy)*111.1111;
				
				mag = getMagForMMI(mmiRef, r);
				
				rateSum += rateModel.get(j) * Math.pow(10, -forecastModel.get_b()*(mag - forecastModel.magComplete));;
			}
			

			// run the timer to see how long this is going to take
			toc = watch.elapsed(TimeUnit.SECONDS);
			if(toc > warnTime){
				warnTime += 10;

				long count = i;
				long total = smoothGridData.size();
				timeEstimate = toc * (double)total/ (double)count;
				System.out.format(initialMessageString + "Approximately %d seconds remaining...\n", (int) ((timeEstimate - toc)));
				initialMessageString = "...";
				
				if (forecastModel.progress != null){
					//									progress.updateProgress(count, total, String.format("%d%% complete. %d seconds remaining", (int) (((double) count)/((double) total) * 100), (int) ((timeEstimate - toc)/1000)));
					forecastModel.progress.setProgressMessage(String.format("%d%% complete. %d seconds remaining", (int) (((double) count)/((double) total) * 100),(int) ((timeEstimate - toc))));
					forecastModel.progress.pack();
				}
			}
			smoothGridData.set(i, rateSum);
		}

		return smoothGridData;
	}
	
	private double[][] mmiModel;
	private double maxDist = Double.POSITIVE_INFINITY;
	
	private double getMagForMMI(double mmiRef, double D){
		
		double mag;
		if (mmiModel == null){
			int npts = 100;
			double[] distances = ETAS_StatsCalc.logspace(1e-5, 1e3, npts);
			distances[0] = 0;
			
			mmiModel = new double[2][npts];
			
			for (int i = 0; i < distances.length; i++){
				mmiModel[0][i] = distances[i]; 
				
				D = distances[i];

				// this one uses Atkinson and Wald 2007, not corrected for sigma
				String loc = "notCEUS";
				double[] c;
				double h, Rt;
				if (loc.equals("CEUS")){
					// CEUS values
					c = new double[]{11.72, 2.36, 0.1155, -0.44, -0.002044, 2.31, -0.479};
					h = 17.0;
					Rt = 80.0;
				} else {
					// call it california
					c = new double[]{12.27, 2.270, 0.1304, -1.30, -0.0007070, 1.95, -0.577};
					h = 14.0;
					Rt = 30.0;
				}

				double R = Math.sqrt(D*D + h*h);
				double B;
				if (R < Rt)
					B = 0.0;
				else
					B = Math.log10(R/Rt);

				double a = c[2];
				double b = c[1] - 12.0*c[2] + c[6]*Math.log10(R);
				double d = c[0] - 6.0*c[1] + 36.0*c[2] + c[3]*Math.log10(R) + c[4]*R + c[5]*B - mmiRef;

				mag = (-b + Math.sqrt(b*b - 4.0*a*d))/2.0/a;

				if (d > b*b/4.0/a){
					// it can never happen
					mag = Double.POSITIVE_INFINITY;
					// update the maxDist -- the distance beyond which it can never happen
					if(D < maxDist)
						maxDist = D;
				}
				
				mmiModel[1][i] = mag;
			}
		}
		
//		return mmiModel.getInterpolatedY_inLogXLogYDomain(D);
		if(D >= maxDist)
			return Double.POSITIVE_INFINITY;
		else{
			int pt = mmiModel[0].length-1;
			while ( D < mmiModel[0][pt] && pt > 0){
				pt--;
			}
			return mmiModel[1][pt];
		}
		
//			if (pt == mmiModel[0].length-1)
//				return mmiModel[1][pt];
//			else{
//				double baseMag = mmiModel[1][pt];
//				double topMag = mmiModel[1][pt+1];
//				
//				double baseDist = mmiModel[0][pt];
//				double topDist = mmiModel[0][pt];
//				
//				return baseMag + (topMag-baseMag)/(topDist-baseDist) * D;
//			}
//
//		}
	}
	
	
	/**	compute rate at one target point for one point source
	 * 
	 */
	private double rateXY(double lon, double lat, double t0, double mag0, double lon0, double lat0, double stressDrop, double ts, double te, double H){
		// quick convert to x,y
		double dx = (lon-lon0)*Math.cos(Math.toRadians(lat0))*111.111;	//using grid adjusted to mainshock lat. will fail near poles
		double dy = (lat-lat0)*111.111;

		// constants
		double r = Math.sqrt(dx*dx + dy*dy);
		double d = ETAS_StatsCalc.magnitude2radius(mag0, stressDrop);

		// compute productivity for weighting this event
		double productivity = Math.pow(10d, forecastModel.getMaxLikelihood_a() + 1d*(mag0 - forecastModel.magComplete));
		// compute rate, integrated over seismogenic depth H
		double spatialDecay = H / (d*d + r*r) / Math.pow(H*H/4d + r*r + d*d, 1d/2d) * d/(2*Math.PI);

		double p = forecastModel.getMaxLikelihood_p();
		double c = forecastModel.getMaxLikelihood_c();
		double timeIntegral = 1d/(1d - p) * ( Math.pow(te - t0 + c, 1d-p) - Math.pow(ts - t0 + c, 1d-p) );  
		//        
		double rate = productivity *  timeIntegral * spatialDecay;

		return rate;
	}

	/** Compute rate at one target point for one finite mainshock source
	 * 
	 */
	private double rateXY(double lon, double lat, double t0, ETASEqkRupture equivalentMainshock, double ts, double te, double H){
		FaultTrace trace = equivalentMainshock.getFaultTrace();

		// quick convert to x,y
		double dx, dy, r, d;
		double productivity, spatialDecay;
		double p, c;
		double timeIntegral;
		double rate = 0;

		Location loc;
		double mag0 = equivalentMainshock.getMag();
		double reach = equivalentMainshock.getSpatialKernelDistance();

		// compute time integral
		p = forecastModel.getMaxLikelihood_p();
		c = forecastModel.getMaxLikelihood_c();
		timeIntegral = 1d/(1d - p) * ( Math.pow(te - t0 + c, 1d-p) - Math.pow(ts - t0 + c, 1d-p) );  

		// compute productivity for this event
		productivity = Math.pow(10d, forecastModel.getMaxLikelihood_a() + 1d*(mag0 - forecastModel.magComplete));
		productivity /= (double) trace.size();

		for (int i = 0; i < trace.size(); i++){
			loc = trace.get(i);

			double lon0 = loc.getLongitude();
			double lat0 = loc.getLatitude();

			// quick convert to x,y
			dx = (lon-lon0)*Math.cos(Math.toRadians(lat0))*111.111;
			dy = (lat-lat0)*111.111;

			// constants
			r = Math.sqrt(dx*dx + dy*dy);
			d = reach;

			// compute rate, integrated over seismogenic depth H
			spatialDecay = H / (d*d + r*r) / Math.pow(H*H/4d + r*r + d*d, 1d/2d) * d/(2*Math.PI);

			//        
			rate += productivity *  timeIntegral * spatialDecay;
		}
		return rate;
	}
	
	public static List<PolyLine> getContours(GriddedGeoDataSet gridData, int nc){
		double[] contourLevels = ETAS_StatsCalc.linspace(0,gridData.getMaxZ(),nc);
		return getContours(gridData, contourLevels);
	}
	
	public static List<PolyLine> getContours(GriddedGeoDataSet gridData, double[] contourLevels){
		double spacing = gridData.getRegion().getSpacing();
		
		// construct contours
		double latmin = gridData.getRegion().getMinLat();
		double latmax = gridData.getRegion().getMaxLat();
		double lonmin = gridData.getRegion().getMinLon();
		double lonmax = gridData.getRegion().getMaxLon();
		int nlat = (int) ((latmax-latmin)/spacing + 0.5)+1;
		int nlon = (int) ((lonmax-lonmin)/spacing + 0.5)+1;
		double[] X = ETAS_StatsCalc.linspace(lonmin,lonmax,nlon);
		double[] Y = ETAS_StatsCalc.linspace(latmin,latmax,nlat);
		double[][] S0 = new double[nlat][nlon];
		int[][] S1 = new int[nlat][nlon];

		// figure out which way the gridded dataset is gridded
		double lat, lat0;
		double lon, lon0;
		lat0 = gridData.getLocation(0).getLatitude();
		lon0 = gridData.getLocation(0).getLatitude();
		int x = -1, y = 0;

		for (int i = 0; i < gridData.size(); i++){
			// starts with minLat, minLon, incrementes lon, then lat
			lat = gridData.getLocation(i).getLatitude();
			lon = gridData.getLocation(i).getLongitude();


			if (Math.abs(lat - lat0) < 1e-6){
				//same y index, iterate x index
				x++;
			} else {
				//new 
				x=0;
				y++;
				lat0 = lat;
			}
			
			S0[y][x] = gridData.get(i);
			S1[y][x] = 1;
			X[x] = lon;
			Y[y] = lat;
		}

		double undefData = 0d;
		int nc = contourLevels.length;
		List<Border> borders = Contour.tracingBorders(S0, X, Y, S1, undefData);
		List<PolyLine> cntr = Contour.tracingContourLines(S0, X, Y, nc, contourLevels, undefData, borders, S1);
		
		return cntr;
	}

	/* Write a contour as KML document 
	 * 
	 */
	public static void writeContoursAsKML(List<PolyLine> contours, String name, String units, File outputFile, CPT cpt){
		double labelSize = 0.8;
		boolean D = false; //debug
		Color color = new Color(0,0,0);
		
		
		int width = 2;
		int labelCount = 3;
		int labelLimit = 8;
		int labelSpaceMin = 4;
		
		StringBuilder outputString = new StringBuilder();
		
		// start building the file
		// insert header
		outputString.append(header(name, labelSize, color, width));
		if(D) System.out.println("Printing " + contours.size() + " contours to file.");
	
		double prevLevel = 0.0;
		boolean isNewLevel = true;
		for (int i = 0; i < contours.size(); i++){
			PolyLine contour = contours.get(i);
			int contourLength = contour.PointList.size();
			int labelSpace = Math.max(labelSpaceMin, contour.PointList.size()/(labelCount));
			int startPt = (int) (Math.random()*labelSpace);
			
			double level = contours.get(i).Value;
			if (Math.abs(level - prevLevel) < 0.001) {
				isNewLevel = false;
			} else {
				isNewLevel = true;
			}
			prevLevel = level;
				
			double lat, lon, depth;
			color = cpt.getColor((float) level);
			
			// start the contour line and label line
			StringBuilder contourString = new StringBuilder();
			StringBuilder labelString = new StringBuilder();
			contourString.append(beginLine(level, color));
			for(int j = 0; j < contour.PointList.size(); j++){
				lon = unwrap(contour.PointList.get(j).X);
				lat = contour.PointList.get(j).Y;
				
				depth = 0.0;
				contourString.append(lon + "," + lat + "," + depth + " ");
				
				// add another label if enough points have been passed
				if ( (contourLength > labelLimit && Math.floorMod(j, labelSpace) == startPt
						&& j <= contour.PointList.size() - labelSpace) 
							|| isNewLevel){
					labelString.append(label(level, lon, lat, depth, units));
				};
				isNewLevel = false;
			}
			// end the contour line
			contourString.append(endLine());
			
			outputString.append(labelString);
			outputString.append(contourString);
			
		}
		
		//add aftershocks so far as pts?
		
		
		outputString.append(endFile());
		
//		if (D) System.out.println(outputString);
		
		// write file
		FileWriter fw;
		try {
			fw = new FileWriter(outputFile, false);
		} catch (IOException e1) {
			//				e1.printStackTrace();
			System.err.println("Couldn't save to file " + outputFile.getAbsolutePath());
			return;
		}

		try {
			fw.append(outputString);
		} catch (IOException e) {
			//					e.printStackTrace();
			System.err.println("Couldn't save to file " + outputFile.getAbsolutePath());
		}
		
		try {
			fw.close();
		} catch (IOException e) {
			//				e.printStackTrace();
			System.err.println("Problem closing file.");
		}
	
	}
	
	private static String header(String name, double labelSize, Color color, int width){
		return "" 
			+"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
		    +"<kml xmlns=\"http://www.opengis.net/kml/2.2\"\n"
		    +" xmlns:gx=\"http://www.google.com/kml/ext/2.2\"\n"
		    +" xmlns:kml=\"http://www.opengis.net/kml/2.2\"\n"
		    +" xmlns:atom=\"http://www.w3.org/2005/Atom\">\n"
		    +"<Document>\n"
		    +"	<name>" + name + ".kml</name>\n"
		    +"	<Style id=\"sn_noicon\">\n"
		    +"		<IconStyle>\n"
		    +"			<Icon>\n"
		    +"			</Icon>\n"
		    +"		</IconStyle>\n"
		    +"		<LabelStyle>\n"
		    +"			<scale>" + labelSize + "</scale>\n"
		    +"		</LabelStyle>\n"
		    +"	</Style>\n"
		    +"	<Style id=\"linestyle\">\n"
		    +"			<LineStyle>\n"
		    +"				<width>" + width + "</width>\n"
		    +"			</LineStyle>\n"
		    +"	</Style>\n";
	}
	
	private static String endFile(){
		return ""
				+"</Document>\n"
				+"</kml>\n";
	}
			
	private static String beginLine(double level, Color color){
		return ""
				+"	<Placemark>\n"
				+"		<name>" + String.format("%.2f", level) + "</name>\n"
				+"		<styleUrl>#linestyle</styleUrl>\n"
				+"		<Style>\n"
				+" 			<LineStyle>\n"
				+"				<color>#ff" + color2hex(color) + "</color>\n"
				+" 			</LineStyle>\n"
				+" 		</Style>\n"
				+"		<LineString>\n"
				+"			<tessellate>1</tessellate>\n"
				+"			<altitudeMode>clampToSeaFloor</altitudeMode>\n"
				+"			<gx:altitudeMode>clampToSeaFloor</gx:altitudeMode>\n"
				+"			<coordinates>\n";
	}
	
	private static String endLine(){
		return "\n"
		    +"			</coordinates>\n"
		    +"		</LineString>\n"
		    +"	</Placemark>\n";
	}		

	private static String label(double level, double lon, double lat, double depth, String units){
		
		String levelStr;
		if (level < 10)
			levelStr = String.format("%.1f", Math.floor(level/0.1)*0.1) + units; //rounds by ones
		else
			levelStr = String.format("%.0f", Math.floor(level/1)*1) + units; //rounds by ones
		
		
//		if (level < 1.001) levelStr = "<1" + units;
//		if (level > 98.99) levelStr = ">99" + units;
		
		return ""
				+"	<Placemark>\n"
				+"		<name>" +  levelStr + "</name>\n"
				+"		<styleUrl>#sn_noicon</styleUrl>\n"
				+"		<Point>\n"
				+"			<altitudeMode>clampToSeaFloor</altitudeMode>\n"
				+"			<gx:altitudeMode>clampToSeaFloor</gx:altitudeMode>\n"
				+"			<coordinates>" + unwrap(lon) + "," + lat + "," + depth +"</coordinates>\n"
				+"		</Point>\n"
				+"	</Placemark>\n";
	}

	private static double unwrap(double lon){
		if (lon > 180) lon -= 360;
		
		return lon;
	}
	
	private static String color2hex(Color color){
		int[] rbg = new int[]{color.getBlue(), color.getGreen(), color.getRed()};
		
		StringBuilder hexString = new StringBuilder();
		hexString.append("");
		for (int i = 0; i< rbg.length; i++){
			hexString.append(String.format("%02x", rbg[i]));
		}
		
		return hexString.toString();
	}
	
	
}
