package scratch.aftershockStatisticsETAS.griddedInterpGMPE;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.TimeUnit;

import javax.swing.JFrame;
import javax.swing.JOptionPane;

import org.jfree.data.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.siteData.impl.WaldAllenGlobalVs30;
import org.opensha.commons.data.xyz.ArbDiscrGeoDataSet;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.exceptions.WarningException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.nshmp2.erf.source.PointSource13b;
import org.opensha.nshmp2.util.FocalMech;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGV_Param;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;

import scratch.aftershockStatisticsETAS.ETAS_ShakingForecastCalc;

public class GriddedInterpGMPE_Calc {
	
	private Boolean D = false; //debug
	private ScalarIMR gmpe;
	private Boolean promptForLongCalc = false;
	
	private DistanceInterpolator distInterp; // index 0
	private IntensityMeasureLevelInterpolator imlInterp; // index N-1
	
	private List<AbstractGMPEInterpolation<?>> allInterps;
	
	private DiscretizedFunc xVals;
	private DiscretizedFunc logXVals;
	
	private IncrementalMagFreqDist inputMFD;
	
	private NDimArrayCalc arrayCalc;
	private double[] allExceedRates; // source non-exceedance rates
	private NDimensionalLinearInterpolation interpolator;
	
	public GriddedInterpGMPE_Calc(ScalarIMR gmpe, DiscretizedFunc xVals, double b, double minMag, double maxMag, int numMag,
			DistanceInterpolator distInterp, AbstractGMPEInterpolation<?>... otherInterps) {
		this(gmpe, xVals, new GutenbergRichterMagFreqDist(b, 1d, minMag, maxMag, numMag), distInterp, otherInterps);
	}
	
	public GriddedInterpGMPE_Calc(ScalarIMR gmpe, DiscretizedFunc xVals, IncrementalMagFreqDist inputMFD,
			DistanceInterpolator distInterp, AbstractGMPEInterpolation<?>... otherInterps) {
		this.gmpe = gmpe;
		
		this.xVals = xVals;
		logXVals = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : xVals)
			logXVals.set(Math.log(pt.getX()), 1d);
		
		this.distInterp = distInterp;
		this.imlInterp = new IntensityMeasureLevelInterpolator("IML", xVals, true);
		
		allInterps = new ArrayList<>();
		allInterps.add(distInterp);
		for (AbstractGMPEInterpolation<?> o : otherInterps)
			if (o != null)
				allInterps.add(o);
		allInterps.add(imlInterp);
		
		this.inputMFD = inputMFD;
		
		inputMFD.scaleToCumRate(0, 1d);
//		System.out.print("Input MFD");
//		for (int i=0; i<inputMFD.size(); i++)
//			System.out.println("\t"+(float)+inputMFD.getX(i)+"\t"+(float)inputMFD.getY(i));
		
		int[] dimensions = new int[allInterps.size()];
		for (int i=0; i<dimensions.length; i++)
			dimensions[i] = allInterps.get(i).getNumBins();
		arrayCalc = new NDimArrayCalc(dimensions);
		interpolator = new NDimensionalLinearInterpolation(dimensions.length);
		
		allExceedRates = new double[arrayCalc.rawArraySize()];
		for (int i=0; i<allExceedRates.length; i++)
			allExceedRates[i] = Double.NaN;
	}
	
	public void precalc(double duration, double[] depths, Map<FocalMech, Double> mechWtMap) {
		Location loc = new Location(0d, 0d);
		
//		Map<FocalMech, Double> mechWtMap = new HashMap<>();
//		double wtEach = 1d/FocalMech.values().length;
//		for (FocalMech mech : FocalMech.values())
//			mechWtMap.put(mech, wtEach);
		PointSource13b source = new PointSource13b(loc, inputMFD, duration, depths, mechWtMap);
		
		if(D) System.out.println("Precalculating GMPE for "+allInterps.size()+" dimensions, "+arrayCalc.rawArraySize()+" values");
		
		Set<EqkRupture> ruptures = new HashSet<>();
		for (ProbEqkRupture rup : source)
			ruptures.add((EqkRupture)rup.clone());
		
		Site site = new Site(loc);
		site.addParameterList(gmpe.getSiteParams());
		gmpe.setSite(site);
		
		try {
			gmpe.setEqkRupture(ruptures.iterator().next());
		} catch (WarningException e) {
			System.err.println(e.getMessage());
			//do nothing then
		}
		
		precalcRecursive(new int[0], source);
		
//		numRupsTrack = new MinMaxAveTracker();
		
//		for (FocalMech mech : FocalMech.values()) {
//			double[] exceedProbs = new double[arrayCalc.rawArraySize()];
//			
//			HashSet<EqkRupture> mechRups = new HashSet<>();
//			for (EqkRupture rup : ruptures) {
//				if (FocalMechInterpolator.forRake(rup.getAveRake()) == mech)
//					mechRups.add(rup);
//			}
//			
//			precalcRecursive(new int[0], source, mechRups, exceedProbs);
//			int numNonZero = 0;
//			for (double val : exceedProbs)
//				if (val > 0)
//					numNonZero++;
//			System.out.println(numNonZero+"/"+exceedProbs.length+" are non-zero for "+mech);
//			this.exceedProbs.put(mech, exceedProbs);
//		}
//		System.out.println("Done precalculating!");
//		System.out.println("Num Rups tracker: "+numRupsTrack);
	}
	
//	private MinMaxAveTracker numRupsTrack;
	
	private void precalcRecursive(int[] upstreamIndexes, PointSource13b source) {
		int curIndex = upstreamIndexes.length;
		AbstractGMPEInterpolation<?> interp = allInterps.get(curIndex);
		int[] indexes = Arrays.copyOf(upstreamIndexes, curIndex+1);
		
		if (interp instanceof IntensityMeasureLevelInterpolator) {
			// we're at the IML level, time to actually calculate
			Preconditions.checkState(indexes.length == allInterps.size(), "IML interpolator must be last");
			
			// all parameters except for the actual site should be set now
			
			Preconditions.checkState(logXVals.size() == interp.getNumBins());
			double[] sourceExceedRates = new double[logXVals.size()];
			
			for (ProbEqkRupture rup : source) {
				try {
					gmpe.setEqkRupture(rup);
				} catch (WarningException e) {
					System.err.println(e.getMessage());
					// do nothing then
				}
				
				gmpe.getExceedProbabilities(logXVals);
				
				double rupProb = rup.getProbability();
				double rupRate = -Math.log(1 - rupProb);
				
				for (int i=0; i<sourceExceedRates.length; i++)
					sourceExceedRates[i] = sourceExceedRates[i] + rupRate * logXVals.getY(i);
			}
			
			// now fold my values into the global array
			for (int i=0; i<sourceExceedRates.length; i++) {
				indexes[indexes.length-1] = i;
				int arrayIndex = arrayCalc.getIndex(indexes);
				Preconditions.checkState(Double.isNaN(allExceedRates[arrayIndex]), "Value already set?");
				allExceedRates[arrayIndex] = sourceExceedRates[i];
			}
			
//			indexes[0] = i; // IML index
//			double exceedProb = interpolator.interpolate(exceedProbs, arrayCalc, indexes);
////			if (Math.random() < 0.00001)
////				System.out.println("Indexes: "+(float)indexes[0]+" "+(float)indexes[1]
////					+" "+(float)indexes[2]+" "+(float)indexes[3]+" ==> "+exceedProb);
//			curves[s].set(i, curves[s].getY(i)*Math.pow(1-rup.getProbability(), exceedProb));
		} else {
			Preconditions.checkState(indexes.length < allInterps.size(), "We're at the end but not an IML interpolator");
			
			for (int i=0; i<interp.getNumBins(); i++) {
				indexes[indexes.length-1] = i;
				interp.setGMPE_Params(gmpe, source, i);
				
				precalcRecursive(indexes, source);
			}
		}
	}
	
	private volatile boolean stopRequested = false;
	public DiscretizedFunc[] calc(GeoDataSet griddedTotCumRates, List<Site> sites) {
		DiscretizedFunc[] curves = new DiscretizedFunc[sites.size()];
		
		for (int i=0; i<sites.size(); i++) {
			curves[i] = new LightFixedXFunc(xVals);
			// initialize the hazard function to 1.0
			for (int j=0; j<curves[i].size(); j++)
				curves[i].set(j, 1d);
		}
		
		double[] indexes = new double[allInterps.size()];
		
		double inputMFD_totCumRate = inputMFD.getCumRate(0);
		
		// set up timer/time estimator
		double toc, timeEstimate, n;
		Stopwatch watch = Stopwatch.createStarted();
		int warnTime = 3;
		boolean userWarned = false;
		double deltaT = 0; //this will record the time spent waiting for the dialog box.
		String initialMessageString = "Calculating shaking map. ";
		
		for (int g=0; g<griddedTotCumRates.size(); g++) {
			Location sourceLoc = griddedTotCumRates.getLocation(g);
			double totCumRate = griddedTotCumRates.get(g);
			double rateScalar = totCumRate / inputMFD_totCumRate;
			
			for (int s=0; s<sites.size(); s++) {
				Site site = sites.get(s);
				double dist = LocationUtils.horzDistanceFast(sourceLoc, site.getLocation());
				if (dist > distInterp.getMax())
					continue;
				if (dist == 0 )
					dist = distInterp.getMin();

				indexes[0] = distInterp.getInterpolatedBinIndex(dist);
				
				for (int j=1; j<allInterps.size()-1; j++)
					indexes[j] = allInterps.get(j).detectInterpolatedBinIndex(null, site);
				
				for (int i=0; i<xVals.size(); i++) {
					indexes[indexes.length-1] = i;
					double sourceExceedRate = interpolator.interpolate(allExceedRates, arrayCalc, indexes);
					
					// now we scale to the actual rate of this source
					sourceExceedRate *= rateScalar;
					
					double sourceExceedProb = 1d - Math.exp(-sourceExceedRate);
					
					double sourceNonExceedProb = 1d - sourceExceedProb;
					
					curves[s].set(i, curves[s].getY(i)*sourceNonExceedProb);
				}
		
				// run the timer to see how long this is going to take
				toc = watch.elapsed(TimeUnit.SECONDS) - deltaT;
				if(toc > warnTime){
					long count = (g)*(sites.size()) + s;
					long total = (sites.size() * griddedTotCumRates.size());
					timeEstimate = toc * total/count;
					System.out.format(initialMessageString + "Approximately %d seconds remaining...\n", (int) ((timeEstimate - toc)));
					initialMessageString = "...";
					
					// if the time estimate is more than 20 seconds, ask if user wants to quit
					if (!userWarned && timeEstimate > 30 && promptForLongCalc) { // only the first time around and if it'll take more than a minute
						userWarned = true;
						// launch a dialog as a new thread
						String message = "It will take approximately " + (int) timeEstimate + " seconds to complete each map at this resolution.\n"
								+ "If plotting MMI multiply this estimate by a factor of 2.\n";
						message += "Are you sure you wish to continue with the current grid spacing of \u0394 (km)?";
						
						String title = "Warning";
						
						double t1 = watch.elapsed(TimeUnit.SECONDS); 
						try {
							int ret = JOptionPane.showConfirmDialog(null, message, title, JOptionPane.OK_CANCEL_OPTION);
							if (ret == JOptionPane.CANCEL_OPTION)
								stopRequested = true;
						} catch (Exception e) {
							System.err.println("Error displaying error message!");
							e.printStackTrace();
						}
						deltaT = watch.elapsed(TimeUnit.SECONDS) - t1;
					}
					warnTime += 10;
				}	
				
				if (stopRequested) {
					System.out.println("Map calculation terminated prematurely");
					return null;
				}
			}
		}
		
//	
//		for (ProbEqkSource source : sources) {
//			// calculate distances
//			Preconditions.checkState(source.getSourceSurface() instanceof PointSurface, "Only point sources supported");
//			Location sourceLoc = ((PointSurface)source.getSourceSurface()).getLocation();
//			double[] distances = new double[sites.size()];
//			for (int s=0; s<sites.size(); s++)
//				distances[s] = LocationUtils.horzDistanceFast(sourceLoc, sites.get(s).getLocation());
//			for (ProbEqkRupture rup : source) {
//				// these 2 are common to all
//				Preconditions.checkState(rup.getMag() >= magInterp.getMin() && rup.getMag() <= magInterp.getMax(),
//						"Rup mag=%s outisde of range [%s %s]", rup.getMag(), magInterp.getMin(), magInterp.getMax());
//				indexes[1] = magInterp.getInterpolatedBinIndex(rup.getMag());
//				FocalMech mech = FocalMechInterpolator.forRake(rup.getAveRake());
//				double[] exceedProbs = this.exceedProbs.get(mech);
//				for (int s=0; s<sites.size(); s++) {
//					Site site = sites.get(s);
//					if (distances[s] > distInterp.getMax())
//						continue;
////					System.out.println("dist: "+distances[s]);
//					// set dist and any other params
//					indexes[2] = distInterp.getInterpolatedBinIndex(distances[s]);
//					for (int j=2; j<allInterps.size(); j++)
//						indexes[j+1] = allInterps.get(j).detectInterpolatedBinIndex(null, site);
//					
//					// multi-dimensional interpolation
////					double mean = interpolator.interpolate(means, arrayCalc, indexes);
////					double stdDev = interpolator.interpolate(stdDevs, arrayCalc, indexes);
////					gmpe.setSite(site);
////					gmpe.setEqkRupture(rup);
////					double mean = gmpe.getMean();
////					double stdDev = gmpe.getStdDev();
//					
////					if (debugMeanScatter != null) {
////						gmpe.setSite(site);
////						gmpe.setEqkRupture(rup);
////						double mean2 = gmpe.getMean();
////						double stdDev2 = gmpe.getStdDev();
////						debugMeanScatter.set(mean2, mean);
////						debugStdDevScatter.set(stdDev2, stdDev);
////					}
//					
//					for (int i=0; i<xVals.size(); i++) {
//						indexes[0] = i; // IML index
//						double exceedProb = interpolator.interpolate(exceedProbs, arrayCalc, indexes);
////						if (Math.random() < 0.00001)
////							System.out.println("Indexes: "+(float)indexes[0]+" "+(float)indexes[1]
////								+" "+(float)indexes[2]+" "+(float)indexes[3]+" ==> "+exceedProb);
//						curves[s].set(i, curves[s].getY(i)*Math.pow(1-rup.getProbability(), exceedProb));
//					}
//				}
//			}
//		}
		
		// convert to exceedance probabilities (currently non-exceedance)
		for (DiscretizedFunc curve : curves)
			for (int i=0; i<curve.size(); i++)
				curve.set(i, 1d - curve.getY(i));
		
		return curves;
	}
	
	public void setPromptForLongCalc(boolean prompt) {
		this.promptForLongCalc = prompt;
	}
	
	public boolean getPromptForLongCalc() {
		return promptForLongCalc;
	}
	
	
	
////	private DefaultXY_DataSet debugMeanScatter = new DefaultXY_DataSet();
////	private DefaultXY_DataSet debugStdDevScatter = new DefaultXY_DataSet();
//	
//	private DiscretizedFunc[] calcTraditional(ERF erf, List<Site> sites, DiscretizedFunc xVals) {
//		HazardCurveCalculator calc = new HazardCurveCalculator();
//		DiscretizedFunc logXVals = new ArbitrarilyDiscretizedFunc();
//		for (Point2D pt : xVals)
//			logXVals.set(Math.log(pt.getX()), 1);
//		
//		DiscretizedFunc[] curves = new DiscretizedFunc[sites.size()];
//		
//		for (int i=0; i<sites.size(); i++) {
//			DiscretizedFunc logCurve = new LightFixedXFunc(logXVals);
//			calc.getHazardCurve(logCurve, sites.get(i), gmpe, erf);
//			curves[i] = new LightFixedXFunc(xVals);
//			// initialize the hazard function to 1.0
//			for (int j=0; j<curves[i].size(); j++)
//				curves[i].set(j, logCurve.getY(j));
//		}
//		
//		return curves;
//	}
//	
//	public static void main(String[] args) throws IOException {
////		ScalarIMR gmpe = AttenRelRef.NGAWest_2014_AVG_NOIDRISS.instance(null);
//		ScalarIMR gmpe = AttenRelRef.BSSA_2014.instance(null);
//		gmpe.setParamDefaults();
////		gmpe.setIntensityMeasure(PGA_Param.NAME);
//		gmpe.setIntensityMeasure(PGV_Param.NAME);
//
//		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(gmpe.getIntensityMeasure());
////		IntensityMeasureLevelInterpolator imlInterp =
////				new IntensityMeasureLevelInterpolator(gmpe.getIntensityMeasure().getName(), xVals, true);
//		
//		double calcSpacing = 0.1;
//		DistanceInterpolator distInterp = new DistanceInterpolator(true, calcSpacing, 500d, 100);
//		
//		double refMag = 5d;
//		double maxMag = 8.5d;
//		double magDelta = 0.1;
//		int numMag = (int)((maxMag - refMag)/magDelta + 0.5) + 1;
//		System.out.println("NumMag: "+numMag);
////		MagnitudeInterpolator magInterp = new MagnitudeInterpolator(refMag, maxMag, numMag);
//		
////		FocalMechInterpolator mechInterp = new FocalMechInterpolator();
//		
////		WaldAllenGlobalVs30 vs30Provider = null;
//		WaldAllenGlobalVs30 vs30Provider = new WaldAllenGlobalVs30();
//		vs30Provider.setActiveCoefficients();
//		
//		DoubleParameterInterpolator vs30Interp = null;
//		if (vs30Provider != null)
//			vs30Interp = new DoubleParameterInterpolator(
//				Vs30_Param.NAME, 180, 760, 20); // matches Wald Allen range
//		
//		double b = 1;
//		
//		double durationYears = 30d/365d;
////		double durationYears = 1d;
//		
//		Map<FocalMech, Double> mechWts = new HashMap<>();
//		mechWts.put(FocalMech.STRIKE_SLIP, 0.5);
//		mechWts.put(FocalMech.NORMAL, 0.25);
//		mechWts.put(FocalMech.REVERSE, 0.25);
////		mechWts.put(FocalMech.STRIKE_SLIP, 0d);
////		mechWts.put(FocalMech.NORMAL, 0.5);
////		mechWts.put(FocalMech.REVERSE, 0.5);
//		List<Map<FocalMech, Double>> mechWtsList = new ArrayList<>();
//		mechWtsList.add(mechWts);
//		
//		GriddedInterpGMPE_Calc calc;
//		if (vs30Provider == null)
//			calc = new GriddedInterpGMPE_Calc(gmpe, xVals, b, refMag, maxMag, numMag, distInterp);
//		else
//			calc = new GriddedInterpGMPE_Calc(gmpe, xVals, b, refMag, maxMag, numMag, distInterp, vs30Interp);
//		double[] depths = { 7, 2 };
//		calc.precalc(durationYears, depths, mechWts);
//		
//		GeoDataSet griddedTotCumRates = ArbDiscrGeoDataSet.loadXYZFile("/home/kevin/OpenSHA/oaf/etas_tests/rateMap.txt", true);
//		
//		System.out.println(gmpe.getIntensityMeasure().getName());
//		
//		ERF erf = new ETAS_ShakingForecastCalc.GriddedForecast(griddedTotCumRates, refMag, maxMag, b, mechWtsList, depths, durationYears);
//		erf.updateForecast();
//		Preconditions.checkState(erf.getTimeSpan().getDuration() == durationYears);
//		
//		GriddedRegion calcRegion = new GriddedRegion(new Region(new Location(griddedTotCumRates.getMaxLat(), griddedTotCumRates.getMaxLon()),
//						new Location(griddedTotCumRates.getMinLat(), griddedTotCumRates.getMinLon())), calcSpacing, null);
//		List<Site> sites = new ArrayList<>();
//		List<Double> vs30s = null;
//		if (vs30Provider != null)
//			vs30s = vs30Provider.getValues(calcRegion.getNodeList());
//		for (int i=0; i<calcRegion.getNodeCount(); i++) {
//			Site site = new Site(calcRegion.locationForIndex(i));
//			for (Parameter<?> param : gmpe.getSiteParams())
//				site.addParameter((Parameter<?>) param.clone());
//			if (vs30s != null)
//				site.getParameter(Double.class, Vs30_Param.NAME).setValue(vs30s.get(i));
////			site.getParameter(Double.class, Vs30_Param.NAME).setValue(325d);
//			sites.add(site);
//		}
//		
//		System.out.println(sites.size()+" sites");
//		
//		System.out.println("Calculating interpolated");
//		Stopwatch watch = Stopwatch.createStarted();
//		DiscretizedFunc[] interpCurves = calc.calc(griddedTotCumRates, sites);
////		DiscretizedFunc[] interpCurves = calc.calc(erf.getSourceList(), sites);
//		watch.stop();
//		System.out.println("Took "+watch.elapsed(TimeUnit.SECONDS)+" seconds");
//		
//		System.out.println("Calculating traditional");
//		watch.reset();
//		watch.start();
//		DiscretizedFunc[] traditionalCurves = calc.calcTraditional(erf, sites, xVals);
//		watch.stop();
//		System.out.println("Took "+watch.elapsed(TimeUnit.SECONDS)+" seconds");
//		
//		GriddedGeoDataSet interpMap = ETAS_ShakingForecastCalc.extractMap(calcRegion, interpCurves, false, 0.1);
//		GriddedGeoDataSet traditionalMap = ETAS_ShakingForecastCalc.extractMap(calcRegion, traditionalCurves, false, 0.1);
//		
//		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, Math.max(interpMap.getMaxZ(), traditionalMap.getMaxZ()));
//		
//		XYZGraphPanel xyzGP = new XYZGraphPanel();
//		
//		XYZPlotSpec interpSpec = new XYZPlotSpec(interpMap, cpt, "Spatial Forecast", "Longitude", "Latitude",
//				gmpe.getIntensityMeasure().getName());
//		XYZPlotSpec traditionalSpec = new XYZPlotSpec(traditionalMap, cpt, "Spatial Forecast", "Longitude", "Latitude",
//				gmpe.getIntensityMeasure().getName());
//		
//		List<XYZPlotSpec> specs = new ArrayList<>();
//		specs.add(interpSpec);
//		specs.add(traditionalSpec);
//		
//		Range xRange = new Range(calcRegion.getMinGridLon()-0.5*calcRegion.getLonSpacing(),
//				calcRegion.getMaxGridLon()+0.5*calcRegion.getLonSpacing());
//		Range yRange = new Range(calcRegion.getMinGridLat()-0.5*calcRegion.getLatSpacing(),
//				calcRegion.getMaxGridLat()+0.5*calcRegion.getLatSpacing());
//		List<Range> xRanges = new ArrayList<>();
//		xRanges.add(xRange);
//		List<Range> yRanges = new ArrayList<>();
//		yRanges.add(yRange);
//		yRanges.add(yRange);
//		
//		xyzGP.drawPlot(specs, false, false, xRanges, yRanges, null);
//		
//		JFrame frame = new JFrame("");
//		frame.setContentPane(xyzGP);
//		frame.pack();
//		frame.setVisible(true);
//		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//	
//		GraphWindow gw = null;
//		
//		int curveMod = interpCurves.length/20;
//		if (curveMod < 1)
//			curveMod = 1;
//		
//		for (int i=0; i<interpCurves.length; i++) {
//			if (i % curveMod > 0)
//				continue;
//			List<DiscretizedFunc> funcs = new ArrayList<>();
//			List<PlotCurveCharacterstics> chars = new ArrayList<>();
//			
//			interpCurves[i].setName("Interpolated");
//			funcs.add(interpCurves[i]);
//			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE));
//			
//			traditionalCurves[i].setName("Traditional");
//			funcs.add(traditionalCurves[i]);
//			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
//			
//			PlotSpec spec = new PlotSpec(funcs, chars, "Site "+i, gmpe.getIntensityMeasure().getName(), "Probability");
//			spec.setLegendVisible(true);
//			
//			if (gw == null) {
//				gw = new GraphWindow(spec);
//				gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
//				gw.setVisible(true);
//			} else {
//				gw.addTab(spec);
//			}
//			gw.setXLog(true);
//			gw.setYLog(true);
//		}
//		
////		if (calc.debugMeanScatter != null) {
////			List<XY_DataSet> funcs = new ArrayList<>();
////			List<PlotCurveCharacterstics> chars = new ArrayList<>();
////			
////			funcs.add(calc.debugMeanScatter);
////			chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.BLACK));
////			
////			PlotSpec spec = new PlotSpec(funcs, chars, "Mean Scatter", "GMPE", "Interpolated");
////			gw = new GraphWindow(spec);
////			Range range = new Range(Math.min(calc.debugMeanScatter.getMinX(), calc.debugMeanScatter.getMinY()),
////					Math.max(calc.debugMeanScatter.getMaxX(), calc.debugMeanScatter.getMaxY()));
////			gw.setX_AxisRange(range);
////			gw.setY_AxisRange(range);
////			
////			funcs = new ArrayList<>();
////			chars = new ArrayList<>();
////			
////			funcs.add(calc.debugStdDevScatter);
////			chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.BLACK));
////			
////			spec = new PlotSpec(funcs, chars, "Std Dev Scatter", "GMPE", "Interpolated");
////			gw = new GraphWindow(spec);
////			range = new Range(Math.min(calc.debugStdDevScatter.getMinX(), calc.debugStdDevScatter.getMinY()),
////					Math.max(calc.debugStdDevScatter.getMaxX(), calc.debugStdDevScatter.getMaxY()));
////			gw.setX_AxisRange(range);
////			gw.setY_AxisRange(range);
////		}
//
//	}

}
