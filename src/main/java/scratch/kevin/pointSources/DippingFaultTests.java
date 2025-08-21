package scratch.kevin.pointSources;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.WeightedList;
import org.opensha.commons.data.WeightedValue;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.param.BackgroundRupType;
import org.opensha.sha.earthquake.rupForecastImpl.PointSourceNshm;
import org.opensha.sha.earthquake.util.GriddedFiniteRuptureSettings;
import org.opensha.sha.earthquake.util.GriddedSeismicitySettings;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.faultSurface.PointSurface.SiteSpecificDistanceCorrected;
import org.opensha.sha.faultSurface.RectangularSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.utils.PointSurfaceBuilder;
import org.opensha.sha.faultSurface.utils.RjbDistributionDistanceCorrection;
import org.opensha.sha.faultSurface.utils.ptSrcCorr.PointSourceDistanceCorrection;
import org.opensha.sha.faultSurface.utils.ptSrcCorr.PointSourceDistanceCorrections;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.BJF_1997_AttenRel;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.util.FocalMech;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import gov.usgs.earthquake.nshmp.gmm.Gmm;
import gov.usgs.earthquake.nshmp.gmm.GmmInput;
import gov.usgs.earthquake.nshmp.gmm.GroundMotionModel;
import gov.usgs.earthquake.nshmp.gmm.Imt;
import net.mahdilamb.colormap.Colors;

public class DippingFaultTests {
	
	enum DistPairScalars {
		RUP_ONLY("Rrup") {
			@Override
			public double getScalar(double rJB, double rRup) {
				return rRup;
			}
		},
		JB_ONLY("Rjb") {
			@Override
			public double getScalar(double rJB, double rRup) {
				return rJB;
			}
		},
		SUM("Sum") {
			@Override
			public double getScalar(double rJB, double rRup) {
				return rRup + rJB;
			}
		},
		SUM_SQ("Sum Sq") {
			@Override
			public double getScalar(double rJB, double rRup) {
				return rRup*rRup + rJB*rJB;
			}
		},
		SQ_RT_SUM_SQ("SRSS") {
			@Override
			public double getScalar(double rJB, double rRup) {
				return Math.sqrt(rRup*rRup + rJB*rJB);
			}
//		},
//		FANCY_TEST("Fancy Test Metric") {
//			@Override
//			public double getScalar(double rJB, double rRup) {
////				return SUM.getScalar(rJB, rRup) + 0.03*rJB;
////				return 1.5*rJB + rRup;
////				return rJB + 3*rRup;
////				double p = 0.5;
////				double a = 1.2;
////				return Math.pow(Math.pow(rRup, p) + Math.pow(a*rJB, p), 1d/p);
////				return (rRup-rJB)/rRup;
////				return 1 + rRup + rJB*(rRup-rJB)/rRup;
////				return SUM.getScalar(rJB, rRup) + SUM_SQ.getScalar(rJB, rRup);
//				double p = 1.5;
//				return Math.pow(Math.pow(rJB, p) + Math.pow(rRup, p), 1d/p);
//			}
//		},
//		SIMPLE_GMPE("Simple GMPE") {
//			GroundMotionModel gmm;
//			@Override
//			public double getScalar(double rJB, double rRup) {
//				if (gmm == null)
////					gmm = Gmm.BJF_97.instance(Imt.PGA);
//					gmm = Gmm.ASK_14_BASE.instance(Imt.PGA);
//				GmmInput input = GmmInput.builder().withDefaults().rX(1d).rJB(rJB).rRup(rRup).build();
//				return -gmm.calc(input).get(0).value().mean();
//			}
		};
		
		private String label;

		private DistPairScalars(String label) {
			this.label = label;
		}
		
		public abstract double getScalar(double rJB, double rRup);
	}

	public static void main(String[] args) throws IOException {
		double mag = 6.55;
//		double mag = 7.05;
		FocalMech mech = FocalMech.REVERSE;
		double[] dists = {0d, 1d, 3d, 5d, 10, 15d, 20d, 24d, 50d};
//		double[] dists = {5d};
//		double[] dists = {15d};
//		double[] dists = {50d};
		double[] periods = {0d, 0.2, 1d, 5d};
		
		double testJB = 41;
		
		boolean printPointRrups = dists.length == 1;
		
//		Boolean[] forceHWs = {null, false, true};
//		Boolean[] forceHWs = {true};
		Boolean[] forceHWs = {null};
		
		Location loc = new Location(0d, 0d);
		
		AttenRelRef gmmRef = AttenRelRef.USGS_NSHM23_ACTIVE;
		TectonicRegionType trt = TectonicRegionType.ACTIVE_SHALLOW;
		
////		PointSourceDistanceCorrection distCorr = PointSourceDistanceCorrections.FIVE_POINT_RJB_DIST.get();
////		PointSourceDistanceCorrection distCorr = PointSourceDistanceCorrections.TWENTY_POINT_RJB_DIST.get();
//		PointSourceDistanceCorrection distCorr = PointSourceDistanceCorrections.FIVE_POINT_SPINNING_DIST.get();
//		boolean randDD = false;
//		boolean randDAS = false;
//		int numRand = 1;
		
		PointSourceDistanceCorrection distCorr = PointSourceDistanceCorrections.FIVE_POINT_SPINNING_DIST_ALONG.get();
		boolean randDD = true;
		boolean randDAS = true;
		int numRand = 11;
		
		int numDD = randDD ? numRand : 1;
		int numDAS = randDAS ? numRand : 1;
		int numSubSurf = numDD * numDAS;
		int numStrikes = Integer.max(180, 5000/numSubSurf);
//		int numStrikes = 5000/numSubSurf;
		
		double[] strikes = buildSpacedSamples(0d, 360d, numStrikes, true);
		double[] dds = randDD ? buildSpacedSamples(0d, 1d, numDD, false) : new double[] {0.5};
		double[] dass = randDAS ? buildSpacedSamples(0d, 1d, numDAS, false) : new double[] {0.5};
		
		CPT subColors = null;
		if (numSubSurf > 1)
			subColors = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(0d, numSubSurf-1);
		
		DistPairScalars[] distPairScalars = DistPairScalars.values();
		
		CPT tab10 = GMT_CPT_Files.CATEGORICAL_TAB10.instance();
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		for (Boolean forceHW : forceHWs) {
			
//			double[] strikes = PointSurfaceBuilder.getEvenlySpanningStrikes(numFinite, 0d);
//			GriddedFiniteRuptureSettings finiteSettings = GriddedFiniteRuptureSettings.DEFAULT_CROSSHAIR.forStrike(0d)
//					.forNumSurfaces(numFinite).forSampleAlongStrike(randDAS).forSampleDownDip(randDD);
			boolean trackHWRanges = !randDAS && !randDD;
			
			File outputDir = new File("/tmp/dipping_fault_tests");
			if (forceHW != null) {
				if (forceHW)
					outputDir = new File(outputDir.getParentFile(), outputDir.getName()+"_hw_only");
				else
					outputDir = new File(outputDir.getParentFile(), outputDir.getName()+"_fw_only");
			}
			Preconditions.checkState(outputDir.exists() || outputDir.mkdirs(),
					"Failed to create output directory: " + outputDir.getAbsolutePath());
			
			File distPairOutDir = null;
			if (forceHW == null || forceHW) {
				distPairOutDir = new File(outputDir, "hw_dist_pair_scalars");
				Preconditions.checkState(distPairOutDir.exists() || distPairOutDir.mkdirs(),
						"Failed to create output directory: " + distPairOutDir.getAbsolutePath());
			}
			
			PointSurfaceBuilder builder = new PointSurfaceBuilder(loc)
					.random(12345l)
					.magnitude(mag);
			
			double dip = mech.dip();
			double dipRad = Math.toRadians(dip);
			double zTop = PointSourceNshm.SURF_BUILDER_DEFAULT.depthForMag(mag);
			double widthDD = PointSourceNshm.SURF_BUILDER_DEFAULT.calcWidth(mag, zTop, dipRad);
			double length = PointSourceNshm.SURF_BUILDER_DEFAULT.calcLength(mag);
			builder.upperDepthWidthAndDip(zTop, widthDD, dip);
			builder.fractionalDAS(0.5).fractionalHypocentralDepth(0.5d);
			builder.length(length);
			
			DecimalFormat oDF = new DecimalFormat("0.##");
			DecimalFormat pDF = new DecimalFormat("0.#%");
			
			PointSurface pointSurf = builder.buildPointSurface();
			
			int numFinite = numStrikes * numSubSurf;
			WeightedList<RectangularSurface> finiteSurfs = new WeightedList<>(numFinite);
			double finiteWeightEach = 1d/(double)numFinite;
			List<Double> finiteRupStrikes = new ArrayList<>(numFinite);
			List<Integer> finiteRupSubSurfIndexes = new ArrayList<>(numFinite);
			int refSubSurfIndex = -1;
			int tmpSubSurfIndex = 0;
			for (int i=0; i<numDD; i++) {
				for (int j=0; j<numDAS; j++) {
					if (dds[i] == 0.5d && dass[j] == 0.5d)
						refSubSurfIndex = tmpSubSurfIndex;
					builder.fractionalHypocentralDepth(dds[i]);
					builder.fractionalDAS(dass[j]);
					for (double strike : strikes) {
						finiteRupStrikes.add(strike);
						finiteRupSubSurfIndexes.add(tmpSubSurfIndex);
						finiteSurfs.add(builder.buildRectSurface(strike), finiteWeightEach);
					}
					tmpSubSurfIndex++;
				}
			}
			
			System.out.println("Fault dimensions: length= "+oDF.format(length)+", zTop="+oDF.format(zTop)+", dip="+oDF.format(dip)
					+", ddw="+oDF.format(widthDD)+", horzWidth="+oDF.format(widthDD*Math.cos(dipRad)));
			
			System.out.println("Built 1 point surfaces and "+finiteSurfs.size()+" finite surfaces");
			
			ScalarIMR gmm = gmmRef.get();
			gmm.setIntensityMeasure(PGA_Param.NAME);
			
			Range imRange = new Range(-5d, 2d);
			Range distRange = new Range(0d, StatUtils.max(dists)+10d);
			EvenlyDiscretizedFunc logIMfunc = new EvenlyDiscretizedFunc(imRange.getLowerBound(), imRange.getUpperBound(), 1000);
			EvenlyDiscretizedFunc distFunc = HistogramFunction.getEncompassingHistogram(distRange.getLowerBound(), distRange.getUpperBound(), 0.5d);
			
			List<List<PlotSpec>> perCombGMPlots = new ArrayList<>();
			for (int p=0; p<periods.length; p++)
				perCombGMPlots.add(new ArrayList<>());
			
			List<List<List<DefaultXY_DataSet>>> hwFiniteDistPairMedScalars = new ArrayList<>();
			for (int i=0; i<distPairScalars.length; i++) {
				List<List<DefaultXY_DataSet>> distScalars = new ArrayList<>();
				for (int d=0; d<dists.length; d++) {
					List<DefaultXY_DataSet> perScalars = new ArrayList<>();
					for (int p=0; p<periods.length; p++)
						perScalars.add(new DefaultXY_DataSet());
					distScalars.add(perScalars);
				}
				hwFiniteDistPairMedScalars.add(distScalars);
			}
			
			for (int d=0; d<dists.length; d++) {
				double dist = dists[d];
				Location siteLoc = dist == 0d ? loc : LocationUtils.location(loc, 0d, dist);
				Site site = new Site(siteLoc);
				site.addParameterList(gmm.getSiteParams());
				gmm.setSite(site);
				
				WeightedList<SiteSpecificDistanceCorrected> myPointSurfs = pointSurf.getForDistanceCorrection(siteLoc, distCorr, trt, mag);
				WeightedList<? extends RuptureSurface> myFiniteSurfs;
				List<Double> myFiniteStrikes = new ArrayList<>(numFinite);
				List<Integer> myFiniteRupSubSurfIndexes = new ArrayList<>(numFinite);
				if (forceHW == null) {
					myFiniteSurfs = finiteSurfs;
					myFiniteStrikes = finiteRupStrikes;
					myFiniteRupSubSurfIndexes = finiteRupSubSurfIndexes;
				} else {
					WeightedList<RuptureSurface> matchingFiniteSurfs = new WeightedList<>();
					for (int i=0; i<finiteSurfs.size(); i++) {
						RuptureSurface surf = finiteSurfs.getValue(i);
						boolean surfHW = surf.getDistanceX(siteLoc) >= 0d;
						if (surfHW == forceHW.booleanValue()) {
							matchingFiniteSurfs.add(surf, finiteWeightEach);
							myFiniteStrikes.add(finiteRupStrikes.get(i));
							myFiniteRupSubSurfIndexes.add(finiteRupSubSurfIndexes.get(i));
						}
					}
					if (matchingFiniteSurfs.isEmpty()) {
						System.out.println("ForceHW="+forceHW+" but no finite surfaces match at rEpi="+oDF.format(dist)+", skipping\n");
						continue;
					} else {
						matchingFiniteSurfs.normalize();
						System.out.println("ForceHW="+forceHW+", found "+matchingFiniteSurfs.size()
								+" matching finite surfaces ("+pDF.format((double)matchingFiniteSurfs.size()/(double)finiteSurfs.size())+")");
					}
					myFiniteSurfs = matchingFiniteSurfs;
					
					WeightedList<SiteSpecificDistanceCorrected> matchingPointSurfs = new WeightedList<>();
					for (WeightedValue<SiteSpecificDistanceCorrected> surfVal : myPointSurfs) {
						boolean surfHW = surfVal.value.getDistanceX(siteLoc) >= 0d;
						if (surfHW == forceHW.booleanValue())
							matchingPointSurfs.add(surfVal.value, surfVal.weight);
					}
					if (matchingPointSurfs.isEmpty()) {
						System.out.println("ForceHW="+forceHW+" but no point surfaces match at rEpi="+oDF.format(dist)+", skipping\n");
						continue;
					} else {
						matchingPointSurfs.normalize();
						System.out.println("ForceHW="+forceHW+", found "+matchingPointSurfs.size()
								+" matching point surfaces ("+pDF.format((double)matchingPointSurfs.size()/(double)myPointSurfs.size())+")");
					}
					myPointSurfs = matchingPointSurfs;
					Preconditions.checkState(myFiniteStrikes.size() == myFiniteSurfs.size());
					Preconditions.checkState(myFiniteRupSubSurfIndexes.size() == myFiniteSurfs.size());
				}
				
				List<XY_DataSet> rJBFuncs = new ArrayList<>();
				List<PlotCurveCharacterstics> rJBChars = new ArrayList<>();
				List<XY_DataSet> rRupFuncs = new ArrayList<>();
				List<PlotCurveCharacterstics> rRupChars = new ArrayList<>();
				List<XY_DataSet> rRupVsJBFuncs = new ArrayList<>();
				List<PlotCurveCharacterstics> rRupVsJBChars = new ArrayList<>();
				
				List<DiscretizedFunc> angleRjbFuncs = new ArrayList<>();
				List<DiscretizedFunc> angleRrupFuncs = new ArrayList<>();
				List<PlotCurveCharacterstics> angleChars = new ArrayList<>();
				
				List<List<XY_DataSet>> gmDistFuncs = new ArrayList<>();
				for (int p=0; p<periods.length; p++)
					gmDistFuncs.add(new ArrayList<>());
				List<PlotCurveCharacterstics> gmDistChars = new ArrayList<>();
				
				double minDist = Double.POSITIVE_INFINITY;
				double maxDist = 0d;
				double histMax = 0d;
				
				String hwAnnAdd = ", HW:";
				
				for (boolean finite : new boolean[] { true, false }) {
					WeightedList<? extends RuptureSurface> surfs = finite ? myFiniteSurfs : myPointSurfs;
					String type = finite ? "finite" : "point";

					System.out.println("Calculating " + type + " surface for rEpi=" + oDF.format(dist) + " km");
					
					EvenlyDiscretizedFunc rJBHist = new EvenlyDiscretizedFunc(distFunc.getMinX(), distFunc.size(), distFunc.getDelta());
					EvenlyDiscretizedFunc rRupHist = new EvenlyDiscretizedFunc(distFunc.getMinX(), distFunc.size(), distFunc.getDelta());
					
					List<List<double[]>> rRupVsJBs = new ArrayList<>();
					List<List<Boolean>> rRupVsJB_hws = new ArrayList<>();
					int numVsFuncs = finite ? numSubSurf : 1;
					for (int i=0; i<numVsFuncs; i++) {
						rRupVsJBs.add(new ArrayList<>());
						rRupVsJB_hws.add(new ArrayList<>());
					}
					
					MinMaxAveTracker rJBtrack = new MinMaxAveTracker();
					MinMaxAveTracker rRuptrack = new MinMaxAveTracker();
					double rJBAvg = 0d;
					double rRupAvg = 0d;
					double fractHW = 0d;
					List<Range> hwRanges = new ArrayList<>();
					double curFirstHWStrike = Double.NaN;
					double curLastHWStrike = Double.NaN;
					
					List<ArbitrarilyDiscretizedFunc> angleRjbFuncsHW = null;
					List<ArbitrarilyDiscretizedFunc> angleRrupFuncsHW = null;
					List<ArbitrarilyDiscretizedFunc> angleRjbFuncsFW = null;
					List<ArbitrarilyDiscretizedFunc> angleRrupFuncsFW = null;
					if (finite) {
						angleRjbFuncsHW = new ArrayList<>();
						angleRrupFuncsHW = new ArrayList<>();
						angleRjbFuncsFW = new ArrayList<>();
						angleRrupFuncsFW = new ArrayList<>();
						for (int i=0; i<numSubSurf; i++) {
							angleRjbFuncsHW.add(new ArbitrarilyDiscretizedFunc());
							angleRrupFuncsHW.add(new ArbitrarilyDiscretizedFunc());
							angleRjbFuncsFW.add(new ArbitrarilyDiscretizedFunc());
							angleRrupFuncsFW.add(new ArbitrarilyDiscretizedFunc());
						}
					}
					
					EvenlyDiscretizedFunc[] exceedFuncs = new EvenlyDiscretizedFunc[periods.length];
					for (int p=0; p<periods.length; p++)
						exceedFuncs[p] = new EvenlyDiscretizedFunc(logIMfunc.getMinX(), logIMfunc.size(), logIMfunc.getDelta());
					
					Preconditions.checkState(surfs.isNormalized());
					
					double closestTestJBDiff = Double.POSITIVE_INFINITY;
					double closestTestJB = Double.NaN;
					double closestRupToTestJB = Double.NaN;
					
					for (int s=0; s<surfs.size(); s++) {
						WeightedValue<? extends RuptureSurface> weightedSurf = surfs.get(s);
						RuptureSurface surf = weightedSurf.value;
						int subIndex = finite ? myFiniteRupSubSurfIndexes.get(s) : 0;
						double weight = weightedSurf.weight;
						EqkRupture rup = new EqkRupture(mag, mech.rake(), surf, null);
						gmm.setEqkRupture(rup);
						
						double rJB = surf.getDistanceJB(siteLoc);
						double rRup = surf.getDistanceRup(siteLoc);
						boolean hw = surf.getDistanceX(siteLoc) >= 0d;
						
						if (finite && hw && Double.isFinite(testJB)) {
							double diff = Math.abs(testJB - rJB);
							if (diff < closestTestJBDiff) {
								closestTestJB = rJB;
								closestTestJBDiff = diff;
								closestRupToTestJB = rRup;
							}
						}
						
						if (finite) {
							double strike = myFiniteStrikes.get(s);
							
							ArbitrarilyDiscretizedFunc angleRjbFuncHW = angleRjbFuncsHW.get(subIndex);
							ArbitrarilyDiscretizedFunc angleRrupFuncHW = angleRrupFuncsHW.get(subIndex);
							ArbitrarilyDiscretizedFunc angleRjbFuncFW = angleRjbFuncsFW.get(subIndex);
							ArbitrarilyDiscretizedFunc angleRrupFuncFW = angleRrupFuncsFW.get(subIndex);
							
							if (hw) {
								angleRjbFuncHW.set(strike, rJB);
								angleRrupFuncHW.set(strike, rRup);
								angleRjbFuncFW.set(strike, Double.NaN);
								angleRrupFuncFW.set(strike, Double.NaN);
							} else {
								angleRjbFuncHW.set(strike, Double.NaN);
								angleRrupFuncHW.set(strike, Double.NaN);
								angleRjbFuncFW.set(strike, rJB);
								angleRrupFuncFW.set(strike, rRup);
							}
						}
						
						if (!finite && printPointRrups && hw)
							System.out.println("\tEstimated hanging-wall rRup="+(float)rRup+" for rJB="+(float)rJB);
						
						rRupVsJBs.get(subIndex).add(new double[] {rJB, rRup});
						rRupVsJB_hws.get(subIndex).add(hw);
						
						minDist = Math.min(minDist, Math.min(rJB, rRup));
						maxDist = Math.max(maxDist, Math.max(rJB, rRup));
						
						rJBHist.add(distFunc.getClosestXIndex(rJB), weight);
						rRupHist.add(distFunc.getClosestXIndex(rRup), weight);
						rJBtrack.addValue(rJB);
						rRuptrack.addValue(rRup);
						rJBAvg += rJB * weight;
						rRupAvg += rRup * weight;
						if (hw) {
							fractHW += weight;
							if (trackHWRanges && finite) {
								// finite
								double strike = myFiniteStrikes.get(s);
								if (Double.isNaN(curFirstHWStrike))
									curFirstHWStrike = strike;
								curLastHWStrike = strike;
							}
						} else {
							if (Double.isFinite(curFirstHWStrike))
								hwRanges.add(new Range(curFirstHWStrike, curLastHWStrike));
							curFirstHWStrike = Double.NaN;
							curLastHWStrike = Double.NaN;
						}
						
						double[] pairScalars = null;
						if (finite) {
							pairScalars = new double[distPairScalars.length];
							for (int i=0; i<distPairScalars.length; i++)
								pairScalars[i] = distPairScalars[i].getScalar(rJB, rRup);
						}

						for (int p=0; p<periods.length; p++) {
							double period = periods[p];
							if (period == 0d) {
								gmm.setIntensityMeasure(PGA_Param.NAME);
							} else {
								gmm.setIntensityMeasure(SA_Param.NAME);
								SA_Param.setPeriodInSA_Param(gmm.getIntensityMeasure(), period);
							}
							
							gmm.getExceedProbabilities(logIMfunc);
							
							for (int i=0; i<logIMfunc.size(); i++)
								exceedFuncs[p].add(i, weight*logIMfunc.getY(i));
							
							if (finite && hw) {
								double median = Math.exp(gmm.getMean());
								for (int i=0; i<pairScalars.length; i++)
									hwFiniteDistPairMedScalars.get(i).get(d).get(p).set(pairScalars[i], median);
							}
						}
					}
					System.out.println("\trJB stats:\t"+rJBtrack);
					System.out.println("\trJB weighted avg:\t"+(float)rJBAvg);
					System.out.println("\trRup stats:\t"+rRuptrack);
					System.out.println("\trRup weighted avg:\t"+(float)rRupAvg);
					System.out.println("\tFract HW:\t"+(float)fractHW);
					
					if (finite && Double.isFinite(testJB) && closestTestJBDiff < 1) {
						System.out.println("Found a close finite FW case to rJB="+(float)testJB);
						System.out.println("\tcloseJB="+(float)closestTestJB);
						System.out.println("\trRup="+(float)closestRupToTestJB);
						System.out.println("Doing pt src equiv");
						double ptRrup = RjbDistributionDistanceCorrection.getCorrDistRup(closestTestJB, dist, zTop,
								pointSurf.getAveRupBottomDepth(), dipRad, length, pointSurf.getAveHorizontalWidth(), false);
						System.out.println("\tPt rRup="+(float)ptRrup);
					}
					
					if (finite)
						hwAnnAdd += " F=";
					else
						hwAnnAdd += " Pt=";
					hwAnnAdd += pDF.format(fractHW);
					
					if (Double.isFinite(curFirstHWStrike))
						hwRanges.add(new Range(curFirstHWStrike, curLastHWStrike));
					if (!hwRanges.isEmpty()) {
						System.out.print("\tHW ranges:");
						for (Range hwRange : hwRanges)
							System.out.print("\t["+oDF.format(hwRange.getLowerBound())+", "+oDF.format(hwRange.getUpperBound())+"]");
						System.out.println();
					}
					
					DefaultXY_DataSet meanRjbLine = new DefaultXY_DataSet();
					meanRjbLine.set(rJBAvg, 0d);
					meanRjbLine.set(rJBAvg, 1d);
					
					DefaultXY_DataSet meanRrupLine = new DefaultXY_DataSet();
					meanRrupLine.set(rRupAvg, 0d);
					meanRrupLine.set(rRupAvg, 1d);
					
					Color histColor;
					Color lineColor;
					PlotSymbol scatterSymbol;
					PlotSymbol thinScatterSymbol;
					Color scatterColorFW;
					Color scatterColorHW;
					PlotLineType gmLineType;
					String name;
					if (finite) {
						name = "Finite";
//						histColor = Colors.tab_lightorange;
						lineColor = Colors.tab_orange;
						histColor = new Color(lineColor.getRed(), lineColor.getGreen(), lineColor.getBlue(), 127);
						scatterColorHW = Colors.tab_orange;
						scatterColorFW = Colors.tab_lightorange;
						scatterSymbol = PlotSymbol.BOLD_CROSS;
						thinScatterSymbol = PlotSymbol.CROSS;
						gmLineType = PlotLineType.SOLID;
					} else {
						name = "Point Sources";
//						histColor = Colors.tab_lightblue;
						lineColor = Colors.tab_blue;
						histColor = new Color(lineColor.getRed(), lineColor.getGreen(), lineColor.getBlue(), 127);
						scatterColorHW = Colors.tab_blue;
						scatterColorFW = Colors.tab_lightblue;
						scatterSymbol = PlotSymbol.BOLD_X;
						thinScatterSymbol = PlotSymbol.X;
						gmLineType = PlotLineType.SHORT_DASHED;
					}
					rJBHist.setName(name+" Distribution");
					rRupHist.setName(name+" Distribution");
					meanRjbLine.setName(name+" Mean");
					meanRrupLine.setName(name+" Mean");
					exceedFuncs[0].setName(name);
					
					rJBFuncs.add(rJBHist);
					rJBChars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, histColor));
					rJBFuncs.add(meanRjbLine);
					rJBChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, lineColor));
					
					rRupFuncs.add(rRupHist);
					rRupChars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, histColor));
					rRupFuncs.add(meanRrupLine);
					rRupChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, lineColor));
					
					if (finite) {
						// re-sort to put the centered one last
						if (numSubSurf > 1) {
							rRupVsJBs.add(rRupVsJBs.remove(refSubSurfIndex));
							rRupVsJB_hws.add(rRupVsJB_hws.remove(refSubSurfIndex));
						}
						
						for (boolean hw : new boolean[] {false,true}) {
							if (forceHW != null && forceHW != hw)
								continue;
							Color color = hw ? scatterColorHW : scatterColorFW;
							for (int i=0; i<numSubSurf; i++) {
								PlotCurveCharacterstics fChar;
								String fName;
								if (i == numSubSurf-1) {
									// centered
									fName = "Centered "+name;
									if (hw)
										fName += " (HW)";
									else
										fName += " (FW)";
									fChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 5f, color);
								} else {
									if (i == 0) {
										fName = "Shifted "+name;
										if (hw)
											fName += " (HW)";
										else
											fName += " (FW)";
									} else {
										fName = null;
									}
									Color myColor = color;
									if (subColors != null)
										myColor = subColors.getColor(i);
									fChar = new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 2f, myColor);
								}
								List<DefaultXY_DataSet> funcs = getSortedVsFuncs(rRupVsJBs.get(i), rRupVsJB_hws.get(i), hw);
								
								for (int f=0; f<funcs.size(); f++) {
									DefaultXY_DataSet func = funcs.get(f);
									if (f == 0)
										func.setName(fName);
									rRupVsJBFuncs.add(func);
									rRupVsJBChars.add(fChar);
								}
							}
						}
					} else {
						List<DefaultXY_DataSet> rRupVsJB_fw = getSortedVsFuncs(rRupVsJBs.get(0), rRupVsJB_hws.get(0), false);
						List<DefaultXY_DataSet> rRupVsJB_hw = getSortedVsFuncs(rRupVsJBs.get(0), rRupVsJB_hws.get(0), true);
						
						if (forceHW == null || !forceHW) {
							for (int f=0; f<rRupVsJB_fw.size(); f++) {
								DefaultXY_DataSet func = rRupVsJB_fw.get(f);
								if (f == 0)
									func.setName(func.size()+" "+name+" (FW)");
								rRupVsJBFuncs.add(func);
								rRupVsJBChars.add(new PlotCurveCharacterstics(scatterSymbol, 5f, scatterColorFW));
							}
						}
						if (forceHW == null || forceHW) {
							for (int f=0; f<rRupVsJB_hw.size(); f++) {
								DefaultXY_DataSet func = rRupVsJB_hw.get(f);
								if (f == 0)
									func.setName(func.size()+" "+name+" (FW)");
								rRupVsJBFuncs.add(func);
								rRupVsJBChars.add(new PlotCurveCharacterstics(scatterSymbol, 5f, scatterColorHW));
							}
						}
					}
					
					if (finite) {
						if (numSubSurf > 1) {
							for (int s=0; s<numSubSurf; s++) {
								ArbitrarilyDiscretizedFunc angleRjbFuncHW = angleRjbFuncsHW.get(s);
								ArbitrarilyDiscretizedFunc angleRrupFuncHW = angleRrupFuncsHW.get(s);
								ArbitrarilyDiscretizedFunc angleRjbFuncFW = angleRjbFuncsFW.get(s);
								ArbitrarilyDiscretizedFunc angleRrupFuncFW = angleRrupFuncsFW.get(s);
								
								if (forceHW == null || !forceHW) {
									angleRjbFuncs.add(angleRjbFuncFW);
									angleRrupFuncs.add(angleRrupFuncFW);
									angleChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, scatterColorFW));
								}
								if (forceHW == null || forceHW) {
									angleRjbFuncs.add(angleRjbFuncHW);
									angleRrupFuncs.add(angleRrupFuncHW);
									angleChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, scatterColorHW));
								}
							}
						}
						ArbitrarilyDiscretizedFunc angleRjbFuncHW = angleRjbFuncsHW.get(refSubSurfIndex);
						ArbitrarilyDiscretizedFunc angleRrupFuncHW = angleRrupFuncsHW.get(refSubSurfIndex);
						ArbitrarilyDiscretizedFunc angleRjbFuncFW = angleRjbFuncsFW.get(refSubSurfIndex);
						ArbitrarilyDiscretizedFunc angleRrupFuncFW = angleRrupFuncsFW.get(refSubSurfIndex);
						if (numSubSurf > 1) {
							// clone to set the name only here
							angleRjbFuncHW = angleRjbFuncHW.deepClone();
							angleRrupFuncHW = angleRrupFuncHW.deepClone();
							angleRjbFuncFW = angleRjbFuncFW.deepClone();
							angleRrupFuncFW = angleRrupFuncFW.deepClone();
						}
						angleRjbFuncFW.setName("Footwall");
						angleRrupFuncFW.setName("Footwall");
						angleRjbFuncHW.setName("Hanging wall");
						angleRrupFuncHW.setName("Hanging wall");
						
						if (forceHW == null || !forceHW) {
							angleRjbFuncs.add(angleRjbFuncFW);
							angleRrupFuncs.add(angleRrupFuncFW);
							angleChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, scatterColorFW));
						}
						if (forceHW == null || forceHW) {
							angleRjbFuncs.add(angleRjbFuncHW);
							angleRrupFuncs.add(angleRrupFuncHW);
							angleChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, scatterColorHW));
						}
					}
					
					for (int p=0; p<periods.length; p++)
						gmDistFuncs.get(p).add(exceedFuncs[p]);
					gmDistChars.add(new PlotCurveCharacterstics(gmLineType, 3f, lineColor));
					
					histMax = Math.max(histMax, Math.max(rJBHist.getMaxY(), rRupHist.getMaxY()));
				}
				System.out.println();
				
				double distBuffer = 0.1*(maxDist - minDist);
				Range plotDistRange = new Range(Math.max(0d, minDist-distBuffer), Math.max(1d, maxDist+distBuffer));
				Range plotHistYRange = new Range(0d, Math.min(1d, histMax+0.1));
				Font annFont = new Font(Font.SANS_SERIF, Font.BOLD, 20);
				
				PlotSpec rJBplot = new PlotSpec(rJBFuncs, rJBChars, null, "Distance (km)", "Fraction");
				rJBplot.setLegendVisible(true);
				XYTextAnnotation rJBAnn = new XYTextAnnotation("rJB"+hwAnnAdd,
						plotDistRange.getLowerBound()+0.9*plotDistRange.getLength(), 0.8*plotHistYRange.getUpperBound());
				rJBAnn.setFont(annFont);
				rJBAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
				rJBplot.addPlotAnnotation(rJBAnn);
				
				PlotSpec rRupPlot = new PlotSpec(rRupFuncs, rRupChars, null, "Distance (km)", "Fraction");
				rRupPlot.setLegendVisible(false);
				XYTextAnnotation rRupAnn = new XYTextAnnotation("rRup"+hwAnnAdd, rJBAnn.getX(), rJBAnn.getY());
				rRupAnn.setFont(annFont);
				rRupAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
				rRupPlot.addPlotAnnotation(rRupAnn);
				
				gp.drawGraphPanel(List.of(rJBplot, rRupPlot), false, false, List.of(plotDistRange), List.of(plotHistYRange, plotHistYRange));
				
				String distPrefix = oDF.format(dist)+"km";
				PlotUtils.writePlots(outputDir, distPrefix+"_distances", gp, 800, 800, true, false, false);
				
				DefaultXY_DataSet oneToOne = new DefaultXY_DataSet();
				oneToOne.set(plotDistRange.getLowerBound(), plotDistRange.getLowerBound());
				oneToOne.set(plotDistRange.getUpperBound(), plotDistRange.getUpperBound());
				
				rRupVsJBFuncs.add(0, oneToOne);
				rRupVsJBChars.add(0, new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.GRAY));
				
				PlotSpec rRupVsJBPlot = new PlotSpec(rRupVsJBFuncs, rRupVsJBChars, null, "Distance JB (km)", "Distance Rup (km)");
				rRupVsJBPlot.setLegendInset(RectangleAnchor.BOTTOM_RIGHT);
				
				gp.drawGraphPanel(rRupVsJBPlot, false, false, plotDistRange, plotDistRange);
				
				PlotUtils.writePlots(outputDir, distPrefix+"_distance_ratios", gp, 800, false, true, false, false);
				
				List<PlotSpec> gmDistPlots = new ArrayList<>();
				List<Range> gmDistYRanges = new ArrayList<>();
				for (int p=0; p<periods.length; p++) {
					PlotSpec plot = new PlotSpec(gmDistFuncs.get(p), gmDistChars, " ", "Ln(IML)", "P(> IML)");
					plot.setLegendVisible(p == 0);
					String text = periods[p] == 0d ? "PGA" : oDF.format(periods[p])+"s SA";
					XYTextAnnotation periodAnn = new XYTextAnnotation(text+hwAnnAdd, imRange.getLowerBound()+imRange.getLength()*0.9, 0.9);
					periodAnn.setFont(annFont);
					periodAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
					plot.addPlotAnnotation(periodAnn);
					gmDistPlots.add(plot);
					gmDistYRanges.add(new Range(0d, 1d));
					
					perCombGMPlots.get(p).add(plot);
				}
				
				gp.drawGraphPanel(gmDistPlots, false, false, List.of(imRange), gmDistYRanges);
				
				PlotUtils.writePlots(outputDir, distPrefix+"_gms", gp, 800, 800, true, false, false);
				
				PlotSpec anglePlotJB = new PlotSpec(angleRjbFuncs, angleChars, " ", "Strike Angle (degrees)", "Distance JB (km)");
				anglePlotJB.setLegendVisible(true);
				PlotSpec anglePlotRup = new PlotSpec(angleRrupFuncs, angleChars, " ", "Strike Angle (degrees)", "Distance Rup (km)");
				
				gp.drawGraphPanel(List.of(anglePlotJB, anglePlotRup), false, false, List.of(new Range(0d, 360d)), List.of(plotDistRange, plotDistRange));
				
				PlotUtils.writePlots(outputDir, distPrefix+"_strikes", gp, 800, 800, true, false, false);
				
				// now distance pair scalars
				for (int i=0; distPairOutDir != null && i<distPairScalars.length; i++) {
					List<DefaultXY_DataSet> funcs = hwFiniteDistPairMedScalars.get(i).get(d);
					List<PlotCurveCharacterstics> chars = new ArrayList<>();
					Preconditions.checkState(funcs.size() == periods.length);
					double minScalar = Double.POSITIVE_INFINITY;
					double maxScalar = Double.NEGATIVE_INFINITY;
					double minGM = Double.POSITIVE_INFINITY;
					double maxGM = Double.NEGATIVE_INFINITY;
					for (int p=0; p<periods.length; p++) {
						Color color = tab10.get(p % tab10.size()).minColor;
						chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 2f, color));
						String text = periods[p] == 0d ? "PGA" : oDF.format(periods[p])+"s SA";;
						funcs.get(p).setName(text);
						minScalar = Math.min(minScalar, funcs.get(p).getMinX());
						maxScalar = Math.max(maxScalar, funcs.get(p).getMaxX());
						minGM = Math.min(minGM, funcs.get(p).getMinY());
						maxGM = Math.max(maxGM, funcs.get(p).getMaxY());
					}
					
//					minGM = Math.pow(10, Math.floor(Math.log10(minGM)));
//					maxGM = Math.pow(10, Math.ceil(Math.log10(maxGM)));
					Range yRange = log10Range(minGM, maxGM);
					
					PlotSpec plot = new PlotSpec(funcs, chars, " ", distPairScalars[i].label, "Median GM (g)");
					plot.setLegendInset(RectangleAnchor.TOP_RIGHT);
					
					gp.drawGraphPanel(plot, false, true, null, yRange);
					
					PlotUtils.writePlots(distPairOutDir, distPrefix+"_"+distPairScalars[i].name(), gp, 800, 800, true, false, false);
					
					gp.drawGraphPanel(plot, true, true, log10Range(minScalar, maxScalar), yRange);
					
					PlotUtils.writePlots(distPairOutDir, distPrefix+"_"+distPairScalars[i].name()+"_log", gp, 800, 800, true, false, false);
				}
			}
			// now distance pair scalars by period with all distances on one plot
			for (int i=0; distPairOutDir != null && i<distPairScalars.length; i++) {
				List<DefaultXY_DataSet> allFuncs = new ArrayList<>();
				List<PlotCurveCharacterstics> allChars = new ArrayList<>();
				double allMinScalar = Double.POSITIVE_INFINITY;
				double allMaxScalar = Double.NEGATIVE_INFINITY;
				double allMinGM = Double.POSITIVE_INFINITY;
				double allMaxGM = Double.NEGATIVE_INFINITY;
				for (int p=0; p<periods.length; p++) {
					List<DefaultXY_DataSet> funcs = new ArrayList<>();
					List<PlotCurveCharacterstics> chars = new ArrayList<>();
					double minScalar = Double.POSITIVE_INFINITY;
					double maxScalar = Double.NEGATIVE_INFINITY;
					double minGM = Double.POSITIVE_INFINITY;
					double maxGM = Double.NEGATIVE_INFINITY;
					String prefix = distPairScalars[i].name()+"_";
					String yLabel = "Median ";
					if (periods[p] == 0d) {
						prefix += "pga";
						yLabel += "PGA (g)";
					} else {
						prefix += oDF.format(periods[p])+"s";
						yLabel += oDF.format(periods[p])+"s SA (g)";
					}
					for (int d=0; d<dists.length; d++) {
						Color color = tab10.get(d % tab10.size()).minColor;
						chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 2f, color));
						String text = oDF.format(dists[d])+" km";
						DefaultXY_DataSet func = hwFiniteDistPairMedScalars.get(i).get(d).get(p);
						func.setName(text);
						funcs.add(func);
						minScalar = Math.min(minScalar, func.getMinX());
						maxScalar = Math.max(maxScalar, func.getMaxX());
						minGM = Math.min(minGM, func.getMinY());
						maxGM = Math.max(maxGM, func.getMaxY());
					}
					
//					minGM = Math.pow(10, Math.floor(Math.log10(minGM)));
//					maxGM = Math.pow(10, Math.ceil(Math.log10(maxGM)));
					Range yRange = log10Range(minGM, maxGM);
					Range logXRange = log10Range(minScalar, maxScalar);
					
					allMinScalar = Math.min(logXRange.getLowerBound(), allMinScalar);
					allMaxScalar = Math.max(logXRange.getUpperBound(), allMaxScalar);
					allMinGM = Math.min(yRange.getLowerBound(), allMinGM);
					allMaxGM = Math.max(yRange.getUpperBound(), allMaxGM);
					allFuncs.addAll(funcs);
					allChars.addAll(chars);
					
					PlotSpec plot = new PlotSpec(funcs, chars, " ", distPairScalars[i].label, yLabel);
					plot.setLegendInset(RectangleAnchor.TOP_RIGHT);
					
					gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
					
					gp.drawGraphPanel(plot, false, true, null, yRange);
					
					PlotUtils.writePlots(distPairOutDir, prefix, gp, 800, 800, true, false, false);
					
					gp.drawGraphPanel(plot, true, true, logXRange, yRange);
					
					PlotUtils.writePlots(distPairOutDir, prefix+"_log", gp, 800, 800, true, false, false);
					
					if (p > 0)
						for (DefaultXY_DataSet func : funcs)
							func.setName(null);
				}
				PlotSpec plot = new PlotSpec(allFuncs, allChars, " ", distPairScalars[i].label, "Median GM (g)");
				plot.setLegendInset(RectangleAnchor.TOP_RIGHT);
				
				gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
				
				gp.drawGraphPanel(plot, false, true, null, new Range(allMinGM, allMaxGM));
				
				PlotUtils.writePlots(distPairOutDir, "all_dists_periods_"+distPairScalars[i].name(), gp, 800, 1200, true, false, false);
				
				gp.drawGraphPanel(plot, true, true, new Range(allMinScalar, allMaxScalar), new Range(allMinGM, allMaxGM));
				
				PlotUtils.writePlots(distPairOutDir, "all_dists_periods_log_"+distPairScalars[i].name(), gp, 800, 1200, true, false, false);
				
				// now by rank
				List<DefaultXY_DataSet> rankFuncs = new ArrayList<>();
				for (DefaultXY_DataSet func : allFuncs) {
					List<Point2D> points = new ArrayList<>(func.size());
					for (int n=0; n<func.size(); n++)
						points.add(func.get(n));
					points.sort(new Comparator<Point2D>() {

						@Override
						public int compare(Point2D o1, Point2D o2) {
							return Double.compare(o1.getX(), o2.getX());
						}
					});
					double[] newX = new double[points.size()];
					double[] newY = new double[points.size()];
					for (int n=0; n<points.size(); n++) {
						newX[n] = n+1;
						newY[n] = points.get(n).getY();
					}
					DefaultXY_DataSet rankFunc = new DefaultXY_DataSet(newX, newY);
					rankFunc.setName(func.getName());
					rankFuncs.add(rankFunc);
				}
				
				plot = new PlotSpec(rankFuncs, allChars, " ", distPairScalars[i].label+" Rank", "Median GM (g)");
				plot.setLegendInset(RectangleAnchor.TOP_RIGHT);
				
				gp.drawGraphPanel(plot, false, true, null, new Range(allMinGM, allMaxGM));
				
				PlotUtils.writePlots(distPairOutDir, "all_dists_periods_rank_"+distPairScalars[i].name(), gp, 800, 1200, true, false, false);
			}
		}
	}
	
	static Range log10Range(double min, double max) {
		double logMin = Math.log10(min <= 0 ? 0.1 : min);
		double logMax = Math.log10(max <= 1d ? 1d : max);
		
		double logDelta = 0.1*Math.max(1d, logMax - logMin);
		return new Range(Math.pow(10, logMin-logDelta), Math.pow(10, logMax+logDelta));
	}
	
	static double[] buildSpacedSamples(double min, double max, int num, boolean sampleEdges) {
		double delta = (max-min)/(double)(sampleEdges ? num-1 : num);
		double ret0 = sampleEdges ? min : min + 0.5*delta;
		double[] ret = new double[num];
		for (int i=0; i<num; i++)
			ret[i] = ret0 + i*delta;
		return ret;
	}
	
	static List<DefaultXY_DataSet> getSortedVsFuncs(List<double[]> data, List<Boolean> hws, boolean hw) {
		List<DefaultXY_DataSet> funcs = new ArrayList<>();
		DefaultXY_DataSet curFunc = null;
		for (int i=0; i<data.size(); i++) {
			double[] pt = data.get(i);
			boolean isHW = hws.get(i);
			if (isHW == hw) {
				if (curFunc == null) {
					curFunc = new DefaultXY_DataSet();
					funcs.add(curFunc);
				}
				curFunc.set(pt[0], pt[1]);
			} else {
				if (curFunc != null)
					curFunc = null;
			}
		}
		return funcs;
	}

}
