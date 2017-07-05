package scratch.kevin.cybershake;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.text.DecimalFormat;
import java.util.List;

import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.jfree.data.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.siteData.OrderedSiteDataProviderList;
import org.opensha.commons.data.siteData.SiteDataValue;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.param.Parameter;
import org.opensha.sha.cybershake.calc.HazardCurveComputation;
import org.opensha.sha.cybershake.db.CybershakeIM;
import org.opensha.sha.cybershake.db.CybershakeSite;
import org.opensha.sha.cybershake.db.Cybershake_OpenSHA_DBApplication;
import org.opensha.sha.cybershake.db.DBAccess;
import org.opensha.sha.cybershake.db.MeanUCERF2_ToDB;
import org.opensha.sha.cybershake.db.PeakAmplitudesFromDB;
import org.opensha.sha.cybershake.db.Runs2DB;
import org.opensha.sha.cybershake.db.CybershakeIM.CyberShakeComponent;
import org.opensha.sha.cybershake.db.CybershakeIM.IMType;
import org.opensha.sha.cybershake.plot.HazardCurvePlotCharacteristics;
import org.opensha.sha.cybershake.plot.HazardCurvePlotter;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.util.SiteTranslator;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public class GMPE_CS_DistributionComparisonPlotter {
	
	private static final DecimalFormat df = new DecimalFormat("0.###");
	
	public static void main(String[] args) {
		int sourceID = 247;
		int rupID = -1; // -1 means all
		
		List<Double> periods = Lists.newArrayList(0.2, 0.5, 2d, 3d, 5d);
		CyberShakeComponent comp = CyberShakeComponent.RotD50;
		
		int datasetID = 61;
		String csSiteName = "PAS";
		
		File outputDir = new File("/home/kevin/CyberShake/gmpe_comparison_distribution");
		
		int csHistEstNum = 15;
		boolean logX = true;
		
		List<Color> gmpeColors = HazardCurvePlotCharacteristics.getDefaultColors();
		
		DBAccess db = null;
		try {
			db = Cybershake_OpenSHA_DBApplication.getDB();
			ERF erf = MeanUCERF2_ToDB.createUCERF2ERF();
			ProbEqkSource source = erf.getSource(sourceID);
			List<Integer> rupIDs = Lists.newArrayList();
			if (rupID < 0)
				for (int r=0; r<source.getNumRuptures(); r++)
					rupIDs.add(r);
			else
				rupIDs.add(rupID);
			
			List<ScalarIMR> gmpes = Lists.newArrayList();
			gmpes.add(AttenRelRef.ASK_2014.instance(null));
			gmpes.add(AttenRelRef.BSSA_2014.instance(null));
			gmpes.add(AttenRelRef.CB_2014.instance(null));
			gmpes.add(AttenRelRef.CY_2014.instance(null));
			
			PeakAmplitudesFromDB amps2db = new PeakAmplitudesFromDB(db);
			List<CybershakeIM> ims = amps2db.getIMs(periods, IMType.SA, comp);
			Runs2DB runs2db = new Runs2DB(db);
			
			RunIDFetcher fetch = new RunIDFetcher(db, datasetID, ims.get(0).getID());
			
			int runID = fetch.getRunID(csSiteName);
			CybershakeSite csSite = fetch.getSite(csSiteName);
			System.out.println("CS Site: "+csSite+", run "+runID);
			Location loc = csSite.createLocation();
			
			int velModelID = runs2db.getRun(runID).getVelModelID();
			Site site = new Site(loc);
			for (ScalarIMR gmpe : gmpes) {
				gmpe.setParamDefaults();
				for (Parameter<?> param : gmpe.getSiteParams())
					if (!site.containsParameter(param))
						site.addParameter((Parameter<?>)param.clone());
				gmpe.setIntensityMeasure(SA_Param.NAME);
			}
			OrderedSiteDataProviderList provs = HazardCurvePlotter.createProviders(velModelID);
			List<SiteDataValue<?>> datas = provs.getAllAvailableData(site.getLocation());
			SiteTranslator trans = new SiteTranslator();
			for (Parameter<?> param : site)
				trans.setParameterValue(param, datas);
			
			for (int p=0; p<periods.size(); p++) {
				double period = periods.get(p);
				String periodStr;
				if (period == (int)period)
					periodStr = (int)period+"";
				else
					periodStr = (float)period+"";
				File subDir = new File(outputDir, source.getName().replaceAll(" ", "_")+"_"+periodStr+"s");
				Preconditions.checkState(subDir.exists() || subDir.mkdir());
				
				CybershakeIM im = ims.get(p);
				
				for (int myRupID : rupIDs) {
					ProbEqkRupture rup = source.getRupture(myRupID);
					List<Double> csAmps = amps2db.getIM_Values(runID, sourceID, myRupID, im);
					System.out.println("Loaded "+csAmps.size()+" amps for s="+sourceID+", r="+myRupID+", p="+periodStr);
					
					ArbDiscrEmpiricalDistFunc csFunc = new ArbDiscrEmpiricalDistFunc();
					for (double val : csAmps)
						csFunc.set(val/HazardCurveComputation.CONVERSION_TO_G,1);
					
					HistogramFunction csHist;
					if (logX) {
						double csHistDelta = (Math.log(csFunc.getMaxX()) - Math.log(csFunc.getMinX())) / csHistEstNum;
						csHist = HistogramFunction.getEncompassingHistogram(
								Math.log(csFunc.getMinX()), Math.log(csFunc.getMaxX()), csHistDelta);
					} else {
						double csHistDelta = (csFunc.getMaxX() - csFunc.getMinX()) / csHistEstNum;
						csHist = HistogramFunction.getEncompassingHistogram(
								csFunc.getMinX(), csFunc.getMaxX(), csHistDelta);
					}
					for (Point2D pt : csFunc)
						if (logX)
							csHist.add(Math.log(pt.getX()), pt.getY());
						else
							csHist.add(pt.getX(), pt.getY());
//					// convert to density
//					csHist.scale(1d/(csHistDelta*csHist.calcSumOfY_Vals()));
					csHist.normalizeBySumOfY_Vals();
					Range yRange = new Range(0, csHist.getMaxY()*1.2);
					
					Range xRange = new Range(csHist.getMinX() - 1.5d, csHist.getMaxX() + 1.5d);
					if (!logX)
						xRange = new Range(0d, Math.exp(xRange.getUpperBound()));
					
					double csMean = NorthridgeCSComparison.logAverage(csAmps)/HazardCurveComputation.CONVERSION_TO_G;
					XY_DataSet csMeanXY = new DefaultXY_DataSet();
					if (logX) {
						csMeanXY.set(Math.log(csMean), 0d);
						csMeanXY.set(Math.log(csMean), yRange.getUpperBound());
					} else {
						csMeanXY.set(csMean, 0d);
						csMeanXY.set(csMean, yRange.getUpperBound());
					}
					
					List<Double> gmpeMeanVals = Lists.newArrayList();
					List<DiscretizedFunc> gmpeDists = Lists.newArrayList();
					for (ScalarIMR gmpe : gmpes) {
						SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), period);
						gmpe.setSite(site);
						gmpe.setEqkRupture(rup);
						double mean = Math.exp(gmpe.getMean());
						double sd = gmpe.getStdDev();
						double shape = sd / mean;
						double scale = Math.log(mean);
						gmpeMeanVals.add(mean);
//						LogNormalDistribution dist = new LogNormalDistribution(scale, shape);
						NormalDistribution normDist = new NormalDistribution(Math.log(mean), sd);
						DiscretizedFunc distFunc = new EvenlyDiscretizedFunc(
								xRange.getLowerBound(), xRange.getUpperBound(), 1000);
						for (int i=0; i<distFunc.size(); i++) {
							if (logX)
								distFunc.set(i, normDist.density(distFunc.getX(i)));
							else
								distFunc.set(i, normDist.density(Math.log(distFunc.getX(i))));
						}
						distFunc.scale(csHist.getMaxY()/distFunc.getMaxY());
//						if (logX) {
//							DiscretizedFunc logDistFunc = new ArbitrarilyDiscretizedFunc();
//							for (Point2D pt : distFunc)
//								if (pt.getX() > 0)
//									logDistFunc.set(Math.log(pt.getX()), pt.getY());
//							distFunc = logDistFunc;
//						}
						gmpeDists.add(distFunc);
						System.out.println(gmpe.getShortName()+" mean="+mean+", sd="+sd);
					}
					
					List<XY_DataSet> funcs = Lists.newArrayList();
					List<PlotCurveCharacterstics> chars = Lists.newArrayList();
					
					funcs.add(csHist);
					chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.DARK_GRAY));
					csHist.setName("CS Distribution");
					
					funcs.add(csMeanXY);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
					csMeanXY.setName("CS Log-Mean: "+df.format(csMean));
					
					double maxGMPEDist = 0;
					
					for (int i=0; i<gmpes.size(); i++) {
						double gmpeMean = gmpeMeanVals.get(i);
						DiscretizedFunc gmpeDist = gmpeDists.get(i);
						String gmpeName = gmpes.get(i).getShortName();
						
						Color c = gmpeColors.get(i % gmpeColors.size());
						
						XY_DataSet gmpeMeanXY = new DefaultXY_DataSet();
						if (logX) {
							gmpeMeanXY.set(Math.log(gmpeMean), 0d);
							gmpeMeanXY.set(Math.log(gmpeMean), yRange.getUpperBound());
						} else {
							gmpeMeanXY.set(gmpeMean, 0d);
							gmpeMeanXY.set(gmpeMean, yRange.getUpperBound());
						}
						
						funcs.add(gmpeMeanXY);
						chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, c));
						gmpeMeanXY.setName(gmpeName+": "+df.format(gmpeMean));
						
						funcs.add(gmpeDist);
						chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, c));
						maxGMPEDist = Math.max(maxGMPEDist, gmpeDist.getMaxX());
					}
					
					String title = source.getName()+" M"+df.format(rup.getMag())+", "+csSiteName;
					
					String xAxisLabel = periodStr+"s SA";
					if (logX)
						xAxisLabel = "Ln("+xAxisLabel+")";
					
					PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisLabel, "");
					spec.setLegendVisible(true);
					
					HeadlessGraphPanel gp = new HeadlessGraphPanel();
					gp.setTickLabelFontSize(18);
					gp.setAxisLabelFontSize(20);
					gp.setPlotLabelFontSize(21);
					gp.setBackgroundColor(Color.WHITE);
					
					gp.drawGraphPanel(spec, false, false, xRange, yRange);
					gp.getChartPanel().setSize(1000, 800);
					String fName = csSiteName+"_"+sourceID+"_"+source.getName().replaceAll(" ", "_")
							+"_rup"+myRupID+"_"+periodStr+"s_"+comp.getShortName();
					gp.saveAsPNG(new File(subDir, fName+".png").getAbsolutePath());
					gp.saveAsPDF(new File(subDir, fName+".pdf").getAbsolutePath());
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (db != null)
				db.destroy();
		}
		
	}

}
