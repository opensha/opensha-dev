package scratch.kevin.cybershake;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.UncertainArbDiscDataset;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.siteData.OrderedSiteDataProviderList;
import org.opensha.commons.data.siteData.SiteDataValue;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.param.Parameter;
import org.opensha.sha.cybershake.calc.HazardCurveComputation;
import org.opensha.sha.cybershake.db.CachedPeakAmplitudesFromDB;
import org.opensha.sha.cybershake.db.CybershakeIM;
import org.opensha.sha.cybershake.db.CybershakeSite;
import org.opensha.sha.cybershake.db.Cybershake_OpenSHA_DBApplication;
import org.opensha.sha.cybershake.db.DBAccess;
import org.opensha.sha.cybershake.db.MeanUCERF2_ToDB;
import org.opensha.sha.cybershake.db.PeakAmplitudesFromDB;
import org.opensha.sha.cybershake.db.Runs2DB;
import org.opensha.sha.cybershake.db.SiteInfo2DB;
import org.opensha.sha.cybershake.db.CybershakeIM.CyberShakeComponent;
import org.opensha.sha.cybershake.db.CybershakeIM.IMType;
import org.opensha.sha.cybershake.db.CybershakeRun;
import org.opensha.sha.cybershake.plot.HazardCurvePlotter;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.util.SiteTranslator;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public class GMPEComparisonScatter {

	public static void main(String[] args) {
		List<Range> magRanges = Lists.newArrayList();
		List<Range> distRanges = Lists.newArrayList();
		
		magRanges.add(new Range(6d, 7d));
		distRanges.add(new Range(0d, 70d));
		
//		File outputDir = new File("/home/kevin/CyberShake/gmpe_comparison_scatter");
//		File outputDir = new File("/home/kevin/CyberShake/gmpe_comparison_scatter_15_12");
		File outputDir = new File("/home/kevin/CyberShake/gmpe_comparison_scatter_GP_2015_vs_2014");
		
//		List<Double> periods = Lists.newArrayList(3d, 5d, 10d);
//		List<Double> periods = Lists.newArrayList(0.2, 0.5, 3d, 5d, 10d);
		List<Double> periods = Lists.newArrayList(3d, 5d, 10d);
		
		CyberShakeComponent comp = CyberShakeComponent.RotD50;
		CyberShakeComponent compForRunID = CyberShakeComponent.GEOM_MEAN;

		
		DBAccess db = null;
		try {
			db = Cybershake_OpenSHA_DBApplication.getDB();
			ERF erf = MeanUCERF2_ToDB.createUCERF2ERF();
			
			ScalarIMR gmpe = AttenRelRef.NGAWest_2014_AVG_NOIDRISS.instance(null);
			
			CachedPeakAmplitudesFromDB amps2db = new CachedPeakAmplitudesFromDB(db, new File("/tmp"), erf);
			List<CybershakeIM> ims = amps2db.getIMs(periods, IMType.SA, comp);
			Runs2DB runs2db = new Runs2DB(db);
			SiteInfo2DB sites2db = new SiteInfo2DB(db);
			
			/* For fetching by site name, single dataset ID */
////			int datasetID = 57;
//			int datasetID = 61;
//			List<String> siteNames = Lists.newArrayList("CCP", "COO", "LADT", "LAPD", "P22", "PAS", "s429",
//					"s603", "s684", "s758", "SBSM", "SMCA", "STNI", "WNGC");
//			RunIDFetcher fetch = new RunIDFetcher(db, datasetID,
//					amps2db.getIMs(periods, IMType.SA, compForRunID).get(0).getID());
//			List<Integer> runIDs = Lists.newArrayList();
//			for (String siteName : siteNames)
//				runIDs.add(fetch.getRunID(siteName));
//			List<Integer> compRunIDs = null;
			
			/* For use with run IDs (multiple datasets */
//			List<Integer> runIDs = Lists.newArrayList(4683, 4686, 4681, 4685);
//			List<Integer> compRunIDs = null;
			List<Integer> runIDs = Lists.newArrayList(4683, 4681);
			List<Integer> compRunIDs = Lists.newArrayList(4686, 4685);
			
			Preconditions.checkState(compRunIDs == null || compRunIDs.size() == runIDs.size());
			
			for (int r=0; r<runIDs.size(); r++) {
				int runID = runIDs.get(r);
				int compRunID = -1;
				if (compRunIDs != null)
					compRunID = compRunIDs.get(r);
				CybershakeRun run = runs2db.getRun(runID);
				CybershakeSite csSite = sites2db.getSiteFromDB(run.getSiteID());
				System.out.println("CS Site: "+csSite+", run "+runID);
				Location loc = csSite.createLocation();
				
				gmpe.setIntensityMeasure(SA_Param.NAME);
				if (compRunIDs == null) {
					int velModelID = run.getVelModelID();
					Site site = new Site(loc);
					gmpe.setParamDefaults();
					for (Parameter<?> param : gmpe.getSiteParams())
						if (!site.containsParameter(param))
							site.addParameter((Parameter<?>)param.clone());
					OrderedSiteDataProviderList provs = HazardCurvePlotter.createProviders(velModelID);
					List<SiteDataValue<?>> datas = provs.getAllAvailableData(site.getLocation());
					SiteTranslator trans = new SiteTranslator();
					System.out.println("Site data:");
					for (Parameter<?> param : site) {
						trans.setParameterValue(param, datas);
						System.out.println("\t"+param.getName()+": "+param.getValue());
					}
					gmpe.setSite(site);
				}
				
				List<Integer> sourceIDs = sites2db.getSrcIdsForSite(run.getSiteID(), run.getERFID());
				
				for (int p=0; p<periods.size(); p++) {
					double period = periods.get(p);
					String periodStr;
					if (period == (int)period)
						periodStr = (int)period+"";
					else
						periodStr = (float)period+"";
					CybershakeIM im = ims.get(p);
					
					SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), period);
					
					for (Range magRange : magRanges) {
						for (Range distRange : distRanges) {
							System.out.println("Calculating for "+period+"s, mag "+rangeFileName(magRange)
										+", dist "+rangeFileName(distRange));
							DefaultXY_DataSet xy = new DefaultXY_DataSet();
							for (int sourceID : sourceIDs) {
								ProbEqkSource source = erf.getSource(sourceID);
								for (int rupID=0; rupID<source.getNumRuptures(); rupID++) {
									ProbEqkRupture rup = source.getRupture(rupID);
									if (!magRange.contains(rup.getMag()))
										continue;
									if (!distRange.contains(rup.getRuptureSurface().getDistanceRup(loc)))
										continue;
									// inside both ranges
									
//									System.out.println("Source "+sourceID+", rup "+rupID);
									
									List<Double> csVals = amps2db.getIM_Values(runID, sourceID, rupID, im);
									for (int i=0; i<csVals.size(); i++)
										// convert to G
										csVals.set(i, csVals.get(i)/HazardCurveComputation.CONVERSION_TO_G);
									double csMean = logAverage(csVals);
									if (compRunIDs == null) {
										gmpe.setEqkRupture(rup);
										double gmpeMean = Math.exp(gmpe.getMean());
										
										xy.set(gmpeMean, csMean);
									} else {
										// CS vs CS
										csVals = amps2db.getIM_Values(compRunID, sourceID, rupID, im);
										for (int i=0; i<csVals.size(); i++)
											// convert to G
											csVals.set(i, csVals.get(i)/HazardCurveComputation.CONVERSION_TO_G);
										double csMean2 = logAverage(csVals);
										
										xy.set(csMean, csMean2);
									}
								}
							}
							System.out.println("Found "+xy.size()+" comparisons");
							Range range = new Range(Math.min(xy.getMinX(), xy.getMinY()),
									Math.max(xy.getMaxX(), xy.getMaxY()));
							
							List<XY_DataSet> funcs = Lists.newArrayList();
							List<PlotCurveCharacterstics> chars = Lists.newArrayList();
							
							// shaded yellow area
							ArbitrarilyDiscretizedFunc oneToOne = new ArbitrarilyDiscretizedFunc();
							oneToOne.set(xy.getMinX(), xy.getMinX());
							oneToOne.set(xy.getMaxX(), xy.getMaxX());
							oneToOne.set(range.getLowerBound(), range.getLowerBound());
							oneToOne.set(range.getUpperBound(), range.getUpperBound());
							oneToOne.set(range.getLowerBound()*2, range.getLowerBound()*2);
							oneToOne.set(range.getLowerBound()/2, range.getLowerBound()/2);
							oneToOne.set(range.getUpperBound()*2, range.getUpperBound()*2);
							oneToOne.set(range.getUpperBound()/2, range.getUpperBound()/2);
							ArbitrarilyDiscretizedFunc upper = new ArbitrarilyDiscretizedFunc();
							for (Point2D pt : oneToOne)
								upper.set(pt.getX(), pt.getY()*2);
							ArbitrarilyDiscretizedFunc lower = new ArbitrarilyDiscretizedFunc();
							for (Point2D pt : oneToOne)
								lower.set(pt.getX(), pt.getY()/2);
							UncertainArbDiscDataset shaded = new UncertainArbDiscDataset(oneToOne, lower, upper);
							
							funcs.add(shaded);
							chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN_TRANS, 1f, Color.YELLOW));
							
							funcs.add(upper);
							chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
							
							funcs.add(oneToOne);
							chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GRAY));
							
							funcs.add(lower);
							chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
							
							funcs.add(xy);
							chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, Color.RED));
							
							// regression in log space
							SimpleRegression regression = new SimpleRegression();
							for (Point2D pt : xy)
								regression.addData(Math.log(pt.getX()), Math.log(pt.getY()));
							double b = regression.getIntercept();
							double m = regression.getSlope();
							DefaultXY_DataSet fit = new DefaultXY_DataSet();
							// use one to one for x values
							for (int i=0; i<oneToOne.size(); i++) {
								double origX = oneToOne.getX(i);
								if (origX < xy.getMinX() || origX > xy.getMaxX())
									continue;
								double x = Math.log(origX);
								double y = m*x + b;
								fit.set(Math.exp(x), Math.exp(y));
							}
							funcs.add(fit);
							chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GREEN.darker()));
							
							String csSiteName = csSite.short_name;
							String prefix = csSiteName+"_mag"+rangeFileName(magRange)+"_dist"+rangeFileName(distRange)
									+"_"+periodStr+"s_run"+runID;
							
							String xAxisLabel, yAxisLabel;
							
							if (compRunID >= 0) {
								prefix += "_vs_"+compRunID;
								xAxisLabel = "CyberShake run "+runID+" "+periodStr+"s SA (g)";
								yAxisLabel = "CyberShake run "+compRunID+" "+periodStr+"s SA (g)";
							} else {
								xAxisLabel = "Empirical "+periodStr+"s SA (g)";
								yAxisLabel = "CyberShake "+periodStr+"s SA (g)";
							}
							
							String title = csSiteName;
							
							PlotSpec spec = new PlotSpec(funcs, chars, title, xAxisLabel, yAxisLabel);
							
							// bottom up
							List<String> annLines = Lists.newArrayList();
							annLines.add(xy.size()+" Ruptures");
							annLines.add(str(distRange.getLowerBound())+" km < Dist < "+str(distRange.getUpperBound())+" km");
							annLines.add(str(magRange.getLowerBound())+" < Mw < "+str(magRange.getUpperBound()));
							
							// this lines up annotations in log space
							EvenlyDiscretizedFunc lineSpacings = new EvenlyDiscretizedFunc(
									Math.log(range.getLowerBound()), Math.log(range.getUpperBound()), 17);
							
							List<XYTextAnnotation> anns = Lists.newArrayList();
							for (int i=0; i<annLines.size(); i++) {
								double x = Math.exp(lineSpacings.getX(lineSpacings.size()-2));
								double y = Math.exp(lineSpacings.getX(1+i));
								String annText = annLines.get(i);
								XYTextAnnotation ann = new XYTextAnnotation(annText, x, y);
								ann.setTextAnchor(TextAnchor.BOTTOM_RIGHT);
								ann.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 24));
								anns.add(ann);
							}
							spec.setPlotAnnotations(anns);
							
							HeadlessGraphPanel gp = new HeadlessGraphPanel();
							gp.setTickLabelFontSize(18);
							gp.setAxisLabelFontSize(20);
							gp.setPlotLabelFontSize(21);
							gp.setBackgroundColor(Color.WHITE);
							
							gp.drawGraphPanel(spec, true, true, range, range);
							gp.getChartPanel().setSize(900, 800);
							gp.saveAsPNG(new File(outputDir, prefix+".png").getAbsolutePath());
							gp.saveAsPDF(new File(outputDir, prefix+".pdf").getAbsolutePath());
						}
					}
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (db != null)
				db.destroy();
		}
	}
	
	private static String str(double num) {
		if (num == (int)num)
			return (int)num+"";
		return (float)num+"";
	}
	
	private static String rangeFileName(Range r) {
		return str(r.getLowerBound())+"to"+str(r.getUpperBound());
	}
	
	private static double logAverage(List<Double> linearVals) {
		double[] vals = new double[linearVals.size()];
		for (int i=0; i<linearVals.size(); i++)
			vals[i] = Math.log(linearVals.get(i));
		return Math.exp(StatUtils.mean(vals));
	}

}
