package scratch.kevin.miscFigures;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.siteData.OrderedSiteDataProviderList;
import org.opensha.commons.data.siteData.SiteDataValue;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.param.Parameter;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.calc.Wald_MMI_Calc;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGV_Param;
import org.opensha.sha.util.SiteTranslator;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO.ETAS_Catalog;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.analysis.ETAS_AbstractPlot;
import scratch.UCERF3.erf.ETAS.launcher.ETAS_Config;
import scratch.UCERF3.erf.mean.MeanUCERF3;
import scratch.UCERF3.erf.mean.MeanUCERF3.Presets;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;

public class HazardCurveTimeScales {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/tmp");
		String prefix = "mmi_hazard_curve";
		
		double[] durations = { 1d/52d, 1d/12d, 1d };
//		Color[] colors = { Color.GREEN.darker(), Color.BLUE.darker(), Color.RED.darker() };
		Color[] colors = { new Color(200, 200, 200), new Color(100, 100, 100), Color.BLACK };
		
		MeanUCERF3 erf = new MeanUCERF3();
		erf.setPreset(Presets.BOTH_FM_BRANCH_AVG);
		ScalarIMR gmpe = AttenRelRef.ASK_2014.instance(null);
		gmpe.setParamDefaults();
		
		Site site = new Site(new Location(34.0192, -118.286)); // USC
		
//		Site site = new Site(new Location(35.052463, -118.173721)); // Mojave
//		prefix = "mojave_"+prefix;
		
		String xLabel = "Modified Mercalli Intensity (MMI)";
		
		OrderedSiteDataProviderList provs = OrderedSiteDataProviderList.createSiteDataProviderDefaults();
		SiteTranslator trans = new SiteTranslator();
		ArrayList<SiteDataValue<?>> datas = provs.getBestAvailableData(site.getLocation());
		
		for (Parameter<?> param : gmpe.getSiteParams()) {
			trans.setParameterValue(param, datas);
			site.addParameter(param);
		}

		EvenlyDiscretizedFunc logPGAVals = new EvenlyDiscretizedFunc(Math.log(1e-6), Math.log(2d), 500);
		EvenlyDiscretizedFunc logPGVVals = new EvenlyDiscretizedFunc(Math.log(1e-2), Math.log(500), 500);
		DiscretizedFunc pgaXVals = new ArbitrarilyDiscretizedFunc();
		DiscretizedFunc pgvXVals = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : logPGAVals)
			pgaXVals.set(Math.exp(pt.getX()), 1d);
		for (Point2D pt : logPGVVals)
			pgvXVals.set(Math.exp(pt.getX()), 1d);
//		DiscretizedFunc pgaXVals = new IMT_Info().getDefaultHazardCurve(PGA_Param.NAME);
//		DiscretizedFunc pgvXVals = new IMT_Info().getDefaultHazardCurve(PGV_Param.NAME);
//		DiscretizedFunc logPGAVals = new ArbitrarilyDiscretizedFunc();
//		for (Point2D pt : pgaXVals)
//			logPGAVals.set(Math.log(pt.getX()), 1d);
		
		double maxY = 1d;
		double minY = 1e-6;
		EvenlyDiscretizedFunc logYVals = new EvenlyDiscretizedFunc(
				Math.log10(minY), Math.log10(maxY), 10000);
		ArbitrarilyDiscretizedFunc yVals = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : logYVals)
			yVals.set(Math.pow(10, pt.getX()), 0d);
		EvenlyDiscretizedFunc discrYVals = new EvenlyDiscretizedFunc(0.5, 1d, 1000);
		for (Point2D pt : discrYVals)
			yVals.set(pt.getX(), 0d);
		
		HazardCurveCalculator calc = new HazardCurveCalculator();
		calc.setMaxSourceDistance(500d);
		
		File etasDir = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2019_09_04-ComCatM7p1_ci38457511_ShakeMapSurfaces");
		ETAS_Config config = ETAS_Config.readJSON(new File(etasDir, "config.json"));
		File binFile = new File(etasDir, "results_m5_preserve_chain.bin");
		FaultSystemSolution fss = config.loadFSS();
		
		List<ETAS_Catalog> catalogs = ETAS_CatalogIO.loadCatalogsBinary(binFile, 5d);
		
		List<DiscretizedFunc> tiFuncs = new ArrayList<>();
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(24);
		gp.setPlotLabelFontSize(24);
		gp.setLegendFontSize(28);
		gp.setBackgroundColor(Color.WHITE);
//		gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
		
		for (boolean etas : new boolean[] {false, true}) {
			
			List<DiscretizedFunc> combFuncs = new ArrayList<>();
			List<PlotCurveCharacterstics> combChars = new ArrayList<>();
			
			if (etas)
				prefix += "_etas";
			else
				tiFuncs = combFuncs;
			for (double duration : durations) {
				erf.getTimeSpan().setDuration(duration);
				erf.updateForecast();
				
				System.out.println("Duration: "+duration+" ETAS? "+etas);
				
				AbstractERF myERF = etas ? new ETAS_ERF(erf, config, fss, catalogs) : erf;
				
				List<DiscretizedFunc> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				System.out.println("Calculating curves");
				
				// PGA curve
				gmpe.setIntensityMeasure(PGA_Param.NAME);
				calc.getHazardCurve(logPGAVals, site, gmpe, myERF);
				
				DiscretizedFunc pgaCurve = pgaXVals.deepClone();
				for (int j=0; j<logPGAVals.size(); j++)
					pgaCurve.set(j, logPGAVals.getY(j));
				System.out.println("\tDone PGA");
				
				// PGV curve
				gmpe.setIntensityMeasure(PGV_Param.NAME);
				calc.getHazardCurve(logPGVVals, site, gmpe, myERF);
				
				DiscretizedFunc pgvCurve = pgvXVals.deepClone();
				for (int j=0; j<logPGVVals.size(); j++)
					pgvCurve.set(j, logPGVVals.getY(j));
				System.out.println("\tDone PGA");
				
//				XY_DataSet mmiCurve = new DefaultXY_DataSet();
				DiscretizedFunc mmiCurve = new ArbitrarilyDiscretizedFunc();
				for (int i=0; i<yVals.size(); i++) {
					double prob = yVals.getX(i);
					if (prob > pgaCurve.getMaxY() || prob > pgvCurve.getMaxY())
						continue;
					if (prob < pgaCurve.getMinY() || prob < pgvCurve.getMinY())
						continue;
					double pga = pgaCurve.getFirstInterpolatedX_inLogXLogYDomain(prob);
					double pgv = pgvCurve.getFirstInterpolatedX_inLogXLogYDomain(prob);
					double mmi = Wald_MMI_Calc.getMMI(pga, pgv);
					if (mmi > 1)
						mmiCurve.set(mmi, prob);
				}
				
				combFuncs.add(mmiCurve);
				combChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, colors[combChars.size()]));
				
				System.out.println("MMI vals: "+mmiCurve.size());
				
				funcs.add(mmiCurve);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
				
				String durLabel = ETAS_AbstractPlot.getTimeLabel(duration, false);
				String durPrefix = ETAS_AbstractPlot.getTimeShortLabel(duration).replaceAll(" ", "");
				
				mmiCurve.setName(durLabel+" ");
				
				PlotSpec spec = new PlotSpec(funcs, chars, " ", xLabel, durLabel+" Probability of Exceedance");
				spec.setLegendVisible(false);
				
				Range xRange = new Range(1, 10);
				Range yRange = new Range(minY, maxY);
				
				gp.drawGraphPanel(spec, false, true, xRange, yRange);
				
				File file = new File(outputDir, prefix+"_"+durPrefix);
				gp.getChartPanel().setSize(800, 600);
				gp.saveAsPNG(file.getAbsolutePath()+".png");
				gp.saveAsPDF(file.getAbsolutePath()+".pdf");
				
				gp.drawGraphPanel(spec, false, false, xRange, new Range(0d, 1d));
				
				file = new File(outputDir, prefix+"_"+durPrefix+"_linear");
				gp.getChartPanel().setSize(800, 600);
				gp.saveAsPNG(file.getAbsolutePath()+".png");
				gp.saveAsPDF(file.getAbsolutePath()+".pdf");
			}
			
			PlotSpec spec = new PlotSpec(combFuncs, combChars, " ", xLabel, "Probability of Exceedance");
			spec.setLegendVisible(true);
			
			Range xRange = new Range(1, 10);
			Range yRange = new Range(minY, maxY);
			
			gp.drawGraphPanel(spec, false, true, xRange, yRange);
			
			File file = new File(outputDir, prefix);
			gp.getChartPanel().setSize(800, 600);
			gp.saveAsPNG(file.getAbsolutePath()+".png");
			gp.saveAsPDF(file.getAbsolutePath()+".pdf");
			
			List<XYTextAnnotation> anns = new ArrayList<>();
			Font annFont = new Font(Font.SANS_SERIF, Font.BOLD, 22);
			DecimalFormat pDF = new DecimalFormat("0.0%");
			double annY = 0.95;
			TextAnchor anchor = TextAnchor.CENTER_RIGHT;
			double deltaY = 0.1;
			for (int d=0; d<durations.length; d++) {
				double prob = combFuncs.get(d).getInterpolatedY_inLogYDomain(6d);
				String label = ETAS_AbstractPlot.getTimeLabel(durations[d], false)+" P(MMI≥6)";
				if (prob*100d < 0.1)
					label += " < 0.1%";
				else
					label += " = "+pDF.format(prob);
				XYTextAnnotation ann = new XYTextAnnotation(label, 9.75, annY);
				ann.setFont(annFont);
				ann.setTextAnchor(anchor);
				anns.add(ann);
				System.out.println("\t"+label);
				annY -= deltaY;
			}
			spec.setPlotAnnotations(anns);
			
			gp.drawGraphPanel(spec, false, false, xRange, new Range(0d, 1d));
			
			file = new File(outputDir, prefix+"_linear");
			gp.getChartPanel().setSize(800, 600);
			gp.saveAsPNG(file.getAbsolutePath()+".png");
			gp.saveAsPDF(file.getAbsolutePath()+".pdf");
			
			if (etas) {
				// now gain
				List<DiscretizedFunc> gainFuncs = new ArrayList<>();
				List<PlotCurveCharacterstics> gainChars = new ArrayList<>();
				
				double maxGain = 1d;
				
				for (int d=0; d<durations.length; d++) {
					ArbitrarilyDiscretizedFunc gainFunc = new ArbitrarilyDiscretizedFunc();
					gainFunc.setName(ETAS_AbstractPlot.getTimeLabel(durations[d], false)+" ");
					
					DiscretizedFunc etasFunc = combFuncs.get(d);
					DiscretizedFunc tiFunc = tiFuncs.get(d);
					
					for (Point2D pt : etasFunc) {
						if (pt.getX() >= tiFunc.getMinX() && pt.getX() <= tiFunc.getMaxX()) {
							double gain = pt.getY()/tiFunc.getInterpolatedY_inLogYDomain(pt.getX());
							gainFunc.set(pt.getX(), gain);
							maxGain = Math.max(maxGain, gain);
						}
					}
					
					gainFuncs.add(gainFunc);
					gainChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, colors[d]));
				}
				System.out.println("Max gain: "+maxGain);
//				maxGain = Math.max(maxGain, 10000);
				
				spec = new PlotSpec(gainFuncs, gainChars, " ", xLabel, "Gain over Time-Independent");
				spec.setLegendVisible(true);
				
				gp.drawGraphPanel(spec, false, true, xRange,
						new Range(1, Math.pow(10, Math.ceil(Math.log10(maxGain)))));
				
				file = new File(outputDir, prefix+"_gains");
				gp.getChartPanel().setSize(800, 600);
				gp.saveAsPNG(file.getAbsolutePath()+".png");
				gp.saveAsPDF(file.getAbsolutePath()+".pdf");
				
				gp.drawGraphPanel(spec, false, false, xRange,
						new Range(1, 10d*Math.ceil(maxGain/10d)));
				
				file = new File(outputDir, prefix+"_gains_linear");
				gp.getChartPanel().setSize(800, 600);
				gp.saveAsPNG(file.getAbsolutePath()+".png");
				gp.saveAsPDF(file.getAbsolutePath()+".pdf");
			}
		}
	}
	
	private static class ETAS_ERF extends AbstractERF {
		
		private List<ProbEqkSource> etasSources;
		private AbstractERF erf;
		
		public ETAS_ERF(AbstractERF erf, ETAS_Config config, FaultSystemSolution fss, List<ETAS_Catalog> catalogs) {
			this.erf = erf;
			etasSources = new ArrayList<>();
			
			long startTime = config.getSimulationStartTimeMillis();
			long endTime = startTime + 
					(long)(erf.getTimeSpan().getDuration()*ProbabilityModelsCalc.MILLISEC_PER_YEAR);
			
			FaultSystemRupSet rupSet = fss.getRupSet();
			
			double probEach = 1d/catalogs.size();
			for (ETAS_Catalog catalog : catalogs) {
				for (ETAS_EqkRupture rup : catalog) {
					if (rup.getOriginTime() > endTime)
						break;
					if (rup.getMag() < 5d)
						continue;
					
					int fssIndex = rup.getFSSIndex();
					RuptureSurface surf;
					double rake;
					if (fssIndex > 0) {
						surf = rupSet.getSurfaceForRupture(fssIndex, 1d);
						rake = rupSet.getAveRakeForRup(fssIndex);
					} else {
						PointSurface psurf = new PointSurface(rup.getHypocenterLocation());
						psurf.setAveDip(90d);
						psurf.setAveStrike(0);
						psurf.setAveWidth(1d);
						rake = 180d;
						surf = psurf;
					}
					ProbEqkRupture probRup = new ProbEqkRupture(rup.getMag(), rake,
							probEach, surf, rup.getHypocenterLocation());
					ProbEqkSource source = new ProbEqkSource() {
						
						@Override
						public RuptureSurface getSourceSurface() {
							return surf;
						}
						
						@Override
						public LocationList getAllSourceLocs() {
							return surf.getEvenlyDiscritizedListOfLocsOnSurface();
						}
						
						@Override
						public ProbEqkRupture getRupture(int nRupture) {
							return probRup;
						}
						
						@Override
						public int getNumRuptures() {
							return 1;
						}
						
						@Override
						public double getMinDistance(Site site) {
							return surf.getDistanceJB(site.getLocation());
						}

						@Override
						public boolean isSourcePoissonian() {
							return false;
						}

						@Override
						public boolean isPoissonianSource() {
							return false;
						}
					};
					etasSources.add(source);
				}
			}
			
			System.out.println("Have "+etasSources.size()+" etas sources");
		}

		@Override
		public int getNumSources() {
			return erf.getNumSources()+etasSources.size();
		}

		@Override
		public ProbEqkSource getSource(int idx) {
			if (idx < erf.getNumSources())
				return erf.getSource(idx);
			return etasSources.get(idx - erf.getNumSources());
		}

		@Override
		public void updateForecast() {
			
		}

		@Override
		public String getName() {
			return null;
		}
		
	}

}
