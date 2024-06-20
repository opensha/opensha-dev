package scratch.kevin.simulators.ruptures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.EqkRuptureParams.MagParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGV_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.utils.RSQSimEqkRupture;

import com.google.common.base.Preconditions;
import com.google.common.collect.BiMap;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.HashBiMap;
import com.google.common.collect.Table;

import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class ResidualStdDevMagUncertPlot {

	public static void main(String[] args) throws IOException {
		RSQSimCatalog catalog = Catalogs.BRUCE_5652.instance();
		File bbpDir = new File("/data/kevin/bbp/parallel/2023_08_29-rundir5652-all-m6.5-skipYears5000-noHF-vmLA_BASIN_500-griddedSites");
		File bbpFile = new File(bbpDir, "results_rotD.zip");
		VelocityModel vm = VelocityModel.LA_BASIN_500;
		AttenRelRef gmmRef = AttenRelRef.ASK_2014;
		
		String prefix = "r5652_ask2014";
		
		File outputDir = new File("/home/kevin/Documents/papers/2023_Bruce_GMM_Multifault/figures/residual_mag_uncert");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		int numSamplesEach = 100;
		EvenlyDiscretizedFunc dmFunc = new EvenlyDiscretizedFunc(0, 0.5, 11);
		System.out.println("Distance bins:");
		List<Double> deltaMags = new ArrayList<>();
		for (int i=0; i<dmFunc.size(); i++) {
			System.out.println("\t"+(float)dmFunc.getX(i));
			deltaMags.add(dmFunc.getX(i));
		}
		
		List<BBP_Site> bbpSites = BBP_Site.readFile(bbpDir);
		BiMap<BBP_Site, Site> bbpToSiteMap = HashBiMap.create();
		List<Site> sites = new ArrayList<>(bbpSites.size());
		for (BBP_Site bbpSite : bbpSites) {
			Site site = bbpSite.buildGMPE_Site(vm);
			bbpToSiteMap.put(bbpSite, site);
			sites.add(site);
		}
		
		List<RSQSimEvent> events = catalog.loader().minMag(6.5d).load();
		Map<Integer, RSQSimEvent> eventsMap = new HashMap<>();
		for (RSQSimEvent e : events)
			eventsMap.put(e.getID(), e);
		BBP_CatalogSimZipLoader bbpLoader = new BBP_CatalogSimZipLoader(bbpFile, bbpSites, bbpToSiteMap, eventsMap);
		
		double[] periods = {-1d, 2, 3, 5, 10};
		
		List<Range> magRanges = new ArrayList<>();
		List<String> magRangeLabels = new ArrayList<>();
		
		magRanges.add(new Range(0d, 6.75));
		magRangeLabels.add("M6.5-6.75");
		
		magRanges.add(new Range(6.75, 7d));
		magRangeLabels.add("M6.75-7");
		
		magRanges.add(new Range(7d, 7.25));
		magRangeLabels.add("M7-7.25");
		
		magRanges.add(new Range(7.25, 7.5));
		magRangeLabels.add("M7.25-7.5");
		
		magRanges.add(new Range(7.5, 7.75));
		magRangeLabels.add("M7.5-7.75");
		
		magRanges.add(new Range(7.75, 10d));
		magRangeLabels.add("M>7.75");
		
		List<Range> distRanges = new ArrayList<>();
		List<String> distRangeLabels = new ArrayList<>();
		
		distRanges.add(new Range(0d, 10d));
		distRangeLabels.add("R<10 km");
		
		distRanges.add(new Range(10d, 20d));
		distRangeLabels.add("10<R<20 km");
		
		distRanges.add(new Range(20d, 40d));
		distRangeLabels.add("20<R<40 km");
		
		distRanges.add(new Range(40d, 80d));
		distRangeLabels.add("40<R<80 km");
		
		distRanges.add(new Range(80d, 160d));
		distRangeLabels.add("80<R<160 km");
		
		distRanges.add(new Range(160d, 1000d));
		distRangeLabels.add("R>160 km");
		
		List<Residuals> allResiduals = initResuduals(deltaMags.size(), periods, magRanges, distRanges);
		
		ExecutorService exec = Executors.newFixedThreadPool(FaultSysTools.defaultNumThreads());
		
		List<Future<?>> siteFutures = new ArrayList<>();
		for (Site site : sites)
			siteFutures.add(exec.submit(new SiteCalcRunnable(catalog, bbpLoader, site, periods, gmmRef, deltaMags,
					numSamplesEach, allResiduals)));
		
		for (Future<?> future : siteFutures) {
			try {
				future.get();
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(1);
			}
		}
		
		exec.shutdown();

		plot(outputDir, prefix+"_residual_std_all_periods", dmFunc, allResiduals, null, false, false, false, periods,
				magRanges, magRangeLabels, distRanges, distRangeLabels);
		plot(outputDir, prefix+"_residual_std_all_periods_by_period", dmFunc, allResiduals, null, false, false, true, periods,
				magRanges, magRangeLabels, distRanges, distRangeLabels);
		for (int p=-1; p<periods.length; p++) {
			String perPrefix;
			Double period;
			if (p < 0) {
				perPrefix = prefix+"_residual_std_all_periods";
				period = null;
			} else if (periods[p] == -1) {
				perPrefix = prefix+"_residual_std_pgv";
				period = -1d;
			} else {
				perPrefix = prefix+"_residual_std_"+oDF.format(periods[p])+"s";
				period = periods[p];
			}
			plot(outputDir, perPrefix+"_by_mag", dmFunc, allResiduals, period, true, false, false, periods,
					magRanges, magRangeLabels, distRanges, distRangeLabels);
			plot(outputDir, perPrefix+"_by_dist", dmFunc, allResiduals, period, false, true, false, periods,
					magRanges, magRangeLabels, distRanges, distRangeLabels);
		}
	}
	
	private static final DecimalFormat oDF = new DecimalFormat("0.##");
	
	private static List<Residuals> initResuduals(int numDMBins, double[] periods, List<Range> magRanges, List<Range> distRanges) {
		List<Residuals> ret = new ArrayList<>(numDMBins);
		for (int i=0; i<numDMBins; i++)
			ret.add(new Residuals(periods, magRanges, distRanges));
		return ret;
	}
	
	private static class Residuals {
		private double[] periods;
		private List<Range> magRanges;
		private List<Range> distRanges;
		
		// <mag range, dist range, period resisuals>
		private Table<Range, Range, List<List<Double>>> residualsTable;
		private Table<Range, Range, List<Double>> detrendTable;
		
		public Residuals(double[] periods, List<Range> magRanges, List<Range> distRanges) {
			this.periods = periods;
			this.magRanges = magRanges;
			this.distRanges = distRanges;
			
			residualsTable = HashBasedTable.create();
			
			for (Range magRange : magRanges) {
				for (Range distRange : distRanges) {
					List<List<Double>> periodResiduals = new ArrayList<>(periods.length);
					for (int p=0; p<periods.length; p++)
						periodResiduals.add(new ArrayList<>());
					residualsTable.put(magRange, distRange, periodResiduals);
				}
			}
		}
		
		public synchronized void add(double mag, double dist, double[] residuals) {
			Preconditions.checkState(residuals.length == periods.length);
			Range magRange = null;
			for (Range test : magRanges) {
				if (test.contains(mag)) {
					magRange = test;
					break;
				}
			}
			Preconditions.checkNotNull(magRange);
			Range distRange = null;
			for (Range test : distRanges) {
				if (test.contains(dist)) {
					distRange = test;
					break;
				}
			}
			Preconditions.checkNotNull(distRange);
			List<List<Double>> periodResiduals = residualsTable.get(magRange, distRange);
			for (int i=0; i<residuals.length; i++)
				periodResiduals.get(i).add(residuals[i]);
		}
		
		private synchronized void checkCalcDetrend() {
			if (detrendTable == null) {
				System.out.println("Calculating detrend averages");
				detrendTable = HashBasedTable.create();
				for (Range magRange : magRanges) {
					for (Range distRange : distRanges) {
						List<List<Double>> perResids = residualsTable.get(magRange, distRange);
						if (perResids.get(0).isEmpty())
							continue;
						List<Double> averages = new ArrayList<>(periods.length);
						for (int p=0; p<periods.length; p++)
							averages.add(perResids.get(p).stream().mapToDouble(D->D).average().getAsDouble());
						System.out.println("Averages for magRange="+magRange+", distRange, N="
								+perResids.get(0).size()+":"+distRange+": "+averages);
						detrendTable.put(magRange, distRange, averages);
						Preconditions.checkNotNull(detrendTable.get(magRange, distRange));
					}
				}
			}
		}
		
		private double calcStdDev(Range restrictMagRange, Range restrictDistRange, Double restrictPeriod) {
			checkCalcDetrend();
			StandardDeviation stdDev = new StandardDeviation();
			int num = 0;
			for (Range magRange : magRanges) {
				if (restrictMagRange != null && restrictMagRange != magRange)
					continue;
				for (Range distRange : distRanges) {
					if (restrictDistRange != null && restrictDistRange != distRange)
						continue;
					List<List<Double>> perResids = residualsTable.get(magRange, distRange);
					if (perResids.get(0).isEmpty())
						continue;
					List<Double> averages = detrendTable.get(magRange, distRange);
					Preconditions.checkNotNull(averages,
							"detrend averages is null for magRange=%s, distRange=%s, but perResids.get(0).size()=%s",
							magRange, distRange, perResids.get(0).size());
					for (int p=0; p<periods.length; p++) {
						if (restrictPeriod != null && restrictPeriod != periods[p])
							continue;
						double average = averages.get(p);
						for (Double residual : perResids.get(p)) {
							stdDev.increment(residual - average);
							num++;
						}
					}
				}
			}
			System.out.println("Calculating stdDev for "+num+" values");
			return stdDev.getResult();
		}
	}
	
	private static class SiteCalcRunnable implements Runnable {
		
		// inputs
		private RSQSimCatalog catalog;
		private BBP_CatalogSimZipLoader bbpLoader;
		private Site site;
		private double[] periods;
		private AttenRelRef gmmRef;
		private List<Double> deltaMags;
		private int numSamplesEach;
		
		// outputs
		private List<Residuals> allResiduals;
		
		private SiteCalcRunnable(RSQSimCatalog catalog, BBP_CatalogSimZipLoader bbpLoader, Site site, double[] periods,
				AttenRelRef gmmRef, List<Double> deltaMags, int numSamplesEach, List<Residuals> allResiduals) {
			super();
			this.catalog = catalog;
			this.bbpLoader = bbpLoader;
			this.site = site;
			this.periods = periods;
			this.gmmRef = gmmRef;
			this.deltaMags = deltaMags;
			this.numSamplesEach = numSamplesEach;
			this.allResiduals = allResiduals;
		}

		@Override
		public void run() {
			ScalarIMR gmm = gmmRef.get();
			Collection<RSQSimEvent> events = bbpLoader.getRupturesForSite(site);
			
			System.out.println("Starting site "+site.getName()+" with "+events.size()+" events");
			Random r = new Random((long)events.size()*(long)numSamplesEach
					+ Double.doubleToLongBits(site.getLocation().lat)
					+ Double.doubleToLongBits(site.getLocation().lon));
			
			gmm.setSite(site);
			
			MagParam magParam = (MagParam) gmm.getParameter(MagParam.NAME);
			
			boolean hasPGV = false;
			for (double period : periods)
				if (period == -1d)
					hasPGV = true;

			try {
				for (RSQSimEvent event : events) {
					RSQSimEqkRupture rup = catalog.getEqkRupture(event);
					
					double mag = rup.getMag();
					double dist = rup.getRuptureSurface().getDistanceRup(site.getLocation());
					
					DiscretizedFunc simRotD = bbpLoader.getRotD50(site, event, 0);
					double pgv = hasPGV ? bbpLoader.getPGV(site, event, 0) : Double.NaN;
					
					gmm.setEqkRupture(rup);
					
					for (int dmI=0; dmI<deltaMags.size(); dmI++) {
						double dm = deltaMags.get(dmI);
						int num = dm > 0 ? numSamplesEach : 1;
						for (int i=0; i<num; i++) {
							double modMag = dm > 0 ? mag + (2d*r.nextDouble()-0.5)*dm : mag; // this is two sided uncertainty
							magParam.setValue(modMag);
							double[] residuals = new double[periods.length];
							for (int p=0; p<periods.length; p++) {
								double simVal;
								if (periods[p] == -1) {
									gmm.setIntensityMeasure(PGV_Param.NAME);
									simVal = pgv;
								} else {
									Preconditions.checkState(periods[p] > 0);
									gmm.setIntensityMeasure(SA_Param.NAME);
									SA_Param.setPeriodInSA_Param(gmm.getIntensityMeasure(), periods[p]);
									simVal = simRotD.getInterpolatedY(periods[p]);
								}
								double gmmMean = gmm.getMean();
								residuals[p] = Math.log(simVal) - gmmMean;
							}
							allResiduals.get(dmI).add(mag, dist, residuals);
						}
					}
				}
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(1);
			}
			System.out.println("Completed site "+site.getName());
		}
		
	}
	
	private static void plot(File outputDir, String prefix, EvenlyDiscretizedFunc dmFunc, List<Residuals> allResiduals, Double period,
			boolean doMagBinned, boolean doDistBinned, boolean doPerBinned,
			double[] periods, List<Range> magRanges, List<String> magRangeLabels, List<Range> distRanges, List<String> distRangeLabels)
					throws IOException {
		System.out.println("Building plot: "+prefix);
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		if (doMagBinned) {
			CPT magCPT = GMT_CPT_Files.SEQUENTIAL_NAVIA_UNIFORM.instance().rescale(0d, magRanges.size()-1);
			
			for (int m=0; m<magRanges.size(); m++) {
				EvenlyDiscretizedFunc sdFunc = calcSDFunc(dmFunc, allResiduals, magRanges.get(m), null, period);
				sdFunc.setName(magRangeLabels.get(m));
				funcs.add(sdFunc);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, magCPT.getColor((float)m)));
			}
		}
		
		if (doDistBinned) {
			CPT distCPT = GMT_CPT_Files.SEQUENTIAL_LAJOLLA_UNIFORM.instance().rescale(0d, distRanges.size()-1);
			
			for (int d=0; d<distRanges.size(); d++) {
				EvenlyDiscretizedFunc sdFunc = calcSDFunc(dmFunc, allResiduals, null, distRanges.get(d), period);
				sdFunc.setName(distRangeLabels.get(d));
				funcs.add(sdFunc);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, distCPT.getColor((float)d)));
			}
		}
		
		if (doPerBinned) {
			CPT perCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(0d, periods.length-1);
			
			for (int p=0; p<periods.length; p++) {
				EvenlyDiscretizedFunc sdFunc = calcSDFunc(dmFunc, allResiduals, null, null, periods[p]);
				if (periods[p] > 0d)
					sdFunc.setName((float)periods[p]+"s SA");
				else
					sdFunc.setName("PGV");
				funcs.add(sdFunc);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, perCPT.getColor((float)p)));
			}
		}
		
		EvenlyDiscretizedFunc sdFunc = calcSDFunc(dmFunc, allResiduals, null, null, period);
		if (!funcs.isEmpty())
			sdFunc.setName("All");
		funcs.add(sdFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 5f, Color.BLACK));
		
		PlotSpec spec = new PlotSpec(funcs, chars, " ", "Magnitude Uncertainty (+/-)", "Residual Standard Deviation");
		spec.setLegendInset(RectangleAnchor.BOTTOM_RIGHT);
		
		Range xRange = new Range(0d, sdFunc.getMaxX());
		Range yRange = new Range(0d, 1d);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(spec, false, false, xRange, yRange);
		
		PlotUtils.writePlots(outputDir, prefix, gp, 800, 700, true, true, true);
	}
	
	private static EvenlyDiscretizedFunc calcSDFunc(EvenlyDiscretizedFunc dmFunc, List<Residuals> allResiduals,
			Range restrictMagRange, Range restrictDistRange, Double restrictPeriod) {
		EvenlyDiscretizedFunc ret = new EvenlyDiscretizedFunc(dmFunc.getMinX(), dmFunc.getMaxX(), dmFunc.size());
		for (int i=0; i<ret.size(); i++)
			ret.set(i, allResiduals.get(i).calcStdDev(restrictMagRange, restrictDistRange, restrictPeriod));
		return ret;
	}

}
