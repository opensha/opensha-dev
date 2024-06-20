package scratch.kevin.simulators.ruptures;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.RectangleEdge;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.ArbDiscrXYZ_DataSet;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.utils.RSQSimUtils;

import com.google.common.base.Preconditions;
import com.google.common.collect.BiMap;
import com.google.common.collect.HashBiMap;
import com.google.common.primitives.Doubles;

import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class BBPSelfSigmaCheckerboardPlot {

	public static void main(String[] args) throws IOException {
		RSQSimCatalog catalog = Catalogs.BRUCE_5652.instance();
		File bbpDir = new File("/data/kevin/bbp/parallel/2023_08_29-rundir5652-all-m6.5-skipYears5000-noHF-vmLA_BASIN_500-griddedSites");
		File bbpFile = new File(bbpDir, "results_rotD.zip");
		VelocityModel vm = VelocityModel.LA_BASIN_500;
		
		String prefix = "r5652_ask2014";
		
		File outputDir = new File("/home/kevin/Documents/papers/2023_Bruce_GMM_Multifault/figures/self_sigma_checkerboards");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
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
		
		double minMag = Double.POSITIVE_INFINITY;
		double maxMag = 0d;
		double minDist = 1;
		double maxDist = 200;
		
		for (RSQSimEvent e : events) {
			minMag = Math.min(minMag, e.getMagnitude());	
			maxMag = Math.max(maxMag, e.getMagnitude());
		}
		minMag = Math.floor(minMag*10d)/10d;
		maxMag = Math.ceil(maxMag*10d)/10d;
		int numDist = 20;
		int numMag = 15;
		double distSpacing = (Math.log10(maxDist) - Math.log10(minDist)) / (double)numDist;
		double magSpacing = (maxMag - minMag) / (double)numMag;
		
		EvenlyDiscrXYZ_DataSet refXYZ = new EvenlyDiscrXYZ_DataSet(numDist, numMag, Math.log10(minDist)+0.5*distSpacing,
				minMag+0.5*magSpacing, distSpacing, magSpacing);
		
		List<SimpleRakeFilter> indvRakeFilters = new ArrayList<>(3);
		List<String> rakeLabels = new ArrayList<>(4);
		List<String> rakeFileLabels = new ArrayList<>(4);
		
		indvRakeFilters.add(new SimpleRakeFilter(new Range(-180.01, -170), new Range(-10, 10), new Range(170, 180.01)));
		rakeLabels.add("Strike-Slip");
		rakeFileLabels.add("rake_ss");
		
		indvRakeFilters.add(new SimpleRakeFilter(new Range(80, 100)));
		rakeLabels.add("Reverse");
		rakeFileLabels.add("rake_reverse");
		
		indvRakeFilters.add(new SimpleRakeFilter(new Range(-100, -80)));
		rakeLabels.add("Normal");
		rakeFileLabels.add("rake_normal");
		
		List<RakeFilter> rakeFilters = new ArrayList<>(4);
		rakeFilters.addAll(indvRakeFilters);
		rakeFilters.add(new ObliqueFilter(indvRakeFilters));
		rakeLabels.add("Oblique");
		rakeFileLabels.add("rake_oblique");
		
		rakeFilters.add(null);
		rakeLabels.add("All");
		rakeFileLabels.add("all");
		
		boolean rJB = false;
		boolean cumDistRup = true;
		
		List<List<Map<Site, List<RSQSimEvent>>>> siteEventMagDists = new ArrayList<>();
		for (int x=0; x<refXYZ.getNumX(); x++) {
			List<Map<Site, List<RSQSimEvent>>> distLists = new ArrayList<>();
			siteEventMagDists.add(distLists);
			for (int y=0; y<refXYZ.getNumY(); y++)
				distLists.add(new HashMap<>());
		}
		
		System.out.println("Mapping sites and events to dist/mag bins");
		sites.parallelStream().forEach(new Consumer<Site>() {

			@Override
			public void accept(Site site) {
				List<RSQSimEvent> rups = new ArrayList<>(bbpLoader.getRupturesForSite(site));
				List<Double> dists = new ArrayList<>(rups.size());
				for (RSQSimEvent rup : rups) {
					EqkRupture eqkRup = cumDistRup ? catalog.getCumDistCalcRupture(rup) : catalog.getMappedSubSectRupture(rup);
					RuptureSurface surf = eqkRup.getRuptureSurface();
					dists.add(rJB ? surf.getDistanceJB(site.getLocation()) : surf.getDistanceRup(site.getLocation()));
				}
				
				synchronized (siteEventMagDists) {
					for (int i=0; i<rups.size(); i++) {
						RSQSimEvent rup = rups.get(i);
						double dist = dists.get(i);
//						if (dist < minDist || dist > maxDist)
//							continue;
//						int xInd = refXYZ.getXIndex(Math.log10(dist));
						int xInd;
						if (dist < minDist)
							xInd = 0;
						else if (dist > maxDist)
							xInd = numDist-1;
						else
							xInd = refXYZ.getXIndex(Math.log10(dist));
						int yInd = refXYZ.getYIndex(rup.getMagnitude());
						
						List<RSQSimEvent> siteBinList = siteEventMagDists.get(xInd).get(yInd).get(site);
						if (siteBinList == null) {
							siteBinList = new ArrayList<>();
							siteEventMagDists.get(xInd).get(yInd).put(site, siteBinList);
						}
						siteBinList.add(rup);
					}
				}
			}
		});
		
		int minNumSims = 2;
		
		DecimalFormat oDF = new DecimalFormat("0.##");
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.getPlotPrefs().setPlotLabelFontSize(36);
		gp.getPlotPrefs().setAxisLabelFontSize(28);
		gp.getPlotPrefs().setTickLabelFontSize(26);
		Font annFont = new Font(Font.SANS_SERIF, Font.BOLD, 36);
		
		CPT stdDevCPT = GMT_CPT_Files.BLACK_RED_YELLOW_UNIFORM.instance().reverse().rescale(0, 1d);
		stdDevCPT.setNanColor(Color.WHITE);
		
		for (int r=0; r<rakeFilters.size(); r++) {
			RakeFilter filter = rakeFilters.get(r);
			String rakePrefix = prefix+"_"+rakeFileLabels.get(r);
			String rakeLabel = rakeLabels.get(r);
			
			System.out.println("Doing "+rakeLabel);
			
			EvenlyDiscrXYZ_DataSet[] stdDevXYZs = new EvenlyDiscrXYZ_DataSet[periods.length];
			for (int p=0; p<periods.length; p++)
				stdDevXYZs[p] = refXYZ.copy();
			
			for (int x=0; x<refXYZ.getNumX(); x++) {
				for (int y=0; y<refXYZ.getNumY(); y++) {
					Map<Site, List<RSQSimEvent>> siteEvents = siteEventMagDists.get(x).get(y);
					
					List<List<Double>> perGMs = new ArrayList<>(periods.length);
					for (int p=0; p<periods.length; p++)
						perGMs.add(new ArrayList<>());
					int numSims = 0;
					for (Site site : siteEvents.keySet()) {
						for (RSQSimEvent rup : siteEvents.get(site)) {
							if (filter != null && !filter.matches(rup))
								continue;
							numSims++;
							DiscretizedFunc rd50 = bbpLoader.getRotD50(site, rup, 0);
							for (int p=0; p<periods.length; p++) {
								double val;
								if (periods[p] < 0)
									val = bbpLoader.getPGV(site, rup, 0);
								else
									val = rd50.getInterpolatedY(periods[p]);
								perGMs.get(p).add(Math.log(val));
							}
						}
					}
					
					if (numSims < minNumSims) {
						for (int p=0; p<periods.length; p++)
							stdDevXYZs[p].set(x, y, Double.NaN);
						continue;
					}
					for (int p=0; p<periods.length; p++) {
						List<Double> gms = perGMs.get(p);
						double[] array = Doubles.toArray(gms);
						stdDevXYZs[p].set(x, y, Math.sqrt(StatUtils.variance(array)));
					}
				}
			}
			
			List<PlotSpec> specs = new ArrayList<>();
			Range xRange = null;
			Range yRange = null;
			
			for (int p=0; p<periods.length; p++) {
				String perPrefix;
				String perLabel;
				if (periods[p] == -1) {
					perPrefix = "pgv";
					perLabel = "PGV";
				} else {
					Preconditions.checkState(periods[p] > 0d);
					perPrefix = oDF.format(periods[p])+"s";
					perLabel = oDF.format(periods[p])+"s SA";
				}
				
//				ArbDiscrXYZ_DataSet xyz = new ArbDiscrXYZ_DataSet();
//				
//				for (int x=0; x<refXYZ.getNumX(); x++) {
//					for (int y=0; y<refXYZ.getNumY(); y++) {
//						double dist = Math.pow(10, refXYZ.getX(x));
//						double mag = refXYZ.getY(y);
//						xyz.set(dist, mag, stdDevXYZs[p].get(x, y));
//					}
//				}
				
				if (xRange == null) {
//					xRange = new Range(Math.log10(minDist), Math.log10(maxDist));
					xRange = new Range(minDist, maxDist);
					yRange = new Range(minMag, maxMag);
				}
				
				String xLabel = rJB ? "Distance JB  (km)" : "Distance Rup  (km)";
				String zLabel = "Residual Std. Dev.";
				
				String title = filter == null ? " " : rakeLabel;
				XYZPlotSpec spec = new XYZPlotSpec(stdDevXYZs[p], stdDevCPT, title, xLabel, "Magnitude", zLabel);
				spec.setCPTPosition(RectangleEdge.LEFT);
				spec.setCPTTickUnit(0.2);
				
				XYTextAnnotation perAnn = new XYTextAnnotation("  "+perLabel, xRange.getLowerBound(),
						yRange.getLowerBound() + 0.98*yRange.getLength());
				perAnn.setFont(annFont);
				perAnn.setTextAnchor(TextAnchor.TOP_LEFT);
				spec.addPlotAnnotation(perAnn);
				
				specs.add(spec);
			}
			
			int width = 400 + 600*specs.size();
			List<Range> xRanges = new ArrayList<>();
			for (int i=0; i<specs.size(); i++)
				xRanges.add(xRange);
			List<Range> yRanges = List.of(yRange);
			gp.drawGraphPanel(specs, true, false, xRanges, yRanges);
//			PlotUtils.setXTick(gp, 0.25d);
			PlotUtils.setYTick(gp, 0.2d);
			PlotUtils.setSubplotGap(gp, 80);
			
			PlotUtils.writePlots(outputDir, rakePrefix+"_combined", gp, width, 600, true, true, false);
		}
	}
	
	private interface RakeFilter {
		public boolean matches(RSQSimEvent event);
	}

	private static class SimpleRakeFilter implements RakeFilter {
		
		private Range[] ranges;

		private SimpleRakeFilter(Range... ranges) {
			this.ranges = ranges;
		}

		@Override
		public boolean matches(RSQSimEvent event) {
			double rake = RSQSimUtils.getElemAvgRake(event, true);
			
			for (Range range : ranges)
				if (range.contains(rake))
					return true;
			return false;
		}
		
	}
	
	private static class ObliqueFilter implements RakeFilter {
		
		private List<Range> ranges;

		private ObliqueFilter(List<SimpleRakeFilter> rakeFilters) {
			ranges = new ArrayList<>();
			for (SimpleRakeFilter filter : rakeFilters)
				for (Range range : filter.ranges)
					ranges.add(range);
		}

		@Override
		public boolean matches(RSQSimEvent event) {
			double rake = RSQSimUtils.getElemAvgRake(event, true);
			
			if (!Double.isFinite(rake))
				return false;
			
			// see if any of our other ranges contain it
			for (Range range : ranges)
				if (range.contains(rake))
					return false;
			return true;
		}
		
	}

}
