package scratch.kevin.simulators.ruptures;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.RectangleEdge;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.MultiIMR_Averaged_AttenRel;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGV_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.utils.RSQSimUtils;

import com.google.common.base.Preconditions;
import com.google.common.collect.BiMap;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.HashBiMap;
import com.google.common.collect.Table;
import com.google.common.collect.Table.Cell;
import com.google.common.primitives.Doubles;

import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class BBPStandaloneCheckerboardPlots {
	
	private enum PlotType {
		MEAN_SIM("Mean Simulation  (ln)", "sim", getLogValCPT()) {
			@Override
			public double[] calc(Table<Site, RSQSimEvent, double[]> simVals,
					Table<Site, RSQSimEvent, GmmComp[]> gmmComps) {
				List<List<Double>> perVals = null;
				for (Cell<Site, RSQSimEvent, double[]> cell : simVals.cellSet()) {
					double[] vals = cell.getValue();
					if (perVals == null) {
						perVals = new ArrayList<>(vals.length);
						for (int p=0; p<vals.length; p++)
							perVals.add(new ArrayList<>());
					}
					for (int p=0; p<vals.length; p++)
						perVals.get(p).add(vals[p]);
				}
				double[] ret = new double[perVals.size()];
				for (int p=0; p<ret.length; p++)
					ret[p] = perVals.get(p).stream().mapToDouble(d->d).average().getAsDouble();
				return ret;
			}
		},
		MEAN_RESIDUAL("Mean Residual  (ln)", "residuals", getResidualCPT()) {
			@Override
			public double[] calc(Table<Site, RSQSimEvent, double[]> simVals,
					Table<Site, RSQSimEvent, GmmComp[]> gmmComps) {
				List<List<Double>> perResiduals = null;
				for (Cell<Site, RSQSimEvent, double[]> cell : simVals.cellSet()) {
					double[] vals = cell.getValue();
					GmmComp[] comps = gmmComps.get(cell.getRowKey(), cell.getColumnKey());
					if (perResiduals == null) {
						perResiduals = new ArrayList<>(vals.length);
						for (int p=0; p<vals.length; p++)
							perResiduals.add(new ArrayList<>());
					}
					for (int p=0; p<vals.length; p++) {
						double residual = vals[p] - comps[p].logMean;
						perResiduals.get(p).add(residual);
					}
				}
				double[] ret = new double[perResiduals.size()];
				for (int p=0; p<ret.length; p++)
					ret[p] = perResiduals.get(p).stream().mapToDouble(d->d).average().getAsDouble();
				return ret;
			}
		},
		RESIDUAL_TO_EPISTEMIC("Mean Residual w.r.t Epistemic (ln)", "epi_residuals", getResidualCPT()) {
			@Override
			public double[] calc(Table<Site, RSQSimEvent, double[]> simVals,
					Table<Site, RSQSimEvent, GmmComp[]> gmmComps) {
				List<List<Double>> perVals = null;
				for (Cell<Site, RSQSimEvent, double[]> cell : simVals.cellSet()) {
					double[] vals = cell.getValue();
					if (perVals == null) {
						perVals = new ArrayList<>(vals.length);
						for (int p=0; p<vals.length; p++)
							perVals.add(new ArrayList<>());
					}
					for (int p=0; p<vals.length; p++)
						perVals.get(p).add(vals[p]);
				}
				double[] perMeans = new double[perVals.size()];
				for (int p=0; p<perMeans.length; p++)
					perMeans[p] = perVals.get(p).stream().mapToDouble(d->d).average().getAsDouble();
				
				List<List<Double>> compMaxs = null;
				List<List<Double>> compMins = null;
				for (Cell<Site, RSQSimEvent, double[]> cell : simVals.cellSet()) {
					GmmComp[] comps = gmmComps.get(cell.getRowKey(), cell.getColumnKey());
					if (compMaxs == null) {
						compMaxs = new ArrayList<>(comps.length);
						compMins = new ArrayList<>(comps.length);
						for (int p=0; p<comps.length; p++) {
							compMaxs.add(new ArrayList<>());
							compMins.add(new ArrayList<>());
						}
					}
					for (int p=0; p<comps.length; p++) {
						compMaxs.get(p).add(comps[p].epistemicMax);
						compMins.get(p).add(comps[p].epistemicMin);
					}
				}
				double[] compMeanMaxs = new double[compMaxs.size()];
				double[] compMeanMins = new double[compMaxs.size()];
				for (int p=0; p<compMeanMaxs.length; p++) {
					compMeanMaxs[p] = compMaxs.get(p).stream().mapToDouble(d->d).average().getAsDouble();
					compMeanMins[p] = compMins.get(p).stream().mapToDouble(d->d).average().getAsDouble();
				}
				
				double[] ret = new double[perMeans.length];
				for (int p=0; p<perMeans.length; p++) {
					double epistemicMax = compMeanMaxs[p];
					double epistemicMin = compMeanMins[p];
					double val;
					if (!Double.isFinite(epistemicMax) || !Double.isFinite(epistemicMin))
						val = Double.NaN;
					else if (perMeans[p] > epistemicMax)
						val = perMeans[p] - epistemicMax;
					else if (perMeans[p] < epistemicMin)
						val = perMeans[p] - epistemicMin;
					else
						val = 0d;
					ret[p] = val;
				}
				return ret;
			}
		},
		RESIDUAL_STD_DEV("Residual Standard Deviation", "residual_sd", getStdDevCPT()) {
			@Override
			public double[] calc(Table<Site, RSQSimEvent, double[]> simVals,
					Table<Site, RSQSimEvent, GmmComp[]> gmmComps) {
				List<List<Double>> perResiduals = null;
				for (Cell<Site, RSQSimEvent, double[]> cell : simVals.cellSet()) {
					double[] vals = cell.getValue();
					GmmComp[] comps = gmmComps.get(cell.getRowKey(), cell.getColumnKey());
					if (perResiduals == null) {
						perResiduals = new ArrayList<>(vals.length);
						for (int p=0; p<vals.length; p++)
							perResiduals.add(new ArrayList<>());
					}
					for (int p=0; p<vals.length; p++)
						perResiduals.get(p).add(vals[p] - comps[p].logMean);
				}
				double[] ret = new double[perResiduals.size()];
				for (int p=0; p<ret.length; p++)
					ret[p] = Math.sqrt(StatUtils.variance(Doubles.toArray(perResiduals.get(p))));
				return ret;
			}
		},
		BIN_STD_DEV("Bin Standard Deviation", "bin_sd", getStdDevCPT()) {
			@Override
			public double[] calc(Table<Site, RSQSimEvent, double[]> simVals,
					Table<Site, RSQSimEvent, GmmComp[]> gmmComps) {
				List<List<Double>> perVals = null;
				for (Cell<Site, RSQSimEvent, double[]> cell : simVals.cellSet()) {
					double[] vals = cell.getValue();
					if (perVals == null) {
						perVals = new ArrayList<>(vals.length);
						for (int p=0; p<vals.length; p++)
							perVals.add(new ArrayList<>());
					}
					for (int p=0; p<vals.length; p++)
						perVals.get(p).add(vals[p]);
				}
				double[] ret = new double[perVals.size()];
				for (int p=0; p<ret.length; p++)
					ret[p] = Math.sqrt(StatUtils.variance(Doubles.toArray(perVals.get(p))));
				return ret;
			}
		};
		
		private String label;
		private String prefix;
		private CPT cpt;

		private PlotType(String label, String prefix, CPT cpt) {
			this.label = label;
			this.prefix = prefix;
			this.cpt = cpt;
		}
		
		public abstract double[] calc(Table<Site, RSQSimEvent, double[]> simVals, Table<Site, RSQSimEvent, GmmComp[]> gmmComps);
	}
	
	private static CPT getStdDevCPT() {
		CPT stdDevCPT;
		try {
			stdDevCPT = GMT_CPT_Files.BLACK_RED_YELLOW_UNIFORM.instance().reverse().rescale(0, 1d);
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		stdDevCPT.setNanColor(Color.WHITE);
		stdDevCPT.setPreferredTickInterval(0.2);
		return stdDevCPT;
	}
	
	private static CPT getResidualCPT() {
		CPT residualCPT;
		try {
//			residualCPT = GMT_CPT_Files.GMT_POLAR.instance().rescale(-1.5, 1.5);
			residualCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-1.5, 1.5);
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		residualCPT.setNanColor(Color.WHITE);
		residualCPT.setPreferredTickInterval(0.5);
		return residualCPT;
	}
	
	private static CPT getLogValCPT() {
		CPT cpt;
		try {
			cpt = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-6, 4);
		} catch (IOException e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		cpt.setNanColor(Color.WHITE);
		cpt.setPreferredTickInterval(1d);
		return cpt;
	}

	public static void main(String[] args) throws IOException {
		RSQSimCatalog catalog = Catalogs.BRUCE_5652.instance();
		File bbpDir = new File("/data/kevin/bbp/parallel/2023_08_29-rundir5652-all-m6.5-skipYears5000-noHF-vmLA_BASIN_500-griddedSites");
		File bbpFile = new File(bbpDir, "results_rotD.zip");
		String prefix = "r5652";
		VelocityModel vm = VelocityModel.LA_BASIN_500;
		
//		PlotType[] types = PlotType.values();
//		AttenRelRef gmmRef = AttenRelRef.NGAWest_2014_AVG_NOIDRISS;
//		prefix += "_ngaw2";
		
		PlotType[] types = {
				PlotType.MEAN_SIM,
				PlotType.MEAN_RESIDUAL,
				PlotType.RESIDUAL_STD_DEV,
				PlotType.BIN_STD_DEV
		};
		AttenRelRef gmmRef = AttenRelRef.ASK_2014;
		prefix += "_ask2014";
		
		File outputDir = new File("/home/kevin/Documents/papers/2023_Bruce_GMM_Multifault/figures/checkerboards");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		outputDir = new File(outputDir, prefix);
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
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
		
		rakeFilters.add(0, null);
		rakeLabels.add(0, "All");
		rakeFileLabels.add(0, "all");
		
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
		sites.parallelStream().forEach(site -> {
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
//					if (dist < minDist || dist > maxDist)
//						continue;
//					int xInd = refXYZ.getXIndex(Math.log10(dist));
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
				System.out.print(".");
			}
		});
		System.out.println();
		
		System.out.println("Calculating GMM values for "+sites.size()+" sites");
		ArrayDeque<ScalarIMR> gmmDeque = new ArrayDeque<>();
		Table<Site, RSQSimEvent, GmmComp[]> gmmCompTable = HashBasedTable.create();
		sites.parallelStream().forEach(site -> {
			ScalarIMR gmm = null;
			synchronized (gmmDeque) {
				if (!gmmDeque.isEmpty())
					gmm = gmmDeque.pop();
			}
			if (gmm == null)
				gmm = gmmRef.get();
			
			gmm.setSite(site);
			
			List<RSQSimEvent> rups = new ArrayList<>(bbpLoader.getRupturesForSite(site));
			List<GmmComp[]> comps = new ArrayList<>(rups.size());
			for (RSQSimEvent rup : rups) {
				EqkRupture eqkRup = cumDistRup ? catalog.getCumDistCalcRupture(rup) : catalog.getMappedSubSectRupture(rup);
				gmm.setEqkRupture(eqkRup);
				
				GmmComp[] comp = new GmmComp[periods.length];
				for (int p=0; p<periods.length; p++) {
					if (periods[p] < 0d) {
						gmm.setIntensityMeasure(PGV_Param.NAME);
					} else if (periods[p] == 0d) {
						gmm.setIntensityMeasure(PGA_Param.NAME);
					} else {
						gmm.setIntensityMeasure(SA_Param.NAME);
						SA_Param.setPeriodInSA_Param(gmm.getIntensityMeasure(), periods[p]);
					}
					
					double mean = gmm.getMean();
					double stdDev = gmm.getStdDev();
					if (gmm instanceof MultiIMR_Averaged_AttenRel) {
						double epiMin = mean;
						double epiMax = mean;
						for (ScalarIMR epiGMM : ((MultiIMR_Averaged_AttenRel)gmm).getIMRs()) {
							double myMean = epiGMM.getMean();
							epiMin = Math.min(epiMin, myMean);
							epiMax = Math.max(epiMax, myMean);
						}
						Preconditions.checkState(epiMax > epiMin);
						comp[p] = new GmmComp(mean, epiMin, epiMax, stdDev);
					} else {
						comp[p] = new GmmComp(mean, stdDev);
					}
				}
				comps.add(comp);
			}
			
			Preconditions.checkState(comps.size() == rups.size());
			
			synchronized (siteEventMagDists) {
				for (int i=0; i<rups.size(); i++)
					gmmCompTable.put(site, rups.get(i), comps.get(i));
				System.out.print(".");
			}
			
			synchronized (gmmDeque) {
				gmmDeque.push(gmm);
			}
		});
		System.out.println();
		
		System.out.println("Calculated "+gmmCompTable.size()+" GMM Comps");
		
		int minNumSimsPerBin = 5;
		
		DecimalFormat oDF = new DecimalFormat("0.##");
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.getPlotPrefs().setPlotLabelFontSize(36);
		gp.getPlotPrefs().setAxisLabelFontSize(28);
		gp.getPlotPrefs().setTickLabelFontSize(26);
		Font annFont = new Font(Font.SANS_SERIF, Font.BOLD, 36);
		
		List<String> lines = new ArrayList<>();
		
		lines.add("# Checkerboard Plots");
		lines.add("");
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		for (int r=0; r<rakeFilters.size(); r++) {
			RakeFilter filter = rakeFilters.get(r);
			String rakePrefix = prefix+"_"+rakeFileLabels.get(r);
			String rakeLabel = rakeLabels.get(r);
			
			System.out.println("Doing "+rakeLabel);
			
			lines.add("## "+rakeLabel);
			lines.add(topLink); lines.add("");
			
			TableBuilder table = MarkdownUtils.tableBuilder();
			
			for (PlotType type : types) {
				System.out.println("\t"+type);
				EvenlyDiscrXYZ_DataSet[] xyzs = new EvenlyDiscrXYZ_DataSet[periods.length];
				for (int p=0; p<periods.length; p++)
					xyzs[p] = refXYZ.copy();
				
				for (int x=0; x<refXYZ.getNumX(); x++) {
					for (int y=0; y<refXYZ.getNumY(); y++) {
						Map<Site, List<RSQSimEvent>> siteEvents = siteEventMagDists.get(x).get(y);
						
						Table<Site, RSQSimEvent, double[]> simVals = HashBasedTable.create();
						
						int numSims = 0;
						for (Site site : siteEvents.keySet()) {
							for (RSQSimEvent rup : siteEvents.get(site)) {
								if (filter != null && !filter.matches(rup))
									continue;
								numSims++;
								Preconditions.checkState(bbpLoader.getNumSimulations(site, rup) == 1);
								DiscretizedFunc rd50 = bbpLoader.getRotD50(site, rup, 0);
								double[] vals = new double[periods.length];
								for (int p=0; p<periods.length; p++) {
									double val;
									if (periods[p] < 0)
										val = bbpLoader.getPGV(site, rup, 0);
									else if (periods[p] == 0)
										val = bbpLoader.getPGA(site, rup, 0);
									else
										val = rd50.getInterpolatedY(periods[p]);
									vals[p] = Math.log(val);
								}
								simVals.put(site, rup, vals);
							}
						}
						
						if (numSims < minNumSimsPerBin) {
							for (int p=0; p<periods.length; p++)
								xyzs[p].set(x, y, Double.NaN);
							continue;
						}
						double[] perVals = type.calc(simVals, gmmCompTable);
						for (int p=0; p<periods.length; p++)
							xyzs[p].set(x, y, perVals[p]);
					}
				}
				
				List<PlotSpec> specs = new ArrayList<>();
				Range xRange = null;
				Range yRange = null;
				
				String plotPrefix = rakePrefix+"_"+type.prefix+"_combined";
				
				String links = "[PDF Link]("+resourcesDir.getName()+"/"+plotPrefix+".pdf), Individual Periods:";
				
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
					
//					ArbDiscrXYZ_DataSet xyz = new ArbDiscrXYZ_DataSet();
//					
//					for (int x=0; x<refXYZ.getNumX(); x++) {
//						for (int y=0; y<refXYZ.getNumY(); y++) {
//							double dist = Math.pow(10, refXYZ.getX(x));
//							double mag = refXYZ.getY(y);
//							xyz.set(dist, mag, stdDevXYZs[p].get(x, y));
//						}
//					}
					
					if (xRange == null) {
//						xRange = new Range(Math.log10(minDist), Math.log10(maxDist));
						xRange = new Range(minDist, maxDist);
						yRange = new Range(minMag, maxMag);
					}
					
					String xLabel = rJB ? "Distance JB  (km)" : "Distance Rup  (km)";
					String zLabel = type.label;
					
					String title = filter == null ? " " : rakeLabel;
					XYZPlotSpec spec = new XYZPlotSpec(xyzs[p], type.cpt, title, xLabel, "Magnitude", zLabel);
					spec.setCPTPosition(RectangleEdge.LEFT);
					
					XYTextAnnotation perAnn = new XYTextAnnotation("  "+perLabel, xRange.getLowerBound(),
							yRange.getLowerBound() + 0.98*yRange.getLength());
					perAnn.setFont(annFont);
					perAnn.setTextAnchor(TextAnchor.TOP_LEFT);
					spec.addPlotAnnotation(perAnn);
					
					gp.drawGraphPanel(spec, true, false, xRange, yRange);
//					PlotUtils.setXTick(gp, 0.25d);
					PlotUtils.setYTick(gp, 0.2d);
					
					PlotUtils.writePlots(resourcesDir, plotPrefix+"_"+perPrefix, gp, 800, 600, true, true, false);
					
					links += " ["+perLabel+"]("+resourcesDir.getName()+"/"+plotPrefix+"_"+perPrefix+".png)";
					
					specs.add(spec);
				}
				
				int width = 400 + 600*specs.size();
				List<Range> xRanges = new ArrayList<>();
				for (int i=0; i<specs.size(); i++)
					xRanges.add(xRange);
				List<Range> yRanges = List.of(yRange);
				gp.drawGraphPanel(specs, true, false, xRanges, yRanges);
//				PlotUtils.setXTick(gp, 0.25d);
				PlotUtils.setYTick(gp, 0.2d);
				PlotUtils.setSubplotGap(gp, 80);
				
				PlotUtils.writePlots(resourcesDir, plotPrefix, gp, width, 600, true, true, false);
				
				table.addLine(MarkdownUtils.boldCentered(type.label));
				table.addLine("![Checkerboards]("+resourcesDir.getName()+"/"+plotPrefix+".png)");
				table.addLine(links);
			}
			lines.addAll(table.build());
			lines.add("");
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
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
	
	private static class GmmComp {
		public final double logMean;
		public final double epistemicMin;
		public final double epistemicMax;
		public final double stdDev;
		
		private GmmComp(double logMean, double stdDev) { 
			this(logMean, Double.NaN, Double.NaN, stdDev);
		}
		
		private GmmComp(double logMean, double epistemicMin, double epistemicMax, double stdDev) { 
			super();
			this.logMean = logMean;
			this.epistemicMin = epistemicMin;
			this.epistemicMax = epistemicMax;
			this.stdDev = stdDev;
		}
	}

}
