package scratch.kevin.ucerf3.etas;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.apache.commons.math3.stat.StatUtils;
import org.dom4j.DocumentException;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.function.UncertainArbDiscDataset;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotElement;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.FileNameComparator;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.nshmp2.erf.source.PointSource13b;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.cybershake.etas.ETASModProbConfig.ETAS_CyberShake_Scenarios;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.griddedSeis.Point2Vert_FaultPoisSource;
import org.opensha.sha.imr.AttenRelRef;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.ETAS_MultiSimAnalysisTools;
import scratch.UCERF3.erf.ETAS.ETAS_SimAnalysisTools;
import scratch.UCERF3.erf.ETAS.ETAS_Simulator.TestScenario;
import scratch.UCERF3.erf.ETAS.ETAS_Utils;
import scratch.UCERF3.erf.mean.TrueMeanBuilder;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;
import scratch.UCERF3.griddedSeismicity.AbstractGridSourceProvider;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.LastEventData;
import scratch.kevin.ucerf3.eal.UCERF3_BranchAvgLossFetcher;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Table;
import com.google.common.io.CharStreams;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;

public class ETAS_CatalogEALCalculator {
	
	static final boolean conf_68 = true; // otherwise 95
	static final String conf_lower_str;
	static final String conf_upper_str;
	static {
		if (conf_68) {
			conf_lower_str = "Lower 68%";
			conf_upper_str = "Upper 68%";
		} else {
			conf_lower_str = "Lower 95%";
			conf_upper_str = "Upper 95%";
		}
	}
	
	private UCERF3_BranchAvgLossFetcher fetcher;
	private List<List<ETAS_EqkRupture>> catalogs;
	private long startTime;
	private FaultModels fm;
	
	// for getting fss index from Ned's "Nth" index
	private FaultSystemSolutionERF erf;
	private FaultSystemSolution meanSol;
	
	private ETAS_EqkRupture triggerRup;
	
	private static boolean rup_mean_loss = true; // otherwise propagate loss distribution
	
	private static final double outside_region_dist_tol = 10d; // km
	
	private boolean triggeredOnly = false;
	
	private static int id_for_scenario = 0;
	
	private static List<List<ETAS_EqkRupture>> loadCatalogs(File resultsBinFile) throws IOException {
		Preconditions.checkArgument(resultsBinFile.exists(), "catalog file doesn't exist");
		
		return loadCatalogs(resultsBinFile, AbstractGridSourceProvider.SOURCE_MIN_MAG_CUTOFF-0.05);
	}
	
	public ETAS_CatalogEALCalculator(UCERF3_BranchAvgLossFetcher fetcher, FaultSystemSolution meanSol,
			FaultModels fm, File resultsBinFile) throws IOException, DocumentException {
		this(fetcher, meanSol, fm, loadCatalogs(resultsBinFile));
	}
	
	public ETAS_CatalogEALCalculator(UCERF3_BranchAvgLossFetcher fetcher, FaultSystemSolution meanSol,
			FaultModels fm, List<List<ETAS_EqkRupture>> catalogs) throws IOException, DocumentException {
		this.fetcher = fetcher;
		this.fm = fm;
		
		this.catalogs = catalogs;
		startTime = Long.MAX_VALUE;
		for (List<ETAS_EqkRupture> catalog : catalogs) {
			if (!catalog.isEmpty()) {
				if (catalog.get(0).getOriginTime() < startTime)
					startTime = catalog.get(0).getOriginTime();
			}
		}
		
		Preconditions.checkState(!catalogs.isEmpty(), "No catalogs loaded!");
		System.out.println("Loaded "+catalogs.size()+" catalogs");
		
		LastEventData.populateSubSects(meanSol.getRupSet().getFaultSectionDataList(), LastEventData.load());
		this.meanSol = meanSol;
		System.out.println("Loading ERF");
		double origMinMag = AbstractGridSourceProvider.SOURCE_MIN_MAG_CUTOFF;
		AbstractGridSourceProvider.SOURCE_MIN_MAG_CUTOFF = 2.55;
		erf = MPJ_ETAS_Simulator.buildERF(meanSol, false, 1d);
		erf.updateForecast();
		AbstractGridSourceProvider.SOURCE_MIN_MAG_CUTOFF = origMinMag;
		System.out.println("Done loading ERF");
	}
	
//	static List<List<ETAS_EqkRupture>> loadCatalogs(File[] etasCatalogsDirs, double minGriddedMag) throws IOException {
//		return loadCatalogs(etasCatalogsDirs, minGriddedMag, null);
//	}
	
	public void setTriggeredOnly(boolean triggeredOnly) {
		this.triggeredOnly = triggeredOnly;
	}
	
	/**
	 * Load ETAS catalogs from file.
	 * @param etasCatalogsDirs
	 * @param catalogDirs if not null, sub directories will be added to this list so that you can keep track
	 * of which directory each catalog came from
	 * @return
	 * @throws IOException
	 */
	static List<List<ETAS_EqkRupture>> loadCatalogs(File resultsBinFile, double minGriddedMag) throws IOException {
		int numEmpty = 0;
		
		List<List<ETAS_EqkRupture>> catalogs = ETAS_CatalogIO.loadCatalogsBinary(resultsBinFile, minGriddedMag);
		for (List<ETAS_EqkRupture> catalog : catalogs)
			if (catalog.isEmpty())
				numEmpty++;
		
		System.out.println(numEmpty+"/"+catalogs.size()+" catalogs are empty "
				+ "(including only fault and gridded above "+minGriddedMag+")");
		return catalogs;
	}
	
	static List<List<ETAS_EqkRupture>> loadCatalogsZip(File zipFile, double minGriddedMag,
			List<String> catalogNames) throws IOException {
		ZipFile zip = new ZipFile(zipFile);
		
		List<List<ETAS_EqkRupture>> catalogs = Lists.newArrayList();
		
		List<? extends ZipEntry> entries = Collections.list(zip.entries());
		// sort for constant ordering
		Collections.sort(entries, new Comparator<ZipEntry>() {

			@Override
			public int compare(ZipEntry o1, ZipEntry o2) {
				return o1.getName().compareTo(o2.getName());
			}
		});
		
		for (ZipEntry entry : entries) {
			if (!entry.isDirectory())
				continue;
//			System.out.println(entry.getName());
			String subEntryName = entry.getName()+"simulatedEvents.txt";
			ZipEntry catEntry = zip.getEntry(subEntryName);
			String infoEntryName = entry.getName()+"infoString.txt";
			ZipEntry infoEntry = zip.getEntry(infoEntryName);
			if (catEntry == null || infoEntry == null)
				continue;
			
			// make sure it's actually done
			BufferedReader reader = new BufferedReader(new InputStreamReader(zip.getInputStream(infoEntry)));
			
			boolean done = false;
			for (String line : CharStreams.readLines(reader)) {
				if (line.contains("Total num ruptures: ")) {
					done = true;
					break;
				}
			}
			if (!done)
				continue;
//			System.out.println("Loading "+catEntry.getName());
			
			List<ETAS_EqkRupture> catalog;
			try {
				catalog = ETAS_CatalogIO.loadCatalog(zip.getInputStream(catEntry), minGriddedMag);
			} catch (Exception e) {
				continue;
			}
			
			if (catalogNames != null)
				catalogNames.add(entry.getName());
			catalogs.add(catalog);
		}
		
		return catalogs;
	}
	
	/**
	 * 
	 * @param attenRelRef
	 * @return loss dists for each catalog. it is a distribution because fault based ruptures have loss distributions
	 * rather than a single value.
	 * @throws IOException
	 */
	public synchronized List<DiscretizedFunc> getLossDists(AttenRelRef attenRelRef, double xAxisScale, double durationYears) throws IOException {
		Map<Double, List<DiscretizedFunc>> vals = getLossDists(attenRelRef, xAxisScale, new double[] {durationYears});
		return vals.get(durationYears);
	}
	
	public synchronized Map<Double, List<DiscretizedFunc>> getLossDists(
			AttenRelRef attenRelRef, double xAxisScale, double[] durations) throws IOException {
		return getLossDists(attenRelRef, xAxisScale, durations, false);
	}
	
	public static final int max_all_sub_durations_n = 1000000;
	
	private double getRoundedMaxCatalogDiration() {
		double maxCatalogDuration = 0d;
		for (List<ETAS_EqkRupture> catalog : catalogs)
			maxCatalogDuration = Math.max(maxCatalogDuration, ETAS_MultiSimAnalysisTools.calcDurationYears(catalog));
		// round to nearest year
		maxCatalogDuration = Math.round(maxCatalogDuration);
		return maxCatalogDuration;
	}
	
	public synchronized Map<Double, List<DiscretizedFunc>> getLossDists(
			AttenRelRef attenRelRef, double xAxisScale, double[] durations, boolean allSubDurations)
					throws IOException {
		// conditional loss distributions (x=loss, y=weight) for each rupture
		DiscretizedFunc[] condLossDists = fetcher.getFaultLosses(attenRelRef, fm, true);
		// mag/loss distributions at each grid node (x=mag, y=loss)
		GriddedRegion region = new CaliforniaRegions.RELM_TESTING_GRIDDED();
		DiscretizedFunc[] griddedMagLossDists = fetcher.getGriddedMagLossDists(
				attenRelRef, region);
		
		double maxCatalogDuration = getRoundedMaxCatalogDiration();
		System.out.println("Max duration: "+maxCatalogDuration+" yrs");
		
		Map<Double, List<DiscretizedFunc>> distsMap = Maps.newHashMap();
		
		for (double durationYears : durations) {
			long maxTime = Long.MAX_VALUE;
			if (durationYears > 00 && !Double.isInfinite(durationYears))
				maxTime = startTime + (long)(durationYears*ProbabilityModelsCalc.MILLISEC_PER_YEAR);
			
			List<DiscretizedFunc> catalogDists = Lists.newArrayList();
			
			List<List<ETAS_EqkRupture>> catalogs = this.catalogs;
			if (allSubDurations && maxCatalogDuration > durationYears*1.1) {
				int numPer = (int)(maxCatalogDuration/durationYears);
				System.out.println(numPer+" sub catalogs for "+durationYears+" yr");
				if (numPer*catalogs.size() > max_all_sub_durations_n) {
					numPer = max_all_sub_durations_n/catalogs.size();
					System.out.println("Too many sub catalogs, trimming to "+numPer+" each");
				}
				Preconditions.checkState(numPer > 1, "bad numPer=%, durationYears=%s, maxCatDuration=%s",
						numPer, durationYears, maxCatalogDuration);
				long millisEach = (long)(durationYears*ProbabilityModelsCalc.MILLISEC_PER_YEAR);
				long[] maxTimes = new long[numPer];
				for (int i=0; i<numPer; i++) {
					if (i == 0)
						maxTimes[i] = startTime + millisEach;
					else
						maxTimes[i] = maxTimes[i-1] + millisEach;
				}
				List<List<ETAS_EqkRupture>> subCatalogs = Lists.newArrayList();
				for (List<ETAS_EqkRupture> catalog : catalogs) {
					int curIndex = 0;
					List<List<ETAS_EqkRupture>> mySubCatalogs = Lists.newArrayList();
					for (int i=0; i<numPer; i++)
						mySubCatalogs.add(null);
					catalogLoop:
					for (ETAS_EqkRupture rup : catalog) {
						while (rup.getOriginTime() > maxTimes[curIndex]) {
							curIndex++;
							if (curIndex == maxTimes.length)
								break catalogLoop;
						}
						if (mySubCatalogs.get(curIndex) == null)
							mySubCatalogs.set(curIndex, new ArrayList<ETAS_EqkRupture>());
						mySubCatalogs.get(curIndex).add(rup);
					}
					subCatalogs.addAll(mySubCatalogs);
				}
				catalogs = subCatalogs;
			}
			
			for (int i = 0; i < catalogs.size(); i++) {
				List<ETAS_EqkRupture> catalog = catalogs.get(i);
				if (catalog == null) {
					catalogDists.add(new LightFixedXFunc(new double[] {0d}, new double[] {1d}));
					continue;
				}
				List<Double> singleLosses = Lists.newArrayList();
				List<DiscretizedFunc> lossDists = Lists.newArrayList();
				
				if (triggeredOnly)
					catalog = ETAS_SimAnalysisTools.getChildrenFromCatalog(catalog, id_for_scenario);
				
				for (ETAS_EqkRupture rup : catalog) {
					if (!allSubDurations && rup.getOriginTime() > maxTime)
						break;
					int fssIndex = getFSSIndex(rup);
					
					double mag = rup.getMag();
					
					if (fssIndex >= 0) {
						// fault based source
						double solMag = meanSol.getRupSet().getMagForRup(fssIndex);
						Preconditions.checkState((float)mag == (float)solMag, "Bad fault mag! %s != %s", mag, solMag);
						if (condLossDists[fssIndex].size() == 0)
							continue;
						lossDists.add(condLossDists[fssIndex]);
						// make sure weights sum to 1
						double sumY = 0;
						for (Point2D pt : condLossDists[fssIndex])
							sumY += pt.getY();
						Preconditions.checkState((float)sumY == 1f, "rup losses don't sum to 1: "+(float)sumY);
					} else {
						// grid source
						double loss = calcGridSourceLoss(rup, region, griddedMagLossDists, "Catalog "+i);
						// single loss value with weight=1
						singleLosses.add(loss);
					}
				}
				
				// first sum up all single losses (easy)
				double totSingleLosses = 0d;
				for (double loss : singleLosses)
					totSingleLosses += loss;
				
				if (rup_mean_loss) {
					for (DiscretizedFunc lossDist : lossDists) {
						double loss = 0;
						double sumWeight = 0;
						for (Point2D pt : lossDist) {
							sumWeight += pt.getY();
							loss += pt.getX()*pt.getY();
						}
						Preconditions.checkState((float)sumWeight == 1f, "Weights don't sum to 1: "+(float)sumWeight);
						totSingleLosses += loss;
					}
				}
				
				ArbitrarilyDiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
				if (lossDists.isEmpty() || rup_mean_loss) {
					// only point sources
					func.set(xAxisScale*totSingleLosses, 1d);
				} else {
					// calculate expected number of loss dists for verification
					int expectedNum = 1;
					for (DiscretizedFunc lossDist : lossDists)
						expectedNum *= lossDist.size();
					
					List<LossChain> lossChains = getLossChains(totSingleLosses, lossDists);
					Preconditions.checkState(lossChains.size() == expectedNum,
							"expected "+expectedNum+" chains, got "+lossChains.size());
					
					double sumWeight = 0d;
					for (LossChain chain : lossChains) {
						double weight = chain.weight;
						double loss = chain.totLoss;
						sumWeight += weight;
						int xInd = UCERF3_BranchAvgLossFetcher.getMatchingXIndexFloatPrecision(loss, func);
						if (xInd < 0)
							func.set(xAxisScale*loss, weight);
						else
							func.set(xAxisScale*loss, weight + func.getY(xInd));
					}
					Preconditions.checkState((float)sumWeight == 1f,
							"chain weights don't sum to 1: "+sumWeight+" ("+lossChains.size()+" chains)");
				}
				
//				double meanLoss = 0d;
//				for (Point2D pt : func)
//					meanLoss += pt.getX()*pt.getY();
				
				catalogDists.add(new LightFixedXFunc(func));
			}
			
			if (triggerRup != null) {
				int fssIndex = getFSSIndex(triggerRup);
				
				double mag = triggerRup.getMag();
				
				double triggerLoss;
				
				if (fssIndex >= 0) {
					// fault based source
					Preconditions.checkState((float)mag == (float)meanSol.getRupSet().getMagForRup(fssIndex));
					// make sure weights sum to 1
					double meanLoss = 0;
					double sumY = 0;
					for (Point2D pt : condLossDists[fssIndex]) {
						sumY += pt.getY();
						meanLoss += pt.getX()*pt.getY();
					}
					Preconditions.checkState((float)sumY == 1f || condLossDists[fssIndex].size()==0,
							"rup losses don't sum to 1: "+(float)sumY+" ("+condLossDists[fssIndex].size()+"");
					triggerLoss = meanLoss;
				} else {
					// grid source
					triggerLoss = calcGridSourceLoss(triggerRup, region, griddedMagLossDists, "TRIGGER");
				}
				triggerLoss *= xAxisScale;
				System.out.println("Trigger M"+(float)mag+" rupture loss: "+triggerLoss);
			}
			
			distsMap.put(durationYears, catalogDists);
		}
		
		return distsMap;
	}
	
	static int calcNodeIndex(ETAS_EqkRupture rup, GriddedRegion region) {
		Location loc = rup.getHypocenterLocation();
		int nodeIndex = region.indexForLocation(loc);
		if (nodeIndex < 0) {
			// find closest
			Location closestLoc = null;
			double closestDist = Double.POSITIVE_INFINITY;
			int closestIndex = -1;
			LocationList nodeList = region.getNodeList();
			for (int j = 0; j < nodeList.size(); j++) {
				Location nodeLoc = nodeList.get(j);
				double dist = LocationUtils.horzDistanceFast(nodeLoc, loc);
				if (dist < closestDist) {
					closestDist = dist;
					closestLoc = nodeLoc;
					closestIndex = j;
				}
			}
			if (closestDist <= outside_region_dist_tol) {
				System.out.println("Location ("+loc+") outside region, mapped to node "
						+(float)closestDist+" km away");
				nodeIndex = closestIndex;
			} else {
				System.out.println("Location ("+loc+") outside region, closest node is "
						+(float)closestDist+" km away. Did not map.");
			}
		}
		return nodeIndex;
	}
	
	double calcGridSourceLoss(ETAS_EqkRupture rup, GriddedRegion region,
			DiscretizedFunc[] griddedMagLossDists, String catName) {
		double mag = rup.getMag();
		if ((float)mag < (float)AbstractGridSourceProvider.SOURCE_MIN_MAG_CUTOFF)
			// below our min mag
			return 0;
		Location loc = rup.getHypocenterLocation();
		int nodeIndex = calcNodeIndex(rup, region);
		Preconditions.checkState(nodeIndex >= 0, "Node not found for loc: "+loc);
		DiscretizedFunc magLossDist = griddedMagLossDists[nodeIndex];
		int magIndex = UCERF3_BranchAvgLossFetcher.getMatchingXIndexFloatPrecision(mag, magLossDist);
		if (magIndex < 0) {
			List<Float> xVals = Lists.newArrayList();
			List<Float> yVals = Lists.newArrayList();
			for (Point2D pt : magLossDist) {
				xVals.add((float)pt.getX());
				yVals.add((float)pt.getY());
			}
			System.out.println("Mag "+mag+" not found in cat "+catName+" for loc="+loc+", gridLoc="
				+region.getNodeList().get(nodeIndex)+"."
				+"\n\tNode mags:\t"+Joiner.on("\t").join(xVals)
				+"\n\tNode loss:\t"+Joiner.on("\t").join(yVals));
			if (rup == triggerRup) {
				// interpolate
				if (mag > magLossDist.getMaxX()) {
					double binWidth;
					if (magLossDist.size() >= 2)
						binWidth = magLossDist.getX(magLossDist.size()-1)
							- magLossDist.getX(magLossDist.size()-2);
					else
						binWidth = 0d;
					if (mag < magLossDist.getMaxX()+binWidth) {
						System.out.println("Above but within last bin, returning that");
						return magLossDist.getY(magLossDist.size()-1);
					} else {
						System.out.println("Above mag loss dist, returning zero");
						return 0d;
					}
				}
				System.out.println("Interpolating losses for trigger rupture");
				return magLossDist.getInterpolatedY(mag);
			}
			if (rup.getNthERF_Index() >= 0) {
				int nth = rup.getNthERF_Index();
				int sourceID = erf.getSrcIndexForNthRup(nth);
				System.out.println("Using Nth ERF Index ("+nth+", src="+sourceID+") instead of catalog location");
				Preconditions.checkState(sourceID >= 0);
				ProbEqkSource source = erf.getSource(sourceID);
				Location sourceLoc = null;
				if (source instanceof PointSource13b)
					sourceLoc = ((PointSource13b)source).getLocation();
				else if (source instanceof Point2Vert_FaultPoisSource)
					sourceLoc = ((Point2Vert_FaultPoisSource)source).getLocation();
				else
					System.out.println("Unkown grid source type, skipping: "+source.getClass().getName());
				if (sourceLoc != null) {
					double dist = LocationUtils.horzDistance(sourceLoc, loc);
					Preconditions.checkState(dist < 15d, "Source location via Nth index is wrong! "
							+dist+" km away, catLoc="+loc+", sourceLoc="+sourceLoc);
					nodeIndex = region.indexForLocation(sourceLoc);
					magLossDist = griddedMagLossDists[nodeIndex];
					magIndex = UCERF3_BranchAvgLossFetcher.getMatchingXIndexFloatPrecision(mag, magLossDist);
					if (magIndex < 0) {
						xVals = Lists.newArrayList();
						yVals = Lists.newArrayList();
						for (Point2D pt : magLossDist) {
							xVals.add((float)pt.getX());
							yVals.add((float)pt.getY());
						}
						System.out.println("Same problem with Nth rup node ("+sourceLoc+"), skipping"
								+"\n\tNode mags:\t"+Joiner.on("\t").join(xVals)
								+"\n\tNode loss:\t"+Joiner.on("\t").join(yVals));
						return 0;
					}
					double loss = magLossDist.getY(magIndex);
					return loss;
				}
			}
//			System.out.println("Dist x vals: "+Joiner.on(",").join(xVals));
//			System.out.println("Node loc: "+region.getNodeList().get(nodeIndex));
			return 0; // TODO figure out what's going on and don't skip!
		}
		Preconditions.checkState(magIndex >= 0, "Mag index not found for grid node. mag="+mag+", loc="+loc);
		double loss = magLossDist.getY(magIndex);
		return loss;
	}
	
	private int getFSSIndex(ETAS_EqkRupture rup) {
		return ETAS_SimAnalysisTools.getFSSIndex(rup, erf);
	}
	
	public void setTriggerGridRup(Location loc, double mag) {
		triggerRup = new ETAS_EqkRupture();
		triggerRup.setPointSurface(loc);
		triggerRup.setHypocenterLocation(loc);
		triggerRup.setMag(mag);
		triggerRup.setNthERF_Index(Integer.MAX_VALUE);
	}
	
	public void setTriggerFaultRup(int fssRupID) {
		triggerRup = new ETAS_EqkRupture();
		triggerRup.setMag(meanSol.getRupSet().getMagForRup(fssRupID));
		triggerRup.setRuptureSurface(meanSol.getRupSet().getSurfaceForRupupture(fssRupID, 1d, false));
		triggerRup.setNthERF_Index(erf.get_nthRupIndicesForSource(erf.getSrcIndexForFltSysRup(fssRupID))[0]);
	}
	
	private List<LossChain> getLossChains(double totSingleLosses, List<DiscretizedFunc> distsToProcess) {
		return getLossChains(new LossChain(totSingleLosses), distsToProcess, 0);
	}
	
	private List<LossChain> getLossChains(LossChain prevChain,
			List<DiscretizedFunc> distsToProcess, int distIndex) {
		if (distIndex == distsToProcess.size())
			// we're done
			return Lists.newArrayList(prevChain);
		DiscretizedFunc lossDist = distsToProcess.get(distIndex);
		
		List<LossChain> ret = Lists.newArrayList();
		
		for (Point2D pt : lossDist) {
			LossChain newChain = prevChain.copy();
			newChain.add(pt.getX(), pt.getY());
			ret.addAll(getLossChains(newChain, distsToProcess, distIndex+1));
		}
		
		return ret;
	}
	
	private class LossChain {
		private double totLoss;
		private double weight;
		
		public LossChain(double totSingleLosses) {
			this(totSingleLosses, 1d);
		}
		
		private LossChain(double totLoss, double weight) {
			this.totLoss = totLoss;
			this.weight = weight;
		}
		
		public void add(double loss, double weight) {
			this.totLoss += loss;
			this.weight *= weight;
		}
		
		public LossChain copy() {
			return new LossChain(totLoss, weight);
		}
	}
	
	public static Map<Double, HistogramFunction> getLossHist(Map<Double, List<DiscretizedFunc>> catalogDists,
			double delta, boolean isLog10) {
		double maxDuration = 0d;
		List<DiscretizedFunc> longestDists = null;
		for (double duration : catalogDists.keySet()) {
			if (duration > maxDuration) {
				maxDuration = duration;
				longestDists = catalogDists.get(duration);
			}
		}
		int num = -1;
		Range range = null;
		while (num == -1 || num == 1) {
			range = calcSmartHistRange(longestDists, delta, isLog10);
			if (!isLog10)
				range = new Range(delta*0.5, range.getUpperBound()); // force zero bin
			num = (int)Math.round((range.getUpperBound()-range.getLowerBound())/delta) + 1;
			if (num == 1)
				delta /= 2d;
		}
		System.out.println("Range: "+range);
		System.out.println("Num: "+num);
		System.out.println("Delta: "+delta);
		Map<Double, HistogramFunction> hists = Maps.newHashMap();
		for (double duration : catalogDists.keySet()) {
			hists.put(duration,getLossHist(catalogDists.get(duration),
					range.getLowerBound(), num, delta, isLog10));
		}
		return hists;
	}
	
	/**
	 * Calculates a nicely rounded range that will exactly cover the data range with the given discretization.
	 * Histogram num points can be calculated as num = (range.getUpperBound()-range.getLowerBound())/delta + 1;
	 * @param delta
	 * @return
	 */
	public static Range calcSmartHistRange(List<DiscretizedFunc> catalogDists, double delta, boolean isLog10) {
		MinMaxAveTracker xTrack = new MinMaxAveTracker();
		for (DiscretizedFunc func : catalogDists) {
			xTrack.addValue(func.getMinX());
			xTrack.addValue(func.getMaxX());
		}
		Range range = new Range(xTrack.getMin(), xTrack.getMax());
		if (isLog10)
			range = new Range(Math.log10(range.getLowerBound()), Math.log10(range.getUpperBound()));
		double min = Math.floor(range.getLowerBound()/delta) * delta;
		min += 0.5*delta;
		// it's possible that this was too conservative and that the bin above could hold the range min val
		// due to the bin width
		if (min + 0.5*delta < range.getLowerBound())
			min += delta;
		double max = min;
		while (max+0.5*delta < range.getUpperBound())
			max += delta;
		
		return new Range(min, max);
	}
	
	public static HistogramFunction getLossHist(List<DiscretizedFunc> catalogDists,
			double minX, int numX, double deltaX, boolean isLog10) {
		HistogramFunction lossHist = new HistogramFunction(minX, numX, deltaX);
		
		double distMin = minX - 0.5*deltaX;
		double distMax = lossHist.getMaxX() + 0.5*deltaX;
		
		double weightMult = 1d/catalogDists.size();
		
		double[] lossVals = new double[catalogDists.size()];
		
		for (int i = 0; i < catalogDists.size(); i++) {
			DiscretizedFunc catalogDist = catalogDists.get(i);
			double meanVal = 0d;
			for (int index=0; index<catalogDist.size(); index++) {
				double loss = catalogDist.getX(index);
				double weight = catalogDist.getY(index);
				
				meanVal += loss*weight;
				
				weight *= weightMult;
				
				if (loss == 0) {
					lossHist.add(0, weight);
					continue;
				}
				
				if (isLog10)
					loss = Math.log10(loss);
				
				if (loss > distMax)
					// above max, put in last bin
					lossHist.add(lossHist.size()-1, weight);
				else if (loss < distMin || lossHist.size() == 1)
					// below min, put in first bin
					lossHist.add(0, weight);
				else
					// put in actual bin
					lossHist.add(loss, weight);
			}
			lossVals[i] = meanVal;
		}
		
		double sumY = lossHist.calcSumOfY_Vals();
//		System.out.println("Sum Y values: "+sumY);
		Preconditions.checkState(DataUtils.getPercentDiff(sumY, 1d) < 1, "loss hist sum of y vals doesn't equal 1: "+(float)sumY);
		lossHist.scale(1d/sumY);
		
		double mean = StatUtils.mean(lossVals);
		double median = DataUtils.median(lossVals);
		double mode = lossHist.getMode();
		System.out.println("Mean loss: "+mean);
		System.out.println("Median loss: "+median);
		String info = "Mean: "+mean+", Median: "+median+", Mode: "+mode;
		lossHist.setInfo(info);
		
		return lossHist;
	}
	
	static String getDurationLabel(double duration) {
		if (duration >= 1d) {
			if (duration == Math.floor(duration))
				return (int)duration+" yr";
			return (float)duration+" yr";
		}
		// check for months
		double oneMonth = 30d/365.25;
		double months = duration/oneMonth;
		if (months == Math.floor(months))
			return (int)(months)+" mo";
		double oneWeek = 7d/365.25;
		if (duration == oneWeek)
			return "1 wk";
		double days = duration*365.25;
		if (days == Math.floor(days))
			return (int)(days)+" day";
		return (float)(days)+" day";
	}
	
	private static final DecimalFormat labelDF = new DecimalFormat("0.000");
	
	public static void writeLossHist(File outputDir, String outputPrefix, Map<Double, HistogramFunction> lossHists, boolean isLogX,
			boolean triggeredOnly, String xAxisLabel, double maxX) throws IOException {
		lossHists = Maps.newHashMap(lossHists);
		if (!outputDir.exists())
			outputDir.mkdir();
		
//		double logMinY = (1d/catalogs.size())/10d; // one order of magnitude below the smallest possible
		double logMinY = 1e-5;
		
		List<PlotElement> elems = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		List<Double> durations = Lists.newArrayList();
		for (Double duration : lossHists.keySet())
			durations.add(duration);
		Collections.sort(durations);
		Collections.reverse(durations); // order from biggest to smallest
		double minDur = durations.get(durations.size()-1);
		double maxDur = durations.get(0);
		
		CPT logDurCPT = null;
		if (lossHists.size() > 1)
			logDurCPT = new CPT(Math.log(minDur), Math.log(maxDur), Color.GRAY, Color.BLACK);
		
		double histsMinX = Double.POSITIVE_INFINITY;
		double histsMaxX = Double.NEGATIVE_INFINITY;
		double histsMinY = Double.POSITIVE_INFINITY;
		double histsMaxY = Double.NEGATIVE_INFINITY;
		
		for (double duration : durations) {
			HistogramFunction lossHist = lossHists.get(duration);
			
//			if (xAxisScale != 1d && xAxisScale != 0d) {
//				HistogramFunction newLossHist = new HistogramFunction(
//						lossHist.getMinX()*xAxisScale, lossHist.size(), lossHist.getDelta()*xAxisScale);
//				System.out.println("Scaled minX: "+newLossHist.getMinX()+" maxX="+newLossHist.getMaxX());
//				for (int i=0; i<lossHist.size(); i++)
//					newLossHist.set(i, lossHist.getY(i));
////				newLossHist.setInfo(info); // TODO
//				lossHist = newLossHist;
//				lossHists.put(duration, lossHist); // we made a copy of the map so this doesn't modify the original
//			}
			
			double mean = lossHist.computeMean();
			
			lossHist.setName(getDurationLabel(duration)+" (mean="+labelDF.format(mean)+")");
			
			elems.add(lossHist);
			Color c;
			if (logDurCPT == null)
				c = Color.BLACK;
			else
				c = logDurCPT.getColor((float)Math.log(duration));
			if (duration == maxDur)
				chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, c));
			else
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, c));
			
			histsMinX = Math.min(histsMinX, lossHist.getMinX()-0.5*lossHist.getDelta());
			histsMaxX = Math.max(histsMaxX, lossHist.getMaxX()+0.5*lossHist.getDelta());
			histsMinY = Math.min(histsMinY, lossHist.getMinY());
			histsMaxY = Math.max(histsMinY, lossHist.getMaxY());
		}
		
		if (isLogX)
			xAxisLabel = "Log10("+xAxisLabel+")";
		
		PlotSpec spec = new PlotSpec(elems, chars, "Simulated Loss Distribution", xAxisLabel, "Frequency");
		if (durations.size() > 1)
			spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		
		gp.setBackgroundColor(Color.WHITE);
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		
//		double minX = lossHist.getMinX()-0.5*lossHist.getDelta();
//		if (!Doubles.isFinite(maxX) || maxX <= 0)
//			maxX = lossHist.getMaxX()+0.5*lossHist.getDelta();
//		double maxY = lossHist.getMaxY()*1.2;
		double minX = histsMinX;
		if (!Doubles.isFinite(maxX) || maxX <= 0)
			maxX = histsMaxX;
		double maxY = histsMaxY*1.2;
		
		for (boolean logY : new boolean[] {false, true}) {
			if (logY)
				gp.setUserBounds(minX, maxX, logMinY, maxY);
			else
				gp.setUserBounds(minX, maxX, 0, maxY);
			gp.drawGraphPanel(spec, false, logY);
			gp.getChartPanel().setSize(1000, 800);
			String myPrefix = outputPrefix;
			if (logY)
				myPrefix += "_logY";
			if (triggeredOnly)
				myPrefix += "_triggered";
			gp.saveAsPNG(new File(outputDir, myPrefix+".png").getAbsolutePath());
			gp.saveAsPDF(new File(outputDir, myPrefix+".pdf").getAbsolutePath());
			gp.saveAsTXT(new File(outputDir, myPrefix+".txt").getAbsolutePath());
		}
	}
	
	static Map<Double, DiscretizedFunc> toExceedFuncs(Map<Double, HistogramFunction> lossHists,
			int nForConf, boolean allSubDurations, double catalogDuration) {
		Map<Double, DiscretizedFunc> exceeds = Maps.newHashMap();
		for (Double duration : lossHists.keySet()) {
			int myNForConf = nForConf;
			if (nForConf > 0 && allSubDurations) {
				Preconditions.checkState(duration > 0);
				// sub durations
				int numPer = (int)(catalogDuration/duration);
				myNForConf = numPer * nForConf;
				if (myNForConf > max_all_sub_durations_n)
					myNForConf = max_all_sub_durations_n;
				System.out.println("All sub durations N for duration "+duration+": "+myNForConf);
			}
			exceeds.put(duration, toExceedFunc(lossHists.get(duration), myNForConf));
		}
		return exceeds;
	}
	
	static DiscretizedFunc toExceedFunc(HistogramFunction lossHist, int nForConf) {
		HistogramFunction cumDist = lossHist.getCumulativeDistFunctionWithHalfBinOffset();
		// convert to exceedance
		for (int i=0; i<cumDist.size(); i++) {
			double y = 1d-cumDist.getY(i);
			if ((float)y == 0f)
				y = 0;
			if (y < 0 || y < 1e-12) {
				Preconditions.checkState(y > -1e-10);
				y = 0;
			}
			cumDist.set(i, y);
		}
		
		if (nForConf > 0) {
			// generate confidence intervals
			HistogramFunction lower = new HistogramFunction(cumDist.getMinX(), cumDist.getMaxX(), cumDist.size());
			HistogramFunction upper = new HistogramFunction(cumDist.getMinX(), cumDist.getMaxX(), cumDist.size());
			for (int i=0; i<cumDist.size(); i++) {
				double p = cumDist.getY(i);
				double[] conf;
				if (conf_68)
					conf = ETAS_Utils.getBinomialProportion68confidenceInterval(p, nForConf);
				else
					conf = ETAS_Utils.getBinomialProportion95confidenceInterval(p, nForConf);
				Preconditions.checkState((float)conf[0] <= (float)p,
						"Bad conf. p=%s, n=%s, [%s %s]", p, nForConf, conf[0], conf[1]);
				Preconditions.checkState((float)conf[1] >= (float)p,
						"Bad conf. p=%s, n=%s, [%s %s]", p, nForConf, conf[0], conf[1]);
				lower.set(i, conf[0]);
				upper.set(i, conf[1]);
			}
			return new UncertainArbDiscDataset(cumDist, lower, upper);
		}
		return cumDist;
	}
	
	private static HistogramFunction getSameNumX(DiscretizedFunc func,
			double minX, double maxX, int num, boolean aboveUseLast) {
		HistogramFunction ret = new HistogramFunction(minX, maxX, num);
		
		for (int i=0; i<num; i++) {
			double x = ret.getX(i);
			double y;
			if (x < func.getMinX()) {
				// if before, use first bin
				y = func.getY(0);
			} else if (x > func.getMaxX()) {
				// if after, set to zero
				if (aboveUseLast)
					// for upper confidence bound
					y = func.getY(func.size()-1);
				else
					y = 0d;
			} else {
				y = func.getY(x);
			}
			ret.set(i, y);
		}
		return ret;
	}
	
	public static void writeLossExceed(File outputDir, String outputPrefix, Map<Double, DiscretizedFunc> exceedFuncs,
			boolean isLogX, boolean triggeredOnly, String xAxisLabel, double maxX) throws IOException {
		Table<String, Double, DiscretizedFunc> table = HashBasedTable.create();
		for (double duration : exceedFuncs.keySet())
			table.put("", duration, exceedFuncs.get(duration));
		writeLossExceed(outputDir, outputPrefix, table, null, isLogX, triggeredOnly, xAxisLabel, maxX, false, true);
	}
	
	public static void writeLossExceed(File outputDir, String outputPrefix,
			Table<String, Double, DiscretizedFunc> exceedFuncs, Map<String, Double> exceedWeightMap,
			boolean isLogX, boolean triggeredOnly, String xAxisLabel, double maxX,
			boolean individualDurations, boolean averageConf) throws IOException {
		if (!outputPrefix.endsWith("exceed"))
			outputPrefix += "_exceed";
		
		List<PlotElement> elems = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();

		double logMinY = 1e-5;

		List<Double> durations = Lists.newArrayList();
		for (Double duration : exceedFuncs.columnKeySet())
			durations.add(duration);
		Collections.sort(durations);
		Collections.reverse(durations); // order from biggest to smallest
		double minDur = durations.get(durations.size()-1);
		double maxDur = durations.get(0);

		List<String> names = Lists.newArrayList();
		for (String row : exceedFuncs.rowKeySet())
			names.add(row);
		if (names.size() > 1)
			// if size == 1, can be null so don't sort
			Collections.sort(names);
		
		boolean hasConf = exceedFuncs.values().iterator().next() instanceof UncertainArbDiscDataset;

		CSVFile<String> csv = new CSVFile<String>(false);
		List<String> header = Lists.newArrayList(xAxisLabel);
		for (double duration : durations) {
			header.add(getDurationLabel(duration));
			if (hasConf) {
				header.add(conf_lower_str);
				header.add(conf_upper_str);
			}
			if (names.size() > 1) {
				for (String name : names) {
					header.add(name);
					if (hasConf) {
						header.add(conf_lower_str);
						header.add(conf_upper_str);
					}
				}
			}
		}
		csv.addLine(header);

		CPT logDurCPT = null;
		if (exceedFuncs.size() > 1)
			logDurCPT = new CPT(Math.log(minDur), Math.log(maxDur), Color.GRAY, Color.BLACK);

		double maxCumY = 0d;
		
		List<PlotSpec> individualSpecs = null;
		if (individualDurations)
			individualSpecs = Lists.newArrayList();

		for (double duration : durations) {
			
			Map<String, DiscretizedFunc> cumDists = exceedFuncs.column(duration);

			DiscretizedFunc cumDist; // mean
			if (cumDists.size() > 1) {
				double distMinX = Double.POSITIVE_INFINITY;
				double distMaxX = 0;
				double distDelta = Double.NaN;

				for (DiscretizedFunc dist : cumDists.values()) {
					double myDelta = dist.getX(1) - dist.getX(0);
					if (Double.isNaN(distDelta))
						distDelta = myDelta;
					else
						Preconditions.checkState((float)distDelta == (float)myDelta);
					distMinX = Math.min(distMinX, dist.getMinX());
					distMaxX = Math.max(distMaxX, dist.getMaxX());
				}
				int num = (int)((distMaxX - distMinX)/distDelta + 0.5) + 1;
				HistogramFunction mean = new HistogramFunction(distMinX, distMaxX, num);
				HistogramFunction lower = null, upper = null;
				if (hasConf) {
					lower = new HistogramFunction(distMinX, distMaxX, num);
					upper = new HistogramFunction(distMinX, distMaxX, num);
					if (!averageConf) {
						// fill in with default vals
						for (int i=0; i<num; i++) {
							lower.set(i, Double.POSITIVE_INFINITY);
							upper.set(i, 0d);
						}
					}
				}
//				cumDist = new HistogramFunction(distMinX, distMaxX, num);

				double rateEven = 1d/(double)cumDists.size();

				for (String name : names) {
					DiscretizedFunc dist = cumDists.get(name);
					// need the same x axis points (one of the sub ones can be smaller, but have same spacing)
//					HistogramFunction newDist = new HistogramFunction(distMinX, distMaxX, num);
					DiscretizedFunc newDist = getSameNumX(dist, distMinX, distMaxX, num, false);
					HistogramFunction myLower = null, myUpper = null;
					if (hasConf) {
						UncertainArbDiscDataset distWithConf = (UncertainArbDiscDataset)dist;
						myLower = getSameNumX(distWithConf.getLower(), distMinX, distMaxX, num, false);
						myUpper = getSameNumX(distWithConf.getUpper(), distMinX, distMaxX, num, true);
					}
					
					double myRate;
					if (exceedWeightMap == null)
						myRate = rateEven;
					else
						myRate = exceedWeightMap.get(name);
					
					for (int i=0; i<newDist.size(); i++) {
						double x = newDist.getX(i);
						double y = newDist.getY(i);
						mean.add(x, y*myRate);
						
						if (hasConf) {
							if (averageConf) {
								lower.add(x, myLower.getY(i)*myRate);
								upper.add(x, myUpper.getY(i)*myRate);
							} else {
								// take min/max
								lower.set(x, Math.min(lower.getY(x), myLower.getY(i)));
								upper.set(x, Math.max(upper.getY(x), myUpper.getY(i)));
							}
						}
					}
					if (hasConf)
						newDist = new UncertainArbDiscDataset(newDist, myLower, myUpper);
					cumDists.put(name, newDist);
				}
				if (hasConf)
					cumDist = new UncertainArbDiscDataset(mean, lower, upper);
				else
					cumDist = mean;
			} else {
				cumDist = cumDists.values().iterator().next();
			}

			if (csv.getNumRows() == 1) {
				// first time through, create empty lines
				for (int i=0; i<cumDist.size(); i++)
					csv.addLine(Lists.newArrayList(cumDist.getX(i)+""));
			} else {
				Preconditions.checkState(cumDist.size() == csv.getNumRows()-1);
			}
			for (int i=0; i<cumDist.size(); i++) {
				List<String> line = csv.getLine(i+1);
				line.add(cumDist.getY(i)+"");
				if (hasConf) {
					line.add(((UncertainArbDiscDataset)cumDist).getLowerY(i)+"");
					line.add(((UncertainArbDiscDataset)cumDist).getUpperY(i)+"");
				}
				if (names.size() > 1) {
					for (String name : names) {
						DiscretizedFunc subDist = cumDists.get(name);
						line.add(subDist.getY(i)+"");
						if (hasConf) {
							line.add(((UncertainArbDiscDataset)subDist).getLowerY(i)+"");
							line.add(((UncertainArbDiscDataset)subDist).getUpperY(i)+"");
						}
					}
				}
			}

			cumDist.setName(getDurationLabel(duration));

			Color c;
			if (logDurCPT == null)
				c = Color.BLACK;
			else
				c = logDurCPT.getColor((float)Math.log(duration));
			
			List<PlotElement> myElems = Lists.newArrayList();
			List<PlotCurveCharacterstics> myChars = Lists.newArrayList();

			Color rangeColor = new Color((c.getRed()+255)/2, (c.getGreen()+255)/2, (c.getBlue()+255)/2);
			if (hasConf) {
//				UncertainArbDiscDataset confRange = null;
//				try {
//					confRange = new UncertainArbDiscDataset(cumDist, conf[0], conf[1]);
//				} catch (RuntimeException e) {
//					System.out.println("x\tlow\tmean\tup");
//					for (int i=0; i<cumDist.size(); i++) {
//						System.out.println((float)cumDist.getX(i)+"\t"+(float)conf[0].getY(i)
//								+"\t"+(float)cumDist.getY(i)+"\t"+(float)conf[1].getY(i));
//					}
//					throw e;
//				}

				myElems.add(cumDist);
				myChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN_TRANS, 1f, rangeColor));
			}

			if (cumDists.size() > 1) {
				HistogramFunction lowerFunc = new HistogramFunction(cumDist.getMinX(), cumDist.getMaxX(), cumDist.size());
				HistogramFunction upperFunc = new HistogramFunction(cumDist.getMinX(), cumDist.getMaxX(), cumDist.size());
				for (int i=0; i<cumDist.size(); i++) {
					double min = Double.POSITIVE_INFINITY;
					double max = 0d;
					for (DiscretizedFunc dist : cumDists.values()) {
						double y = dist.getY(i);
						min = Math.min(min, y);
						max = Math.max(max, y);
					}
					lowerFunc.set(i, min);
					upperFunc.set(i, max);
				}
				UncertainArbDiscDataset meanRange = new UncertainArbDiscDataset(cumDist, lowerFunc, upperFunc);

				myElems.add(meanRange);
				myChars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, rangeColor));
			}

			myElems.add(cumDist);
			myChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, c));

			//			// now confidence intervals
			//			if (conf != null) {
			//				elems.add(conf[0]);
			//				chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, c));
			//				elems.add(conf[1]);
			//				chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, c));
			//			}

			maxCumY = Math.max(maxCumY, cumDist.getMaxY());
			
			elems.addAll(myElems);
			chars.addAll(myChars);
			
			if (individualDurations) {
				// write one for just this duration
				PlotSpec spec = new PlotSpec(myElems, myChars, "Simulated Loss Exceedance Probabilities", xAxisLabel, "Exceedance Prob");
				individualSpecs.add(spec);
			}
		}

		PlotSpec combinedSpec = new PlotSpec(elems, chars, "Simulated Loss Exceedance Probabilities", xAxisLabel, "Exceedance Prob");
		if (durations.size() > 1)
			combinedSpec.setLegendVisible(true);

		double maxY = maxCumY*1.2;
		if (maxY > 1d)
			maxY = 1.05;

		String csvName = outputPrefix;
		if (triggeredOnly)
			csvName += "_triggered";
		csv.writeToFile(new File(outputDir, csvName+".csv"));

		HeadlessGraphPanel gp = new HeadlessGraphPanel();

		gp.setBackgroundColor(Color.WHITE);
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);

		double minX = 0d;
		
		List<PlotSpec> allSpecs = Lists.newArrayList();
		List<String> allPrefixes = Lists.newArrayList();
		
		allSpecs.add(combinedSpec);
		allPrefixes.add(outputPrefix);
		
		if (individualDurations) {
			for (int i=0; i<durations.size(); i++) {
				double duration = durations.get(i);
				allSpecs.add(individualSpecs.get(i));
				allPrefixes.add(outputPrefix+"_"+getDurationLabel(duration).replaceAll(" ", ""));
			}
		}

		for (int i=0; i<allSpecs.size(); i++) {
			PlotSpec spec = allSpecs.get(i);
			outputPrefix = allPrefixes.get(i);
			for (boolean logY : new boolean[] {false, true}) {
				if (logY)
					gp.setUserBounds(minX, maxX, logMinY, maxY);
				else
					gp.setUserBounds(minX, maxX, 0, maxY);
				gp.drawGraphPanel(spec, false, logY);
				gp.getChartPanel().setSize(1000, 800);
				String myPrefix = outputPrefix;
				if (logY)
					myPrefix += "_logY";
				if (triggeredOnly)
					myPrefix += "_triggered";
				gp.saveAsPNG(new File(outputDir, myPrefix+".png").getAbsolutePath());
				gp.saveAsPDF(new File(outputDir, myPrefix+".pdf").getAbsolutePath());
				gp.saveAsTXT(new File(outputDir, myPrefix+".txt").getAbsolutePath());
			}
		}
	}
	
	public void writeLossesToCSV(File dir, String prefix, Map<Double, List<DiscretizedFunc>> lossDists) throws IOException {
//		for (double duration : lossDists.keySet()) {
//			File csvFile = new File(dir, prefix+"_"+getDurationLabel(duration).replaceAll(" ", "")+".csv");
//			writeLossesToCSV(csvFile, lossDists.get(duration));
//		}
		List<Double> durations = Lists.newArrayList(lossDists.keySet());
		Collections.sort(durations);
		
		// write maximum
		double maxDur = durations.get(durations.size()-1);
		List<DiscretizedFunc> maxLosses = lossDists.get(maxDur);
		if (maxLosses.size() == catalogs.size()) {
			// won't if all sub durations
			File csvFile = new File(dir, prefix+"_"+getDurationLabel(maxDur).replaceAll(" ", "")+".csv");
			writeLossesToCSV(csvFile, maxLosses);
		}
		
		// can differ if all sub durations
		int maxCatalogs = 0;
		for (double duration : lossDists.keySet()) {
			int myNum = lossDists.get(duration).size();
			if (myNum > maxCatalogs)
				maxCatalogs = myNum;
		}
		
		// now write combined csv
		CSVFile<String> csv = new CSVFile<String>(true);
		List<String> header = Lists.newArrayList("Index");
		for (double duration : durations)
			header.add(duration+"");
		csv.addLine(header);
		for (int i=0; i<maxCatalogs; i++) {
			List<String> line = Lists.newArrayList(i+"");
			for (double duration : durations) {
				List<DiscretizedFunc> myLosses = lossDists.get(duration);
				if (i >= myLosses.size()) {
					line.add("");
					continue;
				}
				double totLoss = 0d;
				for (Point2D pt : myLosses.get(i))
					totLoss += pt.getX()*pt.getY();
				line.add(totLoss+"");
			}
			csv.addLine(line);
		}
		csv.writeToFile(new File(dir, prefix+"_combined.csv"));
	}
	
	public void writeLossesToCSV(File csvFile, List<DiscretizedFunc> lossDists) throws IOException {
		Preconditions.checkState(lossDists.size() == catalogs.size(), "Have %s dists but %s catalogs!",
				lossDists.size(), catalogs.size());
		
		CSVFile<String> csv = new CSVFile<String>(true);
		
		double cutoffMag = AbstractGridSourceProvider.SOURCE_MIN_MAG_CUTOFF;
		csv.addLine("Index", "Total Mean Loss", "# FSS Ruptures",
				"# M>="+(float)cutoffMag, "Max Mag");
		
		for (int i=0; i<catalogs.size(); i++) {
			double totLoss = 0d;
			for (Point2D pt : lossDists.get(i))
				totLoss += pt.getX()*pt.getY();
			int numFSSRups = 0;
			int numAbove = 0;
			double maxMag = 0;
			for (ETAS_EqkRupture rup : catalogs.get(i)) {
				if (getFSSIndex(rup) >= 0)
					numFSSRups++;
				if (rup.getMag() > maxMag)
					maxMag = rup.getMag();
				if ((float)rup.getMag() >= (float)cutoffMag)
					numAbove++;
			}
			csv.addLine(i+"", totLoss+"", numFSSRups+"", numAbove+"", maxMag+"");
		}
		
		csv.writeToFile(csvFile);
	}
	
	private void scenarioSearch() {
		// this method can be used to search for certain scenarios
		
		// fault based in LA box but not SAF
		Region region = new CaliforniaRegions.LA_BOX();
		HashSet<Integer> sectIDs = new HashSet<Integer>();
		HashSet<Integer> safParents = new HashSet<Integer>(
				FaultModels.FM3_1.getNamedFaultsMapAlt().get("San Andreas"));
		
		FaultSystemRupSet rupSet = meanSol.getRupSet();
		
		long ot = Math.round((2014.0-1970.0)*ProbabilityModelsCalc.MILLISEC_PER_YEAR); // occurs at 2014
		
		sectLoop:
		for (FaultSectionPrefData sect : rupSet.getFaultSectionDataList()) {
			if (safParents.contains(sect.getParentSectionId()))
				continue;
			for (Location loc : sect.getFaultTrace()) {
				if (region.contains(loc)) {
					sectIDs.add(sect.getSectionId());
					continue sectLoop;
				}
			}
		}
		
		for (int i=0; i<catalogs.size(); i++) {
			List<ETAS_EqkRupture> catalog = catalogs.get(i);
			rupLoop:
			for (ETAS_EqkRupture rup : catalog) {
				int fssIndex = getFSSIndex(rup);
				if (fssIndex >= 0) {
					List<FaultSectionPrefData> data = rupSet.getFaultSectionDataForRupture(fssIndex);
					for (int sectID : rupSet.getSectionsIndicesForRup(fssIndex)) {
						if (sectIDs.contains(sectID)) {
							String name = data.size()+" SECTIONS BETWEEN "+data.get(0).getName()
									+" AND "+data.get(data.size()-1).getName();
							float mag = (float)rup.getMag();
							float deltaDays = (float)((rup.getOriginTime()-ot)/1000d/60d/60d/24d);
							System.out.println("catalog "+i+" has a M"+mag+" match "+deltaDays+" days after on: "+name);
							continue rupLoop;
						}
					}
				}
			}
		}
	}
	
	static final double[] durations = { 1d/365.25, 7d/365.25, 30/365.25, 1d, 10d, 30d, 50d, 100d };

	public static void main(String[] args) throws IOException, DocumentException {
//		File parentDir = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/2014_05_28-la_habra/");
//		File parentDir = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/2014_05_28-mojave_7/");
//		File parentDir = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/2014_05_28-spontaneous/");
//		File parentDir = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/2014_06_26-mojave_7/");
		
		// get result dirs
//		List<File> resultDirs = Lists.newArrayList();
//		File[] subDirs = parentDir.listFiles();
//		Arrays.sort(subDirs, new FileNameComparator());
//		for (File subDir : subDirs) {
//			if (subDir.isDirectory() && subDir.getName().startsWith("results"))
//				resultDirs.add(subDir);
//		}
//		File[] etasCatalogsDirs = resultDirs.toArray(new File[0]);
		
//		File zipFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/2014_09_02-mojave_7/results.zip");
//		File zipFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/2014_09_02-napa/results.zip");
//		File zipFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/2014_09_02-spontaneous/results.zip");
		File resultsFile;
		if (args.length >= 1)
			resultsFile = new File(args[0]);
		else
			resultsFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
//				+ "2016_02_19-mojave_m7-10yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14/results_descendents.bin");
//				+ "2016_02_19-mojave_m7-10yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14-combined100k/results_m5_preserve.bin");
//				+ "2016_02_25-surprise_valley_5p0-10yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14/results_descendents.bin");
//				+ "2016_08_24-spontaneous-10yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14-combined/results_m5.bin");
//				+ "2016_02_17-spontaneous-1000yr-scaleMFD1p14-full_td-subSeisSupraNucl-gridSeisCorr/results_m4.bin");
//				+ "2016_08_31-bombay_beach_m4pt8-10yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14-combined/results_m5_preserve.bin");
				+ "2016_02_17-spontaneous-1000yr-scaleMFD1p14-full_td-subSeisSupraNucl-gridSeisCorr/results_m5_preserve.bin");
		
		boolean triggeredOnly = false;
		if (args.length > 1)
			triggeredOnly = Boolean.parseBoolean(args[1]);
		boolean allSubDurations = resultsFile.getParentFile().getName().contains("1000yr");
		if (args.length > 2)
			allSubDurations = Boolean.parseBoolean(args[2]);
		
		// true mean FSS which includes rupture mapping information. this must be the exact file used to calculate EALs
		File trueMeanSolFile = new File("dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_"
				+ "COMPOUND_SOL_TRUE_HAZARD_MEAN_SOL_WITH_MAPPING.zip");

		// directory which contains EAL data
		File ealMainDir = new File("/home/kevin/OpenSHA/UCERF3/eal");
		if (!ealMainDir.exists())
			ealMainDir = new File("/var/www/html/ftp/kmilner/ucerf3/eal_calcs");
		if (!ealMainDir.exists())
			ealMainDir = new File("/home/scec-02/kmilner/ucerf3/eal");
		Preconditions.checkState(ealMainDir.exists(), "Directory doesn't exist: %s", ealMainDir.getAbsolutePath());
		
		// constants for all catalogs
//		double xAxisScale = 1d/1e6; // portfolio units are in thousands (1e3), so convert to billions by dividing by 1e6
		String xAxisLabel = "$ (Billions)";
		double maxX = 200;
//		double deltaX = 1e6;
		double deltaX = 1;
		double thousandsToBillions = 1d/1e6; // portfolio units are in thousands (1e3), so convert to billions by dividing by 1e6
		
		double inflationScalar = 1d/0.9d;
		
		double xAxisScale = thousandsToBillions*inflationScalar;
		
		List<File> dataDirs = Lists.newArrayList();
		
		if (args.length > 3) {
			for (int i=3; i<args.length; i++)
				dataDirs.add(new File(ealMainDir, args[i]));
		} else {
//			dataDirs.add(new File(ealMainDir, "2014_05_28-ucerf3-99percent-wills-smaller"));
//			dataDirs.add(new File(ealMainDir, "2016_06_06-ucerf3-90percent-wald"));
			dataDirs.add(new File(ealMainDir, "2016_10_18-ucerf3-90percent-wald-san-bernardino"));
			dataDirs.add(new File(ealMainDir, "2016_11_28-ucerf3-90percent-wills-san-bernardino"));
//			dataDirs.add(new File(ealMainDir, "2016_10_21-ucerf3-90percent-wald-coachella-valley"));
		}
		
		if (dataDirs.get(0).getName().contains("san-bernardino") || dataDirs.get(0).getName().contains("coachella"))
			maxX = 50;
		
		// CEA proxy wald
//		File dataDir = new File(ealMainDir, "2014_05_05-ucerf3-eal-calc-wald-vs30");
//		String xAxisLabel = "$ (Billions)";
//		double xAxisScale = 1d/1e6; // portfolio units are in thousands (1e3), so convert to billions by dividing by 1e6
//		double maxX = 120;
//		double deltaX = 1e6;
//		String catOutputDirName = "cea_proxy_wald";
		
//		// 99% Wills
//		File dataDir = new File(ealMainDir, "2014_05_28-ucerf3-99percent-wills-smaller");
//		
//		
//		String catOutputDirName = "ca_99_wills";
		
		// Fatality portfolio
//		File dataDir = new File(ealMainDir, "2014_05_28-ucerf3-fatality-smaller");
//		String xAxisLabel = "Fatalities";
//		double xAxisScale = 1d;
//		double maxX = 3000;
//		double deltaX = 1;
//		String catOutputDirName = "fatalities_wills";

		// IMR for which EAL data has already been computed
		Map<AttenRelRef, Double> imrWeightsMap = Maps.newHashMap();
		imrWeightsMap.put(AttenRelRef.CB_2014, 0.22);
		imrWeightsMap.put(AttenRelRef.CY_2014, 0.22);
		imrWeightsMap.put(AttenRelRef.ASK_2014, 0.22);
		imrWeightsMap.put(AttenRelRef.BSSA_2014, 0.22);
		imrWeightsMap.put(AttenRelRef.IDRISS_2014, 0.12);

		// Fault model of interest
		FaultModels fm = FaultModels.FM3_1;

		// Branch averaged FSS
		FaultSystemSolution baSol = FaultSystemIO.loadSol(
				new File("dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_"
				+ "COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));

		// Compound fault system solution
		CompoundFaultSystemSolution cfss = CompoundFaultSystemSolution.fromZipFile(
				new File("dev/scratch/UCERF3/data/scratch/"
						+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip"));
		
		FaultSystemSolution trueMeanSol = FaultSystemIO.loadSol(trueMeanSolFile);
		Map<LogicTreeBranch, List<Integer>> branchMappings = TrueMeanBuilder.loadRuptureMappings(trueMeanSolFile);
		
		System.out.println("Triggered only? "+triggeredOnly);
		System.out.println("All sub durations? "+allSubDurations);
		System.out.println("Data dirs:");
		for (File dataDir : dataDirs)
			System.out.println("\t"+dataDir.getAbsolutePath());
		
		calculate(resultsFile, triggeredOnly, xAxisLabel, maxX, deltaX, xAxisScale, dataDirs, imrWeightsMap, fm, baSol,
				cfss, trueMeanSol, branchMappings, durations, allSubDurations);
	}

	public static void calculate(File resultsFile, boolean triggeredOnly, String xAxisLabel, double maxX,
			double deltaX, double xAxisScale, List<File> dataDirs, Map<AttenRelRef, Double> imrWeightsMap,
			FaultModels fm, FaultSystemSolution baSol, CompoundFaultSystemSolution cfss,
			FaultSystemSolution trueMeanSol, Map<LogicTreeBranch, List<Integer>> branchMappings,
			double[] durations, boolean allSubDurations) throws IOException, DocumentException {
		File lossOutputDir = new File(resultsFile.getParentFile(), "loss_results");
		Preconditions.checkState(lossOutputDir.exists() || lossOutputDir.mkdir());
		
		List<Map<Double, List<DiscretizedFunc>>> lossDistsList = Lists.newArrayList();
		List<Double> lossWeights = Lists.newArrayList();
		Table<String, Double, DiscretizedFunc> allCombLossExceeds = HashBasedTable.create();
		Map<String, Double> allExceedWeightMap = Maps.newHashMap();
		
		TestScenario scenario = ETAS_MultiSimAnalysisTools.detectScenario(resultsFile.getParentFile());
		if (scenario != null && scenario.getFSS_Index() >= 0)
			scenario.updateMag(baSol.getRupSet().getMagForRup(scenario.getFSS_Index()));
		
		if (resultsFile.getParentFile().getName().startsWith("2016_02_19-mojave")) {
			System.out.println("Changing scenario ID");
			id_for_scenario = 9893;
		} else {
			id_for_scenario = 0;
		}
		
		String prefixAdd = "";
		if (allSubDurations)
			prefixAdd += "_all_sub_durations";
		
		String csvPrefixAdd = "";
		if (triggeredOnly)
			csvPrefixAdd = "_triggered";
		
		List<List<ETAS_EqkRupture>> catalogs = null;
		
		ETAS_CatalogEALCalculator calc = null;
		
		boolean isLog10 = false; // x axis
		
		for (int i=0; i<dataDirs.size(); i++) {
			File dataDir = dataDirs.get(i);
			System.out.println("Handling data dir: "+dataDir.getAbsolutePath());
			
			UCERF3_BranchAvgLossFetcher fetcher = new UCERF3_BranchAvgLossFetcher(trueMeanSol, branchMappings, cfss, dataDir);
			
			if (catalogs == null) {
				calc = new ETAS_CatalogEALCalculator(fetcher, baSol, fm, resultsFile);
				catalogs = calc.catalogs;
				// trim durations
				List<Double> myDurations = Lists.newArrayList();
				double maxDuration = 0d;
				for (List<ETAS_EqkRupture> catalog : catalogs)
					if (!catalog.isEmpty())
						maxDuration = Math.max(maxDuration, ETAS_MultiSimAnalysisTools.calcDurationYears(catalog));
				Preconditions.checkState(maxDuration > 0);
				System.out.println("Max catalog direction detected: "+maxDuration);
				// pad max duration by 50%
				maxDuration *= 1.5;
				System.out.println("All durations <= "+maxDuration+" will be considered");
				for (double duration : durations)
					if (duration <= maxDuration)
						myDurations.add(duration);
				Collections.sort(myDurations);
				durations = Doubles.toArray(myDurations);
			} else {
				calc = new ETAS_CatalogEALCalculator(fetcher, baSol, fm, catalogs);
			}
			calc.setTriggeredOnly(triggeredOnly);
			if (scenario != null) {
				if (scenario.getFSS_Index() >= 0)
					calc.setTriggerFaultRup(scenario.getFSS_Index());
				else
					calc.setTriggerGridRup(scenario.getLocation(), scenario.getMagnitude());
			}
			
			List<Map<Double, List<DiscretizedFunc>>> myLossDists = Lists.newArrayList();
			List<Double> myWeights = Lists.newArrayList();
			
			File outputDir = new File(lossOutputDir, dataDir.getName());
			if (!outputDir.exists())
				outputDir.mkdir();
			
			Table<String, Double, DiscretizedFunc> gmpeCombLossExceeds = HashBasedTable.create();
			Map<String, Double> exceedWeightMap = Maps.newHashMap();
			for (AttenRelRef attenRelRef : imrWeightsMap.keySet()) {
				double imrWeight = imrWeightsMap.get(attenRelRef);
				
				System.out.println("Calculating catalog losses");
				Map<Double, List<DiscretizedFunc>> lossDists =
						calc.getLossDists(attenRelRef, xAxisScale, durations, allSubDurations);
				
				myLossDists.add(lossDists);
				myWeights.add(imrWeight);
				
				Map<Double, HistogramFunction> lossHists = getLossHist(lossDists, deltaX, isLog10);
//				writeLossHist(outputDir, attenRelRef.name(), lossHist, isLog10);
				writeLossHist(outputDir, attenRelRef.name()+prefixAdd, lossHists, isLog10, triggeredOnly,
						xAxisLabel, maxX);
				Map<Double, DiscretizedFunc> exceedFuncs = toExceedFuncs(lossHists, calc.catalogs.size(),
						allSubDurations, calc.getRoundedMaxCatalogDiration());
				for (double duration : exceedFuncs.keySet())
					gmpeCombLossExceeds.put(attenRelRef.name(), duration, exceedFuncs.get(duration));
				exceedWeightMap.put(attenRelRef.name(), imrWeight);
				writeLossExceed(outputDir, attenRelRef.name()+prefixAdd, exceedFuncs, isLog10, triggeredOnly, xAxisLabel, maxX);
				
				calc.writeLossesToCSV(outputDir, attenRelRef.name()+prefixAdd+"_losses"+csvPrefixAdd, lossDists);
			}
			
			// combined for all atten rels
			if (imrWeightsMap.size() > 1) {
				Map<Double, List<DiscretizedFunc>> imrCombined = getCombinedLossDists(myLossDists, myWeights);
				Map<Double, HistogramFunction> lossHists = getLossHist(imrCombined, deltaX, isLog10);
				writeLossHist(outputDir, "gmpes_combined"+prefixAdd, lossHists, isLog10, triggeredOnly, xAxisLabel, maxX);
				writeLossExceed(outputDir, "gmpes_combined"+prefixAdd, gmpeCombLossExceeds, exceedWeightMap, isLog10, triggeredOnly, xAxisLabel, maxX, false, true);
//				writeLossExceed(outputDir, "gmpes_combined"+prefixAdd, lossHists, isLog10, triggeredOnly, xAxisLabel, maxX, calc.catalogs.size());
				calc.writeLossesToCSV(outputDir, "gmpes_combined"+prefixAdd+"_losses"+csvPrefixAdd, imrCombined);
			}
			
			lossDistsList.addAll(myLossDists);
			lossWeights.addAll(myWeights);
			for (String name : gmpeCombLossExceeds.rowKeySet()) {
				String combName = dataDir.getName()+"_"+name;
				double weight = exceedWeightMap.get(name)/(double)dataDirs.size();
				for (double duration : gmpeCombLossExceeds.columnKeySet())
					allCombLossExceeds.put(combName, duration, gmpeCombLossExceeds.get(name, duration));
				allExceedWeightMap.put(combName, weight);
			}
		}
		// combine
		if (dataDirs.size() > 1) {
			String combName = "combined";
			if (dataDirs.get(0).getName().contains("bernardino"))
				combName = "combined-san-bernardino";
			else if (dataDirs.get(0).getName().contains("coachella"))
				combName = "combined-coachella";
			File outputDir = new File(lossOutputDir, combName);
			if (!outputDir.exists())
				outputDir.mkdir();
			
			Map<Double, List<DiscretizedFunc>> combined = getCombinedLossDists(lossDistsList, lossWeights);
			Map<Double, HistogramFunction> lossHists = getLossHist(combined, deltaX, isLog10);
			writeLossHist(outputDir, "gmpes_combined"+prefixAdd, lossHists, isLog10, triggeredOnly, xAxisLabel, maxX);
			writeLossExceed(outputDir, "gmpes_combined"+prefixAdd, allCombLossExceeds, allExceedWeightMap, isLog10, triggeredOnly, xAxisLabel, maxX, false, true);
//			writeLossExceed(outputDir, "gmpes_combined"+prefixAdd, lossHists, isLog10, triggeredOnly, xAxisLabel, maxX, calc.catalogs.size());
			calc.writeLossesToCSV(outputDir, "gmpes_combined"+prefixAdd+"_losses"+csvPrefixAdd, combined);
		}
	}
	
	private static Map<Double, List<DiscretizedFunc>> getCombinedLossDists(
			List<Map<Double, List<DiscretizedFunc>>> lossDistsList, List<Double> lossWeights) {
		Map<Double, List<DiscretizedFunc>> combinedLosses = Maps.newHashMap();
		
		for (double duration : lossDistsList.get(0).keySet()) {
			List<DiscretizedFunc> list = Lists.newArrayList();
			for (int i=0; i<lossDistsList.get(0).get(duration).size(); i++) {
				list.add(new ArbitrarilyDiscretizedFunc());
			}
			combinedLosses.put(duration, list);
		}
		
		double totWeight = 0d;
		for (double weight : lossWeights)
			totWeight += weight;
		
		for (int i=0; i<lossDistsList.size(); i++) {
			Map<Double, List<DiscretizedFunc>> lossDistsMap = lossDistsList.get(i);
			double weight = lossWeights.get(i)/totWeight;
			
			for (double duration : lossDistsMap.keySet()) {
				List<DiscretizedFunc> lossDists = lossDistsMap.get(duration);
				Preconditions.checkState(lossDists.size() == combinedLosses.get(duration).size());
				
				for (int n=0; n<combinedLosses.get(duration).size(); n++) {
					DiscretizedFunc combined = combinedLosses.get(duration).get(n);
					for (Point2D pt : lossDists.get(n)) {
						double x = pt.getX();
						double y = pt.getY()*weight;
						int xInd = UCERF3_BranchAvgLossFetcher.getMatchingXIndexFloatPrecision(x, combined);
						if (xInd < 0)
							combined.set(x, y);
						else
							combined.set(x, y + combined.getY(xInd));
					}
				}
			}
		}
		
		return combinedLosses;
	}

}
