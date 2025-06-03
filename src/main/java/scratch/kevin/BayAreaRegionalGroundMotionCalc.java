package scratch.kevin;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.StringWriter;
import java.text.DecimalFormat;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.siteData.impl.ThompsonVs30_2022;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.commons.geo.Region;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.impl.WarningDoubleParameter;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.io.archive.ArchiveOutput;
import org.opensha.commons.util.io.archive.ArchiveOutput.ParallelZipFileOutput;
import org.opensha.sha.calc.params.filters.FixedDistanceCutoffFilter;
import org.opensha.sha.earthquake.DistCachedERFWrapper;
import org.opensha.sha.earthquake.DistCachedERFWrapper.DistCacheWrapperRupture;
import org.opensha.sha.earthquake.PointSource;
import org.opensha.sha.earthquake.PointSource.PoissonPointSource;
import org.opensha.sha.earthquake.PointSource.PoissonPointSourceData;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.erf.BaseFaultSystemSolutionERF;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.param.BackgroundRupType;
import org.opensha.sha.earthquake.param.UseRupMFDsParam;
import org.opensha.sha.earthquake.util.GridCellSuperSamplingPoissonPointSourceData;
import org.opensha.sha.earthquake.util.GridCellSupersamplingSettings;
import org.opensha.sha.earthquake.util.GriddedSeismicitySettings;
import org.opensha.sha.faultSurface.CompoundSurface;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.cache.CustomCacheWrappedSurface;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_WrapperFullParam;
import org.opensha.sha.imr.attenRelImpl.ngaw2.ScalarGroundMotion;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGV_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;

import scratch.kevin.prvi25.figures.MapSourceTypeDisagg;

public class BayAreaRegionalGroundMotionCalc {

	public static void main(String[] args) throws IOException {
		AttenRelRef gmmRef = AttenRelRef.BSSA_2014;
		double period = 0d;
		String perSuffix = "pga";
//		double period = -1d;
//		String perSuffix = "pgv";

		int calcThreads = 16;
		int zipThreads = 16;
		
		NGAW2_WrapperFullParam gmm = (NGAW2_WrapperFullParam)gmmRef.get();
		setPeriod(gmm, period);
		
		BaseFaultSystemSolutionERF erf = new BaseFaultSystemSolutionERF();
		erf.setSolution(FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2024_02_02-nshm23_branches-WUS_FM_v3/"
//				+ "results_WUS_FM_v3_branch_averaged_gridded_simplified.zip")));
				+ "results_WUS_FM_v3_branch_averaged_gridded.zip")));
		double cellGridSpacing = 0.1;
		
//		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.ONLY);
		erf.setParameter(UseRupMFDsParam.NAME, false);
		erf.setCacheGridSources(true);
		
//		GridCellSupersamplingSettings regionalSupersampling = GridCellSupersamplingSettings.QUICK;
		GridCellSupersamplingSettings regionalSupersampling = null;
//		Random relocateRand = null;
		Random relocateRand = new Random(123456);
		GriddedSeismicitySettings gridSettings = erf.getGriddedSeismicitySettings();
		gridSettings = gridSettings.forSupersamplingSettings(null);
		gridSettings = gridSettings.forSurfaceType(BackgroundRupType.FINITE);
		gridSettings = gridSettings.forPointSourceMagCutoff(5d);
		gridSettings = gridSettings.forMinimumMagnitude(3.5d);
		System.out.println("Gridded seismicity settings:\n"+gridSettings);
		erf.setGriddedSeismicitySettings(gridSettings);
		erf.getTimeSpan().setDuration(1d);
		
		if (regionalSupersampling != null && gridSettings.surfaceType != BackgroundRupType.POINT) {
			regionalSupersampling = new GridCellSupersamplingSettings(regionalSupersampling.targetSpacingKM,
					regionalSupersampling.fullDist, regionalSupersampling.borderDist, regionalSupersampling.cornerDist, true);
		}
		
		erf.updateForecast();
		
		System.out.println("Duration: "+erf.getTimeSpan().getDuration());
		System.out.println("ERF has "+erf.getTotNumRupsFromFaultSystem()+" fault ruptures");
		System.out.println("ERF has "+(erf.getTotNumRups()-erf.getTotNumRupsFromFaultSystem())+" gridded ruptures");
		
		DecimalFormat pDF = new DecimalFormat("0.0%");
		
		FixedDistanceCutoffFilter filter = new FixedDistanceCutoffFilter(300d);
		
		File outputDir = new File("/home/kevin/OpenSHA/minson_2025_bay_area_gms");
		
		System.out.println(gmm.getShortName()+" default params (used for phi):");
		System.out.println("ERF:");
		for (Parameter<?> param : gmm.getEqkRuptureParams())
			System.out.println("\t"+param.getName()+":\t"+param.getValue());
		System.out.println("Propogation:");
		for (Parameter<?> param : gmm.getPropagationEffectParams())
			System.out.println("\t"+param.getName()+":\t"+param.getValue());
		
		Region reg = new CaliforniaRegions.CYBERSHAKE_BAY_AREA_MAP_REGION();
		GriddedRegion gridReg = new GriddedRegion(reg, 0.02, GriddedRegion.ANCHOR_0_0);
		System.out.println("Region has "+gridReg.getNodeCount()+" nodes");
		
		List<Site> filterSites = new ArrayList<>();
		for (Location loc : reg.getBorder())
			filterSites.add(new Site(loc));
		
		GriddedGeoDataSet landMask = MapSourceTypeDisagg.buildLandMask(gridReg);
		
//		ThompsonVs30_2020 vs30model = new ThompsonVs30_2020();
		ThompsonVs30_2022 vs30model = new ThompsonVs30_2022(new File(outputDir, "California_vs30_Wills15_hybrid.flt").getAbsolutePath());
		
		ArrayList<Double> vsValues = vs30model.getValues(gridReg.getNodeList());
		
		int numWithBoth = 0;
		int numVs30only = 0;
		int numLandOnly = 0;
		
		List<Site> sites = new ArrayList<>();
		
		for (int i=0; i<gridReg.getNodeCount(); i++) {
			Double vs30 = vsValues.get(i);
			boolean vsValid = vs30model.isValueValid(vs30);
			boolean onLand = landMask.get(i) > 0d;
			if (vsValid && onLand) {
				numWithBoth++;
				Site site = new Site(gridReg.getLocation(i));
				for (Parameter siteParam : gmm.getSiteParams()) {
					if (siteParam.getName().equals(Vs30_Param.NAME)) {
						siteParam = (Parameter) siteParam.clone();
						((WarningDoubleParameter)siteParam).setValueIgnoreWarning(vs30);
					}
					site.addParameter(siteParam);
				}
				sites.add(site);
			} else if (vsValid) {
				numVs30only++;
			} else if (onLand) {
				numLandOnly++;
			}
		}
		
		System.out.println("Have both: "+numWithBoth+" ("+pDF.format((double)numWithBoth/(double)gridReg.getNodeCount())+")");
		System.out.println("Have Vs30 only: "+numVs30only);
		System.out.println("Have land only: "+numLandOnly);
		
		CSVFile<String> sitesCSV = new CSVFile<>(true);
		sitesCSV.addLine("Site index", "Latitude", "Longitude", "Vs30 (m/s)", "phi at r=0");
		for (int i=0; i<sites.size(); i++) {
			Site site = sites.get(i);
			gmm.setSite(site);
			ScalarGroundMotion gm = gmm.getGroundMotion();
			sitesCSV.addLine(i+"", (float)site.getLocation().lat+"", (float)site.getLocation().lon+"",
					((Double)site.getParameter(Vs30_Param.NAME).getValue()).floatValue()+"", (float)gm.phi()+"");
		}
		sitesCSV.writeToFile(new File(outputDir, "sites_and_phi_"+perSuffix+".csv"));
		

		CSVFile<String> eventDetailsCSV = new CSVFile<>(true);
		eventDetailsCSV.addLine("Event ID", "Original Source ID", "Original Rupture ID",
				"Magnitude", "Annual Rate", "Rake", "Dip", "Length (km)", "Down-Dip Width (km)", "Upper Depth (km)",
				"FSS Index", "Grid Node Index", "Gridded Rupture Latitude", "Gridded Rupture Longitude", "Gridded Rupture Strike");
		
		List<ProbEqkRupture> keptEvents = new ArrayList<>();
		int keptFault = 0;
		int keptGridded = 0;
		int numSupersampled = 0;
		int numRelocated = 0;
		GridSourceProvider gridSources = erf.getGridSourceProvider();
		for (int s=0; s<erf.getNumSources(); s++) {
			ProbEqkSource source = erf.getSource(s);
			boolean keepSource = false;
			double minDist = Double.POSITIVE_INFINITY;
			for (Site filterSite : filterSites) {
				double dist = source.getMinDistance(filterSite);
				minDist = Math.min(dist, minDist);
				if (!filter.canSkipSource(source, filterSite, dist))
					keepSource = true;
			}
			if (!keepSource)
				continue;
			boolean gridded = source instanceof PointSource;
			if (regionalSupersampling != null && gridded) {
				PoissonPointSource pointSource = (PoissonPointSource)source;
				PoissonPointSourceData data = pointSource.getData();
				Location center = pointSource.getLocation();
				double halfSpacing = cellGridSpacing*0.5;
				Region gridCell = new Region(new Location(center.lat-halfSpacing, center.lon-halfSpacing),
						new Location(center.lat+halfSpacing, center.lon+halfSpacing));
				GridCellSuperSamplingPoissonPointSourceData samplingData =
						new GridCellSuperSamplingPoissonPointSourceData(data, center, gridCell, regionalSupersampling);
				PoissonPointSourceData sampledData = samplingData.getForDistance(minDist);
				if (sampledData != data) {
					source = new PoissonPointSource(center, source.getTectonicRegionType(),
							pointSource.getDuration(), sampledData, pointSource.getDistCorrs());
					numSupersampled++;
//					System.out.println("Supersampling with minDist="+(float)minDist);
//				} else if (minDist < 50d){
//					System.out.println("Skipping supersampling with minDist="+(float)minDist+", loc="+center);
				}
			}
			List<ProbEqkRupture> rups = source.getRuptureList();
			List<Location> griddedLocs = null;
			if (gridded) {
				griddedLocs = new ArrayList<>(rups.size());
				PoissonPointSource pointSource = (PoissonPointSource)source;
				Location center = pointSource.getLocation();
				if (relocateRand != null) {
					// randomly relocate within the grid cell
					for (int i=0; i<rups.size(); i++) {
						ProbEqkRupture origRup = rups.get(i);
						RuptureSurface surface = origRup.getRuptureSurface();
						double newLat   = center.lat + (relocateRand.nextDouble() - 0.5)*cellGridSpacing;
						double newLon   = center.lon + (relocateRand.nextDouble() - 0.5)*cellGridSpacing;
						Location newLoc = new Location(newLat, newLon);
						griddedLocs.add(newLoc);
						LocationVector vector = LocationUtils.vector(center, newLoc);
						surface = surface.getMoved(vector);
						Location hypo = origRup.getHypocenterLocation();
						if (hypo != null)
							hypo = LocationUtils.location(hypo, vector);
						rups.set(i, new ProbEqkRupture(origRup.getMag(), origRup.getAveRake(), origRup.getProbability(), surface, hypo));
						numRelocated++;
					}
				} else {
					for (int i=0; i<rups.size(); i++)
						griddedLocs.add(center);
				}
			}
			for (int r=0; r<rups.size(); r++) {
				ProbEqkRupture rup = rups.get(r);
				List<String> line = new ArrayList<>(eventDetailsCSV.getNumCols());
//				eventDetailsCSV.addLine("Event ID", "Original Source ID", "Original Rupture ID",
//						"Magnitude", "Annual Rate", "Rake", "Dip", "Length (km)", "Down-Dip Width (km)", "Upper Depth (km)",
//						"FSS Index", "Grid Node Index", "Gridded Rupture Latitude", "Gridded Rupture Longitude", "Gridded Rupture Strike");
				line.add(keptEvents.size()+"");
				line.add(s+"");
				line.add(r+"");
				line.add((float)rup.getMag()+"");
				line.add((float)rup.getMeanAnnualRate(1d)+"");
				line.add((float)rup.getAveRake()+"");
				RuptureSurface surf = rup.getRuptureSurface();
				line.add((float)surf.getAveDip()+"");
				line.add((float)surf.getAveLength()+"");
				line.add((float)surf.getAveWidth()+"");
				line.add((float)surf.getAveRupTopDepth()+"");
				if (gridded) {
					line.add("-1");
					line.add(gridSources.getLocationIndex(((PoissonPointSource)source).getLocation())+"");
					Location loc = griddedLocs.get(r);
					line.add((float)loc.lat+"");
					line.add((float)loc.lon+"");
					if (surf instanceof PointSurface)
						line.add("NaN");
					else
						line.add((float)surf.getAveStrike()+"");
				} else {
					line.add(erf.getFltSysRupIndexForSource(s)+"");
					line.add("-1");
					line.add("");
					line.add("");
					line.add("");
				}
				eventDetailsCSV.addLine(line);
				keptEvents.add(rup);
			}
			if (gridded)
				keptGridded += rups.size();
			else
				keptFault += rups.size();
		}
		eventDetailsCSV.writeToFile(new File(outputDir, "event_details.csv"));
		System.exit(0);
		if (regionalSupersampling != null)
			System.out.println("Supersampled "+numSupersampled+"/"+(erf.getNumSources()-erf.getNumFaultSystemSources())+" gridded sources");
		if (relocateRand != null)
			System.out.println("Relocated "+numRelocated+" gridded ruptures");
		erf = null; // don't use me!
		System.out.println("Kept "+keptEvents.size()+" events");
		System.out.println("Kept "+keptFault+" fault events");
		System.out.println("Kept "+keptGridded+" grid events");
		int bundleSize = 1000;
		List<List<ProbEqkRupture>> eventBundles = new ArrayList<>();
		List<ProbEqkRupture> curBundle = null;
		List<List<Integer>> eventIDBundles = new ArrayList<>();
		List<Integer> curIDBundle = null;
		CSVFile<String> eventsCSV = new CSVFile<>(true);
		eventsCSV.addLine("Event ID", "Magnitude", "Annual Rate", "tau");
		for (int n=0; n<keptEvents.size(); n++) {
			ProbEqkRupture rup = keptEvents.get(n);
			gmm.setEqkRupture(rup);
			ScalarGroundMotion gm = gmm.getGroundMotion();
			eventsCSV.addLine(n+"", (float)rup.getMag()+"", rup.getMeanAnnualRate(1d)+"", (float)gm.tau()+"");
			
			if (curBundle == null) {
				curBundle = new ArrayList<>();
				eventBundles.add(curBundle);
				curIDBundle = new ArrayList<>();
				eventIDBundles.add(curIDBundle);
			}
			curBundle.add(rup);
			curIDBundle.add(n);
			if (curBundle.size() == bundleSize) {
				curBundle = null;
				curIDBundle = null;
			}
		}
		eventsCSV.writeToFile(new File(outputDir, "events_and_tau_"+perSuffix+".csv"));
//		System.out.println("Will keep "+eventsToKeep.cardinality()+"/"+totNumRups
//				+" ("+pDF.format((double)eventsToKeep.cardinality()/(double)totNumRups)+")");
		System.out.println("Will do "+eventBundles.size()+" bundles");
		
		ExecutorService exec = Executors.newFixedThreadPool(calcThreads);
		
		List<CalcCallable> calls = new ArrayList<>();
		for (int i=0; i<calcThreads; i++)
			calls.add(new CalcCallable(gmmRef, period, sites));
		
		CompletableFuture<Void> writeFuture = null;
		
		File outZipFile = new File(outputDir, "event_mu_maps_"+perSuffix+".zip");
		ParallelZipFileOutput outData =
//				new ArchiveOutput.AsynchronousZipFileOutput(outZipFile);
				new ArchiveOutput.ParallelZipFileOutput(outZipFile, zipThreads);
		
		outData.setTrackBlockingTimes(true);
		
		Stopwatch totalWatch = Stopwatch.createStarted();
		Stopwatch calcWatch = Stopwatch.createUnstarted();
		Stopwatch mapWatch = Stopwatch.createUnstarted();
		Stopwatch ioWatch = Stopwatch.createUnstarted();
		
		double millisToMinutes = 1d/(1000d*60d);
		
		for (int b=0; b<eventBundles.size(); b++) {
			List<ProbEqkRupture> bundle = eventBundles.get(b);
			List<Integer> idBundle = eventIDBundles.get(b);
			Preconditions.checkState(idBundle.size() == bundle.size());
			List<Future<List<SiteResult>>> futures = new ArrayList<>(calcThreads);
			
//			if (b == 10)
//				break;
			
			calcWatch.start();
			ArrayDeque<Integer> sitesDeque = new ArrayDeque<>(sites.size());
			for (int i=0; i<sites.size(); i++)
				sitesDeque.add(i);
			
			for (CalcCallable call : calls) {
				call.setEvents(bundle, sitesDeque);
				futures.add(exec.submit(call));
			}
			
			SiteResult[] combResults = new SiteResult[sites.size()];
			for (Future<List<SiteResult>> future : futures) {
				try {
					for (SiteResult result : future.get())
						combResults[result.siteIndex] = result;
				} catch (Exception e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
			}
			calcWatch.stop();
			
			mapWatch.start();
			List<Future<byte[]>> mapFutures = new ArrayList<>(bundle.size());
			for (int e=0; e<bundle.size(); e++) {
				int index = e;
				mapFutures.add(exec.submit(new Callable<byte[]>() {

					@Override
					public byte[] call() throws Exception {
						int size = Integer.max(1000, sites.size()*10);
						StringWriter stringWriter = new StringWriter(size);
						for (int s=0; s<combResults.length; s++) {
							Preconditions.checkNotNull(combResults[s]);
							stringWriter.write((float)combResults[s].muValues[index]+"\n");
						}
						stringWriter.flush();
						return stringWriter.toString().getBytes();
					}
				}));
			}
			
			List<byte[]> mapBuffers = new ArrayList<>(bundle.size());
			for (Future<byte[]> future : mapFutures) {
				try {
                    mapBuffers.add(future.get());
                } catch (Exception e) {
                    throw ExceptionUtils.asRuntimeException(e);
                }
			}
			
			mapWatch.stop();
			
			ioWatch.start();
			if (writeFuture != null)
				writeFuture.join();
			
			writeFuture = CompletableFuture.runAsync(new Runnable() {
				
				@Override
				public void run() {
					try {
						for (int e=0; e<bundle.size(); e++) {
							int eventID = idBundle.get(e);

							outData.putNextEntry(eventID+".txt");
							
							OutputStream stream = outData.getOutputStream();
							stream.write(mapBuffers.get(e));
							stream.flush();
							
							outData.closeEntry();
						}
					} catch (IOException e1) {
						throw ExceptionUtils.asRuntimeException(e1);
					}
				}
			});
			
			ioWatch.stop();
			
			double totalMins = totalWatch.elapsed(TimeUnit.MILLISECONDS)*millisToMinutes;
			double calcMins = calcWatch.elapsed(TimeUnit.MILLISECONDS)*millisToMinutes;
			double mapMins = mapWatch.elapsed(TimeUnit.MILLISECONDS)*millisToMinutes;
			double ioMins = ioWatch.elapsed(TimeUnit.MILLISECONDS)*millisToMinutes;
			
			System.out.println("Done with bundle "+b+"/"+eventBundles.size()+" ("+pDF.format((double)b/(double)eventBundles.size())+")"
					+" in "+(float)totalMins+" m; "
					+pDF.format(calcMins/totalMins)+" calculating; "
					+pDF.format(mapMins/totalMins)+" building maps; "
					+pDF.format(ioMins/totalMins)+" blocking I/O");
			System.out.println("\t"+outData.getBlockingTimeStats());
		}
		
		ioWatch.start();
		writeFuture.join();
		
		outData.close();
		ioWatch.stop();
		
		double totalMins = totalWatch.elapsed(TimeUnit.MILLISECONDS)*millisToMinutes;
		double calcMins = calcWatch.elapsed(TimeUnit.MILLISECONDS)*millisToMinutes;
		double mapMins = mapWatch.elapsed(TimeUnit.MILLISECONDS)*millisToMinutes;
		double ioMins = ioWatch.elapsed(TimeUnit.MILLISECONDS)*millisToMinutes;
		
		System.out.println("DONE in "+(float)totalMins+" m; "
				+pDF.format(calcMins/totalMins)+" calculating; "
				+pDF.format(mapMins/totalMins)+" building maps; "
				+pDF.format(ioMins/totalMins)+" blocking I/O");
		
		exec.shutdown();
	}
	
	private static void setPeriod(ScalarIMR gmm, double period) {
		if (period == 0d) {
			gmm.setIntensityMeasure(PGA_Param.NAME);
		} else if (period == -1) {
			gmm.setIntensityMeasure(PGV_Param.NAME);
		} else {
			Preconditions.checkState(period > 0d);
			gmm.setIntensityMeasure(SA_Param.NAME);
			SA_Param.setPeriodInSA_Param(gmm.getIntensityMeasure(), period);
		}
	}
	
	private static class CalcCallable implements Callable<List<SiteResult>> {
		private NGAW2_WrapperFullParam gmm;
		private List<Site> sites;
		private Map<RuptureSurface, CustomCacheWrappedSurface> wrappedMap;
		
		// inputs for each batch
		private List<ProbEqkRupture> events;
		private ArrayDeque<Integer> sitesDeque;
		
		// outputs for each batch
//		private List
		
		public CalcCallable(AttenRelRef gmmRef, double period, List<Site> sites) {
			gmm = (NGAW2_WrapperFullParam)gmmRef.get();
			setPeriod(gmm, period);
			this.sites = sites;
			this.wrappedMap = new HashMap<>();
		}
		
		public void setEvents(List<ProbEqkRupture> events, ArrayDeque<Integer> sitesDeque) {
			this.events = events;
			this.sitesDeque = sitesDeque;
		}

		@Override
		public List<SiteResult> call() throws Exception {
			Preconditions.checkNotNull(events);
			Preconditions.checkNotNull(sitesDeque);
			
			List<ProbEqkRupture> rups = new ArrayList<>(events.size());
			for (ProbEqkRupture rup : events) {
				RuptureSurface surf = rup.getRuptureSurface();
				if (surf instanceof CompoundSurface) {
					ProbEqkRupture origRup = rup;
					RuptureSurface wrappedSurf = DistCachedERFWrapper.getWrappedSurface(wrappedMap, surf);
					rup = new DistCacheWrapperRupture(rup, wrappedSurf);
					// sanity check
					Preconditions.checkState(origRup.getMag() == rup.getMag());
					Preconditions.checkState(origRup.getProbability() == rup.getProbability());
					Preconditions.checkState(origRup.getAveRake() == rup.getAveRake());
				}
				rups.add(rup);
			}
			
			List<SiteResult> results = new ArrayList<>();
			while (true) {
				int siteIndex;
				synchronized (sitesDeque) {
					if (sitesDeque.isEmpty())
						break;
					siteIndex = sitesDeque.pop();
				}
				
				gmm.setSite(sites.get(siteIndex));
				double[] muValues = new double[rups.size()];
				for (int i=0; i<muValues.length; i++) {
					ProbEqkRupture rup = rups.get(i);
					gmm.setEqkRupture(rup);
					muValues[i] = gmm.getGroundMotion().mean();
				}
				results.add(new SiteResult(siteIndex, muValues));
			}
			
			sitesDeque = null;
			events = null;
			return results;
		}
	}
	
	private static class SiteResult {
		final int siteIndex;
		final double[] muValues;
		private SiteResult(int siteIndex, double[] muValues) {
			super();
			this.siteIndex = siteIndex;
			this.muValues = muValues;
		}
	}

}
