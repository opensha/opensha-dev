package scratch.kevin;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.StringWriter;
import java.text.DecimalFormat;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.util.Precision;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.siteData.impl.ThompsonVs30_2022;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.LocationVector;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.impl.WarningDoubleParameter;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.io.archive.ArchiveOutput;
import org.opensha.commons.util.io.archive.ArchiveOutput.ParallelZipFileOutput;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.params.filters.FixedDistanceCutoffFilter;
import org.opensha.sha.calc.params.filters.SourceFilterManager;
import org.opensha.sha.calc.params.filters.SourceFilters;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.DistCachedERFWrapper;
import org.opensha.sha.earthquake.DistCachedERFWrapper.DistCacheWrapperRupture;
import org.opensha.sha.earthquake.PointSource;
import org.opensha.sha.earthquake.PointSource.PoissonPointSource;
import org.opensha.sha.earthquake.PointSource.PoissonPointSourceData;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.erf.BaseFaultSystemSolutionERF;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultGridAssociations;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.earthquake.observedEarthquake.parsers.UCERF3_CatalogParser;
import org.opensha.sha.earthquake.param.BackgroundRupType;
import org.opensha.sha.earthquake.param.UseRupMFDsParam;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.LocalRegions;
import org.opensha.sha.earthquake.util.GridCellSuperSamplingPoissonPointSourceData;
import org.opensha.sha.earthquake.util.GridCellSupersamplingSettings;
import org.opensha.sha.earthquake.util.GriddedSeismicitySettings;
import org.opensha.sha.faultSurface.ApproxEvenlyGriddedSurface;
import org.opensha.sha.faultSurface.CompoundSurface;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.GriddedSurfaceImpl;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.cache.CustomCacheWrappedSurface;
import org.opensha.sha.gcim.ui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_WrapperFullParam;
import org.opensha.sha.imr.attenRelImpl.ngaw2.ScalarGroundMotion;
import org.opensha.sha.imr.attenRelImpl.nshmp.NSHMP_GMM_Wrapper;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGV_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.SiteParams.DepthTo1pt0kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.DepthTo2pt5kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;

import net.mahdilamb.colormap.Colors;
import scratch.UCERF3.erf.ETAS.association.JeanneFileLoader;
import scratch.kevin.prvi25.figures.MapSourceTypeDisagg;

public class BayAreaRegionalGroundMotionCalc {

	public static void main(String[] args) throws IOException {
		File mainOutputDir = new File("/home/kevin/OpenSHA/minson_2025_bay_area_gms");
		
		AttenRelRef gmmRef = AttenRelRef.BSSA_2014;
		AttenRelRef empGMMRef = gmmRef;
//		empGMMRef = AttenRelRef.WRAPPED_BSSA_2014;
		AttenRelRef hazMapGMMRef = AttenRelRef.USGS_NSHM23_ACTIVE_SF;
		double period = 0d;
		String perSuffix = "pga";
		String perLabel = "PGA (g)";
//		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(PGA_Param.NAME);
//		DiscretizedFunc xVals = new EvenlyDiscretizedFunc(0.1, 1d, 10);
//		DiscretizedFunc xVals = new EvenlyDiscretizedFunc(0.0, 1d, 11);
		DiscretizedFunc xVals = new EvenlyDiscretizedFunc(0.0, 1d, 101);
//		double period = -1d;
//		String perSuffix = "pgv";
		
		File outputDir = new File(mainOutputDir, "nshm23_model");
		File catalogFile = null;
		double catYears = Double.NaN;
		double catalogDefaultRake = Double.NaN;
		double gridMinMag = 3.5d;
		double[] epiMinMags = {3.5d, 4d};
		boolean doEventMapCalc = false;
		boolean doHazardMapCalc = true;
		
//		File outputDir = new File(mainOutputDir, "nshm23_catalog");
//		File catalogFile = new File(mainOutputDir, "nshm23_cat_M4_1967.csv");
//		double catYears = 2023-1967; // ends at start of 2023
//		double catalogDefaultRake = Double.NaN;
////		double catalogDefaultRake = 0;
//		double gridMinMag = Double.NaN;
//		double[] epiMinMags = {4d};
//		boolean doEventMapCalc = true;
//		boolean doHazardMapCalc = false;
		
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());

		int calcThreads = 30;
		int zipThreads = 8;
		
		NGAW2_WrapperFullParam gmm = (NGAW2_WrapperFullParam)gmmRef.get();
		setPeriod(gmm, period);
		
		DecimalFormat pDF = new DecimalFormat("0.0%");
		DecimalFormat oDF = new DecimalFormat("0.##");
		
		System.out.println(gmm.getShortName()+" default params (used for phi):");
		System.out.println("ERF:");
		for (Parameter<?> param : gmm.getEqkRuptureParams())
			System.out.println("\t"+param.getName()+":\t"+param.getValue());
		System.out.println("Propogation:");
		for (Parameter<?> param : gmm.getPropagationEffectParams())
			System.out.println("\t"+param.getName()+":\t"+param.getValue());
		
		SourceFilterManager filters = new SourceFilterManager(SourceFilters.FIXED_DIST_CUTOFF);
		FixedDistanceCutoffFilter filter = filters.getFilterInstance(FixedDistanceCutoffFilter.class);
		filter.setMaxDistance(300d);
		
		Region reg = LocalRegions.CONUS_SF_BAY.load();
		GriddedRegion gridReg = new GriddedRegion(reg, 0.02, GriddedRegion.ANCHOR_0_0);
		System.out.println("Region has "+gridReg.getNodeCount()+" nodes");
		
		List<Site> filterSites = new ArrayList<>();
		for (Location loc : reg.getBorder())
			filterSites.add(new Site(loc));
		
		GriddedGeoDataSet landMask = MapSourceTypeDisagg.buildLandMask(gridReg);
		
//		ThompsonVs30_2020 vs30model = new ThompsonVs30_2020();
		ThompsonVs30_2022 vs30model = new ThompsonVs30_2022(new File(mainOutputDir, "California_vs30_Wills15_hybrid.flt").getAbsolutePath());
		
		ArrayList<Double> vsValues = vs30model.getValues(gridReg.getNodeList());
		
		int numWithBoth = 0;
		int numVs30only = 0;
		int numLandOnly = 0;
		
		List<Site> sites = new ArrayList<>();
		List<Site> sitesAtNodes = new ArrayList<>();
		
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
				sitesAtNodes.add(site);
			} else {
				sitesAtNodes.add(null);
				if (vsValid) {
					numVs30only++;
				} else if (onLand) {
					numLandOnly++;
				}
			}
		}
		Preconditions.checkState(sitesAtNodes.size() == gridReg.getNodeCount());
		
		System.out.println("Have both: "+numWithBoth+" ("+pDF.format((double)numWithBoth/(double)gridReg.getNodeCount())+")");
		System.out.println("Have Vs30 only: "+numVs30only);
		System.out.println("Have land only: "+numLandOnly);
		
		List<Site> basinSites = new ArrayList<>();
		GriddedGeoDataSet[] basins = loadBasinDepths(new File(mainOutputDir, "bay-area.csv"));
		int numNoBasins = 0;
		for (Site site : sites) {
			site = (Site) site.clone(); // this also clones params, e.g., basin
			Location loc = site.getLocation();
			int locIndex = basins[0].indexOf(loc);
			if (locIndex < 0) {
				numNoBasins++;;
			} else {
				double z1 = basins[0].get(locIndex);
				double z25 = basins[1].get(locIndex);
				if (Double.isNaN(z1)) {
					numNoBasins++;
				} else {
					site.getParameter(DepthTo1pt0kmPerSecParam.NAME).setValue(z1);
					site.getParameter(DepthTo2pt5kmPerSecParam.NAME).setValue(z25);
				}
			}
			basinSites.add(site);
		}
		System.out.println(numNoBasins+"/"+basinSites.size()+" sites were missing basin depth data");
		
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
		
		List<ProbEqkRupture> keptEvents;
		List<Double> keptFractsInReg;
		AbstractERF hazardMapERF = null;
		if (catalogFile == null) {
			BaseFaultSystemSolutionERF erf = new BaseFaultSystemSolutionERF();
			erf.setSolution(FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2024_02_02-nshm23_branches-WUS_FM_v3/"
//					+ "results_WUS_FM_v3_branch_averaged_gridded_simplified.zip")));
					+ "results_WUS_FM_v3_branch_averaged_gridded.zip")));
			double cellGridSpacing = 0.1;
			
//			erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.ONLY);
			erf.setParameter(UseRupMFDsParam.NAME, false);
			erf.setCacheGridSources(true);
			
//			GridCellSupersamplingSettings regionalSupersampling = GridCellSupersamplingSettings.QUICK;
			GridCellSupersamplingSettings regionalSupersampling = null;
//			Random relocateRand = null;
			Random relocateRand = new Random(123456);
			GriddedSeismicitySettings origGridSettings = erf.getGriddedSeismicitySettings();
			GriddedSeismicitySettings gridSettings = erf.getGriddedSeismicitySettings();
			gridSettings = gridSettings.forSupersamplingSettings(null);
			gridSettings = gridSettings.forSurfaceType(BackgroundRupType.FINITE);
			gridSettings = gridSettings.forPointSourceMagCutoff(5d);
			gridSettings = gridSettings.forMinimumMagnitude(gridMinMag);
			System.out.println("Gridded seismicity settings:\t"+gridSettings);
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

			CSVFile<String> eventDetailsCSV = new CSVFile<>(true);
			eventDetailsCSV.addLine("Event ID", "Original Source ID", "Original Rupture ID",
					"Magnitude", "Annual Rate", "Rake", "Dip", "Length (km)", "Down-Dip Width (km)", "Upper Depth (km)",
					"FSS Index", "Grid Node Index", "Gridded Rupture Latitude", "Gridded Rupture Longitude", "Gridded Rupture Strike",
					"Fraction in Region");

			FaultSystemRupSet rupSet = erf.getSolution().getRupSet();
			double[] fssRupFracts = rupSet.getFractRupsInsideRegion(reg, false);
			keptEvents = new ArrayList<>();
			keptFractsInReg = new ArrayList<>();
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
						source = new PoissonPointSource(center,
								pointSource.getDuration(), sampledData, pointSource.getDistCorr(), gridSettings.pointSourceMagnitudeCutoff);
						numSupersampled++;
//						System.out.println("Supersampling with minDist="+(float)minDist);
//					} else if (minDist < 50d){
//						System.out.println("Skipping supersampling with minDist="+(float)minDist+", loc="+center);
					}
				}
				List<ProbEqkRupture> rups = source.getRuptureList();
				List<Location> griddedLocs = null;
				List<Double> fracts = new ArrayList<>(rups.size());
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
							fracts.add(reg.contains(newLoc) ? 1d : 0d);
							griddedLocs.add(newLoc);
							LocationVector vector = LocationUtils.vector(center, newLoc);
							surface = surface.getMoved(vector);
							Location hypo = origRup.getHypocenterLocation();
							Preconditions.checkNotNull(hypo);
							hypo = LocationUtils.location(hypo, vector);
							rups.set(i, new ProbEqkRupture(origRup.getMag(), origRup.getAveRake(), origRup.getProbability(), surface, hypo));
							numRelocated++;
						}
					} else {
						boolean inside = reg.contains(center);
						for (int i=0; i<rups.size(); i++) {
							griddedLocs.add(center);
							fracts.add(inside ? 1d : 0d);
						}
					}
				} else {
					int fssIndex = erf.getFltSysRupIndexForSource(s);
					for (int r=0; r<rups.size(); r++)
						fracts.add(fssRupFracts[fssIndex]);
				}
				for (int r=0; r<rups.size(); r++) {
					ProbEqkRupture rup = rups.get(r);
					List<String> line = new ArrayList<>(eventDetailsCSV.getNumCols());
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
					double fract = fracts.get(r);
					line.add((float)fract+"");
					eventDetailsCSV.addLine(line);
					keptEvents.add(rup);
					keptFractsInReg.add(fract);
				}
				if (gridded)
					keptGridded += rups.size();
				else
					keptFault += rups.size();
			}
			eventDetailsCSV.writeToFile(new File(outputDir, "event_details.csv"));
//			System.exit(0);
			if (regionalSupersampling != null)
				System.out.println("Supersampled "+numSupersampled+"/"+(erf.getNumSources()-erf.getNumFaultSystemSources())+" gridded sources");
			if (relocateRand != null)
				System.out.println("Relocated "+numRelocated+" gridded ruptures");
			Preconditions.checkState(keptEvents.size() == keptFractsInReg.size());
			
			System.out.println("Kept "+keptEvents.size()+" events");
			System.out.println("Kept "+keptFault+" fault events");
			System.out.println("Kept "+keptGridded+" grid events");
			
			// now set up ERF for standard NSHM23
			System.out.println("NSHM23 gridded seismicity settings:\t"+origGridSettings);
			erf.setGriddedSeismicitySettings(origGridSettings);
			erf.updateForecast();
			hazardMapERF = erf;
			
			// write nucleation files
			GridSourceProvider gridProv = erf.getGridSourceProvider();
			GriddedRegion nuclRateRegion = new GriddedRegion(reg, 0.1, gridProv.getLocation(0));
			double[] nuclMags = {3.5d, 4d, 4.5, 5d, 5.5, 6d, 6.5, 7d, 7.5d};
			GriddedGeoDataSet[] nuclXYZs = new GriddedGeoDataSet[nuclMags.length];
			for (int m=0; m<nuclMags.length; m++)
				nuclXYZs[m] = new GriddedGeoDataSet(nuclRateRegion);
			
			for (int l=0; l<gridProv.getNumLocations(); l++) {
				Location loc = gridProv.getLocation(l);
				int index = nuclRateRegion.indexForLocation(loc);
				if (index >= 0) {
					IncrementalMagFreqDist mfd = gridProv.getMFD(l);
//					if (mfd == null)
//						continue;
					for (int m=0; m<nuclMags.length; m++) {
						int magIndex = mfd.getClosestXIndex(nuclMags[m]+0.01);
						Preconditions.checkState(mfd.getX(magIndex) > nuclMags[m] && mfd.getX(magIndex) < nuclMags[m]+0.1);
						nuclXYZs[m].add(index, mfd.getCumRate(magIndex));
					}
				}
			}
			
			// now add fault ruptures
			FaultSystemSolution sol = erf.getSolution();
			FaultGridAssociations assoc = rupSet.requireModule(FaultGridAssociations.class);
			
			double[] sectAreas = rupSet.getAreaForAllSections();
			double[] rupAreas = rupSet.getAreaForAllRups();
			for (int rupIndex=0; rupIndex<rupSet.getNumRuptures(); rupIndex++) {
				double sumArea = 0d;
				List<Integer> sects = rupSet.getSectionsIndicesForRup(rupIndex);
				
				double mag = rupSet.getMagForRup(rupIndex);
				double rate = sol.getRateForRup(rupIndex);
				for (int sectIndex : sects) {
					double fractArea = sectAreas[sectIndex]/rupAreas[rupIndex];
					sumArea += sectAreas[sectIndex];
					
					Map<Integer, Double> nodeFracts = assoc.getNodeFractions(rupIndex);
					for (int origNodeIndex : nodeFracts.keySet()) {
						int index = nuclRateRegion.indexForLocation(gridProv.getLocation(origNodeIndex));
						if (index >= 0) {
							double nodeFract = nodeFracts.get(origNodeIndex);
							double fractNucl = rate * fractArea * nodeFract;
							Preconditions.checkState(Double.isFinite(fractNucl));
							for (int m=0; m<nuclMags.length; m++)
								if ((float)mag >= (float)nuclMags[m])
									nuclXYZs[m].add(index, fractNucl);
						}
					}
				}
				Preconditions.checkState(Precision.equals(sumArea, rupAreas[rupIndex], 1e-2));
			}
			
			CSVFile<String> nuclCSV = new CSVFile<>(true);
			List<String> nuclHeader = new ArrayList<>();
			nuclHeader.add("Latitude");
			nuclHeader.add("Longitude");
			for (double nuclMag : nuclMags)
				nuclHeader.add("M>"+(float)nuclMag);
			nuclCSV.addLine(nuclHeader);
			
			GeographicMapMaker mapMaker = new GeographicMapMaker(nuclRateRegion);
			mapMaker.setWriteGeoJSON(false);
			
			for (int m=0; m<nuclMags.length; m++) {
				System.out.println("M"+nuclMags[m]+" nucl range: "+nuclXYZs[m].getMinZ()+", "+nuclXYZs[m].getMaxZ());
				CPT cpt = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(
						Math.floor(Math.log10(nuclXYZs[m].getMinZ())), Math.ceil(Math.log10(nuclXYZs[m].getMaxZ())));
				cpt.setLog10(true);
				mapMaker.plotXYZData(nuclXYZs[m], cpt, "M>"+nuclMags[m]+" nucleation rate");
				
				mapMaker.plot(outputDir, "nucleation_rates_m"+oDF.format(nuclMags[m]), " ");
			}
			
			for (int i=0; i<nuclRateRegion.getNodeCount(); i++) {
				List<String> line = new ArrayList<>(nuclHeader.size());
				Location loc = nuclRateRegion.locationForIndex(i);
				line.add((float)loc.lat+"");
				line.add((float)loc.lon+"");
				for (int m=0; m<nuclMags.length; m++)
					line.add((float)nuclXYZs[m].get(i)+"");
				nuclCSV.addLine(line);
			}
			nuclCSV.writeToFile(new File(outputDir, "nucleation_rates.csv"));
			System.exit(0);
		} else {
			keptEvents = new ArrayList<>();
			keptFractsInReg = new ArrayList<>();
			
			double rateEach = 1d/catYears;
			double probEach = 1d - Math.exp(-rateEach);
			CSVFile<String> catalog = CSVFile.readFile(catalogFile, true);
			for (int row=1; row<catalog.getNumRows(); row++) {
				double lat = catalog.getDouble(row, 1);
				double lon = catalog.getDouble(row, 2);
				double depth = catalog.getDouble(row, 3);
				double mag = catalog.getDouble(row, 4);
				double rake = catalogDefaultRake;
				Location loc = new Location(lat, lon, depth);
				Preconditions.checkState(reg.contains(loc));
				if (mag >= 5d)
					System.out.println("M"+(float)mag+" at "+catalog.get(row, 0));
				keptFractsInReg.add(1d);
				RuptureSurface surf;
				if ((float)mag == 6.89f) {
					System.out.println("Attaching Loma Prieta finite surface");
					File finiteDir = new File("/home/kevin/workspace/opensha/src/main/resources/data/erf/ucerf3/EarthquakeCatalog");
					File finiteFile = new File(finiteDir, "UCERF3_finite.dat");
					ObsEqkRupList inputRups = UCERF3_CatalogParser.loadCatalog(new File(finiteDir, "ofr2013-1165_EarthquakeCat.txt"));
					
					List<ObsEqkRupture> finiteRups = JeanneFileLoader.loadFiniteRups(finiteFile, inputRups);
					surf = null;
					for (ObsEqkRupture rup : finiteRups) {
						if ((float)rup.getMag() == (float)mag && LocationUtils.horzDistanceFast(loc, rup.getHypocenterLocation()) < 1d) {
							surf = rup.getRuptureSurface();
							System.out.println("Found it! Surface type: "+surf.getClass().getName()+", mag: "+(float)rup.getMag()+", id: "+rup.getEventId()+", rake: "+rup.getAveRake());
							if (Double.isFinite(rup.getAveRake()))
								rake = rup.getAveRake();
							if (surf instanceof GriddedSurfaceImpl) {
								// convert to approx
								GriddedSurfaceImpl origSurf = (GriddedSurfaceImpl)surf;
								FaultTrace upper = origSurf.getEvenlyDiscritizedUpperEdge();
								FaultTrace lower = origSurf.getEvenlyDiscritizedLowerEdge();
								surf = new ApproxEvenlyGriddedSurface(upper, lower, 1d);
								System.out.println("\tTranslated to an ApproxEvenlyGriddedSurface");
								System.out.println("\tDip is "+surf.getAveDip());
							}
						}
					}
					Preconditions.checkNotNull(surf, "Didn't find it");
				} else {
//					surf = new PointSurface(lat, lon, depth);
					surf = new PointSurface(new Location(lat, lon, 0d));
					((PointSurface)surf).setAveDip(90d);
				}
				ProbEqkRupture rup = new ProbEqkRupture(mag, rake, probEach, surf, loc);
				keptEvents.add(rup);
			}
		}
		
		int bundleSize = 1000;
		List<List<ProbEqkRupture>> eventBundles = new ArrayList<>();
		List<ProbEqkRupture> curBundle = null;
		List<List<Integer>> eventIDBundles = new ArrayList<>();
		List<Integer> curIDBundle = null;
		CSVFile<String> eventsCSV = new CSVFile<>(true);
		eventsCSV.addLine("Event ID", "Magnitude", "Annual Rate", "Annual Nucleation Rate", "tau");
		double sumRate = 0d;
		double sumNuclRate = 0d;
		for (int n=0; n<keptEvents.size(); n++) {
			ProbEqkRupture rup = keptEvents.get(n);
			gmm.setEqkRupture(rup);
			ScalarGroundMotion gm = gmm.getGroundMotion();
			double rate = rup.getMeanAnnualRate(1d);
			double nuclRate = rate * keptFractsInReg.get(n);
			sumRate += rate;
			sumNuclRate += nuclRate;
			eventsCSV.addLine(n+"", (float)rup.getMag()+"", rate+"", nuclRate+"", (float)gm.tau()+"");
			
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
		System.out.println("Total rate:\t"+(float)sumRate);
		System.out.println("Nucleation rate:\t"+(float)sumNuclRate+" ("+pDF.format(sumNuclRate/sumRate)+")");
//		System.exit(0);
		
		// write out epicenter curves
		DiscretizedFunc[] epiMedianExceeds = new DiscretizedFunc[epiMinMags.length];
		DiscretizedFunc[] epiExceeds = new DiscretizedFunc[epiMinMags.length];
		for (int m=0; m<epiMinMags.length; m++) {
			epiMedianExceeds[m] = xVals.deepClone();
			epiMedianExceeds[m].scale(0d);
			epiExceeds[m] = xVals.deepClone();
			epiExceeds[m].scale(0d);
		}			
		
		ExecutorService exec = Executors.newFixedThreadPool(calcThreads);
		ArrayDeque<Future<EpiCalcCallable>> epiFutures = new ArrayDeque<>(calcThreads);
		
		int numWithHypos = 0;
		double epiSumRate = 0d;
		Map<RuptureSurface, BitSet> subSectInRegionBitSets = new HashMap<>();
		Map<RuptureSurface, Integer> subSectLocsInRegion = new HashMap<>();
		System.out.println("Doing epicentral IM calculations");
		for (ProbEqkRupture rup : keptEvents) {
			Location hypo = rup.getHypocenterLocation();
			double rate = rup.getMeanAnnualRate(1d);
			if (rate == 0d)
				continue;
			List<Site> epiSites = new ArrayList<>();
			double ratePerEpi;
			if (hypo == null) {
				RuptureSurface fullSurf = rup.getRuptureSurface();
				BitSet hypoBitSet = new BitSet(sitesAtNodes.size());
				List<? extends RuptureSurface> surfsToCheck;
				if (fullSurf instanceof CompoundSurface) {
					surfsToCheck = ((CompoundSurface)fullSurf).getSurfaceList();
				} else {
					surfsToCheck = List.of(fullSurf);
				}
				int totSurfLocs = 0;
				int mappedSurfLocs = 0;
				for (RuptureSurface surf : surfsToCheck) {
					BitSet subSurfBitSet = subSectInRegionBitSets.get(surf);
					if (subSurfBitSet != null) {
						// already cached
						hypoBitSet.or(subSurfBitSet);
						mappedSurfLocs += subSectLocsInRegion.get(surf);
						totSurfLocs += surf.getEvenlyDiscretizedNumLocs();
					} else {
						subSurfBitSet = new BitSet(sitesAtNodes.size());
						int myNumLocs = 0;
						int myMappedLocs = 0;
						for (Location loc : surf.getEvenlyDiscritizedListOfLocsOnSurface()) {
							myNumLocs++;
							int gridIndex = gridReg.indexForLocation(loc);
							if (gridIndex >= 0 && sitesAtNodes.get(gridIndex) != null) {
								myMappedLocs++;
								subSurfBitSet.set(gridIndex);
								hypoBitSet.set(gridIndex);
							}
						}
						totSurfLocs += myNumLocs;
						mappedSurfLocs += myMappedLocs;
						if (surfsToCheck.size() > 1) {
							subSectInRegionBitSets.put(surf, subSurfBitSet);
							subSectLocsInRegion.put(surf, myMappedLocs);
						}
					}
				}
				Preconditions.checkState(totSurfLocs > 0);
				Preconditions.checkState(mappedSurfLocs <= totSurfLocs);
				if (hypoBitSet.cardinality() == 0) {
					// no overlap with region
					continue;
				} else {
					double nuclRate = rate * (double)mappedSurfLocs/(double)totSurfLocs;
					ratePerEpi = nuclRate / (double)hypoBitSet.cardinality();
					for (int i=0; i<sitesAtNodes.size(); i++)
						if (hypoBitSet.get(i))
							epiSites.add(sitesAtNodes.get(i));
				}
			} else {
				numWithHypos++;
				int gridIndex = gridReg.indexForLocation(hypo);
				if (gridIndex < 0)
					continue;
				Site site = sitesAtNodes.get(gridIndex);
				if (site == null)
					// offshore
					continue;
				site = (Site)site.clone();
				site.setLocation(new Location(hypo.lat, hypo.lon));
				epiSites.add(site);
				ratePerEpi = rate;
			}
			Preconditions.checkState((float)(ratePerEpi*epiSites.size()) <= (float)rate);
			Preconditions.checkState(ratePerEpi > 0d);
			epiSumRate += ratePerEpi*epiSites.size();
			
			EpiCalcCallable calc;
			if (epiFutures.size() == calcThreads) {
				// wait on a previous future
				try {
					calc = epiFutures.pop().get();
				} catch (InterruptedException | ExecutionException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
				for (int m=0; m<epiMinMags.length; m++) {
					if ((float)calc.event.getMag() >= (float)epiMinMags[m]) {
						for (int i=0; i<xVals.size(); i++) {
							epiMedianExceeds[m].set(i, epiMedianExceeds[m].getY(i) + calc.epiMedianExceeds.getY(i));
							epiExceeds[m].set(i, epiExceeds[m].getY(i) + calc.epiExceeds.getY(i));
						}
					}
				}
			} else {
				calc = new EpiCalcCallable(empGMMRef, period, xVals);
			}
			calc.setEvent(rup, epiSites, ratePerEpi);
			epiFutures.addLast(exec.submit(calc));
		}
		while (!epiFutures.isEmpty()) {
			EpiCalcCallable calc;
			try {
				calc = epiFutures.pop().get();
			} catch (InterruptedException | ExecutionException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
			for (int m=0; m<epiMinMags.length; m++) {
				if ((float)calc.event.getMag() >= (float)epiMinMags[m]) {
					for (int i=0; i<xVals.size(); i++) {
						epiMedianExceeds[m].set(i, epiMedianExceeds[m].getY(i) + calc.epiMedianExceeds.getY(i));
						epiExceeds[m].set(i, epiExceeds[m].getY(i) + calc.epiExceeds.getY(i));
					}
				}
			}
		}
		System.out.println(numWithHypos+"/"+keptEvents.size()+" events had hypos");
		System.out.println("Epi on land rate:\t"+(float)epiSumRate+" ("+pDF.format(epiSumRate/sumNuclRate)
				+" of nucl, "+pDF.format(epiSumRate/sumRate)+" of total)");
		Preconditions.checkState(numWithHypos > 0);
		for (int m=0; m<epiMinMags.length; m++) {
			CSVFile<String> epiFuncsCSV = new CSVFile<>(true);
			epiFuncsCSV.addLine(perLabel, "Epicentral Median Exceedance Rate (1/yr)", "Epicentral Full Exceedance Rate (1/yr)");
			for (int i=0; i<xVals.size(); i++)
				epiFuncsCSV.addLine((float)xVals.getX(i)+"", (float)epiMedianExceeds[m].getY(i)+"", (float)epiExceeds[m].getY(i)+"");
			String epiPrefix = perSuffix+"_epicentral_exceedances_m"+oDF.format(epiMinMags[m]);
			epiFuncsCSV.writeToFile(new File(outputDir, epiPrefix+".csv"));
			epiMedianExceeds[m].setName("Epicenter Median Exceedances");
			epiExceeds[m].setName("Epicenter Exceedances");
			PlotSpec epiPlot = new PlotSpec(List.of(epiMedianExceeds[m], epiExceeds[m]),
					List.of(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_blue),
							new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_orange)),
					" ", perLabel, "Annual Exceedance Rate (1/yr)");
			epiPlot.setLegendVisible(true);
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			gp.drawGraphPanel(epiPlot, false, true, new Range(0d, 1d), new Range(1e-3, 1e1));
			PlotUtils.writePlots(outputDir, epiPrefix, gp, 800, 750, true, true, false);
		}
		
		if (doHazardMapCalc) {
			Preconditions.checkNotNull(hazardMapERF, "ERF is null");
			ArrayDeque<Future<MapCalcCallable>> mapFutures = new ArrayDeque<>(calcThreads);
			
			System.out.println("Calculating hazard maps for "+basinSites.size()+" sites");
			
			List<DiscretizedFunc> siteCurves = new ArrayList<>();
			
			for (Site site : basinSites) {
				MapCalcCallable calc;
				if (mapFutures.size() == calcThreads) {
					// wait on a previous future
					try {
						calc = mapFutures.pop().get();
					} catch (InterruptedException | ExecutionException e) {
						throw ExceptionUtils.asRuntimeException(e);
					}
					siteCurves.add(calc.getCurve());
					System.out.print(".");
					if (siteCurves.size() % 100 == 0)
						System.out.println(" "+siteCurves.size());
				} else {
					calc = new MapCalcCallable(hazMapGMMRef, period, xVals, hazardMapERF, filters);
				}
				calc.setSite(site);
				mapFutures.addLast(exec.submit(calc));
			}
			while (!mapFutures.isEmpty()) {
				MapCalcCallable calc;
				try {
					calc = mapFutures.pop().get();
				} catch (InterruptedException | ExecutionException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
				siteCurves.add(calc.getCurve());
				System.out.print(".");
				if (siteCurves.size() % 100 == 0)
					System.out.println(" "+siteCurves.size());
			}
			System.out.println(" "+siteCurves.size());
			Preconditions.checkState(siteCurves.size() == basinSites.size());
			List<String> header = new ArrayList<>();
			header.addAll(List.of("Index", "Latitude", "Longitude", "Vs30 (m/s)", "Z1.0 (m)", "Z2.5 (km)"));
			for (Point2D pt : xVals)
				header.add((float)pt.getX()+"");
			CSVFile<String> mapCSV = new CSVFile<>(true);
			mapCSV.addLine(header);
			for (int s=0; s<basinSites.size(); s++) {
				Site site = basinSites.get(s);
				List<String> line = new ArrayList<>(header.size());
				line.add(s+"");
				line.add((float)site.getLocation().lat+"");
				line.add((float)site.getLocation().lon+"");
				line.add(site.getParameter(Double.class, Vs30_Param.NAME).getValue().floatValue()+"");
				Double z1 = site.getParameter(Double.class, DepthTo1pt0kmPerSecParam.NAME).getValue();
				Double z25 = site.getParameter(Double.class, DepthTo2pt5kmPerSecParam.NAME).getValue();
				if (z1 == null)
					line.add("NaN");
				else
					line.add(z1.floatValue()+"");
				if (z25 == null)
					line.add("NaN");
				else
					line.add(z25.floatValue()+"");
				DiscretizedFunc curve = siteCurves.get(s);
				for (Point2D pt : curve)
					line.add(pt.getY()+"");
				mapCSV.addLine(line);
			}
			mapCSV.writeToFile(new File(outputDir, "site_hazard_curves.csv"));
			System.out.println("DONE with site curves");
		}
		
		if (doEventMapCalc) {
			List<CalcCallable> calls = new ArrayList<>();
			for (int i=0; i<calcThreads; i++)
				calls.add(new CalcCallable(gmmRef, period, sites));
			
			CompletableFuture<Void> writeFuture = null;
			
			File outZipFile = new File(outputDir, "event_mu_maps_"+perSuffix+".zip");
			ParallelZipFileOutput outData =
//					new ArchiveOutput.AsynchronousZipFileOutput(outZipFile);
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
				
//				if (b == 10)
//					break;
				
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
		}
		
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
	
	private static class EpiCalcCallable implements Callable<EpiCalcCallable> {
		private ScalarIMR gmm;
		private double[] xValsArray;
		private DiscretizedFunc logExceedFunc;
		private Map<RuptureSurface, CustomCacheWrappedSurface> wrappedMap;
		
		// inputs for each batch
		private ProbEqkRupture event;
		private List<Site> sites;
		private double ratePerEpi;
		
		// outputs
		private DiscretizedFunc epiMedianExceeds;
		private DiscretizedFunc epiExceeds;
		
		public EpiCalcCallable(AttenRelRef gmmRef, double period, DiscretizedFunc xVals) {
			gmm = gmmRef.get();
			setPeriod(gmm, period);
			xValsArray = new double[xVals.size()];
			double[] logXValsArray = new double[xVals.size()];
			for (int i=0; i<xVals.size(); i++) {
				xValsArray[i] = xVals.getX(i);
				logXValsArray[i] = Math.log(xValsArray[i]);
			}
			logExceedFunc = new LightFixedXFunc(logXValsArray, new double[xValsArray.length]);
			this.wrappedMap = new HashMap<>();
		}
		
		public void setEvent(ProbEqkRupture event, List<Site> sites, double ratePerEpi) {
			epiMedianExceeds = null;
			epiExceeds = null;
			this.event = event;
			this.sites = sites;
			this.ratePerEpi = ratePerEpi;
		}

		@Override
		public EpiCalcCallable call() throws Exception {
			Preconditions.checkNotNull(event);
			Preconditions.checkNotNull(sites);
			
			RuptureSurface surf = event.getRuptureSurface();
			if (surf instanceof CompoundSurface) {
				RuptureSurface wrappedSurf = DistCachedERFWrapper.getWrappedSurface(wrappedMap, surf);
				event = new DistCacheWrapperRupture(event, wrappedSurf);
			}
			
			epiMedianExceeds = new LightFixedXFunc(xValsArray, new double[xValsArray.length]);
			epiExceeds = new LightFixedXFunc(xValsArray, new double[xValsArray.length]);
			
			gmm.setEqkRupture(event);
			for (Site site : sites) {
				gmm.setSite(site);
				double mean = gmm.getMean();
				gmm.getExceedProbabilities(logExceedFunc);
//				if (site.getLocation().getLatitude() == 37.27 && site.getLocation().lon == -121.67 && event.getMag() == 4.3) {
				if (site.getLocation().getLatitude() == 37.04 && site.getLocation().lon == -121.78 && event.getMag() == 4.7) {
//				if (false) {
					System.out.println("DEBUG!!!");
					System.out.println("\tGM:\tmu="+(float)mean+" ln(g) = "+(float)Math.exp(mean)+" (g), sigma="+(float)gmm.getStdDev());
					if (gmm instanceof NSHMP_GMM_Wrapper) {
						NSHMP_GMM_Wrapper wrappedGMM = (NSHMP_GMM_Wrapper)gmm;
						System.out.println("GM logic tree: "+wrappedGMM.getGroundMotionTree());
						System.out.println("GM input: "+wrappedGMM.getCurrentGmmInput());
					}
					System.out.println(gmm.getShortName()+" params:");
					System.out.println("ERF:");
					for (Parameter<?> param : gmm.getEqkRuptureParams())
						System.out.println("\t"+param.getName()+":\t"+param.getValue());
					System.out.println("Propogation:");
					for (Parameter<?> param : gmm.getPropagationEffectParams())
						System.out.println("\t"+param.getName()+":\t"+param.getValue());
					System.out.println("Other:");
					for (Parameter<?> param : gmm.getOtherParams())
						System.out.println("\t"+param.getName()+":\t"+param.getValue());
					System.out.println("Site:");
					for (Parameter<?> param : gmm.getSiteParams())
						System.out.println("\t"+param.getName()+":\t"+param.getValue());
					System.out.println("\tLocation:\t"+gmm.getSite().getLocation());
					System.out.println("Intensity measure:\t"+gmm.getIntensityMeasure().getName());
				}
				for (int i=0; i<logExceedFunc.size(); i++) {
					if (mean >= logExceedFunc.getX(i))
						epiMedianExceeds.set(i, epiMedianExceeds.getY(i) + ratePerEpi);
					epiExceeds.set(i, epiExceeds.getY(i) + ratePerEpi*logExceedFunc.getY(i));
				}
			}
			sites = null;
			return this;
		}
	}
	
	private static class MapCalcCallable implements Callable<MapCalcCallable> {
		private ScalarIMR gmm;
		private double[] xValsArray;
		private DiscretizedFunc logExceedFunc;
		private DistCachedERFWrapper erf;
		private HazardCurveCalculator calc;
		
		// inputs for each batch
		private Site site;
		
		// outputs
		private DiscretizedFunc curve;
		
		public MapCalcCallable(AttenRelRef gmmRef, double period, DiscretizedFunc xVals, AbstractERF erf, SourceFilterManager filters) {
			gmm = gmmRef.get();
			setPeriod(gmm, period);
			xValsArray = new double[xVals.size()];
			double[] logXValsArray = new double[xVals.size()];
			for (int i=0; i<xVals.size(); i++) {
				xValsArray[i] = xVals.getX(i);
				logXValsArray[i] = Math.log(xValsArray[i]);
			}
			logExceedFunc = new LightFixedXFunc(logXValsArray, new double[xValsArray.length]);
			calc = new HazardCurveCalculator(filters);
			this.erf = new DistCachedERFWrapper(erf);
		}
		
		public void setSite(Site site) {
			this.site = site;
		}

		@Override
		public synchronized MapCalcCallable call() throws Exception {
			Preconditions.checkState(curve == null);
			Preconditions.checkNotNull(site);
			
			calc.getHazardCurve(logExceedFunc, site, gmm, erf);
			curve = new LightFixedXFunc(xValsArray, new double[xValsArray.length]);
			for (int i=0; i<xValsArray.length; i++)
				curve.set(i, logExceedFunc.getY(i));
			
			site = null;
			return this;
		}
		
		public DiscretizedFunc getCurve() {
			DiscretizedFunc curve = this.curve;
			Preconditions.checkNotNull(curve);
			this.curve = null;
			return curve;
		}
	}
	
	private static GriddedGeoDataSet[] loadBasinDepths(File file) throws IOException {
		CSVFile<String> csv = CSVFile.readFile(file, true);
		double minLat = Double.POSITIVE_INFINITY;
		double maxLat = Double.NEGATIVE_INFINITY;
		double minLon = Double.POSITIVE_INFINITY;
		double maxLon = Double.NEGATIVE_INFINITY;
		List<Location> locs = new ArrayList<>(csv.getNumRows()-1);
		List<Double> z1s = new ArrayList<>(csv.getNumRows()-1);
		List<Double> z25s = new ArrayList<>(csv.getNumRows()-1);
		for (int row=1; row<csv.getNumRows(); row++) {
			Location loc = new Location(csv.getDouble(row, 1), csv.getDouble(row, 0));
			minLat = Math.min(minLat, loc.lat);
			maxLat = Math.max(maxLat, loc.lat);
			minLon = Math.min(minLon, loc.lon);
			maxLon = Math.max(maxLon, loc.lon);
			locs.add(loc);
			z1s.add(csv.getDouble(row, 2)*1000d); // km -> m
			z25s.add(csv.getDouble(row, 3)); // km -> m
		}
		System.out.println("Basin depth region: "+minLat+"/"+minLon+" to "+maxLat+"/"+maxLon);
		Region reg = new Region(new Location(minLat-0.0001, minLon-0.0001), new Location(maxLat+0.0001, maxLon+0.0001));
		GriddedRegion gridReg = new GriddedRegion(reg, 0.01, locs.get(0));
		System.out.println("File has "+locs.size()+" locs, grid reg would have "+gridReg.getNodeCount());
		GriddedGeoDataSet[] ret = {
				new GriddedGeoDataSet(gridReg),
				new GriddedGeoDataSet(gridReg)
		};
		// init to NaN for no data
		for (GriddedGeoDataSet xyz : ret)
			for (int i=0; i<gridReg.getNodeCount(); i++)
				xyz.set(i, Double.NaN);
		for (int l=0; l<locs.size(); l++) {
			int ind = gridReg.indexForLocation(locs.get(l));
			Preconditions.checkState(ind >= 0);
			ret[0].set(ind, z1s.get(l));
			ret[1].set(ind, z25s.get(l));
		}
		return ret;
	}

}
