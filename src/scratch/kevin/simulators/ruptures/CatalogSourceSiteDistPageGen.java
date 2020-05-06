package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import org.dom4j.DocumentException;
import org.opensha.commons.data.Site;
import org.opensha.commons.geo.Region;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.iden.LogicalOrRupIden;
import org.opensha.sha.simulators.iden.RegionIden;
import org.opensha.sha.simulators.utils.RSQSimSubSectEqkRupture;

import com.google.common.base.Preconditions;
import com.google.common.collect.BiMap;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.HashBiMap;
import com.google.common.collect.Table;
import com.google.common.primitives.Ints;

import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.simCompare.IMT;
import scratch.kevin.simCompare.RuptureComparison;
import scratch.kevin.simCompare.SimulationRotDProvider;
import scratch.kevin.simCompare.SourceSiteDistPageGen;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.RSQSimCatalog.Loader;

public class CatalogSourceSiteDistPageGen extends SourceSiteDistPageGen<RSQSimEvent> {

	public CatalogSourceSiteDistPageGen(SimulationRotDProvider<RSQSimEvent> simProv, List<Site> sites) {
		super(simProv, sites);
	}
	
	private static Map<AttenRelRef, LinkedList<ScalarIMR>> gmpesInstancesCache = new HashMap<>();
	
	protected static ScalarIMR checkOutGMPE(AttenRelRef gmpeRef) {
		synchronized (gmpesInstancesCache) {
			LinkedList<ScalarIMR> gmpes = gmpesInstancesCache.get(gmpeRef);
			if (gmpes == null) {
				gmpes = new LinkedList<>();
				gmpesInstancesCache.put(gmpeRef, gmpes);
			}
			if (!gmpes.isEmpty())
				return gmpes.pop();
		}
		ScalarIMR gmpe = gmpeRef.instance(null);
		gmpe.setParamDefaults();
		gmpe.setIntensityMeasure(SA_Param.NAME);
		return gmpe;
	}
	
	protected static void checkInGMPE(AttenRelRef gmpeRef, ScalarIMR gmpe) {
		synchronized (gmpesInstancesCache) {
			gmpesInstancesCache.get(gmpeRef).push(gmpe);
		}
	}
	
	private static class EventComparison extends RuptureComparison.Cached<RSQSimEvent> {
		
		private Collection<Site> applicableSites;
		private double catDurationYears;
		private RSQSimSubSectEqkRupture gmpeRup;
		
		public EventComparison(RSQSimEvent event, RSQSimSubSectEqkRupture gmpeRup, double catDurationYears) {
			super(event);
			this.gmpeRup = gmpeRup;
			this.catDurationYears = catDurationYears;
			
			applicableSites = new HashSet<>();
		}
		
		public void addApplicableSite(Site site) {
			applicableSites.add(site);
		}
		
		public boolean isSiteApplicable(Site site) {
			return applicableSites.contains(site);
		}

		@Override
		public Collection<Site> getApplicableSites() {
			return applicableSites;
		}

		@Override
		public EqkRupture getGMPERupture() {
			return gmpeRup;
		}

		@Override
		public double getMagnitude() {
			return getRupture().getMagnitude();
		}

		@Override
		public double getAnnualRate() {
			return 1d/catDurationYears;
		}

		@Override
		public double getRuptureTimeYears() {
			return getRupture().getTimeInYears();
		}
	}
	
	private static class GMPECalcTask implements Runnable {
		

		private RSQSimCatalog catalog;
		private AttenRelRef[] gmpeRefs;
		private IMT[] imts;
		private List<Site> gmpeSites;
		private double catDurationYears;
		private Table<AttenRelRef, String, List<RuptureComparison<RSQSimEvent>>> sourceComps;
		private Map<Integer, String> parentIDtoSourceNameMap;
		private RSQSimEvent event;

		public GMPECalcTask(RSQSimCatalog catalog, AttenRelRef[] gmpeRefs, IMT[] imts, List<Site> gmpeSites,
				double catDurationYears, Table<AttenRelRef, String, List<RuptureComparison<RSQSimEvent>>> sourceComps,
				Map<Integer, String> parentIDtoSourceNameMap, RSQSimEvent event) {
			this.catalog = catalog;
			this.gmpeRefs = gmpeRefs;
			this.imts = imts;
			this.gmpeSites = gmpeSites;
			this.catDurationYears = catDurationYears;
			this.sourceComps = sourceComps;
			this.parentIDtoSourceNameMap = parentIDtoSourceNameMap;
			this.event = event;
		}

		@Override
		public void run() {
			RSQSimSubSectEqkRupture gmpeRup = catalog.getMappedSubSectRupture(event);
			HashSet<String> rupSources = new HashSet<>();
			for (FaultSectionPrefData sect : gmpeRup.getSubSections())
				if (parentIDtoSourceNameMap.containsKey(sect.getParentSectionId()))
					rupSources.add(parentIDtoSourceNameMap.get(sect.getParentSectionId()));
			for (String sourceName : rupSources) {
//				GMPECalcTask task = new GMPECalcTask(comp, gmpeSites, gmpeRefs, periods);
//				gmpeFutures.add(exec.submit(task));
				for (AttenRelRef gmpeRef : gmpeRefs) {
					EventComparison comp = new EventComparison(event, gmpeRup, catDurationYears);
					ScalarIMR gmpe = checkOutGMPE(gmpeRef);
					comp.calculate(gmpe, gmpeSites, imts);
					checkInGMPE(gmpeRef, gmpe);
					synchronized (sourceComps) {
						sourceComps.get(gmpeRef, sourceName).add(comp);
					}
				}
			}
		}
		
	}
	
	public static void main(String[] args) throws ZipException, IOException, DocumentException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		File outputDir = new File("/home/kevin/git/rsqsim-analysis/catalogs");
		File bbpParallelDir = new File("/home/kevin/bbp/parallel");
		
//		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance(baseDir);
////		File bbpDir = new File(bbpParallelDir, "2018_04_13-rundir2585_1myrs-all-m6.5-skipYears5000-noHF-csLASites");
//		File bbpDir = new File(bbpParallelDir, "2019_11_11-rundir2585_1myrs-all-m6.5-skipYears5000-noHF-vmLA_BASIN_500-cs500Sites");
		
//		RSQSimCatalog catalog = Catalogs.BRUCE_2740.instance(baseDir);
//		File bbpDir = new File(bbpParallelDir, "2018_09_10-rundir2740-all-m6.5-skipYears5000-noHF-csLASites");
		
//		RSQSimCatalog catalog = Catalogs.BRUCE_4860_10X.instance(baseDir);
//		File bbpDir = new File(bbpParallelDir, "2020_02_12-rundir4860_multi_combine-all-m6.5-skipYears5000-noHF-vmLA_BASIN_500-cs500Sites");
		
		RSQSimCatalog catalog = Catalogs.BRUCE_4983.instance(baseDir);
		File bbpDir = new File(bbpParallelDir, "2020_04_19-rundir4983-all-m6.5-skipYears5000-noHF-vmLA_BASIN_500-cs500Sites");
		
		VelocityModel vm = RSQSimBBP_Config.detectVM(bbpDir);
		
		List<String> sourceNames = new ArrayList<>();
		List<int[]> parentIDs = new ArrayList<>();
		
		sourceNames.add("San Andreas (Mojave)");
		parentIDs.add(new int[] { 286, 301});
		
		sourceNames.add("Puente Hills");
		parentIDs.add(new int[] { 240});
		
		String[] siteNames = null; // all
		
		double minMag = 6.5d;
		int skipYears = 5000;
		boolean hypoSort = true;
		
		AttenRelRef[] gmpeRefs = { AttenRelRef.ASK_2014, AttenRelRef.BSSA_2014, AttenRelRef.CB_2014, AttenRelRef.CY_2014 };
		IMT[] imts = { IMT.SA3P0, IMT.SA5P0, IMT.SA10P0 };
		
		ZipFile bbpZip = new ZipFile(new File(bbpDir, "results_rotD.zip"));
		
		List<BBP_Site> sites = BBP_Site.readFile(bbpDir);
		
		ArrayList<RegionIden> siteRegIdens = new ArrayList<>();
		HashBiMap<BBP_Site, Site> sitesBBPtoGMPE = HashBiMap.create();
		List<Site> gmpeSites = new ArrayList<>();
		for (BBP_Site site : sites) {
			if (siteNames != null) {
				boolean found = false;
				for (String name : siteNames) {
					if (name.startsWith(site.getName())) {
						found = true;
						break;
					}
				}
				if (!found)
					continue;
			}
			siteRegIdens.add(new RegionIden(new Region(site.getLoc(), MPJ_BBP_CatalogSim.CUTOFF_DIST)));
			Site gmpeSite = site.buildGMPE_Site(vm);
			gmpeSite.setName(RSQSimBBP_Config.siteCleanName(site));
			gmpeSites.add(gmpeSite);
			sitesBBPtoGMPE.put(site, gmpeSite);
		}
		
		Loader loader = catalog.loader();
		if (minMag > 0d)
			loader.minMag(minMag);
		if (skipYears > 0)
			loader.skipYears(skipYears);
		loader.matches(new LogicalOrRupIden(siteRegIdens));
		HashSet<Integer> allParents = new HashSet<>();
		for (int[] parents : parentIDs)
			for (int parent : parents)
				allParents.add(parent);
		loader.forParentSections(true, Ints.toArray(allParents));
		System.out.println("Loading events...");
		List<RSQSimEvent> events = loader.minMag(minMag).load();
		System.out.println("Loaded "+events.size()+" events");
		Map<Integer, RSQSimEvent> eventsMap = new HashMap<>();
		for (RSQSimEvent event : events)
			eventsMap.put(event.getID(), event);
		
		double catDurationYears = catalog.getDurationYears();
		
		BBP_CatalogSimZipLoader bbpZipFile = new BBP_CatalogSimZipLoader(bbpZip, sites, sitesBBPtoGMPE, eventsMap);
		
		Table<AttenRelRef, String, List<RuptureComparison<RSQSimEvent>>> sourceComps = HashBasedTable.create();
		
		Map<Integer, String> parentIDtoSourceNameMap = new HashMap<>();
		for (int i=0; i<parentIDs.size(); i++) {
			String sourceName = sourceNames.get(i);
			for (int parentID : parentIDs.get(i))
				parentIDtoSourceNameMap.put(parentID, sourceName);
			for (AttenRelRef gmpeRef : gmpeRefs) {
				List<RuptureComparison<RSQSimEvent>> comps = new ArrayList<>();
				sourceComps.put(gmpeRef, sourceName, comps);
			}
		}
		
		ExecutorService exec = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
		
		List<Future<?>> gmpeFutures = new ArrayList<>();
		System.out.println("Processing events/building ruptures");
		for (RSQSimEvent event : events)
			gmpeFutures.add(exec.submit(new GMPECalcTask(catalog, gmpeRefs, imts, gmpeSites, catDurationYears,
					sourceComps, parentIDtoSourceNameMap, event)));
		System.out.println("Waiting on "+gmpeFutures.size()+" GMPE futures");
		for (Future<?> future : gmpeFutures) {
			try {
				future.get();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		exec.shutdown();
		
		File catalogOutputDir = new File(outputDir, catalog.getCatalogDir().getName());
		Preconditions.checkState(catalogOutputDir.exists() || catalogOutputDir.mkdir());
		
		File vmDir = new File(catalogOutputDir, "bbp_"+vm.name());
		Preconditions.checkState(vmDir.exists() || vmDir.mkdir());
		
		File sourceOutputDir = new File(vmDir, "source_site_comparisons");
		Preconditions.checkState(sourceOutputDir.exists() || sourceOutputDir.mkdir());
		
		CatalogSourceSiteDistPageGen pageGen = new CatalogSourceSiteDistPageGen(bbpZipFile, gmpeSites);
		
		List<String> headerLines = new ArrayList<>();
		headerLines.add("# "+catalog.getName()+" Source/Site GMPE Comparisons");
		headerLines.add("");
		headerLines.add("**GMPEs:**");
		for (AttenRelRef gmpe : gmpeRefs)
			headerLines.add("* "+gmpe.getName());
		
		pageGen.generatePage(sourceComps, sourceOutputDir, headerLines, imts, hypoSort);
		
		catalog.writeMarkdownSummary(catalogOutputDir, true, false);
		RSQSimCatalog.writeCatalogsIndex(outputDir);
	}

}
