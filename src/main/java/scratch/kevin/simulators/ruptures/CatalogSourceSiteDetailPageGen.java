package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipFile;

import org.dom4j.DocumentException;
import org.opensha.commons.data.Site;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.iden.LogicalOrRupIden;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBiMap;
import com.google.common.primitives.Ints;

import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.simCompare.IMT;
import scratch.kevin.simCompare.SimulationRotDProvider;
import scratch.kevin.simCompare.SourceSiteDetailPageGen;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.RSQSimCatalog.Loader;

public class CatalogSourceSiteDetailPageGen extends SourceSiteDetailPageGen {

	public CatalogSourceSiteDetailPageGen(SimulationRotDProvider<RSQSimEvent> simProv, String sourceName,
			int[] parentSectIDs, RSQSimCatalog catalog, List<RSQSimEvent> events, List<Site> sites) throws IOException {
		super(simProv, sourceName, parentSectIDs, catalog, events, sites);
	}

	public static void main(String[] args) throws IOException, DocumentException {
		File outputDir = new File("/home/kevin/markdown/rsqsim-analysis/catalogs");
		File bbpParallelDir = new File("/home/kevin/bbp/parallel");

		RSQSimCatalog catalog = Catalogs.BRUCE_3062.instance();
//		File bbpDir = new File(bbpParallelDir, "2019_02_09-rundir3062-all-m6.5-skipYears5000-noHF-csLASites");
		File bbpDir = new File(bbpParallelDir, "2020_05_26-rundir3062-all-m6.5-skipYears5000-noHF-vmLA_BASIN_500-cs500Sites");
		
//		RSQSimCatalog catalog = Catalogs.BRUCE_4983_STITCHED.instance();
//		File bbpDir = new File(bbpParallelDir, "2020_05_05-rundir4983_stitched-all-m6.5-skipYears5000-noHF-vmLA_BASIN_500-cs500Sites");
		
		VelocityModel vm = RSQSimBBP_Config.detectVM(bbpDir);
		
//		String sourceName = "San Jacinto (Northern)";
//		String sourcePrefix = "nsjc";
//		int[] parentIDs = { 119, 289, 401, 293, 292 };
		
		String sourceName = "San Jacinto";
		String sourcePrefix = "sjc";
		int[] parentIDs = { 119, 289, 401, 293, 292, 101, 99, 28 };
		
//		String sourceName = "San Andreas (Mojave)";
//		String sourcePrefix = "saf_mojave";
//		int[] parentIDs = { 286, 301};
		
		String[] siteNames = { "SBSM" }; // all
		
		double minMag = 6d;
		int skipYears = 5000;
		
		IMT[] imts = { IMT.SA3P0, IMT.SA5P0, IMT.SA10P0 };
		
		ZipFile bbpZip = new ZipFile(new File(bbpDir, "results_rotD.zip"));
		
		List<BBP_Site> bbpSites = BBP_Site.readFile(bbpDir);
		List<Site> sites = new ArrayList<>();
		HashBiMap<BBP_Site, Site> sitesBBPtoGMPE = HashBiMap.create();
		for (int i=bbpSites.size(); --i>=0;) {
			BBP_Site bbpSite = bbpSites.get(i);
			if (siteNames != null) {
				boolean found = false;
				for (String siteName : siteNames)
					found = found || bbpSite.getName().equals(siteName);
				if (!found) {
					bbpSites.remove(i);
					continue;
				}
			}
			Site site = bbpSite.buildGMPE_Site(vm);
			sites.add(site);
			sitesBBPtoGMPE.put(bbpSite, site);
		}
		
		Loader loader = catalog.loader();
		if (minMag > 0d)
			loader.minMag(minMag);
		if (skipYears > 0)
			loader.skipYears(skipYears);
		HashSet<Integer> allParents = new HashSet<>();
		for (int parent : parentIDs)
			allParents.add(parent);
		loader.forParentSections(true, Ints.toArray(allParents));
		System.out.println("Loading events...");
		List<RSQSimEvent> events = loader.minMag(minMag).load();
		System.out.println("Loaded "+events.size()+" events");
		Map<Integer, RSQSimEvent> eventsMap = new HashMap<>();
		for (RSQSimEvent event : events)
			eventsMap.put(event.getID(), event);
		
		BBP_CatalogSimZipLoader bbpZipFile = new BBP_CatalogSimZipLoader(bbpZip, bbpSites, sitesBBPtoGMPE, eventsMap);
		
		CatalogSourceSiteDetailPageGen pageGen = new CatalogSourceSiteDetailPageGen(
				bbpZipFile, sourceName, parentIDs, catalog, events, sites);
		
		List<String> metadataLines = null;
		
		File catalogOutputDir = new File(outputDir, catalog.getCatalogDir().getName());
		Preconditions.checkState(catalogOutputDir.exists() || catalogOutputDir.mkdir());
		
		File vmDir = new File(catalogOutputDir, "bbp_"+vm.name());
		Preconditions.checkState(vmDir.exists() || vmDir.mkdir());
		
		File sourceOutputDir = new File(vmDir, "source_site_"+sourcePrefix);
		Preconditions.checkState(sourceOutputDir.exists() || sourceOutputDir.mkdir());
		
		pageGen.generatePage(sourceOutputDir, metadataLines, imts);
		
		catalog.writeMarkdownSummary(catalogOutputDir, true, false);
		RSQSimCatalog.writeCatalogsIndex(outputDir);
	}

}
