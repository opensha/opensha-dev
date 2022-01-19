package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipFile;

import org.dom4j.DocumentException;
import org.opensha.commons.data.Site;
import org.opensha.commons.geo.Region;
import org.opensha.commons.util.FileNameComparator;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.iden.LogicalOrRupIden;
import org.opensha.sha.simulators.iden.RegionIden;

import com.google.common.base.Preconditions;
import com.google.common.collect.BiMap;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.HashBiMap;
import com.google.common.collect.Table;

import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.simCompare.IMT;
import scratch.kevin.simCompare.SimulationRotDProvider;
import scratch.kevin.simCompare.SiteHazardCurveComarePageGen;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.RSQSimCatalog.Loader;
import scratch.kevin.simulators.ruptures.CatalogGMPE_Compare.EventComparison;

public class CatalogSiteHazardCurveComparePageGen extends SiteHazardCurveComarePageGen<RSQSimEvent> {

	public CatalogSiteHazardCurveComparePageGen(SimulationRotDProvider<RSQSimEvent> simProv, String simName) {
		super(simProv, simName);
	}
	
	public static void main(String[] args) throws IOException, DocumentException {
		File outputDir = new File("/home/kevin/markdown/rsqsim-analysis/catalogs");
		File bbpParallelDir = new File("/home/kevin/bbp/parallel");
		
		RSQSimCatalog catalog = Catalogs.BRUCE_5310.instance();
//		RSQSimCatalog catalog = Catalogs.BRUCE_4983_STITCHED.instance();
		
		String[] siteNames = { "USC", "OSI", "SBSM" };
//		String[] siteNames = { "SBSM" };
//		String[] siteNames = { "USC", "SMCA", "OSI", "WSS", "SBSM",
//				"LAF", "s022", "STNI", "WNGC", "PDE" };

//		VelocityModel forceVM = VelocityModel.LA_BASIN_500;
//		VelocityModel forceVM = VelocityModel.LA_BASIN_863;
		VelocityModel forceVM = null;
		
		boolean sourceFractional = true;
		
		boolean replotCurves = true;
		boolean replotDisaggs = false;
		
//		AttenRelRef[] gmpeRefs = { AttenRelRef.NGAWest_2014_AVG_NOIDRISS, AttenRelRef.ASK_2014 };
//		AttenRelRef[] gmpeRefs = { AttenRelRef.NGAWest_2014_AVG_NOIDRISS };
		AttenRelRef[] gmpeRefs = { AttenRelRef.ASK_2014 };
		IMT[] imts = IMT.forPeriods(new double[] { 3, 5, 7.5, 10 });
		
		// find BBP parallel dir
		String catalogDirName = catalog.getCatalogDir().getName();
		if (catalogDirName.startsWith("JG_"))
			// I sometimes modify Jacqui's directory names
			catalogDirName = catalogDirName.substring(3);
		File bbpDir = null;
		File bbpZipFile = null;
		File[] allBBPDirs = bbpParallelDir.listFiles();
		Arrays.sort(allBBPDirs, new FileNameComparator());
		VelocityModel bbpVM = null;
		for (File dir : allBBPDirs) {
			String name = dir.getName();
			if (dir.isDirectory() && name.contains(catalogDirName) && name.contains("-all")) {
				if (name.contains("-rg"))
					continue;
				if (name.contains("-gridded"))
					continue;
				if (name.contains("-timeScale"))
					continue;
				File zipFile = new File(dir, "results.zip");
				if (!zipFile.exists())
					zipFile = new File(dir, "results_rotD.zip");
				if (zipFile.exists()) {
					System.out.println("Found candidate zip file: "+zipFile.getAbsolutePath());
					VelocityModel myVM = null;
					List<BBP_Site> sites = BBP_Site.readFile(dir);
					for (VelocityModel v : VelocityModel.values()) {
						if (name.contains(v.name())) {
							myVM = v;
							System.out.println("\tdectected VM: "+myVM);
							break;
						}
						if ((float)sites.get(0).getVs30() == v.getVs30()) {
							myVM = v;
							System.out.println("\tassuming VM from Vs30: "+myVM);
						}
					}
					if (forceVM != null && myVM != forceVM) {
						System.out.println("Skipping dir, wrong VM");
						continue;
					}
					
					bbpDir = dir;
					bbpZipFile = zipFile;
					bbpVM = myVM;
				}
			}
		}
		Preconditions.checkNotNull(bbpDir);
		System.out.println("Located ref BBP dir: "+bbpDir.getAbsolutePath());
		System.out.println("Velocity Model: "+bbpVM);
		
		File gmpeCacheDir = new File(bbpDir, "gmpe_cache");
		Preconditions.checkState(gmpeCacheDir.exists() || gmpeCacheDir.mkdir());
		
		File catalogOutputDir = new File(outputDir, catalog.getCatalogDir().getName());
		Preconditions.checkState(catalogOutputDir.exists() || catalogOutputDir.mkdir());
		
		File vmDir = new File(catalogOutputDir, "bbp_"+bbpVM.name());
		Preconditions.checkState(vmDir.exists() || vmDir.mkdir());
		
		String bbpDirName = bbpDir.getName();
		double minMag;
		if (bbpDirName.contains("-all-m")) {
			String magStr = bbpDirName.substring(bbpDirName.indexOf("-all-m")+"-all-m".length());
			if (magStr.contains("-"))
				magStr = magStr.substring(0, magStr.indexOf("-"));
			minMag = Double.parseDouble(magStr);
			System.out.println("Detected minMag="+minMag);
		} else {
			throw new IllegalStateException("Couldn't detect minMag from "+bbpDirName);
		}
		int skipYears = 0;
		if (bbpDirName.contains("-skipYears")) {
			String yearStr = bbpDirName.substring(bbpDirName.indexOf("-skipYears")+"-skipYears".length());
			if (yearStr.contains("-"))
				yearStr = yearStr.substring(0, yearStr.indexOf("-"));
			skipYears = Integer.parseInt(yearStr);
			System.out.println("Detected skipYears="+skipYears);
		}
		
		List<BBP_Site> sites = BBP_Site.readFile(bbpDir);
		
		System.out.println("Zip file: "+bbpZipFile.getAbsolutePath());
		ZipFile zipFile = new ZipFile(bbpZipFile);
		
		CatalogGMPE_Compare gmpeComp = new CatalogGMPE_Compare(catalog, zipFile, sites, minMag, skipYears,
				gmpeCacheDir, null, bbpVM);
		
		Table<String, RSQSimEvent, Double> sourceContribFracts =
				getSourceContribFracts(catalog, gmpeComp.getEvents(), sourceFractional);
		
		CatalogSiteHazardCurveComparePageGen pageGen = new CatalogSiteHazardCurveComparePageGen(
				gmpeComp.getSimProv(), catalog.getName());
		pageGen.setReplotCurves(replotCurves);
		pageGen.setReplotDisaggs(replotDisaggs);
//		pageGen.setSourceRupContributionFractions(sourceContribFracts, 4e-4, 10);
		pageGen.setSourceRupContributionFractions(sourceContribFracts, 0d, 10); // 0 = RTGM
		
		for (AttenRelRef gmpeRef : gmpeRefs) {
			List<EventComparison> comps = gmpeComp.loadCalcComps(gmpeRef, imts);
			for (String siteName : siteNames) {
				Site site = null;
				for (Site oSite : gmpeComp.getGMPESites())
					if (oSite.getName().equals(siteName))
						site = oSite;
				Preconditions.checkNotNull(site, "Site %s not found", siteName);
				List<EventComparison> myComps = new ArrayList<>();
				for (EventComparison comp : comps)
					if (comp.hasSite(site))
						myComps.add(comp);
				Preconditions.checkState(!myComps.isEmpty());
				List<String> headerLines = null;
				File siteOutputDir = new File(vmDir, "site_hazard_"+siteName+"_"+gmpeRef.getShortName());
				Preconditions.checkState(siteOutputDir.exists() || siteOutputDir.mkdir());
				pageGen.generateSitePage(site, myComps, siteOutputDir, headerLines, imts, gmpeRef);
			}
		}
		
		catalog.writeMarkdownSummary(catalogOutputDir, true, false);
		RSQSimCatalog.writeCatalogsIndex(outputDir);
		
		gmpeComp.shutdown();
		
		System.exit(0);
	}

	private static Table<String, RSQSimEvent, Double> getSourceContribFracts(RSQSimCatalog catalog,
			List<RSQSimEvent> events, boolean sourceFractional) {
		Map<String, List<Integer>> faultNamesToIDsMap = catalog.getFaultModel().getNamedFaultsMapAlt();
		Map<Integer, String> idsToFaultNamesMap = new HashMap<>();
		for (String faultName : faultNamesToIDsMap.keySet()) {
			String name = faultName;
			if (name.startsWith("San Andreas"))
				name = "San Andreas";
			else if (name.startsWith("San Jacinto"))
				name = "San Jacinto";
			for (Integer id : faultNamesToIDsMap.get(faultName))
				idsToFaultNamesMap.put(id, name);
		}
		
		Table<String, RSQSimEvent, Double> table = HashBasedTable.create();
		for (RSQSimEvent event : events) {
			List<? extends FaultSection> sects = catalog.getMappedSubSectRupture(event).getSubSections();
			double totArea = 0d;
			List<Double> areas = new ArrayList<>();
			for (FaultSection sect : sects) {
				double area = sect.getOrigDownDipWidth()*sect.getTraceLength();
				totArea += area;
				areas.add(area);
			}
			Map<String, Double> sourceFracts = new HashMap<>();
			for (int i=0; i<sects.size(); i++) {
				FaultSection sect = sects.get(i);
				int id = sect.getParentSectionId();
				String name = idsToFaultNamesMap.get(id);
				if (name == null) {
					// not a named fault
					name = sect.getParentSectionName();
					if (name.startsWith("San Andreas"))
						name = "San Andreas";
					else if (name.startsWith("San Jacinto"))
						name = "San Jacinto";
				}
				if (sourceFractional) {
					double fract = areas.get(i)/totArea;
					if (sourceFracts.containsKey(name))
						fract += sourceFracts.get(name);
					sourceFracts.put(name, fract);
				} else {
					sourceFracts.put(name, 1d);
				}
			}
			for (String name : sourceFracts.keySet())
				table.put(name, event, sourceFracts.get(name));
		}
		
		return table;
	}

}
