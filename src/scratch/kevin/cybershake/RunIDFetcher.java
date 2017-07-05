package scratch.kevin.cybershake;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.siteData.SiteData;
import org.opensha.commons.data.siteData.impl.CVM4i26BasinDepth;
import org.opensha.commons.data.siteData.impl.CVM_Vs30;
import org.opensha.commons.data.siteData.impl.WillsMap2006;
import org.opensha.commons.geo.Location;
import org.opensha.commons.data.siteData.impl.CVM_Vs30.CVM;
import org.opensha.sha.cybershake.HazardCurveFetcher;
import org.opensha.sha.cybershake.db.CybershakeSite;
import org.opensha.sha.cybershake.db.Cybershake_OpenSHA_DBApplication;
import org.opensha.sha.cybershake.db.DBAccess;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class RunIDFetcher {
	
	private HazardCurveFetcher fetch;
	private List<CybershakeSite> allSites;
	private Map<String, CybershakeSite> siteNameToSiteMap;
	private List<Integer> allRuns;
	private Map<String, Integer> nameToRunMap;
	
	public RunIDFetcher(DBAccess db, int datasetID, int imTypeID) {
		fetch = new HazardCurveFetcher(db, datasetID, imTypeID);
		allSites = fetch.getCurveSites();
		siteNameToSiteMap = Maps.newHashMap();
		allRuns = fetch.getRunIDs();
		nameToRunMap = Maps.newHashMap();
		
		Preconditions.checkState(allSites.size() == allRuns.size());
		for (int i=0; i<allSites.size(); i++) {
			String site = allSites.get(i).short_name;
			int runID = allRuns.get(i);
			Preconditions.checkState(!nameToRunMap.containsKey(site), "duplicate found for site %s", site);
			nameToRunMap.put(site, runID);
			siteNameToSiteMap.put(site, allSites.get(i));
		}
	}
	
	public Integer getRunID(String siteName) {
		return nameToRunMap.get(siteName);
	}
	
	public CybershakeSite getSite(String siteName) {
		return siteNameToSiteMap.get(siteName);
	}
	
	public static void main(String[] args) throws IOException {
//		int datasetID = 61;
		int datasetID = 57;
//		int imTypeID = 21;
		int imTypeID = 146; // 3sec SA, RotD100
		
		boolean writeSiteData = true;
		
		List<List<String>> siteNamesLists = Lists.newArrayList();
		List<String> groupNames = Lists.newArrayList();
		
		groupNames.add("UGMS Orig 14");
		siteNamesLists.add(Lists.newArrayList("CCP", "COO", "LADT", "LAPD", "P22", "PAS", "s429",
				"s603", "s684", "s758", "SBSM", "SMCA", "STNI", "WNGC"));
		
		groupNames.add("LA Basin West Array");
		siteNamesLists.add(Lists.newArrayList("s344", "s345", "s346", "s347", "s348", "s349", "s351", "s353"));
		
		groupNames.add("LA Basin Central Array");
		siteNamesLists.add(Lists.newArrayList("s383", "s385", "s387", "s388", "s389", "s391", "s393", "s395"));
		
		groupNames.add("LA Basin East Array");
		siteNamesLists.add(Lists.newArrayList("s470", "s472", "s474", "s476", "s478", "s480", "s518", "s520"));
		
		groupNames.add("San Fernando Valley");
		siteNamesLists.add(Lists.newArrayList("s311", "s275", "P3", "s238", "P2", "P21"));
		
		groupNames.add("Orange County");
		siteNamesLists.add(Lists.newArrayList("s510", "s512", "s514", "s516", "s550", "s552", "s554", "s556", "s591", "s593", "s595"));
		
		groupNames.add("Inland Empire");
		siteNamesLists.add(Lists.newArrayList("s599", "s601", "s603", "s605", "s644", "s646", "s647", "s688", "s689"));
		
		DBAccess db = Cybershake_OpenSHA_DBApplication.getDB(Cybershake_OpenSHA_DBApplication.ARCHIVE_HOST_NAME);
		RunIDFetcher fetch = new RunIDFetcher(db, datasetID, imTypeID);
		
		Joiner j = Joiner.on(",");
		
		List<String> missingSites = Lists.newArrayList();
		
		List<List<Integer>> runIDsList = Lists.newArrayList();
		
		for (int i=0; i<groupNames.size(); i++) {
			String name = groupNames.get(i);
			System.out.println(name);
			
			List<String> sites = siteNamesLists.get(i);
			List<Integer> runIDs = Lists.newArrayList();
			for (String site : sites) {
				Integer runID = fetch.getRunID(site);
				if (runID == null) {
					missingSites.add(site);
					runID = -1;
				}
//				Preconditions.checkNotNull(runID,
//						"No run ID found for site %s, hazard dataset id=%s, im type id=%s", site, datasetID, imTypeID);
				runIDs.add(runID);
			}
			
			System.out.println("\t"+j.join(sites));
			System.out.println("\t"+j.join(runIDs));
			
			runIDsList.add(runIDs);
		}
		
		if (!missingSites.isEmpty()) {
			System.out.println("WARNING: The following sites were missing for hazard dataset id="+datasetID+", im type id="+imTypeID);
			System.out.println("\t"+j.join(missingSites));
		}
		
		if (writeSiteData) {
			CSVFile<String> csv = new CSVFile<String>(false);
			List<SiteData<Double>> siteDataProvs = Lists.newArrayList();
			List<String> header = Lists.newArrayList("CS Site Name", "CS Run ID", "Latitude", "Longitude");
			siteDataProvs.add(new WillsMap2006());
			header.add("Wills 2006 Vs30 (m/s)");
			siteDataProvs.add(new CVM_Vs30(CVM.CVMS4i26));
			header.add("CVMS-4.26 Vs30 (m/s)");
			siteDataProvs.add(new CVM4i26BasinDepth(SiteData.TYPE_DEPTH_TO_1_0));
			header.add("CVMS-4.26 Z1.0 (km)");
			siteDataProvs.add(new CVM4i26BasinDepth(SiteData.TYPE_DEPTH_TO_2_5));
			header.add("CVMS-4.26 Z2.5 (km)");
			csv.addLine(header);
			
			for (int i=0; i<groupNames.size(); i++) {
				csv.addLine(groupNames.get(i));
				List<String> sites = siteNamesLists.get(i);
				List<Integer> runs = runIDsList.get(i);
				for (int s=0; s<sites.size(); s++) {
					CybershakeSite site = fetch.getSite(sites.get(s));
					List<String> line = Lists.newArrayList(site.short_name, runs.get(s)+"", site.lat+"", site.lon+"");
					Location loc = site.createLocation();
					for (SiteData<Double> prov : siteDataProvs)
						line.add(prov.getValue(loc)+"");
					csv.addLine(line);
				}
			}
			csv.writeToFile(new File("/tmp/ugms_sites.csv"));
		}
		
		db.destroy();
	}

}
