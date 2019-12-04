package scratch.kevin.simulators.ruptures.rotation;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.dom4j.DocumentException;
import org.opensha.commons.data.Site;
import org.opensha.commons.util.FileNameComparator;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_WrapperFullParam;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_Wrappers;
import org.opensha.sha.simulators.RSQSimEvent;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Ints;

import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.simCompare.SimulationRotDProvider;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.ruptures.ASK_EventData;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;
import scratch.kevin.simulators.ruptures.RSQSimBBP_Config;
import scratch.kevin.simulators.ruptures.rotation.RotatedRupVariabilityConfig.RotationSpec;

public class RSQSimRotatedRupVariabilityScenarioPageGen extends RSQSimRotatedRupVariabilityPageGen {

	private Scenario scenario;

	public RSQSimRotatedRupVariabilityScenarioPageGen(RSQSimCatalog catalog, Scenario scenario,
			RSQSimRotatedRupVariabilityConfig config, SimulationRotDProvider<RotationSpec> prov, double[] calcPeriods) {
		super(catalog, config, scenario.getMagnitude(), prov, calcPeriods);
		this.scenario = scenario;
	}

	@Override
	protected String getScenarioName() {
		return scenario.getName();
	}

	@Override
	protected String getScenarioShortName() {
		return scenario.getShortName();
	}

	@Override
	protected String[] getScenarioMatchCriteria() {
		return scenario.getMatchCriteria();
	}

	@Override
	protected Scenario getBBP_PartB_Scenario(RotatedRupVariabilityConfig<RSQSimEvent> config) {
		return scenario;
	}

	@SuppressWarnings("unused")
	public static void main(String[] args) throws IOException, DocumentException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		File outputDir = new File("/home/kevin/git/rsqsim-analysis/catalogs");
		File bbpParallelDir = new File("/home/kevin/bbp/parallel");

//		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance(baseDir);
//		RSQSimCatalog catalog = Catalogs.BRUCE_2740.instance(baseDir);
		RSQSimCatalog catalog = Catalogs.BRUCE_4320.instance(baseDir);

//		File bbpDir = new File(bbpParallelDir,
//				"2019_11_22-rundir2585_1myrs-rotatedRups-m7p2_vert_ss_surface_rnd_mag_0p05"
//				+ "-3dists-18srcAz-1siteSrcAz-100rups-skipYears5000-vmLA_BASIN_500-noHF-1site");
//		File bbpDir = new File(bbpParallelDir,
//				"2019_11_21-rundir2585_1myrs-rotatedRups-2scenarios-3dists-18srcAz-1siteSrcAz"
//				+ "-100rups-skipYears5000-vmLA_BASIN_500-noHF-1site");
//		File bbpDir = new File(bbpParallelDir,
//				"2019_02_27-rundir2585_1myrs-rotatedRups-6scenarios-3dists-18srcAz-1siteSrcAz"
//				+ "-400rups-skipYears5000-vmLA_BASIN_500-noHF-1site");
		File bbpDir = null;
		
//		VelocityModel forceVM = VelocityModel.LA_BASIN_863;
		VelocityModel forceVM = VelocityModel.LA_BASIN_500;
//		VelocityModel forceVM = null;
		
		double timeScale = 1d;
		boolean scaleVelocities = false;
		
		double[] calcPeriods = {1d, 2d, 3d, 4d, 5d, 7.5, 10d};
		double[] periods = {3d, 5d, 7.5, 10d};
		
		Map<Integer, List<ASK_EventData>> realData = ASK_EventData.load(1d);
		
		NGAW2_WrapperFullParam[] refGMPEs = { new NGAW2_Wrappers.ASK_2014_Wrapper(), new NGAW2_Wrappers.BSSA_2014_Wrapper(),
				new NGAW2_Wrappers.CB_2014_Wrapper(), new NGAW2_Wrappers.CY_2014_Wrapper()};
//		NGAW2_WrapperFullParam[] refGMPEs = { new NGAW2_Wrappers.BSSA_2014_Wrapper() };
		
		System.out.println("Catalog: "+catalog.getName());
		// find BBP parallel dir
		String catalogDirName = catalog.getCatalogDir().getName();
		if (catalogDirName.startsWith("JG_"))
			// I sometimes modify Jacqui's directory names
			catalogDirName = catalogDirName.substring(3);
		File bbpZipFile = null;
		File[] allBBPDirs;
		if (bbpDir == null) {
			allBBPDirs = bbpParallelDir.listFiles();
			Arrays.sort(allBBPDirs, new FileNameComparator());
		} else {
			allBBPDirs = new File[] { bbpDir };
		}
		for (File dir : allBBPDirs) {
			String name = dir.getName();
			if (dir.isDirectory() && name.contains(catalogDirName) && name.contains("-rotatedRups-")) {
				File zipFile = new File(dir, "results.zip");
				if (!zipFile.exists())
					zipFile = new File(dir, "results_rotD.zip");
				if (zipFile.exists()) {
					if (forceVM != null && RSQSimBBP_Config.detectVM(dir) != forceVM)
						continue;
					if (timeScale == 1d && name.contains("-timeScale"))
						continue;
					if (timeScale != 1d) {
						if (!name.contains("-timeScale"+(float)timeScale))
							continue;
						if (scaleVelocities && !name.contains("-velScale"))
							continue;
						if (!scaleVelocities && name.contains("-velScale"))
							continue;
					}
					bbpDir = dir;
					bbpZipFile = zipFile;
				}
			}
		}
		Preconditions.checkNotNull(bbpDir);
		System.out.println("Located ref BBP dir: "+bbpDir.getAbsolutePath());
		System.out.println("\tInput file: "+bbpZipFile.getName());
		
		List<BBP_Site> bbpSites = BBP_Site.readFile(bbpDir);
		
		List<Site> sites = new ArrayList<>();
		for (BBP_Site site : bbpSites)
			sites.add(site.buildGMPE_Site());
		
		File catalogOutputDir = new File(outputDir, catalog.getCatalogDir().getName());
		Preconditions.checkState(catalogOutputDir.exists() || catalogOutputDir.mkdir());
		
		VelocityModel vm = RSQSimBBP_Config.detectVM(bbpDir);
		File vmDir = new File(catalogOutputDir, "bbp_"+vm.name());
		Preconditions.checkState(vmDir.exists() || vmDir.mkdir());
		
		Map<Scenario, RSQSimRotatedRupVariabilityScenarioPageGen> pageGensMap = new HashMap<>();
		HashSet<Integer> eventIDsSet = new HashSet<>();
		for (Scenario scenario : Scenario.values()) {
			File rotConfFile = new File(bbpDir, "rotation_config_"+scenario.getPrefix()+".csv");
			if (rotConfFile.exists()) {
				RSQSimRotatedRupVariabilityConfig config = RSQSimRotatedRupVariabilityConfig.loadCSV(
						catalog, rotConfFile, null, sites);
				
				BBP_RotatedRupSimLoader bbpLoader = new BBP_RotatedRupSimLoader(bbpZipFile, bbpSites, scenario.getPrefix());
				
				RSQSimRotatedRupVariabilityScenarioPageGen pageGen =
						new RSQSimRotatedRupVariabilityScenarioPageGen(catalog, scenario, config, bbpLoader, calcPeriods);
				
				pageGen.setGMPEs(refGMPEs);
				
				eventIDsSet.addAll(pageGen.getAllEventIDs());
				
				pageGensMap.put(scenario, pageGen);
			}
		}
		
		Map<Integer, RSQSimEvent> eventsMap = loadEvents(catalog, eventIDsSet);
		
		for (Scenario scenario : Scenario.values()) {
			if (!pageGensMap.containsKey(scenario))
				continue;
			System.out.println("Doing scenario: "+scenario);
			
			RSQSimRotatedRupVariabilityPageGen pageGen = pageGensMap.get(scenario);
			
			pageGen.setEventsMap(eventsMap);
			
			if (realData != null)
				pageGen.setRealEventData(ASK_EventData.getMatches(realData, null, null, scenario.getFaultStyle(), 30d), 100);
			
			String dirName = "rotated_ruptures_"+scenario.getPrefix();
			if (timeScale != 1d) {
				dirName += "_timeScale"+(float)timeScale;
				if (scaleVelocities)
					dirName += "_velScale";
			}
			File rotDir = new File(vmDir, dirName);
			System.out.println("Output dir: "+rotDir.getAbsolutePath());
			Preconditions.checkState(rotDir.exists() || rotDir.mkdir());
			
			List<String> methodSpecificLines = new ArrayList<>();
			
			methodSpecificLines.add("*NOTE: This page uses the SCEC BBP to simulate a 1-dimensional velocity structure. Thus we "
					+ "expect no path variability, and plots of path variabilitiy are included only as verification of the method.*");
			methodSpecificLines.add("");
			methodSpecificLines.add("[RSQSim Catalog Details](../#"+MarkdownUtils.getAnchorName(catalog.getName())+")");
			
			pageGen.generatePage(rotDir, periods, methodSpecificLines);
			pageGen.clearCaches();
		}
		
		catalog.writeMarkdownSummary(catalogOutputDir, true, false);
		RSQSimCatalog.writeCatalogsIndex(outputDir);
	}

}
