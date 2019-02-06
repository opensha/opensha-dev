package scratch.kevin.simulators.ruptures;

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
import org.opensha.sha.simulators.RSQSimEvent;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Ints;

import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.simCompare.SimulationRotDProvider;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;
import scratch.kevin.simulators.ruptures.RotatedRupVariabilityConfig.RotationSpec;

public class RotatedRupVariabilityScenarioPageGen extends RotatedRupVariabilityPageGen {

	private Scenario scenario;

	public RotatedRupVariabilityScenarioPageGen(RSQSimCatalog catalog, Scenario scenario,
			RotatedRupVariabilityConfig config, SimulationRotDProvider<RotationSpec> prov) {
		super(catalog, config, prov);
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

	public static void main(String[] args) throws IOException, DocumentException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		File outputDir = new File("/home/kevin/git/rsqsim-analysis/catalogs");
		File bbpParallelDir = new File("/home/kevin/bbp/parallel");

		RSQSimCatalog catalog = Catalogs.BRUCE_2585.instance(baseDir);
		
		double[] periods = {3d, 5d, 7.5, 10d};
		boolean doExample = true;
		
		System.out.println("Catalog: "+catalog.getName());
		// find BBP parallel dir
		String catalogDirName = catalog.getCatalogDir().getName();
		if (catalogDirName.startsWith("JG_"))
			// I sometimes modify Jacqui's directory names
			catalogDirName = catalogDirName.substring(3);
		File bbpDir = null;
		File bbpZipFile = null;
		File[] allBBPDirs = bbpParallelDir.listFiles();
		VelocityModel vm = VelocityModel.LA_BASIN;
		Arrays.sort(allBBPDirs, new FileNameComparator());
		for (File dir : allBBPDirs) {
			String name = dir.getName();
			if (dir.isDirectory() && name.contains(catalogDirName) && name.contains("-rotatedRups-")) {
				File zipFile = new File(dir, "results.zip");
				if (!zipFile.exists())
					zipFile = new File(dir, "results_rotD.zip");
				if (zipFile.exists()) {
					bbpDir = dir;
					bbpZipFile = zipFile;
					if (name.contains("-vm")) {
						String vmStr = name.substring(name.indexOf("-vm")+3);
						if (vmStr.contains("-"))
							vmStr = vmStr.substring(0, vmStr.indexOf("-"));
						vm = VelocityModel.valueOf(vmStr);
					}
				}
			}
		}
		Preconditions.checkNotNull(bbpDir);
		System.out.println("Located ref BBP dir: "+bbpDir.getAbsolutePath());
		System.out.println("\tInput file: "+bbpZipFile.getName());
		
		List<BBP_Site> bbpSites = BBP_Site.readFile(bbpDir);
		
		List<Site> sites = new ArrayList<>();
		for (BBP_Site site : bbpSites)
			sites.add(site.buildGMPE_Site(vm));
		
		File catalogOutputDir = new File(outputDir, catalog.getCatalogDir().getName());
		Preconditions.checkState(catalogOutputDir.exists() || catalogOutputDir.mkdir());
		
		Map<Scenario, RotatedRupVariabilityScenarioPageGen> pageGensMap = new HashMap<>();
		HashSet<Integer> eventIDsSet = new HashSet<>();
		for (Scenario scenario : Scenario.values()) {
			File rotConfFile = new File(bbpDir, "rotation_config_"+scenario.getPrefix()+".csv");
			if (rotConfFile.exists()) {
				RotatedRupVariabilityConfig config = RotatedRupVariabilityConfig.loadCSV(catalog, rotConfFile, null, sites);
				
				BBP_RotatedRupSimLoader bbpLoader = new BBP_RotatedRupSimLoader(bbpZipFile, bbpSites, scenario.getPrefix());
				
				RotatedRupVariabilityScenarioPageGen pageGen =
						new RotatedRupVariabilityScenarioPageGen(catalog, scenario, config, bbpLoader);
				
				eventIDsSet.addAll(pageGen.getAllEventIDs());
				
				pageGensMap.put(scenario, pageGen);
			}
		}
		
		Map<Integer, RSQSimEvent> eventsMap = loadEvents(catalog, eventIDsSet);
		
		for (Scenario scenario : pageGensMap.keySet()) {
			System.out.println("Doing scenario: "+scenario);
			
			RotatedRupVariabilityPageGen pageGen = pageGensMap.get(scenario);
			
			pageGen.setEventsMap(eventsMap);
			
			File rotDir = new File(catalogOutputDir, "rotated_ruptures_"+scenario.getPrefix());
			Preconditions.checkState(rotDir.exists() || rotDir.mkdir());
			
			List<String> methodSpecificLines = new ArrayList<>();
			
			methodSpecificLines.add("*NOTE: This page uses the SCEC BBP to simulate a 1-dimensional velocity structure. Thus we "
					+ "expect no path variability, and plots of path variabilitiy are included only as verification of the method.*");
			methodSpecificLines.add("");
			methodSpecificLines.add("[RSQSim Catalog Details](../#"+MarkdownUtils.getAnchorName(catalog.getName())+")");
			
			pageGen.generatePage(rotDir, periods, methodSpecificLines);
		}
		
		catalog.writeMarkdownSummary(catalogOutputDir, true, false);
		RSQSimCatalog.writeCatalogsIndex(outputDir);
	}

}
