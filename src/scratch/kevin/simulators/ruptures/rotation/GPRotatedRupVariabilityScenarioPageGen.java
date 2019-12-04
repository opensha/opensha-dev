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

public class GPRotatedRupVariabilityScenarioPageGen extends GPRotatedRupVariabilityPageGen {

	private Scenario scenario;

	public GPRotatedRupVariabilityScenarioPageGen(Scenario scenario,
			GPRotatedRupVariabilityConfig config, SimulationRotDProvider<RotationSpec> prov, double[] calcPeriods) {
		super(config, scenario.getMagnitude(), prov, calcPeriods);
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
	protected Scenario getBBP_PartB_Scenario(RotatedRupVariabilityConfig<GPRotatedRupture> config) {
		return scenario;
	}

	@SuppressWarnings("unused")
	public static void main(String[] args) throws IOException, DocumentException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		File outputDir = new File("/home/kevin/git/rsqsim-analysis/gp_comparisons/");
		File bbpParallelDir = new File("/home/kevin/bbp/parallel");
		
		File bbpDir = new File(bbpParallelDir,
				"2019_12_03-gp-rotatedRups-3scenarios-20.0km-4srcAz-1siteSrcAz-10rups-patchArea1.0-vmLA_BASIN_500-noHF-1site");
		File bbpZipFile = new File(bbpDir, "results_rotD.zip");
		
		double[] calcPeriods = {1d, 2d, 3d, 4d, 5d, 7.5, 10d};
		double[] periods = {3d, 5d, 7.5, 10d};
		
		Map<Integer, List<ASK_EventData>> realData = ASK_EventData.load(1d);
		
		NGAW2_WrapperFullParam[] refGMPEs = { new NGAW2_Wrappers.ASK_2014_Wrapper(), new NGAW2_Wrappers.BSSA_2014_Wrapper(),
				new NGAW2_Wrappers.CB_2014_Wrapper(), new NGAW2_Wrappers.CY_2014_Wrapper()};
//		NGAW2_WrapperFullParam[] refGMPEs = { new NGAW2_Wrappers.BSSA_2014_Wrapper() };
		
		Preconditions.checkNotNull(bbpDir);
		System.out.println("Located ref BBP dir: "+bbpDir.getAbsolutePath());
		System.out.println("\tInput file: "+bbpZipFile.getName());
		
		List<BBP_Site> bbpSites = BBP_Site.readFile(bbpDir);
		
		List<Site> sites = new ArrayList<>();
		for (BBP_Site site : bbpSites)
			sites.add(site.buildGMPE_Site());
		
		File parentDir = new File(outputDir, bbpDir.getName());
		Preconditions.checkState(parentDir.exists() || parentDir.mkdir());
		
		for (Scenario scenario : Scenario.values()) {
			File rotConfFile = new File(bbpDir, "rotation_config_"+scenario.getPrefix()+".csv");
			if (!rotConfFile.exists())
				continue;
			System.out.println("Doing scenario: "+scenario);
			
			GPRotatedRupVariabilityConfig config = GPRotatedRupVariabilityConfig.loadCSV(
					rotConfFile, new File(bbpDir, "rups_"+scenario.getPrefix()), sites);
			
			BBP_RotatedRupSimLoader bbpLoader = new BBP_RotatedRupSimLoader(bbpZipFile, bbpSites, scenario.getPrefix());
			
			GPRotatedRupVariabilityScenarioPageGen pageGen =
					new GPRotatedRupVariabilityScenarioPageGen(scenario, config, bbpLoader, calcPeriods);
			
			pageGen.setGMPEs(refGMPEs);
			
			if (realData != null)
				pageGen.setRealEventData(ASK_EventData.getMatches(realData, null, null, scenario.getFaultStyle(), 30d), 100);
			File rotDir = new File(parentDir, scenario.getPrefix());
			System.out.println("Output dir: "+rotDir.getAbsolutePath());
			Preconditions.checkState(rotDir.exists() || rotDir.mkdir());
			
			List<String> methodSpecificLines = new ArrayList<>();
			
			methodSpecificLines.add("This pages uses the Graves & Pitarka (2016) rupture generator. Rupture surfaces are determined "
					+ "by first computing the Wells & Coppersmith median length for the magnitude, then a down dip width using that "
					+ "length and the area from Somerville (2006). Hypocenters are randomly distributed both down dip and along strike.");
			methodSpecificLines.add("");
			methodSpecificLines.add("*NOTE: This page uses the SCEC BBP to simulate a 1-dimensional velocity structure. Thus we "
					+ "expect no path variability, and plots of path variabilitiy are included only as verification of the method.*");
			
			pageGen.generatePage(rotDir, periods, methodSpecificLines);
			pageGen.clearCaches();
		}
	}

}
