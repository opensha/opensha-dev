package scratch.kevin.simulators.ruptures.rotation;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.dom4j.DocumentException;
import org.opensha.commons.data.Site;
import org.opensha.commons.util.FileNameComparator;
import org.opensha.commons.util.IDPairing;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.ngaw2.FaultStyle;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_WrapperFullParam;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_Wrappers;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.iden.DepthIden;
import org.opensha.sha.simulators.iden.FocalMechIden;
import org.opensha.sha.simulators.iden.LinearRuptureIden;
import org.opensha.sha.simulators.iden.LogicalAndRupIden;
import org.opensha.sha.simulators.iden.MagRangeRuptureIdentifier;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.srf.RSQSimTransValidIden;

import com.google.common.base.Preconditions;
import com.google.common.collect.Range;
import com.google.common.primitives.Ints;

import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.simCompare.SimulationRotDProvider;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.RSQSimCatalog.Loader;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.FilterMethod;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;
import scratch.kevin.simulators.ruptures.RSQSimBBP_Config;
import scratch.kevin.simulators.ruptures.rotation.RotatedRupVariabilityConfig.RotationSpec;

public class RSQSimRotatedRupVariabilityMagDistPageGen extends RSQSimRotatedRupVariabilityPageGen {
	
	public enum RuptureType {
		VERT_STIKE_SLIP("Vertical Strike-Slip with Surface Rupture", "SS", "vert_ss_surface",
				Range.closed(0d, 1d), // Ztor range
				FaultStyle.STRIKE_SLIP, 10, // style, rakeTolerance
				90, 0, // dip, dipTolerance
				0.05, true), // linearFract, linearRelative,
		REVERSE("Reverse, Dip=45", "Reverse", "reverse",
				Range.closed(0d, 5d), // Ztor range
				FaultStyle.REVERSE, 15, // style, rakeTolerance
				45, 10, // dip, dipTolerance
				null, false), // linearFract, linearRelative,
		NORMAL("Normal, Dip=45", "Normal", "normal",
				Range.closed(0d, 5d), // Ztor range
				FaultStyle.NORMAL, 15, // style, rakeTolerance
				45, 10, // dip, dipTolerance
				null, false); // linearFract, linearRelative
		
		private String name;
		private String shortName;
		private String prefix;
		private String[] matchCriteria;
		private Range<Double> zTorRange;
		private FaultStyle style;
		private int rakeTolerance;
		private double dip;
		private int dipTolerance;
		private Double linearFract;
		private boolean linearRelative;

		private RuptureType(String name, String shortName, String prefix,
				Range<Double> zTorRange, FaultStyle style, int rakeTolerance, double dip, int dipTolerance,
				Double linearFract, boolean linearRelative) {
			this.name = name;
			this.shortName = shortName;
			this.prefix = prefix;
			this.zTorRange = zTorRange;
			this.style = style;
			this.rakeTolerance = rakeTolerance;
			this.dip = dip;
			this.dipTolerance = dipTolerance;
			this.linearFract = linearFract;
			this.linearRelative = linearRelative;
			this.matchCriteria = BBP_PartBValidationConfig.buildMatchCriteria(Double.NaN, Double.NaN, zTorRange, style,
					rakeTolerance, dip, dipTolerance, linearFract, linearRelative);
		}
		
		public RuptureIdentifier getCriteria(RSQSimCatalog catalog, double centerMag, double deltaMag) {
			return new LogicalAndRupIden(BBP_PartBValidationConfig.getIdentifiers(
					catalog, centerMag, deltaMag*0.5, zTorRange, style, rakeTolerance, dip, dipTolerance, linearFract, linearRelative));
		}
		
		public List<RSQSimEvent> getMatches(RSQSimCatalog catalog, List<RSQSimEvent> events, double centerMag, double deltaMag)
				throws IOException {
			return getCriteria(catalog, centerMag, deltaMag).getMatches(events);
		}
		
		public List<RSQSimEvent> getMatches(RSQSimCatalog catalog, int skipYears, double centerMag, double deltaMag)
				throws IOException {
			Loader loader = catalog.loader().skipYears(skipYears);
			return loader.matches(getCriteria(catalog, centerMag, deltaMag)).load();
		}

		public String getName() {
			return name;
		}

		public String getShortName() {
			return shortName;
		}
		
		public String getPrefix() {
			return prefix;
		}

		public String getMagPrefix(double mag) {
			String magStr = optionalDigitDF.format(mag).replace('.', 'p');
			return "m"+magStr+"_"+prefix;
		}

		public String[] getMatchCriteria() {
			return matchCriteria;
		}
	}
	
	private RuptureType ruptureType;
	private double minMag;
	private Map<Double, RSQSimRotatedRupVariabilityConfig> magConfigs;

	public RSQSimRotatedRupVariabilityMagDistPageGen(RSQSimCatalog catalog, FilterMethod filter, Map<Double, RSQSimRotatedRupVariabilityConfig> magConfigs,
			Map<Double, SimulationRotDProvider<RotationSpec>> magProvs, RuptureType ruptureType, double[] calcPeriods) {
		super(catalog, filter, magConfigs, magProvs, calcPeriods);
		this.magConfigs = magConfigs;
		minMag = Double.POSITIVE_INFINITY;
		for (Double mag : magConfigs.keySet())
			minMag = Double.min(minMag, mag);
		this.ruptureType = ruptureType;
	}

	@Override
	protected String getScenarioName() {
		return ruptureType.getName();
	}

	@Override
	protected String getScenarioShortName() {
		return ruptureType.getShortName();
	}

	@Override
	protected String[] getScenarioMatchCriteria() {
		return ruptureType.getMatchCriteria();
	}

	@Override
	protected Scenario getBBP_PartB_Scenario(RotatedRupVariabilityConfig config) {
		Double mag = null;
		for (Double m : magConfigs.keySet()) {
			if (config == magConfigs.get(m)) {
				mag = m;
				break;
			}
		}
		Preconditions.checkState(mag != null);
		for (Scenario scenario : Scenario.values())
			if (ruptureType.style == scenario.getFaultStyle() && mag.floatValue() == (float)scenario.getMagnitude())
				return scenario;
		return null;
	}

	public static void main(String[] args) throws IOException, DocumentException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		File outputDir = new File("/home/kevin/git/rsqsim-analysis/catalogs");
		File bbpParallelDir = new File("/home/kevin/bbp/parallel");

		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance(baseDir);
//		RSQSimCatalog catalog = Catalogs.BRUCE_2740.instance(baseDir);
		
		VelocityModel forceVM = null;
		
		NGAW2_WrapperFullParam[] refGMPEs = { new NGAW2_Wrappers.ASK_2014_Wrapper(), new NGAW2_Wrappers.BSSA_2014_Wrapper(),
				new NGAW2_Wrappers.CB_2014_Wrapper(), new NGAW2_Wrappers.CY_2014_Wrapper()};
//		NGAW2_WrapperFullParam[] refGMPEs = { new NGAW2_Wrappers.BSSA_2014_Wrapper() };

		double[] calcPeriods = {1d, 2d, 3d, 4d, 5d, 7.5, 10d};
		double[] periods = {3d, 5d, 7.5, 10d};
		
		double[] highlightMags = {6.5, 7d, 7.5};
		float[] highlightDists = {20f, 60f, 120f};
		
		System.out.println("Catalog: "+catalog.getName());
		// find BBP parallel dir
		String catalogDirName = catalog.getCatalogDir().getName();
		if (catalogDirName.startsWith("JG_"))
			// I sometimes modify Jacqui's directory names
			catalogDirName = catalogDirName.substring(3);
		VelocityModel vm = null;
		File[] allBBPDirs = bbpParallelDir.listFiles();
		Arrays.sort(allBBPDirs, new FileNameComparator());
		List<File> matchingZipFiles = new ArrayList<>();
		FilterMethod filter = null;
		for (File dir : allBBPDirs) {
			String name = dir.getName();
			if (dir.isDirectory() && name.contains(catalogDirName) && name.contains("-rotatedRupsMagDist-")) {
				File zipFile = new File(dir, "results.zip");
				if (!zipFile.exists())
					zipFile = new File(dir, "results_rotD.zip");
				if (zipFile.exists()) {
					VelocityModel newVM = RSQSimBBP_Config.detectVM(dir);
					if (forceVM != null && forceVM != newVM)
						continue;
					FilterMethod myFilter = FilterMethod.fromDirName(name);
					Preconditions.checkState(filter == null || filter == myFilter, "filter mismatch, TODO");
					filter = myFilter;
					if (newVM != vm && !matchingZipFiles.isEmpty()) {
						System.out.println("We switched VMs, only using most recent");
						matchingZipFiles.clear();
						vm = newVM;
					}
					matchingZipFiles.add(zipFile);
				}
			}
		}
		Preconditions.checkState(!matchingZipFiles.isEmpty());
		for (File zipFile : matchingZipFiles) {
			System.out.println("Located ref BBP dir: "+zipFile.getParentFile().getAbsolutePath());
			System.out.println("\tInput file: "+zipFile.getName());
		}
		
		BBP_Site singleSite = RSQSimBBP_Config.getStandardSites().get(0);
		// fix the Vs30 value
		singleSite = new BBP_Site(singleSite.getName(), singleSite.getLoc(), vm.getVs30(),
				singleSite.getLoPassFreq(), singleSite.getHiPassFreq());
		List<BBP_Site> bbpSites = new ArrayList<>();
		bbpSites.add(singleSite);
		
		List<Site> sites = new ArrayList<>();
		for (BBP_Site site : bbpSites)
			sites.add(site.buildGMPE_Site(vm));
		
		File catalogOutputDir = new File(outputDir, catalog.getCatalogDir().getName());
		Preconditions.checkState(catalogOutputDir.exists() || catalogOutputDir.mkdir());

		File vmDir = new File(catalogOutputDir, "bbp_"+vm.name());
		Preconditions.checkState(vmDir.exists() || vmDir.mkdir());
		
		Map<RuptureType, RSQSimRotatedRupVariabilityMagDistPageGen> pageGensMap = new HashMap<>();
		HashSet<Integer> eventIDsSet = new HashSet<>();
		// reverse, latest first
		Collections.reverse(matchingZipFiles);
		for (File bbpZipFile : matchingZipFiles) {
			File bbpDir = bbpZipFile.getParentFile();
			System.out.println("Loading from: "+bbpDir.getName());
			for (RuptureType rupType : RuptureType.values()) {
				if (pageGensMap.containsKey(rupType))
					continue;
				Map<Double, RSQSimRotatedRupVariabilityConfig> configsMap = new HashMap<>();
				Map<Double, SimulationRotDProvider<RotationSpec>> loadersMap = new HashMap<>();
				
				for (File file : bbpDir.listFiles()) {
					String name = file.getName();
					if (!name.startsWith("rotation_config") || !name.contains(rupType.getPrefix()) || !name.endsWith(".csv") || !name.contains("_m"))
						continue;
					String magStr = name.substring(name.indexOf("_m")+2);
					magStr = magStr.replace('p', '.');
					magStr = magStr.substring(0, magStr.indexOf("_"));
					double mag = Double.parseDouble(magStr);
					System.out.println("Located CSV for "+rupType.getName()+", M"+(float)mag+": "+name);

					configsMap.put(mag, RSQSimRotatedRupVariabilityConfig.loadCSV(catalog, file, null, sites));;
					loadersMap.put(mag, new BBP_RotatedRupSimLoader(bbpZipFile, bbpSites, rupType.getMagPrefix(mag)));
				}
				
				if (configsMap.isEmpty())
					continue;
				
				RSQSimRotatedRupVariabilityMagDistPageGen pageGen =
						new RSQSimRotatedRupVariabilityMagDistPageGen(catalog, filter, configsMap, loadersMap, rupType, calcPeriods);
				
				pageGen.setGMPEs(refGMPEs);
				
				eventIDsSet.addAll(pageGen.getAllEventIDs());
				
				pageGensMap.put(rupType, pageGen);
			}
		}
		
		Map<Integer, RSQSimEvent> eventsMap = loadEvents(catalog, eventIDsSet);
		
		for (RuptureType rupType : pageGensMap.keySet()) {
			System.out.println("Doing scenario: "+rupType);
			
			RSQSimRotatedRupVariabilityPageGen pageGen = pageGensMap.get(rupType);
			
			pageGen.setEventsMap(eventsMap);
			
			File rotDir = new File(vmDir, "rotated_ruptures_mag_dist_"+rupType.getPrefix()+"_filter_"+filter.getPrefix());
			Preconditions.checkState(rotDir.exists() || rotDir.mkdir());
			
			List<String> methodSpecificLines = new ArrayList<>();
			
			methodSpecificLines.add("*NOTE: This page uses the SCEC BBP to simulate a 1-dimensional velocity structure. Thus we "
					+ "expect no path variability, and plots of path variabilitiy are included only as verification of the method.*");
			methodSpecificLines.add("");
			methodSpecificLines.add("[RSQSim Catalog Details](../#"+MarkdownUtils.getAnchorName(catalog.getName())+")");
			
			pageGen.generatePage(rotDir, periods, methodSpecificLines, highlightMags, highlightDists);
			pageGen.clearCaches();
		}
		
		catalog.writeMarkdownSummary(catalogOutputDir, true, false);
		RSQSimCatalog.writeCatalogsIndex(outputDir);
	}

}
