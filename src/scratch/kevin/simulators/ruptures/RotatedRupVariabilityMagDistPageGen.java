package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.dom4j.DocumentException;
import org.opensha.commons.data.Site;
import org.opensha.commons.util.FileNameComparator;
import org.opensha.commons.util.IDPairing;
import org.opensha.commons.util.MarkdownUtils;
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
import scratch.kevin.simCompare.SimulationRotDProvider;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.RSQSimCatalog.Loader;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;
import scratch.kevin.simulators.ruptures.RotatedRupVariabilityConfig.RotationSpec;

public class RotatedRupVariabilityMagDistPageGen extends RotatedRupVariabilityPageGen {
	
	
	
	public enum RuptureType {
		VERT_STIKE_SLIP("Vertical Strike-Slip with Surface Rupture", "SS", "vert_ss_surface",
				new String[] { "Ztor=[0,1]", "Rake=[-180,-170] or [-10,10] or [170,180]",
						"Dip=90", "Linear rupture (max 5% deviation from ideal)"}) {
			
			private RSQSimCatalog prevCatalog;
			private Map<IDPairing, Double> horzDistanceCache;
			
			@Override
			public synchronized RuptureIdentifier getCriteria(RSQSimCatalog catalog, double centerMag, double deltaMag) {
				List<RuptureIdentifier> idens = new ArrayList<>();
				idens.add(new MagRangeRuptureIdentifier(centerMag-0.5*deltaMag, centerMag+0.5*deltaMag));
				idens.add(new DepthIden(Range.closed(0d, 1d), null));
				idens.add(FocalMechIden.builder().strikeSlip(10).forDip(90).build());
				if (catalog != prevCatalog)
					horzDistanceCache = null;
				if (horzDistanceCache == null) {
					horzDistanceCache = new HashMap<>();
					prevCatalog = catalog;
				}
				idens.add(new LinearRuptureIden(0.05d, true, horzDistanceCache));
				try {
					idens.add(new RSQSimTransValidIden(catalog.getTransitions(), catalog.getSlipVelocities()));
				} catch (Exception e) {
					System.out.println("Warning, couldn't force events with transitions. Missing trans file? "+e.getMessage());
				}
				return new LogicalAndRupIden(idens);
			}
		},
		REVERSE("Reverse, Dip=45", "Reverse", "reverse",
				new String[] { "Ztor=[0,5]", "Rake=[75,105]", "Dip=[35,55]"}) {
			@Override
			public RuptureIdentifier getCriteria(RSQSimCatalog catalog, double centerMag, double deltaMag) {
				List<RuptureIdentifier> idens = new ArrayList<>();
				idens.add(new MagRangeRuptureIdentifier(centerMag-0.5*deltaMag, centerMag+0.5*deltaMag));
				idens.add(new DepthIden(Range.closed(0d, 5d), null));
				idens.add(FocalMechIden.builder().forRake(75, 105).forDip(35, 55).build());
				try {
					idens.add(new RSQSimTransValidIden(catalog.getTransitions(), catalog.getSlipVelocities()));
				} catch (Exception e) {
					System.out.println("Warning, couldn't force events with transitions. Missing trans file? "+e.getMessage());
				}
				return new LogicalAndRupIden(idens);
			}
		};
		
		private String name;
		private String shortName;
		private String prefix;
		private String[] matchCriteria;

		private RuptureType(String name, String shortName, String prefix, String[] matchCriteria) {
			this.name = name;
			this.shortName = shortName;
			this.prefix = prefix;
			this.matchCriteria = matchCriteria;
		}
		
		public abstract RuptureIdentifier getCriteria(RSQSimCatalog catalog, double centerMag, double deltaMag);
		
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

	public RotatedRupVariabilityMagDistPageGen(RSQSimCatalog catalog, Map<Double, RotatedRupVariabilityConfig> magConfigs,
			Map<Double, SimulationRotDProvider<RotationSpec>> magProvs, RuptureType ruptureType) {
		super(catalog, magConfigs, magProvs);
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

	public static void main(String[] args) throws IOException, DocumentException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		File outputDir = new File("/home/kevin/git/rsqsim-analysis/catalogs");
		File bbpParallelDir = new File("/home/kevin/bbp/parallel");

		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance(baseDir);
		
		double[] periods = {3d, 5d, 7.5, 10d};
		boolean doExample = true;
		
		double[] highlightMags = {6.5, 7d, 7.5};
		float[] highlightDists = {20f, 40f, 80f, 160f};
		
		System.out.println("Catalog: "+catalog.getName());
		// find BBP parallel dir
		String catalogDirName = catalog.getCatalogDir().getName();
		if (catalogDirName.startsWith("JG_"))
			// I sometimes modify Jacqui's directory names
			catalogDirName = catalogDirName.substring(3);
		File bbpDir = null;
		File bbpZipFile = null;
		File[] allBBPDirs = bbpParallelDir.listFiles();
		Arrays.sort(allBBPDirs, new FileNameComparator());
		for (File dir : allBBPDirs) {
			String name = dir.getName();
			if (dir.isDirectory() && name.contains(catalogDirName) && name.contains("-rotatedRupsMagDist-")) {
				File zipFile = new File(dir, "results.zip");
				if (!zipFile.exists())
					zipFile = new File(dir, "results_rotD.zip");
				if (zipFile.exists()) {
					bbpDir = dir;
					bbpZipFile = zipFile;
				}
			}
		}
		Preconditions.checkNotNull(bbpDir);
		System.out.println("Located ref BBP dir: "+bbpDir.getAbsolutePath());
		System.out.println("\tInput file: "+bbpZipFile.getName());
		
		List<BBP_Site> bbpSites = RSQSimBBP_Config.getStandardSites().subList(0, 1);
		
		List<Site> sites = new ArrayList<>();
		for (BBP_Site site : bbpSites)
			sites.add(site.buildGMPE_Site(RSQSimBBP_Config.VM));
		
		File catalogOutputDir = new File(outputDir, catalog.getCatalogDir().getName());
		Preconditions.checkState(catalogOutputDir.exists() || catalogOutputDir.mkdir());
		
		Map<RuptureType, RotatedRupVariabilityMagDistPageGen> pageGensMap = new HashMap<>();
		for (RuptureType rupType : RuptureType.values()) {
			Map<Double, RotatedRupVariabilityConfig> configsMap = new HashMap<>();
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

				configsMap.put(mag, RotatedRupVariabilityConfig.loadCSV(catalog, file, null, sites));;
				loadersMap.put(mag, new BBP_RotatedRupSimLoader(bbpZipFile, bbpSites, rupType.getMagPrefix(mag)));
			}
			
			if (configsMap.isEmpty())
				continue;
			
			RotatedRupVariabilityMagDistPageGen pageGen =
					new RotatedRupVariabilityMagDistPageGen(catalog, configsMap, loadersMap, rupType);
			
			pageGensMap.put(rupType, pageGen);
		}
		
		Map<RuptureType, RSQSimEvent> exampleEventsMap = null;
		if (doExample) {
			Map<RuptureType, Integer> exampleIDs = new HashMap<>();
			for (RotatedRupVariabilityMagDistPageGen pageGen : pageGensMap.values())
				exampleIDs.put(pageGen.ruptureType, pageGen.getEventIDs(pageGen.minMag).get(0));
			System.out.println("Loading "+exampleIDs.size()+" example ruptures");
			List<RSQSimEvent> events = catalog.loader().byIDs(Ints.toArray(exampleIDs.values()));
			exampleEventsMap = new HashMap<>();
			for (RuptureType rupType : exampleIDs.keySet()) {
				int exampleID = exampleIDs.get(rupType);
				for (RSQSimEvent event : events)
					if (event.getID() == exampleID)
						exampleEventsMap.put(rupType, event);
			}
		}
		
		for (RuptureType rupType : pageGensMap.keySet()) {
			System.out.println("Doing scenario: "+rupType);
			
			RotatedRupVariabilityPageGen pageGen = pageGensMap.get(rupType);
			
			if (doExample) {
				RSQSimEvent exampleRup = exampleEventsMap.get(rupType);
				pageGen.setExampleRupture(exampleRup, pageGen.getSites().get(0), 5);
			}
			
			File rotDir = new File(catalogOutputDir, "rotated_ruptures_mag_dist_"+rupType.getPrefix());
			Preconditions.checkState(rotDir.exists() || rotDir.mkdir());
			
			List<String> methodSpecificLines = new ArrayList<>();
			
			methodSpecificLines.add("*NOTE: This page uses the SCEC BBP to simulate a 1-dimensional velocity structure. Thus we "
					+ "expect no path variability, and plots of path variabilitiy are included only as verification of the method.*");
			methodSpecificLines.add("");
			methodSpecificLines.add("[RSQSim Catalog Details](../#"+MarkdownUtils.getAnchorName(catalog.getName())+")");
			
			pageGen.generatePage(rotDir, periods, methodSpecificLines, highlightMags, highlightDists);
		}
		
		catalog.writeMarkdownSummary(catalogOutputDir, true, false);
		RSQSimCatalog.writeCatalogsIndex(outputDir);
	}

}
