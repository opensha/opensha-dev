package scratch.kevin.simulators.ruptures.rotation;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;

import scratch.kevin.bbp.BBP_SimZipLoader;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.simCompare.SimulationRotDProvider;
import scratch.kevin.simulators.ruptures.rotation.RotatedRupVariabilityConfig.RotationSpec;

public class BBP_RotatedRupSimLoader extends BBP_SimZipLoader implements SimulationRotDProvider<RotationSpec> {
	
	private Map<String, DiscretizedFunc> rdMap;
	private Map<String, BBP_Site> bbpSiteMap;
	private String scenarioPrefix;

	public BBP_RotatedRupSimLoader(File file, List<BBP_Site> sites, String scenarioPrefix) throws ZipException, IOException {
		this(new ZipFile(file), sites, scenarioPrefix);
	}

	public BBP_RotatedRupSimLoader(ZipFile zip, List<BBP_Site> sites, String scenarioPrefix) {
		super(zip, sites);
		rdMap = new HashMap<>();
		this.scenarioPrefix = scenarioPrefix;
		bbpSiteMap = new HashMap<>();
		for (BBP_Site site : sites)
			bbpSiteMap.put(site.getName(), site);
	}
	
	public DiscretizedFunc getRotD50(RotationSpec rotation) throws IOException {
		String dirName = scenarioPrefix+"_"+rotation.getPrefix();
		if (rdMap.containsKey(dirName))
			return rdMap.get(dirName);
		DiscretizedFunc rd = readRotD50(bbpSiteMap.get(rotation.site.getName()), dirName);
		rdMap.putIfAbsent(dirName, rd);
		return rd;
	}

	@Override
	public String getName() {
		return "RSQSim/BBP";
	}

	@Override
	public DiscretizedFunc getRotD50(Site site, RotationSpec rupture, int index) throws IOException {
		return getRotD50(rupture);
	}

	@Override
	public DiscretizedFunc getRotD100(Site site, RotationSpec rupture, int index) throws IOException {
		throw new UnsupportedOperationException("not implemented");
	}

	@Override
	public DiscretizedFunc[] getRotD(Site site, RotationSpec rupture, int index) throws IOException {
		throw new UnsupportedOperationException("not implemented");
	}

	@Override
	public DiscretizedFunc getRotDRatio(Site site, RotationSpec rupture, int index) throws IOException {
		throw new UnsupportedOperationException("not implemented");
	}

	@Override
	public int getNumSimulations(Site site, RotationSpec rupture) {
		return 1;
	}

	@Override
	public Location getHypocenter(RotationSpec rupture, int index) {
		throw new UnsupportedOperationException("not implemented");
	}

	@Override
	public Collection<RotationSpec> getRupturesForSite(Site site) {
		throw new UnsupportedOperationException("not implemented");
	}

	@Override
	public double getAnnualRate(RotationSpec rupture) {
		throw new UnsupportedOperationException("not implemented");
	}

	@Override
	public double getMinimumCurvePlotRate(Site site) {
		throw new UnsupportedOperationException("not implemented");
	}

	@Override
	public double getMagnitude(RotationSpec rupture) {
		throw new UnsupportedOperationException("not implemented");
	}

}
