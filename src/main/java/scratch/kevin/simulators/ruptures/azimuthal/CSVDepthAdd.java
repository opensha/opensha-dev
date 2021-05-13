package scratch.kevin.simulators.ruptures.azimuthal;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.sha.simulators.RSQSimEvent;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Ints;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;

class CSVDepthAdd {

	public static void main(String[] args) throws IOException {
		File bbpDir = new File("/data/kevin/bbp/parallel");
		for (File dir : bbpDir.listFiles()) {
			String name = dir.getName();
			if (!name.contains("-azimuthal-"))
				continue;
			double spacing = name.contains("spacing5") ? 5d : 10d;
			double buffer = 100d;
			boolean positiveX = false;
			boolean gp = name.contains("-gp-");
			if (gp) {
				System.out.println("Found GP dir: "+name);
				for (Scenario scenario : Scenario.values()) {
					File locsFile = new File(dir, scenario.getPrefix()+"_locs.csv");
					if (!locsFile.exists())
						continue;
					File rupsFile = new File(dir, scenario.getPrefix()+"_rups.csv");
					AzimuthalSiteConfig<Integer> config = AzimuthalSiteConfig.loadGeneric(scenario, locsFile, rupsFile);
					if (Double.isNaN(config.getHypocenterDDW(config.getRuptures().get(0)))
							|| config.getWidth(config.getRuptures().get(0)) == 0d) {
						System.out.println("\tFixing for "+scenario.getName());
						double elemArea = 1d;
						GPAzimuthalSiteConfig newConfig = new GPAzimuthalSiteConfig(scenario, config.getRuptures().size(),
								elemArea, buffer, spacing, positiveX);
						Preconditions.checkState(newConfig.getGC2XYZ().size() == config.getGC2XYZ().size(),
								"Size mismatch. Old=%s, new=%s", config.getGC2XYZ().size(), newConfig.getGC2XYZ().size());
						for (Integer rupture : config.getRuptures()) {
							double das1 = config.getHypocenterDAS(rupture);
							double das2 = newConfig.getHypocenterDAS(rupture);
							Preconditions.checkState((float)das1 == (float)das2);
						}
						newConfig.writeRupturesCSV(rupsFile);
					}
				}
			} else {
				String catalogName = name.substring(name.indexOf("rundir"));
				catalogName = catalogName.substring(0, catalogName.indexOf("-"));
				if (catalogName.equals("rundir2585_1myrs"))
					catalogName = "rundir2585_1myr";
				System.out.println("Found catalog for: "+catalogName);
				File catalogDir = RSQSimCatalog.locateCatalog(catalogName);
				System.out.println("\tCatalog dir: "+catalogDir.getAbsolutePath());
				RSQSimCatalog catalog = null;
				for (Scenario scenario : Scenario.values()) {
					File locsFile = new File(dir, scenario.getPrefix()+"_locs.csv");
					if (!locsFile.exists())
						continue;
					File rupsFile = new File(dir, scenario.getPrefix()+"_rups.csv");
					AzimuthalSiteConfig<Integer> config = AzimuthalSiteConfig.loadGeneric(scenario, locsFile, rupsFile);
					if (Double.isNaN(config.getHypocenterDDW(config.getRuptures().get(0)))) {
						System.out.println("\tFixing for "+scenario.getName());
						if (catalog == null)
							catalog = new RSQSimCatalog(catalogDir, catalogName, FaultModels.FM3_1, DeformationModels.GEOLOGIC);
						List<RSQSimEvent> matches = catalog.loader().byIDs(Ints.toArray(config.getRuptures()));
						Map<Integer, RSQSimEvent> idRupMap = new HashMap<>();
						for (RSQSimEvent event : matches)
							idRupMap.put(event.getID(), event);
						List<RSQSimEvent> events = new ArrayList<>();
						for (Integer id : config.getRuptures())
							events.add(idRupMap.get(id));
						RSQSimAzimuthalSiteConfig newConfig = new RSQSimAzimuthalSiteConfig(
								catalog, scenario, events, buffer, spacing, positiveX);
						Preconditions.checkState(newConfig.getGC2XYZ().size() == config.getGC2XYZ().size(),
								"Size mismatch. Old=%s, new=%s", config.getGC2XYZ().size(), newConfig.getGC2XYZ().size());
						for (Integer rupture : config.getRuptures()) {
							double das1 = config.getHypocenterDAS(rupture);
							double das2 = newConfig.getHypocenterDAS(idRupMap.get(rupture));
							Preconditions.checkState((float)das1 == (float)das2);
						}
						newConfig.writeRupturesCSV(rupsFile);
					}
				}
			}
		}
	}

}
