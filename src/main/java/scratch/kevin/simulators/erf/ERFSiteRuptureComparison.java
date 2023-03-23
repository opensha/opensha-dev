package scratch.kevin.simulators.erf;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.Site;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.simulators.SimulatorElement;

import com.google.common.base.Preconditions;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.erf.RSQSimSectBundledERF.RSQSimProbEqkRup;
import scratch.kevin.simulators.erf.RSQSimSectBundledERF.RSQSimSectBundledSource;

public class ERFSiteRuptureComparison {

	public static void main(String[] args) throws IOException {
		File localBaseDir = new File("/home/kevin/Simulators/catalogs");
		
		RSQSimCatalog catalog1 = Catalogs.BRUCE_2585_1MYR.instance(localBaseDir);
		RSQSimCatalog catalog2 = Catalogs.BRUCE_2740.instance(localBaseDir);
		
		Preconditions.checkState(catalog1.getElements().size() == catalog2.getElements().size());
		
		File mappingFile1 = new File(catalog1.getCatalogDir(), "erf_mappings.bin");
		RSQSimSectBundledERF erf1 = new RSQSimSectBundledERF(mappingFile1, null,
				catalog1.getFaultModel(), catalog1.getDeformationModel(), catalog1.getSubSects(), catalog1.getElements());
		erf1.updateForecast();
		
		File mappingFile2 = new File(catalog2.getCatalogDir(), "erf_mappings.bin");
		RSQSimSectBundledERF erf2 = new RSQSimSectBundledERF(mappingFile2, null,
				catalog2.getFaultModel(), catalog2.getDeformationModel(), catalog2.getSubSects(), catalog2.getElements());
		erf2.updateForecast();
		
		double distCutoff = 200d;
		Map<String, Location> sitesMap = new HashMap<>();
		sitesMap.put("USC", new Location(34.0192, -118.286));
		sitesMap.put("PAS", new Location(34.148426, -118.17119));
		sitesMap.put("SBSM", new Location(34.064986, -117.29201));
		sitesMap.put("WNGC", new Location(34.041824, -118.0653));
		sitesMap.put("STNI", new Location(33.93088, -118.17881));
		sitesMap.put("LAPD", new Location(34.557, -118.125));
		sitesMap.put("s119", new Location(34.55314, -118.72826));
		sitesMap.put("s279", new Location(34.37809, -118.34757));
		sitesMap.put("s480", new Location(34.15755, -117.87389));
		
		for (String siteName : sitesMap.keySet()) {
			System.out.println(siteName);
			Location siteLoc = sitesMap.get(siteName);
			
			List<RSQSimProbEqkRup> rups1 = getSiteRuptures(erf1, siteLoc, distCutoff);
			HashSet<Integer> elems1 = getAllElementIDs(rups1);
			
			List<RSQSimProbEqkRup> rups2 = getSiteRuptures(erf2, siteLoc, distCutoff);
			HashSet<Integer> elems2 = getAllElementIDs(rups2);
			
			System.out.println("\t"+catalog1.getName()+":");
			System.out.println("\t\t"+rups1.size()+" ruptures");
			System.out.println("\t\t"+elems1.size()+" elements (of "+catalog1.getElements().size()+")");
			System.out.println("\t\t"+getNumUnique(elems1, elems2)+" unique elements");
			System.out.println("\tfurthest elem dist: "+(float)getFurthestElemDist(siteLoc, elems1, catalog1)+" km");
			System.out.println("\tfurthest possible dist: "+(float)getFurthestElemDist(siteLoc, null, catalog1)+" km");
			
			System.out.println("\t"+catalog2.getName()+":");
			System.out.println("\t\t"+rups2.size()+" ruptures");
			System.out.println("\t\t"+elems2.size()+" elements (of "+catalog2.getElements().size()+")");
			System.out.println("\t\t"+getNumUnique(elems2, elems1)+" unique elements");
			System.out.println("\tfurthest elem dist: "+(float)getFurthestElemDist(siteLoc, elems2, catalog2)+" km");
			System.out.println("\tfurthest possible dist: "+(float)getFurthestElemDist(siteLoc, null, catalog2)+" km");
		}
	}
	
	private static List<RSQSimProbEqkRup> getSiteRuptures(RSQSimSectBundledERF erf, Location location, double distCutoff) {
		List<RSQSimProbEqkRup> rups = new ArrayList<>();
		
		Site site = new Site(location);
		
		for (int sourceID=0; sourceID<erf.getNumSources(); sourceID++) {
			RSQSimSectBundledSource source = erf.getSource(sourceID);
			double sourceDist = source.getMinDistance(site);
			if (sourceDist > distCutoff)
				continue;
			for (int rupID=0; rupID<source.getNumRuptures(); rupID++)
				rups.add(source.getRupture(rupID));
		}
		
		return rups;
	}
	
	private static HashSet<Integer> getAllElementIDs(List<RSQSimProbEqkRup> rups) {
		HashSet<Integer> elements = new HashSet<>();
		
		for (RSQSimProbEqkRup rup : rups)
			for (SimulatorElement elem : rup.getElements())
				elements.add(elem.getID());
		
		return elements;
	}
	
	private static int getNumUnique(HashSet<Integer> elems, HashSet<Integer> otherElems) {
		int count = 0;
		
		for (Integer elem : elems)
			if (!otherElems.contains(elem))
				count++;
		
		return count;
	}
	
	private static double getFurthestElemDist(Location loc, HashSet<Integer> ids, RSQSimCatalog catalog) throws IOException {
		double maxDist = 0d;
		
		List<SimulatorElement> elems = catalog.getElements();
		
		if (ids == null) {
			for (SimulatorElement elem : elems) {
				maxDist = Math.max(maxDist, LocationUtils.horzDistanceFast(loc, elem.getCenterLocation()));
			}
		} else {
			for (int id : ids) {
				SimulatorElement elem = elems.get(id-1);
				Preconditions.checkState(elem.getID() == id);
				maxDist = Math.max(maxDist, LocationUtils.horzDistanceFast(loc, elem.getCenterLocation()));
			}
		}
		
		return maxDist;
	}

}
