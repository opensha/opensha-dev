package scratch.kevin.simulators.multiFault;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.FaultSubsectionCluster;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Jump;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.PlausibilityFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.PlausibilityResult;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RuptureConnectionSearch;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SectionDistanceAzimuthCalculator;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.stiffness.SubSectStiffnessCalculator;
import org.opensha.sha.simulators.stiffness.SubSectStiffnessCalculator.StiffnessType;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SubSectionMapping;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class PlausibilityEventIDCompare {

	public static void main(String[] args) throws IOException {
		RSQSimCatalog catalog = Catalogs.BRUCE_4983_STITCHED.instance();
		catalog.setFractForInclusion(0.5);
		
		int[] eventIDs = { 66240 };
		
//		SubSectStiffnessCalculator subSectCalc = new SubSectStiffnessCalculator(catalog.getU3SubSects(),
//				2d, 3e4, 3e4, 0.5);
//		subSectCalc.loadCacheFile(new File("/home/kevin/OpenSHA/UCERF4/rup_sets/"
//				+ "cff_cache_2606sects_2km_lambda30000_mu30000_coeff0.5.csv"), StiffnessType.CFF);
		
//		PlausibilityFilter filter = new ClusterPathCoulombCompatibilityFilter(
//				subSectCalc, StiffnessAggregationMethod.MEDIAN, 0f);
		PlausibilityFilter filter = null;
		
		List<RSQSimEvent> events = catalog.loader().byIDs(eventIDs);
		
		RuptureConnectionSearch connSearch = new RuptureConnectionSearch(
				null, new SectionDistanceAzimuthCalculator(catalog.getU3SubSects()), 100d, false);
		
		RSQSimSubSectionMapper subSectMapper = catalog.getSubSectMapper();
		
		for (RSQSimEvent event : events) {
			List<List<SubSectionMapping>> mappings = subSectMapper.getFilteredSubSectionMappings(event);
			List<FaultSection> subSects = new ArrayList<>();
			for (List<SubSectionMapping> parentMappings : mappings)
				for (SubSectionMapping mapping : parentMappings)
					subSects.add(mapping.getSubSect());
			System.out.println("================================");
			System.out.println("Mapped "+subSects.size()+" sections for rupture "
					+event.getID()+", M"+(float)event.getMagnitude());
			
			List<FaultSubsectionCluster> clusters = connSearch.calcClusters(subSects, false);
			List<Jump> jumps = connSearch.calcRuptureJumps(clusters, false);
			ClusterRupture rupture = connSearch.buildClusterRupture(clusters, jumps, false);
			System.out.println(rupture);
			System.out.println("--------------------------------");
			
			PlausibilityResult result = filter.apply(rupture, true);
			System.out.println("--------------------------------");
			System.out.println("result: "+result);
			System.out.println("================================");
		}
	}

}
