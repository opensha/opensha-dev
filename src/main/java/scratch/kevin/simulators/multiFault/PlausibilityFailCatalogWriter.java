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
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.iden.AbstractRuptureIdentifier;
import org.opensha.sha.simulators.parsers.RSQSimFileWriter;
import org.opensha.sha.simulators.stiffness.SubSectStiffnessCalculator;
import org.opensha.sha.simulators.stiffness.SubSectStiffnessCalculator.StiffnessType;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper.SubSectionMapping;

import com.google.common.base.Preconditions;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.RSQSimCatalog.Loader;

public class PlausibilityFailCatalogWriter {

	public static void main(String[] args) throws IOException {
		RSQSimCatalog catalog = Catalogs.BRUCE_4983_STITCHED.instance();
		catalog.setFractForInclusion(0.5);
		
//		SubSectStiffnessCalculator subSectCalc = new SubSectStiffnessCalculator(catalog.getU3SubSects(),
//				2d, 3e4, 3e4, 0.5);
//		subSectCalc.loadCacheFile(new File("/home/kevin/OpenSHA/UCERF4/rup_sets/"
//				+ "cff_cache_2606sects_2km_lambda30000_mu30000_coeff0.5.csv"), StiffnessType.CFF);
//		
//		PlausibilityFilter filter = new ClusterPathCoulombCompatibilityFilter(
//				subSectCalc, StiffnessAggregationMethod.MEDIAN, 0f);
		PlausibilityFilter filter = null;
		
		File outputDir = new File(catalog.getCatalogDir(), "plausibility_filtered");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		String prefix = "filtered";
		
		RuptureConnectionSearch connSearch = new RuptureConnectionSearch(
				null, new SectionDistanceAzimuthCalculator(catalog.getU3SubSects()), 100d, false);
		
		FilterIden iden = new FilterIden(filter, catalog.getSubSectMapper(), true, connSearch);
		
		Loader loader = catalog.loader().minMag(6.5d);
		
		RSQSimFileWriter.writeFilteredCatalog(loader.iterable(), iden, outputDir, prefix,
				false, catalog.getTransitions());
	}
	
	private static class FilterIden extends AbstractRuptureIdentifier {
		
		private PlausibilityFilter filter;
		private RSQSimSubSectionMapper subSectMapper;
		private boolean failures;
		private RuptureConnectionSearch connSearch;

		public FilterIden(PlausibilityFilter filter, RSQSimSubSectionMapper subSectMapper, boolean failures,
				RuptureConnectionSearch connSearch) {
			this.filter = filter;
			this.subSectMapper = subSectMapper;
			this.failures = failures;
			this.connSearch = connSearch;
		}

		@Override
		public boolean isMatch(SimulatorEvent event) {
			List<List<SubSectionMapping>> mappings = subSectMapper.getFilteredSubSectionMappings(event);
			List<FaultSection> subSects = new ArrayList<>();
			for (List<SubSectionMapping> parentMappings : mappings)
				for (SubSectionMapping mapping : parentMappings)
					subSects.add(mapping.getSubSect());
			
			if (subSects.size() < 2)
				return false;
			
			List<FaultSubsectionCluster> clusters = connSearch.calcClusters(subSects, false);
			List<Jump> jumps = connSearch.calcRuptureJumps(clusters, false);
			ClusterRupture rupture = connSearch.buildClusterRupture(clusters, jumps, false);
			
			PlausibilityResult result = filter.apply(rupture, false);
			return failures && !result.isPass() || !failures && result.isPass();
		}

		@Override
		public String getName() {
			return filter.getName();
		}
		
	}

}
