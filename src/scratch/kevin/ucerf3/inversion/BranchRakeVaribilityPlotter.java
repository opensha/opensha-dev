package scratch.kevin.ucerf3.inversion;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.primitives.Doubles;

import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.utils.DeformationModelFileParser;
import scratch.UCERF3.utils.DeformationModelFileParser.DeformationSection;

public class BranchRakeVaribilityPlotter {

	/**
	 * @param args
	 * @throws IOException 
	 * @throws RuntimeException 
	 * @throws GMT_MapException 
	 */
	public static void main(String[] args) throws IOException, GMT_MapException, RuntimeException {
		List<Double> dmWeights = Lists.newArrayList();
		for (DeformationModels dm : DeformationModels.values()) 
			if (dm.getRelativeWeight(null) > 0)
				dmWeights.add(dm.getRelativeWeight(null));
		
		Region region = new CaliforniaRegions.RELM_TESTING();
		
		File saveDir = new File("/tmp");
		
		for (FaultModels fm : FaultModels.values()) {
			if (fm == FaultModels.FM2_1)
				continue;
			List<FaultSectionPrefData> sects = fm.fetchFaultSections();
			
			List<Map<Integer, DeformationSection>> dmSectsMaps = Lists.newArrayList();
			
			Map<Integer, DeformationSection> geolSects = null;
			
			for (DeformationModels dm : DeformationModels.values()) {
				if (dm.getRelativeWeight(null) > 0) {
					Map<Integer, DeformationSection> dmSects = DeformationModelFileParser.load(dm.getDataFileURL(fm));
					dmSectsMaps.add(dmSects);
					if (dm == DeformationModels.GEOLOGIC)
						geolSects = dmSects;
				}
			}
			Preconditions.checkNotNull(geolSects);
			
			List<LocationList> miniFaults = Lists.newArrayList();
			List<Double> values = Lists.newArrayList();
			
			double maxMaxDiff = 0;
			
			for (FaultSectionPrefData sect : sects) {
				Integer sectID = sect.getSectionId();
				
				int minis = sect.getFaultTrace().size()-1;
				
				for (int mini=0; mini<minis; mini++) {
					double geolRake = geolSects.get(sectID).getRakes().get(mini);
					double maxDiff = 0d;
					for (Map<Integer, DeformationSection> dmSects : dmSectsMaps) {
						double diff = Math.abs(dmSects.get(sectID).getRakes().get(mini) - geolRake);
						if (diff > 180)
							diff = 360-diff;
						if (diff > maxDiff)
							maxDiff = diff;
					}
					
					if (maxDiff > maxMaxDiff)
						maxMaxDiff = maxDiff;
					
					Location loc1 = sect.getFaultTrace().get(mini);
					Location loc2 = sect.getFaultTrace().get(mini+1);
					
					LocationList miniFault = new LocationList();
					miniFault.add(loc1);
					miniFault.add(loc2);
					miniFaults.add(miniFault);
					values.add(maxDiff);
				}
			}
			
			CPT cpt = FaultBasedMapGen.getSlipRateCPT().rescale(0d, maxMaxDiff);
			String prefix = fm.getShortName()+"_rake_max_diff";
			String label = fm.getName()+" Max Rake Diff from Geologic";
			
			double[] valArray = Doubles.toArray(values);
			
			FaultBasedMapGen.makeFaultPlot(cpt, miniFaults, valArray, region, saveDir, prefix, false, false, label);
		}
	}

}
