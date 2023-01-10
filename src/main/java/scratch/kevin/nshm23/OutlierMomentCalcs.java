package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SingleStates;
import org.opensha.sha.faultSurface.FaultSection;

public class OutlierMomentCalcs {

	public static void main(String[] args) throws IOException {
		CSVFile<String> csv = new CSVFile<>(false);
		
//		NSHM23_DeformationModels.OUTLIER_SUB_USE_BOUND = false;
//		File outputFile = new File("/tmp/dm_outlier_substitution_moments.csv");
		NSHM23_DeformationModels.OUTLIER_SUB_USE_BOUND = true;
		File outputFile = new File("/tmp/dm_outlier_bounded_substitution_moments.csv");
		
		double[] ycs = { 2d, 3.5d, 5d };
		
		List<String> header = new ArrayList<>();
		
		header.add("Deformation Model");
		header.add("State");
		header.add("Original Moment Rate (N-m/yr)");
		for (double yc : ycs) {
			header.add("Yc="+(float)yc+" Moment Rate (N-m/yr)");
			header.add("% Change");
		}
		
		csv.addLine(header);
		
		List<String> stateNames = new ArrayList<>();
		List<Region> stateRegions = new ArrayList<>();
		stateNames.add("Full Model");
		stateRegions.add(null);
		for (NSHM23_SingleStates state : NSHM23_SingleStates.values()) {
			stateNames.add(state.getName());
			stateRegions.add(state.loadRegion());
		}
		
		NSHM23_FaultModels fm = NSHM23_FaultModels.NSHM23_v2;
		
//		DecimalFormat pDF = new DecimalFormat("0.00%");
		DecimalFormat eDF = new DecimalFormat("0.000E0");
		DecimalFormat twoDigitsDF = new DecimalFormat("0.00");
		
		for (NSHM23_DeformationModels dm : NSHM23_DeformationModels.values()) {
			if (dm.getNodeWeight(null) == 0d)
				continue;
			
			NSHM23_DeformationModels.OUTLIER_SUB_YC = null;
			List<? extends FaultSection> origSubSects = dm.build(fm);
			
			List<List<? extends FaultSection>> modSubSectsList = new ArrayList<>();
			for (double yc : ycs) {
				NSHM23_DeformationModels.OUTLIER_SUB_YC = yc;
				modSubSectsList.add(dm.build(fm));
			}
			
			for (int s=0; s<stateNames.size(); s++) {
				List<String> line = new ArrayList<>();
				
				line.add(dm.getShortName());
				line.add(stateNames.get(s));
				
				Region region = stateRegions.get(s);
				
				double origMoment = calcMomentRate(origSubSects, region);
				
				line.add(eDF.format(origMoment));
				
				for (List<? extends FaultSection> modSubSects : modSubSectsList) {
					double modMoment = calcMomentRate(modSubSects, region);;
					line.add(eDF.format(modMoment));
					line.add(twoDigitsDF.format(100d*(modMoment-origMoment)/origMoment));
				}
				
				csv.addLine(line);
			}
		}
		
		csv.writeToFile(outputFile);
	}
	
	private static double calcMomentRate(List<? extends FaultSection> sects, Region region) {
		double mo = 0d;
		
		for (FaultSection sect : sects) {
			if (region != null) {
				boolean inside = false;
				for (Location loc : sect.getFaultTrace()) {
					if (region.contains(loc)) {
						inside = true;
						break;
					}
				}
				if (!inside)
					continue;
			}
			mo += sect.calcMomentRate(false);
		}
		
		return mo;
	}

}
