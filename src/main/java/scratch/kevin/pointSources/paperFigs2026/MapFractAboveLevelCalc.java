package scratch.kevin.pointSources.paperFigs2026;

import static scratch.kevin.pointSources.paperFigs2026.ConstantsAndSettings.*;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;

public class MapFractAboveLevelCalc {

	public static void main(String[] args) throws ZipException, IOException {
		ReturnPeriods rp = ReturnPeriods.TWO_IN_50;
		double[] thresholds = {0.1, 0.25};
		
		String mapZipName = "results_hazard_INCLUDE.zip";
		
		String entryName = "map_pga_"+rp.name()+".txt";
		
		File mapFile = new File(Models.AS_PUBLISHED.getMapDir(), mapZipName);
		ZipFile zip = new ZipFile(mapFile);
		
		GriddedGeoDataSet mask = buildLandMask(FULL_GRID_REG);
		GriddedGeoDataSet map = HazardMapFigures.loadXYZ(zip, FULL_GRID_REG, entryName, mask);
		int[] numBelows = new int[thresholds.length];
		int totNum = 0;
		for (int i=0; i<map.size(); i++) {
			double val = map.get(i);
			if (Double.isFinite(val)) {
				totNum++;
				for (int t=0; t<thresholds.length; t++)
					if (val < thresholds[t])
						numBelows[t]++;
			}
		}
		DecimalFormat pDF = new DecimalFormat("0.0%");
		for (int t=0; t<thresholds.length; t++)
			System.out.println(numBelows[t]+"/"+totNum+" ("+pDF.format((double)numBelows[t]/(double)totNum)+") < "+(float)thresholds[t]+" g");
	}

}
