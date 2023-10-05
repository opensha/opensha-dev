package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SingleStates;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

public class RayOffORCoastCalcs {

	public static void main(String[] args) throws IOException {
		Region rect = new Region(new Location(42, -126.5), new Location(46.26, -123));
		Region or = NSHM23_SingleStates.OR.loadRegion();
		
		Region[] offshores = Region.subtract(rect, or);
		Region offshore = null;
		for (Region reg : offshores)
			if (offshore == null || reg.getExtent() > offshore.getExtent())
				offshore = reg;
		
		Feature.write(offshore.toFeature(), new File("/tmp/offshore.geojson"));
		
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2023_06_23-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip"));
		
		IncrementalMagFreqDist refMFD = FaultSysTools.initEmptyMFD(6.05d, 7.95d);
		
		IncrementalMagFreqDist solNucl = sol.calcNucleationMFD_forRegion(
				offshore, refMFD.getMinX(), refMFD.getMaxX(), refMFD.getDelta(), false);
		FaultSystemRupSet rupSet = sol.getRupSet();
		double[] fractRups = rupSet.getFractRupsInsideRegion(offshore, false);
		IncrementalMagFreqDist solPartic = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		for (int r=0; r<fractRups.length; r++)
			if (fractRups[r] > 0d)
				solPartic.add(solPartic.getClosestXIndex(rupSet.getMagForRup(r)), sol.getRateForRup(r));
		
		IncrementalMagFreqDist griddedNucl = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		
		GridSourceProvider gridProv = sol.getGridSourceProvider();
		GriddedRegion gridReg = gridProv.getGriddedRegion();
		
		float minMag = (float)(refMFD.getMinX() - 0.5*refMFD.getDelta());
		float maxMag = (float)(refMFD.getMaxX() + 0.5*refMFD.getDelta());
		for (int i=0; i<gridReg.getNodeCount(); i++) {
			if (offshore.contains(gridReg.getLocation(i))) {
				IncrementalMagFreqDist gridMFD = gridProv.getMFD(i);
				for (int j=0; j<gridMFD.size(); j++) {
					double mag = gridMFD.getX(j);
					if ((float)mag >= minMag && (float)mag <= maxMag)
						griddedNucl.add(griddedNucl.getClosestXIndex(mag), gridMFD.getY(j));
				}
			}
		}
		
		File outputDir = new File("/tmp");
		String prefix = "offshore_or_rates";
		for (boolean cumulative : new boolean[] {false,true}) {
			CSVFile<String> csv = new CSVFile<>(true);
			
			csv.addLine("Magnitude", "Fault Nucleation Rate", "Fault Participation Rate", "Gridded Rate");
			
			EvenlyDiscretizedFunc faultNucl, faultPartic, gridded;
			if (cumulative) {
				faultNucl = solNucl.getCumRateDistWithOffset();
				faultPartic = solPartic.getCumRateDistWithOffset();
				gridded = griddedNucl.getCumRateDistWithOffset();
			} else {
				faultNucl = solNucl;
				faultPartic = solPartic;
				gridded = griddedNucl;
			}
			
			for (int i=0; i<faultNucl.size(); i++) {
				List<String> line = new ArrayList<>();
				if (cumulative)
					line.add("M≥"+(float)faultNucl.getX(i));
				else
					line.add("M"+(float)faultNucl.getX(i));
				line.add((float)faultNucl.getY(i)+"");
				line.add((float)faultPartic.getY(i)+"");
				line.add((float)gridded.getY(i)+"");
				csv.addLine(line);
			}
			File outputFile;
			if (cumulative)
				outputFile = new File(outputDir, prefix+"_cml.csv");
			else
				outputFile = new File(outputDir, prefix+".csv");
			csv.writeToFile(outputFile);
		}
		
	}

}
