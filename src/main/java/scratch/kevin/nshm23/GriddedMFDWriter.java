package scratch.kevin.nshm23;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultGridAssociations;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupMFDsModule;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;

import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.griddedSeismicity.FaultPolyMgr;

public class GriddedMFDWriter {

	public static void main(String[] args) throws IOException {
		File solDir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_02_02-nshm23_branches-WUS_FM_v3/");
		File baSolDir = new File(solDir, "node_branch_averaged");
		File baPairSolDir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2024_02_02-nshm23_branches-WUS_FM_v3-gridded_rebuild/node_branch_averaged_pairs");
		
		Region reg = NSHM23_RegionLoader.loadFullConterminousWUS();
		FaultSystemSolution refSol = FaultSystemSolution.load(new File(solDir, "results_WUS_FM_v3_branch_averaged_gridded.zip"));
		
//		File outputFile = new File("/tmp/nshm23_wus_gridded_mfds.csv");
//		FaultSystemSolution sol = refSol;
		
//		File outputFile = new File("/tmp/nshm23_wus_gridded_mfds_b0.csv");
//		FaultSystemSolution sol = FaultSystemSolution.load(new File(baSolDir, "SupraB_SupraB0.0.zip"));
		
//		File outputFile = new File("/tmp/nshm23_wus_gridded_mfds_b0.25.csv");
//		FaultSystemSolution sol = FaultSystemSolution.load(new File(baSolDir, "SupraB_SupraB0.25.zip"));
		
//		File outputFile = new File("/tmp/nshm23_wus_gridded_mfds_b0.5.csv");
//		FaultSystemSolution sol = FaultSystemSolution.load(new File(baSolDir, "SupraB_SupraB0.5.zip"));
		
//		File outputFile = new File("/tmp/nshm23_wus_gridded_mfds_b0.75.csv");
//		FaultSystemSolution sol = FaultSystemSolution.load(new File(baSolDir, "SupraB_SupraB0.75.zip"));
		
//		File outputFile = new File("/tmp/nshm23_wus_gridded_mfds_b1.csv");
//		FaultSystemSolution sol = FaultSystemSolution.load(new File(baSolDir, "SupraB_SupraB1.0.zip"));
		
//		File outputFile = new File("/tmp/nshm23_wus_gridded_mfds_classic.csv");
//		FaultSystemSolution sol = FaultSystemSolution.load(new File(baSolDir, "SegModel_Classic.zip"));
		
//		File outputFile = new File("/tmp/nshm23_wus_gridded_mfds_seg_low.csv");
//		FaultSystemSolution sol = FaultSystemSolution.load(new File(baSolDir, "SegModel_Low.zip"));
		
//		File outputFile = new File("/tmp/nshm23_wus_gridded_mfds_seg_middle.csv");
//		FaultSystemSolution sol = FaultSystemSolution.load(new File(baSolDir, "SegModel_Middle.zip"));
		
//		File outputFile = new File("/tmp/nshm23_wus_gridded_mfds_seg_high.csv");
//		FaultSystemSolution sol = FaultSystemSolution.load(new File(baSolDir, "SegModel_High.zip"));
		
//		File outputFile = new File("/tmp/nshm23_wus_gridded_mfds_seg_none.csv");
//		FaultSystemSolution sol = FaultSystemSolution.load(new File(baSolDir, "SegModel_None.zip"));
		
//		File outputFile = new File("/tmp/nshm23_wus_gridded_mfds_b0_seg_none.csv");
//		FaultSystemSolution sol = FaultSystemSolution.load(new File(baPairSolDir, "SupraB_SupraB0.0_SegModel_None.zip"));
		
		File outputFile = new File("/tmp/nshm23_wus_gridded_mfds_b1_seg_classic.csv");
		FaultSystemSolution sol = FaultSystemSolution.load(new File(baPairSolDir, "SupraB_SupraB1.0_SegModel_Classic.zip"));
		
//		File outputFile = new File("/tmp/ucerf3_gridded_mfds.csv");
//		refSol = null;
//		reg = new CaliforniaRegions.RELM_TESTING();
//		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2021_11_30-u3_branches-orig_calcs-5h/results_FM3_1_branch_averaged.zip"));
		
		GridSourceProvider gridSources = sol.getGridSourceProvider();
		
		GriddedRegion gridReg = new GriddedRegion(reg, 0.1, GriddedRegion.ANCHOR_0_0);
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(5.01, 8.55);
		
		IncrementalMagFreqDist[] mfds = new IncrementalMagFreqDist[gridReg.getNodeCount()];
		
		// do gridded
		for (int s=0; s<gridReg.getNodeCount(); s++) {
			int gridLocIndex = gridSources.getLocationIndex(gridReg.getLocation(s));
			IncrementalMagFreqDist gridMFD = gridSources.getMFD(gridLocIndex, 5.05);
			mfds[s] = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
			Preconditions.checkState((float)mfds[s].getMinX() == (float)gridMFD.getMinX(), "%s != %s", (float)mfds[s].getMinX(), (float)gridMFD.getMinX());
			for (int i=0; i<gridMFD.size(); i++)
				mfds[s].set(i, gridMFD.getY(i));
		}
		
		// add faults
		FaultSystemRupSet rupSet = sol.getRupSet();
		FaultGridAssociations assoc = rupSet.getModule(FaultGridAssociations.class);
		if (assoc == null) {
			if (refSol != null)
				assoc = refSol.getRupSet().requireModule(FaultGridAssociations.class);
			else if (outputFile.getName().contains("ucerf3"))
				assoc = FaultPolyMgr.loadSerializedUCERF3(FaultModels.FM3_1);
		}
		RupMFDsModule rupMFDs = sol.requireModule(RupMFDsModule.class);
		double totalMappedRate = 0d;
		for (int rupIndex=0; rupIndex<rupSet.getNumRuptures(); rupIndex++) {
			DiscretizedFunc rupMFD = rupMFDs.getRuptureMFD(rupIndex);
			if (rupMFD == null)
				rupMFD = new LightFixedXFunc(new double[] {rupSet.getMagForRup(rupIndex)}, new double[] {sol.getRateForRup(rupIndex)});
			List<Integer> rupSects = rupSet.getSectionsIndicesForRup(rupIndex);
			List<Double> sectAreas = new ArrayList<>(rupSects.size());
			double totArea = 0d;
			for (int s : rupSects) {
				double area = rupSet.getAreaForSection(s);
				sectAreas.add(area);
				totArea += area;
			}
			for (int s=0; s<rupSects.size(); s++) {
				int sectIndex = rupSects.get(s);
				double nuclFract = sectAreas.get(s)/totArea;
				Map<Integer, Double> sectAssoc = assoc.getNodeFractions(sectIndex);
				for (int gridIndex : sectAssoc.keySet()) {
					double fract = sectAssoc.get(gridIndex)*nuclFract;
					int mappedIndex = gridReg.indexForLocation(gridSources.getLocation(gridIndex));
					if (mappedIndex < 0)
						// outside of region
						continue;
					for (Point2D pt : rupMFD) {
						double fractRate = pt.getY()*fract;
						totalMappedRate += fractRate;
						if (pt.getX() < 5d)
							continue;
						mfds[mappedIndex].add(refMFD.getClosestXIndex(pt.getX()), fractRate);
					}
				}
			}
		}
		System.out.println("Mapped "+(float)totalMappedRate+"/"+(float)sol.getTotalRateForAllFaultSystemRups()
			+" = "+(float)(totalMappedRate/sol.getTotalRateForAllFaultSystemRups())+" fract of supra-seis in region");
		
		CSVFile<String> csv = new CSVFile<>(true);
		List<String> header = new ArrayList<>();
		header.add("Grid Index");
		header.add("Latitude");
		header.add("Longitude");
		for (int i=0; i<refMFD.size(); i++)
			header.add((float)refMFD.getX(i)+"");
		csv.addLine(header);
		
		SummedMagFreqDist summedMFD = new SummedMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		
		for (int s=0; s<mfds.length; s++) {
			List<String> line = new ArrayList<>(header.size());
			line.add(s+"");
			Location loc = gridReg.locationForIndex(s);
			line.add((float)loc.lat+"");
			line.add((float)loc.lon+"");
			for (int i=0; i<refMFD.size(); i++)
				line.add((float)mfds[s].getY(i)+"");
			summedMFD.addIncrementalMagFreqDist(mfds[s]);
			csv.addLine(line);
		}
		
		System.out.println(summedMFD);
		
		csv.writeToFile(outputFile);
	}

}
