package scratch.kevin.ucerf3.eal;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.dom4j.DocumentException;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sra.calc.parallel.MPJ_CondLossCalc;

import com.google.common.base.Preconditions;

import scratch.UCERF3.U3CompoundFaultSystemSolution;
import scratch.UCERF3.U3FaultSystemSolution;
import scratch.UCERF3.erf.mean.TrueMeanBuilder;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;
import scratch.UCERF3.utils.U3FaultSystemIO;

public class UCERF3_EAL_MagDistWriter {

	public static void main(String[] args) throws IOException, DocumentException {
		File invSolDir = new File("../opensha-ucerf3/src/scratch/UCERF3/data/scratch/InversionSolutions/");
		File trueMeanSolFile = new File(invSolDir,
				"2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_TRUE_HAZARD_MEAN_SOL_WITH_MAPPING.zip");
		
		File[] ealDirs = {
				new File("/home/kevin/OpenSHA/UCERF3/eal/2014_05_28-ucerf3-99percent-wills-smaller"),
				new File("/home/kevin/OpenSHA/UCERF3/eal/2016_06_06-ucerf3-90percent-wald")};
		
		Map<AttenRelRef, Double> imrWeightsMap = new HashMap<>();
		imrWeightsMap.put(AttenRelRef.CB_2014, 0.22);
		imrWeightsMap.put(AttenRelRef.CY_2014, 0.22);
		imrWeightsMap.put(AttenRelRef.ASK_2014, 0.22);
		imrWeightsMap.put(AttenRelRef.BSSA_2014, 0.22);
		imrWeightsMap.put(AttenRelRef.IDRISS_2014, 0.12);
		
		File outputDir = new File("/tmp");
		
		double thousandsToBillions = 1d/1e6; // portfolio units are in thousands (1e3), so convert to billions by dividing by 1e6
		
		double inflationScalar = 1d/0.9d;
		
		double lossScale = thousandsToBillions*inflationScalar;
		String lossUnits = "$ (Billions)";
		
		System.out.println("Loading true mean");
		U3FaultSystemSolution trueMeanSol = U3FaultSystemIO.loadSol(trueMeanSolFile);
		
		EvenlyDiscretizedFunc magDist = null;
		for (File ealDir : ealDirs) {
			System.out.println("Processing dir: "+ealDir.getName());
			for (AttenRelRef imr : imrWeightsMap.keySet()) {
				System.out.println("Processing GMPE: "+imr.name());
				double myWeight = imrWeightsMap.get(imr)/(double)ealDirs.length;
				
				File rupLossesFile = new File(ealDir, imr.name()+"_fss_index.bin");
				File rupGriddedFile = new File(ealDir, imr.name()+"_fss_gridded.bin");
				
				System.out.println("Loading losses");
				double[][] expectedLosses = MPJ_CondLossCalc.loadResults(rupLossesFile);
				DiscretizedFunc[] griddedFuncs = null;
				if (rupGriddedFile != null) {
					griddedFuncs = MPJ_CondLossCalc.loadGridSourcesFile(rupGriddedFile,
							trueMeanSol.getGridSourceProvider().getGriddedRegion());
					double totCond = 0;
					int gridNonNull = 0;
					for (DiscretizedFunc func : griddedFuncs) {
						if (func != null) {
							gridNonNull++;
							for (Point2D pt : func)
								totCond += pt.getY();
						}
					}
					System.out.println("Tot grid conditional "+totCond+" ("+gridNonNull+" non null)");
				}
				
				EvenlyDiscretizedFunc dist = UCERF3_EAL_Combiner.getTotalEALMagDist(
						trueMeanSol, expectedLosses, griddedFuncs, lossScale);
				if (magDist == null)
					magDist = new EvenlyDiscretizedFunc(dist.getMinX(), dist.getMaxX(), dist.size());
				else
					Preconditions.checkState(magDist.size() == dist.size() && magDist.getMinX() == dist.getMinX());
				for (int i=0; i<dist.size(); i++)
					magDist.add(i, dist.getY(i)*myWeight);
			}
		}
		
		// write it out
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("Magnitude", "EAL, "+lossUnits);
		for (int i=0; i<magDist.size(); i++)
			csv.addLine((float)magDist.getX(i)+"", (float)magDist.getY(i)+"");
		csv.writeToFile(new File(outputDir, "mag_dist_losses.csv"));
		
		System.out.println("Sum y: "+(float)magDist.calcSumOfY_Vals());
	}

}
