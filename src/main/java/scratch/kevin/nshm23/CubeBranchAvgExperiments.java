package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultCubeAssociations;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultCubeAssociations.CubeToGridNodeAggregator;

public class CubeBranchAvgExperiments {

	public static void main(String[] args) throws IOException {
		FaultSystemSolution baSol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2023_04_11-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip"));
		FaultSystemRupSet rupSet = baSol.getRupSet();
		
		FaultCubeAssociations baCubeAssoc = rupSet.requireModule(FaultCubeAssociations.class);
		
		// now aggregate ba gridded to cubes
		CubeToGridNodeAggregator agg = new CubeToGridNodeAggregator(baCubeAssoc);
		
		// are they the same?
		DiffTrack nodeFractTrack = new DiffTrack();
		DiffTrack sectOnNodeFractTrack = new DiffTrack();
		DiffTrack scaledSectOnNodeFractTrack = new DiffTrack();
		for (int n=0; n<agg.getRegion().getNodeCount(); n++) {
			double nodeFrac1 = agg.getNodeFraction(n);
			double nodeFrac2 = baCubeAssoc.getNodeFraction(n);
			nodeFractTrack.add(nodeFrac1, nodeFrac2);
			
			Map<Integer, Double> sectFracts1 = agg.getSectionFracsOnNode(n);
			Map<Integer, Double> sectFracts2 = baCubeAssoc.getSectionFracsOnNode(n);
			
			if (!sectFracts1.keySet().equals(sectFracts2.keySet()))
				System.err.println("WARNING: section mapping mismatch for cube "+n+", only testing union");
			for (int sectIndex : sectFracts1.keySet()) {
				if (sectFracts2.containsKey(sectIndex))
					sectOnNodeFractTrack.add(sectFracts1.get(sectIndex), sectFracts2.get(sectIndex));
			}
			sectFracts1 = agg.getScaledSectFracsOnNode(n);
			sectFracts2 = baCubeAssoc.getScaledSectFracsOnNode(n);
			for (int sectIndex : sectFracts1.keySet()) {
				if (sectFracts2.containsKey(sectIndex))
					scaledSectOnNodeFractTrack.add(sectFracts1.get(sectIndex), sectFracts2.get(sectIndex));
			}
		}
		
		System.out.println("Node fractions\n"+nodeFractTrack);
		System.out.println("Sect on node fractions\n"+sectOnNodeFractTrack);
		System.out.println("Scaled sect on node fractions\n"+scaledSectOnNodeFractTrack);
	}
	
	private static class DiffTrack {
		private int num;
		private double sumAbsDiff;
		private double sumAbsPDiff;
		private double maxAbsDiff;
		private double maxAbsPDiff;
		
		public void add(double val1, double val2) {
			double absDiff = Math.abs(val1 - val2);
			double absPDiff;
			if ((float)val2 == 0f && (float)val1 == 0f)
				absPDiff = 0;
			else
				absPDiff = 100d*absDiff/val2;
			
			num++;
			sumAbsDiff += absDiff;
			sumAbsPDiff += absPDiff;
			maxAbsDiff = Math.max(absDiff, maxAbsDiff);
			maxAbsPDiff = Math.max(absPDiff, maxAbsPDiff);
		}
		
		@Override
		public String toString() {
			double avgDiff = sumAbsDiff/(double)num;
			double avgPDiff = sumAbsPDiff/(double)num;
			return "\t|diff|: agv="+(float)avgDiff+", max="+(float)maxAbsDiff
					+"\n\t|%diff|: agv="+(float)avgPDiff+" %, max="+(float)maxAbsPDiff+" %";
		}
	}

}
