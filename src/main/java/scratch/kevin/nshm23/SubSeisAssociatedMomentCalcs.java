package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Map;

import org.opensha.commons.geo.Region;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultGridAssociations;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.ModelRegion;
import org.opensha.sha.earthquake.faultSysSolution.modules.SectSlipRates;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.gridded.NSHM23_FaultCubeAssociations;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

public class SubSeisAssociatedMomentCalcs {

	public static void main(String[] args) throws IOException {
		File baSolFile = new File("/data/kevin/nshm23/batch_inversions/"
				+ "2023_04_11-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip");
		
		FaultSystemSolution sol = FaultSystemSolution.load(baSolFile);
		FaultSystemRupSet rupSet = sol.getRupSet();
		
//		testCubeAssoc(rupSet);
//		System.exit(0);
		
		GridSourceProvider gridProv = sol.getGridSourceProvider();
		FaultGridAssociations assoc = rupSet.requireModule(FaultGridAssociations.class);
		
		double totAssocSubSeisMo = 0d;
		double[] sectAssocSubSeisMos = new double[rupSet.getNumSections()];
		
		for (int i=0; i<gridProv.size(); i++) {
			Map<Integer, Double> sectFracts = assoc.getScaledSectFracsOnNode(i);
			if (sectFracts != null && !sectFracts.isEmpty()) {
				double sumFract = 0d;
				for (double fract : sectFracts.values())
					sumFract += fract;
				if ((float)sumFract > 1f)
					System.err.println("WARNING: sumFract="+(float)sumFract+" for node "
							+i+" with "+sectFracts.size()+" sections");
				Preconditions.checkState((float)sumFract <= 1f, "Unexpected sumFract=%s for node %s",
						(float)sumFract+"", i);
				
				IncrementalMagFreqDist mfd = gridProv.getMFD_SubSeisOnFault(i);
				Preconditions.checkNotNull(mfd);
				
				double mfdMoment = mfd.getTotalMomentRate();
				totAssocSubSeisMo += mfdMoment;
				
				for (int sectIndex : sectFracts.keySet()) {
					double fract = sectFracts.get(sectIndex);
					sectAssocSubSeisMos[sectIndex] += fract*mfdMoment/sumFract;
				}
			}
		}
		
		MinMaxAveTracker fractAssocSubSiesMoTrack = new MinMaxAveTracker();
		double[] slipBinEdges = { 1d, 5d, 10d, 20d };
		MinMaxAveTracker[] slipBinTracks = new MinMaxAveTracker[slipBinEdges.length+1];
		for (int i=0; i<slipBinTracks.length; i++)
			slipBinTracks[i] = new MinMaxAveTracker();
		double totDefModMo = 0d;
		
		SectSlipRates slips = rupSet.getSectSlipRates();
		
		int numAbove = 0;
		
		for (int s=0; s<sectAssocSubSeisMos.length; s++) {
			double sectMo = slips.calcMomentRate(s);
			totDefModMo += sectMo;
			double fract = sectAssocSubSeisMos[s]/sectMo;
			System.out.println(s+". "+rupSet.getFaultSectionData(s).getSectionName()+":\t"+(float)fract);
			fractAssocSubSiesMoTrack.addValue(fract);
			
			double slipRate = slips.getSlipRate(s)*1e3;
			int binIndex = 0;
			for (binIndex=0; binIndex<slipBinEdges.length; binIndex++) {
				if (slipRate < slipBinEdges[binIndex])
					break;
			}
			slipBinTracks[binIndex].addValue(fract);
			
			if (fract > 1)
				numAbove++;
		}
		
		System.out.println("Total Sub-Seis Moment: "+(float)totAssocSubSeisMo+" / "+(float)totDefModMo
				+" = "+pDF.format(totAssocSubSeisMo/totDefModMo));
		System.out.println("Sect Sub-Seis Stats: "+trackStr(fractAssocSubSiesMoTrack));
		System.out.println("Exceedances: "+numAbove+"/"+rupSet.getNumSections()+" ("
				+pDF.format((double)numAbove/(double)rupSet.getNumSections())+")");
		System.out.println("Slip Rate Binned Sect Stats");
		for (int i=0; i<slipBinTracks.length; i++) {
			String slipStr;
			if (i == 0)
				slipStr = "slip < "+slipDF.format(slipBinEdges[i])+" mm/yr";
			else if (i<slipBinEdges.length)
				slipStr = slipDF.format(slipBinEdges[i-1])+" < slip < "+slipDF.format(slipBinEdges[i])+" mm/yr";
			else
				slipStr = "slip > "+slipDF.format(slipBinEdges[slipBinEdges.length-1])+" mm/yr";
			System.out.println("\t"+slipStr+":\t"+trackStr(slipBinTracks[i]));
		}
	}
	
	private static final DecimalFormat slipDF = new DecimalFormat("0.#");
	private static final DecimalFormat pDF = new DecimalFormat("0.00%");
	
	private static String trackStr(MinMaxAveTracker track) {
		return "avg: "+pDF.format(track.getAverage())+"; range=["+pDF.format(track.getMin())+", "+pDF.format(track.getMax())+"]";
	}
	
	private static void testCubeAssoc(FaultSystemRupSet rupSet) throws IOException {
		FaultGridAssociations storedAssoc = rupSet.requireModule(FaultGridAssociations.class);
		
		Region region = rupSet.requireModule(ModelRegion.class).getRegion();
		NSHM23_FaultCubeAssociations calcAssoc = NSHM23_InvConfigFactory.buildFaultCubeAssociations(rupSet, null, region);
		
		Preconditions.checkState(calcAssoc.getRegion().equalsRegion(storedAssoc.getRegion()));
		int nodeCount = calcAssoc.getRegion().getNodeCount();
		
		for (int n=0; n<nodeCount; n++) {
			Map<Integer, Double> storedFracts = storedAssoc.getScaledNodeFractions(n);
			Map<Integer, Double> calcFracts = calcAssoc.getScaledNodeFractions(n);
			
			Preconditions.checkState((storedFracts == null) == (calcFracts == null));
			
			if (storedFracts != null) {
				if (!storedFracts.keySet().equals(calcFracts.keySet())) {
					// can happen if aseismicity varies on sections affecting this node
					System.err.println("WARNING: keyset mismatch for node "+n+", skipping");
					continue;
				}
				double storedSum = 0d;
				double calcSum = 0d;
				for (int sectIndex : storedFracts.keySet()) {
					double storedFract = storedFracts.get(sectIndex);
					double calcFract = calcFracts.get(sectIndex);
					
					storedSum += storedFract;
					calcSum += calcFract;
					
//					if ((float)storedFract != (float)calcFract)
//						System.err.println("MISMATCH for node "+n+", sect "+sectIndex
//									+":\tcalc="+(float)calcFract+"\tstored="+(float)storedFract);
				}
				if ((float)storedSum != (float)calcSum)
					System.err.println("MISMATCHED sum for node "+n
								+":\tcalc="+(float)calcSum+"\tstored="+(float)storedSum);
			}
		}
	}

}
