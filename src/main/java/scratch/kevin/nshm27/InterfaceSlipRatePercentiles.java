package scratch.kevin.nshm27;

import java.io.IOException;
import java.text.DecimalFormat;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.sha.earthquake.faultSysSolution.logicTree.dmSampling.DeformationModelDistSampler.AverageSampler;
import org.opensha.sha.earthquake.faultSysSolution.logicTree.dmSampling.DeformationModelDistSampler.FixedFractileSampler;
import org.opensha.sha.earthquake.faultSysSolution.logicTree.dmSampling.DeformationModelDistSampler.FixedSampler;

import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_InterfaceDeformationModels;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_InterfaceDeformationModels.DeformationFront;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_InterfaceFaultModels;

public class InterfaceSlipRatePercentiles {

	public static void main(String[] args) throws IOException {
		double[] fractiles = {
				0.00000001,
				0.0000001,
				0.000001,
				0.00001,
				0.0001,
				0.001,
				0.01,
				0.025,
				0.16,
				0.5,
				0.84,
				0.975,
				0.99,
				0.999,
				0.9999,
				0.99999,
				0.999999,
				0.9999999,
				0.99999999,
				1d,
				2d
		};
		
		DecimalFormat pDF = new DecimalFormat("0.#########%");
		DecimalFormat slipDF = new DecimalFormat("0.###");
		DecimalFormat groupedDF = new DecimalFormat("0");
		groupedDF.setGroupingSize(3);
		groupedDF.setGroupingUsed(true);
		
		for (NSHM27_InterfaceFaultModels fm : NSHM27_InterfaceFaultModels.values()) {
			DeformationFront df = NSHM27_InterfaceDeformationModels.getDeformationFront(fm);
			System.out.println(fm.getName());
			
			String[] labels = new String[fractiles.length];
			String[] results = new String[fractiles.length];
			int longestLabel = 0;
			for (int f=0; f < fractiles.length; f++) {
				double fractile = fractiles[f];
				String label;
				FixedSampler sampler;
				if (fractile == 2d) {
					label = "Distribution-average";
					sampler = new AverageSampler();
				} else {
					double rarity = Math.min(1d-fractile, fractile);
					
					label = pDF.format(fractile);
					if (rarity < 0.1)
						label += " (1 in "+groupedDF.format(1d/rarity)+")";
					sampler = new FixedFractileSampler(fractile);
				}
				
				double[] slipRates = NSHM27_InterfaceDeformationModels.getCoupledSlipRates(df, sampler);
				
				longestLabel = Integer.max(longestLabel, label.length());
				labels[f] = label+":";
				results[f] = "mean="+slipDF.format(StatUtils.mean(slipRates))
					+"\trange=["+slipDF.format(StatUtils.min(slipRates))+", "+slipDF.format(StatUtils.max(slipRates))+"]";
			}
			for (int f=0; f<fractiles.length; f++)
				while (labels[f].length() < longestLabel)
					labels[f] += " ";
			for (int f=0; f<fractiles.length; f++)
				System.out.println("\t"+labels[f]+"\t"+results[f]);
		}
	}

}
