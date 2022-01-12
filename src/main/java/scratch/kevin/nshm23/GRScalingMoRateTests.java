package scratch.kevin.nshm23;

import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;

import scratch.UCERF3.enumTreeBranches.ScalingRelationships;

public class GRScalingMoRateTests {

	public static void main(String[] args) {
		
		double minMag = 6.5d;
		double rake = 180d;
		EvenlyDiscretizedFunc xVals = HistogramFunction.getEncompassingHistogram(minMag+0.0001, 8.5d, 0.1d);
		
		ScalingRelationships[] scales = {
				ScalingRelationships.ELLSWORTH_B,
				ScalingRelationships.ELLB_SQRT_LENGTH,
				ScalingRelationships.SHAW_2009_MOD,
				ScalingRelationships.SHAW_CONST_STRESS_DROP,
				ScalingRelationships.HANKS_BAKUN_08
		};
		
		double slipRate = 30d*1e-3;
		double length = 7.5d*1e3;
		double width = 14d*1e3;
		double area = length*width;
		
		double b = 0.8d;
		
		double moRate = FaultMomentCalc.getMoment(area, slipRate);
		
		System.out.println("Orig moment rate: "+moRate);
		GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(
				xVals.getMinX(), xVals.size(), xVals.getDelta(), moRate, b);
		System.out.println("GR tot rate="+(float)gr.calcSumOfY_Vals());
		
		for (ScalingRelationships scale : scales) {
			double estSlipRate = 0d;
			for (int i=0; i<gr.size(); i++) {
				double mag = gr.getX(i);
				double nuclRate = gr.getY(i);
				
				double rupArea = scale.getArea(mag, width);
				double rupWidth = width;
				double rupLength = rupArea/width;
				if (rupLength < rupWidth) {
					// make it a square
					rupLength = Math.sqrt(rupLength*rupWidth);
					rupWidth = rupLength;
				}
				double aveSlip = scale.getAveSlip(rupArea, rupLength, rupWidth, rake);
//				System.out.println("mag="+mag+"\trate="+rate+"\twidth="+rupWidth+"\tlen="+rupLength+"\tslip="+aveSlip);
				
				double particRate;
				if (rupArea < area)
					particRate = nuclRate;
				else
					particRate = nuclRate*rupArea/area;
				estSlipRate += aveSlip*particRate;
			}
			
			System.out.println();
			System.out.println(scale.getName());
			System.out.println("\tCalcSlipRate: "+(float)estSlipRate+" (vs "+(float)slipRate+")");
		}
	}

}
