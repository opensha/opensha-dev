package scratch.ned;

import org.opensha.sha.earthquake.calc.recurInterval.LognormalDistCalc;

public class ConfFromCOV_ForLognormalDist {

	public static void main(String[] args) {
		
		LognormalDistCalc dist = new LognormalDistCalc();
		double mean = 2.0;
		double[] confArray = {0.68,0.95, 0.99};
		double[] last_fuArray = {0,0,0};
		
		System.out.println("COV\tf68%\tf95%\tf99%\tP≥10%");
		for(double cov=0.05;cov<1.41;cov+=0.05) {
			dist.setAll(mean, cov, mean/1000, (int)(mean*10000));
			String tableLine = String.format("%.2f", cov); 
			
			for(int i=0; i<confArray.length; i++) {
				double conf = confArray[i];
				double fu = dist.getFractionalUncertaintyForConfBounds(conf);
				if(Double.isNaN(last_fuArray[i])|| fu<last_fuArray[i])
					fu = Double.NaN;
				tableLine += "\t"+String.format("%.2f", fu);
				last_fuArray[i] = fu;
			}
			double fr10perc = 1.0 - (dist.getCDF().getInterpolatedY(mean*1.1) - dist.getCDF().getInterpolatedY(mean/1.1));
			tableLine += "\t"+String.format("%.2f", fr10perc);
			System.out.println(tableLine);
		}			
	}
}
