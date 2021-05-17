package scratch.kevin.ucerf3.etas;

import org.opensha.commons.util.Interpolate;

public class MFD_GR_Extrap {

	public static void main(String[] args) {
		// Bombay Beach M4.8
//		double rateM3 = 5.98355;
//		double rateM4 = 0.60505;
//		double rateM7 = 0.00168;
		
		// Bombay Beach M4.8
		double rateM3 = 3.2597587;
		double rateM4 = 0.32895;
		double rateM7 = 7.375E-5;
		
		double targetMag = 7d;
		double extrap = Interpolate.findLogY(new double[] {3d, 4d}, new double[] {rateM3, rateM4}, targetMag);
		System.out.println("Extrapolated: "+extrap);
		double f7 = rateM7 / extrap;
		
		System.out.println("F7 = "+rateM7+" / "+extrap+" = "+f7);
	}

}
