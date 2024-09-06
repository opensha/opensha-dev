package scratch.kevin.prvi25;

import java.io.File;
import java.io.IOException;

import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.eq.MagUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.modules.AveSlipModule;

public class SolAreaValidate {

	public static void main(String[] args) throws IOException {
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(new File(args[0]));
		AveSlipModule slips = rupSet.getModule(AveSlipModule.class);
		
		int worstSectMismatch = -1;
		double worstSectDiff = 0d;
		double worstSectSumArea = 0d;
		
		int worstMagMismatch = -1;
		double worstMagDiff = 0d;
		double worstMomentToMag = 0d;
		for (int r=0; r<rupSet.getNumRuptures(); r++) {
			double rupArea = rupSet.getAreaForRup(r);
			double sumSectArea = 0d;
			for (int s : rupSet.getSectionsIndicesForRup(r))
				sumSectArea += rupSet.getAreaForSection(s);
			
			if (slips != null) {
				double mag = rupSet.getMagForRup(r);
				double slip = slips.getAveSlip(r);
//				double moment = MagUtils.magToMoment(mag);
				double moment = FaultMomentCalc.getMoment(rupArea, slip);
				double momentToMag = MagUtils.momentToMag(moment);
				double magDiff = Math.abs(momentToMag - mag);
				if (magDiff > worstMagDiff) {
					worstMagDiff = magDiff;
					worstMagMismatch = r;
					worstMomentToMag = momentToMag;
				}
			}
			
			double sectDiff = Math.abs(sumSectArea - rupArea);
			if (sectDiff > worstSectDiff) {
				worstSectMismatch = r;
				worstSectDiff = sectDiff;
				worstSectSumArea = sumSectArea;
			}
		}
		
		if (worstSectMismatch >= 0) {
			System.out.println("Worst sum[sectAreas] vs rupArea mismatch: rupture "+worstSectMismatch);
			System.out.println("\tsum[sectAreas]:\t"+(float)worstSectSumArea);
			System.out.println("\trupArea:\t"+(float)rupSet.getAreaForRup(worstSectMismatch));
			System.out.println("\tdiff:\t"+(float)(rupSet.getAreaForRup(worstSectMismatch) - worstSectSumArea));
		}
		if (worstMagMismatch >= 0) {
			System.out.println("Worst mag vs area&slip->moment->mag: rupture "+worstMagMismatch);
			System.out.println("\tmag:\t"+(float)rupSet.getMagForRup(worstMagMismatch));
			System.out.println("\tmomentToMag:\t"+(float)worstMomentToMag);
			System.out.println("\trupArea:\t"+(float)rupSet.getAreaForRup(worstMagMismatch));
		}
	}

}
