package scratch.peter.tmp;

import org.opensha.commons.eq.MagUtils;
import org.opensha.nshmp2.util.NSHMP_Utils;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

/**
 * Add comments here
 *
 * 
 * @author Peter Powers
 * @version $Id:$
 */
public class t3 {

	public static void main(String[] args) {
		double mMin = 5.3;
		double mMax = 7.2;
		int mNum = 20;
		double aVal = 0.1556;
		double bVal = 1.23;
		
		IncrementalMagFreqDist mfd = incrMFD(mMin, mMax, mNum, aVal, bVal);
		System.out.println(mfd);
		double tmr = 0.0;
		for (int i=0; i<mfd.size(); i++) {
			tmr += MagUtils.magToMoment(mfd.getX(i)) * mfd.getY(i);
		}
		System.out.println("TotMoRate1: " + tmr);
		System.out.println("TotMoRate2: " + mfd.getTotalMomentRate());
		
		
		double totCumRate = getRate(aVal, bVal, 5.25);
		mfd = new GutenbergRichterMagFreqDist(bVal, totCumRate, mMin, mMax, mNum);
		System.out.println(mfd);
		tmr = 0.0;
		for (int i=0; i<mfd.size(); i++) {
			System.out.println("m & rate: " + mfd.getX(i) + " " + mfd.getY(i));
			tmr += MagUtils.magToMoment(mfd.getX(i)) * mfd.getY(i);
		}
		System.out.println("TotMoRate1: " + tmr);
		System.out.println("TotMoRate2: " + mfd.getTotalMomentRate());
		System.out.println("TotMoRate3: " + NSHMP_Utils.totalMoRate(mMin, mNum, 0.1, aVal, bVal));
		
		
//		POW(10,a-b*mag). Do this for each bin until mag bin 6.95.
//		Using my approach and yours should give the same GR distribution, whereas when I use for one of my examples where Mmin = 5.25, Mmax = 7.2, a-val = 0.1556, b-val = 1.23. My approach gives the distribution as-
//
//		           
//		Data[x,y]:
//		5.3      1.2309289E-7
//		5.4      9.273271E-8
//		5.5      6.9860704E-8
//		5.6      5.262995E-8
//		5.7      3.9649066E-8
//		5.8      2.9869845E-8
//		5.9      2.2502613E-8
//		6.0      1.6952468E-8
//		6.1      1.27712365E-8
//		6.2      9.6212815E-9
//		6.3      7.2482464E-9
//		6.4      5.4605067E-9
//		6.5      4.113703E-9
//		6.6      3.099081E-9
//		6.7      2.33471E-9
//		6.8      1.7588667E-9
//		6.9      1.3250521E-9
//		7.0      9.982353E-10
//		7.1      7.520261E-10
//		7.2      3.0329497E-10
//
//		Your approach give the result as -
//		Data[x,y]:
//		5.3      4.3317326E-7
//		5.4      3.2633346E-7
//		5.5      2.4584514E-7
//		5.6      1.852088E-7
//		5.7      1.3952808E-7
//		5.8      1.0511425E-7
//		5.9      7.918841E-8
//		6.0      5.9657026E-8
//		6.1      4.4942954E-8
//		6.2      3.3858026E-8
//		6.3      2.550713E-8
//		6.4      1.9215939E-8
//		6.5      1.4476434E-8
//		6.6      1.0905902E-8
//		6.7      8.2160225E-9
//		6.8      6.189586E-9
//		6.9      4.662959E-9
//		7.0      3.512866E-9
//		7.1      2.6464373E-9
//		7.2      1.993708E-9
//
	}
	
	public static IncrementalMagFreqDist incrMFD(double min, double max, int num, double aVal, double bVal) {
		IncrementalMagFreqDist mfd = new IncrementalMagFreqDist(min, max, num);
		for (int i=0; i<num; i++) {
			double magLo = min + i * 0.1 - 0.05;
			double magHi = min + i * 0.1 + 0.05;
			double rLo = getRate(aVal, bVal, magLo);
			double rHi = getRate(aVal, bVal, magHi);
			mfd.set(i, rLo - rHi);
		}
		return mfd;
	}
	
	public static double getRate(double a, double b, double M) {
		return Math.pow(10,a-b*M);
	}
	
	
}
