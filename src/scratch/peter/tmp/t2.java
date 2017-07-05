package scratch.peter.tmp;

import java.util.concurrent.TimeUnit;

import org.opensha.commons.data.function.EvenlyDiscretizedFunc;

/**
 * Add comments here
 *
 * 
 * @author Peter Powers
 * @version $Id:$
 * 
 */
public class t2 {

	
	public static void main(String[] args) {
		
		
		EvenlyDiscretizedFunc ff = new EvenlyDiscretizedFunc(5.0, 1, 0.0);
		System.out.println(ff);
		
		double ppp = 1e-18;
		double scale = 1 + 1e-16;
		double yy = 0.99999999999999999;
//		System.out.println(ppp * scale);
//		System.out.println(yy);
//		for(double x:testVals) {
//			int index = (int)Math.round((x-minX)/delta);
////			int index = (delta == 0) ? 0 : (int) Math.round((x-minX)/delta);
//			double binCenter = minX+delta*index;
//			System.out.println(x+"\tindex = "+index+"\tbinCenter="+binCenter);
//		}

//		for(double x:testVals) {
//			
//			int index = getXIndex(x);
//			double binCenter = minX+delta*index;
//			System.out.println(x+"\tindex = "+index+"\tbinCenter="+binCenter);
//		}
//		double yy = 0.99999999999999999;
//		System.out.println((int) yy);
		EvenlyDiscretizedFunc f = new EvenlyDiscretizedFunc(minX, num, delta);
		f.setTolerance(0.15);
//		System.out.println(f);
		for(double x:testVals) {
			int index = f.getXIndex(x);
			int cindex = f.getClosestXIndex(x);
//			System.out.println(x + " " + getClosestXIndex(x));
			double binCenter = minX+delta*index;
			
//			int index = (int)Math.round((x-minX)/delta);
//			double binCenter = minX+delta*index;
//			int alt = (int) Math.rint(((x+delta/2)-minX)/delta);
//			System.out.println(alt);
			System.out.println(x+"\tci = "+cindex+"\tindex = "+index+"\tbinCenter="+binCenter);
		}
		
		
		

	}
	
	static double minX=5.23;
	static double delta = 0.42;
	static int num = 6;

//	static double minX=6.05;
//	static double delta = 0.1;
//	static int num = 10;

	static double[] testVals = {5.02,6,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,7, 7.54,6.499999999999999, 6.49999999999999};
	static double tolerance = 0.05;
	
//	public static int getXIndex( double x) throws Point2DException{
//		if (x < (minX-tolerance) || x > (maxX+tolerance))
//			return -1;
//		// single value functions with delta=0 causes problems
//		System.out.println(((x-minX)/delta));
//		int i = (delta == 0) ? 0 : (int) Math.round((x-minX)/delta);
//		//int i = (int)Math.round((double)((x-minX)/delta));
//		double closeX = funcVals[i]; //getX(i);
//		if( withinTolerance(x, closeX))
//			return i;
//		else
//			return -1;
//	}
//
//	protected static boolean withinTolerance(double x, double xx){
//		return ( Math.abs( x - xx)  <= tolerance);
//	}
	
	public static int getClosestXIndex( double x) {
		double tmp = x / delta;
		double rtmp = Math.rint(tmp);
		int itmp = (int) rtmp;
		double minTmp = minX / delta;
		double rMinTmp = Math.rint(minTmp);
		int iMinTmp = (int) rMinTmp;
		int iDelta = (int) (delta * Math.rint(1 / delta));
		int index = (itmp - iMinTmp) / iDelta;
		
		System.out.println(x+delta/2);
		System.out.println((x+delta/2)-minX);
		System.out.println(((x+delta/2)-minX)/delta);
		int alt = (int) Math.rint(((x+delta/2)-minX)/delta); //(int) Math.rint((x-minX)/delta);
		
		
		
//		System.out.println(tmp + " " + rtmp + " "+itmp+ " "+minTmp+ " "+rMinTmp+ " "+iMinTmp+ " " + pp);
		System.out.println(itmp+ " "+iMinTmp + " " + iDelta + " " + index + " " + alt);
		
		// single value functions with delta=0 causes problems
		int i = (delta == 0) ? 0 : (int) Math.round((x-minX)/delta);
		//int i = (int)Math.round((double)((x-minX)/delta));
		return (i<0) ? 0 : (i>=num) ? num-1 : i;
//		if(i<0) i=0;
//		else if(i>=num) i = num-1;
//		return i;
		
	}



}
