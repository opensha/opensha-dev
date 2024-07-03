package scratch.kevin.pointSources;

import org.opensha.commons.calc.GaussianDistCalc;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;

import com.google.common.base.Preconditions;

public class InvCDF_RJBCorrPointSurface extends FiniteApproxPointSurface {
	
	public static class RJBCorrInvPDFs {
		
		public final EvenlyDiscretizedFunc[][] invCmlPDFs;
		public final EvenlyDiscretizedFunc distFunc;
		public final EvenlyDiscretizedFunc magFunc;
		public RJBCorrInvPDFs(CSVFile<String> csv) {
			// first figure out distance and magnitude spacing
			double minDist = csv.getDouble(1, 0);
			double maxDist = csv.getDouble(csv.getNumRows()-1, 0);
			double distDelta = 0d;
			for (int row=2; row<csv.getNumRows(); row++) {
				double dist = csv.getDouble(row, 0);
				if ((float)dist > (float)minDist) {
					distDelta = dist - minDist;
					break;
				}
			}
			int numDist = (int)((maxDist - minDist)/distDelta + 0.5)+1;
			distFunc = new EvenlyDiscretizedFunc(minDist, maxDist, numDist);
			System.out.println("Detected dist range: ["+minDist+", "+maxDist+"]; size="+numDist);
			
			double minMag = csv.getDouble(1, 1);
			double magDelta = csv.getDouble(2, 1) - minMag;
			double maxMag = csv.getDouble(csv.getNumRows()-1, 1);
			int numMag = (int)((maxMag - minMag)/magDelta + 0.5)+1;
			magFunc = new EvenlyDiscretizedFunc(minMag, maxMag, numMag);
			System.out.println("Detected mag range: ["+minMag+", "+maxMag+"]; size="+numMag);
			
			int expectedRows = numMag*numDist + 1;
			Preconditions.checkState(expectedRows == csv.getNumRows(),
					"Bad row count; expected %s, have %s", expectedRows, csv.getNumRows());
			
			int numProbs = csv.getNumCols()-2;
			EvenlyDiscretizedFunc probVals = new EvenlyDiscretizedFunc(0d, 1d, numProbs);
			for (int i=0; i<numProbs; i++)
				Preconditions.checkState((float)probVals.getX(i) == (float)csv.getDouble(0, i+2));
			
			invCmlPDFs = new EvenlyDiscretizedFunc[numDist][numMag];
			
			int row = 1;
			for (int d=0; d<numDist; d++) {
				for (int m=0; m<numMag; m++) {
					double dist = csv.getDouble(row, 0);
					double mag = csv.getDouble(row, 1);
					Preconditions.checkState((float)dist == (float)distFunc.getX(d));
					Preconditions.checkState((float)mag == (float)magFunc.getX(m));
					
					invCmlPDFs[d][m] = probVals.deepClone();
					for (int i=0; i<numProbs; i++)
						invCmlPDFs[d][m].set(i, csv.getDouble(row, i+2));
					
					row++;
				}
			}
		}
		
		public double calcRJB(double dist, double mag, double prob) {
			// snap within bounds
			if (dist < distFunc.getMinX())
				dist = distFunc.getMinX();
			else if (dist > distFunc.getMaxX())
				dist = distFunc.getMaxX();
			if (mag < magFunc.getMinX())
				mag = magFunc.getMinX();
			else if (mag > magFunc.getMaxX())
				mag = magFunc.getMaxX();
			
			int d = distFunc.getClosestXIndex(dist);
			int m = magFunc.getClosestXIndex(mag);
			double ret = invCmlPDFs[d][m].getInterpolatedY(prob);
//			if (d == 10 && (float)mag == 7.15f) {
//				synchronized (InvCDF_RJBCorrPointSurface.class) {
//					System.out.println("Debug for "+(float)dist+" km, M"+(float)mag+", p="+(float)prob);
//					System.out.println("\trJB["+(float)invCmlPDFs[d][m].getMinX()+"]: "+(float)invCmlPDFs[d][m].getY(0));
//					System.out.println("\trJB["+(float)prob+"]: "+(float)ret);
//					System.out.println("\trJB["+(float)invCmlPDFs[d][m].getMaxX()+"]: "+(float)invCmlPDFs[d][m].getY(invCmlPDFs[d][m].size()-1));
//				}
//			}
			return ret;
//			// this seems to be buggy, and would be slow anyway 
//			int distIndex1 = getXIndexBefore(distFunc, dist);
//			Preconditions.checkState(distIndex1 >= 0);
//			int distIndex2;
//			if (distIndex1 == distFunc.size()-1 || (float)dist == (float)distFunc.getX(distIndex1))
//				distIndex2 = distIndex1;
//			else
//				distIndex2 = distIndex1+1;
//			
//			int magIndex1 = getXIndexBefore(magFunc, mag);
//			Preconditions.checkState(magIndex1 >= 0);
//			int magIndex2;
//			if (magIndex1 == magFunc.size()-1 || (float)mag == (float)magFunc.getX(magIndex1))
//				magIndex2 = magIndex1;
//			else
//				magIndex2 = magIndex1+1;
//			
//			double val11 = invCmlPDFs[distIndex1][magIndex1].getInterpolatedY(prob);
//			double val21;
//			if (distIndex2 == distIndex1)
//				val21 = val11;
//			else
//				val21 = invCmlPDFs[distIndex2][magIndex1].getInterpolatedY(prob);
//			double val12, val22;
//			if (magIndex2 == magIndex1) {
//				val12 = val11;
//				val22 = val21;
//			} else {
//				val12 = invCmlPDFs[distIndex1][magIndex2].getInterpolatedY(prob);
//				if (distIndex2 == distIndex1)
//					val22 = val12;
//				else
//					val22 = invCmlPDFs[distIndex2][magIndex2].getInterpolatedY(prob);
//			}
//			
//			double distFrac = (dist - distFunc.getX(distIndex1))/distFunc.getDelta();
//			double magFrac = (mag - magFunc.getX(magIndex1))/magFunc.getDelta();
//			
//			return (1 - magFrac) * ((1 - distFrac)*val11 + distFrac*val12) + 
//				    magFrac * ((1 - distFrac)*val21 + distFrac*val22);
		}
		
		protected int getXIndexBefore(EvenlyDiscretizedFunc func, double x) {
			return (int)Math.floor((x-func.getMinX())/func.getDelta());
		}
	}
	private RJBCorrInvPDFs invCmlPDFs;
	private Location loc;
	private double mag;
	private double rJB_prob;

	public InvCDF_RJBCorrPointSurface(RJBCorrInvPDFs invCmlPDFs, Location loc, double dip, double zTop, double zBot,
			boolean footwall, double length, double mag, double rJB_prob) {
		super(loc, dip, zTop, zBot, footwall, length);
		this.invCmlPDFs = invCmlPDFs;
		this.loc = loc;
		this.mag = mag;
		this.rJB_prob = rJB_prob;
	}
	
	@Override
	public double getDistanceJB(Location siteLoc) {
		double dist = LocationUtils.horzDistanceFast(loc, siteLoc);
		return invCmlPDFs.calcRJB(dist, mag, rJB_prob);
	}
	
	public static DiscretizedFunc getEvenlySpacedProbs(int num) {
		EvenlyDiscretizedFunc edgeFunc = new EvenlyDiscretizedFunc(0d, 1d, num+1);
		EvenlyDiscretizedFunc centerFunc = new EvenlyDiscretizedFunc(0.5*edgeFunc.getDelta(), 1d-0.5*edgeFunc.getDelta(), num);
		for (int i=0; i<num; i++)
			centerFunc.set(i, 1d/(double)num);
		return centerFunc;
	}
	
	public static DiscretizedFunc getSigmaSpacedProbs(int numSigma) {
		// put bin edges at each sigma
		ArbitrarilyDiscretizedFunc edgeFunc = new ArbitrarilyDiscretizedFunc();
		edgeFunc.set(0d, 0d);
		for (int i=0; i<=numSigma; i++) {
			double cdf = GaussianDistCalc.getCDF(i);
			edgeFunc.set(cdf, 0d);
			if (i > 0)
				edgeFunc.set(1d-cdf, 0d);
		}
		edgeFunc.set(1d, 0d);
		ArbitrarilyDiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
		// now calculate weights
		for (int i=1; i<edgeFunc.size(); i++) {
			double lower = edgeFunc.getX(i-1);
			double upper = edgeFunc.getX(i);
			double center = 0.5*(lower+upper);
			func.set(center, upper-lower);
		}
		double sumWeights = func.calcSumOfY_Vals();
		Preconditions.checkState((float)sumWeights == 1f, "SumWeights should be 1, is %s", (float)sumWeights);
		return func;
	}

}
