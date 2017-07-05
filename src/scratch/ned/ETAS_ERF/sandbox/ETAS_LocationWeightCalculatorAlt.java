package scratch.ned.ETAS_ERF.sandbox;

import java.util.ArrayList;

import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.gui.plot.GraphWindow;

import scratch.UCERF3.erf.ETAS.ETAS_Utils;
import scratch.UCERF3.erf.ETAS.IntegerPDF_FunctionSampler;

public class ETAS_LocationWeightCalculatorAlt {
	
	final static boolean D = true;
	
	int numLatLon, numDepth;
	double maxLatLonDeg, maxDepthKm, latLonDiscrDeg, depthDiscr, midLat, maxDistKm;
	
	double distDecay, minDist;
	
	double cosMidLat;
	
	double[] totWtAtDepth;	
	double[][][] pointWt;
	
	double histLogMin=-2.0;	// log10 distance
	double histLogMax = 4.0;	// log10 distance
	int histNum = 31;
	
	LocationList[][][] subLocsArray;
	IntegerPDF_FunctionSampler[][][] subLocSamplerArray;

	int[] numSubDistances = {100,20,10,5,2,2};
	
	EvenlyDiscretizedFunc logTargetDecay;
	EvenlyDiscretizedFunc logDistWeightHist;
	
	/**
	 * 
	 * @param maxDistKm - the maximum distance for sampling in km
	 * @param maxDepthKm - the max seismogenic thickness
	 * @param latLonDiscrDeg - the lat and lon discretization in degrees (0.2 is recommented)
	 * @param depthDiscr - the depth discretization in km (2.0 is recommended)
	 * @param midLat - the mid latitude used to compute bin widths (since widths decrease with latitude)
	 * @param distDecay - the ETAS distance decay parameter
	 * @param minDist - the ETAS min distance
	 */
	public ETAS_LocationWeightCalculatorAlt(double maxDistKm, double maxDepthKm, double latLonDiscrDeg, double depthDiscr, 
			double midLat, double distDecay, double minDist) {
		
		cosMidLat = Math.cos(midLat*Math.PI/180);
		double aveLatLonDiscrKm = (latLonDiscrDeg+cosMidLat*latLonDiscrDeg)*111/2.0;
		this.maxDistKm = maxDistKm;
		this.maxLatLonDeg = maxDistKm/(111*cosMidLat);	// degrees
		
		this.maxDepthKm = maxDepthKm;
		this.latLonDiscrDeg = latLonDiscrDeg;
		this.depthDiscr = depthDiscr;
		this.midLat = midLat;
		this.distDecay=distDecay;
		this.minDist=minDist;
				
		// the number of points in each direction
		numLatLon = (int)Math.round(maxLatLonDeg/latLonDiscrDeg);
		numDepth = (int)Math.round(maxDepthKm/depthDiscr);
	
		if (D) System.out.println("aveLatLonDiscrKm="+aveLatLonDiscrKm+
				"\nmaxLatLonDeg="+maxLatLonDeg+
				"\ncosMidLat="+cosMidLat+
				"\nnumLatLon="+numLatLon+
				"\nnumDepth="+numDepth);
		
		// the following is info for the close points that are subdivided
		int maxNumPtsWithSubLocs = numSubDistances.length;
		subLocsArray = new LocationList[maxNumPtsWithSubLocs][maxNumPtsWithSubLocs][maxNumPtsWithSubLocs];
		subLocSamplerArray = new IntegerPDF_FunctionSampler[maxNumPtsWithSubLocs][maxNumPtsWithSubLocs][maxNumPtsWithSubLocs];

		pointWt = new double[numLatLon][numLatLon][numDepth];
		
//getSubDistances(0, 1, 2, 4);
//System.exit(0);

		// make distances weight histogram for the various log distances
		logDistWeightHist = new EvenlyDiscretizedFunc(histLogMin,histLogMax,histNum);
		logDistWeightHist.setTolerance(logDistWeightHist.getDelta());
		double[] distances=null;
		for(int iDep=0;iDep<numDepth; iDep++) {
			if (D) System.out.println("Working on depth "+iDep);
			for(int iLat=0;iLat<numLatLon; iLat++) {
				for(int iLon=0;iLon<numLatLon; iLon++) {
					// find the largest index) proxy for farthest distance
					int maxIndex = Math.max(iDep, Math.max(iLat, iLon));
					if(maxIndex<numSubDistances.length) {
						distances = getSubDistances(iLat, iLon, iDep, numSubDistances[maxIndex]);
						for(int i=0;i<distances.length;i++) {
							double dist = distances[i];
							double wt = ETAS_Utils.getDistDecayValue(dist, minDist, -distDecay);
							double logDist = Math.log10(dist);
							if(logDist<logDistWeightHist.getX(0))	// in case it's below the first bin
								logDistWeightHist.add(0, wt);
							else if (dist<maxDistKm)
								logDistWeightHist.add(logDist,wt);
//							if(maxIndex<7)
//							System.out.println(maxIndex+"\t"+distances.length);
						}
					}
					else {
						double dist = getDistance(iLat, iLon, iDep);
						if(dist<maxDistKm) {
							double wt = ETAS_Utils.getDistDecayValue(dist, minDist, -distDecay);;
							logDistWeightHist.add(Math.log10(dist),wt);							
						}
					}
				}
			}
		}
		
		// plot to check for any zero bins
		if (D) {
			GraphWindow graph = new GraphWindow(logDistWeightHist, "test hist"); 
		}

		// make target distances decay histogram (this is what we will match_
		logTargetDecay = new EvenlyDiscretizedFunc(histLogMin,histLogMax,histNum);
		logTargetDecay.setTolerance(logDistWeightHist.getDelta());
		double logBinHalfWidth = logTargetDecay.getDelta()/2;
		double upperBinEdge = Math.pow(10,logTargetDecay.getX(0)+logBinHalfWidth);
		double lowerBinEdge;
		double binWt = ETAS_Utils.getDecayFractionInsideDistance(distDecay, minDist, upperBinEdge);	// everything within the upper edge of first bin
		logTargetDecay.set(0,binWt);
		for(int i=1;i<logTargetDecay.size();i++) {
			double logLowerEdge = logTargetDecay.getX(i)-logBinHalfWidth;
			lowerBinEdge = Math.pow(10,logLowerEdge);
			double logUpperEdge = logTargetDecay.getX(i)+logBinHalfWidth;
			upperBinEdge = Math.pow(10,logUpperEdge);
			double wtLower = ETAS_Utils.getDecayFractionInsideDistance(distDecay, minDist, lowerBinEdge);
			double wtUpper = ETAS_Utils.getDecayFractionInsideDistance(distDecay, minDist, upperBinEdge);
			binWt = wtUpper-wtLower;
			logTargetDecay.set(i,binWt);
		}
		
		// normalize
		double tot = logTargetDecay.calcSumOfY_Vals();
		System.out.println("logWtHistogram.calcSumOfY_Vals()="+tot);
//		logTargetDecay.scale(1.0/tot);
//		System.out.println("logWtHistogram.calcSumOfY_Vals()="+logTargetDecay.calcSumOfY_Vals());
		
//		GraphWindow graph2 = new GraphWindow(logTargetDecay, "logTargetDecay"); 
//		GraphWindow graph3 = new GraphWindow(logWtHistogram2, "logWtHistogram2"); 

		// now fill in weights for each point
		for(int iDep=0;iDep<numDepth; iDep++) {
			System.out.println("Working on depth "+iDep);
			for(int iLat=0;iLat<numLatLon; iLat++) {
				for(int iLon=0;iLon<numLatLon; iLon++) {
					// find the largest index) proxy for farthest distance
					int maxIndex = Math.max(iDep, Math.max(iLat, iLon));
					if(maxIndex<numSubDistances.length) {
						distances = getSubDistances(iLat, iLon, iDep, numSubDistances[maxIndex]);
						for(int i=0;i<distances.length;i++) {
							double dist = distances[i];
							double wt = ETAS_Utils.getDistDecayValue(dist, minDist, -distDecay);
							double logDist = Math.log10(dist);
							if(logDist<logDistWeightHist.getX(0))	{// in case it's below the first bin
								pointWt[iLat][iLon][iDep] += wt*logTargetDecay.getY(0)/logDistWeightHist.getY(0);
							}
							else if (dist<maxDistKm) {
								pointWt[iLat][iLon][iDep] += wt*logTargetDecay.getY(logDist)/logDistWeightHist.getY(logDist);
							}
						}
					}
					else {
						double dist = getDistance(iLat, iLon, iDep);
						if(dist<maxDistKm) {
							double wt = ETAS_Utils.getDistDecayValue(dist, minDist, -distDecay);
							double logDist = Math.log10(dist);
							pointWt[iLat][iLon][iDep] += wt*logTargetDecay.getY(logDist)/logDistWeightHist.getY(logDist);							
						}
					}
				}
			}
		}
		
		// test total weight
		double totWtTest=0;
		for(int iDep=0;iDep<numDepth; iDep++) {
			for(int iLat=0;iLat<numLatLon; iLat++) {
				for(int iLon=0;iLon<numLatLon; iLon++) {
					totWtTest += pointWt[iLat][iLon][iDep];
				}
			}
		}
		System.out.println("totWtTest = "+ totWtTest);
		
		// APPLY THE ABOVE NORMALIZATION?
		
		
		// THE FOLLOWING NEEDS WORK?
		totWtAtDepth = new double[numDepth];
		double testTot=0;
		for(int iDep=0;iDep<numDepth; iDep++) {
			double wtAtDep=0;
			for(int iLat=0;iLat<numLatLon; iLat++) {
				for(int iLon=0;iLon<numLatLon; iLon++) {
					wtAtDep += pointWt[iLat][iLon][iDep];
				}
			}
			totWtAtDepth[iDep]=wtAtDep;
//			System.out.println("totWtAtDepth\t"+iDep+"\t"+wtAtDep);
			testTot += wtAtDep;
		}

		
		
	}
	
	
	/**
	 * This returns a location containing delta lat, lon, and depth based on distance decay
	 * @param relLat
	 * @param relLon
	 * @param relDep
	 * @return
	 */
	public Location getRandomDeltaLoc(double relLat, double relLon, double relDep) {
		int iLat = getLatIndex(relLat);
		int iLon = getLonIndex(relLon);
		int iDep = getDepthIndex(relDep);
		Location loc;	// the location before some added randomness
		double deltaSubLatLon;
		double deltaDepth;
		
		int maxIndex = Math.max(iDep, Math.max(iLat, iLon));
		if(maxIndex<numSubDistances.length) {
			int numSubLoc = numSubDistances[maxIndex];
			deltaSubLatLon = latLonDiscrDeg/numSubLoc;
			deltaDepth = depthDiscr/numSubLoc;
	

			if(subLocsArray[iLat][iLon][iDep] == null) {
				double midLat = getLat(iLat);
				double midLon = getLon(iLon);
				double midDepth = getDepth(iDep);
				
//if(iLat==0 && iLon==1 && iDep==2) {
//					System.out.println("midLat="+midLat+"\tmidLon="+midLon+"\tmidDepth="+midDepth);
//					System.out.println("relLat\trelLon\trelDepth\tdist");
//				}
				LocationList locList = new LocationList();
				IntegerPDF_FunctionSampler newSampler = new IntegerPDF_FunctionSampler(numSubLoc*numSubLoc*numSubLoc);
				int index = 0;
				for(int iSubLat = 0; iSubLat < numSubLoc; iSubLat++) {
					double lat = (midLat-latLonDiscrDeg/2) + iSubLat*deltaSubLatLon + deltaSubLatLon/2;
					for(int iSubLon = 0; iSubLon < numSubLoc; iSubLon++) {
						double lon = (midLon-latLonDiscrDeg/2) + iSubLon*deltaSubLatLon + deltaSubLatLon/2;
						for(int iSubDep = 0; iSubDep < numSubLoc; iSubDep++) {
							double dep = (midDepth-depthDiscr/2) + iSubDep*deltaDepth + deltaDepth/2;
							locList.add(new Location(lat-midLat,lon-midLon,dep-midDepth));	// add the deltaLoc to list
							double dist = getDistance(lat, lon, dep);
							double logDist = Math.log10(dist);
							double wt = ETAS_Utils.getDistDecayValue(dist,minDist, -distDecay);
							double normWt;
							if(logDist<logDistWeightHist.getX(0))
								normWt = logTargetDecay.getY(0)/logDistWeightHist.getY(0);
							else
								normWt = logTargetDecay.getY(logDist)/logDistWeightHist.getY(logDist);
							newSampler.add(index, wt*normWt);		// add the sampler
//if(iLat==0 && iLon==1 && iDep==2) {
//								System.out.println((float)lat+"\t"+(float)lon+"\t"+(float)dep+"\t"+(float)dist+"\t"+(float)wt+"\t"+(float)normWt);
//							}
							index ++;
						}
					}
				}
				subLocsArray[iLat][iLon][iDep] = locList;
				subLocSamplerArray[iLat][iLon][iDep] = newSampler;
			}
			
			int randLocIndex = subLocSamplerArray[iLat][iLon][iDep].getRandomInt();
			loc = subLocsArray[iLat][iLon][iDep].get(randLocIndex);			
		}
		else {
			deltaSubLatLon = latLonDiscrDeg;
			deltaDepth = depthDiscr;
			loc = new Location(0, 0, 0);	// no delta
		}
		// ADD A RANDOM ELEMENT
		return new Location(loc.getLatitude()+deltaSubLatLon*(Math.random()-0.5)*0.999,
				loc.getLongitude()+deltaSubLatLon*(Math.random()-0.5)*0.999,
				loc.getDepth()+deltaDepth*(Math.random()-0.5)*0.999);
//		return loc;
		
	}
	
	/**
	 * Get the distance (km) to the given point
	 * @param iLat
	 * @param iLon
	 * @param iDep
	 * @return
	 */
	private double getDistance(int iLat, int iLon, int iDep) {
		return getDistance(getLat(iLat), getLon(iLon), getDepth(iDep));
	}
	
	/**
	 * Get the distance (km) to the given location (approx distance calculation is applied)
	 * @param relLat
	 * @param relLon
	 * @param relDep
	 * @return
	 */
	private double getDistance(double relLat, double relLon, double relDep) {
		double latDistKm = relLat*111;
		double lonDistKm = relLon*111*cosMidLat;
		return Math.sqrt(latDistKm*latDistKm+lonDistKm*lonDistKm+relDep*relDep);
	}

	
	/**
	 * This give the probability of an event at the given point, and for the 
	 * given main shock hypocenter depth
	 * @param relLat
	 * @param relLon
	 * @param relDep
	 * @param hypoDep
	 * @return
	 */
	public double getProbAtPoint(double relLat, double relLon, double relDep, double hypoDep) {
		int iLat = getLatIndex(relLat);
		int iLon = getLonIndex(relLon);
		int iDep = getDepthIndex(relDep);
		int iHypoDep = getDepthIndex(hypoDep);
		
		// solve for the total weight for the associated layers
		double normWt=0;
		// sum those at same depth and below	// if at surface (iHypoDepth=0), should include all; if at bottom (iHypoDepth=numDepth-1), should include just 0th
		for(int d=0; d<numDepth-iHypoDep;d++)
			normWt += totWtAtDepth[d];
		// sum those above; none if iHypoDepth=0; those above if at bottom (iHypoDepth=numDepth-1)
		if(iHypoDep > 0)
			for(int d=1; d<=iHypoDep;d++)
				normWt += totWtAtDepth[d];
		
		if(iLat >= numLatLon || iLon >= numLatLon || iDep >= numDepth) {
//			System.out.println("relLat="+relLat+"\tiLat="+iLat);
//			System.out.println("relLon="+relLon+"\tiLon="+iLon);
//			System.out.println("relDep="+relDep+"\tiDep="+iDep);
			return 0;
		}

		// factor of four below is to account for the other 3 quadrants
		return pointWt[iLat][iLon][iDep]/(normWt*4);
	}
	
	private double getLat(int iLat) {
		return iLat*latLonDiscrDeg+latLonDiscrDeg/2.0;
	}
	
	private int getLatIndex(double  relLat) {
		return (int) Math.round((relLat-latLonDiscrDeg/2.0)/latLonDiscrDeg);
	}

	
	private double getLon(int iLon) {
		return iLon*latLonDiscrDeg+latLonDiscrDeg/2.0;
	}
	
	private int getLonIndex(double  relLon) {
		return (int) Math.round((relLon-latLonDiscrDeg/2.0)/latLonDiscrDeg);
	}

	private double getDepth(int iDep) {
		return iDep*depthDiscr+depthDiscr/2.0;
	}
	
	private int getDepthIndex(double relDepth) {
		return (int)Math.round((relDepth-depthDiscr/2.0)/depthDiscr);
	}



	public void testRandomSamples(int numSamples) {
		
		//test
//		getRandomDeltaLoc(this.getLat(0), this.getLon(1), this.getDepth(2));
//		System.exit(0);
		
		
		IntegerPDF_FunctionSampler sampler;
		int totNumPts = numLatLon*numLatLon*numDepth;
		sampler = new IntegerPDF_FunctionSampler(totNumPts);
		int[] iLatArray = new int[totNumPts];
		int[] iLonArray = new int[totNumPts];
		int[] iDepArray = new int[totNumPts];
		int index=0;
		for(int iDep=0;iDep<numDepth; iDep++) {
			for(int iLat=0;iLat<numLatLon; iLat++) {
				for(int iLon=0;iLon<numLatLon; iLon++) {
					sampler.set(index,pointWt[iLat][iLon][iDep]);
					iLatArray[index]=iLat;
					iLonArray[index]=iLon;
					iDepArray[index]=iDep;
					index +=1;
				}
			}
		}
		
		// create histogram
		EvenlyDiscretizedFunc testLogHistogram = new EvenlyDiscretizedFunc(histLogMin,histLogMax,histNum);
		testLogHistogram.setTolerance(testLogHistogram.getDelta());
		
		EvenlyDiscretizedFunc testHistogram = new EvenlyDiscretizedFunc(0.5 , 1009.5, 1010);
		testHistogram.setTolerance(testHistogram.getDelta());
		
		// make target histogram
		EvenlyDiscretizedFunc targetHist = new EvenlyDiscretizedFunc(0.5 , 999.5, 1000);
		double halfDelta=targetHist.getDelta()/2;
		for(int i=0;i<targetHist.size();i++) {
			double upper = ETAS_Utils.getDecayFractionInsideDistance(distDecay, minDist, targetHist.getX(i)+halfDelta);
			double lower = ETAS_Utils.getDecayFractionInsideDistance(distDecay, minDist, targetHist.getX(i)-halfDelta);
			targetHist.set(i,upper-lower);
		}

		for(int i=0;i<numSamples;i++) {
			int sampIndex = sampler.getRandomInt();
			double relLat = getLat(iLatArray[sampIndex]);
			double relLon = getLon(iLonArray[sampIndex]);
			double relDep = getDepth(iDepArray[sampIndex]);
			Location deltaLoc=getRandomDeltaLoc(relLat, relLon, relDep);
			double dist = getDistance(relLat+deltaLoc.getLatitude(), relLon+deltaLoc.getLongitude(), relDep+deltaLoc.getDepth());
			if(dist<this.maxDistKm) {
				testHistogram.add(dist, 1.0/numSamples);
				double logDist = Math.log10(dist);
				if(logDist<testLogHistogram.getX(0))
					testLogHistogram.add(0, 1.0/numSamples);
				else if (logDist<histLogMax)
					testLogHistogram.add(logDist,1.0/numSamples);
			}
		}
		
		ArrayList funcs1 = new ArrayList();
		funcs1.add(testLogHistogram);
		funcs1.add(logTargetDecay);

		GraphWindow graph = new GraphWindow(funcs1, "testLogHistogram"); 
		graph.setAxisRange(-2, 3, 1e-6, 1);
		graph.setYLog(true);

		
		ArrayList funcs2 = new ArrayList();
		funcs2.add(testHistogram);
		funcs2.add(targetHist);
		GraphWindow graph2 = new GraphWindow(funcs2, "testHistogram"); 
		graph2.setAxisRange(0.1, 1000, 1e-6, 1);
		graph2.setYLog(true);
		graph2.setXLog(true);
	}
	
	
	
	private double[] getSubDistances(int iLat, int iLon, int iDep, int numDiscr) {
		double[] distances = new double[numDiscr*numDiscr*numDiscr];
		double midLat = getLat(iLat);
		double midLon = getLon(iLon);
		double midDepth = getDepth(iDep);
		double deltaSubLatLon = latLonDiscrDeg/numDiscr;
		double deltaDepth = depthDiscr/numDiscr;
		int index=0;
//System.out.println("midLat="+midLat+"\tmidLon="+midLon+"\tmidDepth="+midDepth);
//System.out.println("relLat\trelLon\trelDepth\tdist");

		for(int latIndex = 0; latIndex < numDiscr; latIndex++) {
			double relLat = (midLat-latLonDiscrDeg/2) + latIndex*deltaSubLatLon + deltaSubLatLon/2;
			for(int lonIndex = 0; lonIndex < numDiscr; lonIndex++) {
				double relLon = (midLon-latLonDiscrDeg/2) + lonIndex*deltaSubLatLon + deltaSubLatLon/2;
				for(int depIndex = 0; depIndex < numDiscr; depIndex++) {
					double relDep = (midDepth-depthDiscr/2) + depIndex*deltaDepth + deltaDepth/2;
					distances[index] = getDistance(relLat, relLon, relDep);
//System.out.println((float)relLat+"\t"+(float)relLon+"\t"+(float)relDep+"\t"+(float)distances[index]);

					index+=1;
				}
			}
		}
		return distances;
	}




	/**
	 * @param args
	 */
	public static void main(String[] args) {
//		ETAS_LocationWeightCalculator calc = new ETAS_LocationWeightCalculator(1000.0, 24.0, 0.01, 1.0, 38.0, 2, 0.3);
		
		double maxDistKm=1000.0;
		double maxDepthKm=24;
//		double latLonDiscrDeg=0.005;
//		double depthDiscr=0.5;
		double latLonDiscrDeg=0.02;
		double depthDiscr=2.0;
		double midLat=38;
		double distDecay=2;
		double minDist=0.3;
		ETAS_LocationWeightCalculatorAlt calc = new ETAS_LocationWeightCalculatorAlt(maxDistKm, maxDepthKm, latLonDiscrDeg, 
				depthDiscr, midLat, distDecay, minDist);
		
//		calc.getRandomDeltaLoc(.01, .01, 0.1);
		
//		for(int i=0;i<100;i++) {
//			Location loc = calc.getRandomDeltaLoc(.01, .01, 0.1);
//			System.out.println(loc.getLatitude()+"\t"+loc.getLongitude()+"\t"+loc.getDepth());
//		}
		
		System.out.println("Testing random samples...");
//		calc.getEquivDistanceFast(0, 0, 0, 100);
		calc.testRandomSamples(1000000);

	}

}
