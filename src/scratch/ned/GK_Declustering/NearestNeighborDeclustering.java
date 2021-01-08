package scratch.ned.GK_Declustering;

import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.magdist.GaussianMagFreqDist;
import org.sparkproject.guava.primitives.Doubles;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.mllib.clustering.GaussianMixture;
import org.apache.spark.mllib.clustering.GaussianMixtureModel;
import org.apache.spark.mllib.clustering.KMeans;
import org.apache.spark.mllib.clustering.KMeansModel;
import org.apache.spark.mllib.linalg.DenseVector;
import org.apache.spark.mllib.linalg.Vector;
import org.apache.spark.mllib.linalg.Vectors;
import org.apache.spark.sql.Dataset;
import org.apache.spark.sql.Encoders;
import org.apache.spark.sql.SparkSession;

import java.util.ArrayList;

import org.apache.commons.math3.analysis.solvers.LaguerreSolver;
import org.apache.commons.math3.complex.Complex;

//import smile.stat.distribution.GaussianMixture;




/**
 * This implements ...
 * 
 * @author field
 *
 */
public class NearestNeighborDeclustering {
	
	static boolean D = true; // debugging flag
	
	double b = 1.0; // GR b-value 
	double d = 1.6; // fractal dimension of catalog
	double q = 0.5; // weighting of temporal (vs spatial) distance in rescaled distances

	double fitWt1;
	double fitWt2;
	double fitMean1;
	double fitMean2;
	double fitSigma1;
	double fitSigma2;
	
	double minNND = -15; //Double.MAX_VALUE;
	double maxNND = 0; //-Double.MAX_VALUE;
	
	ObsEqkRupList fullCatalog, mainshockList, aftershockList, forshockList, bothAftershockAndForeshockList;
	
	int[] indexOfParentArray;
	double[] logNormDistToParentArray, logNormTimeToParentArray, logNNDistanceToParentArray;

	
	int[] numAftershocksForRupArray;
	boolean[] rupIsAfteshockArray, rupIsForeshockArray;

	public NearestNeighborDeclustering(ObsEqkRupList catalog, double minMag) {
		
		// remove M<minMag events
		fullCatalog = new ObsEqkRupList();
		for(ObsEqkRupture rup : catalog)
			if(rup.getMag()>=minMag)
				fullCatalog.add(rup);
		
		if(D) System.out.println("Num Events (M≥"+minMag+"): " + fullCatalog.size());
		
		numAftershocksForRupArray = new int[fullCatalog.size()];
		rupIsAfteshockArray = new boolean[fullCatalog.size()]; // these are initialized to false
		rupIsForeshockArray = new boolean[fullCatalog.size()];

		indexOfParentArray = new int[fullCatalog.size()];
		logNormDistToParentArray = new double[fullCatalog.size()];
		logNormTimeToParentArray = new double[fullCatalog.size()];
		logNNDistanceToParentArray = new double[fullCatalog.size()];
		for(int i=0; i<logNNDistanceToParentArray.length;i++)
			logNNDistanceToParentArray[i] = Double.MAX_VALUE;
		
		
		decluster();
	}
	
	private void decluster() {
		
		// search for aftershocks
		int counter = 0;
		for(int i=fullCatalog.size()-1;i>0;i--) {
			counter += 1;
			if(D && counter == 1000) {
				System.out.println(i+" left to process");
				counter = 0;
			}
			ObsEqkRupture rup = fullCatalog.get(i);
			for(int j=i-1;j>=0;j--) {
				ObsEqkRupture candidateParent = fullCatalog.get(j);
				long timeDiffMills = rup.getOriginTime() - candidateParent.getOriginTime();
				if(timeDiffMills < 0)
					throw new RuntimeException("Error: catalog is not in chronological order");	
				if(timeDiffMills == 0)
					timeDiffMills = 1000;

				double timeDiffYrs = (double)timeDiffMills/(double)(1e3*60*60*24*365.25);
//				double distKm = LocationUtils.linearDistance(rup.getHypocenterLocation(), candidateAftershock.getHypocenterLocation());
				double distKm = LocationUtils.horzDistanceFast(rup.getHypocenterLocation(), candidateParent.getHypocenterLocation());
				if(distKm==0.0) {
					double newLat = candidateParent.getHypocenterLocation().getLatitude() + (-0.05+(0.1*Math.random()));
					double newLon = candidateParent.getHypocenterLocation().getLongitude() + (-0.05+(0.1*Math.random()));
					Location newLoc = new Location(newLat,newLon);
					distKm = LocationUtils.horzDistanceFast(rup.getHypocenterLocation(), newLoc);
//					throw new RuntimeException("Zero Distance");
				}
				double mag = candidateParent.getMag();
				double normDist = Math.pow(distKm, d)*Math.pow(10, -(1-q)*b*mag);
				double normTime = timeDiffYrs*Math.pow(10, -q*b*mag);
				double nnDist = normDist*normTime;
				if(nnDist < logNNDistanceToParentArray[i]) {
					indexOfParentArray[i] = j;
					logNormDistToParentArray[i] = Math.log10(normDist);
					logNormTimeToParentArray[i] = Math.log10(normTime);
					logNNDistanceToParentArray[i] = Math.log10(nnDist);
				}
				if(nnDist==0.0) {
					System.out.println("nnDist="+nnDist);
					System.out.println("normDist="+normDist);
					System.out.println("normTime="+normTime);
					System.out.println("i, j="+i+", "+j);
					System.out.println("timeDiffYrs="+timeDiffYrs);
					System.out.println("distKm="+distKm);
					throw new RuntimeException("Problem related to nnDist being zero");
				}
			}
		}
		
		double[] shorterArray = new double[logNNDistanceToParentArray.length-1];
		for(int i=1; i<logNNDistanceToParentArray.length;i++)
			shorterArray[i-1] = logNNDistanceToParentArray[i];
		
		// Spark gaussian mixture model
		SparkSession spark = SparkSession.builder().master("local").getOrCreate();
		Dataset<Double> dataset = spark.createDataset(Doubles.asList(shorterArray), Encoders.DOUBLE());
		JavaRDD<Vector> parsedData = dataset.toJavaRDD().map(s -> Vectors.dense(s));
		parsedData.cache();
		// default tolerance = 0.01 and maxIterations=100
//		GaussianMixtureModel gmm = new GaussianMixture().setK(2).run(parsedData.rdd());
		
		// this does not do any better
		GaussianMixtureModel gmm = new GaussianMixture().setK(2).setConvergenceTol(0.001).setMaxIterations(1000).run(parsedData.rdd());
		// Output the parameters of the mixture model
		for (int j = 0; j < gmm.k(); j++) {
		  System.out.printf("weight=%f\nmu=%s\nsigma=\n%s\n",
		    gmm.weights()[j], gmm.gaussians()[j].mu(), gmm.gaussians()[j].sigma());
		}
		fitWt1 = gmm.weights()[0];
		fitWt2 = gmm.weights()[1];
		fitMean1 = gmm.gaussians()[0].mu().toArray()[0];
		fitMean2 = gmm.gaussians()[1].mu().toArray()[0];
		fitSigma1 = gmm.gaussians()[0].sigma().toArray()[0];
		fitSigma2 = gmm.gaussians()[1].sigma().toArray()[0];
		
		System.out.println("fitMean1="+fitMean1+"\nfitMean2="+fitMean2+
				"\nfitSigma1="+fitSigma1+"\nfitSigma2="+fitSigma2);
		
		// My attempt to find the roots as in Andrea's matlab code
//		double a = fitSigma1*fitSigma1 - fitSigma2*fitSigma2;
//		double bb = 2*fitMean1*fitSigma2*fitSigma2 - 2*fitMean2*fitSigma1*fitSigma1;
//		double c = fitSigma1*fitSigma1*fitMean2*fitMean2 - fitSigma2*fitSigma2*fitMean1*fitMean1 - 2*fitSigma1*fitSigma1*fitSigma2*fitSigma2*Math.log(fitSigma1/fitSigma2);
//		double[] coefficients = new double[3];
//		coefficients[0] = a;
//		coefficients[1] = bb;
//		coefficients[2] = c;
//		LaguerreSolver solver = new LaguerreSolver();
//		Complex[] complexArray = solver.solveAllComplex(coefficients, -6.0);
//		for(Complex complex : complexArray)
//			System.out.println(complex);
		
		
		// My attemt to use the Smile GaussianMixture model; the means don't look good at all
//		GaussianMixture gm = GaussianMixture.fit(2, shorterArray);
//		System.out.println("mean1: "+gm.components[0].distribution.mean());
//		System.out.println("mean2: "+gm.components[1].distribution.mean());
//		System.out.println("sd1: "+gm.components[0].distribution.sd());
//		System.out.println("sd2: "+gm.components[1].distribution.sd());
//	
//		
//		System.out.println("gm.mean()="+gm.mean());
//		System.out.println("gm.sd()="+gm.sd());
		
		
//		mainshockList = new ObsEqkRupList();
//		aftershockList = new ObsEqkRupList();
//		forshockList = new ObsEqkRupList();
//
//		if(D) {
//			System.out.println("Num Main Shocks: " + mainshockList.size());
//			System.out.println("Num Aftershocks: " + aftershockList.size());
//			System.out.println("Num Foreshocks: " + forshockList.size());
//			
//		}
		
	}
	
	public HistogramFunction getNNDistHistogram() {
		int num = 108;
//		for(int i=1; i<nnDistanceToParentArray.length;i++) {	// first rupture has no parent
//			if(min > Math.log10(nnDistanceToParentArray[i])) min = Math.log10(nnDistanceToParentArray[i]);
//			if(max < Math.log10(nnDistanceToParentArray[i])) max = Math.log10(nnDistanceToParentArray[i]);
//		}	
//		System.out.println("min="+min+"\nmax="+max);
		HistogramFunction hist = new HistogramFunction(minNND, maxNND, num);
		for(int i=1; i<logNNDistanceToParentArray.length;i++) {	// first rupture has no parent
			hist.add(logNNDistanceToParentArray[i], 1.0);
		}
		hist.normalizeBySumOfY_Vals();
		hist.scale(1.0/hist.getDelta());
		System.out.println("Test Hist: "+hist.calcSumOfY_Vals()*hist.getDelta());

		return hist;
	}
	
	
	public ArrayList<HistogramFunction> getFitGaussianFunctions() {
		int num = 1000;
		GaussianMagFreqDist g1 = new GaussianMagFreqDist(minNND,maxNND,num,fitMean1,fitSigma1,1.0);
		g1.normalizeToPDF(); // make a pdf
		HistogramFunction hist1 = new HistogramFunction(minNND,maxNND,num);
		for(int i=0;i<hist1.size();i++)
			hist1.set(i,g1.getY(i)*fitWt1);
		hist1.setName(g1.getName());
		hist1.setInfo(g1.getInfo());
		
		GaussianMagFreqDist g2 = new GaussianMagFreqDist(minNND,maxNND,num,fitMean2,fitSigma2,1.0);
		g2.normalizeToPDF(); // make a pdf
		HistogramFunction hist2 = new HistogramFunction(minNND,maxNND,num);
		for(int i=0;i<hist2.size();i++)
			hist2.set(i,g2.getY(i)*fitWt2);
		hist2.setName(g2.getName());
		hist2.setInfo(g2.getInfo());
		
		// make sum of the two
		HistogramFunction histSum = new HistogramFunction(minNND,maxNND,num);
		for(int i=0;i<histSum.size();i++)
			histSum.set(i,hist1.getY(i)+hist2.getY(i));
		histSum.setName("Sum of two Gaussian fits");

		
		ArrayList<HistogramFunction> list = new ArrayList<HistogramFunction>();
		list.add(hist1);
		list.add(hist2);
		list.add(histSum);
		return list;
	}
	
	
	
	/**
	 * This returns the declustered catalog
	 * @return
	 */
	public ObsEqkRupList getDeclusteredCatalog() {
		return mainshockList;
	}
	
	/**
	 * This is a convenience method to avoid having to instantiate this class if 
	 * only the declusted catalog is desired.
	 * @param rupList - must be in chronological order
	 * @return
	 */
	public static ObsEqkRupList getDeclusteredCatalog(ObsEqkRupList rupList, double minMag) {
		NearestNeighborDeclustering gk_decluster = new NearestNeighborDeclustering(rupList, minMag);
		return gk_decluster.getDeclusteredCatalog();
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}
	
	
}
