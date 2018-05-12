package scratch.aftershockStatisticsETAS;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;

import com.google.common.base.Stopwatch;

/**
 * This class is a stochastic ETAS catalog that represents an extension of the observed catalog,
 * supplied as mainshock and aftershocks objects. Time is relative to the mainshock. If nSims is supplied,
 * a suite of catalogs will be built. Only the last catalog to be built is stored for retrieval.
 * 
 * @param genericETAS_parameters or a,a_sigma,b,p,c,alpha,refMag
 * @param mainshock	
 * @param aftershocks
 * @param startTime
 * @param endTime
 * @param magMag
 * @param maxGeneration
 * @param [nSims]
 * 
 * @author Nicholas van der Elst
 *
 */

public class ETAScatalog {

	private final static boolean D = false; //debug
	
	double[] ams_vec, a_vec, p_vec, c_vec;
	double[][][][] likelihood;
	double alpha;
	double b;
	double refMag;
	double forecastStart;
	double forecastEnd;
	double Mc;
	double maxMagLimit;
	int nSims;
	int maxGenerations; //simulation depth
	private double[] maxMags;	
	private int[] numEventsFinal;
	private int[] numGenerations;
	
	private List<List<float[]>> catalogList;	//list of catalogs
	
	
	public ETAScatalog(double[] ams_vec, double[] a_vec, double[] p_vec, double[] c_vec, double[][][][] likelihood, double alpha, double b, double refMag,
			ObsEqkRupture mainshock, ObsEqkRupList aftershocks,
			double dataStart, double dataEnd, double forecastStart, double forecastEnd, double Mc, double maxMag, int maxGenerations, int nSims){
	
		this.ams_vec = ams_vec;
		this.a_vec = a_vec;
		this.p_vec = p_vec;
		this.c_vec = c_vec;
		this.likelihood = likelihood;
		this.alpha = alpha;
		this.b = b;
		this.refMag = refMag;
		this.forecastStart = forecastStart;
		this.forecastEnd = forecastEnd;
		this.Mc = Mc;
		this.maxMagLimit = maxMag;
		this.maxGenerations = maxGenerations;
		this.nSims = nSims;
		
		if(D) System.out.println("ETAS simulation params: alpha=" + alpha + " b=" + b + " Mref=" + refMag + " Mc=" + Mc + " Mmax=" + maxMag + " nSims=" + nSims); 
				
		
		List<float[]> newEqList = new ArrayList<float[]>();	//catalog containing list of {time, mag, gen}
		List<List<float[]>> catalogList = new ArrayList<List<float[]>>(); //list of catalogs
		
		double[] maxMags = new double[nSims];
		int[] nEvents = new int[nSims];
		int[] nGens = new int[nSims];
		
		if(D) System.out.println("Calculating " + nSims + " " + (int)(forecastEnd - forecastStart) + "-day ETAS catalogs...");
		
		double[][] paramList;
		if(nSims>0){
			//get the list of parameters to supply to each simulation
			paramList = sampleParams(nSims, maxMag);
			
			Stopwatch watch = Stopwatch.createStarted();
			int warnTime = 3;
			long toc;
			double timeEstimate;
			String initialMessageString = "Calculating " + nSims + " " + (int)(forecastEnd - forecastStart) + "-day ETAS catalogs. ";
			
			for(int i = 0; i < nSims ; i++){
				toc = watch.elapsed(TimeUnit.SECONDS);
				if (toc > warnTime){
					warnTime += 10;
					timeEstimate = (double)toc * (double)(nSims)/(double)i;
					System.out.format(initialMessageString + "Approximately %d seconds remaining...\n", (int) ((timeEstimate - toc)));
					initialMessageString = "...";
				}

				double[] params = paramList[i];
				double ams_sample, a_sample, p_sample, c_sample;
				ams_sample = params[0];
				a_sample = params[1];
				p_sample = params[2];
				c_sample = params[3];
				
				if (D && Math.floorMod(i, nSims/10) == 0) System.out.println("Parameter set " + i + ": " + ams_sample + " " + a_sample + " " + p_sample + " " + c_sample);

				// Currently sets the first event as mainshock and adjusts magnitude
				// todo step1: change magnitude of LARGEST earthquake
				// todo step2: depending on the total number of vents, adjust N-largest magnitudes
				ObsEqkRupture simulationMainshock = (ObsEqkRupture) mainshock.clone();
				simulationMainshock.setMag(mainshock.getMag() + (ams_sample - a_sample));

				newEqList = getNewETAScatalog(simulationMainshock, aftershocks, a_sample, p_sample, c_sample, i);
				maxMags[i] = get_maxMag(newEqList);
				nEvents[i] = get_nEvents(newEqList);
				nGens[i] = get_nGenerations(newEqList);
				catalogList.add(i, newEqList);
			}
			toc = watch.elapsed(TimeUnit.SECONDS);
			if(D) System.out.println("It took " + toc + " seconds to generate stochastic catalogs.");
			watch.stop();
		}
		
		//this.eqList = getLastETAScatalog();
		this.catalogList = catalogList;
		this.maxMags = maxMags;
		this.numEventsFinal = nEvents;
		this.numGenerations = nGens;
	}
		
	public List<float[]> getNewETAScatalog(ObsEqkRupture mainshock, ObsEqkRupList aftershocks, double a_sample, double p_sample, double c_sample, int simNumber){
		
		//extract magnitudes and times from supplied Eqk rupture objects to make catalog (combine MS and AS's)
		List<float[]> newEqList = new ArrayList<float[]>();
		List<float[]> finalEqList = new ArrayList<float[]>();
		
		//double[] event = {0, mainshock.getMag(), 0};	//utility variable : {relative time, magnitude, generation number}

		double t0 = mainshock.getOriginTime(); //in milliseconds
		
		//combine lists
		ObsEqkRupList seedQuakes = new ObsEqkRupList();
		seedQuakes.add(mainshock);
		Collections.reverse(aftershocks);
		seedQuakes.addAll(aftershocks);

		//int counter = 0;
		//go through seed (observed) earthquake list, add each event to a pared-down eventList and add simulated children
		for(ObsEqkRupture rup : seedQuakes){
			float[] event = new float[3];
			event[0] = (float) ((rup.getOriginTime() - t0)/ETAS_StatsCalc.MILLISEC_PER_DAY);	//elapsed time in days
			event[1] = (float) rup.getMag();	
			event[2] = 0;	//generation number
			
			//check whether event is prior to forecast start, and larger than Mc
			if( event[0] <= forecastStart && event[0] >= 0 && event[1] >= Mc){
				//System.out.println("Seed "+counter++);
				//newEqList.add(event); 

				//add children
				newEqList = getChildren(newEqList, event[0], event[1], (int)event[2], a_sample, p_sample, c_sample, simNumber);
						
			}else{
				//System.out.println("Skipping Seed "+counter++);
			}
			
		}
		
		// sort catalog
		Collections.sort(newEqList, new java.util.Comparator<float[]>() {
		    public int compare(float[] a, float[] b) {
		        return Double.compare(a[0], b[0]);
		    }
		});
		
		// remove events under mc
		for(float[] eq : newEqList){
			if(eq[1] >= Mc);
				finalEqList.add(eq);
		}
		
		//this.eqList = newEqList;
		return finalEqList;
	}
	
	private List<float[]> getChildren(List<float[]> newEqList, float t, float mag, int ngen, 
			double a_sample, double p_sample, double c_sample, int simNumber){//, double forecastStart, double forecastEnd,
			//double a, double b, double p, double c, double alpha, double refMag, double maxMag, int maxGen){
		
		float newMag;
		float newTime;
			
		//calculate productivity of this quake
		double prod = calculateProductivity(t, mag, forecastStart, forecastEnd, a_sample, b, p_sample, c_sample, alpha, Mc);
		long numNew = assignNumberOfOffspring(prod); 
		
//		if(D) System.out.format("Parent Mag: %.2f Time: %5.2f Generation: %d Number of offspring: %d %n", mag, t, (int)ngen, (int)numNew);
		if(numNew > 0 && ngen < maxGenerations){
			//for each new child, assign a magnitude and time
			for(long i=0; i<numNew; i++){
				float[] event = new float[3];		//this must be declared within for block, in order to generate a new address
				
				// assign a magnitude
				newMag = (float) assignMagnitude(b, Mc, maxMagLimit);
				// assign a time
				newTime = (float) assignTime(t, forecastStart, forecastEnd, p_sample, c_sample);

				// add new child to the list
				event[0] = newTime;
				event[1] = newMag;
				event[2] = ngen + 1;
				
				newEqList.add(event);	
			
				// recursively get children of new child
				newEqList = getChildren(newEqList, newTime, newMag, ngen + 1, a_sample, p_sample, c_sample, simNumber);//, forecastStart, forecastEnd, a, b, p, c, alpha, refMag, maxMag, maxGen);
				
			}
		} else if(ngen == maxGenerations) {
			if(D) System.out.println("Sim=" + simNumber + " t=" + t + " has reached " + maxGenerations + " generations. Cutting it short.");
			if(D) System.out.println("n = " + ETAS_StatsCalc.calculateBranchingRatio(a_sample, p_sample, c_sample, alpha, b, forecastEnd, Mc, maxMagLimit)
					+ " a=" + a_sample + " p=" + p_sample + " c=" + c_sample + " al=" + alpha + " b=" + b + " T=" + forecastEnd + " Mc=" + Mc + " Mmax=" + maxMagLimit);
		}
		return newEqList;
		
	}
	
	private double calculateProductivity(float t, float mag, double forecastStart, double forecastEnd,
			double a_sample, double b, double p, double c, double alpha, double refMag){
		
		double unscaledProductivity;
		
		if(t < forecastStart)
			unscaledProductivity = Math.pow(10,a_sample)/(1-p)*( Math.pow(forecastEnd - t + c, 1-p) - Math.pow(forecastStart - t + c, 1-p) );
		else if(t < forecastEnd)
			unscaledProductivity = Math.pow(10,a_sample)/(1-p)*( Math.pow(forecastEnd - t + c, 1-p) - Math.pow(c, 1-p) );
		else
			unscaledProductivity = 0;
		
		double prod = unscaledProductivity * Math.pow(10,(alpha*(mag-refMag)));
		return prod;
	}
	
	
	/*
	 * Returns a sample of a,p,c using the likelihood array provided. This method inverts the cumulative likelihood function, 
	 * so the sum of the likelihood array needs to be normalized to 1, which it should be if using the likelihood calculator 
	 * in this package. 
	 * 
	 * @author Nicholas van der Elst
	 */
	private double[][] sampleParams(int nsamples, double maxMag){
			
		int h = 0, i = 0, j = 0, k = 0;
		int num_ams = ams_vec.length, num_a = a_vec.length, num_p = p_vec.length, num_c = c_vec.length;
		
		double[][] params = new double[nsamples][4];
		
		// generate vector of random numbers
		double[] uRand = new double[nsamples];
		for(int n = 0; n<nsamples; n++){
			uRand[n] = Math.random();
		}
		// sort vector
		Arrays.sort(uRand);
		
		double nbranch;
		double [][][][] likelihoodTrunc = likelihood.clone();
		
		// set up timer/time estimator
		long toc, timeEstimate;
		Stopwatch watch = Stopwatch.createStarted();
		int warnTime = 3;
		String initialMessageString = "Generating parameter sets. ";
		//truncate likelihood based on criticality
		double cumSum = 0;
		for(h = 0; h < num_ams; h++ ){
			for(i = 0; i < num_a; i++ ){
				for(j = 0; j < num_p; j++ ){
					for(k = 0; k < num_c; k++ ){
						nbranch = ETAS_StatsCalc.calculateBranchingRatio(a_vec[i], p_vec[j], c_vec[k], alpha, b, forecastEnd, Mc, maxMag);
						if(nbranch < 1)
							cumSum += likelihoodTrunc[h][i][j][k];
						else
							likelihoodTrunc[h][i][j][k] = 0;

						// run the timer to see how long this is going to take
						toc = watch.elapsed(TimeUnit.SECONDS);
						if(toc > warnTime){
							warnTime += 10;

							timeEstimate = toc * (num_p*num_c*num_ams*num_a)/((h)*(num_p*num_c*num_a) + (i)*(num_c*num_p) + (j)*(num_c) + k);
							System.out.format(initialMessageString + "Approximately %d seconds remaining...\n", (int) ((timeEstimate - toc)));
							initialMessageString = "...";
						}

					}
				}
			}
		}	// else {
		watch.stop();
		//renormalize the random vector to match the likelihood sum
		for(int n = 0; n < nsamples; n++) uRand[n] /= cumSum;
		
		//get cumulative likelihood array
		int n = 0;
		cumSum = 0;
		for(h = 0; h < ams_vec.length; h++ ){
			for(i = 0; i < a_vec.length; i++ ){
				for(j = 0; j < p_vec.length; j++ ){
					for(k = 0; k < c_vec.length; k++ ){
						cumSum += likelihoodTrunc[h][i][j][k];
						while(n < nsamples && cumSum > uRand[n]){
							//found a hit
							params[n] = new double[]{ams_vec[h],a_vec[i], p_vec[j], c_vec[k]};
							n++;
						}
					}
				}
			}
		}	
		
		//shuffle those parameters for a more accurate duration estimate
		java.util.Collections.shuffle(Arrays.asList(params));
		return params;
	}
	
	
	private long assignNumberOfOffspring(double lambda){
		//return Math.round(lambda); //replace with Poisson random number
		return cern.jet.random.tdouble.Poisson.staticNextInt(lambda);
	}
	
	private double assignMagnitude(double b, double minMag, double Mmax){
		double u = Math.random();
		double mag = minMag - Math.log10(1.0 - u*(1.0 - Math.pow(10, -b*(Mmax-minMag))))/b;
		return mag;
	}
	
	private double assignTime(double t0, double tmin, double tmax, double p, double c){
		
		 double u=Math.random();
		 double a1, a2, a3;
		 double t;
		 
		 if(t0 < tmin){
			 a1= Math.pow(tmax - t0 + c, 1d-p);
			 a2= Math.pow(tmin - t0 + c, 1d-p);
		 } else if(t0 < tmax) {
			 a1= Math.pow(tmax - t0 + c, 1d-p);
			 a2= Math.pow(c, 1d-p);
		 } else {
			 a1= Double.NaN;
			 a2= Double.NaN;
		 }
			 
		 a3 = u*a1 + (1d-u)*a2;
		 t = Math.pow(a3, 1d/(1d-p)) - c + t0;

		 return t;
	}
	
	public List<float[]> getETAScatalog(int index){
		return catalogList.get(index); 
		// return eqList;
	}
	
	public int[] get_nEvents(){
		return this.numEventsFinal;
	}
	
	public int get_nEvents(List<float[]> eqList){
		return eqList.size();
	}
	
	public double[] get_maxMag(){
		return this.maxMags;
	}
	
	public double get_maxMag(List<float[]> eqList){ 
		double maxMag = Double.NEGATIVE_INFINITY;
		double mag;

		for(float[] ev : eqList){
			mag = ev[1];
			if( mag > maxMag )
				maxMag = mag;
		}
		return maxMag; 
	}

	public int[] get_nGenerations(){
		return this.numGenerations;
	}
	
	public int get_nGenerations(List<float[]> eqList){
		double maxGen = 0;
		double ngen;

		for(float[] ev : eqList){
			ngen = ev[2];
			if( ngen > maxGen )
				maxGen = ngen;
		}
		return (int)maxGen; 
	}
	
	public double[][][][] getLikelihood(){
		return likelihood;
	}
	
	
	public String printCatalog(int index){
		List<float[]> eqList = getETAScatalog(index);
		
		StringBuffer paragraph = new StringBuffer("Time Mag Gen\n");
		for(float[] eq: eqList){
			 paragraph.append(String.format("%5.2f %5.2f %d %n", eq[0], eq[1], (int)eq[2]));
		}
		return paragraph.toString();
	}
	
	
}
