package scratch.aftershockStatisticsETAS;


import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;

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

	double[] a_vec, p_vec, c_vec;
	double[][][] likelihood;
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
	private long[] numEventsFinal;
	private int[] numGenerations;
	private boolean D = false;	//debug flag
	private List<List<double[]>> catalogList;	//list of catalogs
	
	public ETAScatalog(double[] a_vec, double[] p_vec, double[] c_vec, double[][][] likelihood, double alpha, double b, double refMag,
			ObsEqkRupture mainshock, ObsEqkRupList aftershocks,
			double dataStart, double dataEnd, double forecastStart, double forecastEnd, double Mc, double maxMag, int maxGenerations, int nSims ){

		this(Double.NaN, a_vec, p_vec, c_vec, likelihood, alpha, b, refMag, 
				mainshock, aftershocks, dataStart, dataEnd, forecastStart, forecastEnd, Mc, 9.5, 100, nSims);
	}
	
	public ETAScatalog(double a_base, double[] a_vec, double[] p_vec, double[] c_vec, double[][][] likelihood, double alpha, double b, double refMag,
			ObsEqkRupture mainshock, ObsEqkRupList aftershocks,
			double dataStart, double dataEnd, double forecastStart, double forecastEnd, double Mc, double maxMag, int maxGenerations, int nSims) {
	
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
		
		List<double[]> newEqList = new ArrayList<double[]>();	//catalog containing list of {time, mag, gen}
		List<List<double[]>> catalogList = new ArrayList<List<double[]>>(); //list of catalogs
		
		double[] maxMags = new double[nSims];
		long[] nEvents = new long[nSims];
		int[] nGens = new int[nSims];
		
		double p_base;
		if(p_vec.length % 2 == 1)
			p_base = p_vec[(int)(p_vec.length - 1)/2];
		else{
			p_base = 0.5 * (p_vec[(int)(p_vec.length/2)] + p_vec[(int)(p_vec.length/2)]);
		}
		
		double c_base;
		if(c_vec.length % 2 == 1)
			c_base = c_vec[(int)(c_vec.length - 1)/2];
		else{
			c_base = 0.5 * (c_vec[(int)(c_vec.length/2)] + c_vec[(int)(c_vec.length/2)]);
		}
		
		
		long tic = System.currentTimeMillis();
		long toc;
		boolean timeWarning = false;
		System.out.println("Calculating " + nSims + " " + (int)(forecastEnd - forecastStart) + "-day ETAS catalogs...");
		for(int i = 0; i < nSims ; i++){
			toc = System.currentTimeMillis();
			if (timeWarning == false && toc - tic > 3000){
				System.out.println("This might take a while. Probably around " + (int)((double)(toc-tic)/i * nSims/1000) + " seconds.");
				timeWarning = true;
			}
			
			int count = 0;
			double n = 2;
			double[] params = new double[3];
			double prodCorr;
			if(Mc != refMag)
				prodCorr = Math.log10((maxMag - refMag)/(maxMag - Mc));
			else
				prodCorr = 0;
			
			while(n > 1){
				//sample parameter values from parameter/likelihood grids
				params = sampleParams();
				params[0] += prodCorr;
				
				//check for supercritical parameters
				n = b * Math.log(10) * (maxMag - Mc) * Math.pow(10, params[0])/(1-params[1]) * ( Math.pow(forecastEnd + params[2], 1-params[1]) - Math.pow(params[2], 1-params[1]) );
				if(++count > 100){
					System.out.println("Found 100 combinations of supercritical params in a row. There's a problem with the generic model.");
				}
				
			}
			
			double a_sample, p_sample, c_sample;
			ObsEqkRupture simulationMainshock = (ObsEqkRupture) mainshock.clone();
			
			if(Double.isNaN(a_base) || Math.abs(a_base - 0) < 1E-6) {
				//not a generic model, use a randomly sampled a for mainshock and aftershock productivity
				a_sample = params[0];
				p_sample = params[1];
				c_sample = params[2];
				
			} else {
				// generic forecast is generating a random mainshock magnitude but leaving all parameters constant...
				// a_base has been supplied. probably a generic model, use a_base for a_sample, and adjust mainshock magnitude. leave p and c alone. Too easy to get supercritical sequence.
				a_sample = a_base + prodCorr;
				p_sample = p_base;
				c_sample = c_base;
				
				simulationMainshock.setMag(mainshock.getMag() + (params[0] - (a_base + prodCorr)));
			}
			//check for supercritical parameters
			n = Math.log(10) * (maxMag - Mc) * Math.pow(10, a_sample)/(1-p_sample) * ( Math.pow(forecastEnd + c_sample, 1-p_sample) - Math.pow(c_sample, 1-p_sample) );
			
			// when running the generic model, we want to get variability in mainshock productivity, without transmitting that 
			// value to the offspring. How do we do this for the two model types?
			
			
			// if the magnitude of completeness (min simulated mag) is not the same as magref, the statistics will be wrong. Make an approximate correction
			// by adjusting k, a by the difference in magnitude range
			
//			System.out.println(params[0] +" "+ a_sample + " " + p_sample + " " + c_sample + " "+ n);
				
			
			
			newEqList = getNewETAScatalog(simulationMainshock, aftershocks, a_sample, p_sample, c_sample);
			maxMags[i] = get_maxMag(newEqList);
			nEvents[i] = get_nEvents(newEqList);
			nGens[i] = get_nGenerations(newEqList);
			catalogList.add(i, newEqList);
		}
		toc = System.currentTimeMillis();
		System.out.println("Finished. It took " + (toc-tic)/1000 + " seconds.");
		
		//this.eqList = getLastETAScatalog();
		this.catalogList = catalogList;
		this.maxMags = maxMags;
		this.numEventsFinal = nEvents;
		this.numGenerations = nGens;
	}
		
	public List<double[]> getNewETAScatalog(ObsEqkRupture mainshock, ObsEqkRupList aftershocks, double a_sample, double p_sample, double c_sample){
		
		//extract magnitudes and times from supplied Eqk rupture objects to make catalog (combine MS and AS's)
		List<double[]> newEqList = new ArrayList<double[]>();
		List<double[]> finalEqList = new ArrayList<double[]>();
		
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
			double[] event = new double[3];
			event[0] = (rup.getOriginTime() - t0)/ETAS_StatsCalc.MILLISEC_PER_DAY;	//elapsed time in days
			event[1] = rup.getMag();	
			event[2] = 0;	//generation number
			
			//check whether event is prior to forecast start, and larger than Mc
			if( event[0] <= forecastStart && event[0] >= 0 && event[1] >= Mc){
				//System.out.println("Seed "+counter++);
				//newEqList.add(event); 

				//add children
				newEqList = getChildren(newEqList, event[0], event[1], (int)event[2], a_sample, p_sample, c_sample);
						//forecastStart, forecastEnd, a, b, p, c, alpha, refMag, maxMagLimit, maxGenerations);
			}else{
				//System.out.println("Skipping Seed "+counter++);
			}
			
		}
		
		// sort catalog
		Collections.sort(newEqList, new java.util.Comparator<double[]>() {
		    public int compare(double[] a, double[] b) {
		        return Double.compare(a[0], b[0]);
		    }
		});
		
		// remove events under mc
		for(double[] eq : newEqList){
			if(eq[1] >= Mc);
				finalEqList.add(eq);
		}
		
		//this.eqList = newEqList;
		return finalEqList;
	}
	
	private List<double[]> getChildren(List<double[]> newEqList, double t, double mag, int ngen, 
			double a_sample, double p_sample, double c_sample){//, double forecastStart, double forecastEnd,
			//double a, double b, double p, double c, double alpha, double refMag, double maxMag, int maxGen){
		
		double newMag;
		double newTime;
		
		
		//calculate productivity of this quake
		double prod = calculateProductivity(t, mag, forecastStart, forecastEnd, a_sample, b, p_sample, c_sample, alpha, Mc);
		long numNew = assignNumberOfOffspring(prod); 
		
		if(D) System.out.format("Parent Mag: %.2f Time: %5.2f Generation: %d Number of offspring: %d %n", mag, t, (int)ngen, (int)numNew);
		if(numNew > 0 && ngen < maxGenerations){
			//for each new child, assign a magnitude and time
			for(long i=0; i<numNew; i++){
				double[] event = new double[3];		//this SOB must be declared within for block, in order to generate a new address
				
				// assign a magnitude
				newMag = assignMagnitude(b, Mc, maxMagLimit);
				// assign a time
				newTime = assignTime(t, forecastStart, forecastEnd, p_sample, c_sample);

				// add new child to the list
				event[0] = newTime;
				event[1] = newMag;
				event[2] = ngen + 1;
				
				newEqList.add(event);	
			
				// recursively get children of new child
//				newEqList.addAll(getChildren(newTime, newMag, ngen + 1, forecastStart, forecastEnd, a, b, p, c, alpha, refMag, maxMag, maxGen));
				newEqList = getChildren(newEqList, newTime, newMag, ngen + 1, a_sample, p_sample, c_sample);//, forecastStart, forecastEnd, a, b, p, c, alpha, refMag, maxMag, maxGen);
				
			}
		} else if(ngen == maxGenerations) {
			
			System.out.println("Simulation has reached " + maxGenerations + " generations. Cutting it short.");
			
		}
		return newEqList;
		
	}
	
	private double calculateProductivity(double t, double mag, double forecastStart, double forecastEnd,
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
	
	
//	private double sampleA_value(double a,double a_sigma){
//		double a_sample;
//		
//		if(a_sigma > 0)
//			a_sample = cern.jet.random.tdouble.Normal.staticNextDouble(a, a_sigma);
//		else
//			a_sample = a;
//		
//		return a_sample;
//	}
	
	/*
	 * Returns a sample of a,p,c using the likelihood array provided. This method inverts the cumulative likelihood function, 
	 * so the sum of the likelihood array needs to be normalized to 1, which it should be if using the likelihood calculator 
	 * in this package. 
	 * 
	 * @author Nicholas van der Elst
	 */
	private double[] sampleParams(){
			
		//double a_sample = cern.jet.random.tdouble.Normal.staticNextDouble(a, a_sigma);
		double cumSum;
		int i = 0, j = 0, k = 0;
		
		double uRand = Math.random();
		
		if(uRand < 0.5){
			cumSum = 0;
			//get cumulative likelihood array
			outerloop:
			for(i = 0; i < a_vec.length; i++ ){
				for(j = 0; j < p_vec.length; j++ ){
					for(k = 0; k < c_vec.length; k++ ){
						cumSum += likelihood[i][j][k];
						if(cumSum > uRand){
//							System.out.println(i +" "+ j +" "+ k +" "+ uRand +" "+ cumSum);
							break outerloop;
						}
					}
				}
			}
		} else {
			cumSum = 1;
			//get cumulative likelihood array
			outerloop:
			for(i = a_vec.length-1; i >= 0; i-- ){
				for(j = p_vec.length-1; j >= 0; j-- ){
					for(k = c_vec.length-1; k >= 0; k-- ){
						cumSum -= likelihood[i][j][k];
						if(cumSum <= uRand){
//							System.out.println(i +" "+ j +" "+ k +" "+ uRand +" "+ cumSum);
							break outerloop;
						}
					}
				}
			}
		}
		
//		double a_sample = a_vec[(int) a_vec.length/2];
//		double p_sample = p_vec[(int) p_vec.length/2];
//		double c_sample = c_vec[(int) c_vec.length/2];;
//		double[] params = new double[]{a_sample, p_sample, c_sample};
		
//		System.out.println(i + " "+ j +" "+ k +" "+ cumSum);
		double[] params = new double[]{a_vec[i], p_vec[j], c_vec[k]};
		
		return params;
	}
	
	
	private long assignNumberOfOffspring(double lambda){
		//return Math.round(lambda); //replace with Poisson random number
		return cern.jet.random.tdouble.Poisson.staticNextInt(lambda);
	}
	
	private double assignMagnitude(double b, double minMag, double Mmax){
		double u=Math.random();
		double mag = 1/b*(minMag-Math.log10(1-u*(1-Math.pow(10, -b*(Mmax-minMag)))));
		return mag;
	}
	
	private double assignTime(double t0, double tmin, double tmax, double p, double c){
		
		 double u=Math.random();
		 double a1, a2, a3;
		 double t;
		 
		 if(t0 < tmin){
			 a1= Math.pow(tmax - t0 + c, 1-p);
			 a2= Math.pow(tmin - t0 + c, 1-p);
		 } else if(t0 < tmax) {
			 a1= Math.pow(tmax - t0 + c, 1-p);
			 a2= Math.pow(c, 1-p);
		 } else {
			 a1= Double.NaN;
			 a2= Double.NaN;
		 }
			 
		 a3 = u*a1 + (1-u)*a2;
		 t = Math.pow(a3, 1./(1-p)) - c + t0;

		 return t;
	}
	
	public List<double[]> getETAScatalog(int index){
		return catalogList.get(index); 
		// return eqList;
	}
	
	public long[] get_nEvents(){
		return this.numEventsFinal;
	}
	
	public long get_nEvents(List<double[]> eqList){
		return eqList.size();
	}
	
	public double[] get_maxMag(){
		return this.maxMags;
	}
	
	public double get_maxMag(List<double[]> eqList){ 
		double maxMag = Double.NEGATIVE_INFINITY;
		double mag;

		for(double[] ev : eqList){
			mag = ev[1];
			if( mag > maxMag )
				maxMag = mag;
		}
		return maxMag; 
	}

	public int[] get_nGenerations(){
		return this.numGenerations;
	}
	
	public int get_nGenerations(List<double[]> eqList){
		double maxGen = 0;
		double ngen;

		for(double[] ev : eqList){
			ngen = ev[2];
			if( ngen > maxGen )
				maxGen = ngen;
		}
		return (int)maxGen; 
	}
	
	public double[][][] getLikelihood(){
		return likelihood;
	}
	
	
	public String printCatalog(int index){
		List<double[]> eqList = getETAScatalog(index);
		
		StringBuffer paragraph = new StringBuffer("Time Mag Gen\n");
		for(double[] eq: eqList){
			 paragraph.append(String.format("%5.2f %5.2f %d %n", eq[0], eq[1], (int)eq[2]));
		}
		return paragraph.toString();
	}
	
	
}
