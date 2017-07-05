package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.dom4j.DocumentException;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.param.HistoricOpenIntervalParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.faultSurface.RuptureSurface;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.analysis.FaultSysSolutionERF_Calc;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.LastEventData;

import com.google.common.collect.Lists;

public class VAHospitalProbCalc {

	/**
	 * This calculates the UCERF3 branch averaged (FM3.1) probabilities of ruptures within
	 * the given radius of a list of California VA hospitals
	 * @param args
	 * @throws IOException
	 * @throws DocumentException
	 */
	public static void main(String[] args) throws IOException, DocumentException {
		File baSol = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/"
				+ "scratch/InversionSolutions/2013_05_10-ucerf3p3-production-10runs_"
				+ "COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip");
		FaultSystemSolution sol = FaultSystemIO.loadSol(baSol);
		LastEventData.populateSubSects(sol.getRupSet().getFaultSectionDataList(), LastEventData.load());
		
		// set ERF params
		double duration = 10d;
		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_PREF_BLEND);
		erf.getTimeSpan().setDuration(duration);
		erf.getTimeSpan().setStartTime(2014);
		erf.setParameter(HistoricOpenIntervalParam.NAME, (double)(2014-1875));
		
		System.out.println("ERF Params:");
		System.out.println("\t"+erf.getAdjustableParameterList().toString().replaceAll(", ", "\n\t"));
		erf.updateForecast();
		
		// list of hospital locations
		List<Location> hospitalLocs = Lists.newArrayList();
		
		double radius = 20d;
		// we use the fault trace as a rough check first - if no site is within this buffer of the
		// trace then we can safely skip the full check
		double traceCheckBuffer = radius*5d;
		
		hospitalLocs.add(new Location(34.05266, -118.45256)); // LA
		hospitalLocs.add(new Location(33.77662, -118.11898)); // Long Beach
		hospitalLocs.add(new Location(32.87495, -117.23228)); // San Diego
		hospitalLocs.add(new Location(37.40511, -122.14018)); // Palo Alto
		hospitalLocs.add(new Location(37.78218, -122.50432)); // San Francisco
		
		// find all sources within radius of at least one of these
		List<ProbEqkSource> sourcesWithin = Lists.newArrayList();
		
		int faultSources = 0;
		int griddedSources = 0;
		
		for (ProbEqkSource source : erf) {
			RuptureSurface surf = source.getSourceSurface();
			LocationList sourceLocs = surf.getEvenlyDiscritizedListOfLocsOnSurface();
			LocationList traceLocs = surf.getPerimeter();
			boolean inside = false;
			hospitalLoop:
			for (Location loc : hospitalLocs) {
				// first check traceLocs
				boolean doFull = false;
				for (Location traceLoc : traceLocs) {
					double dist = LocationUtils.horzDistanceFast(loc, traceLoc);
					if (dist <= radius) {
						inside = true;
						break hospitalLoop;
					} else if (dist <= traceCheckBuffer) {
						doFull = true;
					}
				}
				// now check all source locations if applicable
				if (doFull) {
					for (Location sourceLoc : sourceLocs) {
						double dist = LocationUtils.horzDistanceFast(loc, sourceLoc);
						if (dist <= radius) {
							inside = true;
							break hospitalLoop;
						}
					}
				}
			}
			if (inside) {
				sourcesWithin.add(source);
				if (surf.isPointSurface())
					griddedSources++;
				else
					faultSources++;
			}
		}
		
		System.out.println(sourcesWithin.size()+"/"+erf.getNumSources()+" sources are within "+radius);
		System.out.println("Fault sources: "+faultSources);
		System.out.println("Gridded sources: "+griddedSources);
		
		// start at 6.05 so that cumulative probabilities are correct >=6
		EvenlyDiscretizedFunc mpd = new EvenlyDiscretizedFunc(6.05, 30, 0.1d);
		// we will calculate all probabilities via: totProb = 1 - (1 - prob1)*(1 - prob2)*...*(1 - probN)
		// must first gather all probabilities
		List<List<Double>> mfdProbs = Lists.newArrayList();
		for (int i=0; i<mpd.size(); i++)
			mfdProbs.add(new ArrayList<Double>());
		double minMag = mpd.getMinX()-0.5*mpd.getDelta();
		for (ProbEqkSource source : sourcesWithin) {
			for (ProbEqkRupture rup : source) {
				double mag = rup.getMag();
				if (mag < minMag)
					continue;
				// this is safe because the maximum magnitude is higher than all UCERF3 mags, and we already
				// checked if it was too small
				int index = mpd.getClosestXIndex(mag);
				// cumulative, so add to each mag up to the actual mag bin
				for (int i=0; i<=index; i++)
					mfdProbs.get(i).add(rup.getProbability());
			}
		}
		// now sum the probabilities
		for (int i=0; i<mpd.size(); i++)
			mpd.set(i, FaultSysSolutionERF_Calc.calcSummedProbs(mfdProbs.get(i)));
		
		System.out.println((int)duration+"yr Probabilities:");
		for (int i=0; i<mpd.size(); i++) {
			// adjust mag for cumulative
			double mag = mpd.getX(i)-0.5*mpd.getDelta();
			double prob = mpd.getY(i);
			System.out.println("\tM>="+(float)mag+": "+(float)prob);
		}
	}

}
