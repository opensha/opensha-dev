package scratch.kevin;

import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CompletableFuture;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.param.Parameter;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.sourceFilters.SourceFilterManager;
import org.opensha.sha.calc.sourceFilters.SourceFilters;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.erf.NSHM23_WUS_BranchAveragedERF;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;

public class ThreadedHazardCurveCalcExample {

	public static void main(String[] args) throws IOException {
		// ERF
		NSHM23_WUS_BranchAveragedERF erf = new NSHM23_WUS_BranchAveragedERF();
		erf.getTimeSpan().setDuration(1d);
		erf.updateForecast();
		
		// GMM, can't share an instance across multiple threads so we keep the reference to build instances
		AttenRelRef gmmRef = AttenRelRef.ASK_2014;
		
		// sites, I'll use a gridded region to build the site list for this example
		GriddedRegion sitesGridded = new GriddedRegion(NSHM23_RegionLoader.loadFullConterminousWUS(), 1d, GriddedRegion.ANCHOR_0_0);
		// create the site list
		List<Site> sites = new ArrayList<>();
		// need one GMM in order to get the site parameter list
		ScalarIMR gmm0 = gmmRef.get();
		for (Location loc : sitesGridded.getNodeList()) {
			Site site = new Site(loc);
			for (Parameter<?> param : gmm0.getSiteParams()) {
				// need to clone the site parameter so that it's not shared across multiple sites (in case you want
				// different site data per site)
				site.addParameter((Parameter<?>)param.clone());
			}
			sites.add(site);
		}
		
		// source filters to use for the calculator
		SourceFilterManager sourceFilters = new SourceFilterManager(SourceFilters.FIXED_DIST_CUTOFF);
		
		// this will cache and keep track of our GMM instances
		ArrayDeque<ScalarIMR> gmmDeque = new ArrayDeque<>();
		// and our hazard curve calculator instances
		ArrayDeque<HazardCurveCalculator> calcDeque = new ArrayDeque<>();
		
		// need x values
		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(PGA_Param.NAME);
		// need them in ln spacing
		DiscretizedFunc logXVals = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<xVals.size(); i++)
			logXVals.set(Math.log(xVals.getX(i)), 0d);
		
		// now spin it up in parallel
		List<CompletableFuture<DiscretizedFunc>> futures = new ArrayList<>(sites.size());
		
		System.out.println("Calculating "+sites.size()+" curves");
		for (Site site : sites) {
			futures.add(CompletableFuture.supplyAsync(()-> {
				// this runs in parallel
				
				// get GMM instance for this thread
				ScalarIMR gmm = null;
				synchronized (gmmDeque) {
					if (!gmmDeque.isEmpty())
						gmm = gmmDeque.pop();
				}
				if (gmm == null) {
					gmm = gmmRef.get();
					// could set any custom params here
				}
				// get hazard curve calculator instance for this thread
				HazardCurveCalculator calc = null;
				synchronized (calcDeque) {
					if (!calcDeque.isEmpty())
						calc = calcDeque.pop();
				}
				if (calc == null) {
					calc = new HazardCurveCalculator(sourceFilters);
				}
				
				DiscretizedFunc logCurve = logXVals.deepClone();
				calc.getHazardCurve(logCurve, site, gmm, erf);
				
				// return those instances for reuse by later threads
				synchronized (gmmDeque) {
					gmmDeque.push(gmm);
				}
				synchronized (calcDeque) {
					calcDeque.push(calc);
				}
				
				// return linear curve
				DiscretizedFunc linearCurve = xVals.deepClone();
				for (int i=0; i<linearCurve.size(); i++)
					linearCurve.set(i, logCurve.getY(i));
				return logCurve;
			}));
		}
		
		// now all of the tasks have been submitted, need to "join" them (wait for them to finish)
		for (int s=0; s<sites.size(); s++) {
			Site site = sites.get(s);
			DiscretizedFunc curve = futures.get(s).join();
			// do what you need to do here
			
			// this print helps keep track of progress
			System.out.print(".");
			if (s > 0 && (s % 100 == 0 || s == sites.size()-1))
				System.out.println(" "+s+"/"+sites.size());
		}
	}

}
