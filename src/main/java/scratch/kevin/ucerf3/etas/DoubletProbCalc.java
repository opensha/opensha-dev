package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.faultSurface.RuptureSurface;

import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO.ETAS_Catalog;
import scratch.UCERF3.erf.ETAS.launcher.ETAS_Config;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;

public class DoubletProbCalc {

	public static void main(String[] args) throws IOException {
		File dir = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2021_05_17-Start1919_100yr_kCOV1p5_MaxPtSrcM6_Scale1p29_Spontaneous_HistCatalog");
		File resultsFile = new File(dir, "results_m5.bin");
		File configFile = new File(dir, "config.json");
		ETAS_Config config = ETAS_Config.readJSON(configFile);
		
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(config.getFSS_File());
		
		double minPrimaryMag = 7d;
		double minAftershockMag = 7d;
		double maxSurfDist = 200d;
		double maxDaysBetween = 7d;
		double probDuration = 30d;
		
		long maxMillisBetween = (long)(maxDaysBetween*ProbabilityModelsCalc.MILLISEC_PER_DAY);
		
		long startMillis = config.getSimulationStartTimeMillis();
		long endMillis = startMillis + (long)(config.getDuration()*ProbabilityModelsCalc.MILLISEC_PER_YEAR);
		long latestPossibleFirstEvent = endMillis - maxMillisBetween;
		
		Region[] regions = {
				null,
				new CaliforniaRegions.RELM_NOCAL(),
				new CaliforniaRegions.RELM_SOCAL()
		};
		
		double minMag = Math.min(minPrimaryMag, minAftershockMag);
		
		List<ETAS_Catalog> catalogs = ETAS_CatalogIO.loadCatalogsBinary(resultsFile, minMag);
		
		double totCatalogTime = catalogs.size()*(latestPossibleFirstEvent - startMillis)/ProbabilityModelsCalc.MILLISEC_PER_YEAR;
		
		int[] regionalNumPrimaries = new int[regions.length];
		int[] regionalNumDublets = new int[regions.length];
		
		int numPrimaries = 0;
		int numDoublets = 0;
		
		DecimalFormat pDF = new DecimalFormat("0.00%");
		
		for (ETAS_Catalog catalog : catalogs) {
			for (int r=0; r<catalog.size(); r++) {
				ETAS_EqkRupture primary = catalog.get(r);
				if (primary.getOriginTime() > latestPossibleFirstEvent)
					break;
				if (primary.getMag() < minPrimaryMag)
					continue;
				RuptureSurface primarySurf = primary.getFSSIndex() >= 0 ?
						rupSet.getSurfaceForRupture(primary.getFSSIndex(), 2d) : null;
				// this is a matching primary event
				boolean doublet = false;
				long maxTime = primary.getOriginTime() + maxMillisBetween;
				boolean[] regMatches = new boolean[regions.length];
				for (int j=0; j<regions.length; j++) {
					if (regions[j] == null) {
						regMatches[j] = true;
					} else {
						// match if either events match 
						List<Location> testLocs = new ArrayList<>();
						testLocs.add(primary.getHypocenterLocation());
						if (primarySurf != null)
							testLocs.addAll(primarySurf.getEvenlyDiscritizedListOfLocsOnSurface());
						for (Location loc : testLocs) {
							if (regions[j].contains(loc)) {
								regMatches[j] = true;
								break;
							}
						}
					}
				}
				for (int i=r+1; i<catalog.size(); i++) {
					ETAS_EqkRupture secondary = catalog.get(i);
					if (secondary.getOriginTime() > maxTime)
						break;
					if (secondary.getMag() >= minAftershockMag) {
						// within time range and magnitude range, check distance
						RuptureSurface secondarySurf = secondary.getFSSIndex() >= 0 ?
								rupSet.getSurfaceForRupture(secondary.getFSSIndex(), 2d) : null;
						double dist = LocationUtils.horzDistance(primary.getHypocenterLocation(), secondary.getHypocenterLocation());
						if (primarySurf != null && secondarySurf != null)
							dist = Math.min(dist, primarySurf.getMinDistance(secondarySurf));
						if (primarySurf != null)
							dist = Math.min(dist, primarySurf.getDistanceJB(secondary.getHypocenterLocation()));
						if (secondarySurf != null)
							dist = Math.min(dist, secondarySurf.getDistanceJB(primary.getHypocenterLocation()));
						if (dist < maxSurfDist) {
							// doublet!
							doublet = true;
							for (int j=0; j<regions.length; j++) {
								if (!regMatches[j]) {
									// match if either events match
									List<Location> testLocs = new ArrayList<>();
									testLocs.add(secondary.getHypocenterLocation());
									if (secondarySurf != null)
										testLocs.addAll(secondarySurf.getEvenlyDiscritizedListOfLocsOnSurface());
									for (Location loc : testLocs) {
										if (regions[j].contains(loc)) {
											regMatches[j] = true;
											break;
										}
									}
								}
							}
						}
					}
				}
				numPrimaries++;
				if (doublet) {
					numDoublets++;
					System.out.println("Doublet "+numDoublets+" found:\t"+numDoublets+"/"+numPrimaries
							+" ("+pDF.format((double)numDoublets/(double)numPrimaries)+")");
				}
				for (int i=0; i<regMatches.length; i++) {
					if (regMatches[i]) {
						regionalNumPrimaries[i]++;
						if (doublet)
							regionalNumDublets[i]++;
					}
				}
			}
		}
		
		for (int r=0; r<regions.length; r++) {
			System.out.println();
			if (regions[r] == null)
				System.out.println("Model-wide:");
			else
				System.out.println(regions[r].getName()+":");
			double primaryRate = regionalNumPrimaries[r]/totCatalogTime;
			double primaryProb = 1d-Math.exp(-primaryRate*probDuration);
			System.out.println("\tFound "+regionalNumPrimaries[r]+" primary events;\t"
					+(float)primaryRate+" /yr;\t"+pDF.format(primaryProb)+" in "+(int)probDuration+" yrs");
			double doubletConditionalProb = (double)regionalNumDublets[r]/(double)regionalNumPrimaries[r];
			double doubletRate = regionalNumDublets[r]/totCatalogTime;
			double doubletProb = 1d-Math.exp(-doubletRate*probDuration);
			System.out.println("\tFound "+regionalNumDublets[r]+" doublets sequences;\t"
					+(float)doubletRate+" /yr;\t"+pDF.format(doubletProb)+" in "+(int)probDuration+" yrs;\t"
					+pDF.format(doubletConditionalProb)+" conditional prob"); 
		}
	}

}
