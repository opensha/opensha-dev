package scratch.kevin;

import org.opensha.commons.geo.Location;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.RuptureSurface;

public class ChristineBroadbandRupInfoFetcher {

	public static void main(String[] args) {
		// 85      S. San Andreas;PK+CH+CC+BB+NM+SM
		// 158	Elysian Park (Upper)
		// 112	San Jacinto;SBV+SJV+A+C
		// 39	N. San Andreas;SAO+SAN+SAP+SAS
		// 28	Hayward-Rodgers Creek;HN+HS
		int[] sourceIDs = { 85, 158, 112, 39, 28 };
		
		MeanUCERF2 erf = new MeanUCERF2();
		erf.updateForecast();
		
		for (int sourceID : sourceIDs) {
			ProbEqkSource source = erf.getSource(sourceID);
			RuptureSurface surf = source.getSourceSurface();
			FaultTrace trace = surf.getUpperEdge();
			
			double modalMag = 0d;
			double maxProb = 0d;
			for (ProbEqkRupture rup : source) {
				if (rup.getProbability() > maxProb) {
					maxProb = rup.getProbability();
					modalMag = rup.getMag();
				}
			}
			
			System.out.println("Source Name: "+source.getName());
			System.out.println("Modal Mag: "+(float)modalMag);
			System.out.println("Upper Depth (km): "+(float)surf.getAveRupTopDepth());
			System.out.println("Width (km): "+(float)surf.getAveWidth());
			System.out.println("Dip: "+(float)surf.getAveDip());
			System.out.println("Upper Edge:");
			for (Location loc : trace)
				System.out.println("\t"+loc);
			System.out.println();
		}
	}

}
