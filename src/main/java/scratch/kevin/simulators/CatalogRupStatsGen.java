package scratch.kevin.simulators;

import java.io.IOException;

import org.opensha.sha.simulators.RSQSimEvent;

import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class CatalogRupStatsGen {

	public static void main(String[] args) throws IOException {
		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance();
		
		int rupCount = 0;
		double minMag = Double.POSITIVE_INFINITY;
		double maxMag = Double.NEGATIVE_INFINITY;
		
		double[] mags = { 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8 };
		int[] magCounts = new int[mags.length];
		
		double firstTimeYears = Double.NaN;
		double lastTimeYears = Double.NaN;
		
		for (RSQSimEvent event : catalog.loader().iterable()) {
			if (rupCount == 0)
				firstTimeYears = event.getTimeInYears();
			rupCount++;
			if (rupCount % 100000 == 0)
				System.out.println("Processed "+rupCount+" events ("+(float)(event.getTimeInYears()-firstTimeYears)+" years)");
			double mag = event.getMagnitude();
			minMag = Math.min(minMag, mag);
			maxMag = Math.max(maxMag, mag);
			for (int m=0; m<mags.length; m++)
				if (mag >= mags[m])
					magCounts[m]++;
			
			lastTimeYears = event.getTimeInYears();
		}
		
		System.out.println("Time range: "+(float)firstTimeYears+" => "+(float)lastTimeYears+" years");
		double duration = lastTimeYears - firstTimeYears;
		System.out.println("Duration: "+(float)(duration)+" years");
		System.out.println("Num Events: "+rupCount);
		System.out.println("Overall rate: "+(float)(rupCount/duration)+" events/yr");
		System.out.println("Minimum magnitude: "+minMag);
		System.out.println("Maximum magnitude: "+maxMag);
		for (int m=0; m<mags.length; m++) {
			System.out.println("M>="+(float)mags[m]);
			System.out.println("\t"+magCounts[m]+" events");
			System.out.println("\t"+(float)(magCounts[m]/duration)+" events/yr");
		}
	}

}
