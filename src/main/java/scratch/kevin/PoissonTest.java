package scratch.kevin;

import org.apache.commons.math3.distribution.PoissonDistribution;

public class PoissonTest {
	
	public static void main(String[] args) {
		double mean = 3;
		double testMax = 8;
		
		PoissonDistribution p = new PoissonDistribution(mean);
		
		System.out.println("PROBABLITIES!");
		for (int i=0; i<=testMax; i++)
			System.out.println(i+": "+p.probability(i));
		
		System.out.println("\nCUMULATIVE PROBABLITIES!");
		for (int i=0; i<=testMax; i++)
			System.out.println(i+": "+p.cumulativeProbability(i));
		
		System.out.println("\nMY PROBABILITIES!");
		for (int i=0; i<=testMax; i++) {
			double prob;
			if (i <= mean)
				prob = p.cumulativeProbability(i);
			else
				prob = 1 - p.cumulativeProbability(i); 
			System.out.println(i+": "+prob);
		}
	}

}
