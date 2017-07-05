package scratch.kevin.simulators;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.data.CSVFile;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.iden.ElementMagRangeDescription;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.parsers.EQSIMv06FileReader;
import org.opensha.sha.simulators.utils.General_EQSIM_Tools;

import scratch.kevin.simulators.catBuild.RandomCatalogBuilder;
import scratch.kevin.simulators.dists.RandomDistType;

import com.google.common.collect.Lists;

public class Within5YrsStatisticsGen {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		File dir = new File("/home/kevin/Simulators");
		File geomFile = new File(dir, "ALLCAL2_1-7-11_Geometry.dat");
		System.out.println("Loading geometry...");
		General_EQSIM_Tools tools = new General_EQSIM_Tools(geomFile);
//		File eventFile = new File(dir, "eqs.ALLCAL2_RSQSim_sigma0.5-5_b=0.015.barall");
		File eventFile = new File(dir, "eqs.ALLCAL2_RSQSim_sigma0.5-5_b=0.015.long.barall");
		System.out.println("Loading events...");
		List<? extends SimulatorEvent> events = EQSIMv06FileReader.readEventsFile(eventFile, tools.getElementsList());
		
		List<RuptureIdentifier> rupIdens = Lists.newArrayList();
		
//		ElementMagRangeDescription cholameIden = new ElementMagRangeDescription(
//				ElementMagRangeDescription.SAF_CHOLAME_ELEMENT_ID, 7d, 10d);
		
		ElementMagRangeDescription carrizoIden = new ElementMagRangeDescription("SAF Carrizo 7+",
				ElementMagRangeDescription.SAF_CARRIZO_ELEMENT_ID, 7d, 10d);
		rupIdens.add(carrizoIden);
		
		ElementMagRangeDescription garlockIden = new ElementMagRangeDescription("SAF Garlock 7+", 
				ElementMagRangeDescription.GARLOCK_WEST_ELEMENT_ID, 7d, 10d);
		rupIdens.add(garlockIden);
		
		ElementMagRangeDescription mojaveIden = new ElementMagRangeDescription("SAF Mojave 7+",
				ElementMagRangeDescription.SAF_MOJAVE_ELEMENT_ID, 7d, 10d);
		rupIdens.add(mojaveIden);
		
		ElementMagRangeDescription coachellaIden = new ElementMagRangeDescription("SAF Coachella 7+",
				ElementMagRangeDescription.SAF_COACHELLA_ELEMENT_ID, 7d, 10d);
		rupIdens.add(coachellaIden);
		
		ElementMagRangeDescription sanJacintoIden = new ElementMagRangeDescription("San Jacinto 7+",
				ElementMagRangeDescription.SAN_JACINTO__ELEMENT_ID, 7d, 10d);
		rupIdens.add(sanJacintoIden);
		
		boolean[] randoms = { false, true };
		
		double[] times = { 1, 5, 10 };
		
		for (boolean randomized : randoms) {
			List<? extends SimulatorEvent> theEvents;
			if (randomized) {
				System.out.println("\n\n******* RANDOMIZED ******\n");
				theEvents = RandomCatalogBuilder.getRandomResampledCatalog(events, rupIdens, RandomDistType.ACTUAL, true);
			} else {
				theEvents = events;
			}
			for (double years : times) {
				System.out.println("****** "+(float)years+" years ******");
//				System.out.println("Section\t\t")
				CSVFile<String> probCSV = new CSVFile<String>(true);
				CSVFile<String> percentCSV = new CSVFile<String>(true);
				probCSV.addLine("Section", "Corupture", "Before", "After", "Total", "Indep Prob");
				percentCSV.addLine("Section", "Corupture", "Before", "After", "Total");
				
				addStats(theEvents, coachellaIden, "Coachella", mojaveIden, "Mojave", years, probCSV, percentCSV);
				addStats(theEvents, carrizoIden, "Carrizo", mojaveIden, "Mojave", years, probCSV, percentCSV);
				addStats(theEvents, garlockIden, "Garlock", mojaveIden, "Mojave", years, probCSV, percentCSV);
				addStats(theEvents, sanJacintoIden, "San Jacinto", mojaveIden, "Mojave", years, probCSV, percentCSV);
				
				System.out.println("Probability given Mojave (within "+(float)years+" years)");
				probCSV.printPretty(" | ");
				System.out.println("Percent within "+(float)years+" years of Mojave");
				percentCSV.printPretty(" | ");
				System.out.println("***********************\n");
			}
		}
	}
	
	private static void addStats(List<? extends SimulatorEvent> events,
			RuptureIdentifier targetIden, String targetName,
			RuptureIdentifier givenIden, String givenName, double years,
			CSVFile<String> probCSV, CSVFile<String> percentCSV) {
		List<? extends SimulatorEvent> targetMatches = targetIden.getMatches(events);
		List<? extends SimulatorEvent> givenMatches = givenIden.getMatches(events);
		
		System.out.println("target="+targetMatches.size()+"\tgiven="+givenMatches.size()
				+"\tduration="+General_EQSIM_Tools.getSimulationDuration(events));
		
		int numBefore = 0;
		int numCorupture = 0;
		int numAfter = 0;
		int numTotal = 0;
		
		// find prob that target occurs before/after given
		// only first match in each category is counted
		for (SimulatorEvent given : givenMatches) {
			double givenTime = given.getTimeInYears();
			double givenBeforeTime = givenTime - years;
			double givenAfterTime = givenTime + years;
			
			boolean hasBefore = false;
			boolean hasCorupture = false;
			boolean hasAfter = false;
			
			for (SimulatorEvent target : targetMatches) {
				double targetTime = target.getTimeInYears();
				if (targetTime < givenBeforeTime)
					continue;
				if (targetTime > givenAfterTime)
					break;
				if (targetTime < givenTime)
					hasBefore = true;
				else if (targetTime > givenTime)
					hasAfter = true;
				else
					// equal
					hasCorupture = true;
			}
			if (hasBefore)
				numBefore++;
			if (hasCorupture)
				numCorupture++;
			if (hasAfter)
				numAfter++;
			
			if (hasBefore || hasCorupture || hasAfter)
				numTotal++;
		}
		
//		System.out.println("Probability of "+targetName+" given "+givenName
//				+" (within "+(float)years+" years)");
//		System.out.println("\tBefore: "+getProb(numBefore, givenMatches.size()));
//		System.out.println("\tCorupture: "+getProb(numCorupture, givenMatches.size()));
//		System.out.println("\tAfter: "+getProb(numAfter, givenMatches.size()));
//		System.out.println("\tTotal: "
//				+getProb(numTotal, givenMatches.size()));
		probCSV.addLine(targetName,
				getProb(numCorupture, givenMatches.size()),
				getProb(numBefore, givenMatches.size()),
				getProb(numAfter, givenMatches.size()),
				getProb(numTotal, givenMatches.size()),
				getProb(targetMatches.size(), General_EQSIM_Tools.getSimulationDurationYears(events)));
		
		// now do within calcs
		int withinNumCorupture = 0;
		int withinTargetBeforeNum = 0;
		int withinTargetAfterNum = 0;
		int withinTotal = 0;
		for (SimulatorEvent target : targetMatches) {
			double targetTime = target.getTimeInYears();
			double minTime = targetTime - years;
			double maxTime = targetTime + years;
			
			boolean hasBefore = false;
			boolean hasCorupture = false;
			boolean hasAfter = false;
			
			for (SimulatorEvent given : givenMatches) {
				double givenTime = given.getTimeInYears();
				
				if (givenTime > maxTime)
					break;
				if (givenTime < minTime)
					continue;
				
				// we have a match
				if (targetTime < givenTime)
					hasBefore = true;
				else if (targetTime > givenTime)
					hasAfter = true;
				else
					// equal
					hasCorupture = true;
			}
			

			if (hasBefore)
				withinTargetBeforeNum++;
			if (hasCorupture)
				withinNumCorupture++;
			if (hasAfter)
				withinTargetAfterNum++;
			
			if (hasBefore || hasCorupture || hasAfter)
				withinTotal++;
		}
		
//		System.out.println("Percent of "+targetName+" within "+(float)years+" of "+givenName);
//		System.out.println("\tBefore: "
//				+getProb(withinTargetBeforeNum, targetMatches.size()));
//		System.out.println("\tCorupture: "
//				+getProb(withinNumCorupture, targetMatches.size()));
//		System.out.println("\tAfter: "
//				+getProb(withinTargetAfterNum, targetMatches.size()));
//		System.out.println("\tTotal: "
//				+getProb(withinTotal, targetMatches.size()));
		percentCSV.addLine(targetName,
				getProb(withinNumCorupture, targetMatches.size()),
				getProb(withinTargetBeforeNum, targetMatches.size()),
				getProb(withinTargetAfterNum, targetMatches.size()),
				getProb(withinTotal, targetMatches.size()));
	}
	
	private static String getProb(double occurances, double total) {
		double p = occurances / total;
		return (float)(p*100d)+" %";
	}
}
