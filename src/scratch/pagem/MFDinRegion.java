package scratch.pagem;

import java.io.File;
import java.io.IOException;

import org.dom4j.DocumentException;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.earthquake.calc.ERF_Calculator;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.utils.FaultSystemIO;


// Code from Kevin to find incremental or participation MFDs within a rectangular or circular region
// for UCERF3 T-I or T-D
// See e-mail dated March 21, 2016
public class MFDinRegion{
	
	public static void main(String[] args) throws IOException, DocumentException {

	   FaultSystemSolution sol;
	   
	   sol = FaultSystemIO.loadSol(
		           new File("/Users/pagem/Desktop/"
		           + "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
	   FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
	   
	   // Poisson
	   erf.getParameter(ProbabilityModelParam.NAME).setValue(ProbabilityModelOptions.POISSON);
	   // UCERF3 TD
//	   erf.getParameter(ProbabilityModelParam.NAME).setValue(ProbabilityModelOptions.U3_PREF_BLEND);
//	   erf.getTimeSpan().setStartTime(2014); // start year
//	   erf.setParameter(HistoricOpenIntervalParam.NAME, // historical open interval
//	   erf.getTimeSpan().getStartTimeYear()-1875d);

	   // duration - only really necessary for TD calculations, as MFD later is annualized
	   erf.getTimeSpan().setDuration(30d);
	   erf.updateForecast();

	   // circular region with given center and radius
	   // 1933 Long Beach eq location: (33.665N, 117.975W)
	   Region reg = new Region(new Location(33.665, -117.975), 5d);
	   // rectangular region with these corners
//	   Region reg = new Region(new Location(34, -118), new Location(35, -120));

	   double minMag = 5d;
	   int numMag = 51;
	   double deltaMag = 0.1;

	   // Participation MFD
	   IncrementalMagFreqDist mfd = ERF_Calculator.getParticipationMagFreqDistInRegion(erf, reg, minMag, numMag, deltaMag, true);
	   // Nucleation MFD - this will be slower
//	   IncrementalMagFreqDist mfd = ERF_Calculator.getMagFreqDistInRegion(erf, reg, minMag, numMag, deltaMag, true);

	   System.out.println(mfd);
       
	}
}