package scratch.pagem;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;

import org.dom4j.DocumentException;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.param.HistoricOpenIntervalParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.analysis.FaultSysSolutionERF_Calc;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.enumTreeBranches.MaxMagOffFault;
import scratch.UCERF3.enumTreeBranches.MomentRateFixes;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.enumTreeBranches.SlipAlongRuptureModels;
import scratch.UCERF3.enumTreeBranches.SpatialSeisPDF;
import scratch.UCERF3.enumTreeBranches.TotalMag5Rate;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.inversion.laughTest.LaughTestFilter;
import scratch.UCERF3.simulatedAnnealing.SerialSimulatedAnnealing;
import scratch.UCERF3.simulatedAnnealing.SimulatedAnnealing;
import scratch.UCERF3.simulatedAnnealing.ThreadedSimulatedAnnealing;
import scratch.UCERF3.simulatedAnnealing.completion.CompletionCriteria;
import scratch.UCERF3.simulatedAnnealing.completion.IterationCompletionCriteria;
import scratch.UCERF3.simulatedAnnealing.completion.ProgressTrackingCompletionCriteria;
import scratch.UCERF3.simulatedAnnealing.completion.TimeCompletionCriteria;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.UCERF3_DataUtils;
import scratch.UCERF3.utils.aveSlip.AveSlipConstraint;
import scratch.UCERF3.utils.paleoRateConstraints.PaleoFitPlotter;
import scratch.UCERF3.utils.paleoRateConstraints.PaleoProbabilityModel;
import scratch.UCERF3.utils.paleoRateConstraints.PaleoRateConstraint;
import scratch.UCERF3.utils.paleoRateConstraints.UCERF2_PaleoRateConstraintFetcher;
import scratch.UCERF3.utils.paleoRateConstraints.UCERF3_PaleoRateConstraintFetcher;
import cern.colt.matrix.tdouble.DoubleMatrix2D;


/**
 * This class generates UCERF3-TD participation probs for multiple parent section IDs, without double counting.
 * See Skype message with Kevin Oct 14, 2016
 *
 */

public class U3TD_Participation_Probs {


	public static void main(String[] args) throws IOException, DocumentException {
		
		FaultSystemSolution sol;
		   
		// Make sure this is the updated version that includes the date of last event data in it
		// From: http://opensha.usc.edu/ftp/kmilner/ucerf3/2013_05_10-ucerf3p3-production-10runs/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip
		// Last updated: 10/15/2015
		sol = FaultSystemIO.loadSol(
//			           new File("/Users/pagem/Desktop/"
//			           + "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
						new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions"
								+ "/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		List<FaultSectionPrefData> faultSectionData = sol.getRupSet().getFaultSectionDataList();
		
		// Print list of all parent section numbers and names
		// ListIterator<FaultSectionPrefData> faultSectionDataIterator = faultSectionData.listIterator();
		// int prevIndex=0;
		// while (faultSectionDataIterator.hasNext()) {			
		//	int nextIndex = faultSectionDataIterator.nextIndex();
		//	if (faultSectionData.get(nextIndex).getParentSectionId() != prevIndex) {
		//		System.out.println(faultSectionData.get(nextIndex).getParentSectionId()+"\t"+faultSectionData.get(nextIndex).getParentSectionName());
		//	}
		//	faultSectionDataIterator.next();	
		//	prevIndex = faultSectionData.get(nextIndex).getParentSectionId();
		// }

		
		erf.getParameter(ProbabilityModelParam.NAME).setValue(ProbabilityModelOptions.U3_PREF_BLEND);
		erf.getTimeSpan().setStartTime(2017); // start year
		erf.setParameter(HistoricOpenIntervalParam.NAME, // historical open interval
		erf.getTimeSpan().getStartTimeYear()-1875d);
//		erf.getParameter(ProbabilityModelParam.NAME).setValue(ProbabilityModelOptions.POISSON);
		
		erf.getTimeSpan().setDuration(30d);
		
		erf.updateForecast();
		
		// Add up participation probs for multiple parent sections
		// Arguments: erf, minimum magnitude, parent section IDs (one or more)
		double minMag = 6;
//		int[] parentIDs = {32, 285, 300, 287};
		int[] parentIDs = {170};
		double prob = FaultSysSolutionERF_Calc.calcParticipationProbForParentSects(erf, minMag, parentIDs);
		
		// Print parent section IDs, names, and total participation probability
		for (int i=0; i<parentIDs.length; i++) {
			int parID=parentIDs[i];
			ListIterator<FaultSectionPrefData> faultSectionDataIterator = faultSectionData.listIterator();
			boolean parentSectFound=false;
			while (!parentSectFound && faultSectionDataIterator.hasNext()) {				
				int nextIndex = faultSectionDataIterator.nextIndex();
				if (faultSectionData.get(nextIndex).getParentSectionId()==parID) {
					parentSectFound=true;
					System.out.println(faultSectionData.get(nextIndex).getParentSectionId()+"\t"+faultSectionData.get(nextIndex).getParentSectionName());
				}
				faultSectionDataIterator.next();	
			}
			
		}
		System.out.println("\nTime span = "+erf.getTimeSpan().getDuration()+" "+erf.getTimeSpan().getDurationUnits());
		System.out.println("Minimum magnitude = " + minMag);
		System.out.println("Participation probability = " + prob);
		
			
	}


}

