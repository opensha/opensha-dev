package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.dom4j.Document;
import org.dom4j.DocumentException;
import org.dom4j.Element;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.util.XMLUtils;
import org.opensha.refFaultParamDb.vo.DeformationModel;
import org.opensha.sha.earthquake.param.BPTAveragingTypeParam;
import org.opensha.sha.earthquake.param.HistoricOpenIntervalParam;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityOptions;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.FaultSystemSolutionFetcher;
import scratch.UCERF3.analysis.CompoundFSSPlots.ERFBasedRegionalMagProbPlot;
import scratch.UCERF3.analysis.MPJ_ERF_ProbGainCalc;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.mean.MeanUCERF3;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.UCERF3_DataUtils;

public class ERFFaultProbDebug {
	
	public static void main(String[] args) throws IOException, DocumentException {
		Document doc = XMLUtils.createDocumentWithRoot();
		Element root = doc.getRootElement();
		root.addAttribute("test1", "a\nb\nc");
		String s = System.getProperty("line.separator");
		root.addAttribute("test2", "a"+s+"b"+s+"c");
		s = String.format("%n");
		root.addAttribute("test3", "a"+s+"b"+s+"c");
		s = "&#10;";
		root.addAttribute("test4", "a"+s+"b"+s+"c");
		s = "<br>";
		root.addAttribute("test5", "a"+s+"b"+s+"c");
		XMLUtils.writeDocumentToFile(new File("/tmp/test.xml"), doc);
		System.exit(0);
		// mean ERF
//		InversionFaultSystemSolution meanSol = FaultSystemIO.loadInvSol(
//				new File(new File(UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, "InversionSolutions"),
//						"2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		
		CompoundFaultSystemSolution cfss = CompoundFaultSystemSolution.fromZipFile(
				new File(new File(UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, "InversionSolutions"),
						"2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip"));
//		FaultSystemSolution sol = cfss.getSolution(LogicTreeBranch.fromFileName(
//				"FM3_1_ZENGBB_ShConStrDrp_DsrTap_CharConst_M5Rate9.6_MMaxOff7.6_NoFix_SpatSeisU2"));
		FaultSystemSolution sol = cfss.getSolution(LogicTreeBranch.fromFileName(
				"FM3_1_NEOK_EllBsqrtLen_DsrUni_CharConst_M5Rate6.5_MMaxOff7.3_NoFix_SpatSeisU2"));
		
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		
//		MeanUCERF3 erf = new MeanUCERF3();
//		erf.setMeanParams(0d, false, 0d, MeanUCERF3.RAKE_BASIS_NONE);
//		FaultSystemRupSet rupSet = erf.getSolution().getRupSet();
		
		Map<MagDependentAperiodicityOptions, Double> covMap = Maps.newHashMap();
		covMap.put(MagDependentAperiodicityOptions.LOW_VALUES, FaultSystemSolutionERF.PREF_BLEND_COV_LOW_WEIGHT);
		covMap.put(MagDependentAperiodicityOptions.MID_VALUES, FaultSystemSolutionERF.PREF_BLEND_COV_MID_WEIGHT);
		covMap.put(MagDependentAperiodicityOptions.HIGH_VALUES, FaultSystemSolutionERF.PREF_BLEND_COV_HIGH_WEIGHT);
		covMap.put(null, FaultSystemSolutionERF.PREF_BLEND_POISSON_WEIGHT);
		
		Map<MagDependentAperiodicityOptions, Double> resultsMap = Maps.newHashMap();
		double blendVal = 0;
		
		for (MagDependentAperiodicityOptions cov : covMap.keySet()) {
			boolean[] meanIndeads;
			if (cov == null) {
				erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
				meanIndeads = new boolean [] {false, true}; // also do pref blend
			} else {
				erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_BPT);
				erf.setParameter(MagDependentAperiodicityParam.NAME, cov);
				erf.setParameter(HistoricOpenIntervalParam.NAME,
						(double)(FaultSystemSolutionERF.START_TIME_DEFAULT - 1875));
				meanIndeads = new boolean [] {false};
			}
			for (boolean meanInstead : meanIndeads) {
				if (meanInstead) {
					erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_PREF_BLEND);
					erf.setParameter(HistoricOpenIntervalParam.NAME,
							(double)(FaultSystemSolutionERF.START_TIME_DEFAULT - 1875));
				}
				erf.getTimeSpan().setDuration(30d);
				erf.updateForecast();
				
//				for (Parameter<?> param : erf.getAdjustableParameterList()) {
//					System.out.println("ERF Param: "+param.getName()+" value: "+param.getValue());
//				}
//				System.out.println("ERF Sol: "+((InversionFaultSystemSolution)erf.getSolution()).getLogicTreeBranch());
				
//				System.out.println("Ave type: "+erf.getParameter(BPTAveragingTypeParam.NAME).getValue());
				
				double minMag = 6.7d;
				
				String debugFault = "Calaveras";
				
				Map<String, List<Integer>> mainFaultsMap = FaultModels.parseNamedFaultsAltFile(
						UCERF3_DataUtils.getReader("FaultModels",
						"MainFaultsForTimeDepComparison.txt"));
				ArrayList<String> mainFaultsSorted = Lists.newArrayList(mainFaultsMap.keySet());
				Collections.sort(mainFaultsSorted);
				Map<String, Collection<Integer>> mappings = Maps.newHashMap();
				
				// MPJ calc style
				for (String fault : mainFaultsSorted) {
					HashSet<Integer> rups = new HashSet<Integer>();
					for (Integer parentID : mainFaultsMap.get(fault)) {
						List<Integer> parentRups = rupSet.getRupturesForParentSection(parentID);
						if (parentRups != null)
							rups.addAll(parentRups);
					}
					mappings.put(fault, rups);
				}
				
				double mpjProb = MPJ_ERF_ProbGainCalc.calcFaultProb(erf, rupSet, mappings.get(debugFault), minMag);
				
				EvenlyDiscretizedFunc probFunc = new EvenlyDiscretizedFunc(6.7d, 1, 0.1d);
				EvenlyDiscretizedFunc rateFunc = new EvenlyDiscretizedFunc(6.7d, 1, 0.1d);
				ERFBasedRegionalMagProbPlot.calcFaultProbs(probFunc, rateFunc, erf, rupSet, mappings.get(debugFault));
				double cfssProb = probFunc.getY(0);
				
				if (meanInstead)
					System.out.println("PREF_BLEND");
				else
					System.out.println("COV: "+cov);
				System.out.println("MPJ calc prob: "+mpjProb);
				System.out.println("CFSS calc prob: "+cfssProb);
				
				if (meanInstead)
					blendVal = cfssProb;
				else
					resultsMap.put(cov, cfssProb);
				
			}
		}
		
		double calcBlend = 0;
		for (MagDependentAperiodicityOptions cov : covMap.keySet()) {
			double weight = covMap.get(cov);
			calcBlend += weight*resultsMap.get(cov);
		}
		System.out.println("Blend:\t"+blendVal);
		System.out.println("Calc Blend:\t"+calcBlend);
	}

}
