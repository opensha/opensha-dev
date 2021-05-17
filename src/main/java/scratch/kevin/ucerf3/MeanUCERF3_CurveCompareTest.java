package scratch.kevin.ucerf3;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.dom4j.DocumentException;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.util.DataUtils;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.param.ApplyGardnerKnopoffAftershockFilterParam;
import org.opensha.sha.earthquake.param.BPTAveragingTypeOptions;
import org.opensha.sha.earthquake.param.BPTAveragingTypeParam;
import org.opensha.sha.earthquake.param.HistoricOpenIntervalParam;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;

import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.mean.MeanUCERF3;
import scratch.UCERF3.erf.mean.RuptureCombiner;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.LastEventData;
import scratch.UCERF3.utils.UCERF3_DataUtils;

public class MeanUCERF3_CurveCompareTest {
	
	public static void main(String[] args) throws IOException, DocumentException {
		File dir = new File("/home/kevin/OpenSHA/UCERF3/MeanUCERF3-timedep-tests");
		
		FaultModels fm = FaultModels.FM3_1;
		double udTol = 100d;
		boolean branchAveragedSol = false;
		boolean meanAsAvg = false;
		
		boolean clearCache = true;
		
		FaultSystemSolution baSol = FaultSystemIO.loadSol(
				new File(new File(UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, "InversionSolutions"),
						"2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		
		FaultSystemSolutionERF erf;
		if (branchAveragedSol) {
			udTol = 0d;
			fm = FaultModels.FM3_1;
			baSol = new FaultSystemSolution(baSol.getRupSet(), baSol.getRateForAllRups());
			baSol = getRandomOrder(baSol);
			erf = new FaultSystemSolutionERF(baSol);
		} else {
			erf = new MeanUCERF3();
			((MeanUCERF3)erf).setMeanParams(0d, false, 0d, MeanUCERF3.RAKE_BASIS_MEAN);
//			((MeanUCERF3)erf).setMeanParams(udTol, false, 1d, MeanUCERF3.RAKE_BASIS_MEAN);
			if (fm != null)
				erf.getParameter(MeanUCERF3.FAULT_MODEL_PARAM_NAME).setValue(fm.name());
			if (meanAsAvg) {
				erf.updateForecast();
//				compareSolutions(baSol, erf.getSolution());
				FaultSystemSolution sol = erf.getSolution();
//				sol = copyMagsFrom(sol, baSol);
				sol = getRandomOrder(sol);
				getBASubsetSol(sol, baSol);
				// this uses the BA sol but only rups that are in the mean sol
				Map<Integer, List<LastEventData>> lastEventData = LastEventData.load();
				LastEventData.populateSubSects(sol.getRupSet().getFaultSectionDataList(), lastEventData);
				erf = new FaultSystemSolutionERF(sol);
			}
		}
		erf.getParameter(IncludeBackgroundParam.NAME).setValue(IncludeBackgroundOption.EXCLUDE);
		
		ScalarIMR imr = AttenRelRef.CB_2008.instance(null);
		imr.setParamDefaults();
		imr.setIntensityMeasure(PGA_Param.NAME);
		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(PGA_Param.NAME);
		Site site = new Site(new Location(34.055, -118.2467)); // downtown LA
		site.addParameterList(imr.getSiteParams());
		
		HazardCurveCalculator calc = new HazardCurveCalculator();
		
		String fmAdd;
		if (fm == null)
			fmAdd = "";
		else
			fmAdd = fm.name()+"_";
		
		String udAdd;
		if (udTol == 0)
			udAdd = "";
		else
			udAdd = "ud"+(float)udTol+"_";
		
		String baAdd;
		if (branchAveragedSol)
			baAdd = "baSol_";
		else
			baAdd = "";
		
		File meanPoissonFile = new File(dir, fmAdd+baAdd+udAdd+"meanucerf3_poisson.txt");
		File meanTimeDepFile = new File(dir, fmAdd+baAdd+udAdd+"meanucerf3_timedep.txt");
		
		DiscretizedFunc meanTimeDepCurve, meanPoissonCurve;
		
		erf.setParameter(ApplyGardnerKnopoffAftershockFilterParam.NAME, false);
		if (meanTimeDepFile.exists() && !clearCache) {
			meanTimeDepCurve = ArbitrarilyDiscretizedFunc.loadFuncFromSimpleFile(meanTimeDepFile.getAbsolutePath());
		} else {
			erf.getParameter(ProbabilityModelParam.NAME).setValue(
					ProbabilityModelOptions.U3_PREF_BLEND);
			erf.setParameter(HistoricOpenIntervalParam.NAME,
					(double)(FaultSystemSolutionERF.START_TIME_DEFAULT-1875));
			erf.setParameter(BPTAveragingTypeParam.NAME,
					BPTAveragingTypeOptions.AVE_RI_AVE_NORM_TIME_SINCE);
			
			erf.getTimeSpan().setDuration(30d);
			erf.updateForecast();
			FaultSystemSolution sol = erf.getSolution();
			FaultSystemRupSet rupSet = sol.getRupSet();
			System.out.println("Num Sol Ruptures: "+rupSet.getNumRuptures());
			HashSet<String> uniqueRups = new HashSet<String>();
			int numNonzero = 0;
			List<List<Integer>> sectionIndicesForAllRups = rupSet.getSectionIndicesForAllRups();
			for (int i = 0; i < sectionIndicesForAllRups.size(); i++) {
				List<Integer> rup = sectionIndicesForAllRups.get(i);
				rup = Lists.newArrayList(rup);
				Collections.sort(rup);
				String str = Joiner.on(",").join(rup);
				uniqueRups.add(str);
				if (sol.getRateForRup(i) > 0)
					numNonzero++;
			}
			System.out.println("Num Unique Rups: "+uniqueRups.size());
			System.out.println("Num Nonzero Rups: "+numNonzero);
			
			System.out.println("Calculating Time Dep Curve with "+erf.getNumSources()+" sources");
			int totNumERFRups = 0;
			for (int sourceID=0; sourceID<erf.getNumSources(); sourceID++)
				totNumERFRups += erf.getNumRuptures(sourceID);
			System.out.println("Calculating Time Dep Curve with "+erf.getNumSources()+" sources ("+totNumERFRups+" tot rups)");
			
			meanTimeDepCurve = xVals.deepClone();
			calc.getHazardCurve(meanTimeDepCurve, site, imr, erf);
			
			ArbitrarilyDiscretizedFunc.writeSimpleFuncFile(meanTimeDepCurve, meanTimeDepFile);
		}
		
		if (meanPoissonFile.exists() && !clearCache) {
			meanPoissonCurve = ArbitrarilyDiscretizedFunc.loadFuncFromSimpleFile(meanPoissonFile.getAbsolutePath());
		} else {
			erf.getParameter(ProbabilityModelParam.NAME).setValue(
					ProbabilityModelOptions.POISSON);
			
			erf.getTimeSpan().setDuration(30d);
			erf.updateForecast();
			System.out.println("Num Sources: "+erf.getNumSources());
			
			System.out.println("Calculating Time Indep Curve with "+erf.getNumSources()+" sources");
			
			meanPoissonCurve = xVals.deepClone();
			calc.getHazardCurve(meanPoissonCurve, site, imr, erf);
			
			ArbitrarilyDiscretizedFunc.writeSimpleFuncFile(meanPoissonCurve, meanPoissonFile);
		}
		
		DiscretizedFunc averagedTimeDepCurve = calcMean(new File(dir, "curves_timedep"), fm);
		DiscretizedFunc averagedPoissonCurve = calcMean(new File(dir, "curves_poisson"), fm);
		
		System.out.println("Poisson Diagnostics:");
		writeDiagnostics(averagedPoissonCurve, meanPoissonCurve);
		
		System.out.println("\n\nTime Dep Diagnostics:");
		writeDiagnostics(averagedTimeDepCurve, meanTimeDepCurve);
	}
	
	private static DiscretizedFunc calcMean(File curvesDir, FaultModels fm) throws FileNotFoundException, IOException {
		DiscretizedFunc curve = null;
		double totWeight = 0;
		int numAvg = 0;
		for (File file : curvesDir.listFiles()) {
			String name = file.getName();
			if (!name.endsWith(".txt"))
				continue;
			LogicTreeBranch branch = LogicTreeBranch.fromFileName(name);
			FaultModels fileFM = branch.getValue(FaultModels.class);
			if (fm != null && fileFM != fm)
				continue;
			DiscretizedFunc func = ArbitrarilyDiscretizedFunc.loadFuncFromSimpleFile(file.getAbsolutePath());
			if (curve == null) {
				curve = new ArbitrarilyDiscretizedFunc();
				for (int i=0; i<func.size(); i++)
					curve.set(func.getX(i), 0d);
			}
			
			double fmWeight = 1;
			if (fm == null)
				fmWeight = LogicTreeBranch.getNormalizedWt(fileFM,
					InversionModels.CHAR_CONSTRAINED);
			double dmWeight = LogicTreeBranch.getNormalizedWt(branch.getValue(DeformationModels.class),
					InversionModels.CHAR_CONSTRAINED);
			double scaleWeight = LogicTreeBranch.getNormalizedWt(branch.getValue(ScalingRelationships.class),
					InversionModels.CHAR_CONSTRAINED);
			
			double weight = fmWeight * dmWeight * scaleWeight;
			totWeight += weight;
			
			numAvg++;
			
			for (int i=0; i<curve.size(); i++)
				curve.set(i, curve.getY(i) + weight*func.getY(i));
		}
		
		System.out.println("Tot weight for "+curvesDir.getName()+": "+totWeight+" ("+numAvg+" curves)");
		
		return curve;
	}
	
	private static void writeDiagnostics(DiscretizedFunc avgCurve, DiscretizedFunc meanU3Curve) {
		double maxDiscrep = 0;
		double maxDiscrepPercent = 0;
		for (int i=0; i<avgCurve.size(); i++) {
			double x = avgCurve.getX(i);
			double avgY = avgCurve.getY(i);
			double meanU3Y = meanU3Curve.getY(i);
			double discrep = Math.abs(avgY - meanU3Y);
			double pDiff = DataUtils.getPercentDiff(meanU3Y, avgY);
			
			if (discrep > maxDiscrep)
				maxDiscrep = discrep;
			if (pDiff > maxDiscrepPercent)
				maxDiscrepPercent = pDiff;
			
			System.out.println("X: "+(float)x+"\tavg: "+(float)avgY+"\tmeanU3: "+(float)meanU3Y+"\tdiff: "
					+(float)discrep+"\tpDiff: "+(float)pDiff+" %");
		}
		
		System.out.println("Max discrep: "+(float)+maxDiscrep);
		System.out.println("Max pDiff: "+(float)+maxDiscrepPercent+" %");
	}
	
	public static FaultSystemSolution getRandomOrder(final FaultSystemSolution other) {
		final FaultSystemRupSet origRupSet = other.getRupSet();
		List<Integer> rupIDs = Lists.newArrayList();
		for (int i=0; i<origRupSet.getNumRuptures(); i++)
			rupIDs.add(i);
		Collections.shuffle(rupIDs);
		
		List<List<Integer>> sectionForRups = Lists.newArrayList();
		for (int r : rupIDs)
			sectionForRups.add(origRupSet.getSectionsIndicesForRup(r));
		
		FaultSystemRupSet newRupSet = new FaultSystemRupSet(origRupSet.getFaultSectionDataList(),
				origRupSet.getSlipRateForAllSections(), origRupSet.getSlipRateStdDevForAllSections(),
				origRupSet.getAreaForAllSections(), sectionForRups,
				getResortedArray(origRupSet.getMagForAllRups(), rupIDs),
				getResortedArray(origRupSet.getAveRakeForAllRups(), rupIDs),
				getResortedArray(origRupSet.getAreaForAllRups(), rupIDs),
				getResortedArray(origRupSet.getLengthForAllRups(), rupIDs), origRupSet.getInfoString());
		
		return new FaultSystemSolution(newRupSet, getResortedArray(other.getRateForAllRups(), rupIDs));
	}
	
	private static double[] getResortedArray(double[] orig, List<Integer> order) {
		if (orig == null)
			return null;
		double[] ret = new double[orig.length];
		for (int i=0; i<orig.length; i++)
			ret[i] = orig[order.get(i)];
		return ret;
	}
	
	private static class Rupture implements Comparable<Rupture> {
		private HashSet<Integer> sects;
		private double mag;
		private double rate;
		private double area;
		private double rake;
		
		public Rupture(HashSet<Integer> sects, double mag, double rate,
				double area, double rake) {
			super();
			this.sects = sects;
			this.mag = mag;
			this.rate = rate;
			this.area = area;
			this.rake = rake;
		}
		
		@Override
		public int compareTo(Rupture o) {
			// sort by rate, decreasing
			return -Double.compare(rate, o.rate);
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + ((sects == null) ? 0 : sects.hashCode());
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			Rupture other = (Rupture) obj;
			if (sects == null) {
				if (other.sects != null)
					return false;
			} else if (!sects.equals(other.sects))
				return false;
			return true;
		}
	}
	
	private static void compareSolutions(FaultSystemSolution avgSol, FaultSystemSolution meanSol) {
		Map<String, Integer> avgMap = getRupsMap(avgSol);
		Map<String, Integer> meanMap = getRupsMap(meanSol);
		
		for (String rup : avgMap.keySet()) {
			int avgIndex = avgMap.get(rup);
			printStats(avgSol, avgIndex);
			Integer meanIndex = meanMap.get(rup);
			if (meanIndex == null)
				System.out.println("\t -- NONE --");
			else
				printStats(meanSol, meanIndex);
			System.out.println();
			try {
				Thread.sleep(1000);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	private static void printStats(FaultSystemSolution sol, int index) {
		FaultSystemRupSet rupSet = sol.getRupSet();
		System.out.println("\t"+index+".\t"+(float)rupSet.getMagForRup(index)+"\t"+(float)rupSet.getAveRakeForRup(index)
				+"\t"+(float)rupSet.getAreaForRup(index)+"\t"+(float)sol.getRateForRup(index));
	}
	
	private static Map<String, Integer> getRupsMap(FaultSystemSolution sol) {
		Map<String, Integer> map = Maps.newHashMap();
		
		FaultSystemRupSet rupSet = sol.getRupSet();
		for (int i=0; i<rupSet.getNumRuptures(); i++) {
			List<String> subSects = Lists.newArrayList();
			for (int sect : rupSet.getSectionsIndicesForRup(i)) {
				String name = rupSet.getFaultSectionData(sect).getName();
				if (name.contains(" (instance"))
					name = name.substring(0, name.indexOf(" (instance"));
				subSects.add(name);
			}
			Collections.sort(subSects);
			map.put(Joiner.on(",").join(subSects), i);
		}
		
		return map;
	}
	
	private static FaultSystemSolution copyMagsFrom(FaultSystemSolution meanSol, FaultSystemSolution baSol) {
		Map<String, Integer> avgMap = getRupsMap(baSol);
		Map<String, Integer> meanMap = getRupsMap(meanSol);
		
		FaultSystemRupSet meanRupSet = meanSol.getRupSet();
		
		double[] baMags = baSol.getRupSet().getMagForAllRups();
		double[] meanNewMags = new double[meanRupSet.getNumRuptures()];
		
		for (String rup : meanMap.keySet()) {
			int meanIndex = meanMap.get(rup);
			int avgIndex = avgMap.get(rup);
			
			meanNewMags[meanIndex] = baMags[avgIndex];
		}
		
		return new FaultSystemSolution(new FaultSystemRupSet(
				meanRupSet.getFaultSectionDataList(), meanRupSet.getSlipRateForAllSections(),
				meanRupSet.getSlipRateStdDevForAllSections(), meanRupSet.getAreaForAllSections(),
				meanRupSet.getSectionIndicesForAllRups(), meanNewMags, meanRupSet.getAveRakeForAllRups(),
				meanRupSet.getAreaForAllRups(), meanRupSet.getLengthForAllRups(), meanSol.getInfoString()),
				meanSol.getRateForAllRups());
	}
	
	private static FaultSystemSolution getBASubsetSol(FaultSystemSolution meanSol, FaultSystemSolution baSol) {
		Map<String, Integer> avgMap = getRupsMap(baSol);
		Map<String, Integer> meanMap = getRupsMap(meanSol);
		
		List<Integer> indexesToKeep = Lists.newArrayList();
		
		for (String rup : avgMap.keySet()) {
			int avgIndex = avgMap.get(rup);
			Integer meanIndex = meanMap.get(rup);
			
			if (meanIndex != null)
				indexesToKeep.add(avgIndex);
		}
		
		FaultSystemRupSet rupSet = new RuptureCombiner.SubsetRupSet(baSol.getRupSet(), indexesToKeep);
		double[] newRates = new double[indexesToKeep.size()];
		for (int i=0; i<newRates.length; i++)
			newRates[i] = baSol.getRateForRup(indexesToKeep.get(i));
		
		return new FaultSystemSolution(rupSet, newRates);
	}

}
