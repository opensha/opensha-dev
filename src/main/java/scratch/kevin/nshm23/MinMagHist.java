package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.ModSectMinMags;
import org.opensha.sha.earthquake.faultSysSolution.modules.NamedFaults;
import org.opensha.sha.earthquake.faultSysSolution.modules.SectSlipRates;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.JumpProbabilityCalc;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.RuptureProbabilityCalc.BinaryRuptureProbabilityCalc;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_ConstraintBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs.Builder;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs.SubSeisMoRateReduction;
import org.opensha.sha.faultSurface.FaultSection;

public class MinMagHist {

	public static void main(String[] args) throws IOException {
		File dir = new File("C:\\Users\\Kevin Milner\\Downloads");
		FaultSystemSolution sol = FaultSystemSolution.load(new File(dir,
				"results_NSHM23_v2_CoulombRupSet_branch_averaged.zip"));
		FaultSystemRupSet rupSet = sol.getRupSet();
		rupSet.removeModuleInstances(ModSectMinMags.class);
		
		double supraSeisBValue = 0.5;
		
//		BinaryRuptureProbabilityCalc segModel = null;
//		File csvFile = new File(dir, "supra_b_0p5_no_seg.csv");
//		File subsetCSVFile = new File(dir, "supra_b_0p5_no_seg_subset.csv");
		
		BinaryRuptureProbabilityCalc segModel = new BinaryRuptureProbabilityCalc.LogicalAnd(
				new NSHM23_SegmentationModels.NamedFaultSegmentationModel(rupSet.requireModule(NamedFaults.class)),
				new NSHM23_SegmentationModels.ExcludeRupsThroughCreepingSegmentationModel(NSHM23_ConstraintBuilder.findCreepingSection(rupSet)));
		File csvFile = new File(dir, "supra_b_0p5_classic_seg.csv");
		File subsetCSVFile = new File(dir, "supra_b_0p5_classic_seg_subset.csv");
		
		CSVFile<String> csv = new CSVFile<String>(true);
		
		csv.addLine("Subsection Index", "Subsection Name", "Mmin", "Mmax",
				"Creep-Reduced Slip Rate (mm/yr)", "Sub b=1 MoReductFract", "Sub M>=6.5 MoRedFract");
		CSVFile<String> subsetCSV = new CSVFile<String>(true);
		subsetCSV.addLine(csv.getLine(0));
		
		if (segModel != null) {
			List<Integer> retainedSectIDs = new ArrayList<>(rupSet.getNumSections());
			for (int s=0; s<rupSet.getNumSections(); s++)
				retainedSectIDs.add(s);
			rupSet = rupSet.getForSectionSubSet(retainedSectIDs, segModel);
		}
		
		Builder mfdBuilder = new SupraSeisBValInversionTargetMFDs.Builder(rupSet, supraSeisBValue);
		if (segModel != null)
			mfdBuilder.forBinaryRupProbModel(segModel);

		SectSlipRates subB1slips = mfdBuilder.subSeisMoRateReduction(
				SubSeisMoRateReduction.SUB_SEIS_B_1).buildSlipRatesOnly();
		SectSlipRates subM6p5slips = mfdBuilder.subSeisMoRateReduction(
				SubSeisMoRateReduction.SUPRA_B_TO_M6p5).buildSlipRatesOnly();
		
		MinMaxAveTracker mMinTrack = new MinMaxAveTracker();
		MinMaxAveTracker mMaxTrack = new MinMaxAveTracker();
		MinMaxAveTracker slipRateTrack = new MinMaxAveTracker();
		MinMaxAveTracker subB1Track = new MinMaxAveTracker();
		MinMaxAveTracker subM6p5Track = new MinMaxAveTracker();
		
		for (int s=0; s<rupSet.getNumSections(); s++) {
			List<String> line = new ArrayList<>(csv.getNumCols());
			line.add(s+"");
			FaultSection sect = rupSet.getFaultSectionData(s);
			line.add(sect.getSectionName());
			double mmin = rupSet.getMinMagForSection(s);
			mMinTrack.addValue(mmin);
			line.add((float)mmin+"");
			double mmax = rupSet.getMaxMagForSection(s);
			mMaxTrack.addValue(mmax);
			line.add((float)mmax+"");
			double slip = sect.getReducedAveSlipRate();
			slipRateTrack.addValue(slip);
			line.add((float)slip+"");
			double subB1target = subB1slips.getSlipRate(s)*1e3;
			double subB1red = slip == 0d ? 0d : (slip-subB1target)/slip;
			subB1Track.addValue(subB1red);
			line.add((float)subB1red+"");
			double subM6p5target = subM6p5slips.getSlipRate(s)*1e3;
			double subM6p5red;
			if ((float)subM6p5target == (float)slip || slip == 0d || Math.abs(subM6p5target-slip) < 1e-16)
				subM6p5red = 0;
			else
				subM6p5red = (slip-subM6p5target)/slip;
			subM6p5Track.addValue(subM6p5red);
			line.add((float)subM6p5red+"");
			csv.addLine(line);
			if (subM6p5red > 0)
				subsetCSV.addLine(line);
		}
		// now in min/max/avg line
		List<String> minLine = new ArrayList<>();
		List<String> maxLine = new ArrayList<>();
		List<String> avgLine = new ArrayList<>();
		minLine.add("");
		maxLine.add("");
		avgLine.add("");
		minLine.add("MINIMUM");
		maxLine.add("MAXIMUM");
		avgLine.add("AVERAGE");
		minLine.add((float)mMinTrack.getMin()+"");
		maxLine.add((float)mMinTrack.getMax()+"");
		avgLine.add((float)mMinTrack.getAverage()+"");
		minLine.add((float)mMaxTrack.getMin()+"");
		maxLine.add((float)mMaxTrack.getMax()+"");
		avgLine.add((float)mMaxTrack.getAverage()+"");
		minLine.add((float)slipRateTrack.getMin()+"");
		maxLine.add((float)slipRateTrack.getMax()+"");
		avgLine.add((float)slipRateTrack.getAverage()+"");
		minLine.add((float)subB1Track.getMin()+"");
		maxLine.add((float)subB1Track.getMax()+"");
		avgLine.add((float)subB1Track.getAverage()+"");
		minLine.add((float)subM6p5Track.getMin()+"");
		maxLine.add((float)subM6p5Track.getMax()+"");
		avgLine.add((float)subM6p5Track.getAverage()+"");
		
		int[] numNonZeroes = new int[csv.getNumCols()];
		for (int row=1; row<csv.getNumRows(); row++) {
			for (int col=2; col<numNonZeroes.length; col++)
				if (csv.getDouble(row, col) > 0d)
					numNonZeroes[col]++;
		}
		List<String> numNonZeroLine = new ArrayList<>();
		List<String> fractNonZeroLine = new ArrayList<>();
		numNonZeroLine.add("");
		fractNonZeroLine.add("");
		numNonZeroLine.add("Num > 0");
		fractNonZeroLine.add("Fract > 0");
		for (int col=2; col<numNonZeroes.length; col++) {
			numNonZeroLine.add(numNonZeroes[col]+"");
			double fract = (double)numNonZeroes[col]/(double)rupSet.getNumSections();
			fractNonZeroLine.add((float)fract+"");
		}

		csv.addLine(minLine);
		csv.addLine(maxLine);
		csv.addLine(avgLine);
		csv.addLine(numNonZeroLine);
		csv.addLine(fractNonZeroLine);
		subsetCSV.addLine(minLine);
		subsetCSV.addLine(maxLine);
		subsetCSV.addLine(avgLine);
		subsetCSV.addLine(numNonZeroLine);
		subsetCSV.addLine(fractNonZeroLine);
		
		double[] minMags = new double[rupSet.getNumSections()];
		for (int s=0; s<minMags.length; s++)
			minMags[s] = rupSet.getMinMagForSection(s);
		
		double maxMin = StatUtils.max(minMags);
		double minMin = StatUtils.min(minMags);
		HistogramFunction hist = HistogramFunction.getEncompassingHistogram(minMin, maxMin, 0.1);
		for (double minMag : minMags)
			hist.add(hist.getClosestXIndex(minMag), 1d);
		
		for (int s=0; s<minMags.length; s++)
			if (minMags[s] == minMin)
				System.out.println("Smallest: M"+(float)minMin+" for "+rupSet.getFaultSectionData(s).getSectionName());
		for (int s=0; s<minMags.length; s++)
			if (minMags[s] == maxMin)
				System.out.println("Largest: M"+(float)maxMin+" for "+rupSet.getFaultSectionData(s).getSectionName());
		System.out.println(hist);
		
		csv.writeToFile(csvFile);
		subsetCSV.writeToFile(subsetCSVFile);
	}

}
