package scratch.kevin.nshm23;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertainIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertaintyBoundType;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportMetadata;
import org.opensha.sha.earthquake.faultSysSolution.reports.RupSetMetadata;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.SectBySectDetailPlots;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

public class SectMFDCompare {

	public static void main(String[] args) throws IOException {
		File solFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_08_12-nshm23_branches-NSHM23_v2-CoulombRupSet-NSHM23_Avg-TotNuclRate-SubB1-ThreshAvgIterRelGR/"
				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged.zip");
		
//		int sectIndex = 1474; // Gallina
		int sectIndex = 443; // Broadmouth Canyon - James Peak, Subsection 1
		
		boolean sectParentPage = true;
		
		FaultSystemSolution sol = FaultSystemSolution.load(solFile);
		
		File outputDir = solFile.getParentFile();
		mfdComparison(sol, sectIndex, outputDir);
		
		if (sectParentPage) {
			ReportMetadata meta = new ReportMetadata(new RupSetMetadata("Solution", sol));
			SectBySectDetailPlots plot = new SectBySectDetailPlots();
			plot.plotSingleParent(outputDir, meta, sol.getRupSet().getFaultSectionData(sectIndex).getParentSectionId());
		}
	}
	
	private static Range MFD_Y_RANGE = new Range(1e-12, 1e-4);
	
	private static void mfdComparison(FaultSystemSolution sol, int sectIndex, File outputDir) throws IOException {
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		InversionTargetMFDs targets = rupSet.requireModule(InversionTargetMFDs.class);
		List<? extends IncrementalMagFreqDist> sectTargets = targets.getOnFaultSupraSeisNucleationMFDs();
		
		IncrementalMagFreqDist sectTarget = sectTargets.get(sectIndex);
		IncrementalMagFreqDist solNucl = sol.calcNucleationMFD_forSect(sectIndex,
				sectTarget.getMinX(), sectTarget.getMaxX(), sectTarget.size());
		UncertainIncrMagFreqDist uncertTarget = null;
		if (sectTarget instanceof UncertainIncrMagFreqDist)
			uncertTarget = (UncertainIncrMagFreqDist)sectTarget;
		
		int minNonZeroIndex = sectTarget.size();
		int maxNonZeroIndex = 0;
		for (int i=0; i<sectTarget.size(); i++) {
			if (sectTarget.getY(i) > 0 || solNucl.getY(i) > 0) {
				minNonZeroIndex = Integer.min(i, minNonZeroIndex);
				maxNonZeroIndex = Integer.max(i, maxNonZeroIndex);
			}
		}
		
		if (uncertTarget == null)
			System.out.println("Mag\tTarget\tSolution\t% change");
		else
			System.out.println("Mag\tTarget\tUncert\tSolution\t% change\tz-score");
		for (int i=minNonZeroIndex; i<=maxNonZeroIndex; i++) {
			double mag = sectTarget.getX(i);
			double targetVal = sectTarget.getY(i);
			double solVal = solNucl.getY(i);
			
			double fractDiff = (solVal - targetVal)/targetVal;
			
			if (uncertTarget == null) {
				System.out.println((float)mag+"\t"+eDF.format(targetVal)+"\t"+eDF.format(solVal)+"\t"+pDF.format(fractDiff));
			} else {
				double sd = uncertTarget.getStdDev(i);
				double z = (solVal - targetVal)/sd;
				System.out.println((float)mag+"\t"+eDF.format(targetVal)+"\t"+eDF.format(sd)+"\t"+eDF.format(solVal)
					+"\t"+pDF.format(fractDiff)+"\t"+(float)z);
			}
		}
		System.out.println("TOTAL:\t"+eDF.format(sectTarget.calcSumOfY_Vals())+"\t"+(uncertTarget == null ? "" : "\t")
				+eDF.format(solNucl.calcSumOfY_Vals()));
		
		if (outputDir != null) {
			List<IncrementalMagFreqDist> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			if (uncertTarget != null) {
				if (uncertTarget instanceof UncertainBoundedIncrMagFreqDist) {
					uncertTarget = uncertTarget.deepClone();
					uncertTarget.setName("Bounds");
					funcs.add(uncertTarget);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN_TRANS, 1f, Color.ORANGE));
				}
				UncertainBoundedIncrMagFreqDist plusMinus = uncertTarget.estimateBounds(UncertaintyBoundType.ONE_SIGMA);
				plusMinus.setName("+/- sd");
				funcs.add(plusMinus);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN_TRANS, 1f, Color.CYAN));
			}
			
			sectTarget.setName("Target");
			funcs.add(sectTarget);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.CYAN.darker()));
			
			solNucl.setName("Solution");
			funcs.add(solNucl);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
			
			List<DiscretizedFunc> modFuncs = new ArrayList<>();
			List<PlotCurveCharacterstics> modChars = new ArrayList<>();
			
			for (int i=0; i<funcs.size(); i++) {
				IncrementalMagFreqDist func = funcs.get(i);
				PlotCurveCharacterstics pChar = chars.get(i);
				UncertainBoundedIncrMagFreqDist bounds = func instanceof UncertainBoundedIncrMagFreqDist ?
						(UncertainBoundedIncrMagFreqDist)func : null;
				double halfDelta = func.getDelta()*0.5;
				boolean first = true;
				for (int j=0; j<func.size(); j++) {
					double x = func.getX(j);
					double y = func.getY(j);
					if (y > 0 || bounds != null && bounds.getUpperY(j) > 0d) {
						DiscretizedFunc portion = new ArbitrarilyDiscretizedFunc();
						if (first) {
							portion.setName(func.getName());
							first = false;
						}
						portion.set(x-halfDelta, y);
						portion.set(x+halfDelta, y);
						
						if (bounds != null) {
							ArbitrarilyDiscretizedFunc upper = new ArbitrarilyDiscretizedFunc();
							ArbitrarilyDiscretizedFunc lower = new ArbitrarilyDiscretizedFunc();
							double upperY = bounds.getUpperY(j);
							double lowerY = bounds.getLowerY(j);
							upper.set(x-halfDelta, upperY);
							upper.set(x+halfDelta, upperY);
							lower.set(x-halfDelta, lowerY);
							lower.set(x+halfDelta, lowerY);
							portion = new UncertainArbDiscFunc(portion, lower, upper);
						}
						
						modFuncs.add(portion);
						modChars.add(pChar);
					}
				}
			}
			
			PlotSpec spec = new PlotSpec(modFuncs, modChars, sol.getRupSet().getFaultSectionData(sectIndex).getName(),
					"Magnitude", "Incremental Rate");
			spec.setLegendInset(true);
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			double maxVal = 0d;
			double minNonZero = Double.POSITIVE_INFINITY;
			for (DiscretizedFunc func : funcs) {
				for (Point2D pt : func) {
					if (pt.getY() > 1e-12) {
						maxVal = Math.max(maxVal, pt.getY());
						minNonZero = Math.min(minNonZero, pt.getY());
					}
				}
			}
			Range xRange = new Range(0.5*Math.floor(sectTarget.getX(minNonZeroIndex)*2d), 0.5*Math.ceil(sectTarget.getX(maxNonZeroIndex)*2d));
			Range yRange = new Range(Math.pow(10, Math.floor(Math.log10(minNonZero))), Math.pow(10, Math.ceil(Math.log10(maxVal))));
			if (MFD_Y_RANGE != null)
				yRange = MFD_Y_RANGE;
			
			gp.drawGraphPanel(spec, false, true, xRange, yRange);
			
			PlotUtils.writePlots(outputDir, "sect_nucl_mfds_"+sectIndex, gp, 1000, 900, true, true, true);
		}
	}
	
	private static final DecimalFormat eDF = new DecimalFormat("0.00E0");
	private static final DecimalFormat pDF = new DecimalFormat("0.00%");

}
