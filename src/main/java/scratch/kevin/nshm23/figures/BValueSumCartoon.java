package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.data.Range;
import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertainIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertaintyBoundType;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SupraSeisBValues;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;

class BValueSumCartoon {
	
	private static boolean CARTOON = false;

	public static void main(String[] args) throws IOException {
		Range magRange = CARTOON ? new Range(6d, 8d) : new Range (6d, 8.5d);
		
		double[] bVals = { 0d, 0.5d, 1d };
		
		EvenlyDiscretizedFunc refMFD = SupraSeisBValInversionTargetMFDs.buildRefXValues(magRange.getUpperBound()+0.1);
		
		List<List<IncrementalMagFreqDist>> bValSectGRs = new ArrayList<>();
		List<IncrementalMagFreqDist> bValSums = new ArrayList<>();
		List<EvenlyDiscretizedFunc> bValSectCounts = new ArrayList<>();
		for (int i=0; i<bVals.length; i++) {
			bValSectGRs.add(new ArrayList<>());
			bValSums.add(new SummedMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta()));
			bValSectCounts.add(refMFD.deepClone());
		}
		
		boolean isParents = false;
		
		double minNonZero = Double.POSITIVE_INFINITY;
		double maxVal = 0;
		
		File outputDir;
		if (CARTOON) {
			outputDir = new File("/home/kevin/Documents/papers/2023_NSHM23_Inversion/figures/b_val_sums_cartoon");
			
			int numSects = 100;
			Range slipRateRange = new Range(1d, 15d);
			Range minMagRange = new Range(6d, 6.5d);
			Range maxMagRange = new Range(7.2d, 8d);
			
			Random rand = new Random(numSects);
			
			double areaKM = 14d*8d;
			double area = areaKM*1e6;
			
			for (int s=0; s<numSects; s++) {
				// randomly sample mag range
				double minMag = minMagRange.getLowerBound() + rand.nextDouble()*minMagRange.getLength();
				double maxMag = maxMagRange.getLowerBound() + rand.nextDouble()*maxMagRange.getLength();
				
				// round them
				minMag = refMFD.getX(refMFD.getClosestXIndex(minMag));
				maxMag = refMFD.getX(refMFD.getClosestXIndex(maxMag));
				System.out.println("Section "+s);
				System.out.println("\tMag range: ["+(float)minMag+","+(float)maxMag+"]");
				
				// randomly sample slip rate
				double slipRate = slipRateRange.getLowerBound() + rand.nextDouble()*slipRateRange.getLength();
				System.out.println("\tSlip rate: "+(float)slipRate);
				
				double moRate = FaultMomentCalc.getMoment(area, slipRate*1e-3);
				System.out.println("\tMoment rate: "+(float)moRate);
				
				for (int b=0; b<bVals.length; b++) {
					GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
					gr.setAllButTotCumRate(minMag, maxMag, moRate, bVals[b]);
					minNonZero = Math.min(minNonZero, gr.getY(gr.getClosestXIndex(maxMag)));
					
					bValSectGRs.get(b).add(gr);
					((SummedMagFreqDist)bValSums.get(b)).addIncrementalMagFreqDist(gr);
					maxVal = Math.max(maxVal, bValSums.get(b).getMaxY());
					
					for (int i=0; i<refMFD.size(); i++)
						if (gr.getY(i) > 0)
							bValSectCounts.get(b).add(i, 1d);
				}
			}
		} else {
			outputDir = new File("/home/kevin/Documents/papers/2023_NSHM23_Inversion/figures/b_val_sums");
			
			isParents = true;
			
			LogicTreeBranch<LogicTreeNode> branch = NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT.copy();
			
//			NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory.NoIncompatibleDataAdjust();
			NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
			factory.setCacheDir(new File("/home/kevin/OpenSHA/nshm23/rup_sets/cache"));
			
			for (int b=0; b<bVals.length; b++) {
				SupraSeisBValues bVal = null;
				for (SupraSeisBValues test : SupraSeisBValues.values()) {
					if (test.bValue == bVals[b]) {
						bVal = test;
						break;
					}
				}
				Preconditions.checkNotNull(bVal);
				branch.setValue(bVal);
				FaultSystemRupSet rupSet = factory.buildRuptureSet(branch, 32);
				SupraSeisBValInversionTargetMFDs targetMFDs = rupSet.requireModule(SupraSeisBValInversionTargetMFDs.class);
				
				for (IncrementalMagFreqDist mfd : targetMFDs.getOnFaultSupraSeisNucleationMFDs()) {
					for (int i=0; i<mfd.size(); i++)
						if (mfd.getY(i) > 0d)
							bValSectCounts.get(b).add(i, 1d);
				}
				
				List<? extends IncrementalMagFreqDist> mfds;
				if (isParents) {
					Map<Integer, SummedMagFreqDist> parentMFDs = new HashMap<>();
					
					List<UncertainIncrMagFreqDist> sectMFDs = targetMFDs.getOnFaultSupraSeisNucleationMFDs();
					for (int s=0; s<sectMFDs.size(); s++) {
						int parentID = rupSet.getFaultSectionData(s).getParentSectionId();
						
						IncrementalMagFreqDist mfd = sectMFDs.get(s);
						
						SummedMagFreqDist parentMFD = parentMFDs.get(parentID);
						if (parentMFD == null) {
							parentMFD = new SummedMagFreqDist(mfd.getMinX(), mfd.size(), mfd.getDelta());
							parentMFDs.put(parentID, parentMFD);
						}
						parentMFD.addIncrementalMagFreqDist(mfd);
					}
					mfds = new ArrayList<>(parentMFDs.values());
				} else {
					mfds = targetMFDs.getOnFaultSupraSeisNucleationMFDs();
				}
				for (IncrementalMagFreqDist mfd : mfds) {
					for (Point2D pt : mfd) {
						if (pt.getY() > 0) {
							minNonZero = Math.min(minNonZero, pt.getY());
							maxVal = Math.max(maxVal, pt.getY());
						}
					}
					
					bValSectGRs.get(b).add(mfd);
//					bValSums.get(b).addIncrementalMagFreqDist(mfd);
					maxVal = Math.max(maxVal, bValSums.get(b).getMaxY());
				}
				bValSums.set(b, targetMFDs.getTotalOnFaultSupraSeisMFD());
			}
		}
		
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		Range yRange = new Range(Math.pow(10, Math.floor(Math.log10(minNonZero))),
				Math.pow(10, Math.ceil(Math.log10(maxVal))));
		
		double maxCount = 0d;
		for (EvenlyDiscretizedFunc countFunc : bValSectCounts)
			maxCount = Math.max(maxCount, countFunc.getMaxY());
		System.out.println("Max Count: "+maxCount);
		maxCount += 1;
		if (maxCount > 500)
			maxCount = Math.ceil(maxCount/100d)*100d;
		else
			maxCount = Math.ceil(maxCount/10d)*10d;
		
		Range countRange = new Range(0d, maxCount);
		
		for (int b=0; b<bVals.length; b++) {
			double bVal = bVals[b];
			
			List<DiscretizedFunc> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			IncrementalMagFreqDist bValSum = bValSums.get(b);
			bValSum.setName("Regional Sum");
			funcs.add(bValSum);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
			
			if (bValSum instanceof UncertainIncrMagFreqDist) {
				UncertainBoundedIncrMagFreqDist bounded;
				if (bValSum instanceof UncertainBoundedIncrMagFreqDist)
					bounded = (UncertainBoundedIncrMagFreqDist)bValSum;
				else
					bounded = ((UncertainIncrMagFreqDist)bValSum).estimateBounds(UncertaintyBoundType.ONE_SIGMA);
				System.out.println("b="+(float)bVal+" MFD & bounds:");
				boolean first = true;
				for (int i=0; i<bounded.size(); i++) {
					double x = bounded.getX(i);
					double y = bounded.getY(i);
					if (y > 0 || !first) {
						double up = bounded.getUpperY(i);
						double low = bounded.getLowerY(i);
						first = false;
						System.out.println(i+".\tx="+(float)x+"\ty="+(float)y+"\tbounds=["+(float)low+", "+(float)up+"]");
					}
				}
				// adjust upper to not fall off a cliff above max val
				IncrementalMagFreqDist upper = bounded.getUpper();
				int lowestZeroIndex = -1;
				for (int i=upper.size(); --i>=0;) {
					if (upper.getY(i) <= yRange.getLowerBound())
						lowestZeroIndex = i;
					else
						break;
				}
				if (lowestZeroIndex > 0) {
					upper = upper.deepClone();
					upper.set(lowestZeroIndex, yRange.getLowerBound());
					System.out.println("Setting upper for M"+(float)upper.getX(lowestZeroIndex)+" to "+(float)yRange.getLowerBound());
				}
				// remove the name
				bounded = new UncertainBoundedIncrMagFreqDist(bValSum, bounded.getLower(),
						upper, bounded.getBoundType(), bounded.getStdDevs());
				bounded.setName(bounded.getBoundName());
				funcs.add(bounded);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(0, 0, 255, 60)));
			}
			
			int numSects = bValSectGRs.get(b).size();
			int sectAlpha;
			if (numSects > 5000)
				sectAlpha = 60;
			else if (numSects > 1000)
				sectAlpha = 120;
			else if (numSects > 500)
				sectAlpha = 180;
			else
				sectAlpha = 255;
			Color sectColor = new Color(127, 127, 127, sectAlpha);
			
			for (IncrementalMagFreqDist gr : bValSectGRs.get(b))
				addMFDFuncs(funcs, chars, gr, new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, sectColor), bVal);
			
			funcs.get(funcs.size()-1).setName(isParents ? "Section MFDs" : "Subsection MFDs");
			
			// copy it on top without a name
			IncrementalMagFreqDist copy = bValSum.deepClone();
			copy.setName(null);
			funcs.add(copy);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
			
			PlotSpec spec = new PlotSpec(funcs, chars, "b="+(float)bVal, "Magnitude", "Incremental Nucleation Rate (1/yr)");
			spec.setLegendInset(true);
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			gp.setTickLabelFontSize(22);
			gp.setAxisLabelFontSize(26);
			
			gp.drawGraphPanel(spec, false, true, magRange, yRange);
			
			String prefix = "supra_seis_b_"+(float)bVal;
			
			PlotUtils.writePlots(outputDir, prefix, gp, 800, 750, true, true, false);
			
			// now write combined plot with section counts
			funcs = new ArrayList<>();
			chars = new ArrayList<>();
			
			funcs.add(bValSectCounts.get(b));
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GRAY));
			
			PlotSpec countSpec = new PlotSpec(funcs, chars, spec.getTitle(), spec.getXAxisLabel(), "Subsection Count");
			
			gp.drawGraphPanel(List.of(spec, countSpec), List.of(false), List.of(true, false),
					List.of(magRange), List.of(yRange, countRange));
			
			PlotUtils.setSubPlotWeights(gp, 10, 3);
			PlotUtils.writePlots(outputDir, prefix+"_count", gp, 800, 1100, true, true, false);
		}
	}
	
	private static void addMFDFuncs(List<? super ArbitrarilyDiscretizedFunc> funcs, List<PlotCurveCharacterstics> chars,
			EvenlyDiscretizedFunc mfd, PlotCurveCharacterstics pChar, double bVal) {
		double halfDelta = 0.5*mfd.getDelta();
		double halfDeltaScalar = 1d/Math.pow(10, -bVal*halfDelta);
		
		double prevX = Double.NaN, prevY = Double.NaN;
		ArbitrarilyDiscretizedFunc prevFunc = null;
		for (Point2D pt : mfd) {
			if (pt.getY() > 0) {
				double x1 = pt.getX() - halfDelta;
				double x2 = pt.getX() + halfDelta;
				double y1 = pt.getY()*halfDeltaScalar;
				double y2 = pt.getY()/halfDeltaScalar;
				ArbitrarilyDiscretizedFunc func;
				if (prevFunc != null && (float)x1 == (float)prevX && (float)y1 == (float)prevY) {
					func = prevFunc;
					func.set(x2, y2);
				} else {
					func = new ArbitrarilyDiscretizedFunc();
					func.set(x1, y1);
					func.set(x2, y2);
					funcs.add(func);
					chars.add(pChar);
				}
				
				prevFunc = func;
				prevX = x2;
				prevY = y2;
			} else {
				prevFunc = null;
			}
		}
	}

}
