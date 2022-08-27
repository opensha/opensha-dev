package scratch.kevin.nshm23;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.jfree.data.Range;
import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.base.Preconditions;

public class BValueSumCartoon {

	public static void main(String[] args) throws IOException {
		Range magRange = new Range(6d, 8d);
		EvenlyDiscretizedFunc refMFD = SupraSeisBValInversionTargetMFDs.buildRefXValues(magRange.getUpperBound());
		
		int numSects = 100;
		Range slipRateRange = new Range(1d, 15d);
		Range minMagRange = new Range(6d, 6.5d);
		Range maxMagRange = new Range(7.2d, 8d);
		
		Random rand = new Random(numSects);
		
		double areaKM = 14d*8d;
		double area = areaKM*1e6;
		
		double[] bVals = { 0d, 0.5d, 1d };
		
		File outputDir = new File("/tmp/b_val_sums");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		List<List<GutenbergRichterMagFreqDist>> bValSectGRs = new ArrayList<>();
		List<SummedMagFreqDist> bValSums = new ArrayList<>();
		for (int i=0; i<bVals.length; i++) {
			bValSectGRs.add(new ArrayList<>());
			bValSums.add(new SummedMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta()));
		}
		
		double minNonZero = Double.POSITIVE_INFINITY;
		double maxVal = 0;
		
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
				bValSums.get(b).addIncrementalMagFreqDist(gr);
				maxVal = Math.max(maxVal, bValSums.get(b).getMaxY());
			}
		}
		
		Range yRange = new Range(Math.pow(10, Math.floor(Math.log10(minNonZero))),
				Math.pow(10, Math.ceil(Math.log10(maxVal))));
		
		for (int b=0; b<bVals.length; b++) {
			double bVal = bVals[b];
			
			List<DiscretizedFunc> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			for (GutenbergRichterMagFreqDist gr : bValSectGRs.get(b))
				addMFDFuncs(funcs, chars, gr, new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GRAY), bVal);
			
			funcs.get(0).setName("Subsection MFDs");
			
			addMFDFuncs(funcs, chars, bValSums.get(b), new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK), bVal);
			
			PlotSpec spec = new PlotSpec(funcs, chars, "b="+(float)bVal, "Magnitude", "Incremental Nucleation Rate (1/yr)");
			spec.setLegendVisible(true);
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			gp.drawGraphPanel(spec, false, true, magRange, yRange);
			
			PlotUtils.writePlots(outputDir, "supra_seis_b_"+(float)bVal, gp, 1000, 800, true, true, false);
		}
	}
	
	private static void addMFDFuncs(List<? super ArbitrarilyDiscretizedFunc> funcs, List<PlotCurveCharacterstics> chars,
			EvenlyDiscretizedFunc mfd, PlotCurveCharacterstics pChar, double bVal) {
		double halfDelta = 0.5*mfd.getDelta();
		double halfDeltaScalar = 1d/Math.pow(10, -bVal*halfDelta);
		
		for (Point2D pt : mfd) {
			if (pt.getY() > 0) {
				double x1 = pt.getX() - 0.5*mfd.getDelta();
				double x2 = pt.getX() + 0.5*mfd.getDelta();
				double y1 = pt.getY()*halfDeltaScalar;
				double y2 = pt.getY()/halfDeltaScalar;
				ArbitrarilyDiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
				func.set(x1, y1);
				func.set(x2, y2);
				funcs.add(func);
				chars.add(pChar);
			}
		}
	}

}
