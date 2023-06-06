package scratch.kevin.nshm23;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.util.Interpolate;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

public class SubSeisReductionTests {

	public static void main(String[] args) {
		double supraMMin = 7.25;
//		double supraMMax = 7.95;
		double supraMMax = 7.35;
//		double supraMMax = supraMMin+0.1;
		double subSeisBValue = 1;
		double supraSeisBValue = 0.5;
		
		double plotMinX = 4.5d;
		
		// don't reduce for anything below this magnitude
		double minReductionMag = 6.05;
		// magnitude at which we transition from sub-b to supra-b, if less than supra-min mag
		double supraBtransMag = 6.55;
		
		// don't let the sub-seismogenic portion exceed this fraction of the total moment
		double maxSubSeisReduction = 1;
		
		// fraction of the sub-seis moment that's actually available from the raw gridded model
		double fractAvailGridded = 1.5;
//		double fractAvailGridded = 0.7;
		
		double offFaultMMax = 7.55;
		double offFaultB = 1d;
		
		double slipRate = 10e-3; // 10 mm/yr
		double area = 15d*7.5d*1e6; // m^2
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(supraMMax+0.2);
		
		int minSubIndex = refMFD.getClosestXIndex(minReductionMag);
		int minMagIndex = refMFD.getClosestXIndex(supraMMin);
		int maxMagIndex = refMFD.getClosestXIndex(supraMMax);
		// it'll be supra-b at and above this mag, sub-b below
		int bTransMagIndex = refMFD.getClosestXIndex(Math.min(supraMMin, supraBtransMag));
		
		double targetMoRate = FaultMomentCalc.getMoment(area, slipRate);
		
		double MIN_MAG = refMFD.getMinX();
		double DELTA_MAG = refMFD.getDelta();
		int NUM_MAG = refMFD.size();
		
		IncrementalMagFreqDist supraSeisMFD, subSeisMFD;
		if (minSubIndex < minMagIndex) {
			// we have a sub-seismogenic reduction
			// start with a full G-R with the supra b-value up to the maximum magnitude
			GutenbergRichterMagFreqDist fullSupraB = new GutenbergRichterMagFreqDist(
					MIN_MAG, maxMagIndex+1, DELTA_MAG, targetMoRate, supraSeisBValue);

			// copy it to a regular MFD:
			IncrementalMagFreqDist sectFullMFD = new IncrementalMagFreqDist(MIN_MAG, maxMagIndex+1, DELTA_MAG);
			for (int i=0; i<fullSupraB.size(); i++)
				sectFullMFD.set(i, fullSupraB.getY(i));

			// now correct the sub-seis portion to have the sub-seis b-value

			// first create a full MFD with the sub b-value. this will only be used in a relative sense
			GutenbergRichterMagFreqDist fullSubB = new GutenbergRichterMagFreqDist(
					MIN_MAG, bTransMagIndex+1, DELTA_MAG);
			fullSubB.setAllButTotCumRate(refMFD.getX(minSubIndex), refMFD.getX(bTransMagIndex), targetMoRate, subSeisBValue);

			double targetFirstSupra = fullSupraB.getY(bTransMagIndex);
			double subFirstSupra = fullSubB.getY(bTransMagIndex);
			for (int i=0; i<bTransMagIndex; i++) {
				double targetRatio = fullSubB.getY(i)/subFirstSupra;
				sectFullMFD.set(i, targetFirstSupra*targetRatio);
			}

			// rescale to match the original moment rate
			sectFullMFD.scaleToTotalMomentRate(targetMoRate);

			// split the target G-R into sub-seismo and supra-seismo parts
			subSeisMFD = new IncrementalMagFreqDist(MIN_MAG, NUM_MAG, DELTA_MAG);
			for (int i=0; i<minMagIndex; i++)
				subSeisMFD.set(i, sectFullMFD.getY(i));

			double subMoRate = subSeisMFD.getTotalMomentRate();
			
			double calcSubReduction = subMoRate/targetMoRate;
			if (calcSubReduction > maxSubSeisReduction) {
				// cap it at the maximum
				subMoRate = targetMoRate*maxSubSeisReduction;
			}
			subSeisMFD.scaleToTotalMomentRate(subMoRate);

			supraSeisMFD = new IncrementalMagFreqDist(MIN_MAG, NUM_MAG, DELTA_MAG);
			for (int i=minMagIndex; i<=maxMagIndex; i++)
				supraSeisMFD.set(i, sectFullMFD.getY(i));
			supraSeisMFD.scaleToTotalMomentRate(targetMoRate - subMoRate);
		} else {
			// no reduction
			GutenbergRichterMagFreqDist fullSupraB = new GutenbergRichterMagFreqDist(
					MIN_MAG, maxMagIndex+1, DELTA_MAG);
			fullSupraB.setAllButTotCumRate(refMFD.getX(minMagIndex), refMFD.getX(maxMagIndex),
					targetMoRate, supraSeisBValue);
			supraSeisMFD = fullSupraB;
			
			subSeisMFD = new IncrementalMagFreqDist(MIN_MAG, maxMagIndex+1, DELTA_MAG);
		}
		double subMoRate = subSeisMFD.getTotalMomentRate();
		double supraMoRate = targetMoRate - subMoRate;
		
		List<IncrementalMagFreqDist> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		DecimalFormat pDF = new DecimalFormat("0.0%");
		subSeisMFD.setName("Target Sub-Seis ("+pDF.format(subMoRate/targetMoRate)+" of fault moment)");
		supraSeisMFD.setName("Target Supra-Seis ("+pDF.format(supraMoRate/targetMoRate)+" of fault moment)");
		
		// original gridded seis MFD, before subtracting out any for fault sub-seis
		double totalGriddedMoment = subMoRate*fractAvailGridded;
		GutenbergRichterMagFreqDist rawGriddedMFD = new GutenbergRichterMagFreqDist(
				MIN_MAG, refMFD.getClosestXIndex(offFaultMMax)+1, DELTA_MAG, totalGriddedMoment, offFaultB);
		rawGriddedMFD.setName("Original Total Gridded");
		
		double momentAvailForSubSeis = 0d;
		for (int i=0; i<rawGriddedMFD.size(); i++)
			if (i<maxMagIndex)
				momentAvailForSubSeis += rawGriddedMFD.getMomentRate(i);
		IncrementalMagFreqDist offFaultMFD, associatedMFD;
		if (momentAvailForSubSeis >= subMoRate) {
			// we have enough moment to cover the sub-seis portion
			associatedMFD = subSeisMFD.deepClone();
			associatedMFD.setName("Associated Sub-Seis (fully satisfied)");
			// now figure out what's leftover for off fault
			double extraMomentBelowSupraMax = momentAvailForSubSeis - subMoRate;
			GutenbergRichterMagFreqDist griddedBelowSupraMax = new GutenbergRichterMagFreqDist(
					MIN_MAG, maxMagIndex+1, DELTA_MAG, extraMomentBelowSupraMax, offFaultB);
			offFaultMFD = new IncrementalMagFreqDist(MIN_MAG, NUM_MAG, DELTA_MAG);
			for (int i=0; i<griddedBelowSupraMax.size(); i++)
				offFaultMFD.set(i, griddedBelowSupraMax.getY(i));
			// copy over any off fault above supra mmax
			for (int i=maxMagIndex+1; i<rawGriddedMFD.size(); i++)
				offFaultMFD.set(i, rawGriddedMFD.getY(i));
			offFaultMFD.setName("Truly Off-Fault MFD");
		} else {
			// we don't have enough moment to cover it, need to exceed the original off fault
			// moment rate to cover sub-seis
			
			Preconditions.checkState(minSubIndex < minMagIndex);
			associatedMFD = new IncrementalMagFreqDist(MIN_MAG, NUM_MAG, DELTA_MAG);
			// fill in below the first interpolation point with the raw sub seis MFD
			for (int i=0; i<minSubIndex; i++)
				associatedMFD.set(i, rawGriddedMFD.getY(i));
			
			double interpX1 = refMFD.getX(minSubIndex);
			double interpX2 = refMFD.getX(bTransMagIndex);
			double interpY1 = rawGriddedMFD.getY(minSubIndex);
			double interpY2 = subSeisMFD.getY(bTransMagIndex) + supraSeisMFD.getY(bTransMagIndex);
			
//			double momentBelow = associatedMFD.getTotalMomentRate();
//			double momentLeft = subMoRate - momentBelow;
//			Preconditions.checkState(momentLeft > 0d);
//			// figure out b-value needed to satisfy that, with it pinned to the incremental rate at the lower bin
//			double prevWedgeB = Double.NaN;
//			double wedgeB = 0d;
//			double fractWedgeMomentDiscrep = 0d;
//			GutenbergRichterMagFreqDist wedgeGR = null;
//			System.out.println("Wedge target moment: "+(float)momentLeft);
//			double wedgeMaxX = refMFD.getX(minMagIndex-1);
//			for (int i=0; i<50; i++) {
//				wedgeGR = new GutenbergRichterMagFreqDist(MIN_MAG, NUM_MAG, DELTA_MAG);
//				wedgeGR.setAllButTotCumRate(interpX1, wedgeMaxX, 1e16, wedgeB); // moment placeholder
//				wedgeGR.scaleToIncrRate(interpIndex, interpY1);
//				double wedgeMoment = wedgeGR.getTotalMomentRate();
//				fractWedgeMomentDiscrep = Math.abs(wedgeMoment - momentLeft)/momentLeft;
//				System.out.println("iter "+i+", wedgeB="+(float)wedgeB+", prevWedgeB="+(float)prevWedgeB+", wedgeMo="+(float)wedgeMoment+", discrep="+(float)fractWedgeMomentDiscrep);
//				double curWedgeB = wedgeB;
//				if (wedgeMoment > momentLeft) {
//					// b-value is too low
//					if (Double.isNaN(prevWedgeB))
//						// first time
//						wedgeB += 0.5;
//					else if (wedgeB > prevWedgeB)
//						// we didn't go far enough, keep going
//						wedgeB += 0.5*(wedgeB - prevWedgeB);
//					else
//						// we overcorrected, split the difference
//						wedgeB = 0.5*(wedgeB + prevWedgeB);
//				} else {
//					// b-value is too high
//					if (Double.isNaN(prevWedgeB))
//						// first time
//						wedgeB -= 0.5;
//					else if (wedgeB < prevWedgeB)
//						// we didn't go far enough, keep going
//						wedgeB -= 0.5*(prevWedgeB - wedgeB);
//					else
//						// we overcorrected, split the difference
//						wedgeB = 0.5*(wedgeB + prevWedgeB);
//				}
//				if ((float)wedgeB == (float)curWedgeB)
//					// we have converged
//					break;
//				prevWedgeB = curWedgeB;
//			}
//			System.out.println(wedgeGR);
//			System.out.println("Found wedgeB="+wedgeB+", fractional discrepancy is "+(float)fractWedgeMomentDiscrep);
//			// now extend the wedge to find out what the rate would be at the min supra mag
//			wedgeGR = new GutenbergRichterMagFreqDist(MIN_MAG, NUM_MAG, DELTA_MAG);
//			wedgeGR.setAllButTotCumRate(interpX1, interpX2, 1e16, wedgeB); // moment placeholder
//			wedgeGR.scaleToIncrRate(interpIndex, interpY1);
//			double interpY2 = Math.min(wedgeGR.getY(minMagIndex), supraSeisMFD.getY(minMagIndex));
			
			// interpolate between interp min and supra min
			double logInterpY1 = Math.log10(interpY1);
//			double interpY2 = supraSeisMFD.getY(minMagIndex);
			double logInterpY2 = Math.log10(interpY2);
			for (int i=minSubIndex; i<bTransMagIndex; i++) {
				double logY = Interpolate.findY(interpX1, logInterpY1, interpX2, logInterpY2, refMFD.getX(i));
				
				associatedMFD.set(i, Math.pow(10, logY));
			}
			// copy over any extra above the b-transition magnitude up to the actual min mag
			for (int i=bTransMagIndex; i<minMagIndex; i++)
				associatedMFD.set(i, subSeisMFD.getY(i));
			double associatedMoment = associatedMFD.getTotalMomentRate();
			associatedMFD.setName("Associated Sub-Seis ("+pDF.format(associatedMoment/subMoRate)+" satisfied)");
			
			// truly off fault MFD will be empty unless off fault mmax > supra mmax
			// copy over any off fault above supra mmax
			offFaultMFD = new IncrementalMagFreqDist(MIN_MAG, NUM_MAG, DELTA_MAG);
			for (int i=maxMagIndex+1; i<rawGriddedMFD.size(); i++)
				offFaultMFD.set(i, rawGriddedMFD.getY(i));
			offFaultMFD.setName("Truly Off-Fault MFD");
		}
		
		funcs.add(rawGriddedMFD);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_SQUARE, 3f, Color.GRAY));
		funcs.add(subSeisMFD);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_SQUARE, 3f, Color.BLUE));
		funcs.add(offFaultMFD);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_SQUARE, 3f, Color.ORANGE));
		funcs.add(associatedMFD);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 1.5f, Color.CYAN));
		funcs.add(supraSeisMFD);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_SQUARE, 3f, Color.RED));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "MFD Test", "Magnitude", "Incremental Nucleation Rate");
		spec.setLegendInset(true);
		
		double plotMaxX = Math.ceil(Math.max(supraMMax, offFaultMMax)*5d)/5d;
		Range xRange = new Range(plotMinX, plotMaxX);
		double plotMaxY = 0d;
		double plotMinY = Double.POSITIVE_INFINITY;
		for (IncrementalMagFreqDist func : funcs) {
			for (Point2D pt : func) {
				if (pt.getY() > 0 && xRange.contains(pt.getX())) {
					plotMaxY = Math.max(plotMaxY, pt.getY());
					plotMinY = Math.min(plotMinY, pt.getY());
				}
			}
		}
		plotMaxY = Math.pow(10, Math.ceil(Math.log10(plotMaxY)));
		plotMinY = Math.pow(10, Math.floor(Math.log10(plotMinY)));
		
		Range yRange = new Range(plotMinY, plotMaxY);
		
		GraphWindow gw = new GraphWindow(spec);
		gw.setAxisRange(xRange, yRange);
		gw.setYLog(true);
		gw.setVisible(true);
		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
	}

}
