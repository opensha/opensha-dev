package scratch.kevin.prvi25.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertaintyBoundType;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateFileLoader;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateFileLoader.Direct;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateFileLoader.RateType;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateModel;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SeismicityRateEpoch;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader.PRVI25_SeismicityRegions;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;

public class RateModelComparison {

	public static void main(String[] args) throws IOException {
		File ratesDir = new File("/home/kevin/workspace/opensha/src/main/resources/data/erf/prvi25/seismicity/rates");
		File prevDir = new File(ratesDir, "2025_03_26");
		File newDir1900 = new File(ratesDir, PRVI25_CrustalSeismicityRate.RATE_DATE+"/"+PRVI25_SeismicityRateEpoch.FULL.getRateSubDirName());
//		File newDir1900 = new File(ratesDir, PRVI25_CrustalSeismicityRate.RATE_DATE+"/1973_scaled_to_1900");
		File newDir1973 = new File(ratesDir, PRVI25_CrustalSeismicityRate.RATE_DATE+"/"+PRVI25_SeismicityRateEpoch.RECENT.getRateSubDirName());
		File directDir = new File(ratesDir, PRVI25_CrustalSeismicityRate.RATE_DATE+"/direct");
		
		File outputDir = new File("/tmp/prvi_seis_rates");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		Range magRange = new Range(5d, 8d);
		Range yRange = new Range(1e-5, 1e1);
		
		double mfdMmax = 10d;
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(5.01, mfdMmax-0.01);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.getPlotPrefs().setLegendFontSize(16);
		
		List<PlotSpec> allIncrs = new ArrayList<>();
		List<PlotSpec> allCmls = new ArrayList<>();
		
		DecimalFormat pDF = new DecimalFormat("0%");
		
		List<Double> prevRates = new ArrayList<>();
		List<Double> newRates1900 = new ArrayList<>();
		List<Double> newRates1973 = new ArrayList<>();
		List<Double> directRates1900 = new ArrayList<>();
		List<Double> directRates1973 = new ArrayList<>();
		double directSum1900 = 0d;
		double directSum1973 = 0d;
		
		PRVI25_SeismicityRegions[] regions = PRVI25_SeismicityRegions.values();
		List<String> regionNames = new ArrayList<>();
		for (PRVI25_SeismicityRegions reg : regions) {
			CSVFile<String> origCSV = CSVFile.readFile(new File(prevDir, reg.name()+".csv"), false);
			SeismicityRateModel prevModel = new SeismicityRateModel(origCSV, RateType.M1_TO_MMAX, UncertaintyBoundType.CONF_95);
			
			CSVFile<String> newCSV1900 = CSVFile.readFile(new File(newDir1900, reg.name()+".csv"), false);
			SeismicityRateModel newModel1900 = new SeismicityRateModel(newCSV1900, RateType.M1_TO_MMAX, UncertaintyBoundType.CONF_95);
			
			CSVFile<String> newCSV1973 = CSVFile.readFile(new File(newDir1973, reg.name()+".csv"), false);
			SeismicityRateModel newModel1973 = new SeismicityRateModel(newCSV1973, RateType.M1_TO_MMAX, UncertaintyBoundType.CONF_95);
			
			List<IncrementalMagFreqDist> incrFuncs = new ArrayList<>();
			List<PlotCurveCharacterstics> incrChars = new ArrayList<>();
			List<EvenlyDiscretizedFunc> cmlFuncs = new ArrayList<>();
			List<PlotCurveCharacterstics> cmlChars = new ArrayList<>();
			
			Color color;
			String name, regName;
			switch (reg) {
			case CAR_INTERFACE:
				color = Colors.tab_orange;
				name = "Caribbean Interface";
				regName = "CAR Interface";
				break;
			case CAR_INTRASLAB:
				color = Colors.tab_lightorange;
				name = "Caribbean Intraslab";
				regName = "CAR Intraslab";
				break;
			case CRUSTAL:
				color = Colors.tab_blue;
				name = "Crustal";
				regName = "Crustal";
				break;
			case MUE_INTERFACE:
				color = Colors.tab_green;
				name = "Muertos Interface";
				regName = "MUE Interface";
				break;
			case MUE_INTRASLAB:
				color = Colors.tab_lightgreen;
				name = "Muertos Intraslab";
				regName = "MUE Intraslab";
				break;

			default:
				throw new IllegalStateException("Unknown region: "+reg);
			}
			
			regionNames.add(name);
			
			for (boolean is1973 : new boolean[] {false,true}) {
				PlotCurveCharacterstics boundsChar = new PlotCurveCharacterstics(is1973 ? PlotLineType.DASHED : PlotLineType.SOLID, 2f, color);
				PlotCurveCharacterstics prefChar = new PlotCurveCharacterstics(is1973 ? PlotLineType.DASHED : PlotLineType.SOLID, 5f, color);
				
				SeismicityRateModel newModel = is1973 ? newModel1973 : newModel1900;
				
				String years = is1973 ? "1973-2023" : "1900-2023";
				
				IncrementalMagFreqDist lower = newModel.buildLower(refMFD, mfdMmax);
				IncrementalMagFreqDist pref = newModel.buildPreferred(refMFD, mfdMmax);
				IncrementalMagFreqDist upper = newModel.buildUpper(refMFD, mfdMmax);
				
				pref.setName(years+" "+pref.getName().replace("Observed", "Preferred"));
				incrFuncs.add(pref);
				incrChars.add(prefChar);
				cmlFuncs.add(pref.getCumRateDistWithOffset());
				cmlChars.add(incrChars.get(incrChars.size()-1));
				
				lower.setName("95% M1-Mmax Bounds");
				incrFuncs.add(lower);
				incrChars.add(boundsChar);
				cmlFuncs.add(lower.getCumRateDistWithOffset());
				cmlChars.add(incrChars.get(incrChars.size()-1));
				
				upper.setName(null);
				incrFuncs.add(upper);
				incrChars.add(boundsChar);
				cmlFuncs.add(upper.getCumRateDistWithOffset());
				cmlChars.add(incrChars.get(incrChars.size()-1));
				
//				EvenlyDiscretizedFunc lowerExact = ((Exact)newExactModel.getLowerRecord()).cumulativeDist;
//				lowerExact.setName("95% Exact Bounds");
////				lowerExact.setName(null);
//				cmlFuncs.add(lowerExact);
//				cmlChars.add(exactBoundsChar);
//				
//				EvenlyDiscretizedFunc upperExact = ((Exact)newExactModel.getUpperRecord()).cumulativeDist;
//				upperExact.setName(null);
////				upperExact.setName("95% Exact Bounds");
//				cmlFuncs.add(upperExact);
//				cmlChars.add(exactBoundsChar);
				
				if (is1973)
					newRates1973.add(newModel.getMeanRecord().rateAboveM1);
				else
					newRates1900.add(newModel.getMeanRecord().rateAboveM1);
			}
			
			IncrementalMagFreqDist lower = prevModel.buildLower(refMFD, mfdMmax);
			IncrementalMagFreqDist pref = prevModel.buildPreferred(refMFD, mfdMmax);
			IncrementalMagFreqDist upper = prevModel.buildUpper(refMFD, mfdMmax);
			
			prevRates.add(prevModel.getMeanRecord().rateAboveM1);
			
			color = Color.GRAY;
			PlotCurveCharacterstics boundsChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, color);
			PlotCurveCharacterstics prefChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, color);
			
			pref.setName("Prior "+pref.getName().replace("Observed", "Preferred"));
			incrFuncs.add(pref);
			incrChars.add(prefChar);
			cmlFuncs.add(pref.getCumRateDistWithOffset());
			cmlChars.add(incrChars.get(incrChars.size()-1));
			
//			lower.setName(null);
//			incrFuncs.add(lower);
//			incrChars.add(boundsChar);
//			cmlFuncs.add(lower.getCumRateDistWithOffset());
//			cmlChars.add(incrChars.get(incrChars.size()-1));
//			
//			upper.setName(null);
//			incrFuncs.add(upper);
//			incrChars.add(boundsChar);
//			cmlFuncs.add(upper.getCumRateDistWithOffset());
//			cmlChars.add(incrChars.get(incrChars.size()-1));
			
			for (boolean is1973 : new boolean[] {false,true}) {
				File ratesFile = new File(directDir, PRVI25_CrustalSeismicityRate.getDirectRateFileName(reg,
						is1973 ? PRVI25_SeismicityRateEpoch.RECENT : PRVI25_SeismicityRateEpoch.FULL));
				
				CSVFile<String> totalRateCSV = CSVFile.readFile(ratesFile, false);
				List<Direct> directs = SeismicityRateFileLoader.loadDirectBranches(totalRateCSV);
				
				double minFuncMag = is1973 ? 5.01 : 6.01;
				double maxFuncMag = directs.get(0).maxObsIncrMag-0.01;
				if (maxFuncMag <= 0d) {
					// not found
//					maxFuncMag = directs.get(0).cumulativeDist.getMaxX()+0.01;
					maxFuncMag = directs.get(0).cumulativeDist.getMaxX()-0.01;
				}

				IncrementalMagFreqDist obsRefMFD = FaultSysTools.initEmptyMFD(minFuncMag, maxFuncMag);
				EvenlyDiscretizedFunc refCml = obsRefMFD.getCumRateDistWithOffset();
				
				IncrementalMagFreqDist meanIncrObs = new IncrementalMagFreqDist(obsRefMFD.getMinX(), obsRefMFD.getMaxX(), obsRefMFD.size());
				IncrementalMagFreqDist incr2p5 = new IncrementalMagFreqDist(obsRefMFD.getMinX(), obsRefMFD.getMaxX(), obsRefMFD.size());
				IncrementalMagFreqDist incr16 = new IncrementalMagFreqDist(obsRefMFD.getMinX(), obsRefMFD.getMaxX(), obsRefMFD.size());
				IncrementalMagFreqDist incr84 = new IncrementalMagFreqDist(obsRefMFD.getMinX(), obsRefMFD.getMaxX(), obsRefMFD.size());
				IncrementalMagFreqDist incr97p5 = new IncrementalMagFreqDist(obsRefMFD.getMinX(), obsRefMFD.getMaxX(), obsRefMFD.size());
				
				IncrementalMagFreqDist[] dataIncrs = {meanIncrObs, incr2p5, incr16, incr84, incr97p5};
				
				EvenlyDiscretizedFunc meanObsCml = new EvenlyDiscretizedFunc(refCml.getMinX(), refCml.getMaxX(), refCml.size());
				EvenlyDiscretizedFunc cml2p5 = new EvenlyDiscretizedFunc(refCml.getMinX(), refCml.getMaxX(), refCml.size());
				EvenlyDiscretizedFunc cml16 = new EvenlyDiscretizedFunc(refCml.getMinX(), refCml.getMaxX(), refCml.size());
				EvenlyDiscretizedFunc cml84 = new EvenlyDiscretizedFunc(refCml.getMinX(), refCml.getMaxX(), refCml.size());
				EvenlyDiscretizedFunc cml97p5 = new EvenlyDiscretizedFunc(refCml.getMinX(), refCml.getMaxX(), refCml.size());
				
				EvenlyDiscretizedFunc[] dataCmls = {meanObsCml, cml2p5, cml16, cml84, cml97p5};
				
				for (int i=0; i<meanIncrObs.size(); i++) {
					double mag = meanIncrObs.getX(i);
					Preconditions.checkState((float)mag >= (float)directs.get(0).M1, "mag=%s, M1=%s", (float)mag, (float)directs.get(0).M1);
					if ((float)mag > (float)directs.get(0).incrementalDist.getMaxX()) {
						break;
					}
					Preconditions.checkState(directs.size() == dataIncrs.length);
					for (int j=0; j<directs.size(); j++) {
						EvenlyDiscretizedFunc directIncr = directs.get(j).incrementalDist;
						EvenlyDiscretizedFunc directCml = directs.get(j).cumulativeDist;
						
						dataIncrs[j].set(i, directIncr.getY(directIncr.getClosestXIndex(mag)));
						dataCmls[j].set(i, directCml.getY(directCml.getClosestXIndex(mag-0.05)));
					}
				}
				
				PlotCurveCharacterstics obsChar;
				String obsName;
//				Color darkishGray = new Color(80, 80, 80);
				Color color1900 = Colors.tab_brown;
				Color color1973 = Color.BLACK;
				if (is1973) {
					obsName = "Observed (1973-2023)";
					obsChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, color1973);
				} else {
					obsName = "Observed (1900-2023)";
					obsChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, color1900);
				}
				meanIncrObs.setName(obsName);
				incrFuncs.add(meanIncrObs);
				incrChars.add(obsChar);
				meanObsCml.setName(meanIncrObs.getName());
				cmlFuncs.add(meanObsCml);
				cmlChars.add(obsChar);
				
				if (is1973)
					directRates1973.add(meanObsCml.getY(meanObsCml.getClosestXIndex(5.01)));
				else
					directRates1900.add(meanObsCml.getY(meanObsCml.getClosestXIndex(5.01)));
			}
			
			Preconditions.checkState(incrFuncs.size() == incrChars.size());
			Preconditions.checkState(cmlFuncs.size() == cmlChars.size());
			
//			PlotSpec incrPlot = new PlotSpec(incrFuncs, incrChars, name, "Magnitude", "Incremental Rate (1/yr)");
//			cmlPlot.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
////			cmlPlot.setLegendInset(RectangleAnchor.TOP_RIGHT);
//			allIncrs.add(incrPlot);
//			
//			gp.drawGraphPanel(incrPlot, false, true, magRange, yRange);
//			
//			PlotUtils.writePlots(outputDir, reg.name(), gp, 700, 600, true, true, false);
			
			PlotSpec cmlPlot = new PlotSpec(cmlFuncs, cmlChars, name, "Magnitude", "Cumulative Rate (1/yr)");
			cmlPlot.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
//			cmlPlot.setLegendInset(RectangleAnchor.TOP_RIGHT);
//			cmlPlot.setLegendInset(RectangleAnchor.TOP, 0.5, 0.95, 0.9, false);
			allCmls.add(cmlPlot);
			
			gp.drawGraphPanel(cmlPlot, false, true, magRange, yRange);
			
			PlotUtils.writePlots(outputDir, reg.name()+"_cml", gp, 700, 600, true, true, false);
		}
		
		for (boolean is1973 : new boolean[] {false,true}) {
			if (is1973)
				System.out.println("Rate model, updated 1973-2023 vs previous 1900-2023");
			else
				System.out.println("Rate model, updated 1900-2023 vs previous 1900-2023");
			List<Double> newRates = is1973 ? newRates1973 : newRates1900;
			for (int r=0; r<regions.length; r++) {
				String name = regionNames.get(r);
				
				double prevRate = prevRates.get(r);
				double newRate = newRates.get(r);
				String changeStr = pDF.format((newRate-prevRate)/prevRate);
				if (newRate > prevRate)
					changeStr = "+"+changeStr;
				System.out.println("\t"+name+" M>5 rate:\t"+prevRate+" -> "+newRate+" ("+changeStr+")");
			}
			System.out.println();
		}
		
		for (boolean is1973 : new boolean[] {false,true}) {
			if (is1973)
				System.out.println("Direct rates, 1973-2023");
			else
				System.out.println("Direct rates, 1900-2023");
			List<Double> rates = is1973 ? directRates1973 : directRates1900;
			double sum = rates.stream().mapToDouble(D->D).sum();
			Preconditions.checkState(rates.size()==regions.length);
			System.out.println("\tSum: "+(float)sum);
			for (int r=0; r<regions.length; r++) {
				String name = regionNames.get(r);
				
				double rate = rates.get(r);
				System.out.println("\t"+name+" M>5 rate:\t"+(float)rate+" ("+pDF.format(rate/sum)+")");
			}
			System.out.println();
		}
	}

}
