package scratch.kevin.prvi25.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateFileLoader;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateFileLoader.Direct;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateModel;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SeismicityRateEpoch;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionCaribbeanSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionMuertosSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader.PRVI25_SeismicityRegions;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;

public class RateEpochComparison {

	public static void main(String[] args) throws IOException {
		File ratesDir = new File("/home/kevin/workspace/opensha/src/main/resources/data/erf/prvi25/seismicity/rates");
		File directDir = new File(ratesDir, PRVI25_CrustalSeismicityRate.RATE_DATE+"/direct");
		
		File outputDir = new File("/tmp/prvi_seis_epochs");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		Range magRange = new Range(5d, 8d);
		Range yRange = new Range(1e-5, 1e1);
		
		double mfdMmax = 10d;
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(5.01, mfdMmax-0.01);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.getPlotPrefs().setLegendFontSize(16);
		
		PRVI25_SeismicityRateEpoch[] epochs =  { PRVI25_SeismicityRateEpoch.FULL,
				PRVI25_SeismicityRateEpoch.RECENT, PRVI25_SeismicityRateEpoch.RECENT_SCALED };
		PlotLineType[] epochLines = { PlotLineType.SOLID, PlotLineType.DASHED, PlotLineType.DOTTED };
		
		DecimalFormat pDF = new DecimalFormat("0%");
		
		double sumNobs1900 = 0d;
		double sumNobs1973 = 0d;
		
		PRVI25_SeismicityRegions[] regions = PRVI25_SeismicityRegions.values();
		for (PRVI25_SeismicityRegions reg : regions) {
			List<IncrementalMagFreqDist> incrFuncs = new ArrayList<>();
			List<PlotCurveCharacterstics> incrChars = new ArrayList<>();
			List<DiscretizedFunc> cmlFuncs = new ArrayList<>();
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
			
			List<UncertainBoundedIncrMagFreqDist> epochBounds = new ArrayList<>();
			List<UncertainArbDiscFunc> epochCmlBounds = new ArrayList<>();
			List<Double> epochWeights = new ArrayList<>();
			
			for (int e=0; e<epochs.length; e++) {
				SeismicityRateModel model;
				
				switch (reg) {
				case CAR_INTERFACE:
					model = PRVI25_SubductionCaribbeanSeismicityRate.loadRateModel(epochs[e], false);
					break;
				case CAR_INTRASLAB:
					model = PRVI25_SubductionCaribbeanSeismicityRate.loadRateModel(epochs[e], true);
					break;
				case CRUSTAL:
					model = PRVI25_CrustalSeismicityRate.loadRateModel(epochs[e]);
					break;
				case MUE_INTERFACE:
					model = PRVI25_SubductionMuertosSeismicityRate.loadRateModel(epochs[e], false);
					break;
				case MUE_INTRASLAB:
					model = PRVI25_SubductionMuertosSeismicityRate.loadRateModel(epochs[e], true);
					break;

				default:
					throw new IllegalStateException("Unknown region: "+reg);
				}
				
				UncertainBoundedIncrMagFreqDist bounded = model.getBounded(refMFD, mfdMmax);
				
				IncrementalMagFreqDist lower = bounded.getLower();
				lower.setName(null);
				incrFuncs.add(lower);
				incrChars.add(new PlotCurveCharacterstics(epochLines[e], 2f, color));
				EvenlyDiscretizedFunc cmlLower = lower.getCumRateDistWithOffset();
				cmlFuncs.add(cmlLower);
				cmlChars.add(incrChars.get(incrChars.size()-1));
				
				IncrementalMagFreqDist upper = bounded.getUpper();
				upper.setName(null);
				incrFuncs.add(upper);
				incrChars.add(new PlotCurveCharacterstics(epochLines[e], 2f, color));
				EvenlyDiscretizedFunc cmlUpper = upper.getCumRateDistWithOffset();
				cmlFuncs.add(cmlUpper);
				cmlChars.add(incrChars.get(incrChars.size()-1));
				
				bounded.setName(epochs[e].getShortName()+": "+bounded.getName());;
				incrFuncs.add(bounded);
				incrChars.add(new PlotCurveCharacterstics(epochLines[e], 5f, color));
				EvenlyDiscretizedFunc cmlPref = bounded.getCumRateDistWithOffset();
				cmlFuncs.add(cmlPref);
				cmlChars.add(incrChars.get(incrChars.size()-1));
				
				epochBounds.add(bounded);
				epochCmlBounds.add(new UncertainArbDiscFunc(cmlPref, cmlLower, cmlUpper, bounded.getBoundType(), null));
				epochWeights.add(epochs[e].getNodeWeight(null));
			}
			
			UncertainBoundedIncrMagFreqDist averageBounded = PRVI25_SeismicityRateEpoch.averageUncert(epochBounds, epochWeights);
			UncertainArbDiscFunc averageCmlBounded = PRVI25_SeismicityRateEpoch.averageUncertCml(epochCmlBounds, epochWeights);
			
//			Color avgColor = Color.GRAY;
			Color avgColor = Colors.tab_purple;
			UncertainBoundedIncrMagFreqDist boundsCopy = averageBounded.deepClone();
			boundsCopy.setName("Mixture "+boundsCopy.getBoundName());
			incrFuncs.add(0, boundsCopy);
			incrChars.add(0, new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(0, 0, 0, 60)));
			averageBounded.setName("Weighted average");
			incrFuncs.add(averageBounded);
			incrChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, avgColor));
			
			UncertainArbDiscFunc boundsCmlCopy = averageCmlBounded.deepClone();
			boundsCmlCopy.setName("Mixture "+boundsCmlCopy.getBoundName());
			cmlFuncs.add(0, boundsCmlCopy);
			cmlChars.add(0, new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(0, 0, 0, 60)));
			averageCmlBounded.setName("Weighted average");
			cmlFuncs.add(averageCmlBounded);
			cmlChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, avgColor));
			
			for (boolean is1973 : new boolean[] {false,true}) {
				File ratesFile = new File(directDir, PRVI25_CrustalSeismicityRate.getDirectRateFileName(reg,
						is1973 ? PRVI25_SeismicityRateEpoch.RECENT : PRVI25_SeismicityRateEpoch.FULL));
				
				CSVFile<String> totalRateCSV = CSVFile.readFile(ratesFile, false);
				List<Direct> directs = SeismicityRateFileLoader.loadDirectBranches(totalRateCSV);
				
				double minFuncMag = is1973 ? 5.01 : 6.01;
				double maxFuncMag = directs.get(0).maxObsIncrMag - 0.01;
				if (maxFuncMag <= 0d) {
					// not found
//					maxFuncMag = directs.get(0).cumulativeDist.getMaxX()+0.01;
					maxFuncMag = directs.get(0).cumulativeDist.getMaxX()-0.01;
				}
				
				Direct meanDirect = SeismicityRateFileLoader.locateMean(directs);
				if (is1973)
					sumNobs1973 += meanDirect.nObs;
				else
					sumNobs1900 += meanDirect.nObs;
				
				System.out.println("Direct "+regName+", "+(is1973 ? "1973" : "1900")+" Nobs="+(float)meanDirect.nObs);

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
//				meanIncrObs.setName(obsName);
//				incrFuncs.add(meanIncrObs);
//				incrChars.add(obsChar);
				meanObsCml.setName(meanIncrObs.getName());
				cmlFuncs.add(meanObsCml);
				cmlChars.add(obsChar);
			}
			
			Preconditions.checkState(incrFuncs.size() == incrChars.size());
			Preconditions.checkState(cmlFuncs.size() == cmlChars.size());
			
			PlotSpec incrPlot = new PlotSpec(incrFuncs, incrChars, name, "Magnitude", "Incremental Rate (1/yr)");
			incrPlot.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
//			incrPlot.setLegendInset(RectangleAnchor.TOP_RIGHT);
			
			gp.drawGraphPanel(incrPlot, false, true, magRange, yRange);
			
			PlotUtils.writePlots(outputDir, reg.name(), gp, 700, 600, true, true, false);
			
			PlotSpec cmlPlot = new PlotSpec(cmlFuncs, cmlChars, name, "Magnitude", "Cumulative Rate (1/yr)");
			cmlPlot.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
//			cmlPlot.setLegendInset(RectangleAnchor.TOP_RIGHT);
//			cmlPlot.setLegendInset(RectangleAnchor.TOP, 0.5, 0.95, 0.9, false);
			
			gp.drawGraphPanel(cmlPlot, false, true, magRange, yRange);
			
			PlotUtils.writePlots(outputDir, reg.name()+"_cml", gp, 700, 600, true, true, false);
		}
		
		System.out.println("Direct total 1900 Nobs="+(float)sumNobs1900);
		System.out.println("Direct total 1973 Nobs="+(float)sumNobs1973);
	}

}
