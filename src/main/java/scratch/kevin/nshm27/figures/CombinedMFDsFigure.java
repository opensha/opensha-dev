package scratch.kevin.nshm27.figures;

import static scratch.kevin.nshm27.figures.NSHM27_PaperPaths.*;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchRegionalMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchRegionalMFDs.MFDType;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.RegionsOfInterest;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_SeisClassificationMethod;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_SeisRateModelBranch;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_SeisRateModelSamples;
import gov.usgs.earthquake.nshmp.erf.nshm27.util.NSHM27_RegionLoader;
import gov.usgs.earthquake.nshmp.erf.nshm27.util.NSHM27_RegionLoader.NSHM27_SeismicityRegions;
import gov.usgs.earthquake.nshmp.erf.seismicity.SeismicityRateFileLoader.PureGR;
import net.mahdilamb.colormap.Colors;
import oracle.net.aso.i;

public class CombinedMFDsFigure {

	public static void main(String[] args) throws IOException {
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(4.99, 11.99);
//		Range xRange = new Range(5d, 9.5d);
		Range xRange = new Range(6d, 9.5d);
		Range yRange = new Range (1e-6, 1e1);
		double fullWidth = PlotUtils.DEFAULT_USABLE_PAGE_WIDTH;
		double fullHeight = fullWidth*0.75;
		double halfWidth = PlotUtils.DEFAULT_USABLE_PAGE_WIDTH/2d;
		double halfHeight = halfWidth*0.75;
		
		File outputDir = new File(FIGURES_DIR, "combined_mfds");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		Color fractileBase = Color.DARK_GRAY;
		PlotCurveCharacterstics extremaChar = new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f,
				new Color(fractileBase.getRed(), fractileBase.getGreen(), fractileBase.getBlue(), 50));
		PlotCurveCharacterstics bounds95Char = new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f,
				new Color(fractileBase.getRed(), fractileBase.getGreen(), fractileBase.getBlue(), 70));
		PlotCurveCharacterstics bounds68Char = new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f,
				new Color(fractileBase.getRed(), fractileBase.getGreen(), fractileBase.getBlue(), 100));
		double[] fractiles = {0d, 0.025, 0.16, 0.5, 0.84, 0.975, 1d};
		String fractileLabel = "p[0, 2.5, 16, 84, 97.5, 100]";
		
//		double[] csvMags = { 6d, 7d, 7.5d, 8d, 8.5d, 9d };
		double[] csvMags = { 6d, 7d, 8d, 9d };
		List<String> csvRateHeader = new ArrayList<>();
		List<String> csvRIHeader = new ArrayList<>();
		csvRateHeader.add("Model Component");
		csvRIHeader.add("Model Component");
		for (double mag : csvMags) {
			csvRateHeader.add("M>"+oDF.format(mag)+" rate");
			csvRIHeader.add("M>"+oDF.format(mag)+" RI");
		}
		
		for (NSHM27_SeismicityRegions seisReg : NSHM27_SeismicityRegions.values()) {
			FaultSystemSolution sol = getSolution(seisReg);
			GridSourceList gridList = sol.requireModule(GridSourceList.class);
			Region reg = seisReg.load();
			
			for (boolean interfaceOnly : new boolean[] {false, true}) {
				CSVFile<String> rateCSV = null;
				CSVFile<String> riCSV = null;
				if (!interfaceOnly) {
					rateCSV = new CSVFile<>(true);
					rateCSV.addLine(csvRateHeader);
					riCSV = new CSVFile<>(true);
					riCSV.addLine(csvRIHeader);
				}
				TectonicRegionType targetTRT = interfaceOnly ? TectonicRegionType.SUBDUCTION_INTERFACE : null;
				List<IncrementalMagFreqDist> indvIncrFuncs = new ArrayList<>();
				List<PlotCurveCharacterstics> indvChars = new ArrayList<>();
				
				SummedMagFreqDist gridSum = new SummedMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
				SummedMagFreqDist faultSum = new SummedMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
				
				for (TectonicRegionType trt : TRTs) {
					if (interfaceOnly && trt != TectonicRegionType.SUBDUCTION_INTERFACE)
						continue;
					Color gridColor = getColor(trt, IncludeBackgroundOption.ONLY);
					Color faultColor = getColor(trt, IncludeBackgroundOption.EXCLUDE);
					Color sumColor = getColor(trt, IncludeBackgroundOption.INCLUDE);
					String trtName = NSHM27_RegionLoader.getNameForTRT(trt);
					
					IncrementalMagFreqDist gridMFD = loadGridded(gridList, trt, refMFD);
					IncrementalMagFreqDist faultMFD = loadFault(sol, trt, refMFD, reg);
					
					processCSVs(faultMFD, trtName+" (on-fault)", csvMags, rateCSV, riCSV);
					processCSVs(gridMFD, trtName+" (gridded)", csvMags, rateCSV, riCSV);
					
					if (faultMFD != null) {
						if (!interfaceOnly) {
							SummedMagFreqDist sum = new SummedMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
							sum.addIncrementalMagFreqDist(gridMFD);
							sum.addIncrementalMagFreqDist(faultMFD);
							
							sum.setName(trtName+" (combined)");
							indvIncrFuncs.add(sum);
							indvChars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 2f, sumColor));
							
							processCSVs(sum, trtName+" (combined)", csvMags, rateCSV, riCSV);
						}
						
						faultMFD.setName(trtName+" (on-fault)");
						indvIncrFuncs.add(faultMFD);
						indvChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, faultColor));
						faultSum.addIncrementalMagFreqDist(faultMFD);
					}
					
					gridMFD.setName(trtName+" (gridded)");
					indvIncrFuncs.add(gridMFD);
					indvChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, gridColor));
					gridSum.addIncrementalMagFreqDist(gridMFD);
				}
				
				List<IncrementalMagFreqDist> incrFuncs = new ArrayList<>();
				List<DiscretizedFunc> cmlFuncs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				SummedMagFreqDist mfdSum = new SummedMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
				mfdSum.addIncrementalMagFreqDist(gridSum);
				mfdSum.addIncrementalMagFreqDist(faultSum);
				if (interfaceOnly)
					mfdSum.setName("Total interface");
				else
					mfdSum.setName("Total");
				incrFuncs.add(mfdSum);
				cmlFuncs.add(mfdSum.getCumRateDistWithOffset());
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
				
				processCSVs(mfdSum, "Total", csvMags, rateCSV, riCSV);
				
				// uncertainty for legend
				IncrementalMagFreqDist emptyFunc = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
				UncertainBoundedIncrMagFreqDist emptyUnceratain = new UncertainBoundedIncrMagFreqDist(emptyFunc, emptyFunc, emptyFunc, null);
				emptyUnceratain.setName(fractileLabel);
				incrFuncs.add(emptyUnceratain);
				cmlFuncs.add(emptyUnceratain);
				chars.add(bounds68Char);
				
				// individual components
				for (int i=0; i<indvIncrFuncs.size(); i++) {
					IncrementalMagFreqDist mfd = indvIncrFuncs.get(i);
					incrFuncs.add(mfd);
					cmlFuncs.add(mfd.getCumRateDistWithOffset());
					chars.add(indvChars.get(i));
				}
				
				List<PureGR> samples = new ArrayList<>();
				samples.addAll(NSHM27_SeisRateModelSamples.loadOrigSamples(seisReg, NSHM27_SeisClassificationMethod.PROFACE, targetTRT));
				int numProface = samples.size();
				samples.addAll(NSHM27_SeisRateModelSamples.loadOrigSamples(seisReg, NSHM27_SeisClassificationMethod.PROSLAB, targetTRT));
				int numProslab = samples.size() - numProface;
				double weightProface = NSHM27_SeisClassificationMethod.PROFACE.getNodeWeight();
				double weightProslab = NSHM27_SeisClassificationMethod.PROSLAB.getNodeWeight();
				if ((float)(weightProface+weightProslab) != 1f) {
					double sum = weightProface + weightProslab;
					weightProface /= sum;
					weightProslab /= sum;
				}
				double[] weights = new double[samples.size()];
				double weightProfaceEach = weightProface/(double)numProface;
				double weightProslabEach = weightProslab/(double)numProslab;
				for (int i=0; i<numProface; i++)
					weights[i] = weightProfaceEach;
				for (int i=numProface; i<samples.size(); i++)
					weights[i] = weightProslabEach;
				double weightSum = StatUtils.sum(weights);
				Preconditions.checkState((float)weightSum == 1f, "Bad weightSum=%s for proface=%s (%s), proslab=%s (%s)",
						(float)weightSum, weightProfaceEach, numProface, weightProslabEach, numProslab);
				double[][] incrVals = new double[refMFD.size()][samples.size()];
				double[][] cmlVals = new double[refMFD.size()][samples.size()];
				double obsB = 0d;
				IncrementalMagFreqDist obsIncrMean = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
				EvenlyDiscretizedFunc obsCmlMean = obsIncrMean.getCumRateDistWithOffset();
				double obsProfaceB = 0d;
				IncrementalMagFreqDist obsProfaceIncrMean = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
				EvenlyDiscretizedFunc obsProfaceCmlMean = obsIncrMean.getCumRateDistWithOffset();
				double obsProslabB = 0d;
				IncrementalMagFreqDist obsProslabIncrMean = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
				EvenlyDiscretizedFunc obsProslabCmlMean = obsIncrMean.getCumRateDistWithOffset();
				for (int s=0; s<samples.size(); s++) {
					PureGR sample = samples.get(s);
					GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
					// this sets shape, min/max
					// subtract a tiny amount from mMax so that if it's exactly at a bin edge, e.g. 7.9, it rounds down, e.g. to 7.85
					gr.setAllButTotCumRate(refMFD.getMinX(), refMFD.getMaxX(), 1e16, sample.b);
					// this scales it to match
					// similarly, add a tiny amount to M1 so that if it's exactly at a bin edge (which it should be as it's determined
					// using cumulative binning), it rounds up to the incremental bin for that cumulative edge
					gr.scaleToCumRate(refMFD.getClosestXIndex(sample.M1+0.001), sample.rateAboveM1);
					
					EvenlyDiscretizedFunc cml = gr.getCumRateDistWithOffset();
					
					obsB += sample.b * weights[s];
					if (s < numProface)
						obsProfaceB += sample.b * weights[s];
					else
						obsProslabB += sample.b * weights[s];
					
					for (int m=0; m<refMFD.size(); m++) {
						double incrRate = gr.getY(m);
						double cmlRate = cml.getY(m);
						obsIncrMean.set(m, Math.fma(weights[s], incrRate, obsIncrMean.getY(m)));
						obsCmlMean.set(m, Math.fma(weights[s], cmlRate, obsCmlMean.getY(m)));
						if (s < numProface) {
							obsProfaceIncrMean.set(m, Math.fma(weights[s], incrRate, obsProfaceIncrMean.getY(m)));
							obsProfaceCmlMean.set(m, Math.fma(weights[s], cmlRate, obsProfaceCmlMean.getY(m)));							
						} else {
							obsProslabIncrMean.set(m, Math.fma(weights[s], incrRate, obsProslabIncrMean.getY(m)));
							obsProslabCmlMean.set(m, Math.fma(weights[s], cmlRate, obsProslabCmlMean.getY(m)));
						}
						incrVals[m][s] = incrRate;
						cmlVals[m][s] = cmlRate;
					}
				}
				obsProfaceB /= weightProface;
				obsProfaceIncrMean.scale(1d/weightProface);
				obsProfaceCmlMean.scale(1d/weightProface);
				obsProslabB /= weightProslab;
				obsProslabIncrMean.scale(1d/weightProslab);
				obsProslabCmlMean.scale(1d/weightProslab);
				IncrementalMagFreqDist obsIncrLow = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
				IncrementalMagFreqDist obsIncrHigh = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
				EvenlyDiscretizedFunc obsCmlLow = obsIncrLow.getCumRateDistWithOffset();
				EvenlyDiscretizedFunc obsCmlHigh = obsIncrHigh.getCumRateDistWithOffset();
				for (int m=0; m<refMFD.size(); m++) {
					LightFixedXFunc incrPDF = ArbDiscrEmpiricalDistFunc.calcQuickNormCDF(incrVals[m], weights);
					LightFixedXFunc cmlPDF = ArbDiscrEmpiricalDistFunc.calcQuickNormCDF(cmlVals[m], weights);
					obsIncrLow.set(m, ArbDiscrEmpiricalDistFunc.calcFractileFromNormCDF(incrPDF, 0.025));
					obsIncrHigh.set(m, ArbDiscrEmpiricalDistFunc.calcFractileFromNormCDF(incrPDF, 0.975));
					obsCmlLow.set(m, ArbDiscrEmpiricalDistFunc.calcFractileFromNormCDF(cmlPDF, 0.025));
					obsCmlHigh.set(m, ArbDiscrEmpiricalDistFunc.calcFractileFromNormCDF(cmlPDF, 0.975));
				}
				
				// observed seismicity
				
				// add a b=1 line
				GutenbergRichterMagFreqDist gr1incr = new GutenbergRichterMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
				GutenbergRichterMagFreqDist gr1cml = new GutenbergRichterMagFreqDist(refMFD.getMinX()-0.5*refMFD.getDelta(), refMFD.size(), refMFD.getDelta());
				gr1incr.setAllButTotCumRate(gr1incr.getMinX(), gr1incr.getMaxX(), 1e16, 1d);
				gr1cml.setAllButTotCumRate(gr1cml.getMinX(), gr1cml.getMaxX(), 1e16, 1d);
				int scaleIndex = refMFD.getClosestXIndex(xRange.getLowerBound()+0.01);
				gr1incr.scaleToIncrRate(scaleIndex, mfdSum.getY(scaleIndex));
				gr1cml.scaleToIncrRate(scaleIndex, mfdSum.getCumRate(scaleIndex));
				gr1incr.setName("GR b=1");
				gr1cml.setName(gr1incr.getName());
//				IncrementalMagFreqDist grClone = gr1incr.deepClone();
//				grClone.setName(null);
//				incrFuncs.add(grClone);
//				grClone = gr1cml.deepClone();
//				grClone.setName(null);
//				cmlFuncs.add(grClone);
//				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.WHITE));
//				incrFuncs.add(gr1incr);
//				cmlFuncs.add(gr1cml);
//				chars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 1f, Color.DARK_GRAY));
				incrFuncs.add(gr1incr);
				cmlFuncs.add(gr1cml);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 1f, Color.DARK_GRAY));
				
//				double obsB = ((PureGR)NSHM27_SeisRateModelBranch.PREFFERRED.getRateRecord(seisReg, NSHM27_SeisClassificationMethod.PROFACE, null)).b;
//				IncrementalMagFreqDist obsIncrMean = NSHM27_SeisRateModelBranch.PREFFERRED.build(
//						seisReg, NSHM27_SeisClassificationMethod.PROFACE, null, refMFD, refMFD.getMaxX());
//				EvenlyDiscretizedFunc obsCmlMean = obsIncrMean.getCumRateDistWithOffset();
				
				obsIncrLow.setName(null);
				incrFuncs.add(obsIncrLow);
				obsCmlLow.setName(null);
				cmlFuncs.add(obsCmlLow);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, OBS_RATE_COLOR));
				
				obsIncrHigh.setName(null);
				incrFuncs.add(obsIncrHigh);
				obsCmlHigh.setName(null);
				cmlFuncs.add(obsCmlHigh);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, OBS_RATE_COLOR));
				
//				obsIncrMean.setName(obsIncrMean.getName()+", b="+(float)obsB);
				obsIncrMean.setName("Observed mean & 95% (b="+twoDF.format(obsB)+")");
				incrFuncs.add(obsIncrMean);
				obsCmlMean.setName(obsIncrMean.getName());
				cmlFuncs.add(obsCmlMean);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, OBS_RATE_COLOR));
				
				if (interfaceOnly) {
					// proface and proslab means separately
//					IncrementalMagFreqDist obsProfaceIncrMean = NSHM27_SeisRateModelBranch.PREFFERRED.build(
//							seisReg, NSHM27_SeisClassificationMethod.PROFACE, targetTRT, refMFD, refMFD.getMaxX());
//					IncrementalMagFreqDist obsProslabIncrMean = NSHM27_SeisRateModelBranch.PREFFERRED.build(
//							seisReg, NSHM27_SeisClassificationMethod.PROSLAB, targetTRT, refMFD, refMFD.getMaxX());
					
					obsProfaceIncrMean.setName("Observed pro-face mean");
					incrFuncs.add(obsProfaceIncrMean);
					obsProfaceCmlMean.setName(obsProfaceIncrMean.getName());
					cmlFuncs.add(obsProfaceCmlMean);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1.5f, Colors.tab_lightpurple));
					
					obsProslabIncrMean.setName("Observed pro-slab mean");
					incrFuncs.add(obsProslabIncrMean);
					obsProslabCmlMean.setName(obsProslabIncrMean.getName());
					cmlFuncs.add(obsProslabCmlMean);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1.5f, Colors.tab_lightaqua));
				}
				
				// locate uncertainties
				BranchRegionalMFDs regMFDs = sol.requireModule(BranchRegionalMFDs.class);
				RegionsOfInterest roi = sol.getRupSet().requireModule(RegionsOfInterest.class);
				Preconditions.checkArgument(roi.getRegions().size() == regMFDs.getNumRegions());
				int regionIndex = -1;
				for (int r=0; r<regMFDs.getNumRegions(); r++) {
					Region tmpReg = roi.getRegions().get(r);
					if (roi.getTRTs().get(r) == targetTRT && tmpReg.equalsRegion(reg)) {
						regionIndex = r;
						break;
					}
				}
				Preconditions.checkState(regionIndex >= 0);
				IncrementalMagFreqDist[] incrFractiles = regMFDs.calcRegionalIncrementalFractiles(MFDType.SUM, regionIndex, fractiles);
				UncertainBoundedIncrMagFreqDist incrExtrema = new UncertainBoundedIncrMagFreqDist(incrFractiles[3], incrFractiles[0], incrFractiles[6], null);
				UncertainBoundedIncrMagFreqDist incr95 = new UncertainBoundedIncrMagFreqDist(incrFractiles[3], incrFractiles[1], incrFractiles[5], null);
				UncertainBoundedIncrMagFreqDist incr68 = new UncertainBoundedIncrMagFreqDist(incrFractiles[3], incrFractiles[2], incrFractiles[4], null);
				EvenlyDiscretizedFunc[] cmlFractiles = regMFDs.calcRegionalCumulativeFractiles(MFDType.SUM, regionIndex, fractiles);
				UncertainArbDiscFunc cmlExtrema = new UncertainArbDiscFunc(cmlFractiles[3], cmlFractiles[0], cmlFractiles[6]);
				UncertainArbDiscFunc cml95 = new UncertainArbDiscFunc(cmlFractiles[3], cmlFractiles[1], cmlFractiles[5]);
				UncertainArbDiscFunc cml68 = new UncertainArbDiscFunc(cmlFractiles[3], cmlFractiles[2], cmlFractiles[4]);
				
				processCSVsWithCml(cmlFractiles[0], "Total (minimum)", csvMags, rateCSV, riCSV);
				processCSVsWithCml(cmlFractiles[1], "Total (2.5 %-ile)", csvMags, rateCSV, riCSV);
				processCSVsWithCml(cmlFractiles[2], "Total (16 %-ile)", csvMags, rateCSV, riCSV);
				processCSVsWithCml(cmlFractiles[5], "Total (50 %-ile)", csvMags, rateCSV, riCSV);
				processCSVsWithCml(cmlFractiles[4], "Total (84 %-ile)", csvMags, rateCSV, riCSV);
				processCSVsWithCml(cmlFractiles[5], "Total (97.5 %-ile)", csvMags, rateCSV, riCSV);
				processCSVsWithCml(cmlFractiles[6], "Total (maximum)", csvMags, rateCSV, riCSV);
				
				processCSVsWithCml(obsCmlMean, "Observed", csvMags, rateCSV, riCSV);
				processCSVsWithCml(obsCmlLow, "Observed (2.5 %-ile)", csvMags, rateCSV, riCSV);
				processCSVsWithCml(obsCmlHigh, "Observed (97.5 %-ile)", csvMags, rateCSV, riCSV);

				incrExtrema.setName(null);
				incrFuncs.add(incrExtrema);
				cmlExtrema.setName(null);
				cmlFuncs.add(cmlExtrema);
				chars.add(extremaChar);

				incr95.setName(null);
				incrFuncs.add(incr95);
				cml95.setName(null);
				cmlFuncs.add(cml95);
				chars.add(bounds95Char);

				incr68.setName(null);
				incrFuncs.add(incr68);
				cml68.setName(null);
				cmlFuncs.add(cml68);
				chars.add(bounds68Char);
				
				RectangleAnchor anchor = RectangleAnchor.TOP_RIGHT;
//				RectangleAnchor anchor = RectangleAnchor.BOTTOM_LEFT;
				
				String title = seisReg.getTitleCaseAcronym();
				if (interfaceOnly)
					title += " (interface)";
//				else
//					title += " overall";
				
				PlotSpec incrPlot = new PlotSpec(incrFuncs, chars, title, "Magnitude", "Incremental rate (1/yr)");
				incrPlot.setLegendInset(anchor);
				
				PlotSpec cmlPlot = new PlotSpec(cmlFuncs, chars, title, "Magnitude", "Cumulative rate (1/yr)");
				cmlPlot.setLegendInset(anchor);
				
				HeadlessGraphPanel gp = PlotUtils.initPrintHeadless();
				
				gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
				
				String prefix = seisReg.name();
				if (interfaceOnly)
					prefix += "_interface";
				
				gp.drawGraphPanel(incrPlot, false, true, xRange, yRange);
				
//				PlotUtils.writePrintPlots(outputDir, prefix+"_incr_half", gp, halfWidth, halfHeight, 150, true, true, false);
				PlotUtils.writePrintPlots(outputDir, prefix+"_incr", gp, fullWidth, fullHeight, 150, true, true, false);
				
				gp.drawGraphPanel(cmlPlot, false, true, xRange, yRange);
				
//				PlotUtils.writePrintPlots(outputDir, prefix+"_cml_half", gp, halfWidth, halfHeight, 150, true, true, false);
				PlotUtils.writePrintPlots(outputDir, prefix+"_cml", gp, fullWidth, fullHeight, 150, true, true, false);
				
				if (!interfaceOnly) {
					rateCSV.writeToFile(new File(outputDir, prefix+"_cml_rates.csv"));
					riCSV.writeToFile(new File(outputDir, prefix+"_cml_ris.csv"));
				}
			}
		}
	}
	
	private static IncrementalMagFreqDist loadGridded(GridSourceList gridList,
			TectonicRegionType trt, EvenlyDiscretizedFunc refMFD) {
		SummedMagFreqDist mfdSum = new SummedMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
		for (int l=0; l<gridList.getNumLocations(); l++) {
			IncrementalMagFreqDist mfd = gridList.getMFD(trt, l);
			if (mfd != null)
				mfdSum.addIncrementalMagFreqDist(mfd);
		}
		return mfdSum;
	}
	
	private static IncrementalMagFreqDist loadFault(FaultSystemSolution sol,
			TectonicRegionType trt, EvenlyDiscretizedFunc refMFD, Region reg) {
		SummedMagFreqDist mfdSum = new SummedMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
		double[] fractsInReg = sol.getRupSet().getFractSectsInsideRegion(reg, false);
		
		for (int s=0; s<sol.getRupSet().getNumSections(); s++) {
			FaultSection sect = sol.getRupSet().getFaultSectionData(s);
			if (sect.getTectonicRegionType() == trt) {
				IncrementalMagFreqDist mfd = sol.calcNucleationMFD_forSect(
						s, refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
				mfdSum.addIncrementalMagFreqDist(mfd, fractsInReg[s]);
			}
		}
		if (mfdSum.calcSumOfY_Vals() > 0d)
			return mfdSum;
		return null;
	}
	
	private static void processCSVs(IncrementalMagFreqDist mfd, String name, double[] mags,
			CSVFile<String> rateCSV, CSVFile<String> riCSV) {
		if (mfd == null || rateCSV == null)
			return;
		EvenlyDiscretizedFunc cml = mfd.getCumRateDistWithOffset();
		processCSVsWithCml(cml, name, mags, rateCSV, riCSV);
	}
	
	private static void processCSVsWithCml(EvenlyDiscretizedFunc cmlMFD, String name, double[] mags,
			CSVFile<String> rateCSV, CSVFile<String> riCSV) {
		if (cmlMFD == null || rateCSV == null)
			return;
		List<String> rateLine = new ArrayList<>();
		List<String> riLine = new ArrayList<>();
		rateLine.add(name);
		riLine.add(name);
		for (double mag : mags) {
			int index = cmlMFD.getClosestXIndex(mag);
			double rate = cmlMFD.getY(index);
			rateLine.add((float)rate+"");
			double ri = 1d/rate;
			if (ri > 10d)
				riLine.add((int)Math.round(ri)+"");
			else if (ri > 1d)
				riLine.add(oDF.format(ri));
			else if (ri > 0.1)
				riLine.add(twoDF.format(ri));
			else
				riLine.add((float)ri+"");
		}
		rateCSV.addLine(rateLine);
		riCSV.addLine(riLine);
	}
}
