package scratch.kevin.prvi25.figures;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.apache.commons.math3.util.Precision;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateFileLoader;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateFileLoader.Exact;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateFileLoader.PureGR;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateFileLoader.RateRecord;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateFileLoader.RateType;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SeismicityRateEpoch;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionCaribbeanSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionMuertosSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader.PRVI25_SeismicityRegions;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;

import com.google.common.base.Preconditions;

import net.mahdilamb.colormap.Colors;
import scratch.kevin.prvi25.GriddedRateDistributionSolutionWriter;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

public class ObsUncertaintyBoundsFigure {

	public static void main(String[] args) throws IOException {
		File outputDir = new File(FIGURES_DIR, "obs_mfd_bounds");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		boolean incremental = false;
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(3.05, 10.95);
		
		Color[] colors = {
				Colors.tab_blue,
				Colors.tab_orange,
				Colors.tab_green
		};
		
		PlotLineType[] lineTypes = {
				PlotLineType.DASHED,
				PlotLineType.SOLID,
				PlotLineType.DOTTED
		};
		
		RateType[] types = {
			RateType.M1,
			RateType.M1_TO_MMAX,
			RateType.EXACT
		};
		
		List<PRVI25_SeismicityRegions> seisRegions = new ArrayList<>(List.of(PRVI25_SeismicityRegions.values()));
		seisRegions.add(null); // reference flag
		PRVI25_SeismicityRegions refRef = PRVI25_SeismicityRegions.CRUSTAL;
		PRVI25_SeismicityRateEpoch epoch = PRVI25_SeismicityRateEpoch.DEFAULT;
//		File refRatePairsFile = null;
//		File refRatePairsFile = new File("/home/kevin/OpenSHA/nshm23/prvi/rate_raw_data/2025_03_26/rbpairs-Crustal-Full-v3.csv");
		File refRatePairsFile = new File("/home/kevin/OpenSHA/nshm23/prvi/rate_raw_data/2025_07_17/1900_2023/rbpairs-Crustal-Prob-v9.csv");
		
		double weightLow = PRVI25_CrustalSeismicityRate.LOW.getNodeWeight(null);
		double weightPref = PRVI25_CrustalSeismicityRate.PREFFERRED.getNodeWeight(null);
		double weightHigh = PRVI25_CrustalSeismicityRate.HIGH.getNodeWeight(null);
		
		for (PRVI25_SeismicityRegions seisReg : seisRegions) {
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			boolean ref = seisReg == null;
			String prefix, title;
			if (ref) {
				seisReg = refRef;
				prefix = "reference";
				title = " ";
			} else {
				prefix = seisReg.name();
				switch (seisReg) {
				case CRUSTAL:
					title = "Active crustal";
					break;
				case CAR_INTERFACE:
					title = "Caribbean interface";
					break;
				case CAR_INTRASLAB:
					title = "Caribbean intraslab";
					break;
				case MUE_INTERFACE:
					title = "Muertos interface";
					break;
				case MUE_INTRASLAB:
					title = "Muertos intraslab";
					break;

				default:
					throw new IllegalStateException();
				}
			}
			
			DecimalFormat oDF = new DecimalFormat("0.#");

			Double m1 = null;
			Double mMax = null;
			Map<RateType, EvenlyDiscretizedFunc[]> typeMagFuncs = new HashMap<>();
			EvenlyDiscretizedFunc overallMean = null;
			for (RateType type : types) {
				List<? extends RateRecord> rates;
				
				switch (seisReg) {
				case CRUSTAL:
					rates = PRVI25_CrustalSeismicityRate.loadRates(epoch, type);
					break;
				case CAR_INTERFACE:
					rates = PRVI25_SubductionCaribbeanSeismicityRate.loadRates(epoch, type, false);
					break;
				case CAR_INTRASLAB:
					rates = PRVI25_SubductionCaribbeanSeismicityRate.loadRates(epoch, type, true);
					break;
				case MUE_INTERFACE:
					rates = PRVI25_SubductionMuertosSeismicityRate.loadRates(epoch, type, false);
					break;
				case MUE_INTRASLAB:
					rates = PRVI25_SubductionMuertosSeismicityRate.loadRates(epoch, type, true);
					break;

				default:
					throw new IllegalStateException();
				}
				
				RateRecord meanRec = SeismicityRateFileLoader.locateMean(rates);
				if (m1 == null)
					m1 = meanRec.M1;
				if (mMax == null && meanRec instanceof PureGR)
					mMax = ((PureGR)meanRec).Mmax;
				
				EvenlyDiscretizedFunc meanMFD;
				if (incremental)
					meanMFD = SeismicityRateFileLoader.buildIncrementalMFD(meanRec, refMFD, refMFD.getMaxX());
				else
					meanMFD = cmlMFD(meanRec, refMFD);
				
				if (funcs.isEmpty()) {
					meanMFD.setName("Mean");
					funcs.add(meanMFD);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
					overallMean = meanMFD;
				}
				
				Color color = colors[type.ordinal()];
				PlotLineType plt = lineTypes[type.ordinal()];
				RateRecord low = SeismicityRateFileLoader.locateQuantile(rates, 0.025);
				RateRecord high = SeismicityRateFileLoader.locateQuantile(rates, 0.975);

				EvenlyDiscretizedFunc lowMFD;
				EvenlyDiscretizedFunc highMFD;
				if (incremental) {
					lowMFD = SeismicityRateFileLoader.buildIncrementalMFD(low, refMFD, refMFD.getMaxX());
					highMFD = SeismicityRateFileLoader.buildIncrementalMFD(high, refMFD, refMFD.getMaxX());
				} else {
					lowMFD = cmlMFD(low, refMFD);
					highMFD = cmlMFD(high, refMFD);
				}
				
				typeMagFuncs.put(type, new EvenlyDiscretizedFunc[] {lowMFD, highMFD});
				
				lowMFD.setName(type.toString().replace("Branches", "branches"));
				funcs.add(lowMFD);
				chars.add(new PlotCurveCharacterstics(plt, 3f, color));
				highMFD.setName(null);
				funcs.add(highMFD);
				chars.add(new PlotCurveCharacterstics(plt, 3f, color));
				
				if (!ref && type == RateType.M1_TO_MMAX) {
					EvenlyDiscretizedFunc weightAvg = new EvenlyDiscretizedFunc(meanMFD.getMinX(), meanMFD.size(), meanMFD.getDelta());
					Preconditions.checkState((float)meanMFD.getMinX() == (float)lowMFD.getMinX());
					Preconditions.checkState((float)meanMFD.getMinX() == (float)highMFD.getMinX());
					for (int i=0; i<weightAvg.size(); i++)
						weightAvg.set(i, weightLow*lowMFD.getY(i) + weightPref*meanMFD.getY(i) + weightHigh*highMFD.getY(i));
					
					if (!incremental) {
						System.out.println(title);
						System.out.println("\tM>5: "+(float)weightAvg.getY(weightAvg.getClosestXIndex(5.01)));
						System.out.println("\tM>6: "+(float)weightAvg.getY(weightAvg.getClosestXIndex(6.01)));
						System.out.println("\tM>6 snapped: "+(float)weightAvg.getX(weightAvg.getClosestXIndex(6.01)));
						System.out.println("\tM1="+m1.floatValue());
						System.out.println("\tMmax="+mMax.floatValue());
					}
					
					weightAvg.setName(type.toString()+" Average");
					funcs.add(weightAvg);
					chars.add(new PlotCurveCharacterstics(plt, 3f, Color.DARK_GRAY));
				}
			}
			
			for (XY_DataSet func : funcs) {
				if (func.getName() != null && func.getName().contains("M1"))
					func.setName(func.getName().replace("M1", "M₁"));
				if (func.getName() != null && func.getName().contains("Mmax"))
					func.setName(func.getName().replace("Mmax", "Mₘₐₓ"));
			}
			
			Range xRange = new Range(4d, 8d);
			Range yRange = incremental ? new Range(1e-5, 1e1) : new Range(1e-4, 1e2);
			
			List<XYTextAnnotation> anns = new ArrayList<>();
			Font annFont = new Font(Font.SANS_SERIF, Font.PLAIN, 22);
			
			DefaultXY_DataSet m1Line = new DefaultXY_DataSet();
			m1Line.set(m1, yRange.getLowerBound());
			m1Line.set(m1, yRange.getUpperBound());
			funcs.add(m1Line);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, Color.DARK_GRAY));
			
			DefaultXY_DataSet mMaxLine = new DefaultXY_DataSet();
			mMaxLine.set(mMax, yRange.getLowerBound());
			mMaxLine.set(mMax, yRange.getUpperBound());
			funcs.add(mMaxLine);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, Color.DARK_GRAY));
			
			DefaultXY_DataSet mc2Line = new DefaultXY_DataSet();
			mc2Line.set(6d, yRange.getLowerBound());
			mc2Line.set(6d, yRange.getUpperBound());
			funcs.add(mc2Line);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, Color.DARK_GRAY));
			
			XYTextAnnotation m1Ann = new XYTextAnnotation(" M₁=5", m1, yRange.getUpperBound());
			m1Ann.setFont(annFont);
			m1Ann.setTextAnchor(TextAnchor.TOP_LEFT);
			anns.add(m1Ann);
			
			XYTextAnnotation mMaxAnn = new XYTextAnnotation("Mₘₐₓ="+mMax.floatValue()+" ", mMax, yRange.getLowerBound());
			mMaxAnn.setFont(annFont);
			mMaxAnn.setTextAnchor(TextAnchor.BOTTOM_RIGHT);
			anns.add(mMaxAnn);
			
			XYTextAnnotation mc1Ann = new XYTextAnnotation(" Mc₁₉₇₃=5", m1, yRange.getLowerBound());
			mc1Ann.setFont(annFont);
			mc1Ann.setTextAnchor(TextAnchor.BOTTOM_LEFT);
			anns.add(mc1Ann);
			
			XYTextAnnotation mc2Ann = new XYTextAnnotation(" Mc₁₉₀₀=6", 6d, yRange.getLowerBound());
			mc2Ann.setFont(annFont);
			mc2Ann.setTextAnchor(TextAnchor.BOTTOM_LEFT);
			anns.add(mc2Ann);
			
			PlotSpec plot = new PlotSpec(funcs, chars, title, "Magnitude", incremental ? "Incremental Rate (1/yr)" : "Cumulative Rate (1/yr)");
			plot.setLegendInset(true);
			plot.setPlotAnnotations(anns);
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			gp.drawGraphPanel(plot, false, true, xRange, yRange);
			
			PlotUtils.writePlots(outputDir, prefix, gp, 700, 650, true, true, false);
			
			if (ref && refRatePairsFile != null && !incremental) {
				List<double[]> ratePairs = GriddedRateDistributionSolutionWriter.loadRates(refRatePairsFile);
				Collections.shuffle(ratePairs, new Random(ratePairs.size()));
				
//				int c = 200;
//				int a = 127;
//				int c = 150;
//				int a = 80;
				int c = 180;
				int a = 60;
				PlotCurveCharacterstics indvChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(c, c, c, a));
				int maxNumRates = 1000;
				int numRates = Integer.min(maxNumRates, ratePairs.size());
				GutenbergRichterMagFreqDist[] rateMFDs = new GutenbergRichterMagFreqDist[ratePairs.size()];
				for (int i=0; i<rateMFDs.length; i++) {
					double[] pair = ratePairs.get(i);
					GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(refMFD.getMinX()-0.5*refMFD.getDelta(), refMFD.size(), refMFD.getDelta());
					// this sets shape, min/max
					// subtract a tiny amount from mMax so that if it's exactly at a bin edge, e.g. 7.9, it rounds down, e.g. to 7.85
					gr.setAllButTotCumRate(gr.getX(0), gr.getMaxX(), 1e16, pair[1]);
					// this scales it to match
					// similarly, add a tiny amount to M1 so that if it's exactly at a bin edge (which it should be as it's determined
					// using cumulative binning), it rounds up to the incremental bin for that cumulative edge
//					gr.scaleToCumRate(refMFD.getClosestXIndex(m1+0.001), pair[0]);
					// this is actually cumulative, but the GR dist object works in incremental space
					gr.scaleToIncrRate(gr.getClosestXIndex(m1+0.001), pair[0]);
					gr.setName(null);
					
					rateMFDs[i] = gr;
					
					if (i < numRates) {
						funcs.add(0, gr);
						chars.add(0, indvChar);
						
						if (i == 0) {
							DiscretizedFunc forLegend = gr.deepClone();
							forLegend.scale(1e-10);
							forLegend.setName("Individual samples");
							funcs.add(0, forLegend);
							chars.add(0, new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, new Color(c, c, c))); // alpha disabled
						}
					}
				}
				
				gp.drawGraphPanel(plot, false, true, xRange, yRange);
				
				PlotUtils.writePlots(outputDir, prefix+"_with_indv", gp, 700, 650, true, true, false);
				
				// write histogram
				double[] histMags = {5d, 6d, 7d};
				
				for (double histMag : histMags) {
					double minRate = Double.POSITIVE_INFINITY;
					double maxRate = Double.NEGATIVE_INFINITY;
					
					double[] typeLowers = new double[types.length];
					double[] typeUppers = new double[types.length];
					
					for (int t=0; t<types.length; t++) {
						RateType type = types[t];
						EvenlyDiscretizedFunc[] typeFuncs = typeMagFuncs.get(type);
						int index = typeFuncs[0].getClosestXIndex(histMag);
						Preconditions.checkState(Precision.equals(histMag, typeFuncs[0].getX(index), 1e-4));
						typeLowers[t] = typeFuncs[0].getY(index);
						typeUppers[t] = typeFuncs[1].getY(index);
						minRate = Math.min(minRate, typeLowers[t]);
						maxRate = Math.max(maxRate, typeUppers[t]);
					}
					
					double[] histValues = new double[ratePairs.size()];
					for (int r=0; r<histValues.length; r++) {
						int index = rateMFDs[r].getClosestXIndex(histMag);
						Preconditions.checkState(Precision.equals(histMag, rateMFDs[r].getX(index), 1e-4),
								"Bad sample mag[%s]=%s for histMag=%s", index, rateMFDs[r].getX(index), histMag);
						histValues[r] = rateMFDs[r].getY(index);
						minRate = Math.min(minRate, histValues[r]);
						maxRate = Math.max(maxRate, histValues[r]);
					}
					
					int index = overallMean.getClosestXIndex(histMag);
					Preconditions.checkState(Precision.equals(histMag, overallMean.getX(index), 1e-4));
					double meanValue = overallMean.getY(index);
					
					
					double logMinRate = Math.log10(minRate) - 0.1;
					double logMaxRate = Math.log10(maxRate) + 0.1;
					double logMeanValue = Math.log10(meanValue);
					// center around mean value
					double maxDiff = Math.max(Math.abs(logMeanValue - logMaxRate), Math.abs(logMeanValue - logMinRate));
					logMinRate = logMeanValue - maxDiff;
					logMaxRate = logMeanValue + maxDiff;
					
					HistogramFunction logHist = HistogramFunction.getEncompassingHistogram(logMinRate, logMaxRate, 0.05);
					
					for (double histValue : histValues)
						logHist.add(logHist.getClosestXIndex(Math.log10(histValue)), 1d);
					
					double maxY = logHist.getMaxY()*1.1d;
					DiscretizedFunc linearHist = new ArbitrarilyDiscretizedFunc();
					for (Point2D pt : logHist)
						linearHist.set(Math.pow(10, pt.getX()), pt.getY());
					
					funcs = new ArrayList<>();
					chars = new ArrayList<>();
					
					linearHist.setName("Sampled Distribution");
					funcs.add(linearHist);
					chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GRAY));
					
					
					DefaultXY_DataSet meanXY = new DefaultXY_DataSet();
					meanXY.set(meanValue, 0d);
					meanXY.set(meanValue, maxY);
					meanXY.setName("Mean");
					funcs.add(meanXY);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));

					for (int t=0; t<types.length; t++) {
						RateType type = types[t];
						Color color = colors[type.ordinal()];
						PlotLineType plt = lineTypes[type.ordinal()];
						
						DefaultXY_DataSet lowerXY = new DefaultXY_DataSet();
						lowerXY.set(typeLowers[t], 0d);
						lowerXY.set(typeLowers[t], maxY);
						lowerXY.setName(type.toString());
						
						DefaultXY_DataSet upperXY = new DefaultXY_DataSet();
						upperXY.set(typeUppers[t], 0d);
						upperXY.set(typeUppers[t], maxY);
						
						funcs.add(lowerXY);
						chars.add(new PlotCurveCharacterstics(plt, 3f, color));
						
						funcs.add(upperXY);
						chars.add(new PlotCurveCharacterstics(plt, 3f, color));
					}
					
					plot = new PlotSpec(funcs, chars, " ", "M>"+oDF.format(histMag)+" Rate", "Sample Count");
					plot.setLegendInset(true);
					
					gp.drawGraphPanel(plot, true, false, new Range(Math.pow(10, logMinRate), Math.pow(10, logMaxRate)), new Range(0d, maxY));
					
					PlotUtils.writePlots(outputDir, prefix+"_hist_m"+oDF.format(histMag), gp, 700, 650, true, true, false);
				}
			}
		}
	}
	
	private static EvenlyDiscretizedFunc cmlMFD(RateRecord record, EvenlyDiscretizedFunc refMFD) {
		if (record.type == RateType.EXACT)
			return ((Exact)record).cumulativeDist;
		Preconditions.checkState(record instanceof PureGR);
		PureGR grRec = (PureGR)record;
		// fake a cml GR
		GutenbergRichterMagFreqDist grMFD = new GutenbergRichterMagFreqDist(
				grRec.b, 1d, refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
		grMFD.scaleToIncrRate(grMFD.getX(grMFD.getClosestXIndex(grRec.M1+0.01)), grRec.rateAboveM1);

		EvenlyDiscretizedFunc cmlGR = new EvenlyDiscretizedFunc(
				refMFD.getMinX()-0.5*refMFD.getDelta(), refMFD.size(), refMFD.getDelta());
		for (int i=0; i<cmlGR.size(); i++)
			cmlGR.set(i, grMFD.getY(i));
		return cmlGR;
	}

}
