package scratch.kevin.prvi25.figures;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
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

			Double m1 = null;
			Double mMax = null;
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
				
				lowMFD.setName(type.toString());
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
