package scratch.kevin.simCompare;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import javax.imageio.ImageIO;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.calc.GaussianDistCalc;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.function.UncertainArbDiscDataset;
import org.opensha.commons.gui.plot.AnimatedGIFRenderer;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.commons.util.FileUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGV_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SignificantDurationParam;
import org.opensha.sra.rtgm.RTGM;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;

import scratch.kevin.simulators.hazard.HazardMapComparePlotter;

public class SimulationHazardPlotter<E> {
	
	// plotting constants
	private static PlotLineType[] gmpe_alt_line_types = { PlotLineType.DASHED,
			PlotLineType.DOTTED, PlotLineType.DOTTED_AND_DASHED };

	private static PlotCurveCharacterstics gmpeCurveChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE);
	private static PlotCurveCharacterstics gmpeSpectraChar = new PlotCurveCharacterstics(
			PlotLineType.SOLID, 3f, PlotSymbol.FILLED_CIRCLE, 4.5f, Color.BLUE);
	
	private static PlotCurveCharacterstics simCurveChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK);
	private static PlotCurveCharacterstics simSpectraChar = new PlotCurveCharacterstics(
			PlotLineType.SOLID, 4f, PlotSymbol.FILLED_CIRCLE, 5.5f, Color.BLACK);
	private static PlotCurveCharacterstics simUncertainChar =
			new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, new Color(120, 120, 120, 120));
	private static PlotCurveCharacterstics[] simCompCurveChars = {
			new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.DARK_GRAY),
			new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.DARK_GRAY),
			new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, Color.DARK_GRAY),
			new PlotCurveCharacterstics(PlotLineType.DOTTED_AND_DASHED, 2f, Color.DARK_GRAY)
	};
	
	// plotting settings
	// these should be empty arrays if disabled, not null
	private int[] hazardCurveRPs = { 1000, 2500, 10000 };
	private double[] gmpeTruncs = { 3d, 2d, 1d };
//	private double[] gmpe_fixed_sigmas = { 0d, 0.3 };
	private double[] gmpeFixedSigmas = { 0d };
	private Map<SimulationHazardCurveCalc<?>, PlotCurveCharacterstics> customPlotChars = new HashMap<>();
	private Map<String, Range> curveXRanges;
	private Range curveYRange = new Range(1e-8, 1e0);

	// inputs
	private SimulationHazardCurveCalc<E> simCalc;
	private List<SimulationHazardCurveCalc<?>> compCalcs;
	private List<? extends RuptureComparison<E>> comps;
	private Site site;
	private double curveDuration;
	private AttenRelRef gmpeRef;
	private Map<String, DiscretizedFunc> xValsMap;
	
	private Table<String, E, Double> sourceRupContribFracts;
	private boolean sourceContribMutuallyExclusive;
	
	// results so far
	private Table<SimulationHazardCurveCalc<?>, IMT, DiscretizedFunc> simCurves = HashBasedTable.create();
	private Map<IMT, DiscretizedFunc> gmpeCurves = new HashMap<>();
	/**
	 * IMT, Trucation Level, Curve
	 */
	private Table<IMT, Double, DiscretizedFunc> gmpeTruncCurves = HashBasedTable.create();
	/**
	 * IMT, Sigma, Curve
	 */
	private Table<IMT, Double, DiscretizedFunc> gmpeFixedSigmaCurves = HashBasedTable.create();
	private Table<IMT, String, DiscretizedFunc> simSourceCurves = HashBasedTable.create();
	private Table<IMT, String, DiscretizedFunc> gmpeSourceCurves = HashBasedTable.create();
	
	public SimulationHazardPlotter(SimulationHazardCurveCalc<E> simCalc, List<? extends RuptureComparison<E>> comps,
			Site site, double curveDuration, AttenRelRef gmpeRef) {
		this(simCalc, null, comps, site, curveDuration, gmpeRef);
	}
	
	public SimulationHazardPlotter(SimulationHazardCurveCalc<E> simCalc, List<SimulationHazardCurveCalc<?>> compCalcs,
			List<? extends RuptureComparison<E>> comps, Site site, double curveDuration, AttenRelRef gmpeRef) {
		this.simCalc = simCalc;
		if (compCalcs == null)
			compCalcs = new ArrayList<>();
		this.compCalcs = compCalcs;
		this.comps = comps;
		this.site = site;
		this.curveDuration = curveDuration;
		this.gmpeRef = gmpeRef;

		this.xValsMap = simCalc.getXValsMap();
		curveXRanges = new HashMap<>();
		curveXRanges.put(SA_Param.NAME, new Range(1e-3, 1e1));
		curveXRanges.put(PGA_Param.NAME, new Range(1e-3, 1e1));
		curveXRanges.put(PGV_Param.NAME, new Range(1e-1, 1e3));
		curveXRanges.put(SignificantDurationParam.NAME, new Range(1e-1, 1e3));
	}
	
	public Site getSite() {
		return site;
	}
	
	public AttenRelRef getGMPE() {
		return gmpeRef;
	}
	
	public void setHazardCurveRPs(int[] hazardCurveRPs) {
		if (hazardCurveRPs == null)
			hazardCurveRPs = new int[0];
		this.hazardCurveRPs = hazardCurveRPs;
	}
	
	public void setGMPE_TruncationLevels(double[] gmpeTruncs) {
		if (gmpeTruncs == null)
			gmpeTruncs = new double[0];
		this.gmpeTruncs = gmpeTruncs;
	}
	
	public void setGMPE_FixedSigmas(double[] gmpeFixedSigmas) {
		if (gmpeFixedSigmas == null)
			gmpeFixedSigmas = new double[0];
		this.gmpeFixedSigmas = gmpeFixedSigmas;
	}
	
	public void setCustomPlotColors(SimulationHazardCurveCalc<?> simCalc, PlotCurveCharacterstics plotChar) {
		customPlotChars.put(simCalc, plotChar);
	}
	
	public void setSourceContributionFractions(Table<String, E, Double> sourceRupContribFracts, boolean mutuallyExclusive) {
		this.sourceRupContribFracts = sourceRupContribFracts;
		this.sourceContribMutuallyExclusive = mutuallyExclusive;
	}
	
	public synchronized DiscretizedFunc getCalcSimCurve(SimulationHazardCurveCalc<?> simCalc, IMT imt) throws IOException {
		if (simCurves.contains(simCalc, imt))
			return simCurves.get(simCalc, imt);
		DiscretizedFunc curve = simCalc.calc(site, imt, curveDuration);
		simCurves.put(simCalc, imt, curve);
		return curve;
	}
	
	private DiscretizedFunc getGMPECurve(IMT imt, double sigmaTruncation, double fixedSigma) {
		Preconditions.checkState(sigmaTruncation < 0 || fixedSigma < 0);
		if (sigmaTruncation < 0 && fixedSigma < 0 && gmpeCurves.containsKey(imt))
			// regular gmpe curve
			return gmpeCurves.get(imt);
		if (sigmaTruncation > 0 && gmpeTruncCurves.contains(imt, sigmaTruncation))
			// truncated curve
			return gmpeTruncCurves.get(imt, sigmaTruncation);
		if (fixedSigma >= 0 && gmpeFixedSigmaCurves.contains(imt, fixedSigma))
			// fixed sigma curve
			return gmpeFixedSigmaCurves.get(imt, fixedSigma);
		return null; // not calculated
	}
	
	public synchronized DiscretizedFunc getCalcGMPECurve(IMT imt) {
		return getCalcGMPECurve(imt, -1, -1);
	}
	
	public synchronized DiscretizedFunc getCalcGMPECurve(IMT imt, double sigmaTruncation, double fixedSigma) {
		DiscretizedFunc curve = getGMPECurve(imt, sigmaTruncation, fixedSigma);
		if (curve != null)
			return curve;
		
		// we don't have the curve we need, calculate them all for this period
		
		DiscretizedFunc gmpeCurve = xValsMap.get(imt.getParamName()).deepClone();
		// init to 1, non-exceedance curves
		for (int i=0; i<gmpeCurve.size(); i++)
			gmpeCurve.set(i, 1d);
		
		DiscretizedFunc[] calcTruncCurves = new DiscretizedFunc[gmpeTruncs.length];
		for (int t=0; t<calcTruncCurves.length; t++) {
			calcTruncCurves[t] = gmpeCurve.deepClone();
			Preconditions.checkState(gmpeTruncs[t] > 0, "Truncation level not positive: %s", gmpeTruncs[t]);
			// init to 1, non-exceedance curves
			for (int i=0; i<calcTruncCurves[t].size(); i++)
				calcTruncCurves[t].set(i, 1d);
		}
		
		DiscretizedFunc[] calcFixedSigmaCurves = new DiscretizedFunc[gmpeFixedSigmas.length];
		for (int t=0; t<calcFixedSigmaCurves.length; t++) {
			calcFixedSigmaCurves[t] = gmpeCurve.deepClone();
			Preconditions.checkState(gmpeFixedSigmas[t] >= 0, "Fixed sigma level is negative: %s", gmpeFixedSigmas[t]);
			if (gmpeFixedSigmas[t] > 0) {
				// init to 1, non-exceedance curves
				for (int i=0; i<calcFixedSigmaCurves[t].size(); i++)
					calcFixedSigmaCurves[t].set(i, 1d);
			} else {
				// init to 0
				for (int i=0; i<calcFixedSigmaCurves[t].size(); i++)
					calcFixedSigmaCurves[t].set(i, 0);
			}
		}
		
		DiscretizedFunc logXVals = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : gmpeCurve)
			logXVals.set(Math.log(pt.getX()), 1d);
		logXVals = new LightFixedXFunc(logXVals);
		
		for (RuptureComparison<E> comp : new RuptureComparisonFilter.SiteFilter<E>().getMatches(comps, site)) {
			double rupRate = comp.getAnnualRate();
			if (rupRate == 0)
				continue;
			double rupProb = 1d - Math.exp(-rupRate*curveDuration);
			
			double mean = comp.getLogMean(site, imt);
			double stdDev = comp.getStdDev(site, imt);
			
			// regular curve
			for(int k=0; k<logXVals.size(); k++) {
				double stRndVar = (logXVals.getX(k) - mean) / stdDev;
				double exceedProb = GaussianDistCalc.getExceedProb(stRndVar);
				gmpeCurve.set(k, gmpeCurve.getY(k)*Math.pow(1d-rupProb, exceedProb));
			}
			
			// truncated curves
			for (int i=0; gmpeTruncs != null && i<calcTruncCurves.length; i++) {
				double trunc = gmpeTruncs[i];
				for(int k=0; k<logXVals.size(); k++) {
					double stRndVar = (logXVals.getX(k) - mean) / stdDev;
					double exceedProb = GaussianDistCalc.getExceedProb(stRndVar, 1, trunc);
					calcTruncCurves[i].set(k, calcTruncCurves[i].getY(k)*Math.pow(1d-rupProb, exceedProb));
				}
			}
			
			// fixed sigma curves
			for (int i=0; gmpeFixedSigmas != null && i<calcFixedSigmaCurves.length; i++) {
				double sigma = gmpeFixedSigmas[i];
				if (sigma > 0) {
					for(int k=0; k<logXVals.size(); k++) {
						double stRndVar = (logXVals.getX(k) - mean) / sigma;
						double exceedProb = GaussianDistCalc.getExceedProb(stRndVar);
						calcFixedSigmaCurves[i].set(k, calcFixedSigmaCurves[i].getY(k)*Math.pow(1d-rupProb, exceedProb));
					}
				} else {
					Preconditions.checkState(Double.isFinite(mean));
					Preconditions.checkState(sigma == 0d);
					for(int k=0; k<logXVals.size(); k++)
						if (mean >= logXVals.getX(k))
							calcFixedSigmaCurves[i].set(k, calcFixedSigmaCurves[i].getY(k) + rupRate);
				}
			}
		}
		// convert to exceedance probabilities
		for (int i=0; i<gmpeCurve.size(); i++)
			gmpeCurve.set(i, 1d-gmpeCurve.getY(i));
		gmpeCurves.put(imt, gmpeCurve);
		for (int t=0; t<calcTruncCurves.length; t++) {
			for (int i=0; i<calcTruncCurves[t].size(); i++)
				calcTruncCurves[t].set(i, 1d-calcTruncCurves[t].getY(i));
			gmpeTruncCurves.put(imt, gmpeTruncs[t], calcTruncCurves[t]);
		}
		for (int t=0; t< calcFixedSigmaCurves.length; t++) {
			if (gmpeFixedSigmas[t] > 0) {
				for (int i=0; i<calcFixedSigmaCurves[t].size(); i++)
					calcFixedSigmaCurves[t].set(i, 1d-calcFixedSigmaCurves[t].getY(i));
			} else {
				// now mean only curve -> probabilities
				for (int i=0; i<calcFixedSigmaCurves[t].size(); i++) {
					double rate = calcFixedSigmaCurves[t].getY(i);
					double prob = 1d - Math.exp(-rate*curveDuration);
					calcFixedSigmaCurves[t].set(i, prob);
				}
			}
			gmpeFixedSigmaCurves.put(imt, gmpeFixedSigmas[t], calcFixedSigmaCurves[t]);
		}
		
		// now return it
		return getGMPECurve(imt, sigmaTruncation, fixedSigma);
	}
	
	public File plotHazardCurves(File outputDir, String prefix, double period) throws IOException {
		return plotHazardCurves(outputDir, prefix, IMT.forPeriod(period));
	}
	
	public File plotHazardCurves(File outputDir, String prefix, IMT imt) throws IOException {
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		DiscretizedFunc gmpeCurve = getCalcGMPECurve(imt);
		
		// primary GMPE curve
		gmpeCurve.setName(gmpeRef.getShortName());
		funcs.add(gmpeCurve);
		chars.add(gmpeCurveChar);
		
		// GMPE truncated curves
		for (int i=0; i<gmpeTruncs.length; i++) {
			DiscretizedFunc gmpeTruncCurve = getCalcGMPECurve(imt, gmpeTruncs[i], -1);
			String name = optionalDigitDF.format(gmpeTruncs[i])+" σₜ";
			float thickness = 2f;
			gmpeTruncCurve.setName(name);
			funcs.add(gmpeTruncCurve);
			chars.add(new PlotCurveCharacterstics(gmpe_alt_line_types[i % gmpe_alt_line_types.length], thickness, Color.BLUE));
		}
		
		// GMPE fixed sigma curves
		for (int i=0; i<gmpeFixedSigmas.length; i++) {
			DiscretizedFunc gmpeFixedSigmaCurve = getCalcGMPECurve(imt, -1, gmpeFixedSigmas[i]);
			String name = "σ="+optionalDigitDF.format(gmpeFixedSigmas[i]);
			if (i == 0)
				name = "Fixed "+name;
			float thickness = 2f;
			gmpeFixedSigmaCurve.setName(name);
			funcs.add(gmpeFixedSigmaCurve);
			chars.add(new PlotCurveCharacterstics(gmpe_alt_line_types[i % gmpe_alt_line_types.length],
					thickness, Color.GREEN.darker()));
		}
		
		// primary simulation curve
		addSimulationCurves(funcs, chars, imt, true, true);
		
		return plotHazardCurves(outputDir, prefix, site, imt, curveDuration, funcs, chars, null, true);
	}
	
	public File plotGMPE_SimHazardCurves(File outputDir, String prefix, IMT imt,
			GMPESimulationBasedProvider<E> gmpeSimProv, int numGMPESims) throws IOException {
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		DiscretizedFunc[] gmpeSimCurves = new DiscretizedFunc[numGMPESims];
		SimulationHazardCurveCalc<E> gmpeSimCalc = new SimulationHazardCurveCalc<>(
				gmpeSimProv, xValsMap);
		for (int i=0; i<gmpeSimCurves.length; i++) {
			gmpeSimProv.clearCache();
			gmpeSimCurves[i] = gmpeSimCalc.calc(site, imt, curveDuration);
		}
		
		DiscretizedFunc meanCurve = gmpeSimCalc.getXValsMap().get(imt.getParamName()).deepClone();
		for (int i=0; i<meanCurve.size(); i++)
			meanCurve.set(i, 0);
		meanCurve.setName(gmpeRef.getShortName()+" mean, N="+gmpeSimCurves.length);
		Color singleSimColor = new Color(255, 120, 120);
		Color meanSimColor = Color.RED;
		
		for (int n=0; n<gmpeSimCurves.length; n++) {
			DiscretizedFunc simCurve = gmpeSimCurves[n];
			Preconditions.checkState(simCurve.size() > 0, "Empty sim curve?");
			if (n == 0) {
				if (gmpeSimCurves.length == 0)
					simCurve.setName(gmpeRef.getShortName()+" Simulated");
				else
					simCurve.setName(gmpeRef.getShortName()+" Simulations");
			} else 
				simCurve.setName(null);
			Preconditions.checkState(meanCurve.size() >= simCurve.size(),
					"Size mismatch for n=%s, sim=%s, mean=%s", n, simCurve.size(), meanCurve.size());
			for (int i=0; i<simCurve.size(); i++)
				meanCurve.set(i, meanCurve.getY(i)+simCurve.getY(i));
			funcs.add(simCurve);
			if (gmpeSimCurves.length == 0)
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, meanSimColor));
			else
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, singleSimColor));
		}
		if (gmpeSimCurves.length > 1) {
			meanCurve.scale(1d/gmpeSimCurves.length);
			ArbitrarilyDiscretizedFunc trimmedMeanCurve = new ArbitrarilyDiscretizedFunc();
			for (Point2D pt : meanCurve)
				if (pt.getY() > 0)
					trimmedMeanCurve.set(pt);
			Preconditions.checkState(trimmedMeanCurve.size() > 0, "Empty trimmed mean curve? Orig was %s", meanCurve);
			trimmedMeanCurve.setName(meanCurve.getName());
			funcs.add(trimmedMeanCurve);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, meanSimColor));
		}
		
		DiscretizedFunc gmpeCurve = getCalcGMPECurve(imt);
		
		// primary GMPE curve
		gmpeCurve.setName(gmpeRef.getShortName());
		funcs.add(gmpeCurve);
		chars.add(gmpeCurveChar);
		
		// primary simulation curve
		addSimulationCurves(funcs, chars, imt, true, false);
		
		return plotHazardCurves(outputDir, prefix, site, imt, curveDuration, funcs, chars, null, false);
	}
	
	private synchronized Map<String, DiscretizedFunc> getCalcSimSourceCurves(IMT imt) throws IOException {
		if (simSourceCurves.containsRow(imt))
			return simSourceCurves.row(imt);
		Map<String, DiscretizedFunc> curves = simCalc.calcSourceContributionCurves(
				site, imt, curveDuration, sourceRupContribFracts);
		for (String sourceName : curves.keySet())
			simSourceCurves.put(imt, sourceName, curves.get(sourceName));
		return curves;
	}
	
	private synchronized Map<String, DiscretizedFunc> getCalcGMPESourceCurves(IMT imt)
			throws IOException {
		if (gmpeSourceCurves.containsRow(imt))
			return gmpeSourceCurves.row(imt);
		
		Map<String, DiscretizedFunc> curves = new HashMap<>();
		for (String sourceName : sourceRupContribFracts.rowKeySet()) {
			Map<E, Double> rupContribFracts = sourceRupContribFracts.row(sourceName);
			DiscretizedFunc simSourceCurve = getCalcSimSourceCurves(imt).get(sourceName);
			if (simSourceCurve == null || simSourceCurve.size() == 0) {
				System.out.println("Skipping "+sourceName+" (null ? "+(simSourceCurve == null)+")");
				continue;
			}
//			if (contribIML > 0) {
//				if (simSourceCurve.getMaxX() < contribIML)
//					continue;
//				double simSourceVal = simSourceCurve.getInterpolatedY_inLogXLogYDomain(contribIML);
//				if (!Double.isFinite(simSourceVal) || simSourceVal <= 0)
//					continue;
//			}
			List<ModRateRuptureComparison<E>> modComps = new ArrayList<>();
			for (RuptureComparison<E> comp : new RuptureComparisonFilter.SiteFilter<E>().getMatches(comps, site)) {
				double rupRate = comp.getAnnualRate();
				Double fract = rupContribFracts.get(comp.getRupture());
				if (fract == null || fract <= 0 || rupRate == 0)
					continue;
				rupRate *= fract;
				modComps.add(new ModRateRuptureComparison<>(comp, rupRate));
			}
			
			SimulationHazardPlotter<E> modCalc = new SimulationHazardPlotter<>(simCalc, modComps, site, curveDuration, gmpeRef);
			modCalc.setGMPE_FixedSigmas(null);
			modCalc.setGMPE_TruncationLevels(null);
			DiscretizedFunc gmpeCurve = modCalc.getCalcGMPECurve(imt);
			
			curves.put(sourceName, gmpeCurve);
		}
		for (String sourceName : curves.keySet())
			gmpeSourceCurves.put(imt, sourceName, curves.get(sourceName));
		return curves;
	}
	
	public static double calcRTGM(DiscretizedFunc curve, double curveDuration) {
		// first discretize evenly
		DiscretizedFunc discrCurve;
		if (curve instanceof EvenlyDiscretizedFunc) {
			discrCurve = curve;
		} else {
			discrCurve = new EvenlyDiscretizedFunc(
					curve.getMinX(), curve.getMaxX(), 100);
			for (int i=0; i<discrCurve.size(); i++) {
				double x = discrCurve.getX(i);
				double y;
				if ((float)x == 0f)
					y = curve.getY(0);
				else if ((float)x == (float)curve.getMaxX())
					y = curve.getY(curve.size()-1);
				else
					y = curve.getInterpolatedY(x);
				y = - Math.log(1-y)/curveDuration;
				discrCurve.set(i, y);
			}
		}
		if (RTGM.AFE4UHGM < discrCurve.getMinY() || RTGM.AFE4UHGM > discrCurve.getMaxY()) {
			return 0d;
		}
		RTGM calc = RTGM.create(discrCurve, null, null);
		try {
			calc.call();
		} catch (RuntimeException e) {
			System.err.println("RTGM Calc failed for Hazard Curve:\n"+discrCurve);
			System.err.flush();
			throw e;
		}
		double rtgm = calc.get();
		if (!Double.isFinite(rtgm)) {
			System.err.println("WARNING: could not calculate RTGM for "+curve.getName()+": "+rtgm);
			rtgm = 0d;
		}
//		Preconditions.checkState(rtgm >= 0, "RTGM is not >=0: %s\n\n%s", rtgm, discrCurve);
		return rtgm;
	}
	
	public File plotSourceContributionHazardCurves(File outputDir, String prefix, IMT imt,
			double sourceContribSortProb, int numSourceCurves, boolean gmpeSources) throws IOException {
		Preconditions.checkState(sourceRupContribFracts != null, "Source rup contribution fractions not set!");
		Map<String, Double> simSourceSortVals = new HashMap<>();
		
		DiscretizedFunc refCurve = getCalcSimCurve(simCalc, imt);
		double contribIML;
		if (sourceContribSortProb > 0) {
			if (sourceContribSortProb < refCurve.getMinY() || sourceContribSortProb > refCurve.getMaxY()) {
				System.out.println("Source contribution sort probability ("+(float)sourceContribSortProb
						+") is not contained in simulation curve, can't determine ground motion");
				return null;
			}
			contribIML = refCurve.getFirstInterpolatedX_inLogXLogYDomain(sourceContribSortProb);
			System.out.println("Sorting source contirubtions for IML="+(float)+contribIML
					+" at p="+(float)sourceContribSortProb);
		} else {
			// use RTGM
			contribIML = 0d;
		}
		
		Map<String, DiscretizedFunc> simSourceCurves = getCalcSimSourceCurves(imt);
		Map<String, DiscretizedFunc> gmpeSourceCurves =
				gmpeSources ? getCalcGMPESourceCurves(imt) : null;
		
		for (String sourceName : simSourceCurves.keySet()) {
			DiscretizedFunc simSourceCurve = simSourceCurves.get(sourceName);
			if (simSourceCurve == null || simSourceCurve.size() == 0)
				continue;
			double simSourceVal;
			if (sourceContribSortProb > 0) {
				if (simSourceCurve.getMaxX() < contribIML)
					continue;
				simSourceVal = simSourceCurve.getInterpolatedY_inLogXLogYDomain(contribIML);
			} else {
				// calculate RTGM
				simSourceVal = calcRTGM(simSourceCurve, curveDuration);
//				System.out.println(sourceName+" RTGM: "+simSourceVal+", y(0): "+discrCurve.getY(0));
			}
			
			simSourceSortVals.put(sourceName, simSourceVal);
		}
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		// primary GMPE curve
		if (gmpeSources) {
			DiscretizedFunc gmpeCurve = getCalcGMPECurve(imt, -1, -1);
			gmpeCurve.setName(gmpeRef.getShortName());
			funcs.add(gmpeCurve);
			chars.add(gmpeCurveChar);
		} else {
			// primary simulation curve
			addSimulationCurves(funcs, chars, imt, true, false);
		}
		
		addSourceContributionFuncs(funcs, chars, simSourceCurves, gmpeSourceCurves, simSourceSortVals,
				numSourceCurves, sourceContribMutuallyExclusive, gmpeSources, false, imt);
		
		return plotHazardCurves(outputDir, prefix, site, imt, curveDuration, funcs, chars, null, true);
	}
	
	private void addSourceContributionFuncs(List<DiscretizedFunc> funcs, List<PlotCurveCharacterstics> chars,
			Map<String, DiscretizedFunc> simSourceCurves, Map<String, DiscretizedFunc> gmpeSourceCurves,
			Map<String, Double> simSourceSortVals, int numSourceCurves, boolean sourceContribMutuallyExclusive,
			boolean gmpeSources, boolean symbols, IMT imt) throws IOException {
		List<String> sortedSourceNames = ComparablePairing.getSortedData(simSourceSortVals);
		Collections.reverse(sortedSourceNames);
		List<DiscretizedFunc> otherCurves = new ArrayList<>();
		
		int numSoFar = 0;
		List<DiscretizedFunc> curvesToAdd = new ArrayList<>();
		for (String sourceName : sortedSourceNames) {
			numSoFar++;
			DiscretizedFunc curve;
			if (gmpeSources)
				curve = gmpeSourceCurves.get(sourceName);
			else
				curve = simSourceCurves.get(sourceName);
			Preconditions.checkNotNull(curve, "No curve for %s, gmpeSources=%s", sourceName, gmpeSources);
			curve.setName(sourceName);
			double val = simSourceSortVals.get(sourceName);
			if (numSoFar >= numSourceCurves || val == 0)
				otherCurves.add(curve);
			else
				curvesToAdd.add(curve);
		}
		Preconditions.checkState(!curvesToAdd.isEmpty());
		if (curvesToAdd.size() == 1) {
			funcs.add(curvesToAdd.get(0));
			Color color = GMT_CPT_Files.MAX_SPECTRUM.instance().getMaxColor().darker();
			if (symbols)
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, PlotSymbol.FILLED_CIRCLE, 3f, color));
			else
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, color));
		} else {
			CPT colorCPT = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0, curvesToAdd.size()-1);
			colorCPT = colorCPT.reverse();
			colorCPT.setBelowMinColor(colorCPT.getMinColor());
			colorCPT.setAboveMaxColor(colorCPT.getMaxColor());
			for (int i=0; i <curvesToAdd.size(); i++) {
				funcs.add(curvesToAdd.get(i));
				Color color = colorCPT.getColor((float)i).darker();
				if (symbols)
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, PlotSymbol.FILLED_CIRCLE, 3f, color));
				else
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, color));
			}
		}
		
		if (sourceContribMutuallyExclusive) {
			// now build other sources curves
			if (!otherCurves.isEmpty()) {
				DiscretizedFunc otherCurve = xValsMap.get(imt.getParamName()).deepClone();
				// init to 1, do in non-exceed space
				for (int i=0; i<otherCurve.size(); i++)
					otherCurve.set(i, 1d);
				for (DiscretizedFunc curve : otherCurves) {
					for (int i=0; i<curve.size(); i++)
						otherCurve.set(i, otherCurve.getY(i)*(1d-curve.getY(i)));
				}
				// convert to exceed probs
				for (int i=0; i<otherCurve.size(); i++)
					otherCurve.set(i, 1d-otherCurve.getY(i));
				otherCurve.setName(otherCurves.size()+" Other Sources");
				funcs.add(otherCurve);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GRAY));
			}
		}
	}
	
	private void addSimulationCurves(List<DiscretizedFunc> funcs, List<PlotCurveCharacterstics> chars,
			IMT imt, boolean includeUncertainty, boolean includeComps) throws IOException {
		DiscretizedFunc simCurve = getCalcSimCurve(simCalc, imt);
		if (includeUncertainty && simCurve instanceof UncertainArbDiscDataset) {
			// this is so that it is first in the legend
			funcs.add(simCurve);
			chars.add(simCurveChar);
			// now clone so that we can change the name
			UncertainArbDiscDataset uncertain = (UncertainArbDiscDataset)simCurve;
			uncertain = new UncertainArbDiscDataset(simCurve, uncertain.getLower(), uncertain.getUpper());
			uncertain.setName("95% conf");
			funcs.add(uncertain);
			chars.add(simUncertainChar);
		} else {
			funcs.add(simCurve);
			chars.add(simCurveChar);
		}
		
		int charIndex = 0;
		for (int i=0; i<compCalcs.size() && includeComps; i++) {
			SimulationHazardCurveCalc<?> compCalc = compCalcs.get(i);
			DiscretizedFunc curve = getCalcSimCurve(compCalc, imt);
			funcs.add(curve);
			if (customPlotChars.containsKey(compCalc))
				chars.add(customPlotChars.get(compCalc));
			else
				chars.add(simCompCurveChars[charIndex++ % simCompCurveChars.length]);
		}
		
		// put final simulation curve back on top, without an extra legend
		DiscretizedFunc curveNoLegend = simCurve.deepClone();
		curveNoLegend.setName(null);
		funcs.add(curveNoLegend);
		chars.add(simCurveChar);
	}
	
	public void plotCurveAnimation(File outputFile, double[] timeRange, double delta, IMT imt) throws IOException {
		Preconditions.checkState(delta > 0);
		
		List<DiscretizedFunc> prevSimCurves = new ArrayList<>();
		List<DiscretizedFunc> prevGMPECurves = new ArrayList<>();
		List<Double> prevRelLengths = new ArrayList<>();
		
		double totalRange = timeRange[1] - timeRange[0];
		int numFrames = (int)Math.ceil(totalRange / delta);
		
		double fps = 1;
		
		System.out.println("Creating animation with "+numFrames+" frames, "+optionalDigitDF.format(delta)+" years each");

		CPT prevSimCPT = new CPT(0, numFrames, new Color(100, 100, 100, 180), new Color(100, 100, 100, 20));
		CPT finalSimCPT = new CPT(0, numFrames, Color.DARK_GRAY, Color.LIGHT_GRAY);
		CPT prevGMPECPT = new CPT(0, numFrames, new Color(100, 100, 255, 180), new Color(100, 100, 255, 20));
		
		File tempDir = Files.createTempDir();
		
		List<File> imageFiles = new ArrayList<>();
		
		for (int frameIndex=0; frameIndex<numFrames; frameIndex++) {
			double windowEnd = timeRange[0] + (frameIndex + 1)*delta;
			double relativeLength = delta*(frameIndex+1);
			if (windowEnd > timeRange[1])
				relativeLength -= (windowEnd - timeRange[1]);
			System.out.println("Animation frame "+frameIndex+", "+optionalDigitDF.format(relativeLength)
				+" yr ("+optionalDigitDF.format(timeRange[0])+" => "+optionalDigitDF.format(windowEnd)+")");
			
			double rateScalar = (timeRange[1]-timeRange[0])/(windowEnd - timeRange[0]);
			System.out.println("\tRate scale: "+(float)rateScalar);
			
			List<RuptureComparison<E>> subComps = new ArrayList<>();
			HashSet<E> ruptures = new HashSet<>();
			for (RuptureComparison<E> comp : comps) {
				double time = comp.getRuptureTimeYears();
				if (time < windowEnd) {
					subComps.add(new ModRateRuptureComparison<>(comp, comp.getAnnualRate()*rateScalar));
					ruptures.add(comp.getRupture());
				}
			}
			Preconditions.checkState(!subComps.isEmpty(), "No ruptures found in window!");
			
			SubSetSimulationRotDProvider<E> subSetProv = new SubSetSimulationRotDProvider<>(simCalc.getSimProv(), ruptures, rateScalar);
			SimulationHazardCurveCalc<E> subSetCalc = new SimulationHazardCurveCalc<>(subSetProv);
			
			DiscretizedFunc simCurve = subSetCalc.calc(site, imt, curveDuration);
			SimulationHazardPlotter<E> subSetPlotter = new SimulationHazardPlotter<>(simCalc, subComps, site, curveDuration, gmpeRef);
			subSetPlotter.setGMPE_FixedSigmas(null);
			subSetPlotter.setGMPE_TruncationLevels(null);
			DiscretizedFunc gmpeCurve = subSetPlotter.getCalcGMPECurve(imt);
			
			List<DiscretizedFunc> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			if (frameIndex > 0) {
				for (int i=0; i<prevSimCurves.size(); i++) {
					int generation = prevSimCurves.size() - i;
					
					DiscretizedFunc prevSimCurve = prevSimCurves.get(i);
					prevSimCurve.setName(null);
					funcs.add(prevSimCurve);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, prevSimCPT.getColor((float)generation)));
					
					DiscretizedFunc prevGMPECurve = prevGMPECurves.get(i);
					prevGMPECurve.setName(null);
					funcs.add(prevGMPECurve);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1.5f, prevGMPECPT.getColor((float)generation)));
				}
			}
			
			funcs.add(gmpeCurve);
			chars.add(gmpeCurveChar);
			
			funcs.add(simCurve);
			chars.add(simCurveChar);
			
			List<XYAnnotation> anns = new ArrayList<>();
			XYTextAnnotation durationAnn = new XYTextAnnotation(groupedIntDF.format(relativeLength)+" years", 2e-3, 3e-8);
			durationAnn.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 26));
			durationAnn.setTextAnchor(TextAnchor.BASELINE_LEFT);
			anns.add(durationAnn);
			
			imageFiles.add(plotHazardCurves(tempDir, "frame_"+frameIndex, site, imt, curveDuration, funcs, chars, anns, false));
			
			if (frameIndex == numFrames-1 && numFrames > 1) {
				funcs = new ArrayList<>();
				chars = new ArrayList<>();
				
				for (int i=0; i<prevSimCurves.size(); i++) {
					int generation = prevSimCurves.size() - i;
					
					DiscretizedFunc prevSimCurve = prevSimCurves.get(i);
//					prevSimCurve.setName(groupedIntDF.format(prevRelLengths.get(i))+" yrs");
					if (i == 0)
						prevSimCurve.setName("Subset Curves");
					else
						prevSimCurve.setName(null);
					funcs.add(prevSimCurve);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, finalSimCPT.getColor((float)generation)));
//					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.DARK_GRAY));
				}
				
//				funcs.add(gmpeCurve);
//				chars.add(gmpeCurveChar);
				
				simCurve.setName("Complete Model");
				funcs.add(simCurve);
				chars.add(simCurveChar);
				
				double minY = 1d;
				double maxX = 0d;
				for (Point2D pt : simCurve) {
					if (pt.getY() > 0) {
						minY = Math.min(minY, pt.getY());
						maxX = Math.max(maxX, pt.getX());
					}
				}
				double maxY = 1e-2;
				minY = Math.pow(10, Math.floor(Math.log10(minY)-0.2));
				double minX = simCurve.getFirstInterpolatedX(maxY);
				maxX = Math.pow(10, Math.ceil(Math.log10(maxX)+0.1));
				
				String prefix = outputFile.getName().replaceAll(".gif", "")+"_final";
				plotHazardCurves(outputFile.getParentFile(), prefix, site, imt, curveDuration,
						funcs, chars, null, true, new Range(minX, maxX), new Range(minY, maxY));
				
				// now a subset
				funcs = new ArrayList<>();
				chars = new ArrayList<>();
				
				int mod = 1;
				
				for (int i=0; i<prevSimCurves.size(); i++) {
					if (i == 1)
						mod = 2;
					else if (i == 2)
						mod = 4;
					if ((i+1) % mod > 0)
						continue;
					
					DiscretizedFunc prevSimCurve = prevSimCurves.get(i);
					prevSimCurve.setName(groupedIntDF.format(prevRelLengths.get(i))+" yrs");
					funcs.add(prevSimCurve);
//					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, finalSimCPT.getColor((float)generation)));
//					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.DARK_GRAY));
				}
				
				CPT summaryCPT = new CPT(0, Integer.max(1, funcs.size()-1), Color.LIGHT_GRAY, Color.DARK_GRAY);
				for (int i=0; i<funcs.size(); i++)
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, summaryCPT.getColor((float)i)));
				
//				funcs.add(gmpeCurve);
//				chars.add(gmpeCurveChar);
				
				simCurve.setName("Complete Model");
				funcs.add(simCurve);
				chars.add(simCurveChar);
				
				prefix = outputFile.getName().replaceAll(".gif", "")+"_summary";
				plotHazardCurves(outputFile.getParentFile(), prefix, site, imt, curveDuration,
						funcs, chars, null, true, new Range(minX, maxX), new Range(minY, maxY));
			}
			
			prevGMPECurves.add(gmpeCurve);
			prevSimCurves.add(simCurve);
			prevRelLengths.add(relativeLength);
		}
		
		AnimatedGIFRenderer gifRender = new AnimatedGIFRenderer(outputFile, fps, true);
		
		for (File imageFile : imageFiles) {
			BufferedImage img = ImageIO.read(imageFile);
			
			gifRender.writeFrame(img);
		}
		
		gifRender.finalizeAnimation();
		
		FileUtils.deleteRecursive(tempDir);
	}
	
	File plotHazardCurves(File outputDir, String prefix, Site site, IMT imt, double curveDuration,
			List<DiscretizedFunc> funcs, List<PlotCurveCharacterstics> chars, List<XYAnnotation> anns,
			boolean pubFig) throws IOException {
		return plotHazardCurves(outputDir, prefix, site, imt, curveDuration, funcs, chars, anns, pubFig,
				curveXRanges.get(imt.getParamName()), curveYRange);
	}
	
	File plotHazardCurves(File outputDir, String prefix, Site site, IMT imt, double curveDuration,
			List<DiscretizedFunc> funcs, List<PlotCurveCharacterstics> chars, List<XYAnnotation> anns,
			boolean pubFig, Range curveXRange, Range curveYRange) throws IOException {
		// now plot
		String xAxisLabel = imt.getDisplayName()+" ("+imt.getUnits()+")";
		String yAxisLabel;
		if (curveDuration == 1d)
			yAxisLabel = "Annual Probability of Exceedance";
		else
			yAxisLabel = optionalDigitDF.format(curveDuration)+"yr Probability of Exceedance";
		
		if (hazardCurveRPs != null && hazardCurveRPs.length > 0) {
			if (anns == null)
				anns = new ArrayList<>();
			CPT rpCPT = HazardMapComparePlotter.getRPlogCPT(hazardCurveRPs);
			Font font = new Font(Font.SANS_SERIF, Font.BOLD, 16);
			for (int rp : hazardCurveRPs) {
				Color color = rpCPT.getColor((float)Math.log10(rp));
				double probLevel = curveDuration/(double)rp;
				DiscretizedFunc probLine = new ArbitrarilyDiscretizedFunc();
				probLine.set(curveXRange.getLowerBound(), probLevel);
				probLine.set(curveXRange.getUpperBound(), probLevel);
				funcs.add(probLine);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, color));
				
				XYTextAnnotation ann = new XYTextAnnotation(" "+rp+"yr", curveXRange.getLowerBound(), probLevel);
				ann.setTextAnchor(TextAnchor.BOTTOM_LEFT);
				ann.setFont(font);
				ann.setPaint(Color.DARK_GRAY);
				anns.add(ann);
			}
		}

		PlotSpec spec = new PlotSpec(funcs, chars, site.getName()+" Hazard Curves", xAxisLabel, yAxisLabel);
		spec.setLegendVisible(true);
		spec.setPlotAnnotations(anns);
		
		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(18);
		plotPrefs.setAxisLabelFontSize(20);
		plotPrefs.setPlotLabelFontSize(21);
		plotPrefs.setBackgroundColor(Color.WHITE);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
		
		gp.drawGraphPanel(spec, true, true, curveXRange, curveYRange);
		gp.getChartPanel().setSize(800, 600);
		File pngFile = new File(outputDir, prefix+".png");
		gp.saveAsPNG(pngFile.getAbsolutePath());
		File pdfFile = new File(outputDir, prefix+".pdf");
		gp.saveAsPDF(pdfFile.getAbsolutePath());
		
		if (pubFig) {
			// write version without title, for publications
			spec.setTitle(" ");
			gp.drawGraphPanel(spec, true, true, curveXRange, curveYRange);
			gp.getChartPanel().setSize(800, 600);
			File pubPNGFile = new File(outputDir, prefix+"_pub.png");
			gp.saveAsPNG(pubPNGFile.getAbsolutePath());
			File pubPDFFile = new File(outputDir, prefix+"_pub.pdf");
			gp.saveAsPDF(pubPDFFile.getAbsolutePath());
		}
		System.out.println("DONE "+prefix+", "+site.getName()+", "+xAxisLabel);
		return pngFile;
	}
	
	public List<String> getCurveLegend(boolean spectra, boolean hasGMPETruncs, boolean hasGMPESigmas, int numGMPESims) {
		List<String> lines = new ArrayList<>();
		String typeStr = spectra ? "Spectra" : "Curves";
		lines.add("**Legend**:");
		if (spectra)
			lines.add("* **Simulation "+typeStr+"**");
		else
			lines.add("* **Simulation "+typeStr+"** *(truncated below lowest possible y-value)*");
		lines.add("  * Black Solid Line: "+simCalc.getSimProv().getName());
		
		int charIndex = 0;
		for (int i=0; i<compCalcs.size(); i++) {
			SimulationHazardCurveCalc<?> compCalc = compCalcs.get(i);
			PlotLineType type;
			String colorStr;
			if (customPlotChars.containsKey(compCalc)) {
				PlotCurveCharacterstics pChar = customPlotChars.get(compCalc);
				type = pChar.getLineType();
				Color color = pChar.getColor();
				if (color == Color.BLACK)
					colorStr = "Black";
				else if (color == Color.BLUE)
					colorStr = "Blue";
				else if (color == Color.CYAN)
					colorStr = "Cyan";
				else if (color == Color.DARK_GRAY)
					colorStr = "Dark Gray";
				else if (color == Color.GRAY)
					colorStr = "Gray";
				else if (color == Color.GREEN)
					colorStr = "Green";
				else if (color == Color.LIGHT_GRAY)
					colorStr = "Light Gray";
				else if (color == Color.MAGENTA)
					colorStr = "Magenta";
				else if (color == Color.ORANGE)
					colorStr = "Orange";
				else if (color == Color.PINK)
					colorStr = "Pint";
				else if (color == Color.RED)
					colorStr = "Red";
				else if (color == Color.YELLOW)
					colorStr = "Yellow";
				else
					colorStr = "("+color.getRed()+","+color.getGreen()+","+color.getBlue()+")";
			} else {
				type = simCompCurveChars[charIndex++ % simCompCurveChars.length].getLineType();
				colorStr = "Gray";
			}
			String lineType = type.name().replaceAll("_", " ");
			lineType = lineType.substring(0, 1).toUpperCase()+lineType.substring(1).toLowerCase();
			lines.add("  * "+colorStr+" "+lineType+" Line: "+compCalc.getSimProv().getName());
		}
		lines.add("* **GMPE "+typeStr+"**");
		String gmpeName = gmpeRef.getShortName();
		if (numGMPESims > 0) {
			if (numGMPESims == 1) {
				lines.add("  * Red Solid Line: "+gmpeName+", simulated (with samples from GMPE log-normal distribution)");
			} else {
				lines.add("  * Light Red Thin Solid Lines: "+gmpeName+" simulations (with samples from GMPE log-normal distribution)");
				lines.add("  * Red Solid Line: "+gmpeName+", mean of "+numGMPESims+" simulations");
			}
		}
		lines.add("  * Blue Solid Line: "+gmpeName+" full "+typeStr.toLowerCase());
		for (int i=0; hasGMPETruncs && i<gmpeTruncs.length; i++) {
			String truncAdd = optionalDigitDF.format(gmpeTruncs[i])+"-sigma truncation";
			PlotLineType type = gmpe_alt_line_types[i % gmpe_alt_line_types.length];
			String lineType = type.name().replaceAll("_", " ");
			lineType = lineType.substring(0, 1).toUpperCase()+lineType.substring(1).toLowerCase();
			lines.add("  * Blue "+lineType+" Line: "+gmpeName+", "+truncAdd);
		}
		for (int i=0; hasGMPESigmas && i<gmpeFixedSigmas.length; i++) {
			String sigmaAdd = "Fixed sigma="+optionalDigitDF.format(gmpeFixedSigmas[i]);
			PlotLineType type = gmpe_alt_line_types[i % gmpe_alt_line_types.length];
			String lineType = type.name().replaceAll("_", " ");
			lineType = lineType.substring(0, 1).toUpperCase()+lineType.substring(1).toLowerCase();
			lines.add("  * Green "+lineType+" Line: "+gmpeName+", "+sigmaAdd);
		}
		if (!spectra && hazardCurveRPs.length > 0)
			lines.add("* Gray Dashed Lines: "+Joiner.on(" yr, ").join(Ints.asList(hazardCurveRPs))+" yr return periods");
		return lines;
	}
	
	private static DiscretizedFunc[] buildEmptyFuncs(int num) {
		DiscretizedFunc[] ret = new DiscretizedFunc[num];
		for (int i=0; i<num; i++)
			ret[i] = new ArbitrarilyDiscretizedFunc();
		return ret;
	}
	
	private static void addToSpectra(DiscretizedFunc spectra, double period, DiscretizedFunc curve, double probLevel) {
		if (probLevel >= curve.getMinY()) {
			double val = curve.getFirstInterpolatedX_inLogXLogYDomain(probLevel);
			spectra.set(period, val);
		}
	}
	
	public File plotHazardSpectra(File outputDir, String prefix, double probLevel, String yAxisLabel, double[] periods) throws IOException {
		System.out.println("Plotting spectrum (and calculating all necessary curves)");
		Preconditions.checkArgument(periods.length > 1);
		
		DiscretizedFunc simSpectra = new ArbitrarilyDiscretizedFunc();
		simSpectra.setName(simCalc.getSimProv().getName());
		DiscretizedFunc gmpeSpectra = new ArbitrarilyDiscretizedFunc();
		gmpeSpectra.setName(gmpeRef.getShortName());
		
		DiscretizedFunc[] gmpeTruncSpectras = buildEmptyFuncs(gmpeTruncs.length);
		DiscretizedFunc[] gmpeFixedSigmaSpectras = buildEmptyFuncs(gmpeTruncs.length);
		
		DiscretizedFunc[] compSimSpectras = buildEmptyFuncs(compCalcs.size());
		
		List<IMT> imts = new ArrayList<>();
		for (double period : periods)
			if (period > 0)
				imts.add(IMT.forPeriod(period));
		Preconditions.checkState(!imts.isEmpty(), "No SA periods supplied");
		
		for (IMT imt : imts) {
			DiscretizedFunc gmpeCurve = getCalcGMPECurve(imt);
			addToSpectra(gmpeSpectra, imt.getPeriod(), gmpeCurve, probLevel);
			
			DiscretizedFunc simCurve = getCalcSimCurve(simCalc, imt);
			addToSpectra(simSpectra, imt.getPeriod(), simCurve, probLevel);
			
			// GMPE truncated curves
			for (int i=0; i<gmpeTruncs.length; i++) {
				DiscretizedFunc gmpeTruncCurve = getCalcGMPECurve(imt, gmpeTruncs[i], -1);
				addToSpectra(gmpeTruncSpectras[i], imt.getPeriod(), gmpeTruncCurve, probLevel);
			}
			
			// GMPE fixed sigma curves
			for (int i=0; i<gmpeFixedSigmas.length; i++) {
				DiscretizedFunc gmpeFixedSigmaCurve = getCalcGMPECurve(imt, -1, gmpeFixedSigmas[i]);
				addToSpectra(gmpeFixedSigmaSpectras[i], imt.getPeriod(), gmpeFixedSigmaCurve, probLevel);
			}
			
			for (int i=0; i<compSimSpectras.length; i++) {
				DiscretizedFunc curve = getCalcSimCurve(compCalcs.get(i), imt);
				addToSpectra(compSimSpectras[i], imt.getPeriod(), curve, probLevel);
			}
		}
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		// primary GMPE curve
		funcs.add(gmpeSpectra);
		chars.add(gmpeSpectraChar);
		
		// GMPE truncated curves
		for (int i=0; i<gmpeTruncs.length; i++) {
			String name = optionalDigitDF.format(gmpeTruncs[i])+" σₜ";
			float thickness = 2f;
			gmpeTruncSpectras[i].setName(name);
			funcs.add(gmpeTruncSpectras[i]);
			chars.add(new PlotCurveCharacterstics(gmpe_alt_line_types[i % gmpe_alt_line_types.length], thickness, Color.BLUE));
		}
		
		// GMPE fixed sigma curves
		for (int i=0; i<gmpeFixedSigmas.length; i++) {
			String name = "σ="+optionalDigitDF.format(gmpeFixedSigmas[i]);
			if (i == 0)
				name = "Fixed "+name;
			float thickness = 2f;
			gmpeFixedSigmaSpectras[i].setName(name);
			funcs.add(gmpeFixedSigmaSpectras[i]);
			chars.add(new PlotCurveCharacterstics(gmpe_alt_line_types[i % gmpe_alt_line_types.length],
					thickness, Color.GREEN.darker()));
		}
		
		// primary simulation curve
		funcs.add(simSpectra);
		chars.add(simSpectraChar);
		
		int charIndex = 0;
		for (int i=0; i<compCalcs.size(); i++) {
			SimulationHazardCurveCalc<?> compCalc = compCalcs.get(i);
			compSimSpectras[i].setName(compCalc.getSimProv().getName());
			funcs.add(compSimSpectras[i]);
			if (customPlotChars.containsKey(compCalc))
				chars.add(customPlotChars.get(compCalc));
			else
				chars.add(simCompCurveChars[charIndex++ % simCompCurveChars.length]);
		}
		
		// again on top simulation curve, no name
		DiscretizedFunc clonedSpectra = simSpectra.deepClone();
		clonedSpectra.setName(null);
		funcs.add(clonedSpectra);
		chars.add(simSpectraChar);
		
		// now plot
		String xAxisLabel = "Period (s)";

		PlotSpec spec = new PlotSpec(funcs, chars, site.getName()+" Hazard Spectra", xAxisLabel, yAxisLabel);
		spec.setLegendVisible(true);

		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(18);
		plotPrefs.setAxisLabelFontSize(20);
		plotPrefs.setPlotLabelFontSize(21);
		plotPrefs.setBackgroundColor(Color.WHITE);

		HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
		
		double minY = Double.POSITIVE_INFINITY;
		double maxY = Double.NEGATIVE_INFINITY;
		for (DiscretizedFunc func : funcs) {
			for (Point2D pt : func) {
				if (pt.getY() > 0) {
					minY = Math.min(minY, pt.getY());
					maxY = Math.max(maxY, pt.getY());
				}
			}
		}
		minY = Math.pow(10, Math.floor(Math.log10(minY)));
		maxY = Math.pow(10, Math.ceil(Math.log10(maxY)));
		
		Range xRange = new Range(StatUtils.min(periods), StatUtils.max(periods));
		Range yRange = new Range(minY, maxY);

		gp.drawGraphPanel(spec, false, true, xRange, yRange);
		gp.getChartPanel().setSize(800, 600);
		File pngFile = new File(outputDir, prefix+".png");
		gp.saveAsPNG(pngFile.getAbsolutePath());
		System.out.println("DONE "+prefix+", "+site.getName()+", "+xAxisLabel);
		return pngFile;
	}
	
	private void addContributionToSpectra(DiscretizedFunc spectra, double period, DiscretizedFunc refCurve,
			DiscretizedFunc sourceCurve, double probLevel) {
		DiscretizedFunc withoutCurve = refCurve.deepClone();
		for (int i=0; i<refCurve.size() && i<sourceCurve.size(); i++) {
			Preconditions.checkState(refCurve.getX(i) == sourceCurve.getX(i));
			double refProb = refCurve.getY(i);
			double sourceProb = sourceCurve.getY(i);
			double refRate = -Math.log(1 - refProb)/curveDuration;
			double sourceRate = -Math.log(1 - sourceProb)/curveDuration;
			double withoutRate = refRate - sourceRate;
			double withoutProb = 1 - Math.exp(-withoutRate*curveDuration);
			withoutCurve.set(i, withoutProb);
		}
		if (probLevel >= refCurve.getMinY() && probLevel >= withoutCurve.getMinY()) {
			double withVal = refCurve.getFirstInterpolatedX_inLogXLogYDomain(probLevel);
			double withoutVal = withoutCurve.getFirstInterpolatedX_inLogXLogYDomain(probLevel);
//			if (period == 3d) {
//				System.out.println("3s debug for "+spectra.getName());
//				System.out.println("\tWith val: "+withVal);
//				System.out.println("\tWithout val: "+withoutVal);
//			}
			spectra.set(period, withVal - withoutVal);
		}
	}
	
	public File plotSourceContributionHazardSpectra(File outputDir, String prefix, double[] periods,
			double probLevel, String yAxisLabel, int numSourceCurves, boolean gmpeSources) throws IOException {
		Preconditions.checkState(sourceRupContribFracts != null, "Source rup contribution fractions not set!");
		System.out.println("Plotting spectrum (and calculating all necessary curves)");
		Preconditions.checkArgument(periods.length > 1);
		
		DiscretizedFunc simSpectra = new ArbitrarilyDiscretizedFunc();
		simSpectra.setName(simCalc.getSimProv().getName());
		DiscretizedFunc gmpeSpectra = new ArbitrarilyDiscretizedFunc();
		gmpeSpectra.setName(gmpeRef.getShortName());
		
		Map<String, DiscretizedFunc> simSourceSpectra = new HashMap<>();
		Map<String, DiscretizedFunc> gmpeSourceSpectra = new HashMap<>();
		Map<String, List<Double>> fractionalContribs = new HashMap<>();
		
		List<IMT> imts = new ArrayList<>();
		for (double period : periods)
			if (period > 0)
				imts.add(IMT.forPeriod(period));
		Preconditions.checkState(!imts.isEmpty(), "No SA periods supplied");
		
		for (IMT imt : imts) {
			DiscretizedFunc gmpeCurve = getCalcGMPECurve(imt);
			addToSpectra(gmpeSpectra, imt.getPeriod(), gmpeCurve, probLevel);
			
			DiscretizedFunc simCurve = getCalcSimCurve(simCalc, imt);
			addToSpectra(simSpectra, imt.getPeriod(), simCurve, probLevel);
			
			if (probLevel < simCurve.getMinY() || probLevel > simCurve.getMaxY())
				continue;
			double contribIML = simCurve.getFirstInterpolatedX_inLogXLogYDomain(probLevel);
//			System.out.println("Sorting source contirubtions for IML="+(float)+contribIML
//					+" at p="+(float)probLevel);
			
			Map<String, DiscretizedFunc> simCurves = getCalcSimSourceCurves(imt);
			Map<String, DiscretizedFunc> gmpeCurves = getCalcGMPESourceCurves(imt);
			
			for (String sourceName : simCurves.keySet()) {
				DiscretizedFunc sourceSimCurve = simCurves.get(sourceName);
				DiscretizedFunc sourceGMPECurve = gmpeCurves.get(sourceName);
				if (sourceSimCurve == null || sourceGMPECurve == null)
					continue;
				if (!simSourceSpectra.containsKey(sourceName)) {
					simSourceSpectra.put(sourceName, new ArbitrarilyDiscretizedFunc(sourceName));
					gmpeSourceSpectra.put(sourceName, new ArbitrarilyDiscretizedFunc(sourceName));
					fractionalContribs.put(sourceName, new ArrayList<>());
				}
				addContributionToSpectra(simSourceSpectra.get(sourceName), imt.getPeriod(), simCurve, sourceSimCurve, probLevel);
				addContributionToSpectra(gmpeSourceSpectra.get(sourceName), imt.getPeriod(), gmpeCurve, sourceGMPECurve, probLevel);
//				addToSpectra(simSourceSpectra.get(sourceName), period, sourceSimCurve, probLevel);
//				addToSpectra(gmpeSourceSpectra.get(sourceName), period, sourceGMPECurve, probLevel);
				if (simSourceSpectra.get(sourceName).hasX(imt.getPeriod()))
					fractionalContribs.get(sourceName).add(simSourceSpectra.get(sourceName).getY(imt.getPeriod())/contribIML);
//					fractionalContribs.get(sourceName).add(simSourceSpectra.get(sourceName).getY(period));
				else
					fractionalContribs.get(sourceName).add(0d);
			}
		}
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		// primary GMPE curve
		if (gmpeSources) {
			funcs.add(gmpeSpectra);
			chars.add(gmpeSpectraChar);
		} else {
			// primary simulation curve
			funcs.add(simSpectra);
			chars.add(simSpectraChar);
		}
		
		Map<String, Double> simSourceSortVals = new HashMap<>();
		for (String sourceName : fractionalContribs.keySet()) {
			List<Double> vals = fractionalContribs.get(sourceName);
			simSourceSortVals.put(sourceName, StatUtils.mean(Doubles.toArray(vals)));
		}
		
		addSourceContributionFuncs(funcs, chars, simSourceSpectra, gmpeSourceSpectra, simSourceSortVals,
				numSourceCurves, sourceContribMutuallyExclusive, gmpeSources, true, IMT.SA1P0); // IMT determines x values only
		
		// now plot
		String xAxisLabel = "Period (s)";

		PlotSpec spec = new PlotSpec(funcs, chars, site.getName()+" Hazard Spectrum", xAxisLabel, yAxisLabel);
		spec.setLegendVisible(true);

		PlotPreferences plotPrefs = PlotPreferences.getDefault();
		plotPrefs.setTickLabelFontSize(18);
		plotPrefs.setAxisLabelFontSize(20);
		plotPrefs.setPlotLabelFontSize(21);
		plotPrefs.setBackgroundColor(Color.WHITE);

		HeadlessGraphPanel gp = new HeadlessGraphPanel(plotPrefs);
		
		double minY = Double.POSITIVE_INFINITY;
		double maxY = Double.NEGATIVE_INFINITY;
		for (DiscretizedFunc func : funcs) {
			for (Point2D pt : func) {
				if (pt.getY() > 0) {
					minY = Math.min(minY, pt.getY());
					maxY = Math.max(maxY, pt.getY());
				}
			}
		}
		minY = Math.pow(10, Math.floor(Math.log10(minY)));
		maxY = Math.pow(10, Math.ceil(Math.log10(maxY)));
		
		Range xRange = new Range(StatUtils.min(periods), StatUtils.max(periods));
		Range yRange = new Range(minY, maxY);

		gp.drawGraphPanel(spec, false, true, xRange, yRange);
		gp.getChartPanel().setSize(800, 600);
		File pngFile = new File(outputDir, prefix+".png");
		gp.saveAsPNG(pngFile.getAbsolutePath());
		System.out.println("DONE "+prefix+", "+site.getName()+", "+xAxisLabel);
		return pngFile;
	}
	
	static final DecimalFormat optionalDigitDF = new DecimalFormat("0.##");
	private static DecimalFormat groupedIntDF = new DecimalFormat("#");
	static {
		groupedIntDF.setGroupingUsed(true);
		groupedIntDF.setGroupingSize(3);
	}
	
}
