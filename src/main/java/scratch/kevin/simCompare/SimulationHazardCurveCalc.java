package scratch.kevin.simCompare;

import java.awt.geom.Point2D;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGV_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SignificantDurationParam;

import com.google.common.base.Preconditions;
import com.google.common.collect.Table;
import com.google.common.primitives.Doubles;

import scratch.UCERF3.erf.ETAS.ETAS_Utils;

public class SimulationHazardCurveCalc<E> {
	
	private SimulationRotDProvider<E> simProv;

	private Map<String, DiscretizedFunc> xValsMap;
	
	public static int DEFAULT_X_VAL_MULT = 4;
	
	public static DiscretizedFunc getDefaultHazardCurve(String paramName, int xValMult) {
		ArbitrarilyDiscretizedFunc xValues = new IMT_Info().getDefaultHazardCurve(paramName);
		if (xValMult > 0) {
			ArbitrarilyDiscretizedFunc newXValues = new ArbitrarilyDiscretizedFunc();
			for (int i=0; i<xValues.size()-1; i++) {
				double x0 = Math.log(xValues.getX(i));
				double x1 = Math.log(xValues.getX(i+1));
				double dx = (x1 - x0)/xValMult;
				for (int j=0; j<xValMult; j++)
					newXValues.set(Math.exp(x0 + j*dx), 1d);
			}
			newXValues.set(xValues.getMaxX(), 1d);
			xValues = newXValues;
		}
		return xValues;
	}

	public SimulationHazardCurveCalc(SimulationRotDProvider<E> simProv) {
		this(simProv, null);
	}

	public SimulationHazardCurveCalc(SimulationRotDProvider<E> simProv, Map<String, DiscretizedFunc> xValsMap) {
		this.simProv = simProv;
		if (xValsMap == null)
			xValsMap = new HashMap<>();
		if (!xValsMap.containsKey(SA_Param.NAME))
			xValsMap.put(SA_Param.NAME, getDefaultHazardCurve(SA_Param.NAME, DEFAULT_X_VAL_MULT));
		if (!xValsMap.containsKey(PGV_Param.NAME))
			xValsMap.put(PGV_Param.NAME, getDefaultHazardCurve(PGV_Param.NAME, DEFAULT_X_VAL_MULT));
		if (!xValsMap.containsKey(PGA_Param.NAME))
			xValsMap.put(PGA_Param.NAME, getDefaultHazardCurve(PGA_Param.NAME, DEFAULT_X_VAL_MULT));
		if (!xValsMap.containsKey(SignificantDurationParam.NAME))
			xValsMap.put(SignificantDurationParam.NAME, getDefaultHazardCurve(SignificantDurationParam.NAME, DEFAULT_X_VAL_MULT));
		this.xValsMap = xValsMap;
	}
	
	public SimulationRotDProvider<E> getSimProv() {
		return simProv;
	}
	
	public Map<String, DiscretizedFunc> getXValsMap() {
		return xValsMap;
	}
	
	public DiscretizedFunc calc(Site site, IMT imt, double curveDuration) throws IOException {
		return calc(site, imt, curveDuration, null, null);
	}
	
	public DiscretizedFunc calcSimDistributionFractileCurve(Site site, IMT imt,
			double curveDuration, double fractile) throws IOException {
		return calc(site, imt, curveDuration, null, fractile);
	}
	
	public Map<String, DiscretizedFunc> calcSourceContributionCurves(Site site, IMT imt, double curveDuration,
			Table<String, E, Double> sourceRupContribFracts) throws IOException {
		Map<String, DiscretizedFunc> ret = new HashMap<>();
		
		for (String sourceName : sourceRupContribFracts.rowKeySet())
			ret.put(sourceName, calc(site, imt, curveDuration, sourceRupContribFracts.row(sourceName), null));
		
		return ret;
	}
	
	private DiscretizedFunc calc(Site site, IMT imt, double curveDuration,
			Map<E, Double> rupRateScalars, Double fractile) throws IOException {
		// annual rate curve
		DiscretizedFunc curve = xValsMap.get(imt.getParamName()).deepClone();
		for (int i=0; i<curve.size(); i++)
			curve.set(i, 0d);
		int[] numExceed = new int[curve.size()];
		int numRuptures = 0;
		double firstRate = -1;
//		double minRate = Double.POSITIVE_INFINITY;
		simProv.getMinimumCurvePlotRate(site);
		boolean allRatesSame = true;
		for (E rupture : simProv.getRupturesForSite(site)) {
			double rupRate = simProv.getAnnualRate(rupture);
			if (rupRateScalars != null) {
				Double scale = rupRateScalars.get(rupture);
				if (scale == null)
					continue;
				rupRate *= scale;
			}
			if (rupRate == 0)
				continue;
			if (firstRate == -1)
				firstRate = rupRate;
			else
				allRatesSame = allRatesSame && firstRate == rupRate;
//			minRate = Math.min(rupRate, minRate);
			List<Double> vals = simProv.getValues(site, rupture, imt);
			if (fractile == null) {
				for (int j=0; j<vals.size(); j++) {
					double simRate = simProv.getIndividualSimulationRate(rupture, rupRate, j, vals.size());
					double val = vals.get(j);
					for (int i=0; i<curve.size(); i++) {
						if (curve.getX(i) <= val) {
							numExceed[i]++;
							curve.set(i, curve.getY(i)+simRate);
						}
					}
					numRuptures++;
				}
			} else {
				Preconditions.checkState(vals.size() > 1,
						"Must have multiple values per rupture for fractile curves");
				double[] array = Doubles.toArray(vals);
				double val = StatUtils.percentile(array, fractile*100d);
				for (int i=0; i<curve.size(); i++) {
					if (curve.getX(i) <= val) {
						numExceed[i]++;
						curve.set(i, curve.getY(i)+rupRate);
					}
				}
				numRuptures++;
			}
		}
		if (firstRate < 0)
			return null;
		
		DiscretizedFunc lowerCurve = null;
		DiscretizedFunc upperCurve = null;
		if (allRatesSame && rupRateScalars == null) {
			lowerCurve = curve.deepClone();
			upperCurve = curve.deepClone();
			
			double scale = firstRate*numRuptures;
			
			for (int i=0; i<curve.size(); i++) {
				double[] conf = ETAS_Utils.getBinomialProportion95confidenceInterval(
						(double)numExceed[i]/(double)numRuptures, numRuptures);
				lowerCurve.set(i, 1d - Math.exp(-conf[0]*scale*curveDuration));
				upperCurve.set(i, 1d - Math.exp(-conf[1]*scale*curveDuration));
//				System.out.println("x="+(float)xVals.getX(i)+"\ty="+(float)curve.getY(i)
//					+"\tc[0]="+(float)conf[0]+"\tc[0]*s="+(float)(conf[0]*scale)
//					+"\tc[1]="+(float)conf[1]+"\tc[1]*s="+(float)(conf[1]*scale));
			}
		}
		
		// now probabilities
		for (int i=0; i<curve.size(); i++) {
			double rate = curve.getY(i);
			double prob = 1d - Math.exp(-rate*curveDuration);
			curve.set(i, prob);
		}
		
		double minPlotRate = simProv.getMinimumCurvePlotRate(site);
		if (minPlotRate > 0 && Double.isFinite(minPlotRate)) {
//		minRate = Math.min(minRate, );
//		if (minRate > 0) {
			double minProb = 1d - Math.exp(-minPlotRate*curveDuration);
			// truncate curve to remove x values never seen in finite catalog
			ArbitrarilyDiscretizedFunc truncatedCurve = new ArbitrarilyDiscretizedFunc();
			ArbitrarilyDiscretizedFunc truncatedLowerCurve = null;
			ArbitrarilyDiscretizedFunc truncatedUpperCurve = null;
			if (lowerCurve != null) {
				truncatedLowerCurve = new ArbitrarilyDiscretizedFunc();
				truncatedUpperCurve = new ArbitrarilyDiscretizedFunc();
			}
			for (int i=0; i<curve.size(); i++) {
				Point2D pt = curve.get(i);
				if (pt.getY() >= minProb) {
					truncatedCurve.set(pt);
					if (truncatedLowerCurve != null) {
						truncatedLowerCurve.set(lowerCurve.get(i));
						truncatedUpperCurve.set(upperCurve.get(i));
					}
				}
			}
			curve = truncatedCurve;
			lowerCurve = truncatedLowerCurve;
			upperCurve = truncatedUpperCurve;
		}
		
		if (lowerCurve != null)
			curve = new UncertainArbDiscFunc(curve, lowerCurve, upperCurve);
		
		curve.setName(simProv.getName());
		
		return curve;
	}

}
