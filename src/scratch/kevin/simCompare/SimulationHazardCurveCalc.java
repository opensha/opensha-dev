package scratch.kevin.simCompare;

import java.awt.geom.Point2D;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.UncertainArbDiscDataset;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;

import com.google.common.collect.Table;

import scratch.UCERF3.erf.ETAS.ETAS_Utils;

public class SimulationHazardCurveCalc<E> {
	
	private SimulationRotDProvider<E> simProv;
	
	private DiscretizedFunc xVals;
	
	private static DiscretizedFunc getDefaultHazardCurve(int xValMult) {
		ArbitrarilyDiscretizedFunc xValues = new IMT_Info().getDefaultHazardCurve(SA_Param.NAME);
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
		this(simProv, getDefaultHazardCurve(4));
	}

	public SimulationHazardCurveCalc(SimulationRotDProvider<E> simProv, DiscretizedFunc xVals) {
		this.simProv = simProv;
		this.xVals = xVals;
	}
	
	public SimulationRotDProvider<E> getSimProv() {
		return simProv;
	}
	
	public DiscretizedFunc getXVals() {
		return xVals;
	}
	
	public DiscretizedFunc calc(Site site, double period, double curveDuration) throws IOException {
		return calc(site, period, curveDuration, null);
	}
	
	public Map<String, DiscretizedFunc> calcSourceContributionCurves(Site site, double period, double curveDuration,
			Table<String, E, Double> sourceRupContribFracts) throws IOException {
		Map<String, DiscretizedFunc> ret = new HashMap<>();
		
		for (String sourceName : sourceRupContribFracts.rowKeySet())
			ret.put(sourceName, calc(site, period, curveDuration, sourceRupContribFracts.row(sourceName)));
		
		return ret;
	}
	
	private DiscretizedFunc calc(Site site, double period, double curveDuration, Map<E, Double> rupRateScalars)
			throws IOException {
		// annual rate curve
		DiscretizedFunc curve = xVals.deepClone();
		for (int i=0; i<curve.size(); i++)
			curve.set(i, 0d);
		int[] numExceed = new int[xVals.size()];
		int numRuptures = 0;
		double firstRate = -1;
		double minRate = Double.POSITIVE_INFINITY;
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
			minRate = Math.min(rupRate, minRate);
			List<DiscretizedFunc> spectras = simProv.getRotD50s(site, rupture);
			rupRate /= spectras.size();
			for (DiscretizedFunc spectra : spectras) {
				double rd50 = spectra.getInterpolatedY(period);
				for (int i=0; i<curve.size(); i++) {
					if (curve.getX(i) <= rd50) {
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
			lowerCurve = xVals.deepClone();
			upperCurve = xVals.deepClone();
			
			double scale = firstRate*numRuptures;
			
			for (int i=0; i<xVals.size(); i++) {
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
		
		minRate = Math.min(minRate, simProv.getMinimumCurvePlotRate());
		if (minRate > 0) {
			double minProb = 1d - Math.exp(-minRate*curveDuration);
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
			curve = new UncertainArbDiscDataset(curve, lowerCurve, upperCurve);
		
		curve.setName(simProv.getName());
		
		return curve;
	}

}
