package scratch.kevin.simCompare;

import java.awt.geom.Point2D;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;

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
	
	public DiscretizedFunc getXVals() {
		return xVals;
	}
	
	public DiscretizedFunc calc(Site site, double period, double curveDuration) throws IOException {
		// annual rate curve
		DiscretizedFunc curve = xVals.deepClone();
		for (int i=0; i<curve.size(); i++)
			curve.set(i, 0d);
		for (E rupture : simProv.getRupturesForSite(site)) {
			double rupRate = simProv.getAnnualRate(rupture);
			List<DiscretizedFunc> spectras = simProv.getRotD50s(site, rupture);
			rupRate /= spectras.size();
			for (DiscretizedFunc spectra : spectras) {
				double rd50 = spectra.getInterpolatedY(period);
				for (int i=0; i<curve.size(); i++)
					if (curve.getX(i) <= rd50)
						curve.set(i, curve.getY(i)+rupRate);
			}
		}
		
		// now probabilities
		for (int i=0; i<curve.size(); i++) {
			double rate = curve.getY(i);
			double prob = 1d - Math.exp(-rate*curveDuration);
			curve.set(i, prob);
		}
		
		double minRate = simProv.getMinimumCurvePlotRate();
		if (minRate > 0) {
			double minProb = 1d - Math.exp(-minRate*curveDuration);
			// truncate curve to remove x values never seen in finite catalog
			ArbitrarilyDiscretizedFunc truncatedCurve = new ArbitrarilyDiscretizedFunc();
			for (Point2D pt : curve)
				if (pt.getY() >= minProb)
					truncatedCurve.set(pt);
			curve = truncatedCurve;
		}
		
		curve.setName(simProv.getName());
		
		return curve;
	}

}
