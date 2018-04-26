package scratch.aftershockStatisticsETAS.griddedInterpGMPE;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.imr.ScalarIMR;

public class IntensityMeasureLevelInterpolator extends AbstractGMPEInterpolation.Discrete<Double> {
	
	private DiscretizedFunc xValsFunc;

	IntensityMeasureLevelInterpolator(String name, DiscretizedFunc xValsFunc, boolean convertToLog) {
		super(name, toList(xValsFunc, convertToLog));
		this.xValsFunc = xValsFunc;
	}
	
	private static List<Double> toList(DiscretizedFunc xValsFunc, boolean convertToLog) {
		List<Double> values = new ArrayList<>();
		for (Point2D pt : xValsFunc) {
			if (convertToLog)
				values.add(Math.log(pt.getX()));
			else
				values.add(pt.getX());
		}
		return values;
	}
	
	/**
	 * 
	 * @return x-values function in linear space
	 */
	public DiscretizedFunc getXValsFunc() {
		return xValsFunc;
	}

	@Override
	public void setGMPE_Params(ScalarIMR gmpe, ProbEqkSource source, int index) {
		gmpe.setIntensityMeasureLevel(getValue(index));
	}

//	@Override
//	public Set<EqkRupture> getViableRuptures(Set<EqkRupture> ruptures, int index) {
//		return ruptures;
//	}

	@Override
	public Double detectCurrentVal(ScalarIMR gmpe, Site site) {
		return (Double)gmpe.getIntensityMeasureLevel();
	}

}
