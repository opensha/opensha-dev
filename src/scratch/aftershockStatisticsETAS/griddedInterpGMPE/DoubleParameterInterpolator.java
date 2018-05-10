package scratch.aftershockStatisticsETAS.griddedInterpGMPE;

import java.util.Set;

import org.opensha.commons.data.Site;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.impl.WarningDoubleParameter;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.imr.ScalarIMR;

import scratch.aftershockStatisticsETAS.griddedInterpGMPE.AbstractGMPEInterpolation.EvenlySpacedDouble;

public class DoubleParameterInterpolator extends EvenlySpacedDouble {

	public DoubleParameterInterpolator(String paramName, double min, double max, int numBins) {
		super(paramName, min, max, numBins);
	}

	@Override
	public void setGMPE_Params(ScalarIMR gmpe, ProbEqkSource source, int index) {
		Parameter<Double> param = gmpe.getParameter(getName());
		if (param instanceof WarningDoubleParameter)
			((WarningDoubleParameter)param).setValueIgnoreWarning(getValue(index));
		else
			param.setValue(getValue(index));
	}

//	@Override
//	public Set<EqkRupture> getViableRuptures(Set<EqkRupture> ruptures, int index) {
//		return ruptures;
//	}

	@Override
	public Double detectCurrentVal(ScalarIMR gmpe, Site site) {
		if (site.containsParameter(getName()))
			return (Double)site.getParameter(Double.class, getName()).getValue();
		return (Double)gmpe.getParameter(getName()).getValue();
	}

}
