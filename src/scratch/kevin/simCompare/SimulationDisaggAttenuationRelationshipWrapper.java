package scratch.kevin.simCompare;

import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.exceptions.IMRException;
import org.opensha.commons.exceptions.ParameterException;
import org.opensha.commons.param.constraint.impl.DoubleDiscreteConstraint;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.imr.AttenuationRelationship;
import org.opensha.sha.imr.param.IntensityMeasureParams.DampingParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.PeriodParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

import scratch.kevin.simCompare.RuptureComparisonERF.CompProbEqkRupture;

public class SimulationDisaggAttenuationRelationshipWrapper<E> extends AttenuationRelationship {
	
	public enum Source {
		SIMULATION,
		GMPE
	}
	
	private SimulationRotDProvider<E> simProv;
	
	private Source probSource;
	private Source meanSource;
	private Source stdDevSource;
	
	private double[] periods;

	public SimulationDisaggAttenuationRelationshipWrapper(SimulationRotDProvider<E> simProv,
			Source probSource, Source meanSource, Source stdDevSource, double[] periods) {
		this.simProv = simProv;
		this.probSource = probSource;
		this.meanSource = meanSource;
		this.stdDevSource = stdDevSource;
		this.periods = periods;
		
		initSupportedIntensityMeasureParams();
	}

	public Source getProbSource() {
		return probSource;
	}

	public void setProbSource(Source probSource) {
		this.probSource = probSource;
	}

	public Source getMeanSource() {
		return meanSource;
	}

	public void setMeanSource(Source meanSource) {
		this.meanSource = meanSource;
	}

	public Source getStdDevSource() {
		return stdDevSource;
	}

	public void setStdDevSource(Source stdDevSource) {
		this.stdDevSource = stdDevSource;
	}
	
	private E prevRup;
	private Site prevSite;
	private double prevPeriod;
	
	private double[] prevSimValues;
	
	/**
	 * @param rupture
	 * @param site
	 * @param period
	 * @return simulation values for the given rupture in log space, using a single valued cache
	 */
	private synchronized double[] getSimValues(E rupture, Site site, double period) {
		if (prevRup != rupture || site != prevSite || prevPeriod != period)
			prevSimValues = null;
		
		if (prevSimValues == null) {
			List<DiscretizedFunc> simRotDs;
			try {
				simRotDs = simProv.getRotD50s(site, rupture);
			} catch (IOException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
			Preconditions.checkState(simRotDs.size() > 0);
			prevSimValues = new double[simRotDs.size()];
			for (int i=0; i<prevSimValues.length; i++)
				prevSimValues[i] = Math.log(simRotDs.get(i).getInterpolatedY(period));
		}
		
		return prevSimValues;
	}

	@Override
	public double getMean() {
		EqkRupture rup = getEqkRupture();
		Preconditions.checkState(rup instanceof CompProbEqkRupture);
		RuptureComparison<E> comp = ((RuptureComparisonERF<E>.CompProbEqkRupture)rup).getRuptureComparison();
		Site site = getSite();
		double period = saPeriodParam.getValue();
		switch (meanSource) {
		case SIMULATION:
			double[] simVals = getSimValues(comp.getRupture(), site, period);
			return StatUtils.mean(simVals);
		case GMPE:
			return comp.getLogMean(site, period);

		default:
			throw new IllegalStateException("Unknown source: "+meanSource);
		}
	}

	@Override
	public double getStdDev() {
		EqkRupture rup = getEqkRupture();
		Preconditions.checkState(rup instanceof CompProbEqkRupture);
		RuptureComparison<E> comp = ((RuptureComparisonERF<E>.CompProbEqkRupture)rup).getRuptureComparison();
		Site site = getSite();
		double period = saPeriodParam.getValue();
		switch (stdDevSource) {
		case SIMULATION:
			double[] simVals = getSimValues(comp.getRupture(), site, period);
			return Math.sqrt(StatUtils.variance(simVals));
		case GMPE:
			return comp.getStdDev(site, period);

		default:
			throw new IllegalStateException("Unknown source: "+meanSource);
		}
	}

//	@Override
//	public double getEpsilon() {
//		double iml = ((Double) im.getValue()).doubleValue();
//		return (iml - getMean())/getStdDev();
//	}

	@Override
	public double getExceedProbability() throws ParameterException, IMRException {
		EqkRupture rup = getEqkRupture();
		Preconditions.checkState(rup instanceof CompProbEqkRupture);
		RuptureComparison<E> comp = ((RuptureComparisonERF<E>.CompProbEqkRupture)rup).getRuptureComparison();
		Site site = getSite();
		double period = saPeriodParam.getValue();
		double iml = (Double)getIntensityMeasureLevel();
		switch (probSource) {
		case SIMULATION:
			int numAbove = 0;
			double[] simVals = getSimValues(comp.getRupture(), site, period);
			for (double simIM : simVals)
				if (simIM > iml)
					numAbove++;
			return (double)numAbove/(double)simVals.length;
		case GMPE:
			double gmpeMean = comp.getLogMean(site, period);
			double gmpeStdDev = comp.getStdDev(site, period);
			return AttenuationRelationship.getExceedProbability(gmpeMean, gmpeStdDev, iml, null, null);

		default:
			throw new IllegalStateException("Unknown source: "+meanSource);
		}
	}

	@Override
	public void setParamDefaults() {}
	
	/**
	 *  Returns name of the IntensityMeasureRelationship.
	 *
	 * @return    The name string
	 */
	public String getName() {
		return "Simulation Disagg Wrapper: "+simProv.getName();
	}

	@Override
	public String getShortName() {
		return "SimDisaggWrap";
	}

	@Override
	protected void setPropagationEffectParams() {}

	@Override
	protected void initSupportedIntensityMeasureParams() {
		saDampingParam = new DampingParam();
		saDampingParam.setNonEditable();
		
		DoubleDiscreteConstraint periodList = new DoubleDiscreteConstraint(Doubles.asList(periods));
		
		saPeriodParam = new PeriodParam(periodList, periods[0], false);
		saPeriodParam.setValueAsDefault();
		
		supportedIMParams.clear();
		saParam = new SA_Param(saPeriodParam, saDampingParam);
		saParam.setNonEditable();
		supportedIMParams.addParameter(saParam);
		
		setIntensityMeasure(saParam);
	}

	@Override
	protected void initSiteParams() {}

	@Override
	protected void initEqkRuptureParams() {}

	@Override
	protected void initPropagationEffectParams() {}

}
