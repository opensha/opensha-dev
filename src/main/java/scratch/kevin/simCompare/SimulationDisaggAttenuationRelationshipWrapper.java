package scratch.kevin.simCompare;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.exceptions.IMRException;
import org.opensha.commons.exceptions.ParameterException;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.constraint.impl.DoubleDiscreteConstraint;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.faultSurface.cache.SurfaceDistances;
import org.opensha.sha.imr.AttenuationRelationship;
import org.opensha.sha.imr.param.IntensityMeasureParams.DampingParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGV_Param;
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
	
	private IMT[] imts;

	public SimulationDisaggAttenuationRelationshipWrapper(SimulationRotDProvider<E> simProv,
			Source probSource, Source meanSource, Source stdDevSource, IMT[] imts) {
		this.simProv = simProv;
		this.probSource = probSource;
		this.meanSource = meanSource;
		this.stdDevSource = stdDevSource;
		this.imts = imts;
		
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
	private IMT prevIMT;
	
	private double[] prevSimValues;
	
	/**
	 * @param rupture
	 * @param site
	 * @param imt
	 * @return simulation values for the given rupture in log space, using a single valued cache
	 */
	private synchronized double[] getSimValues(E rupture, Site site, IMT imt) {
		if (prevRup != rupture || site != prevSite || prevIMT != imt)
			prevSimValues = null;
		
		if (prevSimValues == null) {
			List<Double> simVals;
			try {
				simVals = simProv.getValues(site, rupture, imt);
			} catch (IOException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
			Preconditions.checkState(simVals.size() > 0);
			prevSimValues = new double[simVals.size()];
			for (int i=0; i<prevSimValues.length; i++)
				prevSimValues[i] = Math.log(simVals.get(i));
		}
		
		return prevSimValues;
	}

	@Override
	public double getMean() {
		EqkRupture rup = getEqkRupture();
		Preconditions.checkState(rup instanceof CompProbEqkRupture);
		RuptureComparison<E> comp = ((RuptureComparisonERF<E>.CompProbEqkRupture)rup).getRuptureComparison();
		Site site = getSite();
		IMT imt = getIMT();
		switch (meanSource) {
		case SIMULATION:
			double[] simVals = getSimValues(comp.getRupture(), site, imt);
			return StatUtils.mean(simVals);
		case GMPE:
			return comp.getLogMean(site, imt);

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
		IMT imt = getIMT();
		switch (stdDevSource) {
		case SIMULATION:
			double[] simVals = getSimValues(comp.getRupture(), site, imt);
			return Math.sqrt(StatUtils.variance(simVals));
		case GMPE:
			return comp.getStdDev(site, imt);

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
		IMT imt = getIMT();
		double iml = (Double)getIntensityMeasureLevel();
		switch (probSource) {
		case SIMULATION:
			int numAbove = 0;
			double[] simVals = getSimValues(comp.getRupture(), site, imt);
			for (double simIM : simVals)
				if (simIM > iml)
					numAbove++;
			return (double)numAbove/(double)simVals.length;
		case GMPE:
			double gmpeMean = comp.getLogMean(site, imt);
			double gmpeStdDev = comp.getStdDev(site, imt);
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
	public void setPropagationEffectParams(SurfaceDistances distances) {}

	private IMT getIMT() {
		Parameter<?> imtParam = getIntensityMeasure();
		if (imtParam.getName().equals(SA_Param.NAME))
			return IMT.forPeriod(SA_Param.getPeriodInSA_Param(imtParam));
		else if (imtParam.getName().equals(PGV_Param.NAME))
			return IMT.PGV;
		else if (imtParam.getName().equals(PGA_Param.NAME))
			return IMT.PGA;
		else
			throw new IllegalStateException("IMT not yet supported: "+imtParam.getName());
	}

	@Override
	protected void initSupportedIntensityMeasureParams() {
		List<Double> periodsList = new ArrayList<>();
		for (IMT imt : imts) {
			if (imt.getParamName().equals(SA_Param.NAME)) {
				periodsList.add(imt.getPeriod());
			} else if (imt == IMT.PGV) {
				pgvParam = new PGV_Param();
				supportedIMParams.addParameter(pgvParam);
				
				setIntensityMeasure(pgvParam);
			} else if (imt == IMT.PGA) {
				pgaParam = new PGA_Param();
				supportedIMParams.addParameter(pgaParam);
				
				setIntensityMeasure(pgaParam);
			} else {
				throw new IllegalStateException("IMT not yet supported: "+imt);
			}
		}
		
		if (!periodsList.isEmpty()) {
			saDampingParam = new DampingParam();
			saDampingParam.setNonEditable();
			
			DoubleDiscreteConstraint periodConstraint = new DoubleDiscreteConstraint(periodsList);
			
			saPeriodParam = new PeriodParam(periodConstraint, periodsList.get(0), false);
			saPeriodParam.setValueAsDefault();
			
			supportedIMParams.clear();
			saParam = new SA_Param(saPeriodParam, saDampingParam);
			saParam.setNonEditable();
			supportedIMParams.addParameter(saParam);
			
			setIntensityMeasure(saParam);
		}
	}

	@Override
	protected void initSiteParams() {}

	@Override
	protected void initEqkRuptureParams() {}

	@Override
	protected void initPropagationEffectParams() {}

}
