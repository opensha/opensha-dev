package scratch.kevin.tdProbModelPlayground;

import java.util.EnumSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;

import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.param.ParameterList;
import org.opensha.commons.param.event.ParameterChangeEvent;
import org.opensha.commons.param.event.ParameterChangeListener;
import org.opensha.commons.param.impl.EnumParameter;
import org.opensha.sha.earthquake.calc.recurInterval.EqkProbDistCalc;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;

import com.google.common.base.Preconditions;

public class UCERF3_ProbabilityModel implements FSS_ERF_ProbabilityModel, ParameterChangeListener {
	
	private AperiodicityModels aperiodicityModel;
	private RenewalDistributions renewalDist;
	private boolean aveRecurIntervals;
	private boolean aveNormTimeSinceLast;
	
	// implementation note: if we ever add AperiodicityModels that return continuous values, this cache will explode and
	// will need to be rethought
	private ConcurrentMap<Double, EvenlyDiscretizedFunc> normCDFsCache;
	
	static double max_time_for_normBPT_CDF=5;
	static int num_for_normBPT_CDF=501;
	
	public static final String APERIODICITY_PARAM_NAME = "Aperiodicity Model";
	public static final String RENEWAL_MODEL_PARAM_NAME = "Renewal Model";
	
	private EnumParameter<AperiodicityModels> aperiodicityParam;
	private EnumParameter<RenewalDistributions> renewalDistParam;
	
	private ParameterList params;
	
	public UCERF3_ProbabilityModel(
			AperiodicityModels aperiodicityModel, EnumSet<AperiodicityModels> supportedAperiodicityModels,
			RenewalDistributions renewalDist, EnumSet<RenewalDistributions> supportedRenewalDists) {
		this(aperiodicityModel, supportedAperiodicityModels, renewalDist, supportedRenewalDists, false, true);
	}
	
	public UCERF3_ProbabilityModel(AperiodicityModels aperiodicityModel, EnumSet<AperiodicityModels> supportedAperiodicityModels,
			RenewalDistributions renewalDist, EnumSet<RenewalDistributions> supportedRenewalDists,
			boolean aveRecurIntervals, boolean aveNormTimeSinceLast) {
		params = new ParameterList();
		
		this.aperiodicityModel = aperiodicityModel;
		aperiodicityParam = initModelParam(APERIODICITY_PARAM_NAME, aperiodicityModel, supportedAperiodicityModels);
		aperiodicityParam.addParameterChangeListener(this);
		params.addParameter(aperiodicityParam);
		this.renewalDist = renewalDist;
		renewalDistParam = initModelParam(RENEWAL_MODEL_PARAM_NAME, renewalDist, supportedRenewalDists);
		renewalDistParam.addParameterChangeListener(this);
		params.addParameter(renewalDistParam);
		// TODO: add these to params
		this.aveRecurIntervals = aveRecurIntervals;
		this.aveNormTimeSinceLast = aveNormTimeSinceLast;
		
		normCDFsCache = new ConcurrentHashMap<>();
	}
	
	private static <E extends Enum<E>> EnumParameter<E> initModelParam(String name, E defaultValue, EnumSet<E> options) {
		Preconditions.checkNotNull(defaultValue, "Default value is null");
		if (options == null)
			options = EnumSet.of(defaultValue);
		else
			Preconditions.checkState(options.contains(defaultValue),
					"Allowed set doesn't contain default value: "+defaultValue);
		return new EnumParameter<E>(name, options, defaultValue, null);
	}
	
	private EvenlyDiscretizedFunc getNormCDF(double aperiodicity) {
		if (!normCDFsCache.containsKey(aperiodicity)) {
			EqkProbDistCalc distCalc = renewalDist.instance();
			double delta = max_time_for_normBPT_CDF/(num_for_normBPT_CDF-1);
			distCalc.setAll(1.0, aperiodicity, delta, num_for_normBPT_CDF);
			EvenlyDiscretizedFunc normCDF = distCalc.getCDF();
			// putIfAbsent for thread-safety
			normCDFsCache.putIfAbsent(aperiodicity, normCDF);
		}
		return normCDFsCache.get(aperiodicity);
	}

	@Override
	public double getProbability(FaultSystemSolution fltSysSolution, int fltSysRupIndex, long forecastStartTimeMillis,
			double durationYears) {
		double aperiodicity = aperiodicityModel.getAperiodicity(fltSysSolution, fltSysRupIndex);
		
		EvenlyDiscretizedFunc normCDF = getNormCDF(aperiodicity);
		
		// TODO: do the actual calculation
		
		return Double.NaN;
	}

	@Override
	public double getProbabilityGain(FaultSystemSolution fltSysSolution, int fltSysRupIndex, long forecastStartTimeMillis,
			double durationYears) {
		double aperiodicity = aperiodicityModel.getAperiodicity(fltSysSolution, fltSysRupIndex);
		
		EvenlyDiscretizedFunc normCDF = getNormCDF(aperiodicity);
		
		// TODO: do the actual calculation
		
		return Double.NaN;
	}

	@Override
	public void parameterChange(ParameterChangeEvent event) {
		if (event.getParameter() == aperiodicityParam) {
			aperiodicityModel = aperiodicityParam.getValue();
		} else if (event.getParameter() == renewalDistParam) {
			normCDFsCache.clear();
			renewalDist = renewalDistParam.getValue();
		}
	}

	@Override
	public ParameterList getAdjustableParameters() {
		return params;
	}

}
