package scratch.aftershockStatisticsETAS;

import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;

public class ETAS_forecast {

	public ETAS_forecast(GenericETAS_Parameters genericETAS_Parameters,
			ObsEqkRupture mainshock, ObsEqkRupList aftershocks,
			double startTime, double endTime, double maxMag, int maxGenerations) {
			
		this(genericETAS_Parameters.get_aValueMean(), genericETAS_Parameters.get_bValue(), 
				genericETAS_Parameters.get_pValue(), genericETAS_Parameters.get_cValue(),
				genericETAS_Parameters.get_alpha(), genericETAS_Parameters.get_refMag(),
				mainshock, aftershocks, startTime, endTime, maxMag, maxGenerations);
	}

		
	public ETAS_forecast(double a, double b, double p, double c, double alpha, double refMag, 
			ObsEqkRupture mainshock, ObsEqkRupList aftershocks,
			double forecastStart, double forecastEnd, double maxMag, int maxGenerations) {
	}	
}
