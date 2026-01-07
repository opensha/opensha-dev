package scratch.kevin.tdProbModelPlayground;

import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;

public class WG02_ProbabilityModel implements FSS_ERF_ProbabilityModel {

	@Override
	public double getProbability(FaultSystemSolution fltSysSolution, int fltSysRupIndex, long forecastStartTimeMillis,
			double durationYears) {
		throw new UnsupportedOperationException("Not yet implemented");
	}

	@Override
	public double getProbabilityGain(FaultSystemSolution fltSysSolution, int fltSysRupIndex, long forecastStartTimeMillis,
			double durationYears) {
		throw new UnsupportedOperationException("Not yet implemented");
	}

}
