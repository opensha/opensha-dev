package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;

import org.opensha.sha.earthquake.AbstractNthRupERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;

import scratch.UCERF3.erf.ETAS.launcher.ETAS_Config;
import scratch.UCERF3.erf.ETAS.launcher.ETAS_Launcher;

public class GridSeisCorrRateCalcs {

	public static void main(String[] args) throws IOException {
		ETAS_Config conf = ETAS_Config.readJSON(new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2025_09_15-Start1930_100yr_kCOV1p5_Spontaneous_HistCatalog/config.json"));
		
		conf.setGridSeisCorr(true);
		
		ETAS_Launcher launcher = new ETAS_Launcher(conf, false);
		AbstractNthRupERF erf = launcher.checkOutERF();
//		erf.getTimeSpan().setDuration(1d);
//		erf.updateForecast();
		double duration = erf.getTimeSpan().getDuration();
		
		double m5RateWith = 0d;
		for (ProbEqkSource source : erf)
			for (ProbEqkRupture rup : source)
				if (rup.getMag() >= 5d)
					m5RateWith += rup.getMeanAnnualRate(duration);
		
		conf.setGridSeisCorr(false);
		
		launcher = new ETAS_Launcher(conf, false);
		erf = launcher.checkOutERF();
//		erf.getTimeSpan().setDuration(1d);
//		erf.updateForecast();
		duration = erf.getTimeSpan().getDuration();
		
		double m5RateWithout = 0d;
		for (ProbEqkSource source : erf)
			for (ProbEqkRupture rup : source)
				if (rup.getMag() >= 5d)
					m5RateWithout += rup.getMeanAnnualRate(duration);
		
		System.out.println("M5 rate with corr: "+(float)m5RateWith);
		System.out.println("M5 rate without corr: "+(float)m5RateWithout);
		System.out.println("M5 rate gain from corr: "+(float)(m5RateWith/m5RateWithout));
	}

}
