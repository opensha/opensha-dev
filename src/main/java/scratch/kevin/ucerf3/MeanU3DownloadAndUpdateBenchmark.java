package scratch.kevin.ucerf3;

import java.util.concurrent.TimeUnit;

import com.google.common.base.Stopwatch;

import scratch.UCERF3.erf.mean.MeanUCERF3;

public class MeanU3DownloadAndUpdateBenchmark {

	public static void main(String[] args) {
		Stopwatch watch = Stopwatch.createStarted();
		MeanUCERF3 u3 = new MeanUCERF3();
		
		u3.updateForecast();
		
		watch.stop();
		System.out.println("Took "+watch.elapsed(TimeUnit.SECONDS)+"s");
	}

}
