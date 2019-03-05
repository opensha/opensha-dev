package scratch.kevin.bbp;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.opensha.commons.util.ExceptionUtils;

import com.google.common.base.Preconditions;

import scratch.kevin.bbp.BBP_Module.Method;
import scratch.kevin.bbp.BBP_Module.VelocityModel;

public class BBP_ThreadedRunner {

	public static void main(String[] args) {
		File eventDir = new File("/data/kevin/simulators/catalogs/rundir2194_long/event_srfs");
		File srcFile = new File(eventDir, "event_136704.src");
		File sitesFile = new File("/home/kevin/bbp/bbp_data/run/stations_cs_sites.stl");
		VelocityModel vm = VelocityModel.LA_BASIN_863;
		Method method = Method.GP;
		
		int numRealizations = 100;
		
		File mainOutputDir = new File(eventDir, "event_136704_gp_rg");
		Preconditions.checkState(mainOutputDir.exists() || mainOutputDir.mkdir());
		ExecutorService exec = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
		
		List<Future<?>> futures = new ArrayList<>();
		
		Random r = new Random();
		
		for (int i=0; i<numRealizations; i++) {
			int index = i;
			File outputDir = new File(mainOutputDir, "run_"+i);
			Runnable run = new Runnable() {
				
				@Override
				public void run() {
					System.out.println("Executiong "+index+"/"+numRealizations);
					BBP_Wrapper wrapper = new BBP_Wrapper(vm, method, srcFile,
							(long)r.nextInt(Short.MAX_VALUE), null, sitesFile, outputDir);
					wrapper.run();
					System.out.println("DONE "+index+"/"+numRealizations);
				}
			};
			futures.add(exec.submit(run));
			try {
				Thread.sleep(5000);
			} catch (InterruptedException e) {}
		}
		
		for (Future<?> f : futures) {
			try {
				f.get();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		exec.shutdown();
		System.exit(0);
	}

}
