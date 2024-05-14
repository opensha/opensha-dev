package scratch.kevin.mpj;

import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;

public class MPJ_EnvWrapperTest extends MPJTaskCalculator {

	public MPJ_EnvWrapperTest(CommandLine cmd) {
		super(cmd);
		
		// see if we have the env var
		Map<String, String> env = System.getenv();
		debug("$TEST_ENV: "+env.get("TEST_ENV"));
	}

	@Override
	protected int getNumTasks() {
		return size;
	}

	@Override
	protected void calculateBatch(int[] batch) throws Exception {
		for (int index : batch)
			debug("Calc "+index);
	}

	@Override
	protected void doFinalAssembly() throws Exception {}

	public static void main(String[] args) {
		System.setProperty("java.awt.headless", "true");
		try {
			args = MPJTaskCalculator.initMPJ(args);
			
			Options options = createOptions();
			
			CommandLine cmd = parse(options, args, MPJ_EnvWrapperTest.class);
			
			MPJ_EnvWrapperTest driver = new MPJ_EnvWrapperTest(cmd);
			driver.run();
			
			finalizeMPJ();
			
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}

}
