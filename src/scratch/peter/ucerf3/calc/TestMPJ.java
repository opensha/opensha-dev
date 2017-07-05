package scratch.peter.ucerf3.calc;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.opensha.commons.hpc.mpj.taskDispatch.MPJTaskCalculator;

/**
 * Tests functioning of MPJ lib.
 * 
 * @author Peter Powers
 * @version $Id:$
 */
public class TestMPJ extends MPJTaskCalculator {

	public TestMPJ(CommandLine cmd, String[] args)
			throws IOException, InvocationTargetException, FileNotFoundException {
		super(cmd);
		System.out.println("Initialized");
	}

	@Override
	protected int getNumTasks() {
		return 100;
	}

	@Override
	protected void calculateBatch(int[] batch) throws Exception {
		for (int i : batch) {
			System.out.println("Completed job " + i);
		}
	}

	@Override
	protected void doFinalAssembly() throws Exception {}

	public static void main(String[] args) {
		args = MPJTaskCalculator.initMPJ(args);

		try {
			Options options = createOptions();
			CommandLine cmd = parse(options, args,
				TestMPJ.class);
			args = cmd.getArgs();
			TestMPJ driver = new TestMPJ(cmd,
				args);
			driver.run();
			finalizeMPJ();
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}

}
