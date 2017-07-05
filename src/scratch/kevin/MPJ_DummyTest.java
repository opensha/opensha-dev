package scratch.kevin;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.opensha.commons.hpc.mpj.taskDispatch.MPJTaskCalculator;
import org.opensha.commons.util.ExceptionUtils;

public class MPJ_DummyTest extends MPJTaskCalculator {
	
	private int num = 1000;
	private long sleepDuration = 1000;
	private FileWriter fw;

	public MPJ_DummyTest(CommandLine cmd) {
		super(cmd);
		
		if (cmd.hasOption("num"))
			num = Integer.parseInt(cmd.getOptionValue("num"));
		
		if (cmd.hasOption("time"))
			sleepDuration = (long)(1000d*Double.parseDouble(cmd.getOptionValue("time")));
		
		if (cmd.hasOption("file") && rank == 0) {
			File file = new File(cmd.getOptionValue("file"));
			try {
				fw = new FileWriter(file);
				fw.write(getDebugText("Starting")+"\n");
				fw.flush();
			} catch (IOException e) {
				ExceptionUtils.throwAsRuntimeException(e);
			}
		}
	}

	@Override
	protected int getNumTasks() {
		return num;
	}

	@Override
	protected void calculateBatch(int[] batch) throws Exception {
		for (int i=0; i<batch.length; i++)
			Thread.sleep(sleepDuration);
	}

	@Override
	protected void doFinalAssembly() throws Exception {
		if (fw != null) {
			try {
				fw.write(getDebugText("DONE")+"\n");
				fw.close();
			} catch (IOException e) {
				ExceptionUtils.throwAsRuntimeException(e);
			}
		}
	}
	
	public static Options createOptions() {
		Options ops = MPJTaskCalculator.createOptions();
		
		Option numOp = new Option("n", "num", true, "Number of tasks");
		numOp.setRequired(false);
		ops.addOption(numOp);
		
		Option timeOp = new Option("ts", "time", true, "Time per task in seconds");
		timeOp.setRequired(false);
		ops.addOption(timeOp);
		
		Option statusFile = new Option("f", "file", true, "Status file for writing job status");
		statusFile.setRequired(false);
		ops.addOption(statusFile);
		
		return ops;
	}
	
	public static void main(String[] args) {
		args = MPJTaskCalculator.initMPJ(args);
		
		try {
			Options options = createOptions();
			
			CommandLine cmd = parse(options, args, MPJ_DummyTest.class);
			
			MPJ_DummyTest driver = new MPJ_DummyTest(cmd);
			driver.run();
			
			finalizeMPJ();
			
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}

}
