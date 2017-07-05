package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.dom4j.DocumentException;
import org.opensha.commons.util.ClassUtils;

import com.google.common.base.Preconditions;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.LastEventData;

public class BranchAverageLastEventFixer {

	/**
	 * branch averaged solutions are no longer Inversion FSS's, so they don't have last event data loaded 
	 * by default. therefore we need to add the data to the files themselves. This will
	 * @param args
	 * @throws IOException 
	 * @throws DocumentException 
	 */
	public static void main(String[] args) throws IOException, DocumentException {
		if (args.length != 1) {
			System.err.println("USAGE: "+ClassUtils.getClassNameWithoutPackage(BranchAverageLastEventFixer.class)+" <dir>");
			System.exit(2);
		}
		Map<Integer, List<LastEventData>> data = LastEventData.load();
		File dir = new File(args[0]);
		Preconditions.checkState(dir.exists());
		if (dir.isFile()) {
			Preconditions.checkState(dir.getName().endsWith(".zip"), "If file supplied must be a zip file");
			handleSolFile(dir, data);
		} else {
			handleDir(dir, data);
		}
	}
	
	private static void handleDir(File dir, Map<Integer, List<LastEventData>> data) throws IOException, DocumentException {
		for (File file : dir.listFiles()) {
			String fName = file.getName().toLowerCase();
			if (fName.endsWith(".zip") && fName.contains("branch_avg"))
				handleSolFile(file, data);
			else if (file.isDirectory() && !fName.startsWith("."))
				handleDir(file, data);
		}
	}
	
	private static void handleSolFile(File file, Map<Integer, List<LastEventData>> data) throws IOException, DocumentException {
		System.out.println("Loading "+file.getAbsolutePath());
		FaultSystemSolution sol = FaultSystemIO.loadSol(file);
		LastEventData.populateSubSects(sol.getRupSet().getFaultSectionDataList(), data);
		System.out.println("Writing populated sol");
		FaultSystemIO.writeSol(sol, file);
	}

}
