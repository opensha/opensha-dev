package scratch.kevin.simulators.erf;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.parsers.EQSIMv06FileReader;
import org.opensha.sha.simulators.utils.General_EQSIM_Tools;

public class ScratchJFreeDebug {

	public static void main(String[] args) throws IOException {
		File dataDir = new File("/home/kevin/Simulators");
		File geomFile = new File(dataDir, "ALLCAL2_1-7-11_Geometry.dat");
		System.out.println("Loading geometry from "+geomFile.getAbsolutePath());
		General_EQSIM_Tools tools = new General_EQSIM_Tools(geomFile);
		File eventFile = new File(dataDir, "eqs.ALLCAL2_RSQSim_sigma0.5-5_b=0.015.barall");
//		File eventFile = new File(dataDir, "eqs.ALLCAL2_RSQSim_sigma0.5-5_b=0.015.long.barall");
		System.out.println("Loading events...");
		
		List<? extends SimulatorEvent> events = EQSIMv06FileReader.readEventsFile(eventFile, tools.getElementsList());
		
		tools.testTimePredictability(Double.NaN, false, null, false);
	}

}
