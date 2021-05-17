package scratch.kevin.simulators;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.iden.EventTimeIdentifier;
import org.opensha.sha.simulators.iden.LogicalAndRupIden;
import org.opensha.sha.simulators.iden.MagRangeRuptureIdentifier;
import org.opensha.sha.simulators.parsers.RSQSimFileReader;
import org.opensha.sha.simulators.utils.RSQSimUtils;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import scratch.UCERF3.SlipEnabledSolution;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;

public class RSQSimUCERF3ComparisonPlotGen {

	public static void main(String[] args) throws IOException, GMT_MapException, RuntimeException {
		File geomFile, eventDir, plotDir;
		
		if (args.length > 0) {
			Preconditions.checkState(args.length <= 2, "USAGE: <sim-dir> [<output-dir>]");
			File runDir = new File(args[0]);
			Preconditions.checkState(runDir.exists());
			if (args.length == 2)
				plotDir = new File(args[1], "ucerf3_fss_comparison_plots");
			else
				plotDir = new File(runDir, "ucerf3_fss_comparison_plots");
			geomFile = new File(runDir, "zfault_Deepen.in");
			if (!geomFile.exists()) {
				// not one of Bruce's
				geomFile = null;
				for (File file : runDir.listFiles()) {
					if (file.getName().endsWith(".flt")) {
						geomFile = file;
						break;
					}
				}
				Preconditions.checkNotNull(geomFile);
			}
			System.out.println("Geometry file: "+geomFile.getAbsolutePath());
			System.out.println("Output dir: "+plotDir.getAbsolutePath());
			eventDir = runDir;
		} else {
//			File dir = new File("/home/kevin/Simulators/UCERF3_35kyrs");
//			geomFile = new File(dir, "UCERF3.1km.tri.flt");
//			File dir = new File("/home/kevin/Simulators/UCERF3_125kyrs");
//			geomFile = new File(dir, "UCERF3.D3.1.1km.tri.2.flt");
//			File dir = new File("/home/kevin/Simulators/bruce/rundir1435");
//			geomFile = new File(dir, "zfault_Deepen.in");
			File dir = new File("/home/kevin/Simulators/UCERF3_JG_supraSeisGeo2");
			geomFile = new File(dir, "UCERF3.D3.1.1km.tri.2.flt");
//			for (Location loc : elements.get(0).getVertices())
//				System.out.println(loc);
			eventDir = dir;
			plotDir = new File(eventDir, "ucerf3_fss_comparison_plots");
		}
		
		List<SimulatorElement> elements = RSQSimFileReader.readGeometryFile(geomFile, 11, 'S');
		System.out.println("Loaded "+elements.size()+" elements");
		
		double minMag = 6d;
		List<RSQSimEvent> events = RSQSimFileReader.readEventsFile(eventDir, elements,
				Lists.newArrayList(new LogicalAndRupIden(new EventTimeIdentifier(5000d, Double.POSITIVE_INFINITY, true),
						new MagRangeRuptureIdentifier(minMag, 10d))));
		double duration = events.get(events.size()-1).getTimeInYears() - events.get(0).getTimeInYears();
		System.out.println("First event time: "+events.get(0).getTimeInYears()+", duration: "+duration);
		
		FaultModels fm = FaultModels.FM3_1;
		DeformationModels dm = DeformationModels.GEOLOGIC;
		SlipEnabledSolution sol = RSQSimUtils.buildFaultSystemSolution(RSQSimUtils.getUCERF3SubSectsForComparison(
				fm, dm), elements, events, minMag);
		
		Preconditions.checkState(plotDir.exists() ||  plotDir.mkdir());
//		MFDCalc.writeMFDPlots(elements, events, plotDir, new CaliforniaRegions.RELM_SOCAL(),
//				new CaliforniaRegions.RELM_NOCAL(), new CaliforniaRegions.LA_BOX(), new CaliforniaRegions.NORTHRIDGE_BOX(),
//				new CaliforniaRegions.SF_BOX(), new CaliforniaRegions.RELM_TESTING());
		RSQSimUtils.writeUCERF3ComparisonPlots(sol, fm, dm, plotDir, "rsqsim_comparison");
	}

}
