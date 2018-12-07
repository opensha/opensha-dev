package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.srf.RSQSimEventSlipTimeFunc;
import org.opensha.sha.simulators.srf.RSQSimSRFGenerator;
import org.opensha.sha.simulators.srf.RSQSimSRFGenerator.SRFInterpolationMode;

import com.google.common.base.Preconditions;

import org.opensha.sha.simulators.srf.RSQSimStateTransitionFileReader;
import org.opensha.sha.simulators.srf.SRF_PointData;
import org.opensha.sha.simulators.utils.RupturePlotGenerator;

import scratch.kevin.simulators.RSQSimCatalog;

public class VariableSlipSpeedDebug {

	public static void main(String[] args) throws IOException {
		RSQSimCatalog catalog = new RSQSimCatalog(new File("/home/kevin/Simulators/catalogs/variableSlipRun6"),
				"Variable slip speed test", null, null);
		
		List<RSQSimEvent> events = catalog.loader().load();
		System.out.println("Loaded "+events.size()+" events");
		
		System.out.println("Variable slip speed? "+catalog.isVariableSlipSpeed());
		
		// sort by magnitude, decreasing
		Collections.sort(events, new Comparator<RSQSimEvent>() {

			@Override
			public int compare(RSQSimEvent e1, RSQSimEvent e2) {
				return -Double.compare(e1.getMagnitude(), e2.getMagnitude());
			}
		});
		
		int numToPlot = 5;
		double dt = 0.05;
		SRFInterpolationMode mode = SRFInterpolationMode.ADJ_VEL;
		
		for (int i=0; i<numToPlot && i<events.size(); i++) {
			RSQSimEvent event = events.get(i);
			List<SimulatorElement> patches = event.getAllElements();
			System.out.println("Plotting event "+event.getID()+" (M"+(float)event.getMagnitude()+") with "
					+patches.size()+" patches");
			
			File eventDir = new File(catalog.getCatalogDir(), "event_"+event.getID());
			Preconditions.checkState(eventDir.exists() || eventDir.mkdir());
			
			RSQSimStateTransitionFileReader trans = catalog.getTransitions();
			PrintStream transWriter = new PrintStream(new File(eventDir, "event_"+event.getID()+"_transitions.txt"));
			RSQSimStateTransitionFileReader.printTransitions(event, trans.getTransitions(event), transWriter);
			transWriter.close();
			
			RSQSimEventSlipTimeFunc func = catalog.getSlipTimeFunc(event);
			
			List<SRF_PointData> srf = RSQSimSRFGenerator.buildSRF(func, patches, dt, mode);
			File srfFile = new File(eventDir, "event_"+event.getID()+"_"+(float)dt+"s_"+mode+".srf");
			System.out.println("Writing SRF: "+srfFile.getAbsolutePath());
			SRF_PointData.writeSRF(srfFile, srf, 1d);
			int mod = patches.size() < 30 ? 1 : patches.size() / 20;
			for (int p=0; p<patches.size(); p += mod) {
				SimulatorElement patch = patches.get(p);
				RSQSimSRFGenerator.plotSlip(eventDir, "patch_"+patch.getID(), func, patch, dt, mode);
			}
			System.out.println("Plotting slip animation");
			RupturePlotGenerator.writeSlipAnimation(event, func, new File(eventDir, "event_"+event.getID()+".gif"), 10);
		}
	}

}
