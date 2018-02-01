package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.srf.RSQSimEventSlipTimeFunc;
import org.opensha.sha.simulators.srf.RSQSimSRFGenerator;
import org.opensha.sha.simulators.srf.RSQSimSRFGenerator.SRFInterpolationMode;
import org.opensha.sha.simulators.srf.SRF_PointData;
import org.opensha.sha.simulators.utils.SimulatorUtils;

import com.google.common.base.Preconditions;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class CatalogSRF_Writer {

	public static void main(String[] args) throws IOException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		RSQSimCatalog catalog = Catalogs.BRUCE_2457.instance(baseDir);
		
		File outputDir = new File(catalog.getCatalogDir(), "srfs");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		double minMag = 6.5;
		double skipYears = 5000;
		double maxDuration = 10000;
		
		double dt = RSQSimBBP_Config.SRF_DT;
		SRFInterpolationMode interp = RSQSimBBP_Config.SRF_INTERP_MODE;
		double srfVersion = RSQSimBBP_Config.SRF_VERSION;
		
		List<RSQSimEvent> events = catalog.loader().minMag(minMag).skipYears(skipYears).maxDuration(maxDuration).load();
		System.out.println("Loaded "+events.size()+" events in "+SimulatorUtils.getSimulationDurationYears(events)+" yrs");
		
		int cnt = 0;
		for (RSQSimEvent e : events) {
			if (cnt % 100 == 0)
				System.out.println("Writing event "+cnt+"/"+events.size());
			cnt++;
			RSQSimEventSlipTimeFunc slipTimeFunc = catalog.getSlipTimeFunc(e);
			File srfFile = new File(outputDir, "event_"+e.getID()+".srf");
			List<SRF_PointData> srfPoints = RSQSimSRFGenerator.buildSRF(slipTimeFunc, e.getAllElements(), dt, interp);
			SRF_PointData.writeSRF(srfFile, srfPoints, srfVersion);
		}
	}

}
