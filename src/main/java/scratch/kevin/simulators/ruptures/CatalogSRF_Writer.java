package scratch.kevin.simulators.ruptures;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.text.DecimalFormat;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.srf.RSQSimEventSlipTimeFunc;
import org.opensha.sha.simulators.srf.RSQSimSRFGenerator;
import org.opensha.sha.simulators.srf.RSQSimSRFGenerator.SRFInterpolationMode;
import org.opensha.sha.simulators.srf.SRF_PointData;
import org.opensha.sha.simulators.utils.SimulatorUtils;

import com.google.common.base.Preconditions;

import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.RSQSimCatalog.Loader;

public class CatalogSRF_Writer {

	public static void main(String[] args) throws IOException {
		File baseDir = new File("/data/kevin/simulators/catalogs");
		RSQSimCatalog catalog = Catalogs.BRUCE_4983_STITCHED.instance(baseDir);
		
		File outputFile = new File(catalog.getCatalogDir(), "srfs.zip");
//		File outputDir = new File(catalog.getCatalogDir(), "srfs");
//		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		double minMag = 6.5;
		double skipYears = 65000;
		double maxDuration = 0;
		
//		double dt = RSQSimBBP_Config.SRF_DT;
		double dt = 0.1;
		SRFInterpolationMode interp = RSQSimBBP_Config.SRF_INTERP_MODE;
		double srfVersion = RSQSimBBP_Config.SRF_VERSION;
		
		Loader loader = catalog.loader().minMag(minMag).skipYears(skipYears);
		if (maxDuration > 0 && Double.isFinite(maxDuration))
			loader.maxDuration(maxDuration);
		List<RSQSimEvent> events = loader.load();
		System.out.println("Loaded "+events.size()+" events in "+SimulatorUtils.getSimulationDurationYears(events)+" yrs");
		
		ZipOutputStream zout = new ZipOutputStream(new FileOutputStream(outputFile));
		
		BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(zout));
		
		int cnt = 0;
		DecimalFormat pDF = new DecimalFormat("0.0%");
		for (RSQSimEvent e : events) {
			if (cnt % 100 == 0)
				System.out.println("Writing event "+cnt+"/"+events.size()+" ("+pDF.format((double)cnt/(double)events.size())+")");
			cnt++;
			RSQSimEventSlipTimeFunc slipTimeFunc = catalog.getSlipTimeFunc(e);
//			File srfFile = new File(outputDir, "event_"+e.getID()+".srf");
			zout.putNextEntry(new ZipEntry("event_"+e.getID()+".srf"));
			List<SRF_PointData> srfPoints = RSQSimSRFGenerator.buildSRF(slipTimeFunc, e.getAllElements(), dt, interp);
			SRF_PointData.writeSRF(writer, srfPoints, srfVersion);
			writer.flush();
			zout.flush();
			zout.closeEntry();
		}
		
		zout.close();
	}

}
