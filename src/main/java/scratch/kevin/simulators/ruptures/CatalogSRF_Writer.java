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

import org.opensha.commons.util.io.archive.ArchiveOutput;
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
		RSQSimCatalog catalog;
		double minMag = 6.5;
		double skipYears = 5000;
		if (args.length > 0) {
			catalog = Catalogs.valueOf(args[0]).instance();
			if (args.length > 1)
				minMag = Double.parseDouble(args[1]);
			if (args.length > 2)
				skipYears = Double.parseDouble(args[2]);
		} else {
			catalog = Catalogs.BRUCE_5566.instance();
		}
		
		File outputFile = new File(catalog.getCatalogDir(), "srfs.zip");
//		File outputDir = new File(catalog.getCatalogDir(), "srfs");
//		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		double maxDuration = 0;
		
		double dt = RSQSimBBP_Config.SRF_DT;
//		double dt = 0.1;
		SRFInterpolationMode interp = RSQSimBBP_Config.SRF_INTERP_MODE;
		double srfVersion = RSQSimBBP_Config.SRF_VERSION;
		
		Loader loader = catalog.loader().minMag(minMag).skipYears(skipYears);
		if (maxDuration > 0 && Double.isFinite(maxDuration))
			loader.maxDuration(maxDuration);
		List<RSQSimEvent> events = loader.load();
		System.out.println("Loaded "+events.size()+" events in "+SimulatorUtils.getSimulationDurationYears(events)+" yrs");
		
		ArchiveOutput output = new ArchiveOutput.AsynchronousZipFileOutput(outputFile);
		
		int cnt = 0;
		DecimalFormat pDF = new DecimalFormat("0.0%");
		for (RSQSimEvent e : events) {
			if (cnt % 100 == 0)
				System.out.println("Writing event "+cnt+"/"+events.size()+" ("+pDF.format((double)cnt/(double)events.size())+")");
			cnt++;
			RSQSimEventSlipTimeFunc slipTimeFunc = catalog.getSlipTimeFunc(e);
//			File srfFile = new File(outputDir, "event_"+e.getID()+".srf");
			output.putNextEntry("event_"+e.getID()+".srf");;
			List<SRF_PointData> srfPoints = RSQSimSRFGenerator.buildSRF(slipTimeFunc, e.getAllElements(), dt, interp);
			BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(output.getOutputStream()));
			SRF_PointData.writeSRF(writer, srfPoints, srfVersion);
			writer.flush();
			output.closeEntry();
		}
		
		output.close();
	}

}
