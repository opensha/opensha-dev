package scratch.kevin.prvi25;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.util.ClassUtils;
import org.opensha.commons.util.modules.AverageableModule.AveragingAccumulator;
import org.opensha.commons.util.modules.ModuleArchive;
import org.opensha.commons.util.modules.OpenSHA_Module;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.modules.RegionsOfInterest;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

public class RegionalMFD_Spider {

	public static void main(String[] args) throws IOException {
		if (args.length != 1) {
			System.err.println("USAGE: <results-dir>");
			System.exit(1);
		}
		File dir = new File(args[0]);
		
		AveragingAccumulator<RegionsOfInterest> averager = null;
		
		ModuleArchive.VERBOSE_DEFAULT = false;
		
		for (File subdir : dir.listFiles()) {
			if (!subdir.isDirectory())
				continue;
			File solFile = new File(subdir, "solution.zip");
			if (!solFile.exists())
				continue;
			System.out.println(subdir.getName());
			ModuleArchive<OpenSHA_Module> archive = new ModuleArchive<>(solFile);
			RegionsOfInterest roi = archive.loadUnlistedModule(RegionsOfInterest.class, FaultSystemRupSet.NESTING_PREFIX);
			archive.getInput().close();
			
			List<IncrementalMagFreqDist> mfds = roi.getMFDs();
			System.out.println("\tRegion count: "+roi.getRegions().size());
			System.out.println("\tMFDs: "+(mfds == null ? "NULL" : mfds.size()+""));
			if (mfds != null) {
				for (int i=0; i<mfds.size(); i++) {
					IncrementalMagFreqDist mfd = mfds.get(i);
					System.out.println("\t\tMFD "+i+": "+(mfd == null ? "NULL" :
							mfd.getName()+" ("+ClassUtils.getClassNameWithoutPackage(mfd.getClass())+")"));
				}
			}
			List<TectonicRegionType> trts = roi.getTRTs();
			System.out.println("\tTRTs: "+(trts == null ? "NULL" : trts.size()+""));
			if (trts != null) {
				for (int i=0; i<trts.size(); i++) {
					TectonicRegionType trt = trts.get(i);
					System.out.println("\t\tTRT "+i+": "+(trt == null ? "NULL" : trt));
				}
			}
			
			if (averager == null)
				averager = roi.averagingAccumulator();
			averager.process(roi, 1d);
		}
		
		System.out.println("Building average");
		RegionsOfInterest roi = averager.getAverage();
		List<IncrementalMagFreqDist> mfds = roi.getMFDs();
		System.out.println("\tRegion count: "+roi.getRegions().size());
		System.out.println("\tMFDs: "+(mfds == null ? "NULL" : mfds.size()+""));
		if (mfds != null) {
			for (int i=0; i<mfds.size(); i++) {
				IncrementalMagFreqDist mfd = mfds.get(i);
				System.out.println("\t\tMFD "+i+": "+(mfd == null ? "NULL" :
						mfd.getName()+" ("+ClassUtils.getClassNameWithoutPackage(mfd.getClass())+")"));
			}
		}
		List<TectonicRegionType> trts = roi.getTRTs();
		System.out.println("\tTRTs: "+(trts == null ? "NULL" : trts.size()+""));
		if (trts != null) {
			for (int i=0; i<trts.size(); i++) {
				TectonicRegionType trt = trts.get(i);
				System.out.println("\t\tTRT "+i+": "+(trt == null ? "NULL" : trt));
			}
		}
	}

}
