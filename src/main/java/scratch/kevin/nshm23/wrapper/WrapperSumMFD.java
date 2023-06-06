package scratch.kevin.nshm23.wrapper;

import java.nio.file.Path;
import java.util.EnumSet;
import java.util.Set;

import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import gov.usgs.earthquake.nshmp.model.NshmErf;

public class WrapperSumMFD {

	public static void main(String[] args) {
		Path erfPath = Path.of("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/nshm-conus-5.2.0-wasatch");
		IncludeBackgroundOption griddedOp = IncludeBackgroundOption.EXCLUDE;
		boolean subduction = false;
		
		Set<TectonicRegionType> trts = EnumSet.of(TectonicRegionType.ACTIVE_SHALLOW);
		if (subduction) {
			trts.add(TectonicRegionType.SUBDUCTION_INTERFACE);
			trts.add(TectonicRegionType.SUBDUCTION_SLAB);
		}

		NshmErf erf = new NshmErf(erfPath, trts, griddedOp);
		System.out.println("NSHM ERF size: " + erf.getNumSources());
		erf.getTimeSpan().setDuration(1.0);
		erf.updateForecast();
		
		EvenlyDiscretizedFunc refXVals = FaultSysTools.initEmptyMFD(8.5);
		IncrementalMagFreqDist mfd = new IncrementalMagFreqDist(refXVals.getMinX(), refXVals.getMaxX(), refXVals.size());
		
		for (ProbEqkSource src : erf) {
			for (ProbEqkRupture rup : src) {
				double mag = rup.getMag();
				double rate = rup.getMeanAnnualRate(1d);
				
				mfd.add(mfd.getClosestXIndex(mag), rate);
			}
		}
		
		System.out.println("Incremental MFD:\n\n"+mfd);
		
		System.out.println("\n\nCumulative MFD:\n\n"+mfd.getCumRateDistWithOffset());
	}

}
