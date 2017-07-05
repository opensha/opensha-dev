package scratch.peter.tmp;

import java.awt.Color;
import java.util.Iterator;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.RegionUtils;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.rupForecastImpl.FaultRuptureSource;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.FaultSegmentData;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.UCERF2;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.data.A_FaultsFetcher;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.utils.GriddedSurfaceUtils;

import com.google.common.collect.Iterators;

public class distXtests {

	public static void main(String[] args) {
//		A_FaultsFetcher fetcher = new A_FaultsFetcher();
//		FaultSegmentData fsd = fetcher.getFaultSegmentData("San Andreas (North Coast)", false);
//		System.out.println(fsd);
		
		MeanUCERF2 uc2 = new MeanUCERF2();
		uc2.setParameter(UCERF2.BACK_SEIS_NAME, UCERF2.BACK_SEIS_EXCLUDE);
//		uc2.setParameter(UCERF2.RUP_MODEL_TYPE_NAME, UCERF2.UNSEGMENTED_A_FAULT_MODEL);
		uc2.updateForecast();
		int idx = 0;
//		for (ProbEqkSource source : uc2.getSourceList()) {
//			System.out.println(idx++ + " " + source.getName());
//		}
		
		FaultRuptureSource frs = (FaultRuptureSource) uc2.getSource(33);
		
//		LocationList surfLocs = new LocationList();
//		Iterator<Location> it = frs.getSourceSurface().getLocationsIterator();
//		while (it.hasNext()) {
//			surfLocs.add(it.next());
//		}
//		RegionUtils.locListToKML(surfLocs, "NSA_surface", Color.BLUE);
//		
//		LocationList traceLocs = new LocationList();
//		LocationList trace = frs.getSourceSurface().getEvenlyDiscritizedUpperEdge();
//		RegionUtils.locListToKML(trace, "NSA_trace", Color.ORANGE);
		
		FaultTrace trace = frs.getSourceSurface().getEvenlyDiscritizedUpperEdge();
		GriddedSurfaceUtils.getDistanceX(trace, new Location(38, -123));
		
	}
}
