package scratch.kevin.cybershake.ucerf3.safWallToWallTests;

import java.util.Map;
import java.util.concurrent.ConcurrentMap;

import org.opensha.commons.geo.Location;
import org.opensha.sha.cybershake.db.DBAccess;
import org.opensha.sha.cybershake.db.ERF2DB;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.faultSurface.EvenlyGriddedSurface;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.InterpolatedEvenlyGriddedSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.StirlingGriddedSurface;

import com.google.common.collect.Maps;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;

public class FSS_ERF2DB extends ERF2DB {
	
	private double gridSpacing;
	
	public FSS_ERF2DB(FaultSystemSolution sol, DBAccess db, boolean ucerf2StyleSurfaces, double gridSpacing) {
		super(db);
		if (ucerf2StyleSurfaces) {
			sol = new FaultSystemSolution(new UCERF2StyleSolRupSetWrapper(sol.getRupSet(), gridSpacing), sol.getRateForAllRups());
		} else {
			throw new IllegalStateException("Non UCERF2 style surfaces not yet supported");
		}
		this.gridSpacing = gridSpacing;
		this.eqkRupForecast = getUCERF3_ERF(sol);
	}
	
	public static FaultSystemSolutionERF getUCERF3_ERF(FaultSystemSolution sol) {
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
		erf.getTimeSpan().setDuration(1d);
		
		erf.updateForecast();
		
		return erf;
	}
	
	private static class UCERF2StyleSolRupSetWrapper extends FaultSystemRupSet {
		
		private ConcurrentMap<Integer, RuptureSurface> surfCache = Maps.newConcurrentMap();
		private double gridSpacing;
		private double createGridSpacing;
		
		public UCERF2StyleSolRupSetWrapper(FaultSystemRupSet rupSet, double gridSpacing) {
			init(rupSet);
			this.gridSpacing = gridSpacing;
			if (gridSpacing < 0.5)
				createGridSpacing = 1d;
			else
				createGridSpacing = gridSpacing;
		}

		@Override
		public synchronized RuptureSurface getSurfaceForRupupture(int rupIndex, double gridSpacing,
				boolean quadRupSurface) {
			// use global grid spacing
			gridSpacing = this.createGridSpacing;
			RuptureSurface cached = surfCache.get(rupIndex);
			if (cached != null)
				return cached;
			
			RuptureSurface combSurface = super.getSurfaceForRupupture(rupIndex, gridSpacing, quadRupSurface);
			double aveTopDepth = combSurface.getAveRupTopDepth();
			double aveWidth = combSurface.getAveWidth();
			
			FaultTrace upperEdge = combSurface.getUpperEdge();
			// now make it constant upper depth
			for (int i=0; i<upperEdge.size(); i++) {
				Location loc = upperEdge.get(i);
				if (loc.getDepth() != 0)
					upperEdge.set(i, new Location(loc.getLatitude(), loc.getLongitude(), 0d));
			}
			EvenlyGriddedSurface surf = new StirlingGriddedSurface(upperEdge, combSurface.getAveDip(),
					aveTopDepth, aveTopDepth+aveWidth, gridSpacing, combSurface.getAveDipDirection());
			if (this.createGridSpacing != this.gridSpacing)
				surf = new InterpolatedEvenlyGriddedSurface(surf, this.gridSpacing);
			
			surfCache.putIfAbsent(rupIndex, surf);
			
			return surf;
		}
		
	}

}
