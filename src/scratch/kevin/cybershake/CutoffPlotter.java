package scratch.kevin.cybershake;

import java.awt.Color;
import java.io.File;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.UncertainArbDiscDataset;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.util.DataUtils;
import org.opensha.sha.cybershake.calc.HazardCurveComputation;
import org.opensha.sha.cybershake.db.CachedPeakAmplitudesFromDB;
import org.opensha.sha.cybershake.db.CybershakeIM;
import org.opensha.sha.cybershake.db.CybershakeSiteInfo2DB;
import org.opensha.sha.cybershake.db.Cybershake_OpenSHA_DBApplication;
import org.opensha.sha.cybershake.db.DBAccess;
import org.opensha.sha.cybershake.db.MeanUCERF2_ToDB;
import org.opensha.sha.cybershake.db.Runs2DB;
import org.opensha.sha.cybershake.db.CybershakeIM.CyberShakeComponent;
import org.opensha.sha.cybershake.db.CybershakeIM.IMType;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.primitives.Doubles;

public class CutoffPlotter {
	
	public static void main(String[] args) throws SQLException {
		ERF erf = MeanUCERF2_ToDB.createUCERF2ERF();
		
		DBAccess db = Cybershake_OpenSHA_DBApplication.getDB();
		CachedPeakAmplitudesFromDB amps2db = new CachedPeakAmplitudesFromDB(
				db, new File("/home/kevin/CyberShake/MCER/.amps_cache"), erf);
		Runs2DB runs2db = new Runs2DB(db);
		CybershakeSiteInfo2DB sites2db = new CybershakeSiteInfo2DB(db);
		
		CybershakeIM im = new CybershakeIM(162, IMType.SA, 3d, null, CyberShakeComponent.RotD50);
//		CybershakeIM im = new CybershakeIM(146, IMType.SA, 3d, null, CyberShakeComponent.RotD100);
		
//		List<Integer> runIDs = Lists.newArrayList(2657);
		List<Integer> runIDs = Lists.newArrayList(2657,3037,2722,3022,3030,3027,2636,2638,2660,2703,3504,2988,2965,3007);
		
		double minDist = 55d;
		int numDist = 15;
		double deltaDist = 10d;
		
		double cutoffMin = minDist - 0.5*deltaDist;
		
		EvenlyDiscretizedFunc xVals = new EvenlyDiscretizedFunc(minDist, numDist, deltaDist);
		Preconditions.checkState((float)xVals.getMaxX() == 195f);
		double cutoffMax = xVals.getMaxX() + 0.5*deltaDist;
		
		double minMag = 6.75d;
		int numMag = 4;
		double deltaMag = 0.5;
		
		double cutoffMinMag = minMag - 0.5*deltaDist;
		
		EvenlyDiscretizedFunc magVals = new EvenlyDiscretizedFunc(minMag, numMag, deltaMag);
		Preconditions.checkState((float)magVals.getMaxX() == 8.25f);
		double cutoffMaxMag = magVals.getMaxX() + 0.5*deltaMag;
		
		List<List<List<Double>>> saVals = Lists.newArrayList();
		for (int m=0; m<numMag; m++) {
			List<List<Double>> magSAVals = Lists.newArrayList();
			saVals.add(magSAVals);
			for (int i=0; i<numDist; i++)
				magSAVals.add(new ArrayList<Double>());
		}
		
		for (int runID : runIDs) {
			Location siteLoc = sites2db.getSiteFromDB(runs2db.getSiteID(runID)).createLocation();
			double[][][] vals = amps2db.getAllIM_Values(runID, im);
			
			for (int sourceID=0; sourceID<erf.getNumSources(); sourceID++) {
				if (vals[sourceID] == null)
					continue;
				ProbEqkSource source = erf.getSource(sourceID);
				for (int rupID=0; rupID<source.getNumRuptures(); rupID++) {
					if (vals[sourceID][rupID] == null)
						continue;
					ProbEqkRupture rup = source.getRupture(rupID);
					double mag = rup.getMag();
					if (mag < cutoffMinMag || mag > cutoffMaxMag)
						continue;
					int magIndex = magVals.getClosestXIndex(mag);
					double dist = rup.getRuptureSurface().getDistanceRup(siteLoc);
//					if (dist < cutoffMin)
					if (dist < cutoffMin || dist > cutoffMax)
						continue;
//					Preconditions.checkState((float)dist <= (float)cutoffMax, "dist too big! "+dist+" > "+cutoffMax);
					int index = xVals.getClosestXIndex(dist);
					saVals.get(magIndex).get(index).addAll(Doubles.asList(vals[sourceID][rupID]));
				}
			}
		}
		
		for (int m=0; m<numMag; m++) {
			// now build functions
			ArbitrarilyDiscretizedFunc minFunc = new ArbitrarilyDiscretizedFunc();
			ArbitrarilyDiscretizedFunc lower95Func = new ArbitrarilyDiscretizedFunc();
			ArbitrarilyDiscretizedFunc lower68Func = new ArbitrarilyDiscretizedFunc();
			ArbitrarilyDiscretizedFunc medianFunc = new ArbitrarilyDiscretizedFunc();
			ArbitrarilyDiscretizedFunc upper68Func = new ArbitrarilyDiscretizedFunc();
			ArbitrarilyDiscretizedFunc upper95Func = new ArbitrarilyDiscretizedFunc();
			ArbitrarilyDiscretizedFunc maxFunc = new ArbitrarilyDiscretizedFunc();
			
			for (int index=0; index<numDist; index++) {
				List<Double> saValsList = saVals.get(m).get(index);
				if (saValsList.isEmpty())
					continue;
				double[] saValsArray = Doubles.toArray(saValsList);
				Arrays.sort(saValsArray);
				// convert to g
				for (int i=0; i<saValsArray.length; i++)
					saValsArray[i] /= HazardCurveComputation.CONVERSION_TO_G;
				double x = xVals.getX(index);
				minFunc.set(x, saValsArray[0]);
				lower95Func.set(x, StatUtils.percentile(saValsArray, 2.5d));
				lower68Func.set(x, StatUtils.percentile(saValsArray, 16d));
				medianFunc.set(x, DataUtils.median_sorted(saValsArray));
				upper68Func.set(x, StatUtils.percentile(saValsArray, 84d));
				upper95Func.set(x, StatUtils.percentile(saValsArray, 97.5d));
				maxFunc.set(x, saValsArray[saValsArray.length-1]);
			}
			
			UncertainArbDiscDataset minMaxBounds = new UncertainArbDiscDataset(medianFunc, minFunc, maxFunc);
			UncertainArbDiscDataset percent95Bounds = new UncertainArbDiscDataset(medianFunc, lower95Func, upper95Func);
			UncertainArbDiscDataset percent68Bounds = new UncertainArbDiscDataset(medianFunc, lower68Func, upper68Func);
			
			List<DiscretizedFunc> funcs = Lists.newArrayList();
			List<PlotCurveCharacterstics> chars = Lists.newArrayList();
			
			funcs.add(minMaxBounds);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, Color.GRAY));
			
			funcs.add(percent95Bounds);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, Color.GREEN));
			
			funcs.add(percent68Bounds);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, Color.BLUE));
			
			funcs.add(medianFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
			
			double yMax = 0;
			for (DiscretizedFunc func : funcs) {
				if (func instanceof UncertainArbDiscDataset)
					yMax = Math.max(yMax, ((UncertainArbDiscDataset)func).getUpperMaxY());
				else
					yMax = Math.max(yMax, func.getMaxY());
			}
			
			double magCenter = magVals.getX(m);
			float myMinMag = (float)(magCenter - 0.5*deltaMag);
			float myMaxMag = (float)(magCenter + 0.5*deltaMag);
			
			PlotSpec spec = new PlotSpec(funcs, chars, "M"+myMinMag+" => M"+myMaxMag, "Distance (km)",
					"SA (g), "+im.getComponent().getShortName());
			
			GraphWindow gw = new GraphWindow(spec);
			gw.setAxisRange(cutoffMin, cutoffMax, 0, yMax);
		}
	}
}
