package scratch.kevin.nshm23.wrapper;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Deque;
import java.util.List;
import java.util.concurrent.ConcurrentLinkedDeque;
import java.util.concurrent.TimeUnit;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.util.FileUtils;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_LogicTreeHazardCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;

import gov.usgs.earthquake.nshmp.model.NshmErf;

public class WrapperHazardCalc {

	public static void main(String[] args) throws IOException {
		Path erfPath = Path.of("/home/kevin/git/nshm-conus");
		NshmErf erf = new NshmErf(erfPath, false, true);
		System.out.println("NSHM ERF size: " + erf.getNumSources());
		erf.getTimeSpan().setDuration(1.0);
		erf.updateForecast();

		AttenRelRef gmpeRef = AttenRelRef.ASK_2014;

		Region modelReg = NSHM23_RegionLoader.loadFullConterminousWUS();
		// Region modelReg = new Region(new Location(34, -118), 10d);
		GriddedRegion gridReg = new GriddedRegion(modelReg, 0.2d,
				GriddedRegion.ANCHOR_0_0);
		File outputDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_09_23-nshm18-hazard-ask2014-noSub/results");
//				+ "2022_09_23-nshm18-hazard-ask2014-noSub-faultOnly/results");

//		Region modelReg = new CaliforniaRegions.RELM_TESTING();
//		GriddedRegion gridReg = new GriddedRegion(modelReg, 0.2d, GriddedRegion.ANCHOR_0_0);
//		File outputDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
////				+ "2022_09_21-nshm18-ca-hazard-ask2014/results");
////				+ "2022_09_23-nshm18-ca-hazard-ask2014-noSub-faultOnly/results");
//				+ "2022_09_23-nshm18-ca-hazard-ask2014-noSub/results");
		
		int threads = 8;
		
		Preconditions.checkState(outputDir.getParentFile().exists() || outputDir.getParentFile().mkdir());
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());

		System.out.println("Will calculate for " + gridReg.getNodeCount() + " locations");

		IMT_Info imtInfo = new IMT_Info();

		double[] periods = { 0d, 1d };
		DiscretizedFunc[] linearXVals = {
				imtInfo.getDefaultHazardCurve(PGA_Param.NAME),
				imtInfo.getDefaultHazardCurve(SA_Param.NAME) };

		ReturnPeriods[] rps = ReturnPeriods.values();

		DiscretizedFunc[] logXVals = new DiscretizedFunc[linearXVals.length];
		for (int i = 0; i < logXVals.length; i++) {
			logXVals[i] = new ArbitrarilyDiscretizedFunc();
			for (Point2D pt : linearXVals[i])
				logXVals[i].set(Math.log(pt.getX()), 0d);
			// System.out.println("Log curve for p="+periods[i]+"\n"+logXVals[i]);
		}

		DiscretizedFunc[][] curves = new DiscretizedFunc[periods.length][gridReg.getNodeCount()];

		List<Integer> calcIndexes = new ArrayList<>();
		for (int s = 0; s < gridReg.getNodeCount(); s++)
			calcIndexes.add(s);
		// shuffle them to reduce thread contention
		Collections.shuffle(calcIndexes);
		ConcurrentLinkedDeque<Integer> deque = new ConcurrentLinkedDeque<>(calcIndexes);

		Stopwatch watch = Stopwatch.createStarted();
		List<CalcThread> calcThreads = new ArrayList<>();
		for (int t = 0; t < threads; t++) {
			CalcThread thread =
					new CalcThread(erf, gmpeRef, gridReg, periods, logXVals, linearXVals, curves, deque);
			thread.start();

			calcThreads.add(thread);
		}

		for (CalcThread thread : calcThreads) {
			try {
				thread.join();
			} catch (InterruptedException e) {
				e.printStackTrace();
				System.exit(1);
			}
		}
		watch.stop();
		double mins = watch.elapsed(TimeUnit.SECONDS) / 60d;
		System.out.println(
				"Took " + (float) mins + " minutes to calculate " + gridReg.getNodeCount() + " sites");
		System.out.println("Rate: " + (mins * 60d / gridReg.getNodeCount()) + " curves per second");

		List<String> zipFileNames = new ArrayList<>();

		for (int p = 0; p < periods.length; p++) {
			File curvesFile =
					new File(outputDir, SolHazardMapCalc.getCSV_FileName("curves", periods[p]) + ".gz");
			SolHazardMapCalc.writeCurvesCSV(curvesFile, curves[p], gridReg.getNodeList());
			zipFileNames.add(curvesFile.getName());

			for (ReturnPeriods rp : rps) {
				GriddedGeoDataSet xyz = new GriddedGeoDataSet(gridReg, false);

				double curveLevel = rp.oneYearProb;

				for (int i = 0; i < curves[p].length; i++) {
					DiscretizedFunc curve = curves[p][i];
					Preconditions.checkNotNull(curve, "Curve not calculated at index %s", i);
					double val;
					// curveLevel is a probability, return the IML at that probability
					if (curveLevel > curve.getMaxY())
						val = 0d;
					else if (curveLevel < curve.getMinY())
						// saturated
						val = curve.getMaxX();
					else
						val = curve.getFirstInterpolatedX_inLogXLogYDomain(curveLevel);

					xyz.set(i, val);
				}

				String mapPrefix = MPJ_LogicTreeHazardCalc.mapPrefix(periods[p], rp);
				File mapFile = new File(outputDir, mapPrefix + ".txt");
				zipFileNames.add(mapFile.getName());

				GriddedGeoDataSet.writeXYZFile(xyz, mapFile);
			}
		}

		File gridFile = new File(outputDir, "gridded_region.geojson");
		Feature.write(gridReg.toFeature(), gridFile);
		zipFileNames.add(gridFile.getName());

		File zipFile = new File(outputDir.getAbsolutePath() + ".zip");

		FileUtils.createZipFile(zipFile.getAbsolutePath(), outputDir.getAbsolutePath(), zipFileNames);
	}

	public static class CalcThread extends Thread {

		private AbstractERF erf;
		private AttenRelRef gmpeRef;
		private GriddedRegion gridReg;
		private double[] periods;
		private DiscretizedFunc[] logXVals;
		private DiscretizedFunc[] linearXVals;
		private DiscretizedFunc[][] curves;

		private Deque<Integer> siteIndexDeque;

		public CalcThread(NshmErf erf, AttenRelRef gmpeRef, GriddedRegion gridReg, double[] periods,
				DiscretizedFunc[] logXVals, DiscretizedFunc[] linearXVals, DiscretizedFunc[][] curves,
				Deque<Integer> siteIndexDeque) {
			// this.erf = new DistCachedERFWrapper(erf);
			this.erf = erf;
			this.gmpeRef = gmpeRef;
			this.gridReg = gridReg;
			this.periods = periods;
			this.logXVals = logXVals;
			this.linearXVals = linearXVals;
			this.curves = curves;
			this.siteIndexDeque = siteIndexDeque;
		}

		@Override
		public void run() {
			ScalarIMR gmpe = gmpeRef.instance(null);
			gmpe.setParamDefaults();

			HazardCurveCalculator calc = new HazardCurveCalculator();

			while (true) {
				Integer index = siteIndexDeque.pollFirst();
				if (index == null)
					break;

				Location loc = gridReg.getLocation(index);
				Site site = new Site(loc);
				// site.addParameterList(gmpe.getSiteParams());
				for (Parameter<?> param : gmpe.getSiteParams()) {
					site.addParameter((Parameter<?>) param.clone());
				}
				site.getParameter(Double.class, Vs30_Param.NAME).setValue(760.0);

				Stopwatch watch = Stopwatch.createStarted();

				for (int p = 0; p < periods.length; p++) {
					if (periods[p] == 0d) {
						gmpe.setIntensityMeasure(PGA_Param.NAME);
					} else {
						Preconditions.checkState(periods[p] > 0d);
						gmpe.setIntensityMeasure(SA_Param.NAME);
						SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), periods[p]);
					}

					DiscretizedFunc logCurve = logXVals[p].deepClone();
					// System.out.println("Input curve:\n"+logCurve);

					calc.getHazardCurve(logCurve, site, gmpe, erf);

					// System.out.println("Curve for "+loc+":\n"+logCurve);

					DiscretizedFunc linearCurve = linearXVals[p].deepClone();
					Preconditions.checkState(linearCurve.size() == logCurve.size());
					for (int i = 0; i < logCurve.size(); i++) {
						double y = logCurve.getY(i);
						Preconditions.checkState(Double.isFinite(y), "bad y=%s for x=%s, logX=%s",
								y, linearCurve.getX(i), logCurve.getX(i));
						linearCurve.set(i, y);
					}
					curves[p][index] = linearCurve;
				}
				watch.stop();
				double secs = watch.elapsed(TimeUnit.MILLISECONDS) / 1000d;
				int nodeCount = gridReg.getNodeCount();
				int done = nodeCount - siteIndexDeque.size();
				System.out.println("DONE site " + index + "/" + nodeCount + " in " + (float) secs + " s" +
						"\t(" + (float) (100d * done / nodeCount) + " % completed or running)");
			}
		}

	}

}
