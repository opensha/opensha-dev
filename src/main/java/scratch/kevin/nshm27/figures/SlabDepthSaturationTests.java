package scratch.kevin.nshm27.figures;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.util.ColorUtils;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.rupForecastImpl.PointSourceNshm.DistanceCorrection2013;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.gridded.NSHM23_SingleRegionGridSourceProvider.NSHM23GridRupDimensions;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.gridded.NSHM23_SingleRegionGridSourceProvider.NSHM23_WUS_FiniteRuptureConverter;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.util.TectonicRegionType;

public class SlabDepthSaturationTests {

	private static final boolean CRUSTAL = true;
	private static final double[] MAGNITUDES = { 7d, 7.5d, 8d };
	private static final double RAKE = 0d;
	private static final double[] EPICENTRAL_DISTANCES = { 0d, 10d, 20d, 50d };
	private static final EvenlyDiscretizedFunc DEPTHS = new EvenlyDiscretizedFunc(0d, 301, 1d);
	
	private static final Range Y_RANGE = new Range(1e-2, 1e0);

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/tmp");
		Location sourceHorzLoc = new Location(0d, 0d);
		DistanceCorrection2013 distanceCorrection = new DistanceCorrection2013();
		DecimalFormat magDF = new DecimalFormat("0.#");

		TectonicRegionType trt = CRUSTAL ? TectonicRegionType.ACTIVE_SHALLOW
				: TectonicRegionType.SUBDUCTION_SLAB;
		ScalarIMR gmm = (CRUSTAL ? AttenRelRef.USGS_NSHM23_ACTIVE
				: AttenRelRef.USGS_NSHM27_TEMP_SLAB).get();
		gmm.setParamDefaults();

		for (double magnitude : MAGNITUDES) {
			NSHM23GridRupDimensions dimensions = CRUSTAL
					? NSHM23_WUS_FiniteRuptureConverter.getDimensions(trt, magnitude, 90d)
					: NSHM23_WUS_FiniteRuptureConverter.getSlabDimensions(magnitude, 0d);
			String magLabel = magDF.format(magnitude);
			String magPrefix = "m"+magLabel.replace('.', 'p');

			for (double period : new double[] { 0d, 1d }) {
				String imtName;
				String prefix;
				if (period == 0d) {
					gmm.setIntensityMeasure(PGA_Param.NAME);
					imtName = "PGA";
					prefix = "pga";
				} else {
					gmm.setIntensityMeasure(SA_Param.NAME);
					SA_Param.setPeriodInSA_Param(gmm.getIntensityMeasure(), period);
					imtName = "1s SA";
					prefix = "sa_1s";
				}
				prefix = (CRUSTAL ? "nshm27_crustal_" : "nshm27_slab_")+magPrefix+"_"+prefix;

				List<DiscretizedFunc> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();

				for (int d=0; d<EPICENTRAL_DISTANCES.length; d++) {
					double epicentralDistance = EPICENTRAL_DISTANCES[d];
					Location siteLoc = LocationUtils.location(sourceHorzLoc, 0d, epicentralDistance);
					Site site = new Site(siteLoc);
					site.addParameterList(gmm.getSiteParams());
					gmm.setSite(site);

					EvenlyDiscretizedFunc truePointFunc = DEPTHS.deepClone();
					truePointFunc.setName((int)epicentralDistance+" km, true point");
					EvenlyDiscretizedFunc correctedFunc = DEPTHS.deepClone();
					correctedFunc.setName((int)epicentralDistance+" km, distance corrected");

					for (int i=0; i<DEPTHS.size(); i++) {
						double upperDepth = DEPTHS.getX(i);
						double lowerDepth = upperDepth + dimensions.width();
						Location sourceLoc = new Location(sourceHorzLoc.lat, sourceHorzLoc.lon, upperDepth);

						// Keep all rupture properties identical except for the finite length used by the
						// NSHM13 distance correction. A zero length is the true-point convention.
						PointSurface truePointSurface = new PointSurface(sourceLoc, dimensions.dip(),
								upperDepth, lowerDepth, 0d);
						EqkRupture truePointRupture = new EqkRupture(magnitude, RAKE, truePointSurface, null);
						truePointFunc.set(i, Math.exp(gmm.getMean(truePointRupture)));

						PointSurface correctedPointSurface = new PointSurface(sourceLoc, dimensions.dip(),
								upperDepth, lowerDepth, dimensions.length());
						EqkRupture correctedRupture = new EqkRupture(magnitude, RAKE,
								correctedPointSurface.getForDistanceCorrection(distanceCorrection,
										trt, magnitude), null);
						correctedFunc.set(i, Math.exp(gmm.getMean(correctedRupture)));
					}

					funcs.add(truePointFunc);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SHORT_DASHED, 3f, ColorUtils.TAB_10[d]));
					funcs.add(correctedFunc);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, ColorUtils.TAB_10[d]));
				}

				PlotSpec plot = new PlotSpec(funcs, chars, "M"+magLabel+" "
						+(CRUSTAL ? "Crustal" : "Intraslab")+" Point Sources",
						"Depth to Top of Rupture (km)", "Median "+imtName+" (g)");
				plot.setLegendInset(true);

				HeadlessGraphPanel gp = PlotUtils.initScreenHeadless();
				gp.drawGraphPanel(plot, false, true, new Range(DEPTHS.getMinX(), DEPTHS.getMaxX()), Y_RANGE);
				PlotUtils.writePlots(outputDir, prefix, gp, 900, 750, true, true, false);
			}
		}
	}

}
