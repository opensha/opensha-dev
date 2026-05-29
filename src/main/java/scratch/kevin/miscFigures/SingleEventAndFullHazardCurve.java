package scratch.kevin.miscFigures;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.siteData.OrderedSiteDataProviderList;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.util.modules.ModuleArchive;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.sourceFilters.SourceFilterManager;
import org.opensha.sha.calc.sourceFilters.SourceFilters;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.erf.NSHM23_WUS_BranchAveragedERF;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.util.SiteTranslator;

import com.google.common.base.Joiner;

import net.mahdilamb.colormap.Colors;

public class SingleEventAndFullHazardCurve {

	public static void main(String[] args) throws IOException {
		ModuleArchive.VERBOSE_DEFAULT = false;
		NSHM23_WUS_BranchAveragedERF erf = new NSHM23_WUS_BranchAveragedERF();
		
		File outputDir = new File("/tmp");
		
		int duration = 30;
		erf.getTimeSpan().setDuration((double)duration);
		erf.updateForecast();
		
		FaultSystemSolution sol = erf.getSolution();
		FaultSystemRupSet rupSet = sol.getRupSet();
		double targetMag = 7.5;
		int targetParentID = FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "Mojave", "south");
		String targetName = "If a San Andreas M7.5 occurred";
		String hazardCurveName = "From all possible earthquakes in "+duration+" years";
		
		ScalarIMR gmm = AttenRelRef.NGAWest_2014_AVG_NOIDRISS.get();
		String siteName = "USC";
		Site site = new Site(new Location(34.019265173984806, -118.28634952355786));
		SiteTranslator trans = new SiteTranslator();
		OrderedSiteDataProviderList provs = OrderedSiteDataProviderList.createSiteDataProviderDefaults();
		trans.setAllSiteParams(gmm, provs.getBestAvailableData(site.getLocation()));
		site.addParameterList(gmm.getSiteParams());
		System.out.println("Site params:");
		for (Parameter<?> param : site)
			System.out.println("\t"+param.getName()+":\t"+param.getValue());
		
		String xAxisName = "Peak Ground Acceleration (g)";
		gmm.setIntensityMeasure(PGA_Param.NAME);
//		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(PGA_Param.NAME);
		EvenlyDiscretizedFunc xVals = new EvenlyDiscretizedFunc(0d, 1d, 200);
		double[] logXarray = new double[xVals.size()];
		for (int i=0; i<logXarray.length; i++)
			logXarray[i] = Math.log(xVals.getX(i));
		DiscretizedFunc logXVals = new LightFixedXFunc(logXarray, new double[logXarray.length]);
		HazardCurveCalculator calc = new HazardCurveCalculator(new SourceFilterManager(SourceFilters.TRT_DIST_CUTOFFS));
		Range<Double> magRange = Range.of(targetMag-0.05, targetMag+0.05);
		System.out.println("Mojave Parent: "+targetParentID);
		
		String[] avoidParents = {
				"San Jacinto",
				"Reservoir Canyon",
				"Crafton Hills",
				"San Gabriel",
				"Cucamonga",
				"Pine",
				"Red Hill",
				"Fontana"
		};
		
		List<Integer> matchingRups = new ArrayList<>();
		for (int rupIndex : rupSet.getRupturesForParentSection(targetParentID)) {
			double mag = rupSet.getMagForRup(rupIndex);
			if (magRange.contains(mag)) {
				boolean avoids = true;
				List<String> myParents = new ArrayList<>();
				int prevParent = -1;
				for (FaultSection sect : rupSet.getFaultSectionDataForRupture(rupIndex)) {
					int parent = sect.getParentSectionId();
					if (parent != prevParent) {
						String parentName = sect.getParentSectionName();
						for (String avoid : avoidParents) {
							if (parentName.contains(avoid)) {
								avoids = false;
								break;
							}
						}
						if (!avoids)
							break;
						myParents.add(parentName);
						prevParent = parent;
					}
				}
				if (avoids) {
					matchingRups.add(rupIndex);
					System.out.println("Match: M"+(float)mag+" on: "+Joiner.on("; ").join(myParents));
				}
			}
		}
		System.out.println("Identified "+matchingRups.size()+" matching ruptures");
		
		System.out.println("Calculating conditional curve");
		double sumMatchRate = 0d;
		DiscretizedFunc condFunc = xVals.deepClone();
		for (int i=0; i<condFunc.size(); i++)
			condFunc.set(i, 0d);
		gmm.setSite(site);
		for (int rupIndex : matchingRups) {
			RuptureSurface surf = rupSet.getSurfaceForRupture(rupIndex, 1d);
			EqkRupture rup = new EqkRupture(targetMag, rupSet.getAveRakeForRup(rupIndex), surf, null);
			gmm.setEqkRupture(rup);
			gmm.getExceedProbabilities(logXVals);
			double rate = sol.getRateForRup(rupIndex);
			sumMatchRate += rate;
			for (int i=0; i<condFunc.size(); i++)
				condFunc.set(i, condFunc.getY(i) + rate*logXVals.getY(i));
		}
		condFunc.scale(1d/sumMatchRate);
		System.out.println("DONE calculating conditional curve");
		
		// now full hazard curve
		System.out.println("Calculating full hazard curve");
		calc.getHazardCurve(logXVals, site, gmm, erf);
		System.out.println("DONE calculating full hazard curve");
		
		DiscretizedFunc hazardCurve = xVals.deepClone();
		for (int i=0; i<xVals.size(); i++)
			hazardCurve.set(i, logXVals.getY(i));
		
		// convert to percent
		condFunc.scale(100d);
		hazardCurve.scale(100d);
		
		HeadlessGraphPanel gp = PlotUtils.initScreenHeadless();
		PlotPreferences prefs = gp.getPlotPrefs();
		prefs.scaleFontSizes(1.5);
		
		for (boolean condOnly : new boolean[] {true,false}) {
			List<DiscretizedFunc> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			condFunc.setName(targetName);
			funcs.add(condFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLACK));
			
			if (!condOnly) {
				hazardCurve.setName(hazardCurveName);
				funcs.add(hazardCurve);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Colors.tab_red));
			}
			
			PlotSpec plot = new PlotSpec(funcs, chars, siteName+" Hazard Curves", xAxisName, "Probability (%)");
			plot.setLegendInset(true);
			
			gp.drawGraphPanel(plot, false, false, new org.jfree.data.Range(0d, 1d), new org.jfree.data.Range(0d, 100d));
			
			String prefix = "hazard_curve_example";
			if (condOnly)
				prefix += "_cond_only";
			PlotUtils.writePlots(outputDir, prefix, gp, 1000, 800, true, true, false);
		}
		
		
		
//		for (int rupIndex=0; rupIndex<rupSet.getNumRuptures(); rupIndex++) {
//			for (FaultSectionUtils se)
//		}
	}

}
