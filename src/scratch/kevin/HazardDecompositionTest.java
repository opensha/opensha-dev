package scratch.kevin;

import java.awt.Color;
import java.util.Collections;
import java.util.List;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.siteData.OrderedSiteDataProviderList;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.hazardMap.HazardCurveSetCalculator;
import org.opensha.sha.cybershake.db.MeanUCERF2_ToDB;
import org.opensha.sha.cybershake.plot.HazardCurvePlotter;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ERFTestSubset;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.util.SiteTranslator;
import org.opensha.sra.rtgm.RTGM;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public class HazardDecompositionTest {

	public static void main(String[] args) {
		AbstractERF erf = MeanUCERF2_ToDB.createUCERF2ERF();
		ScalarIMR imr = AttenRelRef.CB_2014.instance(null);
		imr.setParamDefaults();
		imr.setIntensityMeasure(SA_Param.NAME);
		SA_Param.setPeriodInSA_Param(imr.getIntensityMeasure(), 3d);
		String xAxisLabel = "3sec SA";
		
		Location loc = new Location(34.05204, -118.25713); // LADT
		int numToInclude = 10;
		
		Site site = new Site(loc);
		
		int velModelID = 5;
		OrderedSiteDataProviderList providers = HazardCurvePlotter.createProviders(velModelID);
		
		SiteTranslator trans = new SiteTranslator();
		trans.setAllSiteParams(imr, providers.getBestAvailableData(loc));
		site.addParameterList(imr.getSiteParams());
		
		HazardCurveCalculator calc = new HazardCurveCalculator();
		
		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(imr.getIntensityMeasure());
		DiscretizedFunc logXVals = HazardCurveSetCalculator.getLogFunction(xVals.deepClone());
		
		DiscretizedFunc totalHazard = calc.getAnnualizedRates(
				HazardCurveSetCalculator.unLogFunction(xVals,
						calc.getHazardCurve(logXVals.deepClone(), site, imr, erf)), 1d);
		totalHazard.setName("Total Hazard");
		
		List<DiscretizedFunc> sourceFuncs = Lists.newArrayList();
		List<Double> sourceRTGMContributions = Lists.newArrayList();
		
		double totRTGM = calcRTGM(totalHazard);
		
		int numAboveZero = 0;
		
		for (int sourceID=0; sourceID<erf.getNumSources(); sourceID++) {
			ProbEqkSource source = erf.getSource(sourceID);
			if (source.getMinDistance(site) > 200d) {
				sourceFuncs.add(null);
				sourceRTGMContributions.add(0d);
				continue;
			}
			
			System.out.println("Source "+sourceID);
			ERFTestSubset subset = new ERFTestSubset(erf);
			subset.includeSource(sourceID);
			subset.updateForecast();
			Preconditions.checkState(subset.getNumSources() == 1);
			DiscretizedFunc srcHazard = calc.getAnnualizedRates(
					HazardCurveSetCalculator.unLogFunction(xVals,
							calc.getHazardCurve(logXVals.deepClone(), site, imr, subset)), 1d);
			sourceFuncs.add(srcHazard);
			
			// now all sources except this one for disagg
			subset = new ERFTestSubset(erf);
			subset.includeAllExcept(sourceID);
			subset.updateForecast();
			Preconditions.checkState(subset.getNumSources() == erf.getNumSources()-1);
			DiscretizedFunc srcWithoutHazard = calc.getAnnualizedRates(
					HazardCurveSetCalculator.unLogFunction(xVals,
							calc.getHazardCurve(logXVals.deepClone(), site, imr, subset)), 1d);
			
			double withoutRTGM = calcRTGM(srcWithoutHazard);
			double deltaRTGM = totRTGM - withoutRTGM;
			if (deltaRTGM > 0 )
				numAboveZero++;
			Preconditions.checkState(deltaRTGM >= 0);
			sourceRTGMContributions.add(deltaRTGM);
			
			srcHazard.setName(source.getName());
			srcHazard.setInfo("RTGM Contribution: "+(float)totRTGM+" - "+(float)withoutRTGM+" = "+(float)deltaRTGM);
			System.out.println(source.getName()+" "+srcHazard.getInfo());
		}
		
		System.out.println("Num above zero: "+numAboveZero);
		
		List<ComparablePairing<Double, DiscretizedFunc>> pairings =
				ComparablePairing.build(sourceRTGMContributions, sourceFuncs);
		Collections.sort(pairings);
		Collections.reverse(pairings);
		
		List<DiscretizedFunc> funcs = Lists.newArrayList();
		List<Color> colors = GraphWindow.generateDefaultColors();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		funcs.add(totalHazard);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, colors.get((funcs.size()-1) % colors.size())));
		
		for (int i=0; i<numToInclude; i++) {
			funcs.add(pairings.get(i).getData());
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, colors.get((funcs.size()-1) % colors.size())));
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Hazard Decomposition", xAxisLabel, "Annual Rate Of Exceedance");
		spec.setLegendVisible(true);
		
		new GraphWindow(spec);
	}
	
	private static double calcRTGM(DiscretizedFunc curve) {
		RTGM calc = RTGM.create(curve, null, null);
		try {
			calc.call();
		} catch (RuntimeException e) {
			System.err.println("RTGM Calc failed for Hazard Curve:\n"+curve);
			System.err.flush();
			throw e;
		}
		double rtgm = calc.get();
		Preconditions.checkState(rtgm > 0, "RTGM is not positive");
		return rtgm;
	}

}
