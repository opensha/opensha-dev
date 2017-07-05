package scratch.kevin.cybershake.ucerf3;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.dom4j.DocumentException;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.siteData.OrderedSiteDataProviderList;
import org.opensha.commons.data.siteData.SiteDataValue;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.param.Parameter;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.hazardMap.HazardCurveSetCalculator;
import org.opensha.sha.cybershake.calc.HazardCurveComputation;
import org.opensha.sha.cybershake.calc.UCERF2_AleatoryMagVarRemovalMod;
import org.opensha.sha.cybershake.db.CybershakeIM;
import org.opensha.sha.cybershake.db.CybershakeRun;
import org.opensha.sha.cybershake.db.CybershakeSite;
import org.opensha.sha.cybershake.db.CybershakeSiteInfo2DB;
import org.opensha.sha.cybershake.db.Cybershake_OpenSHA_DBApplication;
import org.opensha.sha.cybershake.db.DBAccess;
import org.opensha.sha.cybershake.db.HazardCurve2DB;
import org.opensha.sha.cybershake.db.MeanUCERF2_ToDB;
import org.opensha.sha.cybershake.db.Runs2DB;
import org.opensha.sha.cybershake.plot.HazardCurvePlotter;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.param.AleatoryMagAreaStdDevParam;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.attenRelImpl.NGAWest_2014_Averaged_AttenRel;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.util.SiteTranslator;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.utils.FaultSystemIO;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

/**
 * This will show the influence of Aleatory Variability on UCERF2 results to help inform our
 * decision to include/exclude it for UCERF3 calculations
 * @author kevin
 *
 */
public class AleatoryMagVariabilitiyTest {

	public static void main(String[] args) throws IOException, DocumentException {
		ERF erf = MeanUCERF2_ToDB.createUCERF2ERF();
		UCERF2_AleatoryMagVarRemovalMod mod = new UCERF2_AleatoryMagVarRemovalMod(erf);
		
		File plotDir = new File("/home/kevin/CyberShake/ucerf3/aleatory_test_ucerf2");
		
		int imTypeID = 21;
		int erfID = 35;
		int sgtVarID = 8;
		int rupVarScenID = 4;
		int velModelID = 5;
		
		List<CybershakeSite> curveSites = Lists.newArrayList();
		DBAccess db = Cybershake_OpenSHA_DBApplication.getDB();
		CybershakeSiteInfo2DB sites2db = new CybershakeSiteInfo2DB(db);
		curveSites.add(sites2db.getSiteFromDB("STNI"));
		curveSites.add(sites2db.getSiteFromDB("SBSM"));
		curveSites.add(sites2db.getSiteFromDB("PAS"));
		curveSites.add(sites2db.getSiteFromDB("FIL"));
		Runs2DB runs2db = new Runs2DB(db);
		List<CybershakeRun> runs = Lists.newArrayList();
		
		for (CybershakeSite site : curveSites) {
			CybershakeRun run = runs2db.getLatestRun(site.id, erfID, sgtVarID, rupVarScenID, velModelID, null, null, null, null);
			Preconditions.checkNotNull(run);
			runs.add(run);
		}
		CybershakeIM imType = new HazardCurve2DB(db).getIMFromID(imTypeID);
		
		List<Double> xVals = Lists.newArrayList();
		for (Point2D pt : new IMT_Info().getDefaultHazardCurve(SA_Param.NAME))
			xVals.add(pt.getX());
		
		List<DiscretizedFunc> regularCurves = Lists.newArrayList();
		List<DiscretizedFunc> noAleatoryCurves = Lists.newArrayList();
		
		HazardCurveComputation calc = new HazardCurveComputation(db);
		// first do normal calc
		System.out.println("Calculating original curves");
		for (int i=0; i<curveSites.size(); i++)
			regularCurves.add(calc.computeHazardCurve(xVals, runs.get(i), imType));
		// now no aleatory calc
		System.out.println("Calculating modified curves");
		calc.setRupProbModifier(mod);
		for (int i=0; i<curveSites.size(); i++)
			noAleatoryCurves.add(calc.computeHazardCurve(xVals, runs.get(i), imType));
		db.destroy();
//		System.out.println("Stripped aleatory variability from "+mod.aleatorySources.size()+" sources");
//		System.out.println("Raymond? "+mod.aleatorySources.contains(247));
		
		// now do some sanity checks
		int[] testIDs = { 90, 247, 246, 245 };
		for (int testID : testIDs) {
			ProbEqkSource testSource = erf.getSource(testID);
			System.out.println("Probs for "+testSource.getName());
			for (int rupID=0; rupID<testSource.getNumRuptures(); rupID++) {
				ProbEqkRupture rup = testSource.getRupture(rupID);
				double origProb = rup.getProbability();
				double modProb = mod.getModifiedProb(testID, rupID, origProb);
				System.out.println("\t"+(float)rup.getMag()+"\t"+(float)origProb+"\t"+(float)modProb);
			}
			System.out.println();
		}
		
		for (int i=0; i<curveSites.size(); i++) {
			List<DiscretizedFunc> funcs = Lists.newArrayList();
			List<PlotCurveCharacterstics> chars = Lists.newArrayList();
			
			DiscretizedFunc origCurve = regularCurves.get(i);
			origCurve.setName("Original");
			funcs.add(origCurve);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
			DiscretizedFunc newCurve = noAleatoryCurves.get(i);
			newCurve.setName("No Aleatory Mag Variability");
			funcs.add(newCurve);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
			
			String siteName = curveSites.get(i).short_name;
			
			PlotSpec spec = new PlotSpec(funcs, chars, siteName, "3sec SA (g)", "Exceed. Prob");
			spec.setLegendVisible(true);
			GraphWindow gw = new GraphWindow(spec);
			gw.getGraphWidget().setBackgroundColor(Color.WHITE);
			gw.setXLog(true);
			gw.setYLog(true);
			gw.setX_AxisRange(1e-2, 3);
			gw.setY_AxisRange(1e-8, 1e-1);
			gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
			if (plotDir != null) {
				try {
					gw.saveAsPNG(new File(plotDir, "curve_"+siteName+"_comparison.png").getAbsolutePath());
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		
		// now do UCERF3 GMPE tests
		FaultSystemSolution inputSol = FaultSystemIO.loadSol(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		NGAWest_2014_Averaged_AttenRel gmpe = new NGAWest_2014_Averaged_AttenRel(null, false);
		
		erf = new FaultSystemSolutionERF(inputSol);
		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
		erf.getTimeSpan().setDuration(1d);
		erf.updateForecast();
		
		// set to 3s SA
		gmpe.setParamDefaults();
		gmpe.setIntensityMeasure(SA_Param.NAME);
		SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), 3d);
		HazardCurveCalculator gmpeCalc = new HazardCurveCalculator();
		DiscretizedFunc xValFunc = new IMT_Info().getDefaultHazardCurve(gmpe.getIntensityMeasure());

		List<DiscretizedFunc> origCurves = Lists.newArrayList();
		List<DiscretizedFunc> newCurves = Lists.newArrayList();

		OrderedSiteDataProviderList provs = HazardCurvePlotter.createProviders(velModelID);
		SiteTranslator trans = new SiteTranslator();
		List<Site> sites = Lists.newArrayList();
		// create sites with site data
		for (CybershakeSite csSite : curveSites) {
			Site site = new Site(csSite.createLocation());

			ArrayList<SiteDataValue<?>> datas = provs.getBestAvailableData(site.getLocation());

			for (Parameter<?> param : gmpe.getSiteParams()) {
				param = (Parameter<?>) param.clone();
				trans.setParameterValue(param, datas);
				site.addParameter(param);
			}
			sites.add(site);
		}
		
		System.out.println("Calculating original curves");
		for (Site site : sites) {
			DiscretizedFunc curve = HazardCurveSetCalculator.getLogFunction(xValFunc);
			gmpeCalc.getHazardCurve(curve, site, gmpe, erf);
			origCurves.add(HazardCurveSetCalculator.unLogFunction(xValFunc, curve));
		}
		
		erf = new FaultSystemSolutionERF(inputSol);
		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
		erf.setParameter(AleatoryMagAreaStdDevParam.NAME, 0.12);
		erf.getTimeSpan().setDuration(1d);
		erf.updateForecast();
		
		System.out.println("Calculating new curves");
		for (Site site : sites) {
			DiscretizedFunc curve = HazardCurveSetCalculator.getLogFunction(xValFunc);
			gmpeCalc.getHazardCurve(curve, site, gmpe, erf);
			newCurves.add(HazardCurveSetCalculator.unLogFunction(xValFunc, curve));
		}
		
		for (int i=0; i<sites.size(); i++) {
			List<DiscretizedFunc> funcs = Lists.newArrayList();
			List<PlotCurveCharacterstics> chars = Lists.newArrayList();
			
			DiscretizedFunc origCurve = origCurves.get(i);
			origCurve.setName("Original");
			funcs.add(origCurve);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
			DiscretizedFunc newCurve = newCurves.get(i);
			newCurve.setName("No Aleatory Mag Variability");
			funcs.add(newCurve);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
			
			String siteName = curveSites.get(i).short_name;
			
			PlotSpec spec = new PlotSpec(funcs, chars, siteName+" GMPE Curves", "3sec SA (g)", "Exceed. Prob");
			spec.setLegendVisible(true);
			GraphWindow gw = new GraphWindow(spec);
			gw.getGraphWidget().setBackgroundColor(Color.WHITE);
			gw.setXLog(true);
			gw.setYLog(true);
			gw.setX_AxisRange(1e-2, 3);
			gw.setY_AxisRange(1e-8, 1e-1);
			gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
			if (plotDir != null) {
				try {
					gw.saveAsPNG(new File(plotDir, "gmpe_curve_"+siteName+"_comparison.png").getAbsolutePath());
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}

}
