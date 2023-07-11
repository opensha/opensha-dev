package scratch.kevin.nshm23;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultGridAssociations;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.NucleationRatePlot;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.faultSurface.CompoundSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.util.NEHRP_TestCity;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

import gov.usgs.earthquake.nshmp.model.HazardModel;
import gov.usgs.earthquake.nshmp.model.NshmErf;
import gov.usgs.earthquake.nshmp.model.NshmSurface;
import scratch.UCERF3.erf.FaultSystemSolutionERF;

public class SingleSiteHazardAndDataComparisonPageGen {
	
	private static AttenRelRef GMPE_DEFAULT = AttenRelRef.ASK_2014;
	
//	private static final double[] PERIODS_DEFAULT = { 0d, 0.2d, 1d };
	public static final double[] PERIODS_DEFAULT = { 0d, 1d };

	public static void main(String[] args) throws IOException {
		if (args.length == 1 && args[0].equals("--hardcoded")) {
			String argStr = "--input-file /home/kevin/OpenSHA/nshm23/batch_inversions/"
					+ "2023_01_17-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
					+ "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip";
			argStr += " --comp-nshm-erf /data/kevin/nshm23/nshmp-haz-models/nshm-conus-5.2.0";
			argStr += " --comp-name NSHM18";
//			argStr += " --nehrp-site santabarbara";
//			argStr += " --output-dir /tmp/curves_test_sb";
//			argStr += " --nehrp-site reno";
//			argStr += " --output-dir /tmp/curves_test_reno";
//			argStr += " --nehrp-site saltlakecity";
//			argStr += " --output-dir /tmp/curves_test_slc";
			argStr += " --nehrp-site seattle";
			argStr += " --output-dir /tmp/curves_test_seattle";
			args = argStr.split(" ");
		}
		Options ops = createOptions();
		
		CommandLine cmd = FaultSysTools.parseOptions(ops, args, SingleSiteHazardAndDataComparisonPageGen.class);
		
		File outputDir = new File(cmd.getOptionValue("output-dir"));
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		Location loc;
		String locName;
		
		if (ops.hasOption("nehrp-site")) {
			Preconditions.checkState(!cmd.hasOption("location"), "Can't supply --nehrp-site and --location");
			String rawName = cmd.getOptionValue("nehrp-site");
			String siteName = rawName.toLowerCase().replace(" ", "").replace("_", "");
			NEHRP_TestCity match = null;
			for (NEHRP_TestCity site : NEHRP_TestCity.values()) {
				String test = site.toString().toLowerCase().replace(" ", "").replace("_", "");
				if (test.equals(siteName)) {
					match = site;
					break;
				}
			}
			Preconditions.checkNotNull(match, "No NEHRP site found matching name: %s", rawName);
			
			loc = match.location();
			locName = match.toString();
		} else if (cmd.hasOption("location")) {
			String locStr = cmd.getOptionValue("location");
			Preconditions.checkState(locStr.contains(","), "Location should be specified as lat,lon: %s", locStr);
			String[] split = locStr.trim().split(",");
			Preconditions.checkState(split.length == 2, "Location should be specified as lat,lon: %s", locStr);
			loc = new Location(Double.parseDouble(split[0]), Double.parseDouble(split[1]));
			locName = locStr.trim();
		} else {
			throw new IllegalArgumentException("Must specify --location or --nehrp-site");
		}
		
		if (cmd.hasOption("name"))
			locName = cmd.getOptionValue("name");
		
		File inputFile = new File(cmd.getOptionValue("input-file"));
		Preconditions.checkState(inputFile.exists());
		
		FaultSystemSolution sol = FaultSystemSolution.load(inputFile);
		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		erf.getTimeSpan().setDuration(1d);
		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
		
		NshmErf compGridERF = null;
		NshmErf compFaultERF = null;
		FaultSystemSolution compSol = null;
		FaultSystemSolutionERF compSolERF = null;
		String compName = cmd.hasOption("comp-name") ? cmd.getOptionValue("comp-name") : "Comparison";
		boolean hasComp = false;
		
		double maxMag = sol.getRupSet().getMaxMag();
		
		if (cmd.hasOption("comp-nshm-erf")) {
			Preconditions.checkState(!cmd.hasOption("comp-sol"), "Can't supply both --comp-nshm-erf and --comp-sol");
			Path path = Path.of(cmd.getOptionValue("comp-nshm-erf"));
			HazardModel model = HazardModel.load(path);
			
			// TODO subduction option?
			Set<TectonicRegionType> trts = EnumSet.of(TectonicRegionType.ACTIVE_SHALLOW,
					TectonicRegionType.STABLE_SHALLOW);
			compFaultERF = new NshmErf(model, trts, IncludeBackgroundOption.EXCLUDE);
			compFaultERF.getTimeSpan().setDuration(1d);
			compFaultERF.updateForecast();
			compGridERF = new NshmErf(model, trts, IncludeBackgroundOption.ONLY);
			compGridERF.getTimeSpan().setDuration(1d);
			compGridERF.updateForecast();
			
			for (ProbEqkSource source : compFaultERF)
				for (ProbEqkRupture rup : source)
					maxMag = Math.max(maxMag, rup.getMag());
			for (ProbEqkSource source : compGridERF)
				for (ProbEqkRupture rup : source)
					maxMag = Math.max(maxMag, rup.getMag());
			hasComp = true;
		} else if (cmd.hasOption("comp-sol")) {
			compSol = FaultSystemSolution.load(new File(cmd.getOptionValue("comp-sol")));
			compSolERF = new FaultSystemSolutionERF(sol);
			compSolERF.getTimeSpan().setDuration(1d);
			compSolERF.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
			
			maxMag = Math.max(maxMag, compSol.getRupSet().getMaxMag());
			hasComp = true;
		}
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(maxMag);
		
		HazardCurveCalculator curveCalc = new HazardCurveCalculator();
		
		AttenRelRef gmpeRef = GMPE_DEFAULT;
		double[] periods = PERIODS_DEFAULT;
		if (cmd.hasOption("gmpe"))
			gmpeRef = AttenRelRef.valueOf(cmd.getOptionValue("gmpe"));
		
		if (cmd.hasOption("periods")) {
			List<Double> periodsList = new ArrayList<>();
			String periodsStr = cmd.getOptionValue("periods");
			if (periodsStr.contains(",")) {
				String[] split = periodsStr.split(",");
				for (String str : split)
					periodsList.add(Double.parseDouble(str));
			} else {
				periodsList.add(Double.parseDouble(periodsStr));
			}
			periods = Doubles.toArray(periodsList);
		}
		
		DiscretizedFunc[] perXVals = new DiscretizedFunc[periods.length];
		DiscretizedFunc[] perLogXVals = new DiscretizedFunc[periods.length];
		IMT_Info imtInfo = new IMT_Info();
		for (int p=0; p<periods.length; p++) {
			if (periods[p] > 0d)
				perXVals[p] = imtInfo.getDefaultHazardCurve(SA_Param.NAME);
			else if (periods[p] == 0d)
				perXVals[p] = imtInfo.getDefaultHazardCurve(PGA_Param.NAME);
			else
				throw new IllegalArgumentException("Periods must be >= 0");
			perLogXVals[p] = new ArbitrarilyDiscretizedFunc();
			for (Point2D pt : perXVals[p])
				perLogXVals[p].set(Math.log(pt.getX()), 0d);
		}
		
		double gridSpacing = 0.1d;
		
		double plotBufferDist = 100d;
		Location above = LocationUtils.location(loc, 0, plotBufferDist);
		Location right = LocationUtils.location(loc, 0.5*Math.PI, plotBufferDist);
		Location below = LocationUtils.location(loc, Math.PI, plotBufferDist);
		Location left = LocationUtils.location(loc, 1.5*Math.PI, plotBufferDist);
		Region plotReg = new Region(new Location(below.getLatitude(), left.getLongitude()),
				new Location(above.getLatitude(), right.getLongitude()));
		GriddedRegion plotGridReg = new GriddedRegion(plotReg, gridSpacing, GriddedRegion.ANCHOR_0_0);
		double halfSpacing = gridSpacing/2d;
		plotReg = new Region(new Location(plotGridReg.getMinGridLat()-halfSpacing, plotGridReg.getMinGridLon()-halfSpacing),
				new Location(plotGridReg.getMaxGridLat()+halfSpacing, plotGridReg.getMaxGridLon()+halfSpacing));
		
		GeographicMapMaker mapMaker = new RupSetMapMaker(sol.getRupSet(), plotReg);
		mapMaker.setWritePDFs(false);
		mapMaker.plotScatters(List.of(loc), Color.CYAN);
		mapMaker.setScatterSymbol(PlotSymbol.FILLED_INV_TRIANGLE, 10f, PlotSymbol.INV_TRIANGLE, Color.BLACK);
		mapMaker.setSectOutlineChar(null);
		
		double[] magThresholds = { 5d, 6.5d, 7d, 7.5d };
		double[] dists = { 200d, 100d, 50d };
		GriddedRegion[] distRegions = new GriddedRegion[dists.length];
		for (int d=0; d<dists.length; d++)
			distRegions[d] = new GriddedRegion(new Region(loc, dists[d]), gridSpacing, GriddedRegion.ANCHOR_0_0);
		
		System.out.println("Calculating regional perticipation maps and MFDs for input solution");
		RegionalParticipationResult solFaultPartic = calcFSSFaultPartic(sol, plotGridReg, magThresholds, refMFD);
		RegionalParticipationResult solGridPartic = calcFSSGriddedPartic(sol, plotGridReg, magThresholds, refMFD);
		RegionalParticipationResult[] distFaultPartics = new RegionalParticipationResult[dists.length];
		IncrementalMagFreqDist[] distFaultMFDs = new IncrementalMagFreqDist[dists.length];
		EvenlyDiscretizedFunc[] distFaultCmlMFDs = new EvenlyDiscretizedFunc[dists.length];
		RegionalParticipationResult[] distGridPartics = new RegionalParticipationResult[dists.length];
		IncrementalMagFreqDist[] distGridMFDs = new IncrementalMagFreqDist[dists.length];
		EvenlyDiscretizedFunc[] distGridCmlMFDs = new EvenlyDiscretizedFunc[dists.length];
		for (int d=0; d<dists.length; d++) {
			distFaultPartics[d] = calcFSSFaultPartic(sol, distRegions[d], new double[0], refMFD);
			distFaultMFDs[d] = distFaultPartics[d].particMFD;
			distFaultCmlMFDs[d] = distFaultMFDs[d].getCumRateDistWithOffset();
			distGridPartics[d] = calcFSSGriddedPartic(sol, distRegions[d], new double[0], refMFD);
			distGridMFDs[d] = distGridPartics[d].particMFD;
			distGridCmlMFDs[d] = distGridMFDs[d].getCumRateDistWithOffset();
		}
		
		RegionalParticipationResult compFaultPartic = null;
		RegionalParticipationResult compGridPartic = null;
		RegionalParticipationResult[] compDistFaultPartics = null;
		IncrementalMagFreqDist[] compDistFaultMFDs = null;
		EvenlyDiscretizedFunc[] compDistFaultCmlMFDs = null;
		RegionalParticipationResult[] compDistGridPartics = null;
		IncrementalMagFreqDist[] compDistGridMFDs = null;
		EvenlyDiscretizedFunc[] compDistGridCmlMFDs = null;
		if (hasComp) {
			System.out.println("Calculating comparison participation maps and MFDs");
			compDistFaultPartics = new RegionalParticipationResult[dists.length];
			compDistFaultMFDs = new IncrementalMagFreqDist[dists.length];
			compDistFaultCmlMFDs = new EvenlyDiscretizedFunc[dists.length];
			compDistGridPartics = new RegionalParticipationResult[dists.length];
			compDistGridMFDs = new IncrementalMagFreqDist[dists.length];
			compDistGridCmlMFDs = new EvenlyDiscretizedFunc[dists.length];
			if (compFaultERF != null) {
				compFaultPartic = calcNshmErfPartic(compFaultERF, plotGridReg, magThresholds, refMFD);
				compGridPartic = calcNshmErfPartic(compGridERF, plotGridReg, magThresholds, refMFD);
				for (int d=0; d<dists.length; d++) {
					compDistFaultPartics[d] = calcNshmErfPartic(compFaultERF, distRegions[d], new double[0], refMFD);
					compDistGridPartics[d] = calcNshmErfPartic(compGridERF, distRegions[d], new double[0], refMFD);
				}
			} else if (compSol != null) {
				compFaultPartic = calcFSSFaultPartic(compSol, plotGridReg, magThresholds, refMFD);
				compGridPartic = calcFSSGriddedPartic(compSol, plotGridReg, magThresholds, refMFD);
				for (int d=0; d<dists.length; d++) {
					compDistFaultPartics[d] = calcFSSFaultPartic(compSol, distRegions[d], new double[0], refMFD);
					compDistGridPartics[d] = calcFSSGriddedPartic(compSol, distRegions[d], new double[0], refMFD);
				}
			}
			for (int d=0; d<dists.length; d++) {
				compDistFaultMFDs[d] = compDistFaultPartics[d].particMFD;
				compDistFaultCmlMFDs[d] = compDistFaultMFDs[d].getCumRateDistWithOffset();
				compDistGridMFDs[d] = compDistGridPartics[d].particMFD;
				compDistGridCmlMFDs[d] = compDistGridMFDs[d].getCumRateDistWithOffset();
			}
		}
		
		List<String> lines = new ArrayList<>();
		
		if (hasComp)
			lines.add("# Hazard And Model Comparisons for "+locName);
		else
			lines.add("# Hazard And Model Summary for "+locName);
		lines.add("");
		
		TableBuilder table = MarkdownUtils.tableBuilder();
		
		table.addLine("", "Location and GMPE Parameters");
		table.addLine("__Location__", (float)loc.getLatitude()+", "+(float)loc.getLongitude());
		table.addLine("__GMP Used__", gmpeRef.getName());
		
		ScalarIMR gmpe = gmpeRef.get();
		gmpe.setParamDefaults();
		Site site = new Site(loc);
		for (Parameter<?> param : gmpe.getSiteParams()) {
			site.addParameter(param);
			String valStr;
			if (param.getValue() instanceof Double && param.getValue() != null)
				valStr = ((Double)param.getValue()).floatValue()+"";
			else
				valStr = param.getValue()+"";
			table.addLine("__"+param.getName()+"__", valStr);
		}
		
		lines.addAll(table.build());
		lines.add("");
		
		int tocIndex = lines.size();
		String topLink = "*[(top)](#table-of-contents)*";
		
		lines.add("## Hazard Curves");
		lines.add(topLink); lines.add("");
		
		table = MarkdownUtils.tableBuilder();
		table.initNewLine();
		for (int p=0; p<periods.length; p++)
			table.addColumn(periods[p] == 0d ? "PGA" : (float)periods[p]+"s SA");
		table.finalizeLine();
		
		Color mainColor = Color.RED;
		Color compColor = Color.BLUE;
		
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		DiscretizedFunc[] fullCurves = new DiscretizedFunc[periods.length];
		DiscretizedFunc[] compFullCurves = hasComp ? new DiscretizedFunc[periods.length] : null;
		
		ReturnPeriods[] rps = ReturnPeriods.values();
		
		table.initNewLine();
		for (int p=0; p<periods.length; p++) {
			String xLabel;
			String prefix;
			if (periods[p] == 0d) {
				gmpe.setIntensityMeasure(PGA_Param.NAME);
				xLabel = "Peak Ground Acceleration (g)";
				prefix = "curves_pga";
			} else {
				gmpe.setIntensityMeasure(SA_Param.NAME);
				SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), periods[p]);
				xLabel = (float)periods[p]+"s Spectral Acceleration (g)";
				prefix = "curves_"+(float)periods[p]+"_sa";
			}
			
			DiscretizedFunc logXVals = perLogXVals[p];
			DiscretizedFunc xVals = perXVals[p];
			
			System.out.println("Calculating curves for period="+(float)periods[p]);
			erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
			erf.updateForecast();
			curveCalc.getHazardCurve(logXVals, site, gmpe, erf);
			DiscretizedFunc faultCurve = toLinearCurve(logXVals, xVals);
			erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.ONLY);
			erf.updateForecast();
			curveCalc.getHazardCurve(logXVals, site, gmpe, erf);
			DiscretizedFunc gridCurve = toLinearCurve(logXVals, xVals);
			
			DiscretizedFunc sumCurve = sumCurve(faultCurve, gridCurve);
			fullCurves[p] = sumCurve;
			
			DiscretizedFunc compFaultCurve = null;
			DiscretizedFunc compGridCurve = null;
			DiscretizedFunc compSumCurve = null;
			
			if (compFaultERF != null) {
				System.out.println("Calculating comparison curves for period="+(float)periods[p]);
				curveCalc.getHazardCurve(logXVals, site, gmpe, compFaultERF);
				compFaultCurve = toLinearCurve(logXVals, xVals);
				curveCalc.getHazardCurve(logXVals, site, gmpe, compGridERF);
				compGridCurve = toLinearCurve(logXVals, xVals);
				compSumCurve = sumCurve(compFaultCurve, compGridCurve);
			} else if (compSolERF != null) {
				System.out.println("Calculating comparison curves for period="+(float)periods[p]);
				compSolERF.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
				compSolERF.updateForecast();
				curveCalc.getHazardCurve(logXVals, site, gmpe, compSolERF);
				compFaultCurve = toLinearCurve(logXVals, xVals);
				compSolERF.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.ONLY);
				compSolERF.updateForecast();
				curveCalc.getHazardCurve(logXVals, site, gmpe, compSolERF);
				compGridCurve = toLinearCurve(logXVals, xVals);
				compSumCurve = sumCurve(compFaultCurve, compGridCurve);
			}
			if (hasComp)
				compFullCurves[p] = compSumCurve;
			
			List<DiscretizedFunc> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			sumCurve.setName("Full Hazard Curve");
			funcs.add(sumCurve);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, mainColor));
			
			faultCurve.setName("Faults Only");
			funcs.add(faultCurve);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, mainColor));
			
			gridCurve.setName("Gridded Only");
			funcs.add(gridCurve);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, mainColor));
			
			if (compSumCurve != null) {
				compSumCurve.setName(compName+" Full");
				funcs.add(compSumCurve);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, compColor));
				
				compFaultCurve.setName("Faults Only");
				funcs.add(compFaultCurve);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, compColor));
				
				compGridCurve.setName("Gridded Only");
				funcs.add(compGridCurve);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, compColor));
			}
			
			Range xRange = new Range(1e-2, 1e1);
			
			for (ReturnPeriods rp : rps) {
				DiscretizedFunc rpFunc = new ArbitrarilyDiscretizedFunc();
				rpFunc.set(xRange.getLowerBound(), rp.oneYearProb);
				rpFunc.set(xRange.getUpperBound(), rp.oneYearProb);
				funcs.add(rpFunc);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.GRAY));
			}
			
			PlotSpec spec = new PlotSpec(funcs, chars, locName+" Hazard Curves",
					xLabel, "Annual Probability of Exceedance");
			spec.setLegendInset(true);
			
			for (ReturnPeriods rp : rps) {
				XYTextAnnotation ann = new XYTextAnnotation(rp.label, xRange.getLowerBound(), rp.oneYearProb);
				ann.setTextAnchor(TextAnchor.BASELINE_LEFT);
				ann.setPaint(Color.DARK_GRAY);
				ann.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 18));;
				spec.addPlotAnnotation(ann);
			}
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
			gp.drawGraphPanel(spec, true, true, xRange, new Range(1e-6, 1e0));
			
			PlotUtils.writePlots(resourcesDir, prefix, gp, 1000, 850, true, true, true);
			
			table.addColumn("![Hazard Curves]("+resourcesDir.getName()+"/"+prefix+".png)");
		}
		table.finalizeLine();
		
		lines.addAll(table.wrap(2, 0).build());
		lines.add("");
		
		table = MarkdownUtils.tableBuilder();
		table.initNewLine();
		table.addColumns("Return Period");
		for (double period : periods) {
			if (period == 0d)
				table.addColumn("PGA (g)");
			else
				table.addColumn((float)period+"s SA (g)");
			if (hasComp)
				table.addColumns(compName, "% Change");
		}
		table.finalizeLine();
		
		for (ReturnPeriods rp : rps) {
			table.initNewLine();
			
			table.addColumn(rp.label);
			for (int p=0; p<periods.length; p++) {
				double iml = calcCurveIML(rp, fullCurves[p]);
				table.addColumn((float)iml);
				if (hasComp) {
					double compIML = calcCurveIML(rp, compFullCurves[p]);
					table.addColumn((float)compIML);
					double pDiff = (iml-compIML)/compIML;
					table.addColumn((iml>compIML ? "+" : "")+pDF.format(pDiff));
				}
			}
			
			table.finalizeLine();
		}
		
		lines.addAll(table.build());
		lines.add("");
		
		lines.add("## Model Comparison Summary");
		lines.add(topLink); lines.add("");
		
		lines.add("");
		lines.add("");
		
		if (hasComp) {
			// summary table
			table = MarkdownUtils.tableBuilder();
			
			// will invert this table
			
			table.initNewLine().addColumn("");
			table.addColumn("Total Moment Rate");
			table.addColumn("Fault Moment Rate");
			table.addColumn("Gridded Moment Rate");
			for (double mag : magThresholds) {
				String magStr = oDF.format(mag);
				table.addColumn("Total Rate M&ge;"+magStr);
				table.addColumn("Fault Rate M&ge;"+magStr);
				table.addColumn("Grid Rate M&ge;"+magStr);
			}
			table.finalizeLine();
			
			for (int d=0; d<dists.length; d++) {
				table.addColumn("% Change Within "+oDF.format(dists[d])+" km");
				
				double modelFaultMoment = distFaultPartics[d].momentRateMap.getSumZ();
				double modelGridMoment = distGridPartics[d].momentRateMap.getSumZ();
				double modelTotMoment = modelFaultMoment+modelGridMoment;
				
				double compFaultMoment = compDistFaultPartics[d].momentRateMap.getSumZ();
				double compGridMoment = compDistGridPartics[d].momentRateMap.getSumZ();
				double compTotMoment = compFaultMoment+compGridMoment;
				
				double totMomFractDiff = (modelTotMoment - compTotMoment)/compTotMoment;
				double faultMomFractDiff = (modelFaultMoment - compFaultMoment)/compFaultMoment;
				double gridMomFractDiff = (modelGridMoment - compGridMoment)/compGridMoment;
				
				table.addColumn((totMomFractDiff > 0 ? "+" : "")+pDF.format(totMomFractDiff));
				table.addColumn((faultMomFractDiff > 0 ? "+" : "")+pDF.format(faultMomFractDiff));
				table.addColumn((gridMomFractDiff > 0 ? "+" : "")+pDF.format(gridMomFractDiff));
				
				for (int m=0; m<magThresholds.length; m++) {
					double mag = magThresholds[m];
					IncrementalMagFreqDist modelFaultMFD = distFaultMFDs[d];
					IncrementalMagFreqDist modelGridMFD = distGridMFDs[d];
					IncrementalMagFreqDist modelTotMFD = sumMFD(modelFaultMFD, modelGridMFD);
					
					IncrementalMagFreqDist compFaultMFD = compDistFaultMFDs[d];
					IncrementalMagFreqDist compGridMFD = compDistGridMFDs[d];
					IncrementalMagFreqDist compTotMFD = sumMFD(compFaultMFD, compGridMFD);
					
					double totRatePDiff = cmlMFDFractDiff(modelTotMFD, compTotMFD, mag);
					double faultRatePDiff = cmlMFDFractDiff(modelFaultMFD, compFaultMFD, mag);
					double gridRatePDiff = cmlMFDFractDiff(modelGridMFD, compGridMFD, mag);
					
					table.addColumn((totRatePDiff > 0 ? "+" : "")+pDF.format(totRatePDiff));
					table.addColumn((faultRatePDiff > 0 ? "+" : "")+pDF.format(faultRatePDiff));
					table.addColumn((gridRatePDiff > 0 ? "+" : "")+pDF.format(gridRatePDiff));
				}
				table.finalizeLine();
			}
			
			lines.addAll(table.invert().build());
		}
		
		lines.add("## Regional Participation MFDs");
		lines.add(topLink); lines.add("");
		
		table = MarkdownUtils.tableBuilder();
		table.addLine("Distance Threshold", "Incremental MFDs", "Cumulative MFDs");
		
		for (int d=0; d<dists.length; d++) {
			double roundedMaxMag = Math.ceil(refMFD.getMaxX()*2d)/2d;
			Range xRange = new Range(5d, roundedMaxMag);
			Range yRange = new Range(1e-8, 1e0);
			
			List<DiscretizedFunc> incrFuncs = new ArrayList<>();
			List<DiscretizedFunc> cmlFuncs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			IncrementalMagFreqDist faultMFD = distFaultMFDs[d];
			IncrementalMagFreqDist gridMFD = distGridMFDs[d];
			IncrementalMagFreqDist fullMFD = sumMFD(faultMFD, gridMFD);
			
			fullMFD.setName("Full MFD");
			incrFuncs.add(fullMFD);
			EvenlyDiscretizedFunc fullCml = fullMFD.getCumRateDistWithOffset();
			fullCml.setName(fullMFD.getName());
			cmlFuncs.add(fullCml);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, mainColor));
			
			faultMFD.setName("Fault MFD");
			incrFuncs.add(faultMFD);
			EvenlyDiscretizedFunc faultCml = faultMFD.getCumRateDistWithOffset();
			faultCml.setName(faultMFD.getName());
			cmlFuncs.add(faultCml);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, mainColor));
			
			gridMFD.setName("Gridded MFD");
			incrFuncs.add(gridMFD);
			EvenlyDiscretizedFunc gridCml = gridMFD.getCumRateDistWithOffset();
			gridCml.setName(gridMFD.getName());
			cmlFuncs.add(gridCml);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, mainColor));
			
			if (hasComp) {
				faultMFD = compDistFaultMFDs[d];
				gridMFD = compDistGridMFDs[d];
				fullMFD = sumMFD(faultMFD, gridMFD);
				
				fullMFD.setName(compName+" Full MFD");
				incrFuncs.add(fullMFD);
				fullCml = fullMFD.getCumRateDistWithOffset();
				fullCml.setName(fullMFD.getName());
				cmlFuncs.add(fullCml);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, compColor));
				
				faultMFD.setName("Fault MFD");
				incrFuncs.add(faultMFD);
				faultCml = faultMFD.getCumRateDistWithOffset();
				faultCml.setName(faultMFD.getName());
				cmlFuncs.add(faultCml);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, compColor));
				
				gridMFD.setName("Gridded MFD");
				incrFuncs.add(gridMFD);
				gridCml = gridMFD.getCumRateDistWithOffset();
				gridCml.setName(gridMFD.getName());
				cmlFuncs.add(gridCml);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, compColor));
			}
			
			String title = locName+" MFDs within "+oDF.format(dists[d])+" km";
			
			PlotSpec incrSpec = new PlotSpec(incrFuncs, chars, title, "Magnitude", "Incremental Participation Rate (/yr)");
			incrSpec.setLegendInset(true);
			PlotSpec cmlSpec = new PlotSpec(cmlFuncs, chars, title, "Magnitude", "Cumulative Participation Rate (/yr)");
			cmlSpec.setLegendInset(true);
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			String prefix = "mfds_"+oDF.format(dists[d])+"km";
			
			gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
			gp.drawGraphPanel(incrSpec, false, true, xRange, yRange);
			PlotUtils.writePlots(resourcesDir, prefix+"_incr", gp, 1000, 850, true, true, true);
			
			gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
			gp.drawGraphPanel(cmlSpec, false, true, xRange, yRange);
			PlotUtils.writePlots(resourcesDir, prefix+"_cml", gp, 1000, 850, true, true, true);
			
			table.initNewLine();
			table.addColumn("Within "+oDF.format(dists[d])+" km");
			table.addColumn("![MFDs]("+resourcesDir.getName()+"/"+prefix+"_incr.png)");
			table.addColumn("![MFDs]("+resourcesDir.getName()+"/"+prefix+"_cml.png)");
			table.finalizeLine();
		}
		
		lines.addAll(table.build());
		lines.add("");
		
		lines.add("## Moment Rate Maps");
		lines.add(topLink); lines.add("");
		
		table = MarkdownUtils.tableBuilder();
		
		table.initNewLine();
		if (hasComp)
			table.addColumn("");
		table.addColumns("Fault Sources", "Grid Sources", "Total");
		table.finalizeLine();
		
		Color transparent = new Color(255, 255, 255, 0);
		
		CPT moCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(10d, 18d);
		moCPT.setNanColor(transparent);
		
		CPT moRateDiffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-1e16, 1e16);
		moRateDiffCPT.setNanColor(transparent);

		CPT pDiffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-100d, 100d);
		pDiffCPT.setNanColor(transparent);
		
		table.initNewLine();
		
		if (hasComp)
			table.addColumn("Model");
		
		mapMaker.plotXYZData(asLog10(solFaultPartic.momentRateMap), moCPT, "Fault Moment Rate (N-m)");
		mapMaker.plot(resourcesDir, "fault_moment_map", " ");
		table.addColumn("![Moment Map]("+resourcesDir.getName()+"/fault_moment_map.png)");
		
		mapMaker.plotXYZData(asLog10(solGridPartic.momentRateMap), moCPT, "Gridded Moment Rate (N-m)");
		mapMaker.plot(resourcesDir, "grid_moment_map", " ");
		table.addColumn("![Moment Map]("+resourcesDir.getName()+"/grid_moment_map.png)");
		
		GriddedGeoDataSet momentRateSumMap = sumMap(solFaultPartic.momentRateMap, solGridPartic.momentRateMap);
		mapMaker.plotXYZData(asLog10(momentRateSumMap), moCPT, "Total Moment Rate (N-m)");
		mapMaker.plot(resourcesDir, "total_moment_map", " ");
		table.addColumn("![Moment Map]("+resourcesDir.getName()+"/total_moment_map.png)");
		
		table.finalizeLine();
		
		if (hasComp) {
			table.initNewLine().addColumn(compName);
			
			mapMaker.plotXYZData(asLog10(compFaultPartic.momentRateMap), moCPT, compName+" Fault Moment Rate (N-m)");
			mapMaker.plot(resourcesDir, "comp_fault_moment_map", " ");
			table.addColumn("![Moment Map]("+resourcesDir.getName()+"/comp_fault_moment_map.png)");
			
			mapMaker.plotXYZData(asLog10(compGridPartic.momentRateMap), moCPT, compName+" Gridded Moment Rate (N-m)");
			mapMaker.plot(resourcesDir, "comp_grid_moment_map", " ");
			table.addColumn("![Moment Map]("+resourcesDir.getName()+"/comp_grid_moment_map.png)");
			
			GriddedGeoDataSet compMomentRateSumMap = sumMap(compFaultPartic.momentRateMap, compGridPartic.momentRateMap);
			mapMaker.plotXYZData(asLog10(compMomentRateSumMap), moCPT, compName+" Total Moment Rate (N-m)");
			mapMaker.plot(resourcesDir, "comp_total_moment_map", " ");
			table.addColumn("![Moment Map]("+resourcesDir.getName()+"/comp_total_moment_map.png)");
			
			table.finalizeLine().initNewLine().addColumn("% Change");
			
			mapMaker.plotXYZData(mapPDiff(solFaultPartic.momentRateMap, compFaultPartic.momentRateMap),
					pDiffCPT, "% Change, Fault Moment Rate");
			mapMaker.plot(resourcesDir, "comp_fault_moment_map_pDiff", " ");
			table.addColumn("![Moment Map]("+resourcesDir.getName()+"/comp_fault_moment_map_pDiff.png)");
			
			mapMaker.plotXYZData(mapPDiff(solGridPartic.momentRateMap, compGridPartic.momentRateMap),
					pDiffCPT, "% Change, Gridded Moment Rate");
			mapMaker.plot(resourcesDir, "comp_grid_moment_map_pDiff", " ");
			table.addColumn("![Moment Map]("+resourcesDir.getName()+"/comp_grid_moment_map_pDiff.png)");
			
			mapMaker.plotXYZData(mapPDiff(momentRateSumMap, compMomentRateSumMap),
					pDiffCPT, "% Change, Total Moment Rate");
			mapMaker.plot(resourcesDir, "comp_total_moment_map_pDiff", " ");
			table.addColumn("![Moment Map]("+resourcesDir.getName()+"/comp_total_moment_map_pDiff.png)");
			
			table.finalizeLine().initNewLine().addColumn("Difference");
			
			mapMaker.plotXYZData(mapDiff(solFaultPartic.momentRateMap, compFaultPartic.momentRateMap),
					moRateDiffCPT, "Difference, Fault Moment Rate (N-m)");
			mapMaker.plot(resourcesDir, "comp_fault_moment_map_diff", " ");
			table.addColumn("![Moment Map]("+resourcesDir.getName()+"/comp_fault_moment_map_diff.png)");
			
			mapMaker.plotXYZData(mapDiff(solGridPartic.momentRateMap, compGridPartic.momentRateMap),
					moRateDiffCPT, "Difference, Gridded Moment Rate (N-m)");
			mapMaker.plot(resourcesDir, "comp_grid_moment_map_diff", " ");
			table.addColumn("![Moment Map]("+resourcesDir.getName()+"/comp_grid_moment_map_diff.png)");
			
			mapMaker.plotXYZData(mapDiff(momentRateSumMap, compMomentRateSumMap),
					moRateDiffCPT, "Difference, Total Moment Rate (N-m)");
			mapMaker.plot(resourcesDir, "comp_total_moment_map_diff", " ");
			table.addColumn("![Moment Map]("+resourcesDir.getName()+"/comp_total_moment_map_diff.png)");
			
			table.finalizeLine();
		}
		
		lines.addAll(table.build());
		lines.add("");
		
		lines.add("## Rate Maps");
		lines.add(topLink); lines.add("");
		
		table = MarkdownUtils.tableBuilder();
		table.initNewLine();
		if (hasComp)
			table.addColumn("");
		for (double mag : magThresholds)
			table.addColumn("M&ge;"+oDF.format(mag));
		table.finalizeLine();
		
		GriddedGeoDataSet[] modelSums = new GriddedGeoDataSet[magThresholds.length];
		CPT[] rateCPTs = new CPT[magThresholds.length];
		CPT[] rateDiffCPTs = new CPT[magThresholds.length];
		for (int m=0; m<magThresholds.length; m++) {
			modelSums[m] = sumMap(solFaultPartic.particRateMaps[m], solGridPartic.particRateMaps[m]);
			double maxRate = modelSums[m].getMaxZ();
			
			double logMax = Math.ceil(Math.log10(maxRate));
			
			rateCPTs[m] = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(logMax-6d, logMax);
			rateCPTs[m].setNanColor(transparent);
			
			double diffMax = Math.pow(10, logMax-1);
			
			rateDiffCPTs[m] = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-diffMax, diffMax);
			rateDiffCPTs[m].setNanColor(transparent);
		}
		
		table.initNewLine();
		if (hasComp)
			table.addColumn("Model");
		for (int m=0; m<magThresholds.length; m++) {
			mapMaker.plotXYZData(asLog10(modelSums[m]), rateCPTs[m], "M≥"+oDF.format(magThresholds[m])+" Participation Rate");
			String prefix = "rate_map_m"+oDF.format(magThresholds[m]);
			mapMaker.plot(resourcesDir, prefix, " ");
			table.addColumn("![Rate Map]("+resourcesDir.getName()+"/"+prefix+".png)");
		}
		table.finalizeLine();
		if (hasComp) {
			table.initNewLine().addColumn(compName);
			GriddedGeoDataSet[] compSums = new GriddedGeoDataSet[magThresholds.length];
			for (int m=0; m<magThresholds.length; m++) {
				compSums[m] = sumMap(compFaultPartic.particRateMaps[m], compGridPartic.particRateMaps[m]);
				mapMaker.plotXYZData(asLog10(compSums[m]), rateCPTs[m], compName+" "+"M≥"+oDF.format(magThresholds[m])+" Participation Rate");
				String prefix = "comp_rate_map_m"+oDF.format(magThresholds[m]);
				mapMaker.plot(resourcesDir, prefix, " ");
				table.addColumn("![Rate Map]("+resourcesDir.getName()+"/"+prefix+".png)");
			}
			table.finalizeLine().initNewLine().addColumn("% Change");
			for (int m=0; m<magThresholds.length; m++) {
				GriddedGeoDataSet pDiff = mapPDiff(modelSums[m], compSums[m]);
				mapMaker.plotXYZData(pDiff, pDiffCPT, "% Change in M≥"+oDF.format(magThresholds[m])+" Participation Rate");
				String prefix = "rate_map_pDiff_m"+oDF.format(magThresholds[m]);
				mapMaker.plot(resourcesDir, prefix, " ");
				table.addColumn("![Change Map]("+resourcesDir.getName()+"/"+prefix+".png)");
			}
			table.finalizeLine().initNewLine().addColumn("Difference");
			for (int m=0; m<magThresholds.length; m++) {
				GriddedGeoDataSet diff = mapDiff(modelSums[m], compSums[m]);
				mapMaker.plotXYZData(diff, rateDiffCPTs[m], "Difference in M≥"+oDF.format(magThresholds[m])+" Participation Rate");
				String prefix = "rate_map_diff_m"+oDF.format(magThresholds[m]);
				mapMaker.plot(resourcesDir, prefix, " ");
				table.addColumn("![Change Map]("+resourcesDir.getName()+"/"+prefix+".png)");
			}
			table.finalizeLine();
		}
		lines.addAll(table.build());
		lines.add("");
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}
	
	private static Options createOptions() {
		Options ops = new Options();
		
		ops.addOption("ns", "nehrp-site", true, "Name of a NEHRP site to use. Must supply this or a location with "
				+ "--location <lat,lon>. Will be matched to list of NEHRP sites, ignoring capitalization and spacing.");
		ops.addOption("l", "location", true, "Site location specified as lat,lon. Must supply this or --nehrp-site <name>.");
		ops.addOption("n", "name", true, "Location name to use in plots.");
		ops.addRequiredOption("od", "output-dir", true, "Path to output directory");
		ops.addRequiredOption("if", "input-file", true, "Path to input fault system solution");
		ops.addOption("cne", "comp-nshm-erf", true, "Path to NSHM ERF to use as a comparison");
		ops.addOption("cs", "comp-sol", true, "Path to comparison fault system solution");
		ops.addOption("cn", "comp-name", true, "Name of NSHM ERF (default: 'Comparison')");
		ops.addOption("gm", "gmpe", true, "Sets GMPE. Note that this will be overriden if the Logic Tree "
				+ "supplies GMPE choices. Default: "+GMPE_DEFAULT.name());
		ops.addOption("p", "periods", true, "Calculation period(s). Mutliple can be comma separated");
		
		return ops;
	}
	
	private static DecimalFormat oDF = new DecimalFormat("0.##");
	private static DecimalFormat pDF = new DecimalFormat("0.00%");
	
	public static class RegionalParticipationResult {
		public final double[] magThresholds;
		public final GriddedGeoDataSet[] particRateMaps;
		public final GriddedGeoDataSet[] nuclRateMaps;
		public final GriddedGeoDataSet momentRateMap;
		public final IncrementalMagFreqDist particMFD;
		
		public RegionalParticipationResult(double[] magThresholds, GriddedGeoDataSet[] particRateMaps,
				GriddedGeoDataSet[] nuclRateMaps,
				GriddedGeoDataSet momentRateMap, IncrementalMagFreqDist particMFD) {
			this.magThresholds = magThresholds;
			this.particRateMaps = particRateMaps;
			this.nuclRateMaps = nuclRateMaps;
			this.momentRateMap = momentRateMap;
			this.particMFD = particMFD;
		}
	}
	
	public static RegionalParticipationResult calcNshmErfPartic(NshmErf erf, GriddedRegion gridReg,
			double[] magThresholds, EvenlyDiscretizedFunc refMFD) {
		ExecutorService exec = Executors.newFixedThreadPool(FaultSysTools.defaultNumThreads());
		
		List<Future<?>> futures = new ArrayList<>();
		
		// create a larger region around the gridded region, we'll skip any sources whose centroid is outside of this
		Location botLeft = new Location(gridReg.getMinGridLat(), gridReg.getMinGridLon());
		Location topRight = new Location(gridReg.getMaxGridLat(), gridReg.getMaxGridLon());
		botLeft = LocationUtils.location(botLeft, 5d*Math.PI/4d, 50d);
		topRight = LocationUtils.location(topRight, Math.PI/4d, 50d);
		Region largerTestReg = new Region(botLeft, topRight);
		
		Map<RuptureSurface, int[]> mappedIndexes = new ConcurrentHashMap<>();

		GriddedGeoDataSet[] particRateMaps = new GriddedGeoDataSet[magThresholds.length];
		for (int m=0; m<magThresholds.length; m++)
			particRateMaps[m] = new GriddedGeoDataSet(gridReg);
		GriddedGeoDataSet[] nuclRateMaps = new GriddedGeoDataSet[magThresholds.length];
		for (int m=0; m<magThresholds.length; m++)
			nuclRateMaps[m] = new GriddedGeoDataSet(gridReg);
		GriddedGeoDataSet momentRateMap = new GriddedGeoDataSet(gridReg);
		IncrementalMagFreqDist particMFD = refMFD == null ? null :
			new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		
		System.out.println("Building NshmErf participation rate maps");
		for (ProbEqkSource source : erf) {
			futures.add(exec.submit(new Runnable() {
				
				@Override
				public void run() {
					// see if we can skip
					boolean skip = true;
					for (ProbEqkRupture rup : source) {
						if (!canSkipNshmSurf(largerTestReg, rup.getRuptureSurface())) {
							skip = false;
							break;
						}
					}
					if (skip)
						return;
					double[][] sourceParticRates = new double[gridReg.getNodeCount()][magThresholds.length];
					double[][] sourceNuclRates = new double[gridReg.getNodeCount()][magThresholds.length];
					double[] sourceMoRates = new double[gridReg.getNodeCount()];
					double[] mfdRates = refMFD == null ? null : new double[refMFD.size()];
					for (ProbEqkRupture rup : source) {
						double mag = rup.getMag();
						double rate = rup.getMeanAnnualRate(1d);
						RuptureSurface surf = rup.getRuptureSurface();
						double moment = rate*MagUtils.magToMoment(mag);
						int fullIndexCount = 0;
						List<int[]> fullIndexes = new ArrayList<>();
						if (surf instanceof CompoundSurface) {
							List<? extends RuptureSurface> subSurfs = ((CompoundSurface)surf).getSurfaceList();
							for (RuptureSurface subSurf : subSurfs) {
								int[] indexes;
									indexes = mappedIndexes.get(subSurf);
								if (indexes == null) {
									LocationList locs = subSurf.getEvenlyDiscritizedListOfLocsOnSurface();
									indexes = new int[locs.size()];
									for (int i=0; i<indexes.length; i++)
										indexes[i] = gridReg.indexForLocation(locs.get(i));
										mappedIndexes.put(subSurf, indexes);
								}
								fullIndexCount += indexes.length;
								fullIndexes.add(indexes);
							}
						} else {
							LocationList locs = surf.getEvenlyDiscritizedListOfLocsOnSurface();
							if (locs.size() == 0)
								System.err.println("WARNING: source haze no points? rate="+rate);
							int[] indexes = new int[locs.size()];
							for (int i=0; i<indexes.length; i++)
								indexes[i] = gridReg.indexForLocation(locs.get(i));
							fullIndexes.add(indexes);
							fullIndexCount += indexes.length;
						}
						double momentEach = moment/(double)fullIndexCount;
						double rateEach = rate/(double)fullIndexCount;
						BitSet cellParticipation = new BitSet(sourceMoRates.length);
						boolean any = false;
						for (int[] indexes : fullIndexes) {
							for (int idx : indexes) {
								if (idx >= 0d) {
									cellParticipation.set(idx);
									sourceMoRates[idx] += momentEach;
									for (int m=0; m<magThresholds.length; m++)
										if (mag >= magThresholds[m])
											sourceNuclRates[idx][m] += rateEach;
									any = true;
								}
							}
						}
						for (int i = cellParticipation.nextSetBit(0); i >= 0; i = cellParticipation.nextSetBit(i+1))
							for (int m=0; m<magThresholds.length; m++)
								if (mag >= magThresholds[m])
									sourceParticRates[i][m] += rate;
						if (any && mfdRates != null)
							mfdRates[refMFD.getClosestXIndex(mag)] += rate;
					}
					synchronized (particRateMaps) {
						for (int i=0; i<sourceParticRates.length; i++) {
							if (sourceMoRates[i] > 0) {
								for (int m=0; m<magThresholds.length; m++) {
									if (sourceParticRates[i][m] > 0)
										particRateMaps[m].set(i, particRateMaps[m].get(i)+sourceParticRates[i][m]);
									if (sourceNuclRates[i][m] > 0)
										nuclRateMaps[m].set(i, particRateMaps[m].get(i)+sourceNuclRates[i][m]);
								}
								momentRateMap.set(i, momentRateMap.get(i)+sourceMoRates[i]);
							}
						}
						if (mfdRates != null)
							for (int i=0; i<mfdRates.length; i++)
								particMFD.add(i, mfdRates[i]);
					}
				}
			}));
		}
		
		System.out.println("Waiting on "+futures.size()+" NshmErf source futures");
		for (Future<?> future : futures) {
			try {
				future.get();
			} catch (InterruptedException | ExecutionException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
		}
		
		exec.shutdown();
		
		return new RegionalParticipationResult(magThresholds, particRateMaps, nuclRateMaps, momentRateMap, particMFD);
	}
	
	static boolean canSkipNshmSurf(Region reg, RuptureSurface surf) {
		if (surf instanceof NshmSurface) {
			NshmSurface nshmSurf = (NshmSurface)surf;
			Location centroid;
			try {
				centroid = nshmSurf.centroid();
			} catch (Exception e) {
				return false;
			}
			return !reg.contains(centroid);
		} else if (surf instanceof CompoundSurface) {
			for (RuptureSurface subSurf : ((CompoundSurface)surf).getSurfaceList()) {
				if (!canSkipNshmSurf(reg, subSurf))
					// at least one is false, return false
					return false;
			}
			// all of them are true, so return true
			return true;
		} else {
			return false;
		}
	}
	
	public static RegionalParticipationResult calcFSSFaultPartic(FaultSystemSolution sol, GriddedRegion gridReg,
			double[] magThresholds, EvenlyDiscretizedFunc refMFD) {
		FaultSystemRupSet rupSet = sol.getRupSet();
		FaultGridAssociations assoc = FaultGridAssociations.getIntersectionAssociations(rupSet, gridReg);
		List<IncrementalMagFreqDist> sectNuclMFDs = NucleationRatePlot.calcNuclMFDs(sol);
		GriddedGeoDataSet[] faultParticRates = new GriddedGeoDataSet[magThresholds.length];
		for (int m=0; m<magThresholds.length; m++)
			faultParticRates[m] = new GriddedGeoDataSet(gridReg, false);
		GriddedGeoDataSet[] faultNuclRates = new GriddedGeoDataSet[magThresholds.length];
		for (int m=0; m<magThresholds.length; m++)
//			faultNuclRates[m] = new GriddedGeoDataSet(gridReg, false);
			faultNuclRates[m] = NucleationRatePlot.calcFaultNucleationRates(gridReg, sol, assoc, sectNuclMFDs, magThresholds[m]);
		IncrementalMagFreqDist mfd = refMFD == null ? null :
			new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		// now do fault participation
		BitSet sectBits = new BitSet(rupSet.getNumSections());
		for (int r=0; r<rupSet.getNumRuptures(); r++) {
			sectBits.clear();
			for (int sectIndex : rupSet.getSectionsIndicesForRup(r))
				for (int gridIndex : assoc.getNodeFractions(sectIndex).keySet())
					sectBits.set(gridIndex);
			double mag = rupSet.getMagForRup(r);
			double rate = sol.getRateForRup(r);
			boolean any = false;
			for (int i = sectBits.nextSetBit(0); i >= 0; i = sectBits.nextSetBit(i+1)) {
				any = true;
				for (int m=0; m<magThresholds.length; m++) {
					if (mag >= magThresholds[m])
						faultParticRates[m].set(i, faultParticRates[m].get(i)+rate);
				}
			}
			if (any && mfd != null)
				mfd.add(mfd.getClosestXIndex(mag), rate);
		}
		GriddedGeoDataSet faultMoRates = NucleationRatePlot.calcFaultNucleationMomentRates(
				gridReg, sol, assoc, sectNuclMFDs);
		
		return new RegionalParticipationResult(magThresholds, faultParticRates, faultNuclRates, faultMoRates, mfd);
	}
	
	public static RegionalParticipationResult calcFSSGriddedPartic(FaultSystemSolution sol, GriddedRegion gridReg,
			double[] magThresholds, EvenlyDiscretizedFunc refMFD) {
		GridSourceProvider gridProv = sol.getGridSourceProvider();
		if (gridProv == null) {
			 // return zeros
			GriddedGeoDataSet[] griddedRates = new GriddedGeoDataSet[magThresholds.length];
			for (int m=0; m<magThresholds.length; m++)
				new GriddedGeoDataSet(gridReg);
			GriddedGeoDataSet griddedMoRates = new GriddedGeoDataSet(gridReg);
			IncrementalMagFreqDist mfd = null;
			if (refMFD != null)
				mfd = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
			return new RegionalParticipationResult(magThresholds, griddedRates, griddedRates, griddedMoRates, mfd);
		}
		GriddedGeoDataSet[] griddedRates = new GriddedGeoDataSet[magThresholds.length];
		for (int m=0; m<magThresholds.length; m++)
			// nucleation for gridded
			griddedRates[m] = NucleationRatePlot.calcGriddedNucleationRates(
					sol.getGridSourceProvider(), gridReg, magThresholds[m]);
		GriddedGeoDataSet griddedMoRates = NucleationRatePlot.calcGriddedNucleationMomentRates(gridProv, gridReg);
		IncrementalMagFreqDist mfd = null;
		if (refMFD != null) {
			mfd = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
			for (int i=0; i<gridProv.size(); i++) {
				Location gridLoc = gridProv.getGriddedRegion().getLocation(i);
				int index = gridReg.indexForLocation(gridLoc);
				if (index >= 0) {
					IncrementalMagFreqDist nodeMFD = gridProv.getMFD(i, mfd.getMinX());
					for (int j=0; j<nodeMFD.size(); j++) {
						double mag = nodeMFD.getX(j);
						double rate = nodeMFD.getY(j);
						if (rate > 0d)
							mfd.add(mfd.getClosestXIndex(mag), rate);
					}
				}
			}
		}
		
		return new RegionalParticipationResult(magThresholds, griddedRates, griddedRates, griddedMoRates, mfd);
	}
	
	public static RegionalParticipationResult calcFSSPartic(FaultSystemSolution sol, GriddedRegion gridReg,
			double[] magThresholds, EvenlyDiscretizedFunc refMFD) {
		RegionalParticipationResult faultPartic = calcFSSFaultPartic(sol, gridReg, magThresholds, refMFD);
		RegionalParticipationResult gridPartic = calcFSSGriddedPartic(sol, gridReg, magThresholds, refMFD);
		
		GriddedGeoDataSet[] sumParticRates = new GriddedGeoDataSet[magThresholds.length];
		for (int m=0; m<magThresholds.length; m++)
			sumParticRates[m] = new GriddedGeoDataSet(gridReg);
		GriddedGeoDataSet[] sumNuclRates = new GriddedGeoDataSet[magThresholds.length];
		for (int m=0; m<magThresholds.length; m++)
			sumNuclRates[m] = new GriddedGeoDataSet(gridReg);
		GriddedGeoDataSet sumMomentRate = new GriddedGeoDataSet(gridReg, false);
		
		for (int i=0; i<sumMomentRate.size(); i++) {
			sumMomentRate.set(i, faultPartic.momentRateMap.get(i)+gridPartic.momentRateMap.get(i));
			for (int m=0; m<magThresholds.length; m++) {
				sumParticRates[m].set(i, faultPartic.particRateMaps[m].get(i)+gridPartic.particRateMaps[m].get(i));
				sumNuclRates[m].set(i, faultPartic.nuclRateMaps[m].get(i)+gridPartic.nuclRateMaps[m].get(i));
			}
		}
		
		IncrementalMagFreqDist mfd = null;
		if (refMFD != null) {
			mfd = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
			for (int i=0; i<mfd.size(); i++)
				mfd.set(i, faultPartic.particMFD.getY(i)+gridPartic.particMFD.getY(i));
		}
		
		return new RegionalParticipationResult(magThresholds, sumParticRates, sumNuclRates, sumMomentRate, mfd);
	}
	
	private static DiscretizedFunc toLinearCurve(DiscretizedFunc curve, DiscretizedFunc xVals) {
		DiscretizedFunc ret = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<xVals.size(); i++)
			ret.set(xVals.getX(i), curve.getY(i));
		return ret;
	}
	
	private static DiscretizedFunc sumCurve(DiscretizedFunc curve1, DiscretizedFunc curve2) {
		Preconditions.checkState(curve1.size() == curve2.size());
		DiscretizedFunc ret = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<curve1.size(); i++)
			ret.set(curve1.getX(i), 1d - (1d-curve1.getY(i))*(1d-curve2.getY(i)));
		return ret;
	}
	
	private static IncrementalMagFreqDist sumMFD(IncrementalMagFreqDist mfd1, IncrementalMagFreqDist mfd2) {
		Preconditions.checkState(mfd1.size() == mfd2.size());
		IncrementalMagFreqDist ret = new IncrementalMagFreqDist(mfd1.getMinX(), mfd1.size(), mfd1.getDelta());
		for (int i=0; i<mfd1.size(); i++)
			ret.set(mfd1.getX(i), mfd1.getY(i)+mfd2.getY(i));
		return ret;
	}
	
	private static GriddedGeoDataSet sumMap(GriddedGeoDataSet map1, GriddedGeoDataSet map2) {
		Preconditions.checkState(map1.size() == map2.size());
		GriddedGeoDataSet ret = new GriddedGeoDataSet(map1.getRegion());
		for (int i=0; i<ret.size(); i++)
			ret.set(i, map1.get(i)+map2.get(i));
		return ret;
	}
	
	private static GriddedGeoDataSet mapPDiff(GriddedGeoDataSet map1, GriddedGeoDataSet map2) {
		Preconditions.checkState(map1.size() == map2.size());
		GriddedGeoDataSet ret = new GriddedGeoDataSet(map1.getRegion());
		for (int i=0; i<ret.size(); i++)
			ret.set(i, 100d*(map1.get(i)-map2.get(i))/map2.get(i));
		return ret;
	}
	
	private static GriddedGeoDataSet mapDiff(GriddedGeoDataSet map1, GriddedGeoDataSet map2) {
		Preconditions.checkState(map1.size() == map2.size());
		GriddedGeoDataSet ret = new GriddedGeoDataSet(map1.getRegion());
		for (int i=0; i<ret.size(); i++)
			ret.set(i, map1.get(i)-map2.get(i));
		return ret;
	}
	
	private static GriddedGeoDataSet asLog10(GriddedGeoDataSet map) {
		map = map.copy();
		for (int i=0; i<map.size(); i++)
			if (map.get(i) == 0d)
				map.set(i, Double.NaN);
		map.log10();
		return map;
	}
	
	private static double calcCurveIML(ReturnPeriods rp, DiscretizedFunc curve) {
		double curveLevel = rp.oneYearProb;
		if (curveLevel > curve.getMaxY())
			return 0d;
		else if (curveLevel < curve.getMinY())
			// saturated
			return curve.getMaxX();
		else
			return curve.getFirstInterpolatedX_inLogXLogYDomain(curveLevel);
	}
	
	private static double cmlMFDFractDiff(IncrementalMagFreqDist mfd1, IncrementalMagFreqDist mfd2, double mag) {
		EvenlyDiscretizedFunc cml1 = mfd1.getCumRateDistWithOffset();
		EvenlyDiscretizedFunc cml2 = mfd2.getCumRateDistWithOffset();
		
		int index = cml1.getClosestXIndex(mag);
		return (cml1.getY(index)-cml2.getY(index))/cml2.getY(index);
	}

}
