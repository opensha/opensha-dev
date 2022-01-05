package scratch.kevin.ucerf3;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.net.InetAddress;
import java.net.UnknownHostException;
import java.nio.charset.Charset;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.time.Duration;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.time.temporal.ChronoUnit;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Date;
import java.util.Enumeration;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.TimeZone;
import java.util.concurrent.ExecutionException;
import java.util.regex.Matcher;
import java.util.stream.Collectors;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import javax.swing.JFrame;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.dom4j.DocumentException;
import org.jfree.data.Range;
import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.region.CaliforniaRegions.RELM_SOCAL;
import org.opensha.commons.data.siteData.OrderedSiteDataProviderList;
import org.opensha.commons.data.siteData.SiteDataValue;
import org.opensha.commons.data.siteData.impl.WillsMap2006;
import org.opensha.commons.data.siteData.impl.WillsMap2015;
import org.opensha.commons.data.uncertainty.Uncertainty;
import org.opensha.commons.data.uncertainty.UncertaintyBoundType;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.PlaneUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.editor.AbstractParameterEditor;
import org.opensha.commons.param.editor.impl.NumericTextField;
import org.opensha.commons.param.impl.DoubleParameter;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FaultUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.FileUtils;
import org.opensha.commons.util.IDPairing;
import org.opensha.commons.util.ServerPrefUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.commons.util.modules.OpenSHA_Module;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.hazardMap.HazardCurveSetCalculator;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.calc.ERF_Calculator;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfigurationFactory;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.InversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.UncertainDataConstraint;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupMFDsModule;
import org.opensha.sha.earthquake.faultSysSolution.modules.WaterLevelRates;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.SolMFDPlot;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.FaultSubsectionCluster;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Jump;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.GeoJSONFaultReader;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SectionDistanceAzimuthCalculator;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.UniqueRupture;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.earthquake.observedEarthquake.parsers.UCERF3_CatalogParser;
import org.opensha.sha.earthquake.param.BackgroundRupParam;
import org.opensha.sha.earthquake.param.BackgroundRupType;
import org.opensha.sha.earthquake.param.HistoricOpenIntervalParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.earthquake.rupForecastImpl.PointSource13b;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.StirlingGriddedSurface;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.AttenuationRelationship;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.USGS_Combined_2004_AttenRel;
import org.opensha.sha.imr.attenRelImpl.ngaw2.ASK_2014;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.SiteParams.DepthTo1pt0kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.DepthTo2pt5kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimUtils;
import org.opensha.sha.util.SiteTranslator;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Table;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.UCERF3.analysis.FaultSysSolutionERF_Calc;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.SlipAlongRuptureModels;
import scratch.UCERF3.enumTreeBranches.TotalMag5Rate;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO.BinarayCatalogsIterable;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO.ETAS_Catalog;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.ETAS_MultiSimAnalysisTools;
import scratch.UCERF3.erf.ETAS.ETAS_SimAnalysisTools;
import scratch.UCERF3.erf.ETAS.ETAS_Simulator.TestScenario;
import scratch.UCERF3.erf.ETAS.ETAS_Utils;
import scratch.UCERF3.erf.ETAS.FaultSystemSolutionERF_ETAS;
import scratch.UCERF3.erf.ETAS.association.FiniteFaultMappingData;
import scratch.UCERF3.erf.ETAS.launcher.ETAS_Config;
import scratch.UCERF3.erf.ETAS.launcher.ETAS_Launcher;
import scratch.UCERF3.erf.mean.MeanUCERF3;
import scratch.UCERF3.erf.mean.MeanUCERF3.Presets;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;
import scratch.UCERF3.griddedSeismicity.AbstractGridSourceProvider;
import scratch.UCERF3.inversion.CommandLineInversionRunner;
import scratch.UCERF3.inversion.InversionFaultSystemRupSetFactory;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.logicTree.APrioriBranchWeightProvider;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;
import scratch.UCERF3.utils.DeformationModelFetcher;
import scratch.UCERF3.utils.FaultSectionDataWriter;
import scratch.UCERF3.utils.U3FaultSystemIO;
import scratch.UCERF3.utils.LastEventData;
import scratch.UCERF3.utils.UCERF3_DataUtils;
import scratch.UCERF3.utils.aveSlip.U3AveSlipConstraint;
import scratch.UCERF3.utils.paleoRateConstraints.U3PaleoRateConstraint;
import scratch.UCERF3.utils.paleoRateConstraints.UCERF3_PaleoRateConstraintFetcher;
import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.ruptures.RSQSimBBP_Config;
import smile.stat.distribution.GaussianMixture;

public class PureScratch {

	private static void test1() {
		Location hypoLoc = new Location(34, -118);
		int num = 1000;
		double mag = 5d;

		ETAS_Utils utils = new ETAS_Utils();

		DefaultXY_DataSet xy = new DefaultXY_DataSet();

		int retries = 0;

		for (int i=0; i<num; i++) {
			double radius = ETAS_Utils.getRuptureRadiusFromMag(mag);
			double testRadius = radius+1;
			Location loc = null;
			while(testRadius>radius) {
				double lat = hypoLoc.getLatitude()+(2.0*utils.getRandomDouble()-1.0)*(radius/111.0);
				double lon = hypoLoc.getLongitude()+(2.0*utils.getRandomDouble()-1.0)*(radius/(111*Math.cos(hypoLoc.getLatRad())));
				double depthBottom = hypoLoc.getDepth()+radius;
				if(depthBottom>24.0)
					depthBottom=24.0;
				double depthTop = hypoLoc.getDepth()-radius;
				if(depthTop<0.0)
					depthTop=0.0;
				double depth = depthTop + utils.getRandomDouble()*(depthBottom-depthTop);
				loc = new Location(lat,lon,depth);
				testRadius=LocationUtils.linearDistanceFast(loc, hypoLoc);
				retries++;
			}
			xy.set(loc.getLongitude(), loc.getLatitude());
		}

		System.out.println("Retries: "+retries);
		GraphWindow gw = new GraphWindow(xy, "Circular Random Test", new PlotCurveCharacterstics(PlotSymbol.X, 2f, Color.BLACK));
		gw.setAxisRange(xy.getMinX(), xy.getMaxX(), xy.getMinY(), xy.getMaxY());
		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
	}

	private static void test2() throws IOException, DocumentException {
		FaultSystemSolution sol = U3FaultSystemIO.loadSol(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		FaultSystemRupSet rupSet = sol.getRupSet();

		HashSet<Integer> parents = new HashSet(FaultModels.FM3_1.getNamedFaultsMapAlt().get("San Andreas"));

		Region soCal = new CaliforniaRegions.RELM_SOCAL();

		Map<Integer, Boolean> safSectsInSoCal = Maps.newHashMap();
		for (int sectIndex=0; sectIndex<rupSet.getNumSections(); sectIndex++) {
			FaultSection sect = rupSet.getFaultSectionData(sectIndex);
			if (!parents.contains(sect.getParentSectionId()))
				continue;
			boolean inside = false;
			for (Location loc : sect.getFaultTrace()) {
				if (soCal.contains(loc)) {
					inside = true;
					break;
				}
			}
			safSectsInSoCal.put(sectIndex, inside);
			//			System.out.println(sect.getName()+": "+inside);
		}

		int numSAF = 0;
		int numPartiallySAFSoCal = 0;
		int numOnlySAFSoCal = 0;
		rupLoop:
			for (int rupIndex = 0; rupIndex<rupSet.getNumRuptures(); rupIndex++) {
				for (int parent : rupSet.getParentSectionsForRup(rupIndex)) {
					if (!parents.contains(parent))
						continue rupLoop;
				}
				numSAF++;
				boolean partial = false;
				boolean only = true;
				for (int sectIndex : rupSet.getSectionsIndicesForRup(rupIndex)) {
					boolean inside = safSectsInSoCal.get(sectIndex);
					partial = partial || inside;
					only = only && inside;
				}
				if (partial)
					numPartiallySAFSoCal++;
				if (only)
					numOnlySAFSoCal++;
			}

		System.out.println(numSAF+"/"+rupSet.getNumRuptures()+" ruputures are only on SAF");
		System.out.println(numPartiallySAFSoCal+" of those are at least partially in SoCal");
		System.out.println(numOnlySAFSoCal+" of those are entirely in SoCal");
	}

	private static void test3() throws IOException, DocumentException {
		FaultSystemSolution sol = U3FaultSystemIO.loadSol(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		FaultSystemRupSet rupSet = sol.getRupSet();

		Region region = new CaliforniaRegions.SF_BOX();

		Map<String, Integer> parentsInBox = Maps.newHashMap();

		for (int sectIndex=0; sectIndex<rupSet.getNumSections(); sectIndex++) {
			FaultSection sect = rupSet.getFaultSectionData(sectIndex);
			boolean inside = false;
			for (Location loc : sect.getFaultTrace()) {
				if (region.contains(loc)) {
					inside = true;
					break;
				}
			}
			if (inside)
				parentsInBox.put(sect.getParentSectionName(), sect.getParentSectionId());
			//			System.out.println(sect.getName()+": "+inside);
		}

		List<String> names = Lists.newArrayList(parentsInBox.keySet());
		Collections.sort(names);

		for (String name : names)
			System.out.println(parentsInBox.get(name)+". "+name);
	}

	private static void test4() throws IOException, DocumentException {
		FaultSystemSolution sol_31 = U3FaultSystemIO.loadSol(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		FaultSystemSolution sol_32 = U3FaultSystemIO.loadSol(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_2_MEAN_BRANCH_AVG_SOL.zip"));

		System.out.println("FM3.1: "+sol_31.getRupSet().getNumRuptures());
		System.out.println("FM3.2: "+sol_32.getRupSet().getNumRuptures());
	}

	private static void test5() throws IOException {
		List<? extends FaultSection> subSects = new DeformationModelFetcher(
				FaultModels.FM3_1, DeformationModels.GEOLOGIC,
				UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, 0.1d).getSubSectionList();
		LastEventData.populateSubSects(subSects, LastEventData.load());

		List<LocationList> fmTraces = Lists.newArrayList();
		for (int i=0; i<subSects.size(); i++) {
			fmTraces.add(subSects.get(i).getFaultTrace());
			if (subSects.get(i).getName().contains("Baker"))
				System.out.println(i+". "+subSects.get(i).getName());
			//			fmSubSectIndexMap.put(fm, subSects.get(i).getName(), i);
		}
	}

	private static void test6() throws ZipException, IOException, DocumentException {
		FaultSystemRupSet rupSet = U3FaultSystemIO.loadRupSet(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));

		for (FaultSection sect : rupSet.getFaultSectionDataList()) {
			System.out.println(sect.getSectionId()+". "+sect.getName());
			if (sect.getSectionId() > 50)
				break;
		}

		Map<String, List<Integer>> parentCounts = Maps.newHashMap();

		for (FaultSection sect : rupSet.getFaultSectionDataList()) {
			String parentName = sect.getParentSectionName();
			List<Integer> indexes = parentCounts.get(parentName);
			if (indexes == null) {
				indexes = Lists.newArrayList();
				parentCounts.put(parentName, indexes);
			}
			indexes.add(sect.getSectionId());
		}

		List<String> over10s = Lists.newArrayList();
		for (String parentName : parentCounts.keySet()) {
			if (parentCounts.get(parentName).size() > 10)
				over10s.add(parentName);
		}
		Collections.sort(over10s);
		System.out.println(over10s.size()+"/"+parentCounts.size()+" sections are affected:");
		for (String name : over10s)
			System.out.println("\t"+name);
	}

	private static void test7() throws IOException {
		List<? extends List<ETAS_EqkRupture>> catalogs = ETAS_CatalogIO.loadCatalogsBinary(
				new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
						//				+ "2016_01_05-spontaneous-10000yr-mc10-applyGrGridded-full_td-noApplyLTR/results_first50_m4.bin"));
						+ "2016_01_05-spontaneous-10000yr-mc10-applyGrGridded-full_td-noApplyLTR/results_second50_m4.bin"));
		MinMaxAveTracker track = new MinMaxAveTracker();
		for (List<ETAS_EqkRupture> catalog : catalogs)
			track.addValue(ETAS_MultiSimAnalysisTools.calcDurationYears(catalog));
		System.out.println(track);
	}

	private static void test8() throws Exception {
		ObsEqkRupList loadedRups = UCERF3_CatalogParser.loadCatalog(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/EarthquakeCatalog/"
						+ "ofr2013-1165_EarthquakeCat.txt"));
		File xmlFile = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/EarthquakeCatalog/"
				+ "finite_fault_mappings.xml");
		FaultSystemRupSet rupSet = U3FaultSystemIO.loadRupSet(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		FiniteFaultMappingData.loadRuptureSurfaces(xmlFile, loadedRups, FaultModels.FM3_1, rupSet);
		for (ObsEqkRupture rup : loadedRups) {
			if (rup.getEventId().equals("14607652")) {
				System.out.println("Writing!");
				FileWriter fw = new FileWriter("/tmp/el_mayor.txt");
				for (Location loc : rup.getRuptureSurface().getEvenlyDiscritizedListOfLocsOnSurface())
					fw.write(loc.getLatitude()+" "+loc.getLongitude()+" "+loc.getDepth()+"\n");
				fw.close();
				return;
			}
		}
	}

	private static void test9() {
		SimpleDateFormat catDateFormat = new SimpleDateFormat("yyyy\tMM\tdd\tHH\tmm\tss.SSS");
		TimeZone utc = TimeZone.getTimeZone("UTC");
		catDateFormat.setTimeZone(utc);
		Date date = new Date();
		System.out.println("Date object: "+date);

		System.out.println("UTC Formatted:");
		System.out.println("\tyyyy\tMM\tdd\tHH\tmm\tss.SSS");
		System.out.println("\t"+catDateFormat.format(date));

		int startYear = 2012;
		long time = Math.round((startYear-1970.0)*ProbabilityModelsCalc.MILLISEC_PER_YEAR)+1;
		System.out.println("Start time: "+time);
		long rupTime = Math.round((startYear-1970.0)*ProbabilityModelsCalc.MILLISEC_PER_YEAR);
		System.out.println("Rup time: "+rupTime);

		GregorianCalendar cal = new GregorianCalendar(TimeZone.getTimeZone("UTC"));
		cal.set(2012, 0, 1, 0, 0, 0);
		time = cal.getTimeInMillis();
		System.out.println("New cal test");
		System.out.println("Cal millis: "+cal.getTimeInMillis());
		System.out.println("Formatted:");
		System.out.println("\tyyyy\tMM\tdd\tHH\tmm\tss.SSS");
		System.out.println("\t"+catDateFormat.format(cal.getTime()));

		System.out.println();
		for (startYear=1970; startYear<1980; startYear++) {
			//			time = Math.round((startYear-1970.0)*ProbabilityModelsCalc.MILLISEC_PER_YEAR);
			double inner = (startYear-1970.0)*ProbabilityModelsCalc.MILLISEC_PER_YEAR;
			time = Math.round(inner);
			//			time = (long)((long)startYear - 1970l) * (long)ProbabilityModelsCalc.MILLISEC_PER_YEAR;
			boolean equal = inner == time;

			System.out.println(startYear+" epoch="+time+", inner double to be rounded: "+inner+". Equal? "+equal);
			//			System.out.println(time);
			System.out.println("\t"+catDateFormat.format(new Date(time)));
		}
	}

	private static void test10() throws ZipException, IOException, DocumentException {
		FaultSystemRupSet rupSet = U3FaultSystemIO.loadRupSet(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		int id = 193821;
		System.out.println(Joiner.on(",").join(rupSet.getSectionsIndicesForRup(id)));
	}

	private static void test11() throws IOException, DocumentException {
		FaultSystemSolution sol = U3FaultSystemIO.loadSol(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));

		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		// Poisson
		erf.getParameter(ProbabilityModelParam.NAME).setValue(ProbabilityModelOptions.POISSON);
		// UCERF3 TD
		//		erf.getParameter(ProbabilityModelParam.NAME).setValue(ProbabilityModelOptions.U3_PREF_BLEND);
		//		erf.getTimeSpan().setStartTime(2014); // start year
		//		erf.setParameter(HistoricOpenIntervalParam.NAME, // historical open interval
		//				erf.getTimeSpan().getStartTimeYear()-1875d);

		// duration - only really necessary for TD calculations, as MFD later is annualized
		erf.getTimeSpan().setDuration(30d);

		erf.updateForecast();

		// circular region with given center and radius
		Region reg = new Region(new Location(34, -118), 100d);
		// rectangular region with these corners
		//		Region reg = new Region(new Location(34, -118), new Location(35, -120));

		double minMag = 5d;
		int numMag = 51;
		double deltaMag = 0.1;

		// Participation MFD
		IncrementalMagFreqDist mfd = ERF_Calculator.getParticipationMagFreqDistInRegion(
				erf, reg, minMag, numMag, deltaMag, true);
		// Nucleation MFD - this will be slower
		//		IncrementalMagFreqDist mfd = ERF_Calculator.getMagFreqDistInRegion(erf, reg, minMag, numMag, deltaMag, true);

		System.out.println(mfd);
	}

	private static void test12() {
		AttenuationRelationship attenRel;
		Site site;
		EqkRupture rup;
		FaultTrace faultTrace = new FaultTrace("Seismic");

		rup = new EqkRupture();
		rup.setMag(9.3D);
		rup.setAveRake(90.0D);


		faultTrace.add( new Location(46.04674444, -123.9145741, 0));
		faultTrace.add( new Location(47.97927793, -123.9145741, 0));

		rup.setRuptureSurface(new StirlingGriddedSurface(faultTrace,45D,0D, 10D,1.0D));

		attenRel = new USGS_Combined_2004_AttenRel(null);
		attenRel.setParamDefaults();
		attenRel.setIntensityMeasure("PGA");

		site = new Site();
		Location loc = new Location(46.347669, -119.278117);
		site.setLocation(loc);
		site.addParameter(attenRel.getParameter("Vs30"));

		attenRel.setEqkRupture(rup);
		attenRel.setSite(site);

		//        attenRel.getExceedProbabilities(intensityMeasureLevels)

		DiscretizedFunc imls = IMT_Info.getUSGS_PGA_Function();
		HazardCurveSetCalculator.getLogFunction(imls);
		System.out.println(attenRel.getExceedProbabilities(imls));

		System.out.println( attenRel.getIML_AtExceedProb(0.9));

		//        pga = Math.exp( attenRel.getIML_AtExceedProb(prob) );
	}
	
	private static void test13() {
		double x1 = 0.5;
		double x2 = 1.0;
		
		System.out.println("geo-mean: "+StatUtils.geometricMean(new double[] {x1, x2}));
		System.out.println("log-mean: "+Math.exp(StatUtils.mean(new double[] {Math.log(x1), Math.log(x2)})));
	}
	
	private static void test14() {
		String name = "34.05_-118.57.txt";
//		String latlon = StringUtils.replaceChars(name, '_', ',');
		String latlon = StringUtils.replaceChars(StringUtils.substringBeforeLast(
				name, "."), '_', ',');
		System.out.println(latlon);
	}
	
	private static void test15() {
		MeanUCERF2 erf = new MeanUCERF2();
		erf.updateForecast();
		System.out.println("Orig width: "+erf.getSource(90).getSourceSurface().getAveWidth());
		erf.setParameter(MeanUCERF2.CYBERSHAKE_DDW_CORR_PARAM_NAME, true);
		erf.updateForecast();
		System.out.println("Corr width: "+erf.getSource(90).getSourceSurface().getAveWidth());
	}
	
	private static void test16() {
		System.out.println(MagUtils.magToMoment(7.8)/MagUtils.magToMoment(6.7));
	}
	
	private static void test17() {
		HashMap<double[], String> map = Maps.newHashMap();
		double[] array1 = {1d, 2d};
		double[] array2 = {2d, 3d};
		double[] array3 = {4d, 5d};
		map.put(array1, "array1");
		System.out.println(map.get(array1));
	}
	
	private static void test18() {
		double[] bVals = { -1d, 0d, 0.5d, 1d, 1.5d, 2d };
		
		List<DiscretizedFunc> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		List<Color> colors = GraphWindow.generateDefaultColors();
		
		for (int i=0; i<bVals.length; i++) {
			double b = bVals[i];
			DiscretizedFunc mfd = new GutenbergRichterMagFreqDist(b, 1d, 5d, 9d, 50);
			funcs.add(mfd);
			mfd.setName("b="+(float)b);
			Color c = colors.get(i % colors.size());
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, c));
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, "MFDs", "Mag", "Incr rate");
		spec.setLegendVisible(true);
		GraphWindow gw = new GraphWindow(spec);
		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
	}
	
	private static void test19() throws IOException, DocumentException {
		FaultSystemSolution sol = U3FaultSystemIO.loadSol(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		FaultModels fm = FaultModels.FM3_1;
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		RELM_SOCAL soCalReg = new CaliforniaRegions.RELM_SOCAL();
		
		double statewide = 0d;
		double soCal = 0d;
		double[] fractInSoCal = rupSet.getFractRupsInsideRegion(soCalReg, true);
		for (int r=0; r<rupSet.getNumRuptures(); r++) {
			if (rupSet.getMagForRup(r) >= 7) {
				double rate = sol.getRateForRup(r);
				statewide += rate;
				soCal += rate*fractInSoCal[r];
			}
		}
		
		
		Map<String, List<Integer>> map = fm.getNamedFaultsMapAlt();
		List<Integer> safSects = map.get("San Andreas");
		List<Integer> sjSects = map.get("San Jacinto (SB to C)");
		sjSects.addAll(map.get("San Jacinto (CC to SM)"));
		
		Map<Integer, FaultSection> sectMap = fm.fetchFaultSectionsMap();
		
//		System.out.println("SSAF sects before removal: "+safSects.size());
//		
//		for (int i=safSects.size(); --i>=0;) {
//			FaultSection sect = sectMap.get(safSects.get(i));
//			boolean inside = false;
//			for (Location loc : sect.getFaultTrace()) {
//				if (soCalReg.contains(loc)) {
//					inside = true;
//					break;
//				}
//			}
//			if (!inside)
//				safSects.remove(i);
//		}
//		
//		System.out.println("SSAF sects after removal: "+safSects.size());
		
		System.out.println("Statewide: "+statewide+" (RI: "+(1d/statewide)+")");
		System.out.println("SoCal: "+soCal+" (RI: "+(1d/soCal)+")");
		
		double ssaf = partProbForParents(sol, 7d, safSects, fractInSoCal);
		double sj = partProbForParents(sol, 7d, sjSects, fractInSoCal);
		System.out.println("sSAF: "+ssaf+" (RI: "+(1d/ssaf)+")");
		System.out.println("SJF: "+sj+" (RI: "+(1d/sj)+")");
	}
	
	private static double partProbForParents(FaultSystemSolution sol, double mag, List<Integer> parentSects, double[] fractInSoCal) {
		HashSet<Integer> rups = new HashSet<Integer>();
		FaultSystemRupSet rupSet = sol.getRupSet();
		for (int parent : parentSects)
			rups.addAll(rupSet.getRupturesForParentSection(parent));
		double sum = 0d;
		for (int r : rups)
			if (rupSet.getMagForRup(r) >= mag)
				sum += sol.getRateForRup(r)*fractInSoCal[r];
		return sum;
	}
	
	private static void test20() throws IOException {
//		File file = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
//				+ "2016_02_19-mojave_m7-10yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14-combined100k/"
//				+ "results_descendents_m5_preserve.bin");
//		File file = new File("/tmp/results_descendents.bin");
		File file = new File("/tmp/results_descendents_m5_preserve.bin");
//		File file = new File("/tmp/results_m5_preserve.bin");
		int scenID = 9893;
		
//		List<ETAS_EqkRupture> catalog = ETAS_CatalogIO.getBinaryCatalogsIterable(file, 0d).iterator().next();
//		checkCatalog(catalog, scenID);
//		System.out.println("DONE");
		
		List<? extends List<ETAS_EqkRupture>> catalogs = ETAS_CatalogIO.loadCatalogsBinary(file);
		for (int i=0; i<catalogs.size(); i++) {
			checkCatalog(catalogs.get(0), scenID);
		}
		System.out.println("DONE");
	}
	
	private static void checkCatalog(List<ETAS_EqkRupture> catalog, int scenID) {
		List<ETAS_EqkRupture> children = ETAS_SimAnalysisTools.getChildrenFromCatalog(catalog, scenID);
		if (catalog.size() != children.size()) {
			System.out.println("*** Original ***");
			printCatalog(catalog);
			System.out.println("*** Children ***");
			printCatalog(children);
		}
		Preconditions.checkState(catalog.size() == children.size(),
				"not the same! %s != %s", catalog.size(), children.size());
		// now test filter preserve
		double minMag = 5d;
		List<ETAS_EqkRupture> childrenAbove = ETAS_SimAnalysisTools.getAboveMagPreservingChain(children, minMag);
		List<ETAS_EqkRupture> childrenAbove2 = ETAS_SimAnalysisTools.getChildrenFromCatalog(childrenAbove, scenID);
		Preconditions.checkState(childrenAbove.size() == childrenAbove2.size(),
				"not the same after above mag check! %s != %s", childrenAbove.size(), childrenAbove2.size());
//		printCatalog(childrenAbove2);
	}
	
	private static void printCatalog(List<ETAS_EqkRupture> catalog) {
		System.out.println("ID\tMag\tGen\tParent");
		for (ETAS_EqkRupture rup : catalog)
			System.out.println(rup.getID()+"\t"+(float)rup.getMag()+"\t"+rup.getGeneration()+"\t"+rup.getParentID());
	}
	
	private static void test21() {
		ArbitrarilyDiscretizedFunc func1 = new ArbitrarilyDiscretizedFunc();
		func1.set(1d, 2d);
		List<ArbitrarilyDiscretizedFunc> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		funcs.add(func1);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		PlotSpec spec1 = new PlotSpec(funcs, chars, "Spec 1", "x", "y");
		PlotSpec spec2 = new PlotSpec(funcs, chars, "Spec 2", "x", "y");
		GraphWindow gw = new GraphWindow(spec1, false);
		gw.setXLog(true);
		gw.addTab(spec2);
		gw.setVisible(true);
	}
	
	private static void test22() throws IOException, DocumentException {
//		File solFile = new File("/home/kevin/Documents/2016_SCEC_AM/ucerf3/FM3_1_ref_slip_high.zip");
		File solFile = new File("/home/kevin/Documents/2016_SCEC_AM/ucerf3/FM3_1_ref_paleo_high.zip");
//		File solFile = new File("/home/kevin/Documents/2016_SCEC_AM/ucerf3/FM3_1_ref.zip");
		InversionFaultSystemSolution sol = U3FaultSystemIO.loadInvSol(solFile);
		File outputDir = new File(solFile.getParentFile(), solFile.getName().replace(".zip", ""));
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		ArrayList<U3PaleoRateConstraint> paleoRateConstraints =
				CommandLineInversionRunner.getPaleoConstraints(sol.getRupSet().getFaultModel(), sol.getRupSet());
		System.out.println(paleoRateConstraints.size()+" paleo constraints");
		List<U3AveSlipConstraint> aveSlipConstraints = U3AveSlipConstraint.load(sol.getRupSet().getFaultSectionDataList());
		System.out.println(aveSlipConstraints.size()+" ave slip constraints");
		Map<String, List<Integer>> namedFaultsMap = sol.getRupSet().getFaultModel().getNamedFaultsMapAlt();
		for (String name : Lists.newArrayList(namedFaultsMap.keySet())) {
			if (!name.contains("ndreas"))
				namedFaultsMap.remove(name);
		}
		CommandLineInversionRunner.writePaleoFaultPlots(paleoRateConstraints, aveSlipConstraints, namedFaultsMap, sol, outputDir);
	}
	
	private static void test23() {
		double duration = 30d;
		double prob = 0.217;
		
		// prob = 1 - Math.exp(-rate*duration)
		// Math.exp(-rate*duration) = 1 - prob
		// -rate*duration = Math.log(1 - prob)
		// rate = -duration*Math.log(1 - prob)
		
		double rate = -Math.log(1 - prob)/duration;
		
		System.out.println("Annual rate: "+rate);
		
		double targetDuration = 7d/365.25;
		double targetRate = targetDuration*rate;
		double targetProb = 1 - Math.exp(-rate*targetDuration);
		
		System.out.println("1-week rate: "+targetRate);
		System.out.println("1-week prob: "+targetProb);
		
		
	}
	
	private static void test24() {
//		double probability = 0.21707605530773;
		double probability = 0.113;
		double duration = 30d;
		double annualRate = -Math.log(1 - probability)/duration;
		
		double weeklyRate = annualRate *7d/365.25;
		System.out.println("Annual rate: "+annualRate);
		System.out.println("Weekly rate: "+weeklyRate);
	}
	
//	private static void test25() {
//		DiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
//		func.set(1.0E-4, 3.65255E-4);
//		func.set(1.3E-4,  3.457038E-4);
//		func.set(1.6E-4,  3.2601133E-4);
//		func.set(2.0E-4,  3.011314E-4);
//		func.set(2.5E-4,  2.733292E-4);
//		func.set(3.2E-4,  2.4073821E-4);
//		func.set(4.0E-4,  2.111737E-4);
//		func.set(5.0E-4,  1.8285471E-4);
//		func.set(6.3E-4,  1.5588816E-4);
//		func.set(7.9E-4,  1.3244183E-4);
//		func.set(0.001,   1.11381836E-4);
//		func.set(0.00126, 9.397233E-5);
//		func.set(0.00158, 7.973755E-5);
//		func.set(0.002,   6.7452595E-5);
//		func.set(0.00251, 5.764316E-5);
//		func.set(0.00316, 4.932694E-5);
//		func.set(0.00398, 4.2291384E-5);
//		func.set(0.00501, 3.6282872E-5);
//		func.set(0.00631, 3.1062427E-5);
//		func.set(0.00794, 2.6521913E-5);
//		func.set(0.01,    2.254247E-5);
//		func.set(0.01259, 1.9099985E-5);
//		func.set(0.01585, 1.6150763E-5);
//		func.set(0.01995, 1.3649939E-5);
//		func.set(0.02512, 1.15329385E-5);
//		func.set(0.03162, 9.743127E-6);
//		func.set(0.03981, 8.220153E-6);
//		func.set(0.05012, 6.9222365E-6);
//		func.set(0.0631,  5.814752E-6);
//		func.set(0.07943, 4.8654542E-6);
//		func.set(0.1,     4.03794E-6);
//		func.set(0.12589, 3.2988742E-6);
//		func.set(0.15849, 2.6238315E-6);
//		func.set(0.19953, 2.0062375E-6);
//		func.set(0.25119, 1.4563641E-6);
//		func.set(0.31623, 9.925378E-7);
//		func.set(0.39811, 6.2928785E-7);
//		func.set(0.50119, 3.6850898E-7);
//		func.set(0.63096, 1.9821493E-7);
//		func.set(0.79433, 9.7517116E-8);
//		func.set(1.0,     4.3739632E-8);
//		func.set(1.25893, 1.7841097E-8);
//		func.set(1.58489, 6.6055765E-9);
//		func.set(1.99526, 2.2166786E-9);
//		func.set(2.51189, 6.7364303E-10);
//		func.set(3.16228, 1.8535695E-10);
//		func.set(3.98107, 4.619971E-11);
//		func.set(5.01187, 1.04463105E-11);
//		func.set(6.30957, 2.1446178E-12);
//		func.set(7.94328, 3.9890313E-13);
//		func.set(10.0,    6.550316E-14);
//		AbstractMCErProbabilisticCalc.calcRTGM(func);
//	}
	
	private static void test26() {
//		double p = 1.1E-4;
//		double p = 1.110223E-16;
		double p = 1.2000000000000899E-4;
		int n = 500;
//		double p = 0;
//		int n = 100000;
		double[] ret = ETAS_Utils.getBinomialProportion95confidenceInterval(p, n);
		System.out.println("95: "+ret[0]+" "+ret[1]);
		ret = ETAS_Utils.getBinomialProportion68confidenceInterval(p, n);
		System.out.println("68: "+ret[0]+" "+ret[1]);
		
		long l = (1325419200000l + 7l*ProbabilityModelsCalc.MILLISEC_PER_DAY);
		System.out.println(l);
		ETAS_SimAnalysisTools.writeMemoryUse("test");
	}
	
//	private static void test27() throws IOException {
//		String[] names = { "s688", "s605", "s646" };
//		DBAccess db = Cybershake_OpenSHA_DBApplication.getDB();
//		CVM_Vs30 cvm = new CVM_Vs30(CVM.CVMS4i26);
//		
//		SiteInfo2DB sites2db = new SiteInfo2DB(db);
//		
//		for (String name : names) {
//			Location loc = sites2db.getSiteFromDB(name).createLocation();
//			double vs30 = cvm.getValue(loc);
//			System.out.println(name+" ("+(float)loc.getLatitude()+", "+(float)loc.getLongitude()+"): "+(float)vs30);
//		}
//		
//		db.destroy();
//	}
	
	private static void test28() throws IOException {
		WillsMap2006 wills2006 = new WillsMap2006();
		WillsMap2015 wills2015 = new WillsMap2015();
		
		GriddedRegion reg = new CaliforniaRegions.RELM_TESTING_GRIDDED(0.01);
		
		MinMaxAveTracker track2006 = new MinMaxAveTracker();
		MinMaxAveTracker track2015 = new MinMaxAveTracker();
		
		for (double val : wills2006.getValues(reg.getNodeList()))
			track2006.addValue(val);
		for (double val : wills2015.getValues(reg.getNodeList()))
			track2015.addValue(val);
		
		System.out.println(track2006);
		System.out.println(track2015);
	}
	
	private static void test29() throws IOException {
		File dir = new File("/home/kevin/OpenSHA/UCERF3/etas/hazard/"
				+ "2017_03_03-haywired_m7_fulltd_descendents-NGA2-0.02-site-effects-with-basin");
		for (File file : dir.listFiles()) {
			String name = file.getName();
			if (!file.getName().endsWith(".bin"))
				continue;
			
			name = name.replaceAll("results_results", "results");
			name = name.replaceAll("pga_combined", "pga");
			name = name.replaceAll("pgv_combined", "pgv");
			
			System.out.println(file.getName()+" => "+name);
			if (!name.equals(file.getName()))
				Files.move(file, new File(dir, name));
		}
	}
	
	private static void test30() {
		File catalogFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2017_01_02-haywired_m7-10yr-gridded-only-200kcombined/results_descendents_m5_preserve.bin");
		
		double minMag = 7d;
		
		int numCatsWith = 0;
		int totalNumAbove = 0;
		
		BinarayCatalogsIterable catsIt = ETAS_CatalogIO.getBinaryCatalogsIterable(catalogFile, minMag);
		
		double[] counts = new double[catsIt.getNumCatalogs()];
		int index = 0;
		
		for (List<ETAS_EqkRupture> catalog : catsIt) {
			int count = catalog.size();
			if (count > 0)
				numCatsWith++;
			totalNumAbove += count;
			counts[index++] = count;
		}
		
		double fractCatsWith = (double)numCatsWith/counts.length;
		double mean = (double)totalNumAbove/counts.length;
		double median = DataUtils.median(counts);
		
		System.out.println(numCatsWith+"/"+counts.length+", "+(float)(fractCatsWith*100d)
				+" % of catalogs have at least one M"+(float)minMag);
		System.out.print("Mean: "+mean+", Median: "+median);
	}
	
	private static void test31() throws IOException, DocumentException {
		int bendID = 287;
		int mojaveNID = 286;
		int mojaveSID = 301;
		
		FaultSystemSolution sol = U3FaultSystemIO.loadSol(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		HashSet<Integer> bendRups = new HashSet<Integer>(rupSet.getRupturesForParentSection(bendID));
		HashSet<Integer> mojaveNRups = new HashSet<Integer>(rupSet.getRupturesForParentSection(mojaveNID));
		HashSet<Integer> mojaveSRups = new HashSet<Integer>(rupSet.getRupturesForParentSection(mojaveSID));
		
		List<Integer> bendAndMojaveNRups = Lists.newArrayList();
		List<Integer> bendAndMojaveSRups = Lists.newArrayList();
		for (Integer rup : bendRups) {
			if (mojaveNRups.contains(rup))
				bendAndMojaveNRups.add(rup);
			if (mojaveSRups.contains(rup))
				bendAndMojaveSRups.add(rup);
		}
		
		double minMag = 7.3;
		
		double rateBend = rateForRups(sol, bendRups, minMag);
		double rateBenMojaveN = rateForRups(sol, bendAndMojaveNRups, minMag);
		double rateBenMojaveS = rateForRups(sol, bendAndMojaveSRups, minMag);
		
		System.out.println("Rate Big Bend: "+rateBend);
		System.out.println("Rate Big Bend & Mojave N: "+rateBenMojaveN+" ("+(float)(100d*rateBenMojaveN/rateBend)+" %)");
		System.out.println("Rate Big Bend & Mojave S: "+rateBenMojaveS+" ("+(float)(100d*rateBenMojaveS/rateBend)+" %)");
	}
	
	private static void test32() throws IOException, DocumentException {
		int niOffshoreID = 122;
		int niOnshoreID = 235;
		int rcID = 123;
		
		FaultSystemSolution sol = U3FaultSystemIO.loadSol(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		HashSet<Integer> niOffshoreRups = new HashSet<Integer>(rupSet.getRupturesForParentSection(niOffshoreID));
		HashSet<Integer> niOnshoreRups = new HashSet<Integer>(rupSet.getRupturesForParentSection(niOnshoreID));
		HashSet<Integer> rcRups = new HashSet<Integer>(rupSet.getRupturesForParentSection(rcID));
		
		List<Integer> rcAndOffshoreRups = Lists.newArrayList();
		List<Integer> rcAndOnshoreRups = Lists.newArrayList();
		for (Integer rup : rcRups) {
			if (niOffshoreRups.contains(rup))
				rcAndOffshoreRups.add(rup);
			if (niOnshoreRups.contains(rup))
				rcAndOnshoreRups.add(rup);
		}
		
		double[] minMags = { 6.5, 7d, 7.5d };
		
		for (double minMag : minMags) {
			double rateRC = rateForRups(sol, rcRups, minMag);
			double rateOffshore = rateForRups(sol, niOffshoreRups, minMag);
			double rateOnshore = rateForRups(sol, niOnshoreRups, minMag);
			double rateRCandOffshore = rateForRups(sol, rcAndOffshoreRups, minMag);
			double rateRCandOnshore = rateForRups(sol, rcAndOnshoreRups, minMag);
			
			System.out.println("Min mag: "+minMag);
			System.out.println("\tRate RC: "+rateRC);
			System.out.println("\tRate NI Offshore: "+rateOffshore);
			System.out.println("\tRate NI Onshore: "+rateOnshore);
			System.out.println("\tRate RC and NI Offshore: "+rateRCandOffshore+" ("+(float)(100d*rateRCandOffshore/rateRC)+" %)");
			System.out.println("\tRate RC and NI Onshore: "+rateRCandOnshore+" ("+(float)(100d*rateRCandOnshore/rateRC)+" %)");
		}
	}
	
	private static double rateForRups(FaultSystemSolution sol, Collection<Integer> rups, double minMag) {
		double rate = 0d;
		FaultSystemRupSet rupSet = sol.getRupSet();
		for (int rup : rups)
			if (rupSet.getMagForRup(rup) >= minMag)
				rate += sol.getRateForRup(rup);
		return rate;
	}
	
	private static void test33() throws IOException, DocumentException {
		FaultSystemSolution sol1 = U3FaultSystemIO.loadSol(new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		FaultSystemSolution sol2 = U3FaultSystemIO.loadSol(new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
				+ "UCERF3_ERF/cached_FM3_1_dep100.0_depMean_rakeMean.zip"));
		int parentID = 295;
		List<Integer> rups1 = sol1.getRupSet().getRupturesForParentSection(parentID);
//		Collections.sort(rups1);
		List<Integer> rups2 = sol2.getRupSet().getRupturesForParentSection(parentID);
//		Collections.sort(rups2);
		System.out.println("Hello?");
		System.out.println("Sol 1 "+rups1.size()+" rups: ");
		System.out.println("HI?");
		System.out.println("Sol 2 "+rups2.size()+" rups: ");
	}
	
	private static void test34() throws IOException, DocumentException {
		FaultSystemRupSet rupSet = U3FaultSystemIO.loadRupSet(new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		int index = TestScenario.MOJAVE_M7.getFSS_Index();
		for (FaultSection sect : rupSet.getFaultSectionDataForRupture(index))
			System.out.println(sect.getSectionId()+": "+sect.getSectionName());
	}
	
	private static void test35() throws IOException, DocumentException {
		FaultSystemSolution sol = U3FaultSystemIO.loadSol(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));

		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		erf.setParameter(BackgroundRupParam.NAME, BackgroundRupType.CROSSHAIR);
		erf.updateForecast();
		System.out.println("Total Num Ruptures: "+erf.getTotNumRups());
	}
	
	private static void test36() {
		double rateM5 = TotalMag5Rate.RATE_7p9.getRateMag5();
		double b = 1d;
		GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(b, 1d, 5.05, 9.05, 41);
		EvenlyDiscretizedFunc cumGR = gr.getCumRateDistWithOffset();
		cumGR.scale(rateM5/cumGR.getY(0));
		System.out.println(cumGR);
	}
	
	private static void test37() throws IOException, DocumentException {
		FaultSystemSolution sol = U3FaultSystemIO.loadSol(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_2_MEAN_BRANCH_AVG_SOL.zip"));
//						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_TRUE_HAZARD_MEAN_SOL.zip"));
		System.out.println("Fault ruptures: "+sol.getRupSet().getNumRuptures());
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		erf.setParameter(BackgroundRupParam.NAME, BackgroundRupType.CROSSHAIR);
		erf.updateForecast();
		System.out.println("ERF rups from FSS: "+erf.getTotNumRupsFromFaultSystem());
		System.out.println("Gridded ERF rups: "+(erf.getTotNumRups() - erf.getTotNumRupsFromFaultSystem()));
	}
	
	private static void test38() throws IOException, DocumentException {
		FaultSystemSolution sol = U3FaultSystemIO.loadSol(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));

		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
//		// Poisson
//		erf.getParameter(ProbabilityModelParam.NAME).setValue(ProbabilityModelOptions.POISSON);
		// UCERF3 TD
		erf.getParameter(ProbabilityModelParam.NAME).setValue(ProbabilityModelOptions.U3_PREF_BLEND);
		erf.getTimeSpan().setStartTime(2017); // start year
		erf.setParameter(HistoricOpenIntervalParam.NAME, // historical open interval
				 erf.getTimeSpan().getStartTimeYear()-1875d);

		// duration - only really necessary for TD calculations, as MFD later is annualized
		erf.getTimeSpan().setDuration(1d);

		erf.updateForecast();
		
		double cholameProb = FaultSysSolutionERF_Calc.calcParticipationProbForParentSects(erf, 7.8, 285);
		System.out.println("Cholame prob: "+cholameProb);
		double ssafProb = FaultSysSolutionERF_Calc.calcParticipationProbForParentSects(erf, 7.8, 285, 300, 287, 286, 301, 282, 283, 284, 295);
		System.out.println("S.SAF prob: "+ssafProb);
	}
	
	private static void test39() throws IOException, DocumentException {
		FaultSystemSolution sol = U3FaultSystemIO.loadSol(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
						+ "FM3_1_GEOL_MEAN_BRANCH_AVG_SOL.zip"));

		IncrementalMagFreqDist supraMFD = new IncrementalMagFreqDist(5.05, 41, 0.1);
		supraMFD.setName("Fault Supra");
		IncrementalMagFreqDist subMFD = new IncrementalMagFreqDist(5.05, 51, 0.1);
		subMFD.setName("Fault Sub");
		IncrementalMagFreqDist offMFD = new IncrementalMagFreqDist(5.05, 51, 0.1);
		offMFD.setName("Off Fault");
		
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		for (int r=0; r<rupSet.getNumRuptures(); r++)
			supraMFD.add(supraMFD.getClosestXIndex(rupSet.getMagForRup(r)), sol.getRateForRup(r));
		
		GridSourceProvider gridProv = sol.getGridSourceProvider();
		
		for (int n=0; n<gridProv.size(); n++) {
			IncrementalMagFreqDist nodeOffMFD = gridProv.getNodeUnassociatedMFD(n);
			IncrementalMagFreqDist nodeSubSeisMFD = gridProv.getNodeSubSeisMFD(n);
			
			if (nodeOffMFD != null) {
				for (int i=0; i<nodeOffMFD.size(); i++) {
					if (nodeOffMFD.getX(i) >= offMFD.getMinX())
						offMFD.add(nodeOffMFD.getX(i), nodeOffMFD.getY(i));
				}
			}
			if (nodeSubSeisMFD != null) {
				for (int i=0; i<nodeSubSeisMFD.size(); i++) {
					if (nodeSubSeisMFD.getX(i) >= offMFD.getMinX())
						subMFD.add(nodeSubSeisMFD.getX(i), nodeSubSeisMFD.getY(i));
				}
			}
		}
		
		SummedMagFreqDist totalMFD = new SummedMagFreqDist(supraMFD.getMinX(), supraMFD.getMaxX(), supraMFD.size());
		totalMFD.setName("Total (on+off)");
		SummedMagFreqDist faultMFD = new SummedMagFreqDist(supraMFD.getMinX(), supraMFD.getMaxX(), supraMFD.size());
		faultMFD.setName("Fault Sub+Supra");
		
		totalMFD.addIncrementalMagFreqDist(supraMFD);
		totalMFD.addIncrementalMagFreqDist(subMFD);
		totalMFD.addIncrementalMagFreqDist(offMFD);
		
		faultMFD.addIncrementalMagFreqDist(supraMFD);
		faultMFD.addIncrementalMagFreqDist(subMFD);
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		funcs.add(offMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
		
		funcs.add(subMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GREEN));
		
		funcs.add(supraMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		
		funcs.add(faultMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		
		funcs.add(totalMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "FM3.1 Geol MFDs", "Magnitude", "Incremental Rate (1/yr)");
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setBackgroundColor(Color.WHITE);
		gp.setTickLabelFontSize(18);
		gp.setAxisLabelFontSize(20);
		gp.setPlotLabelFontSize(21);
		
		String prefix = "/tmp/fm3_1_geol_mfds";
		
		gp.drawGraphPanel(spec, false, true, new Range(5d, 9d), new Range(1e-8, 1e1));
		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPNG(prefix+".png");
		gp.saveAsPDF(prefix+".pdf");
		gp.saveAsTXT(prefix+".txt");
	}
	
	private static void test40() {
		int num = 100;
		Random r = new Random();
		for (int i=0; i<num; i++) {
//			long seed = System.nanoTime()*r.nextLong();
			long seed = System.nanoTime()+r.nextLong();
			System.out.println(seed);
		}
	}
	
	private static void test41() {
		double strike = 90;
		double dip = 90;
		double rake = 90;
		double[] vect = PlaneUtils.getSlipVector(new double[] {strike, dip, rake});
		System.out.println(vect[0]+" "+vect[1]+" "+vect[2]);
	}
	
	private static void test42() throws IOException {
		FileUtils.createZipFile(new File("/tmp/bbp_test7.zip"), new File("/tmp/bbp_test7"), false);
	}
	
	private static void test43() throws IOException {
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, 10d);
		double fract = 0.7;
		for (CPTVal cv : cpt) {
			Color c1 = cv.minColor;
			Color c2 = cv.maxColor;
			c1 = new Color((int)(c1.getRed()*fract), (int)(c1.getGreen()*fract), (int)(c1.getBlue()*fract));
			c2 = new Color((int)(c2.getRed()*fract), (int)(c2.getGreen()*fract), (int)(c2.getBlue()*fract));
			cv.minColor = c1;
			cv.maxColor = c2;
		}
		cpt.setBelowMinColor(cpt.getMinColor());
		cpt.setAboveMaxColor(cpt.getMaxColor());
		cpt.setNanColor(Color.LIGHT_GRAY);
		cpt.writeCPTFile(new File("/tmp/dark_cpt.cpt"));
	}
	
	private static void test44() {
		double r1 = 0.02;
		double t1 = 50;
		
		double r1Star = r1*(1 + 0.5*r1);
		System.out.println("r1 = "+r1);
		System.out.println("r1* = "+r1Star);
		System.out.println(-Math.log(1-r1));
		System.out.println("T1 = "+t1);
		
		double[] rps = { 1d, 50d, 500d, 1547.0297, 1733, 2474, 2500 };
		
		for (double t2 : rps) {
			double r2Star = r1Star*t2/t1;
			
			double r2 = Math.sqrt(1 + 2*r2Star) - 1;
			
			System.out.println("T2="+(float)t2+"\tR2*="+(float)r2Star+"\tR2="+(float)r2);
		}
		
		System.out.println("Median test with fixed r2");
		double r2 = 0.5;
		double r2Star = r2*(1 + 0.5*r2);
		System.out.println("r2="+(float)r2+"\tr2*="+r2Star);
		double t2 = t1*r2Star/r1Star;
		System.out.println("t2="+(float)t2);
		
		System.out.println("Median test with fixed r2*");
		r2Star = 0.5;
		r2 = Math.sqrt(1 + 2d*r2Star) - 1;
		System.out.println("r2="+(float)r2+"\tr2*="+r2Star);
		t2 = t1*r2Star/r1Star;
		System.out.println("t2="+(float)t2);
	}
	
	private static void test45() {
//		System.out.println(Double.parseDouble("NaN"));
//		Map<String, String> env = System.getenv();
//		for (String key : env.keySet())
//			System.out.println(key+": "+env.get(key));
		
//		TemporalAccessor date = DateTimeFormatter.ISO_LOCAL_DATE_TIME.parse("2018-02-22T15:04:59");
		LocalDateTime future = LocalDateTime.parse("2018-02-23T15:22:59", DateTimeFormatter.ISO_LOCAL_DATE_TIME);
		LocalDateTime now = LocalDateTime.now();
		Duration duration = Duration.between(now, future);
		long secs = duration.get(ChronoUnit.SECONDS);
		System.out.println("Duration: "+secs+" s");
//		date.
//		date.to
//		date.get(TemporalField)
//		System.out.println(date);
//		List<Integer> list1 = new ArrayList<>();
//		List<Integer> list2 = new ArrayList<>();
//		
//		list1.add(0);
//		list1.add(1);
//		list1.add(4);
//		list1.add(5);
//		list1.add(8);
//		
//		list2.add(0);
//		list2.add(1);
//		list2.add(4);
//		list2.add(5);
//		list2.add(8);
//		
//		System.out.println("Hash1: "+list1.hashCode());
//		System.out.println("Hash2: "+list2.hashCode());
//		
//		System.out.println("Equals? "+list1.equals(list2));
	}
	
	private static void test46() {
		double mean = 1d;
		double sd = 0.3;
		LogNormalDistribution ln = new LogNormalDistribution(Math.log(mean), sd);
		
		SummaryStatistics stats = new SummaryStatistics();
		SummaryStatistics logStats = new SummaryStatistics();
		
		for (int i=0; i<10000000; i++) {
			double sample = ln.sample();
			stats.addValue(sample);
			logStats.addValue(Math.log(sample));
		}
		
		System.out.println("Linear Stats");
		System.out.println(stats);
		System.out.println("\nLog Stats");
		System.out.println(logStats);
	}
	
	private static void test47() throws IOException {
		File file = new File("/tmp/catalog_000.bin");
		List<ETAS_EqkRupture> catalog = ETAS_CatalogIO.loadCatalogBinary(file);
		
		System.out.println("Loaded "+catalog.size()+" events");
		
		int numFSS = 0;
		int numGridded = 0;
		
		for (ETAS_EqkRupture rup : catalog) {
			if (rup.getFSSIndex() >= 0)
				numFSS++;
			if (rup.getGridNodeIndex() >= 0)
				numGridded++;
			if (rup.getFSSIndex() >= 0 && rup.getGridNodeIndex() >= 0)
				System.out.println("Uh oh, has both indexes!");
			if (rup.getFSSIndex() < 0 && rup.getGridNodeIndex() < 0)
				System.out.println("Uh oh, has neither indexes!");
		}
		
		System.out.println(numFSS+" with FSS index");
		System.out.println(numGridded+" with grid node index");
	}
	
	private static void test48() throws IOException {
		Location laHabra = new Location(33.9225, -117.9352);
		Location eq = new Location(33.844, -119.716);
		
		System.out.println("Distance: "+LocationUtils.horzDistance(laHabra, eq));
	}
	
	private static void test49() throws IOException {
		File profileFile = new File(System.getProperty("user.home")+"/.bash_profile");
		System.out.println("exists? "+profileFile.exists());
		System.exit(0);
		int num2 = 0;
		int num3plus = 0;
		int max = 0;
		
		List<? extends FaultSection> sects = RSQSimUtils.getUCERF3SubSectsForComparison(FaultModels.FM3_1, DeformationModels.GEOLOGIC);
		
		for (FaultSection sect : sects) {
			int num = sect.getFaultTrace().size();
			if (num == 2)
				num2++;
			else
				num3plus++;
			if (num > max)
				max = num;
		}
		
		System.out.println(num3plus+"/"+(num2+num3plus)+" have 3 or more points");
		System.out.println("Most points in a sub sect: "+max);
	}
	
	private static void test50() {
		int n = 1000000;
		NormalDistribution stdNorm = new NormalDistribution(0d, 1d);
		double[] offsets = { 0d, 0.5, -0.5, 1d, 2d };
		
		for (double offset : offsets) {
			Variance var = new Variance(true);
			for (int i=0; i<n; i++) {
				double val = stdNorm.sample()+offset;
				var.increment(val);
			}
			
			System.out.println("Offset: "+(float)offset);
			double variance = var.getResult();
			System.out.println("\tVariance: "+(float)variance);
			System.out.println("\tStd. Dev.: "+(float)Math.sqrt(variance));
		}
	}
	
	private static void test51() {
		System.out.println(MagUtils.momentToMag(6.528036E20));
	}
	
	private static void test52() {
		DoubleParameter param = new DoubleParameter("");
		if (param.getEditor() instanceof AbstractParameterEditor) {
			AbstractParameterEditor<?> editor = (AbstractParameterEditor<?>)param.getEditor();
			// update title font
			editor.setTitleFont(new Font(Font.SANS_SERIF, Font.BOLD, 18));
			if (editor.getWidget() instanceof NumericTextField) {
				NumericTextField textField = (NumericTextField)editor.getWidget();
				textField.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 12));
			}
		}
	}
	
	private static void test53() throws IOException, DocumentException {
//		FaultSystemSolution sol = FaultSystemIO.loadSol(
//				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
//						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		MeanUCERF3 erf = new MeanUCERF3();
		erf.setPreset(Presets.FM3_1_BRANCH_AVG);
		erf.updateForecast();
		FaultSystemSolution sol = erf.getSolution();
		int index = 44602;
		for (FaultSection sect : sol.getRupSet().getFaultSectionDataForRupture(index)) {
			System.out.println("Section: "+sect.getName());
		}
		RuptureSurface surf = sol.getRupSet().getSurfaceForRupture(index, 1d);
		System.out.println("Ztor: "+surf.getAveRupTopDepth());
		Location loc = new Location(37.78849, -122.26912);
		System.out.println("Rjb: "+surf.getDistanceJB(loc));
		System.out.println("Rrup: "+surf.getDistanceRup(loc));
		System.out.println("Rx: "+surf.getDistanceX(loc));
		double minDist = Double.POSITIVE_INFINITY;
		double minSurfDist = Double.POSITIVE_INFINITY;
		Location closestLoc = null;
		for (Location pt : surf.getEvenlyDiscritizedListOfLocsOnSurface()) {
			double dist = LocationUtils.linearDistance(pt, loc);
			if (dist < minDist) {
				closestLoc = pt;
				minDist = dist;
			}
			minSurfDist = Double.min(minSurfDist, LocationUtils.horzDistance(loc, pt));
		}
		System.out.println("3-D min dist: "+minDist);
		System.out.println("Closest: "+closestLoc);
		System.out.println("2-D min dist: "+minSurfDist);
	}
	
	private static void test54() throws IOException, DocumentException {
		FaultSystemSolution sol = U3FaultSystemIO.loadSol(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		Location mojave = new Location(34.42295, -117.80177, 5.8);
		Location bombay = new Location(33.3172, -115.72800000000001, 5.96);
		double minMojave = Double.POSITIVE_INFINITY;
		double minBombay = Double.POSITIVE_INFINITY;
		for (FaultSection sect : sol.getRupSet().getFaultSectionDataList()) {
			if (!sect.getName().contains("Mojave") && !sect.getName().contains("Coachella"))
				continue;
			RuptureSurface surf = sect.getFaultSurface(0.1, false, false);
			for (Location loc : surf.getEvenlyDiscritizedListOfLocsOnSurface()) {
				minMojave = Double.min(minMojave, LocationUtils.linearDistanceFast(mojave, loc));
				minBombay = Double.min(minBombay, LocationUtils.linearDistanceFast(bombay, loc));
			}
		}
		System.out.println("Mojave dist: "+minMojave);
		System.out.println("Bombay dist: "+minBombay);
	}
	
	private static void test55() throws IOException {
		List<U3PaleoRateConstraint> biasiScharerSites = new ArrayList<>();
		for (U3PaleoRateConstraint constraint : UCERF3_PaleoRateConstraintFetcher.getConstraints()) {
			String name = constraint.getPaleoSiteName();
			if (name.equals("S. San Andreas - Coachella") || name.equals("San Jacinto - Hog Lake")
					|| name.equals("Frazier Mountian, SSAF") || name.equals("N. San Andreas - Santa Cruz Seg.")
					|| name.equals("Hayward fault - South"))
				biasiScharerSites.add(constraint);
		}
		for (U3PaleoRateConstraint constraint : biasiScharerSites)
			System.out.println(constraint.getPaleoSiteName()+": "+(float)(1d/constraint.getMeanRate()));
		Preconditions.checkState(biasiScharerSites.size() == 5);
	}
	
	private static void test56() {
		GsonBuilder builder = new GsonBuilder();
		builder.setPrettyPrinting();
		Gson gson = builder.create();
		
		Map<String, String> map1 = new HashMap<>();
		map1.put("key1", "val1");
		map1.put("key2", "multi\n\nline\nval2");
		
		String map1Str = gson.toJson(map1);
		System.out.println(map1Str);
		Map<String, String> map1Recreate = gson.fromJson(map1Str, Map.class);
		for (String key : map1Recreate.keySet())
			System.out.println(key+": "+map1Recreate.get(key));
		
		System.out.println();
		Map<String, List<String>> map2 = new HashMap<>();
		map2.put("key1", Lists.newArrayList("val1"));
		map2.put("key2", Lists.newArrayList("multi", "", "line", "val2"));
		
		String map2Str = gson.toJson(map2);
		System.out.println(map2Str);
		Map<String, List<String>> map2Recreate = gson.fromJson(map2Str, Map.class);
		for (String key : map2Recreate.keySet())
			System.out.println(key+": "+map2Recreate.get(key));
	}
	
	private static void test57() throws IOException {
		CPT cpt = GMT_CPT_Files.BLACK_RED_YELLOW_UNIFORM.instance().reverse().rescale(0d, 5d);
		cpt.setNanColor(Color.LIGHT_GRAY);
		cpt.writeCPTFile(new File("/tmp/cpt.cpt"));
	}
	
	private static void test58() throws IOException {
		double mean = 5d;
		
		int num = 10000;
		double[] rnds = { 0.1d, 0.5d, 1d, 2d, 3d, 4d, 5d };
		
		for (double rnd : rnds) {
			double[] vals = new double[num];
			double[] diffs = new double[num];
			double rms = 0d;
			double rmsDiff = 0d;
			for (int i=0; i<num; i++) {
				double r = Math.random() - 0.5;
				vals[i] = mean + rnd*r;
				rms += vals[i]*vals[i];
				rmsDiff += (vals[i]-mean)*(vals[i]-mean);
				diffs[i] = Math.abs(vals[i] - mean);
			}
			rms /= (double)num;
			rms = Math.sqrt(rms);
			rmsDiff /= (double)num;
			rmsDiff = Math.sqrt(rmsDiff);
			
			System.out.println("Rnd="+(float)rnd+"\tmean="+(float)StatUtils.mean(vals)
				+"\tgeoMean="+(float)StatUtils.geometricMean(vals)+"\trms="+(float)rms
				+"\trmsDiff="+(float)rmsDiff+"\taveDiff="+(float)StatUtils.mean(diffs)
				+"\tstdDev="+(float)Math.sqrt(StatUtils.variance(vals)));
		}
		
		double momentNM = 1e20;
		double mag1 = (Math.log10(momentNM) - 9.05) / 1.5;
		System.out.println("mag1: "+mag1);
		double momentDC = momentNM*1e7;
		double mag2 = (2d/3d)*Math.log10(momentDC) - 10.7;
		System.out.println("mag2: "+mag2);
	}
	
	private static void test59() throws IOException {
		List<? extends FaultSection> subSects31 = DeformationModels.loadSubSects(FaultModels.FM3_1, DeformationModels.GEOLOGIC);
		List<? extends FaultSection> subSects32 = DeformationModels.loadSubSects(FaultModels.FM3_2, DeformationModels.GEOLOGIC);
		FaultSectionDataWriter.writeSectionsToFile(subSects31, null, new File("/tmp/fm_3_1.txt"), false);
		FaultSectionDataWriter.writeSectionsToFile(subSects32, null, new File("/tmp/fm_3_2.txt"), false);
	}
	
	private static void test60() throws IOException {
		System.out.println(ASK_2014.calcZ1ref(500));
	}
	
	private static void test61() {
		List<? extends FaultSection> subSects = DeformationModels.loadSubSects(
				FaultModels.FM3_1, DeformationModels.GEOLOGIC);
		
		List<FaultSection> gv07Sects = new ArrayList<>();
		List<FaultSection> gv08Sects = new ArrayList<>();
		List<FaultSection> gv10Sects = new ArrayList<>();
		List<FaultSection> ortigalitaSects = new ArrayList<>(); //Ortigalita
		
		for (FaultSection sect : subSects) {
			if (sect.getParentSectionId() == 136)
				gv07Sects.add(sect);
			if (sect.getParentSectionId() == 137)
				gv08Sects.add(sect);
			if (sect.getParentSectionId() == 138)
				gv10Sects.add(sect);
			if (sect.getParentSectionId() == 921 || sect.getParentSectionId() == 9) 
				ortigalitaSects.add(sect);
		}

		double min07Dist = Double.POSITIVE_INFINITY;
		for (FaultSection gvSect : gv07Sects)
			for (FaultSection oSect : ortigalitaSects)
				min07Dist = Double.min(min07Dist, minDist(gvSect, oSect));
		System.out.println("GV 07 to O dist: "+min07Dist);
		
		double min08Dist = Double.POSITIVE_INFINITY;
		for (FaultSection gvSect : gv08Sects)
			for (FaultSection oSect : ortigalitaSects)
				min08Dist = Double.min(min08Dist, minDist(gvSect, oSect));
		System.out.println("GV 08 to O dist: "+min08Dist);
		
		double min10Dist = Double.POSITIVE_INFINITY;
		for (FaultSection gvSect : gv10Sects)
			for (FaultSection oSect : ortigalitaSects)
				min10Dist = Double.min(min10Dist, minDist(gvSect, oSect));
		System.out.println("GV 10 to O dist: "+min10Dist);
		
		// reset depths
		double ave07Upper = 0d;
		for (FaultSection gvSect : gv07Sects)
			ave07Upper += gvSect.getOrigAveUpperDepth();
		ave07Upper /= gv07Sects.size();
		double ave08Upper = 0d;
		for (FaultSection gvSect : gv08Sects)
			ave08Upper += gvSect.getOrigAveUpperDepth();
		ave08Upper /= gv08Sects.size();
		System.out.println("GV 07 Upper: "+ave07Upper);
		System.out.println("GV 08 Upper: "+ave08Upper);
		
		// mod gv08 depths
		double delta = ave07Upper - ave08Upper;
//		delta *= 3;
		for (FaultSection sect : gv08Sects) {
			((FaultSectionPrefData)sect).setAveLowerDepth(sect.getAveLowerDepth()+delta);
			((FaultSectionPrefData)sect).setAveUpperDepth(sect.getOrigAveUpperDepth()+delta);
			RuptureSurface surf = sect.getFaultSurface(1d+0.05*Math.random()); // to clear cache
			System.out.println("Surf upper: "+surf.getAveRupTopDepth());
//			System.out.println("Surf lower: "+surf.getave);
		}
		
		double mod08Dist = Double.POSITIVE_INFINITY;
		for (FaultSection gvSect : gv08Sects)
			for (FaultSection oSect : ortigalitaSects)
				mod08Dist = Double.min(mod08Dist, minDist(gvSect, oSect));
		System.out.println("Mod 08 to O dist: "+mod08Dist);
	}

	
	private static double minDist(FaultSection s1, FaultSection s2) {
		RuptureSurface surf1 = s1.getFaultSurface(1d, false, false);
		RuptureSurface surf2 = s2.getFaultSurface(1d, false, false);
		
		double minDist = Double.POSITIVE_INFINITY;
		for (Location l1 : surf1.getEvenlyDiscritizedListOfLocsOnSurface())
			for (Location l2 : surf2.getEvenlyDiscritizedListOfLocsOnSurface())
				minDist = Math.min(minDist, LocationUtils.linearDistanceFast(l1, l2));
		
		return minDist;
	}
	
	private static void test62() throws IOException {
		RSQSimCatalog catalog = Catalogs.BRUCE_2585_1MYR.instance();
		RSQSimSubSectionMapper mapper = catalog.getSubSectMapper();
		HashSet<String> names = new HashSet<>();
		for (FaultSection sect : mapper.getSubSections()) {
			if (sect.getName().contains("Mojave")) {
				Collection<SimulatorElement> elems = mapper.getElementsForSection(sect);
				names.add(elems.iterator().next().getSectionName());
			}
		}
		for (String name : names)
			System.out.println(name);
	}
	
	private static void test63() {
		int num = (int)Math.round((9d)/0.1d)+1;
		EvenlyDiscretizedFunc logXVals = new EvenlyDiscretizedFunc(0d, num, 0.1);
		System.out.println(logXVals);
		DiscretizedFunc xVals = new ArbitrarilyDiscretizedFunc();
		xVals.set(0d, 0d);
		for (Point2D pt : logXVals)
			xVals.set(Math.pow(10, pt.getX()), 0d);
		System.out.println(xVals);
	}
	
	private static void test64() throws IOException {
		RSQSimCatalog catalog = Catalogs.BRUCE_4965.instance();
		double firstTime = -1;
		for (RSQSimEvent event : catalog.loader().iterable()) {
			double time = event.getTimeInYears();
			if (firstTime < 0) {
				System.out.println("First event is at "+time+"yr");
				firstTime = time;
			}
			if (event.getMagnitude() > 8.5)
				System.out.println("Event "+event.getID()+" is M"+(float)event.getMagnitude()+" at "+time+"yr");
		}
	}
	
	private static void test65() {
		float start = 1000f;
		float next = Math.nextUp(start);
		float maxDiff = 0;
		for (int i=0; i<1000; i++) {
			float diff = next - start;
			maxDiff = Float.max(maxDiff, diff);
		}
		System.out.println("diff: "+maxDiff);
	}
	
	private static void test66() throws UnknownHostException {
		InetAddress host = InetAddress.getLocalHost();
		System.out.println("My hostname (getCanonicalHostName): "+host.getCanonicalHostName());
		System.out.println("My hostname (getHostName): "+host.getHostName());
		System.out.println("My hostname (getHostAddress): "+host.getHostAddress());
		System.out.println("Servlet URL: "+ServerPrefUtils.SERVER_PREFS.getServletBaseURL());
		System.out.println("Testing servlet:");
		System.out.flush();
		try {
			new WillsMap2015().getValue(new Location(34, -118));
			System.out.println("Success!");
		} catch (Exception e) {
			System.err.println(e.getMessage());
			System.err.flush();
			System.out.println("Fail :(");
		}
		System.out.println("Testing http:");
		System.out.flush();
		try {
			FileUtils.downloadURL("http://opensha.usc.edu/ftp", File.createTempFile("test", ".tmp"));
			System.out.println("Success!");
		} catch (Exception e) {
			System.err.println(e.getMessage());
			System.err.flush();
			System.out.println("Fail :(");
		}
		System.out.println("Testing http IP:");
		System.out.flush();
		try {
			FileUtils.downloadURL("http://68.181.32.140/ftp", File.createTempFile("test", ".tmp"));
			System.out.println("Success!");
		} catch (Exception e) {
			System.err.println(e.getMessage());
			System.err.flush();
			System.out.println("Fail :(");
		}
		System.out.println("Testing localhost:");
		System.out.flush();
		try {
			FileUtils.downloadURL("http://localhost/ftp", File.createTempFile("test", ".tmp"));
			System.out.println("Success!");
		} catch (Exception e) {
			System.err.println(e.getMessage());
			System.err.flush();
			System.out.println("Fail :(");
		}
		System.out.println("Testing 127.0.0.1:");
		System.out.flush();
		try {
			FileUtils.downloadURL("http://127.0.0.1/ftp", File.createTempFile("test", ".tmp"));
			System.out.println("Success!");
		} catch (Exception e) {
			System.err.println(e.getMessage());
			System.err.flush();
			System.out.println("Fail :(");
		}
	}
	
	private static void test67() {
		System.out.println(new CaliforniaRegions.RELM_TESTING().contains(new Location(38.1589, -117.8749)));
	}
	
	private static void test68() throws IOException {
		List<FaultSection> sects31 = FaultModels.FM3_1.fetchFaultSections();
		List<FaultSection> sects32 = FaultModels.FM3_1.fetchFaultSections();
		FaultSectionDataWriter.writeSectionsToFile(sects31, null, new File("/tmp/fm_3_1.txt"), false);
		FaultSectionDataWriter.writeSectionsToFile(sects32, null, new File("/tmp/fm_3_2.txt"), false);
	}
	
	private static void test69() throws ZipException, IOException {
		File zipFile = new File("/home/kevin/OpenSHA/UCERF3/fss_csvs/full_model_csvs.zip");
		ZipFile zip = new ZipFile(zipFile);
		
		Enumeration<? extends ZipEntry> entriesEnum = zip.entries();
		
		while (entriesEnum.hasMoreElements()) {
			ZipEntry entry = entriesEnum.nextElement();
			if (!entry.getName().contains("FM3_1"))
				continue;
			CSVFile<String> csv = CSVFile.readStream(zip.getInputStream(entry), true);
			if (csv.getNumRows() < 253707)
				System.out.println("Bad CSV? "+entry.getName()+" has "+csv.getNumRows());
		}
		zip.close();
	}
	
	private static void test70() throws ZipException, IOException {
		ArrayList<FaultSection> sects = FaultModels.FM3_1.fetchFaultSections();
		for (FaultSection sect : sects)
			System.out.println(sect.getSectionId()+". "+sect.getSectionName());
		System.out.println(sects.size()+" sects");
	}
	
	private static void test71() throws ZipException, IOException {
		ObsEqkRupList loadedRups = UCERF3_CatalogParser.loadCatalog(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/EarthquakeCatalog/"
						+ "ofr2013-1165_EarthquakeCat.txt"));
		for (ObsEqkRupture rup : loadedRups) {
			GregorianCalendar cal = rup.getOriginTimeCal();
			int year = cal.get(GregorianCalendar.YEAR);
			int month = cal.get(GregorianCalendar.MONTH)+1;
			int day = cal.get(GregorianCalendar.DAY_OF_MONTH);
			if (rup.getMag() > 7d)
				System.out.println("M"+(float)rup.getMag()+" on "+year+"/"+month+"/"+day
						+" at "+rup.getHypocenterLocation());
		}
	}
	
	private static void test72() throws IOException, DocumentException {
		Map<FaultModels, FaultSystemRupSet> rupSetMap = new HashMap<>();
		rupSetMap.put(FaultModels.FM3_1, U3FaultSystemIO.loadRupSet(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip")));
		rupSetMap.put(FaultModels.FM3_2, U3FaultSystemIO.loadRupSet(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_2_MEAN_BRANCH_AVG_SOL.zip")));
		
		for (FaultModels fm : rupSetMap.keySet()) {
			FaultSystemRupSet origRupSet = rupSetMap.get(fm);
			System.out.println(fm);
			FaultSystemRupSet newRupSet = InversionFaultSystemRupSetFactory.forBranch(fm);
			System.out.println("new has "+newRupSet.getNumRuptures()+" ruptures");
			System.out.println("orig has "+origRupSet.getNumRuptures()+" ruptures");
			if (newRupSet.getNumRuptures() != origRupSet.getNumRuptures())
				continue;
			int numCountDiffs = 0;
			int numContentsDiffs = 0;
			for (int r=0; r<origRupSet.getNumRuptures(); r++) {
				List<Integer> sects1 = origRupSet.getSectionsIndicesForRup(r);
				List<Integer> sects2 = newRupSet.getSectionsIndicesForRup(r);
				if (sects1.size() != sects2.size()) {
					numCountDiffs++;
					continue;
				}
				for (int s=0; s<sects1.size(); s++) {
					if (!sects1.get(s).equals(sects2.get(s))) {
						numContentsDiffs++;
						break;
					}
				}
				if (r == 0 || r == 10) {
					System.out.println("Orig rup "+r+": "+Joiner.on(",").join(sects1));
					System.out.println("New rup "+r+": "+Joiner.on(",").join(sects2));
				}
			}
			System.out.println("Differences: "+numCountDiffs+" in count");
			System.out.println("Differences: "+numContentsDiffs+" in contents only");
		}
	}
	
	private static void test73() {
		FaultSection sect = FaultModels.FM3_1.fetchFaultSectionsMap().get(541);
		System.out.println(sect.getSectionName());
		RuptureSurface surf = sect.getFaultSurface(1d, false, false);
		for (BBP_Site site : RSQSimBBP_Config.getCyberShakeVs500LASites()) {
			Location loc = site.getLoc();
			double dist = surf.getDistanceJB(loc);
			System.out.println("Site "+site.getName()+" is "+(float)dist+" km away");
		}
	}
	
	private static void test74() {
		ArbitrarilyDiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
		func.set(1e-1, 0.1d);
		func.set(2e-1, 0.2d);
		func.set(3e-1, 0.3d);
		func.set(4e-1, 0.4d);
		func.set(5e-1, 0.5d);
		
		GraphWindow gw = new GraphWindow(func, "Test",
				new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));
		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
		gw.setXLog(true);
	}
	
	private static void test75() {
		int subSize = 	60000;
		int fullSize = 	100000;
		double n = ((double)subSize*fullSize)/((double)subSize + fullSize);
		System.out.println("\tN="+(float)n+" (leaf count: "+subSize+")");
//		System.out.println("\tDn: "+(float)minDn);
		double dnThresh = 1.63/Math.sqrt(n);
		System.out.println("\tThreshold: "+dnThresh);
//		boolean passes = minDn <= 1.63/Math.sqrt(n);
		
		System.out.println("now the other way around");
		n = ((double)subSize + fullSize)/((double)subSize*fullSize);
		System.out.println("\tN="+(float)n+" (leaf count: "+subSize+")");
//		System.out.println("\tDn: "+(float)minDn);
		dnThresh = 1.63*Math.sqrt(n);
		System.out.println("\tThreshold: "+dnThresh);
	}
	
	private static void test76() throws ZipException, IOException, DocumentException {
		File rupSetsDir = new File("/home/kevin/OpenSHA/UCERF4/rup_sets");
		File rupSetFile = new File(rupSetsDir, "fm3_1_ucerf3.zip");
		
		try {
			Thread.sleep(5000);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		FaultSystemRupSet rupSet = U3FaultSystemIO.loadRupSet(rupSetFile);
		SectionDistanceAzimuthCalculator fakeDistAzCalc =
				new SectionDistanceAzimuthCalculator(rupSet.getFaultSectionDataList()) {
			public double getDistance(FaultSection sect1, FaultSection sect2) {
				return 0d;
			}
			
			public double getDistance(int id1, int id2) {
				return 0d;
			}
		};
		System.gc();
		System.out.println("Loaded rupture set");
		try {
			Thread.sleep(20000);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.out.println("Loading cluster ruptures");
		List<ClusterRupture> rups = new ArrayList<>();
//		for (int r=0; r<rupSet.getNumRuptures(); r++)
		for (int r=0; r<20000; r++)
			rups.add(ClusterRupture.forOrderedSingleStrandRupture(
					rupSet.getFaultSectionDataForRupture(r), fakeDistAzCalc));
		System.out.println("DONE");
		System.gc();
		try {
			Thread.sleep(100000000);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private static void test77() throws ZipException, IOException, DocumentException {
		U3FaultSystemIO.loadRupSet(new File("/home/kevin/OpenSHA/UCERF4/rup_sets/fm3_1_reproduce_ucerf3.zip"));
	}
	
//	public void process(Collection<? extends ObsEqkRupture> catalog)
	
//	public <E extends ObsEqkRupture> List<E> process(Collection<E> catalog) {
//		List<E> ret = new ArrayList<>();
//		// do stuff
//		return ret;
//	}
//	
//	public <E>
	
	private static void test78() throws IOException, DocumentException {
		File fssFile = new File("/home/kevin/git/ucerf3-etas-launcher/inputs/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_SpatSeisU3_MEAN_BRANCH_AVG_SOL.zip");
		File resultsFile = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2019_11_05-Start2012_500yr_kCOV1p5_Spontaneous_HistoricalCatalog/"
				+ "results_m5_preserve_chain.bin");
		AbstractGridSourceProvider.SOURCE_MIN_MAG_CUTOFF = 2.55d;
		FaultSystemSolution sol = U3FaultSystemIO.loadSol(fssFile);
		FaultSystemSolutionERF_ETAS erf = ETAS_Launcher.buildERF(sol,
				false, 500d, 2012);
		erf.updateForecast();

		// 5d will filter out below M5
		List<ETAS_Catalog> catalogs =
				ETAS_CatalogIO.loadCatalogsBinary(resultsFile, 5d);

		int numPtMatches = 0;
		int numPtFails = 0;
		for (ETAS_Catalog catalog : catalogs) {
			for (ETAS_EqkRupture rup : catalog) {
				int nth = rup.getNthERF_Index();
//				EqkRupture nthRup = erf.getNthRupture(nth);
				int sourceID = erf.getSrcIndexForNthRup(nth);
				ProbEqkSource source = erf.getSource(sourceID);
				int rupID = erf.getRupIndexInSourceForNthRup(nth);
				EqkRupture nthRup = source.getRupture(rupID);
				if ((float)rup.getMag() != (float)nthRup.getMag()) {
					Preconditions.checkState(source instanceof PointSource13b);
					if (numPtFails == 0) {
						System.out.println("ETAS rup w/ M="+(float)rup.getMag()+", loc="+rup.getHypocenterLocation());
						System.out.println("\tsource "+sourceID+" for nth="+nth+" is type "+source.getClass());
						System.out.println("\trup "+rupID+" for nth has M="+(float)nthRup.getMag()
							+", loc="+nthRup.getHypocenterLocation());
						System.out.println("MAG MISMATCH. Nth rup has M="+(float)nthRup.getMag());
						PointSource13b ptSrc = (PointSource13b)source;
						for (int r=0; r<ptSrc.getNumRuptures(); r++) {
							ProbEqkRupture ptRup = ptSrc.getRupture(r);
							System.out.println("\tr="+r+",\tM"+(float)ptRup.getMag()
								+",\trake="+(float)ptRup.getAveRake()+",\tp="+(float)ptRup.getProbability());
						}
					}
					numPtFails++;
				} else {
					if (source instanceof PointSource13b)
						numPtMatches++;
					else
						Preconditions.checkState((float)rup.getMag() ==
							(float)nthRup.getMag(),
							"Nth rupture mag=%s doesn't match ETAS mag=%s",
							nthRup.getMag(), rup.getMag());
				}
//				Preconditions.checkState((float)rup.getMag() ==
//						(float)nthRup.getMag(),
//						"Nth rupture mag=%s doesn't match ETAS mag=%s",
//						nthRup.getMag(), rup.getMag());
				// you can then either use nthRup, or set the surface in the original rupture as here:
				rup.setRuptureSurface(nthRup.getRuptureSurface());
//				new PointSource13
			}
		}
		System.out.println(numPtMatches+" pt matches and "+numPtFails
				+" fails, "+(numPtMatches+numPtFails)+" tot pt source");
	}
	
	private static void test79() {
//		SparkSession spark = SparkSession.builder().master("local").getOrCreate();
//		double[] data = { 0, 1, 2, 3, 4 };
//		Dataset<Double> dataset = spark.createDataset(Doubles.asList(data), Encoders.DOUBLE());
//		dataset.toJavaRDD().map(s -> Vectors.dense(s));
//		
//		List<FaultSection> sects = null;
//		Map<Integer, List<FaultSection>> parentMap = sects.stream().collect(Collectors.groupingBy(s -> s.getParentSectionId()));
//		sects.removeIf(s -> s.getParentSectionId() == 2);
//		GaussianMixture.fit(data);
	}
	
	private static void test80() throws IOException, GMT_MapException {
		GriddedRegion reg = new GriddedRegion(new CaliforniaRegions.RELM_TESTING(), 0.02, GriddedRegion.ANCHOR_0_0);
		GriddedGeoDataSet xyz = new GriddedGeoDataSet(reg, false);
		for (String line : Files.readLines(new File("/home/kevin/Downloads/map_data.txt"), Charset.defaultCharset())) {
			line = line.trim();
			if (line.isEmpty())
				continue;
			String[] split = line.split("\t");
			Preconditions.checkState(split.length == 3);
			double lat = Double.parseDouble(split[0]);
			double lon = Double.parseDouble(split[1]);
			double val = Double.parseDouble(split[2]);
			xyz.set(new Location(lat, lon), val);
		}
		int lineIndex = 0;
		double prevLat = 0d;
		for (int i=0; i<xyz.size(); i++) {
			Location loc = xyz.getLocation(i);
			if (loc.getLatitude() != prevLat) {
				prevLat = loc.getLatitude();
				lineIndex++;
			}
			if (lineIndex % 2 != 0)
				xyz.set(i, -6d);
			else
				xyz.set(i, -(i % 4));
		}
		System.out.println(lineIndex+" lines");
		Region plotReg = reg;
//		Region plotReg = new Region(new Location(34, -118), new Location(35, -120));
		
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(-6, -2);
		GMT_Map map = new GMT_Map(plotReg, xyz, reg.getSpacing(), cpt);
		map.setUseGMTSmoothing(false);
		map.setCustomLabel("Test map");
//		map.setDpi(600);
//		map.setDpi(300);
//		map.setImageWidth(20d);
		map.setUseGRDView(true);
		
		FaultBasedMapGen.plotMap(new File("/tmp"), "test_map", false, map);
	}

	private static void test81() {
		// Get current size of heap in bytes
		double heapSize = Runtime.getRuntime().totalMemory(); 
		heapSize /= 1024; // KB
		heapSize /= 1024; // MB
		System.out.println("Current heap: "+(float)heapSize+" MB = "+(float)(heapSize/1024d)+" GB");

		// Get maximum size of heap in bytes. The heap cannot grow beyond this size.// Any attempt will result in an OutOfMemoryException.
		double heapMaxSize = Runtime.getRuntime().maxMemory();
		heapMaxSize /= 1024; // KB
		heapMaxSize /= 1024; // MB
		System.out.println("Max heap: "+(float)heapMaxSize+" MB = "+(float)(heapMaxSize/1024d)+" GB");
	}

	private static void test82() {
		String str = "|asdf|";
//		System.out.println(str.replaceAll("|", Matcher.quoteReplacement("\\|")));
		System.out.println(str.replace("|", "\\|"));
		System.exit(0);
		
		// list of faults with 90 degree dip and non-SS rakes
		FaultModels fm = FaultModels.FM3_1;
		List<FaultSection> sects = fm.fetchFaultSections();
		Map<Integer, FaultSection> parentIDsMap = sects.stream().collect(Collectors.toMap(S -> S.getSectionId(), S->S));
		
		Table<DeformationModels, Integer, Double> dmParentRakes = HashBasedTable.create();
		for (DeformationModels dm : DeformationModels.values()) {
			if (dm.getRelativeWeight(null) == 0d)
				continue;
			List<? extends FaultSection> subSects = DeformationModels.loadSubSects(fm, dm);
			Map<Integer, List<FaultSection>> parentSects = subSects.stream().collect(Collectors.groupingBy(S -> S.getParentSectionId()));
			for (Integer parentID : parentIDsMap.keySet()) {
				List<Double> rakes = new ArrayList<>();
				for (FaultSection sect : parentSects.get(parentID))
					rakes.add(sect.getAveRake());
				double aveRake = FaultUtils.getInRakeRange(FaultUtils.getAngleAverage(rakes));
//				if (parentID == 301) {
//					System.out.println("301. rakes="+Joiner.on(",").join(rakes)+". ave="+aveRake);
//					System.exit(0);
//				}
				dmParentRakes.put(dm, parentID, aveRake);
			}
		}
		for (FaultSection sect : sects) {
			if (sect.getAveDip() != 90d)
				continue;
			Map<DeformationModels, Double> dmRakes = dmParentRakes.column(sect.getSectionId());
			double geolRake = dmRakes.get(DeformationModels.GEOLOGIC);
			if ((float)geolRake != -180f && (float)geolRake != 0f && (float)geolRake != 180f) {
				System.out.println(sect.getName()+" (ID="+sect.getSectionId()+") has dip=90 and non-SS geologic rake="+(float)geolRake);
				for (DeformationModels dm : DeformationModels.values()) {
					if (dmRakes.containsKey(dm))
						System.out.println("\t"+dm.getName()+": "+dmRakes.get(dm).floatValue());
				}
			}
		}
	}
	
	private static void test83() throws IOException {
		List<? extends FaultSection> subSects = DeformationModels.loadSubSects(FaultModels.FM3_1, DeformationModels.GEOLOGIC);
		FaultSubsectionCluster cluster = new FaultSubsectionCluster(subSects.subList(0, 3));
		FaultSubsectionCluster splayCluster = new FaultSubsectionCluster(subSects.subList(30, 31));
		Jump jump = new Jump(cluster.subSects.get(1), cluster, splayCluster.subSects.get(0), splayCluster, 1d);
		cluster.addConnection(jump);
		ClusterRupture rup = new ClusterRupture(cluster);
		rup = rup.take(jump);
		System.out.println(rup);
		File file = new File("/tmp/splay.json");
		ClusterRupture.writeJSON(file, Lists.newArrayList(rup), subSects);
		List<ClusterRupture> rups = ClusterRupture.readJSON(file, subSects);
		System.out.println(rups.get(0));
	}
	
	private static void test84() throws IOException {
		// corupture rate calculation for Mike Oskin & Alba Padilla
		HashSet<String> sects1 = new HashSet<>();
		HashSet<String> sects2 = new HashSet<>();
		
		String name1 = "SAF San Bernardino";
		sects1.add("San Andreas (San Bernardino N), Subsection 0");
		sects1.add("San Andreas (San Bernardino N), Subsection 1");
		sects1.add("San Andreas (San Bernardino N), Subsection 2");
		sects1.add("San Andreas (San Bernardino N), Subsection 3");
		sects1.add("San Andreas (San Bernardino N), Subsection 4");
		
		String name2 = "SJF San Bernardino";
		sects2.add("San Jacinto (San Bernardino), Subsection 0");
		sects2.add("San Jacinto (San Bernardino), Subsection 1");
		sects2.add("San Jacinto (San Bernardino), Subsection 2");
		
//		String name1 = "SAF Wrightwood";
//		sects1.add("San Andreas (Mojave S), Subsection 13");
//		
//		String name2 = "SJF Mystic Lake";
//		sects2.add("San Jacinto (Stepovers Combined), Subsection 0");
		
		CompoundFaultSystemSolution cfss = CompoundFaultSystemSolution.fromZipFile(
				new File("/home/kevin/OpenSHA/UCERF3/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip"));
		APrioriBranchWeightProvider weightProv = new APrioriBranchWeightProvider();
		double totWeight = 0d;
		Map<FaultModels, List<Integer>> fmCorupsMap = new HashMap<>();
		Map<FaultModels, List<Integer>> fmRups1Map = new HashMap<>();
		Map<FaultModels, List<Integer>> fmRups2Map = new HashMap<>();
		List<Double> corupVals = new ArrayList<>();
		List<Double> vals1 = new ArrayList<>();
		List<Double> vals2 = new ArrayList<>();
		List<Double> weights = new ArrayList<>();
		
		for (U3LogicTreeBranch branch : cfss.getBranches()) {
			double weight = weightProv.getWeight(branch);
			weights.add(weight);
			totWeight += weight;
			FaultModels fm = branch.getValue(FaultModels.class);
			List<Integer> corupIDs = fmCorupsMap.get(fm);
			if (corupIDs == null) {
				FaultSystemRupSet rupSet = cfss.getSolution(branch).getRupSet();
				corupIDs = new ArrayList<>();
				List<Integer> ids1 = new ArrayList<>();
				List<Integer> ids2 = new ArrayList<>();
				for (int r=0; r<rupSet.getNumRuptures(); r++) {
					boolean has1 = false;
					boolean has2 = false;
					for (FaultSection sect : rupSet.getFaultSectionDataForRupture(r)) {
						has1 = has1 || sects1.contains(sect.getSectionName());
						has2 = has2 || sects2.contains(sect.getSectionName());
					}
					if (has1 && has2)
						corupIDs.add(r);
					if (has1)
						ids1.add(r);
					if (has2)
						ids2.add(r);
				}
				System.out.println("Rup counts for "+fm);
				System.out.println("\tFault 1: "+ids1.size());
				System.out.println("\tFault 2: "+ids2.size());
				System.out.println("\tCoruptures: "+corupIDs.size());
				fmCorupsMap.put(fm, corupIDs);
				fmRups1Map.put(fm, ids1);
				fmRups2Map.put(fm, ids2);
			}
			List<Integer> ids1 = fmRups1Map.get(fm);
			List<Integer> ids2 = fmRups2Map.get(fm);
			double[] rates = cfss.getRates(branch);
			corupVals.add(calcRate(rates, corupIDs));
			vals1.add(calcRate(rates, ids1));
			vals2.add(calcRate(rates, ids2));
		}
		
		System.out.println("Co-rupture rate:");
		printRateStats(corupVals, weights, totWeight);
		System.out.println(name1+" rates");
		printRateStats(vals1, weights, totWeight);
		System.out.println(name2+" rates");
		printRateStats(vals2, weights, totWeight);
	}
	
	private static void printRateStats(List<Double> rates, List<Double> weights, double totWeight) {
		double mean = 0d;
		double min = Double.POSITIVE_INFINITY;
		double max = Double.NEGATIVE_INFINITY;
		for (int i=0; i<rates.size(); i++) {
			double val = rates.get(i);
			mean += val*weights.get(i);
			min = Math.min(min, val);
			max = Math.max(max, val);
		}
		mean /= totWeight;
		
		System.out.println("\tBranch-averaged rate: "+(float)mean);
		System.out.println("\tBranch rate range: ["+(float)min+", "+(float)max+"]");
	}
	
	private static double calcRate(double[] rates, List<Integer> ids) {
		double rate = 0d;
		for (int id : ids)
			rate += rates[id];
		return rate;
	}
	
	private static void test85() {
		for (FaultSection sect : FaultModels.FM3_2.fetchFaultSections()) {
			float dipDir = sect.getDipDirection();
			float calcDipDir = (float)(sect.getFaultTrace().getStrikeDirection()+90d);
			if (calcDipDir >= 360f)
				calcDipDir -=360f;
			if (Math.abs(calcDipDir-dipDir) > 2d)
				System.out.println(sect.getSectionName()+"\tdipDir="+dipDir+"\tcalcDir="+calcDipDir);
		}
	}
	
	private static void test86() throws ZipException, IOException, DocumentException {
//		// do we have landers?
//		FaultSystemRupSet u3 = FaultSystemIO.loadRupSet(new File("/home/kevin/OpenSHA/UCERF4/rup_sets/fm3_1_ucerf3.zip"));
//		FaultSystemRupSet candidate = FaultSystemIO.loadRupSet(new File("/home/kevin/OpenSHA/UCERF4/rup_sets/"
//				+ "fm3_1_plausible10km_direct_slipP0.05incr_cff0.75IntsPos_comb2Paths_cffFavP0.02_cffFavRatioN2P0.5_sectFractPerm0.05.zip"));
//		int landersID = 246711;
//		ClusterRupture landersU3 = ClusterRupture.forOrderedSingleStrandRupture(
//				u3.getFaultSectionDataForRupture(landersID), candidate.getPlausibilityConfiguration().getDistAzCalc());
//		System.out.println("UCERF3 Landers: "+landersU3);
//		for (ClusterRupture rup : candidate.getClusterRuptures()) {
//			if (landersU3.unique.equals(rup.unique)) {
//				System.out.println("Candidate has landers!");
//				System.out.println("\t"+rup);
//				break;
//			}
//		}
	}
	
	private static void test87() {
		Map<Integer, FaultSection> map31 = FaultModels.FM3_1.fetchFaultSectionsMap();
		Map<Integer, FaultSection> map32 = FaultModels.FM3_2.fetchFaultSectionsMap();
		List<String> uniques1 = new ArrayList<>();
		for (Integer id : map31.keySet())
			if (!map32.containsKey(id))
				uniques1.add(map31.get(id).getName());
		Collections.sort(uniques1);
		for (String name : uniques1)
			System.out.println(name);
	}
	
	private static void test88() throws IOException {
		ETAS_Config conf = ETAS_Config.readJSON(new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2021_03_02-Start2012_500yr_CompletenessSTRICT_Spontaneous_HistoricalCatalog/config.json"));
		ETAS_Launcher launcher = new ETAS_Launcher(conf, false);
		AbstractERF erf = launcher.checkOutERF();
		System.out.println("ERF params");
		for (Parameter<?> param : erf.getAdjustableParameterList())
			System.out.println(param.getName()+":\t"+param.getValue());
		System.out.println("ETAS params");
		for (Parameter<?> param : launcher.getETAS_Params())
			System.out.println(param.getName()+":\t"+param.getValue());
		List<ETAS_EqkRupture> hist = launcher.getHistQkList();
		System.out.println("Have "+hist.size()+" hist quakes");
		ETAS_EqkRupture last = hist.get(hist.size()-1);
		System.out.println("\tLast: "+last.getEventId()+": M"+last.getMag());
	}
	
	private static void test89() throws IOException {
		File binFile = new File("/data/kevin/ucerf3/etas/simulations/"
//				+ "2016_02_17-spontaneous-1000yr-scaleMFD1p14-full_td-subSeisSupraNucl-gridSeisCorr/results_m5_preserve.bin");
				+ "2021_03_08-spontaneous-500yr-full_td-subSeisSupraNucl-gridSeisCorr-scale1.14-REDO-2016-JAR/results_m5.bin");
		double duration = 1000d;
		double modDuration = 500d;
		List<ETAS_Catalog> catalogs = ETAS_CatalogIO.loadCatalogsBinary(binFile);
		int[] maxIdenticalEvents = { -1 };
		if (modDuration != duration) {
			long firstEvent = Long.MAX_VALUE;
			for (ETAS_Catalog catalog : catalogs)
				firstEvent = Long.min(firstEvent, catalog.get(0).getOriginTime());
			long maxTime = (long)(firstEvent + ProbabilityModelsCalc.MILLISEC_PER_YEAR*modDuration);
			List<ETAS_Catalog> newCats = new ArrayList<>();
			for (ETAS_Catalog cat : catalogs) {
				ETAS_Catalog newCat = new ETAS_Catalog(null);
				for (ETAS_EqkRupture rup : cat) {
					if (rup.getOriginTime() > maxTime)
						break;
					newCat.add(rup);
				}
				newCats.add(newCat);
			}
			System.out.println("Truncated to "+modDuration+" years");
			catalogs = newCats;
			duration = modDuration;
		}
		
		for (int maxIdentical : maxIdenticalEvents) {
			List<HashSet<Integer>> grouped = new ArrayList<>();
			int numMatches = 0;
			System.out.println("Matching with "+maxIdentical+" events...");
			for (int i=0; i<catalogs.size(); i++) {
				ETAS_Catalog cat1 = catalogs.get(i);
				for (int j=i+1; j<catalogs.size(); j++) {
					ETAS_Catalog cat2 = catalogs.get(j);
					boolean match = true;
					int maxLen = Integer.max(cat1.size(), cat2.size());
					for (int n=0; match && n<maxLen && (maxIdentical <= 0 || n < maxIdentical); n++) {
						if (n == cat1.size() || n == cat2.size()) {
							match = false;
							break;
						}
						ETAS_EqkRupture rup1 = cat1.get(n);
						ETAS_EqkRupture rup2 = cat2.get(n);
						match = rup1.getGridNodeIndex() == rup2.getGridNodeIndex()
								&& rup1.getFSSIndex() == rup2.getFSSIndex()
								&& (float)rup1.getMag() == (float)rup2.getMag();
					}
					if (match) {
//						System.out.println("\t"+i+" matches "+j);
						boolean found = false;
						HashSet<Integer> myGroup = null;
						for (HashSet<Integer> group : grouped) {
							if (group.contains(i)) {
								myGroup = group;
								break;
							}
						}
						if (myGroup == null) {
							myGroup = new HashSet<>();
							grouped.add(myGroup);
						}
						myGroup.add(i);
						myGroup.add(j);
						numMatches++;
					}
				}
			}
			System.out.println("Found "+grouped.size()+" groups:");
			int testNum = 0;
			int numNonUnique = 0;
			HashSet<Integer> nonUniqueIndexes = new HashSet<>();
			for (HashSet<Integer> group : grouped) {
//				testNum += group.size();
				for (int i=0; i<group.size(); i++)
					for (int j=i+1; j<group.size(); j++)
						testNum++;
				System.out.println("\t"+Joiner.on(",").join(group));
				numNonUnique += group.size();
				nonUniqueIndexes.addAll(group);
			}
			Preconditions.checkState(numMatches == testNum, "%s != %s", numMatches, testNum);
			System.out.println("\t"+(catalogs.size()-numNonUnique)+"/"+catalogs.size()+" are unique");
			System.out.println("Total M5 rate: "+annualRate(catalogs, duration, 5d));
			List<ETAS_Catalog> uniqueCatalogs = new ArrayList<>();
			for (int i=0; i<catalogs.size(); i++)
				if (!nonUniqueIndexes.contains(i))
					uniqueCatalogs.add(catalogs.get(i));
			System.out.println("Unique catalog M5 rate: "+annualRate(uniqueCatalogs, duration, 5d));
		}
	}
	
	private static double annualRate(Collection<ETAS_Catalog> catalogs, double duration, double minMag) {
		long count = 0;
		for (ETAS_Catalog catalog : catalogs)
			for (ETAS_EqkRupture rup : catalog)
				if ((float)rup.getMag() >= 5f)
					count++;
		return (double)count/(duration*catalogs.size());
	}
	
	private static void test90() throws ZipException, IOException, DocumentException {
//		File rupSetsDir = new File("/home/kevin/OpenSHA/UCERF4/rup_sets/");
//		File mainFile = new File(rupSetsDir, "fm3_1_plausibleMulti10km_direct_slipP0.05incr_cff0.75IntsPos_comb2Paths_"
//				+ "cffFavP0.02_cffFavRatioN2P0.5_sectFractPerm0.05.zip");
//		File altFile = new File(rupSetsDir, "fm3_1_plausibleMulti10km_direct_slipP0.05incr_cff0.75IntsPos_comb2Paths_"
//				+ "cffFavP0.02_cffFavRatioN2P0.5_sectFractPerm0.05_comp/alt_perm_Bilateral_Adaptive_5SectIncrease_MaintainConnectivity.zip");
//		
//		FaultSystemRupSet mainRupSet = FaultSystemIO.loadRupSet(mainFile);
//		FaultSystemRupSet altRupSet = FaultSystemIO.loadRupSet(altFile);
//		
//		ClusterRupture mainRup = mainRupSet.getClusterRuptures().get(137320);
//		ClusterRupture altRup = altRupSet.getClusterRuptures().get(226038);
//		
//		System.out.println("Main: "+mainRup);
//		System.out.println("\thash: "+mainRup.unique.hashCode());
//		System.out.println("Alt: "+altRup);
//		System.out.println("\thash: "+altRup.unique.hashCode());
//		System.out.println("Unique equals? "+mainRup.unique.equals(altRup.unique));
//		System.out.println("Regular equals? "+mainRup.equals(altRup));
	}
	
	private static void test91() throws IOException, DocumentException {
		FaultSystemSolution sol = U3FaultSystemIO.loadSol(
				new File("/home/kevin/OpenSHA/UCERF4/rup_sets/fm3_1_ucerf3.zip"));
		
		double rsRate = 1d/714516d;
		int count = 0;
		for (double rate : sol.getRateForAllRups())
			if (rate <rsRate)
				count++;
		FaultSystemRupSet rupSet = sol.getRupSet();
		System.out.println(count+"/"+rupSet.getNumRuptures()+" ruptures have rates below "+rsRate+" ("
				+new DecimalFormat("0.00%").format((double)count/(double)rupSet.getNumRuptures())+")");
		
		FaultSystemRupSet rsRupSet = U3FaultSystemIO.loadRupSet(new File("/home/kevin/OpenSHA/UCERF4/rup_sets/"
				+ "rsqsim_4983_stitched_m6.5_skip65000_sectArea0.5.zip"));
		HashSet<UniqueRupture> uniques = new HashSet<>();
		for (int r=0; r<rsRupSet.getNumRuptures(); r++)
			uniques.add(UniqueRupture.forIDs(rsRupSet.getSectionsIndicesForRup(r)));
		System.out.println("RSQSim has "+uniques.size()+"/"+rsRupSet.getNumRuptures()+" unique ruptures ("
				+new DecimalFormat("0.00%").format((double)uniques.size()/(double)rsRupSet.getNumRuptures())+")");
	}
	
	private static void test92() {
		int n = 1001;
		
		int cnt = 0;
		for (int i=0; i<n; i++)
			for (int j=i+1; j<n; j++)
				cnt++;
		System.out.println(cnt);
		System.out.println(n*(n-1)/2);
	}
	
	private static void test93() throws ZipException, IOException {
		Thread printingHook = new Thread(() -> System.out.println("In the middle of a shutdown"));
		Runtime.getRuntime().addShutdownHook(printingHook);
		JFrame frame = new JFrame();
		frame.setSize(200, 200);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setVisible(true);
		
		File file = new File("/home/kevin/OpenSHA/UCERF4/rup_sets/fm3_1_ucerf3.zip");
		ZipFile zip = new ZipFile(file);
		
		Runnable closeRunnable = new Runnable() {
			
			@Override
			public void run() {
				try {
					zip.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		};
		Runtime.getRuntime().addShutdownHook(new Thread(closeRunnable));
	}
	
	private static void test94() throws IOException, DocumentException {
//		FaultSystemRupSet rupSet  = FaultSystemIO.loadRupSet(
//				new File("/home/kevin/OpenSHA/UCERF4/rup_sets/"
//						+ "nshm23_v1p2_all_plausibleMulti15km_adaptive6km_direct_cmlRake360_jumpP0.001_"
//						+ "slipP0.05incrCapDist_cff0.75IntsPos_comb2Paths_cffFavP0.01_cffFavRatioN2P0.5_"
//						+ "sectFractGrow0.1.zip"));
//		System.out.println("Done loading");
//		List<ClusterRupture> rups = rupSet.getClusterRuptures();
//		System.out.println("Have "+rups.size()+" cluster ruptures, now calling get(0)");
//		rups.get(0);
//		System.out.println("Done with get(0)");
	}
	
	private static void test95() {
		double lat = 37.09084333164281;
		double lon = -121.96972173834948;
//		Preconditions.checkState(lat == new Location(lat, lon).getLatitude());
		Location loc1 = new Location(lat, lon);
		Location loc2 = new Location(loc1.getLatitude(), loc1.getLongitude());
		System.out.println("INPUT lat/lon");
		System.out.println(lat+"\t"+lon);
		System.out.println("ORIG loc.get*()");
		System.out.println(loc1.getLatitude()+"\t"+loc1.getLongitude());
		System.out.println("ORIG loc.get*Rad()");
		System.out.println(loc1.getLatRad()+"\t"+loc1.getLonRad());
		System.out.println("RECONSTRUCTED loc.get*()");
		System.out.println(loc2.getLatitude()+"\t"+loc2.getLongitude());
		System.out.println("RECONSTRUCTED loc.get*Rad()");
		System.out.println(loc2.getLatRad()+"\t"+loc2.getLonRad());
	}
	
	private static void test96() {
		double lat = 37.09084333164281;
		System.out.println("INPUT DEG\t"+lat);
		System.out.println("Using Java 9+ Math.toRadians/toDegrees:");
		double rad = Math.toRadians(lat);
		double deg = Math.toDegrees(rad);
		System.out.println("\tRAD:\t"+rad);
		System.out.println("\tDEG:\t"+deg);
		System.out.println("Using /180*PI and *180/PI");
		rad = lat/180d * Math.PI;
		deg = rad*180d / Math.PI;
		System.out.println("\tRAD:\t"+rad);
		System.out.println("\tDEG:\t"+deg);
		System.out.println("Using BigDecimal");
		BigDecimal PI = new BigDecimal(
		        "3.14159265358979323846264338327950288419716939937510" +
		        "5820974944592307816406286208998628034825342117067982");
		rad = new BigDecimal(lat).divide(new BigDecimal(180), 100, RoundingMode.HALF_UP).multiply(PI).doubleValue(); 
		deg = new BigDecimal(rad).multiply(new BigDecimal(180)).divide(PI, 100, RoundingMode.HALF_UP).doubleValue(); 
		System.out.println("\tRAD:\t"+rad);
		System.out.println("\tDEG:\t"+deg);
	}
	
	private static void test97() {
		BigDecimal PI = new BigDecimal(
				"3.14159265358979323846264338327950288419716939937510" +
				"5820974944592307816406286208998628034825342117067982");

		double DEG_TO_RAD = PI.divide(new BigDecimal(180), 100, RoundingMode.HALF_UP).doubleValue();
		double RAD_TO_DEG = new BigDecimal(180).divide(PI, 100, RoundingMode.HALF_UP).doubleValue(); 

		Random random = new Random(0);

		int toRadBetterTechniqueWins = 0;
		int toRadBetterTechniqueLoses = 0;
		int toDegBetterTechniqueWins = 0;
		int toDegBetterTechniqueLoses = 0;

		for (int i = 0; i < 10000; i++) {
			double degrees = random.nextInt(360) + random.nextDouble();
//			double degrees = 37.09084333164281;

			double standard = degrees/180d * Math.PI;
			double better = degrees * DEG_TO_RAD;
			double best = new BigDecimal(degrees)
					.divide(new BigDecimal(180), 100, RoundingMode.HALF_UP)
					.multiply(PI)
					.doubleValue();

			double standardError = Math.abs(best - standard);
			double betterError = Math.abs(best - better);

			if (betterError < standardError) {
				toRadBetterTechniqueWins++;
			} else if (betterError > standardError) {
				toRadBetterTechniqueLoses++;
			}

			if (!(standard == better && better == best)) {
				System.out.println(
						degrees + " => " +
								"Standard: " + standard + "; " +
								"Better: " + better + "; " +
								"Best: " + best);
			}
		}

		for (int i = 0; i < 10000; i++) {
			double radians = random.nextDouble() * (2 * Math.PI);
//			double radians = 0.6473573384785501;

			double standard = radians*180d / Math.PI;

			// seemingly more accurate to divide than to multiply by the reciprocal:
			double better = radians / DEG_TO_RAD;
			//double better = radians * RAD_TO_DEG;

			double best = new BigDecimal(radians)
					.multiply(new BigDecimal(180))
					.divide(PI, 100, RoundingMode.HALF_UP)
					.doubleValue();

			double standardError = Math.abs(best - standard);
			double betterError = Math.abs(best - better);

			if (betterError < standardError) {
				toDegBetterTechniqueWins++;
			} else if (betterError > standardError) {
				toDegBetterTechniqueLoses++;
			}

			if (!(standard == better && better == best)) {
				System.out.println(
						radians + " => " +
								"Standard: " + standard + "; " +
								"Better: " + better + "; " +
								"Best: " + best);
			}
		}

		System.out.println("When converting from degrees to radians:");
		System.out.println("Better technique wins: " + toRadBetterTechniqueWins);
		System.out.println("Better technique loses: " + toRadBetterTechniqueLoses);

		System.out.println("When converting from radians to degrees:");
		System.out.println("Better technique wins: " + toDegBetterTechniqueWins);
		System.out.println("Better technique loses: " + toDegBetterTechniqueLoses);

		//System.out.println("DEG_TO_RAD = " + DEG_TO_RAD);
		//System.out.println("RAD_TO_DEG = " + RAD_TO_DEG); 
	}
	
	private static void test98() {
		double[] vals = {
				0.004557008898170912,
				0.006147943144848745,
				5.52179640057068E-4,
				4.279434356649868E-7,
				0.0026159678814425035,
				1.7629342361803746E-8,
				5.78051085579677E-4,
				6.305592933062744E-8,
				2.865482448602668E-7,
				3.77428890089215E-8,
				8.701321712129223E-7,
				Double.NaN,
				Double.POSITIVE_INFINITY,
				Double.NEGATIVE_INFINITY
		};
		for (double val : vals) {
			System.out.println("Orig:\t"+val);
			double rVal = DataUtils.roundFixed(val, 8, RoundingMode.HALF_EVEN);
			double backVal = Double.parseDouble(rVal+"");
			System.out.println("\t"+rVal+"\t->\t"+backVal);
			rVal = DataUtils.roundSigFigs(val, 8, RoundingMode.HALF_EVEN);
			backVal = Double.parseDouble(rVal+"");
			System.out.println("\t"+rVal+"\t->\t"+backVal);
		}
	}
	
	private static void test99() throws IOException {
//		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/UCERF4/rup_sets/fm3_1_ucerf3.zip"));
//		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/UCERF3/rup_sets/modular/"
//				+ "FM3_1_ZENGBB_Shaw09Mod_DsrUni_CharConst_M5Rate7.9_MMaxOff7.6_NoFix_SpatSeisU3.zip"));
		CompoundFaultSystemSolution cfss = CompoundFaultSystemSolution.fromZipFile(
				new File("/home/kevin/OpenSHA/UCERF3/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip"));
		FaultSystemSolution sol = cfss.getSolution(U3LogicTreeBranch.fromValues(true, SlipAlongRuptureModels.UNIFORM));
		System.out.println(sol.requireModule(GridSourceProvider.class).getNodeUnassociatedMFD(0));
	}
	
	private static void test100() throws IOException {
		File inDir = new File("/data/kevin/markdown/inversions/2021_09_13-coulomb-nshm23_geol_dm_v1-slip_constr-gr_constr-1hr");
//		File inRupSetFile = new File(inDir, "rupture_set.zip");
		File inRupSetFile = new File(inDir, "solution2.zip");
		File solFile = new File(inDir, "solution.zip");
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(inRupSetFile);
		FaultSystemSolution origSol = FaultSystemSolution.load(solFile);
		FaultSystemSolution newSol = new FaultSystemSolution(rupSet, origSol.getRateForAllRups());
		for (OpenSHA_Module module : origSol.getModules())
			newSol.addModule(module);
		newSol.write(new File(inDir, "solution2.zip"));
	}
	
	private static void test101() throws IOException {
		File baseDir = new File("/home/kevin/OpenSHA/UCERF4/fault_models/");
		
		File fmFile = new File(baseDir, "NSHM23_FaultSections_v1p1/NSHM23_FaultSections_v1p1.geojson");
		File dmFile = new File(baseDir, "NSHM23_GeolDefMod_v1/NSHM23_GeolDefMod_v1.csv");
		
		List<GeoJSONFaultSection> sects = GeoJSONFaultReader.readFaultSections(fmFile);
		System.out.println("Loaded "+sects.size()+" sections");
		
		GeoJSONFaultReader.attachGeoDefModel(sects, dmFile);
		
		double confStdDevsAway = 2d;

		MinMaxAveTracker stdDevTrack = new MinMaxAveTracker();
		MinMaxAveTracker fractTrack = new MinMaxAveTracker();
		
		DefaultXY_DataSet xy = new DefaultXY_DataSet();

		List<Double> stdDevs = new ArrayList<>();
		List<Double> fracts = new ArrayList<>();
		
		for (GeoJSONFaultSection sect : sects) {
			double high = sect.getProperty("HighRate", Double.NaN);
			double low = sect.getProperty("LowRate", Double.NaN);
			double mean = sect.getOrigAveSlipRate();
			if (mean == 0d)
				continue;
			
			double upperStdDev = (high-mean)/confStdDevsAway;
			double lowerStdDev = (mean-low)/confStdDevsAway;
//			double comb = 0.5*(upperStdDev + lowerStdDev);
			double comb = (high-low)/(2*confStdDevsAway);
			stdDevTrack.addValue(comb);
			stdDevs.add(comb);
			double fract = comb/mean;
			fractTrack.addValue(fract);
			fracts.add(fract);
			
			System.out.println(sect.getSectionId()+". "+sect.getSectionName());
			System.out.println("\tSlip Rate:\t"+(float)mean+" ["+(float)low+", "+(float)high+"]");
			System.out.println("\tImplied Std Dev:\t"+(float)comb+" ("+(float)lowerStdDev+", "+(float)upperStdDev+")");
			System.out.println("\tImplied Fractional:\t"+(float)(fract));
			
			xy.set(mean, fract);
		}
		System.out.println("Standard Deviations: "+stdDevTrack);
		System.out.println("\tMedian: "+DataUtils.median(Doubles.toArray(stdDevs)));
		System.out.println("Standard Deviation Fractions: "+fractTrack);
		System.out.println("\tMedian: "+DataUtils.median(Doubles.toArray(fracts)));
		
		GraphWindow gw = new GraphWindow(xy, "Fractional Std Dev vs Slip Rate",
				new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, Color.BLACK));
		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
	}
	
	private static void test102() throws IOException {
		File solFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2021_10_19-reproduce-ucerf3-ref_branch-tapered-convergence-u3Iters/mean_solution.zip");
		FaultSystemSolution sol = FaultSystemSolution.load(solFile);
		
		WaterLevelRates wlRates = sol.getModule(WaterLevelRates.class);
		double[] wl = wlRates.get();
		
		for (int r=0; r<wl.length; r++) {
			double rate = sol.getRateForRup(r);
			double noMinRate = rate-wl[r];
			if (noMinRate < 1e-16 && Math.random() < 0.01)
				System.out.println(r+"\ttotRate="+rate+"\twl="+wl[r]+"\tnoMinRate="+noMinRate);
		}
	}
	
	private static void test103() throws IOException {
		NormalDistribution normDist = new NormalDistribution(0d, 1d);
		MinMaxAveTracker track = new MinMaxAveTracker();
		MinMaxAveTracker absTrack = new MinMaxAveTracker();
		for (int r=0; r<100000; r++) {
			double v = normDist.sample();
			track.addValue(v);
			absTrack.addValue(Math.abs(v));
		}
		System.out.println("Values: "+track);
		System.out.println("Abs Values: "+absTrack);
	}
	
	private static void test104() throws IOException {
		File solFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2021_10_19-reproduce-ucerf3-ref_branch-tapered-convergence-u3Iters/mean_solution.zip");
		
		System.out.println("Loading as rupture set...");
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(solFile);
		rupSet.setParent(null);
		
		System.out.println();
		System.out.println("Loading as solution reusing rup set...");
		FaultSystemSolution sol = FaultSystemSolution.load(solFile, rupSet);
		Preconditions.checkState(sol.getRupSet() == rupSet);
	}
	
	private static void test105() throws IOException {
		File solFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2021_10_19-reproduce-ucerf3-ref_branch-tapered-convergence-u3Iters/mean_solution.zip");
		
		System.out.println("Loading as rupture set...");
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(solFile);
		
		rupSet.requireModule(InversionTargetMFDs.class);
	}
	
	private static void test106() {
		GutenbergRichterMagFreqDist refMFD = new GutenbergRichterMagFreqDist(6.05, 30, 0.1);
		double targetB = 1.0;
		double moRate = 1e17;
		refMFD.setAllButTotCumRate(refMFD.getMinX(), refMFD.getMaxX(), moRate, targetB);
		
		double refRate = refMFD.getTotCumRate();
		
		System.out.println("Total rate w/ b="+(float)targetB+":\t"+refRate);
		
		double fakeB = targetB+0.5;
		
		GutenbergRichterMagFreqDist fakeGR = new GutenbergRichterMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		fakeGR.setAllButTotCumRate(refMFD.getMinX(), refMFD.getMaxX(), moRate, fakeB);
		
		double fakeRate = fakeGR.getTotCumRate();
		
		System.out.println("Total rate w/ b="+(float)fakeB+":\t"+fakeRate);
		
		double oppositeB = targetB - (fakeB-targetB);
		
		GutenbergRichterMagFreqDist oppGR = new GutenbergRichterMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		oppGR.setAllButTotCumRate(refMFD.getMinX(), refMFD.getMaxX(), moRate, oppositeB);
		
		double oppRate = oppGR.getTotCumRate();
		
		System.out.println("Opposite rate w/ b="+(float)oppositeB+":\t"+oppRate);
	}
	
	private static void test107() {
		GutenbergRichterMagFreqDist refMFD = new GutenbergRichterMagFreqDist(6.05, 30, 0.1);
		double targetB = 1.0;
		double moRate = 1e17;
		refMFD.setAllButTotCumRate(refMFD.getMinX(), refMFD.getMaxX(), moRate, targetB);
		
		GutenbergRichterMagFreqDist charMFD = new GutenbergRichterMagFreqDist(6.05, 30, 0.1);
		charMFD.setAllButTotCumRate(charMFD.getMaxX(), charMFD.getMaxX(), moRate, targetB);
		
		System.out.println("Ref MFD: "+refMFD);
		System.out.println("Char MFD: "+charMFD);
		System.out.println("Char moment: "+MagUtils.magToMoment(charMFD.getMaxX())*charMFD.getY(charMFD.size()-1));
	}
	
	private static void test108() throws IOException {
//		File inputFile = new File("/home/kevin/.opensha/ucerf3_erf/cached_dep100.0_depMean_rakeMean.zip");
		File inputFile = new File("/home/kevin/.opensha/ucerf3_erf/cached_FM3_1_dep100.0_depMean_rakeMean.zip");
		File outputFile = new File(inputFile.getParentFile(), "modular_"+inputFile.getName());
		FaultSystemSolution.load(inputFile).write(outputFile);
	}
	
	private static void test109() throws IOException {
		UncertainDataConstraint prior = new UncertainDataConstraint("Prior", 0.1,
				UncertaintyBoundType.ONE_SIGMA.estimate(0.1, 0.01));
		
		UncertainDataConstraint paleo = new UncertainDataConstraint("Paleo", 0.5,
				UncertaintyBoundType.ONE_SIGMA.estimate(0.5, 0.1));
		
		System.out.println("Prior:\t"+prior);
		System.out.println("Paleo:\t"+paleo);
		
		double priorSD_2 = Math.pow(prior.getPreferredStdDev(), 2);
		double paleoSD_2 = Math.pow(paleo.getPreferredStdDev(), 2);
		
		
//		double meanEst = (priorSD_2*prior.bestEstimate + paleoSD_2*paleo.bestEstimate)/(priorSD_2 + paleoSD_2);
//		double meanEst = (priorSD_2*paleo.bestEstimate + paleoSD_2*prior.bestEstimate)/(priorSD_2 + paleoSD_2);
//		double sdEst = Math.sqrt((priorSD_2*paleoSD_2)/(priorSD_2 + paleoSD_2));
		
		
		// from https://stats.stackexchange.com/a/15273
//		double nPrior = 10;
//		double nPaleo = 10;
//		double meanEst = (nPrior*prior.bestEstimate/priorSD_2 + nPaleo*paleo.bestEstimate/paleoSD_2)/(nPrior/priorSD_2 + nPaleo/paleoSD_2);
//		double sdEst = Math.sqrt(1d/(nPrior/priorSD_2 + nPaleo/paleoSD_2));
		
//		int nPriorSamples = 1000000;
//		int nPaleoSamples = 1000000;
//		NormalDistribution priorNorm = new NormalDistribution(prior.bestEstimate, prior.getPreferredStdDev());
//		NormalDistribution paleoNorm = new NormalDistribution(paleo.bestEstimate, paleo.getPreferredStdDev());
//		double[] vals = new double[nPriorSamples+nPaleoSamples];
//		for (int i=0; i<vals.length; i++) {
//			if (i < nPriorSamples)
//				vals[i] = priorNorm.sample();
//			else
//				vals[i] = paleoNorm.sample();
//		}
//		double meanEst = StatUtils.mean(vals);
//		double sdEst = Math.sqrt(StatUtils.variance(vals));
		
		double meanEst = 0.5*(prior.bestEstimate+paleo.bestEstimate);
		double sdEst = 0.5*Math.sqrt(priorSD_2 + paleoSD_2);
		
		UncertainDataConstraint updated = new UncertainDataConstraint("Updated",
				meanEst, UncertaintyBoundType.ONE_SIGMA.estimate(meanEst, sdEst));
		System.out.println("Updated:\t"+updated);
	}
	
	private static void test110() {
		double[] rates = new double[10];
		double relStdDev = 0.1;
		for (int i=0; i<rates.length; i++)
			rates[i] = Math.random()*Math.random()*Math.random()*10000;
//			rates[i] = 1;
		
		double[] stdDevs = new double[rates.length];
		for (int i=0; i<rates.length; i++)
			stdDevs[i] = rates[i]*relStdDev;
		
		double sumVar = 0d;
		for (double stdDev : stdDevs)
			sumVar += stdDev*stdDev;
		
		double sqrtSumVar = Math.sqrt(sumVar);
		
		double sumRates = StatUtils.sum(rates);
		
		double sumSqrtRel = sqrtSumVar/sumRates;
		System.out.println("Sum rates: "+sumRates);
		System.out.println("Sum variance: "+sumVar);
		System.out.println("Sqrt(sum var): "+sqrtSumVar);
		System.out.println("\tReltaive: "+sumSqrtRel);
		
		double sumStdDev = StatUtils.sum(stdDevs);
		double sumStdDevRel = sumStdDev/sumRates;
		System.out.println("Sum of standard deviations: "+sumStdDev);
		System.out.println("\tReltaive: "+sumStdDevRel);
	}
	
	private static void test111() {
		double impliedRate = 0.0060084797;
		double constrRate = 3.1561314E-4;
	}
	
	private static void test112() {
		Options ops = new Options();
		ops.addRequiredOption("cl", "class", true, "asdf");
		
		String[] args = {"--class", "scratch.kevin.nshm23.U3InversionConfigFactory$NoPaleoParkfieldSingleReg"};
		CommandLine cmd = FaultSysTools.parseOptions(ops, args, PureScratch.class);
		System.out.println("Class: "+cmd.getOptionValue("class"));
		
		InversionConfigurationFactory factory;
		try {
			@SuppressWarnings("unchecked")
			Class<? extends InversionConfigurationFactory> factoryClass = (Class<? extends InversionConfigurationFactory>)
					Class.forName(cmd.getOptionValue("class"));
			factory = factoryClass.getDeclaredConstructor().newInstance();
		} catch (Exception e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		System.out.println("Factory type: "+factory.getClass().getName());
//		U3InversionConfigFactory factory = new U3InversionConfigFactory.NoPaleoParkfieldSingleReg();
		
		U3LogicTreeBranch branch = U3LogicTreeBranch.DEFAULT;
		FaultSystemRupSet rupSet = factory.buildRuptureSet(branch, 32);
		InversionConfiguration config = factory.buildInversionConfig(rupSet, branch, 32);
		for (InversionConstraint constraint : config.getConstraints())
			System.out.println(constraint.getName()+" has "+constraint.getNumRows()+" rows");
	}
	
	private static void test113() {
		double supraB = 0.5;
		double subB = 1d;
		double totMoRate = 1e17;
		
		EvenlyDiscretizedFunc refXVals = new EvenlyDiscretizedFunc(0.05, 78, 0.1);
		int supraIndex = 72;
		System.out.println("Max m="+refXVals.getMaxX());
		System.out.println("Supra m="+refXVals.getX(supraIndex));
		
		// start with a full G-R with the supra b-value
		GutenbergRichterMagFreqDist fullSupraB = new GutenbergRichterMagFreqDist(
				refXVals.getMinX(), refXVals.size(), refXVals.getDelta(), totMoRate, supraB);
		
		// copy it to a regular MFD:
		IncrementalMagFreqDist fullMFD = new IncrementalMagFreqDist(refXVals.getMinX(), refXVals.size(), refXVals.getDelta());
		for (int i=0; i<fullSupraB.size(); i++)
			fullMFD.set(i, fullSupraB.getY(i));
		
		System.out.println("Orig mo rate: "+fullMFD.getTotalMomentRate());
		
		// now correct the sub-seis portion to have the sub-seis b-value
		
		// first create a full MFD with the sub b-value. this will only be used in a relative sense
		GutenbergRichterMagFreqDist fullSubB = new GutenbergRichterMagFreqDist(
				refXVals.getMinX(), refXVals.size(), refXVals.getDelta(), totMoRate, subB);
		
		double targetFirstSupra = fullSupraB.getY(supraIndex);
		double subFirstSupra = fullSubB.getY(supraIndex);
		for (int i=0; i<supraIndex; i++) {
			double targetRatio = fullSubB.getY(i)/subFirstSupra;
			fullMFD.set(i, targetFirstSupra*targetRatio);
		}
		
		System.out.println("Pre-scaled mo rate: "+fullMFD.getTotalMomentRate());
		
		fullMFD.scaleToTotalMomentRate(totMoRate);
		
		System.out.println("Scaled mo rate: "+fullMFD.getTotalMomentRate());
		
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		funcs.add(fullSupraB);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		funcs.add(fullSubB);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE));
		
		funcs.add(fullMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.RED));
		
		GraphWindow gw = new GraphWindow(funcs, "MFD test", chars);
		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
		
		gw.setYLog(true);
	}
	
	private static void test114() throws IOException {
		System.gc();
		Runtime rt = Runtime.getRuntime();
		long totalMB = rt.totalMemory() / 1024 / 1024;
		long freeMB = rt.freeMemory() / 1024 / 1024;
		long usedMB = totalMB - freeMB;
		System.out.println("mem t/u/f: "+totalMB+"/"+usedMB+"/"+freeMB);
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/tmp/solution.zip"));
		totalMB = rt.totalMemory() / 1024 / 1024;
		freeMB = rt.freeMemory() / 1024 / 1024;
		usedMB = totalMB - freeMB;
		System.out.println("mem t/u/f: "+totalMB+"/"+usedMB+"/"+freeMB);
		System.gc();
		totalMB = rt.totalMemory() / 1024 / 1024;
		freeMB = rt.freeMemory() / 1024 / 1024;
		usedMB = totalMB - freeMB;
		System.out.println("mem t/u/f: "+totalMB+"/"+usedMB+"/"+freeMB);
		System.out.println(sol.getTotalRateForAllFaultSystemRups());
	}
	
	private static void test115() throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2021_12_17-nshm23_draft_branches-FM3_1-CoulombRupSet/node_branch_averaged/SegModel_ShawR0_2.zip"));
		SolMFDPlot plot = new SolMFDPlot();
		File outputDir = new File("/tmp/mfd_test");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		File withDir = new File(outputDir, "with");
		Preconditions.checkState(withDir.exists() || withDir.mkdir());
		plot.writePlot(sol.getRupSet(), sol, "Orig", withDir);
		Preconditions.checkState(sol.hasModule(RupMFDsModule.class));
		sol.removeModuleInstances(RupMFDsModule.class);
		File withoutDir = new File(outputDir, "without");
		Preconditions.checkState(withoutDir.exists() || withoutDir.mkdir());
		plot.writePlot(sol.getRupSet(), sol, "Without", withoutDir);
	}
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		test115();
	}

}
