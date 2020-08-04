package scratch.kevin.ucerf3;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.InetAddress;
import java.net.UnknownHostException;
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
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.dom4j.DocumentException;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.region.CaliforniaRegions.RELM_SOCAL;
import org.opensha.commons.data.siteData.impl.WillsMap2006;
import org.opensha.commons.data.siteData.impl.WillsMap2015;
import org.opensha.commons.eq.MagUtils;
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
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.param.editor.AbstractParameterEditor;
import org.opensha.commons.param.editor.impl.NumericTextField;
import org.opensha.commons.param.impl.DoubleParameter;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.FileUtils;
import org.opensha.commons.util.ServerPrefUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.calc.hazardMap.HazardCurveSetCalculator;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.calc.ERF_Calculator;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.earthquake.observedEarthquake.parsers.UCERF3_CatalogParser;
import org.opensha.sha.earthquake.param.BackgroundRupParam;
import org.opensha.sha.earthquake.param.BackgroundRupType;
import org.opensha.sha.earthquake.param.HistoricOpenIntervalParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.StirlingGriddedSurface;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenuationRelationship;
import org.opensha.sha.imr.attenRelImpl.USGS_Combined_2004_AttenRel;
import org.opensha.sha.imr.attenRelImpl.ngaw2.ASK_2014;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.utils.RSQSimSubSectionMapper;
import org.opensha.sha.simulators.utils.RSQSimUtils;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.Files;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.analysis.FaultSysSolutionERF_Calc;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.TotalMag5Rate;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO.BinarayCatalogsIterable;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.ETAS_MultiSimAnalysisTools;
import scratch.UCERF3.erf.ETAS.ETAS_SimAnalysisTools;
import scratch.UCERF3.erf.ETAS.ETAS_Simulator.TestScenario;
import scratch.UCERF3.erf.ETAS.ETAS_Utils;
import scratch.UCERF3.erf.ETAS.association.FiniteFaultMappingData;
import scratch.UCERF3.erf.mean.MeanUCERF3;
import scratch.UCERF3.erf.mean.MeanUCERF3.Presets;
import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;
import scratch.UCERF3.griddedSeismicity.GridSourceProvider;
import scratch.UCERF3.inversion.CommandLineInversionRunner;
import scratch.UCERF3.inversion.InversionFaultSystemRupSetFactory;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.utils.DeformationModelFetcher;
import scratch.UCERF3.utils.FaultSectionDataWriter;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.utils.LastEventData;
import scratch.UCERF3.utils.UCERF3_DataUtils;
import scratch.UCERF3.utils.aveSlip.AveSlipConstraint;
import scratch.UCERF3.utils.paleoRateConstraints.PaleoRateConstraint;
import scratch.UCERF3.utils.paleoRateConstraints.UCERF3_PaleoRateConstraintFetcher;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

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
		FaultSystemSolution sol = FaultSystemIO.loadSol(
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
		FaultSystemSolution sol = FaultSystemIO.loadSol(
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
		FaultSystemSolution sol_31 = FaultSystemIO.loadSol(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		FaultSystemSolution sol_32 = FaultSystemIO.loadSol(
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
		FaultSystemRupSet rupSet = FaultSystemIO.loadRupSet(
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
		FaultSystemRupSet rupSet = FaultSystemIO.loadRupSet(
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
		FaultSystemRupSet rupSet = FaultSystemIO.loadRupSet(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
						+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		int id = 193821;
		System.out.println(Joiner.on(",").join(rupSet.getSectionsIndicesForRup(id)));
	}

	private static void test11() throws IOException, DocumentException {
		FaultSystemSolution sol = FaultSystemIO.loadSol(
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
		FaultSystemSolution sol = FaultSystemIO.loadSol(
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
		InversionFaultSystemSolution sol = FaultSystemIO.loadInvSol(solFile);
		File outputDir = new File(solFile.getParentFile(), solFile.getName().replace(".zip", ""));
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		ArrayList<PaleoRateConstraint> paleoRateConstraints =
				CommandLineInversionRunner.getPaleoConstraints(sol.getRupSet().getFaultModel(), sol.getRupSet());
		System.out.println(paleoRateConstraints.size()+" paleo constraints");
		List<AveSlipConstraint> aveSlipConstraints = AveSlipConstraint.load(sol.getRupSet().getFaultSectionDataList());
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
		
		FaultSystemSolution sol = FaultSystemIO.loadSol(
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
		
		FaultSystemSolution sol = FaultSystemIO.loadSol(
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
		FaultSystemSolution sol1 = FaultSystemIO.loadSol(new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		FaultSystemSolution sol2 = FaultSystemIO.loadSol(new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
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
		FaultSystemRupSet rupSet = FaultSystemIO.loadRupSet(new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		int index = TestScenario.MOJAVE_M7.getFSS_Index();
		for (FaultSection sect : rupSet.getFaultSectionDataForRupture(index))
			System.out.println(sect.getSectionId()+": "+sect.getSectionName());
	}
	
	private static void test35() throws IOException, DocumentException {
		FaultSystemSolution sol = FaultSystemIO.loadSol(
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
		FaultSystemSolution sol = FaultSystemIO.loadSol(
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
		FaultSystemSolution sol = FaultSystemIO.loadSol(
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
		FaultSystemSolution sol = FaultSystemIO.loadSol(
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
		RuptureSurface surf = sol.getRupSet().getSurfaceForRupupture(index, 1d);
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
		FaultSystemSolution sol = FaultSystemIO.loadSol(
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
		List<PaleoRateConstraint> biasiScharerSites = new ArrayList<>();
		for (PaleoRateConstraint constraint : UCERF3_PaleoRateConstraintFetcher.getConstraints()) {
			String name = constraint.getPaleoSiteName();
			if (name.equals("S. San Andreas - Coachella") || name.equals("San Jacinto - Hog Lake")
					|| name.equals("Frazier Mountian, SSAF") || name.equals("N. San Andreas - Santa Cruz Seg.")
					|| name.equals("Hayward fault - South"))
				biasiScharerSites.add(constraint);
		}
		for (PaleoRateConstraint constraint : biasiScharerSites)
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
		rupSetMap.put(FaultModels.FM3_1, FaultSystemIO.loadRupSet(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip")));
		rupSetMap.put(FaultModels.FM3_2, FaultSystemIO.loadRupSet(
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

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		test72();

		////		FaultSystemSolution sol3 = FaultSystemIO.loadSol(new File("/tmp/avg_SpatSeisU3/"
		////				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		////		System.out.println(sol3.getSubSeismoOnFaultMFD_List().size());
		////		System.exit(0);
		////		CompoundFaultSystemSolution cfss2 = CompoundFaultSystemSolution.fromZipFile(new File(
		////				"/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
		////				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip"));
		////		FaultSystemSolution sol1 = cfss2.getSolution(LogicTreeBranch.DEFAULT);
		//		FaultSystemSolution inputSol = FaultSystemIO.loadSol(
		////				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
		////					+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		//				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/UCERF3_ERF/"
		//						+ "cached_dep10.0_depMean_rakeGEOLOGIC.zip"));
		//		
		//		Region region = new CaliforniaRegions.RELM_SOCAL();
		//		
		//		FaultSystemSolution sol1 = CSDownsampledSolCreator.getDownsampledSol(inputSol, region);
		//		
		//		double[] counts = new double[sol1.getRupSet().getNumSections()];
		//		for (int sectIndex=0; sectIndex<sol1.getRupSet().getNumSections(); sectIndex++)
		//			counts[sectIndex] = sol1.getRupSet().getRupturesForSection(sectIndex).size();
		//		HistogramFunction numHist = HistogramFunction.getEncompassingHistogram(0d, 50000d, 1000d);
		////		HistogramFunction numHist = new HistogramFunction(50d, 10050d, (100));
		//		for (double count : counts)
		//			numHist.add(numHist.getClosestXIndex(count), 1d);
		//		
		//		List<DiscretizedFunc> funcs = Lists.newArrayList();
		//		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		//		
		//		funcs.add(numHist);
		//		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));
		//		
		//		GraphWindow gw = new GraphWindow(funcs, "Ruptures Per Subsection", chars);
		//		gw.setX_AxisLabel("Num Ruptures Including Subsection");
		//		gw.setY_AxisLabel("Num Subsections");
		//		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
		//		gw.saveAsPNG("/tmp/rup_count_hist.png");
		//		
		//		// now fract
		//		EvenlyDiscretizedFunc fractFunc = new EvenlyDiscretizedFunc(1d, 100000, 1d);
		//		Arrays.sort(counts);
		//		for (int i=0; i<fractFunc.size(); i++) {
		//			double x = fractFunc.getX(i);
		//			int index = Arrays.binarySearch(counts, x);
		//			if (index < 0) {
		//				index = -(index + 1);
		//			}
		//			double y = 1d - (double)(index)/(double)(counts.length-1);
		//			fractFunc.set(i, y);
		//		}
		//		
		//		funcs = Lists.newArrayList();
		//		chars = Lists.newArrayList();
		//		
		//		funcs.add(fractFunc);
		//		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		//		
		//		gw = new GraphWindow(funcs, "Ruptures Per Subsection", chars);
		//		gw.setX_AxisLabel("Num Ruptures Including Subsection");
		//		gw.setY_AxisLabel("Fraction of Subsections >= Num");
		//		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
		//		gw.saveAsPNG("/tmp/rup_count_fract_func.png");
		//		
		//		while (" ".length() > 0) {
		//			try {
		//				Thread.sleep(1000);
		//			} catch (InterruptedException e) {
		//				// TODO Auto-generated catch block
		//				e.printStackTrace();
		//			}
		//		}

		//		int subSect = 1267;
		//		System.out.println(sol1.getSubSeismoOnFaultMFD_List().get(subSect));
		//		int numContained = 0;
		//		Region reg = sol1.getGridSourceProvider().getGriddedRegion();
		//		System.out.println("Poly overlap? "+(Region.intersect(reg,
		//				sol1.getRupSet().getInversionTargetMFDs().getGridSeisUtils().getPolyMgr().getPoly(subSect)) != null));
		////		((UCERF3_GridSourceGenerator)sol1.getGridSourceProvider()).get
		//		for (Location loc : sol1.getRupSet().getFaultSectionData(subSect).getFaultTrace()) {
		//			if (reg.contains(loc))
		//				numContained++;
		//		}

		//		System.out.println(numContained+" points contained in "+reg.getName());
		//		File sol1File = new File("/tmp/sol1.zip");
		//		FaultSystemIO.writeSol(sol1, sol1File);
		//		FaultSystemSolution sol2 = FaultSystemIO.loadSol(sol1File);
		//		System.out.println("# sub seismos: "+sol2.getSubSeismoOnFaultMFD_List().size());
		//		System.out.println("Orig first: "+sol1.getSubSeismoOnFaultMFD_List().get(0));
		//		System.out.println("New first: "+sol2.getSubSeismoOnFaultMFD_List().get(0));
		//		System.exit(0);
		//		RegionUtils.regionToKML(new CaliforniaRegions.LA_BOX(), "la_box", Color.BLACK);
		//		RegionUtils.regionToKML(new CaliforniaRegions.SF_BOX(), "sf_box", Color.BLACK);
		//		System.exit(0);
		//		FaultSystemSolution theSol = FaultSystemIO.loadSol(
		//				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
		//						+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		//		FaultSystemRupSet rupSet = theSol.getRupSet();
		//		int numWith = 0;
		//		for (FaultSection sect : rupSet.getFaultSectionDataList())
		//			if (sect.getDateOfLastEvent() > Long.MIN_VALUE)
		//				numWith++;
		//		System.out.println(numWith+"/"+rupSet.getNumSections()+" have last event data");
		//		System.exit(0);
		//		HistogramFunction hist = new HistogramFunction(0.25, 40, 0.5);
		//		for (FaultSection sect : theSol.getRupSet().getFaultSectionDataList()) {
		//			double len = sect.getTraceLength();
		//			hist.add(len, 1d);
		//		}
		//		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		//		List<HistogramFunction> elems = Lists.newArrayList();
		//		elems.add(hist);
		//		List<PlotCurveCharacterstics> chars = Lists.newArrayList(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));
		//		PlotSpec spec = new PlotSpec(elems, chars, "Sub Sect Lenghts", "Length (km)", "Number");
		//		gp.drawGraphPanel(spec);
		//		gp.getCartPanel().setSize(1000, 800);
		//		gp.setBackground(Color.WHITE);
		//		gp.saveAsPNG("/tmp/sub_sect_length_hist.png");
		//		System.exit(0);
		//		CompoundFaultSystemSolution cfss = CompoundFaultSystemSolution.fromZipFile(new File("/tmp/compound/2013_05_10-ucerf3p3-production-10runs_run1_COMPOUND_SOL.zip"));
		//		List<LogicTreeBranch> branches = Lists.newArrayList(cfss.getBranches());
		//		Collections.shuffle(branches);
		//		HashSet<FaultModels> models = new HashSet<FaultModels>();
		//		for (int i=0; i<branches.size(); i++) {
		//			LogicTreeBranch branch = branches.get(i);
		//			FaultModels model = branch.getValue(FaultModels.class);
		//			if (models.contains(model))
		//				continue;
		//			InversionFaultSystemSolution sol = cfss.getSolution(branch);
		//			System.out.println(model.getShortName()+": "+sol.getRupSet().getNumRuptures());
		//			System.out.println(sol.getClass().getName());
		//			if (sol instanceof AverageFaultSystemSolution)
		//				System.out.println(((AverageFaultSystemSolution)sol).getNumSolutions()+" sols");
		//			models.add(model);
		//		}
		////		InversionFaultSystemSolution sol = cfss.getSolution(cfss.getBranches().iterator().next());
		////		CommandLineInversionRunner.writeParentSectionMFDPlots(sol, new File("/tmp/newmfd"));
		//		System.exit(0);
		//		File f = new File("/tmp/FM3_2_ZENGBB_EllBsqrtLen_DsrTap_CharConst_M5Rate9.6_MMaxOff7.6_NoFix_SpatSeisU2_run0_sol.zip");
		//		InversionFaultSystemSolution invSol = FaultSystemIO.loadInvSol(f);
		//		System.out.println(invSol.getClass());
		//		System.out.println(invSol.getLogicTreeBranch());
		//		System.out.println(invSol.getInversionConfiguration());
		//		System.out.println(invSol.getInvModel());
		//		System.out.println(invSol.getMisfits().size());
		//		invSol.getGridSourceProvider();
		//		
		//		cfss = CompoundFaultSystemSolution.fromZipFile(new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/2013_05_03-ucerf3p3-production-first-five_MEAN_COMPOUND_SOL.zip"));
		//		invSol = cfss.getSolution(invSol.getLogicTreeBranch());
		//		invSol.getGridSourceProvider();
		//		System.out.println("Got it from a grid source provider");
		//		System.exit(0);
		//		
		//		System.out.println("LEGACY MEAN SOLUTION");
		//		f = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/2013_01_14-stampede_3p2_production_runs_combined_FM3_1_MEAN_BRANCH_AVG_SOL.zip");
		//		invSol = FaultSystemIO.loadInvSol(f);
		//		System.out.println(invSol.getLogicTreeBranch());
		//		System.out.println(invSol.getInversionConfiguration());
		//		System.out.println(invSol.getInvModel());
		//		System.out.println(invSol.getMisfits().size());
		//		
		//		System.out.println("LEGACY NORMAL SOLUTION");
		//		f = new File("/tmp/FM3_1_ZENGBB_EllB_DsrUni_CharConst_M5Rate7.9_MMaxOff7.6_NoFix_SpatSeisU3_VarSlipTypeBOTH_VarSlipWtUnNorm100_sol.zip.1");
		//		invSol = FaultSystemIO.loadInvSol(f);
		//		System.out.println(invSol.getLogicTreeBranch());
		//		System.out.println(invSol.getInversionConfiguration());
		//		System.out.println(invSol.getInvModel());
		//		System.out.println(invSol.getMisfits().size());
		//		
		//		// loading in AVFSS
		//		f = new File("/tmp/compound_tests_data/subset/FM3_1_NEOK_Shaw09Mod_DsrUni_CharConst_M5Rate7.9_MMaxOff7.6_NoFix_SpatSeisU3_mean/FM3_1_NEOK_Shaw09Mod_DsrUni_CharConst_M5Rate7.9_MMaxOff7.6_NoFix_SpatSeisU3_mean_sol.zip");
		//		AverageFaultSystemSolution avgSol = FaultSystemIO.loadAvgInvSol(f);
		//		System.out.println("Average sol with "+avgSol.getNumSolutions()+" sols");
		//		System.exit(0);
		//		
		//		double a1 = 8.62552e32;
		//		double a2 = 1.67242e34;
		//		double a3 = 1.77448e20;
		//		double a4 = 9.05759e20;
		//		double c = 8.0021909e37;
		//		double me = 5.9742e24;
		//		double re = 6371000;
		//		
		//		double p1 = (c - a1*me/a3) / (a2-a1*a4/a3);
		//		double p2 = (me - a4*p1)/a3;
		//		
		//		System.out.println("p1: "+p1);
		//		System.out.println("p2: "+p2);
		//		
		//		System.exit(0);
		//		
		//		
		//		CB_2008_AttenRel imr = new CB_2008_AttenRel(null);
		//		imr.setParamDefaults();
		//		MeanUCERF2 ucerf = new MeanUCERF2();
		//		ucerf.updateForecast();
		//		ProbEqkSource src = ucerf.getSource(3337);
		//		System.out.println(src.getName());
		//		ProbEqkRupture theRup = src.getRupture(77);
		//		System.out.println("Rup 77 pts: "+theRup.getRuptureSurface().getEvenlyDiscritizedListOfLocsOnSurface().size());
		//		imr.setEqkRupture(theRup);
		//		imr.setIntensityMeasure(PGA_Param.NAME);
		//		Site site = new Site(new Location(34.10688, -118.22060));
		//		site.addParameterList(imr.getSiteParams());
		//		imr.setAll(theRup, site, imr.getIntensityMeasure());
		//		imr.getMean();
		//		System.exit(0);
		//		
		//		UCERF2_TimeIndependentEpistemicList ti_ep = new UCERF2_TimeIndependentEpistemicList();
		//		UCERF2_TimeDependentEpistemicList td_ep = new UCERF2_TimeDependentEpistemicList();
		////		ep.updateForecast();
		//		System.out.println(ti_ep.getNumERFs());
		//		System.out.println(td_ep.getNumERFs());
		//		
		//		System.exit(0);
		//		
		//		
		//		// this get's the DB accessor (version 3)
		//		DB_AccessAPI db = DB_ConnectionPool.getDB3ReadOnlyConn();
		//
		//		PrefFaultSectionDataDB_DAO faultSectionDB_DAO = new PrefFaultSectionDataDB_DAO(db);
		//
		//		List<FaultSection> sections = faultSectionDB_DAO.getAllFaultSection(); 
		//		for (FaultSection data : sections)
		//			System.out.println(data);
		//		System.exit(0);
		//		
		//		
		//		double minX = -9d;
		//		double maxX = 0d;
		//		int num = 200;
		//		EvenlyDiscretizedFunc ucerf2Func = new EvenlyDiscretizedFunc(minX, maxX, num);
		//		double delta = ucerf2Func.getDelta();
		//		ucerf2Func.setName("UCERF2");
		//		
		//		boolean doUCERF2 = true;
		//		
		//		if (doUCERF2) {
		//			System.out.println("Creating UCERF2");
		//			ERF erf = new MeanUCERF2();
		//			erf.setParameter(UCERF2.PROB_MODEL_PARAM_NAME, UCERF2.PROB_MODEL_POISSON);
		//			erf.setParameter(UCERF2.BACK_SEIS_NAME, UCERF2.BACK_SEIS_EXCLUDE);
		//			erf.getTimeSpan().setDuration(1, TimeSpan.YEARS);
		//			erf.updateForecast();
		//			
		//			System.out.println("Setting UCERF2 rates");
		//			for (ProbEqkSource source : erf) {
		//				for (ProbEqkRupture rup : source) {
		//					if (Math.random() > 0.2d)
		//						continue;
		//					double prob = rup.getProbability();
		//					double log10prob = Math.log10(prob);
		//					if (log10prob < minX || log10prob > maxX) {
		//						System.out.println("Prob outside of bounds: "+prob + " (log10: "+log10prob+")");
		//					}
		//					int ind = (int)Math.round((log10prob-minX)/delta);
		//					ucerf2Func.set(ind, ucerf2Func.getY(ind)+1);
		//				}
		//			}
		//		}
		//		
		//		File dir = new File("/home/kevin/OpenSHA/UCERF3/test_inversion/bench/2011_09_08-morgan-CS_fixed");
		//		File binFile = new File(dir, "run1.mat");
		//		
		//		double[] rupRateSolution = MatrixIO.doubleArrayFromFile(binFile);
		//		int numNonZero = 0;
		//		for (double rate : rupRateSolution)
		//			if (rate > 0)
		//				 numNonZero++;
		//		double[] nonZeros = new double[numNonZero];
		//		int cnt = 0;
		//		for (double rate : rupRateSolution) {
		//			if (rate > 0)
		//				nonZeros[cnt++] = rate;
		//		}
		//		EvenlyDiscretizedFunc inversionFunc = new EvenlyDiscretizedFunc(minX, maxX, num);
		//		inversionFunc.setName("UCERF3 Inversion");
		//		
		//		
		//		System.out.println("Setting inversion rates");
		//		for (int i=0; i<nonZeros.length; i++) {
		//			double log10rate = Math.log10(nonZeros[i]);
		//			if (log10rate < minX || log10rate > maxX) {
		//				System.out.println("Prob outside of bounds: "+nonZeros[i] + " (log10: "+log10rate+")");
		//			}
		//			int ind = (int)Math.round((log10rate-minX)/delta);
		//			inversionFunc.set(ind, inversionFunc.getY(ind)+1);
		//		}
		//		
		//		ArrayList<DiscretizedFunc> funcs = new ArrayList<DiscretizedFunc>();
		//		funcs.add(ucerf2Func);
		//		funcs.add(inversionFunc);
		//		
		//		chars = new ArrayList<PlotCurveCharacterstics>();
		//		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
		//		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLUE));
		//		
		//		System.out.println("Displaying graph!");
		//		
		//		new GraphWindow(funcs, "Rupture Rates", chars);
		//		System.out.println("DONE");
	}

}
