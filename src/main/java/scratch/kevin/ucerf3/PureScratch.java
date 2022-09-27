package scratch.kevin.ucerf3;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.lang.reflect.Type;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.RealDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.Ellsworth_B_WG02_MagAreaRel;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.IntegerSampler.ExclusionIntegerSampler;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.IntegerPDF_FunctionSampler;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.uncertainty.BoundedUncertainty;
import org.opensha.commons.data.uncertainty.UncertainIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.Uncertainty;
import org.opensha.commons.data.uncertainty.UncertaintyBoundType;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.geo.json.FeatureCollection;
import org.opensha.commons.geo.json.FeatureProperties;
import org.opensha.commons.geo.json.GeoJSON_Type;
import org.opensha.commons.geo.json.Geometry;
import org.opensha.commons.geo.json.Geometry.GeometryCollection;
import org.opensha.commons.geo.json.Geometry.MultiPoint;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.BranchWeightProvider;
import org.opensha.commons.logicTree.BranchWeightProvider.OriginalWeights;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FaultUtils;
import org.opensha.commons.util.modules.OpenSHA_Module;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.RupSetScalingRelationship;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfigurationFactory;
import org.opensha.sha.earthquake.faultSysSolution.inversion.Inversions;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.InversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.U3MFDSubSectNuclInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.UncertainDataConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.UncertainDataConstraint.SectMappedUncertainDataConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.ConstraintRange;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.TimeCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.modules.AveSlipModule;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionMisfits;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.NamedFaults;
import org.opensha.sha.earthquake.faultSysSolution.modules.PaleoseismicConstraintData;
import org.opensha.sha.earthquake.faultSysSolution.modules.RegionsOfInterest;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupMFDsModule;
import org.opensha.sha.earthquake.faultSysSolution.modules.SectSlipRates;
import org.opensha.sha.earthquake.faultSysSolution.modules.SlipAlongRuptureModel;
import org.opensha.sha.earthquake.faultSysSolution.modules.SlipAlongRuptureModel.Tapered;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionSlipRates;
import org.opensha.sha.earthquake.faultSysSolution.modules.WaterLevelRates;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen;
import org.opensha.sha.earthquake.faultSysSolution.reports.ReportPageGen.PlotLevel;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.SectBySectDetailPlots;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.SolMFDPlot;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.GeoJSONFaultReader;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.earthquake.faultSysSolution.util.AverageSolutionCreator;
import org.opensha.sha.earthquake.faultSysSolution.util.BranchAverageSolutionCreator;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_ConstraintBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.data.NSHM23_PaleoDataLoader;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_ConstraintBuilder.ParkfieldSelectionCriteria;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.DistDependSegShift;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SingleStates;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_U3_HybridLogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.RupturePlausibilityModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SegmentationMFD_Adjustment;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.ShawSegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SupraSeisBValues;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.U3_UncertAddDeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.SeismicityRegions;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.simulators.RSQSimEvent;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.reflect.TypeToken;

import scratch.UCERF3.analysis.FaultSystemRupSetCalc;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.enumTreeBranches.TotalMag5Rate;
import scratch.UCERF3.inversion.U3InversionTargetMFDs;
import scratch.UCERF3.inversion.UCERF3InversionConfiguration;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;
import scratch.UCERF3.utils.FaultSectionDataWriter;
import scratch.UCERF3.utils.U3SectionMFD_constraint;
import scratch.UCERF3.utils.UCERF2_A_FaultMapper;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class PureScratch {
	
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
	
	private static void test112() throws IOException {
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
	
	private static void test116() throws IOException {
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(new File("/home/kevin/OpenSHA/UCERF4/rup_sets/"
				+ "nshm23_geo_dm_v1p2_all_plausibleMulti15km_adaptive6km_direct_cmlRake360_jumpP0.001"
				+ "_slipP0.05incrCapDist_cff0.75IntsPos_comb2Paths_cffFavP0.01_cffFavRatioN2P0.5_sectFractGrow0.1.zip"));
//		PaleoseismicConstraintData.loadUCERF3(rupSet);
		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
		
		LogicTreeBranch<?> branch = NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT;
		
		rupSet = factory.updateRuptureSetForBranch(rupSet, branch);
		
		InversionConfiguration config = factory.buildInversionConfig(rupSet, branch, 16);
		config = InversionConfiguration.builder(config).completion(TimeCompletionCriteria.getInMinutes(1)).build();
		
		System.out.println("Inversion Constraints");
		for (InversionConstraint constraint : config.getConstraints())
			System.out.println("\t"+constraint.getName()+", weight="+constraint.getWeight());
	
		Inversions.run(rupSet, config).write(new File("/tmp/scrach_sol.zip"));
	}
	
	private static void test117() throws IOException {
		File inputDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_01_18-reproduce-ucerf3-ref_branch-long_reweight_test-NEOK-EllB-DsrTap-SupraB0.0-"
				+ "NuclMFD-reweight_sqrt_conserve_phased-initial_parkfield/");
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(new File(inputDir, "rupture_set.zip"));
		InversionConfiguration config = InversionConfiguration.readJSON(
				new File(inputDir, "sub_0.5_per-avg_50_per/config.json"), rupSet);
		IntegerPDF_FunctionSampler sampler = (IntegerPDF_FunctionSampler)config.getSampler();
		
		InversionConfiguration.writeJSON(config, new File("/tmp/config1.json"));
		
		HashSet<Integer> excludes = new HashSet<>();
		for (int i=0; i<sampler.size(); i++)
			if (sampler.getY(i) == 0d)
				excludes.add(i);
		
		ExclusionIntegerSampler exclusionSampler = new ExclusionIntegerSampler(0, rupSet.getNumRuptures(), excludes);
		config = InversionConfiguration.builder(config).sampler(exclusionSampler).build();
		
		InversionConfiguration.writeJSON(config, new File("/tmp/config2.json"));
		
		System.out.println(InversionConfiguration.readJSON(new File("/tmp/config1.json"), rupSet).getSampler().getClass());
		System.out.println(InversionConfiguration.readJSON(new File("/tmp/config2.json"), rupSet).getSampler().getClass());
	}
	
	private static void test118() throws IOException {
		EvenlyDiscretizedFunc refMFD = new EvenlyDiscretizedFunc(0.05, 78, 0.1);
		
		double totMoRate = 1e16;
		double supraB = 0;
		
		// start with a full G-R with the supra b-value
		GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(
				refMFD.getMinX(), refMFD.size(), refMFD.getDelta(), totMoRate, supraB);
		
		System.out.println(gr.getY(1)/gr.getY(0));
		System.out.println(gr.getY(6)/gr.getY(5));
		
		double delta = refMFD.getDelta();
		double relFract = (Math.pow(10, -supraB*delta) - Math.pow(10, -supraB*delta*2))
				/(Math.pow(10, 0) - Math.pow(10, -supraB*delta));
		System.out.println(relFract);
	}
	
	private static void test119() throws IOException {
		FaultModels fm = FaultModels.FM3_1;
		for (DeformationModels dm : DeformationModels.values()) {
			if (dm.getNodeWeight(null) > 0d) {
				MinMaxAveTracker slipTrack = new MinMaxAveTracker();
				for (FaultSection sect : dm.build(fm))
					slipTrack.addValue(sect.getReducedAveSlipRate());
				System.out.println(dm.getShortName()+" slip rates: "+slipTrack);
			}
		}
	}
	
	private static void test120() throws IOException {
		// fix improperly attached modules in an inversion run
		File mainDir = new File("/data/kevin/ucerf4/batch_inversions/"
				+ "2022_02_15-coulomb-fm31-ref_branch-seg_model_adjustments-U3_ZENG-Shaw09Mod-DsrUni-SupraB0.8-TotNuclRate-ShawR0_3");
		for (File subDir : mainDir.listFiles()) {
			File solFile = new File(subDir, "solution.zip");
			if (!solFile.exists())
				continue;
			File rupSetFile = new File(subDir, "rupture_set.zip");
			FaultSystemRupSet rupSet = FaultSystemRupSet.load(rupSetFile);
			FaultSystemSolution sol = FaultSystemSolution.load(solFile);
			
			sol = sol.copy(rupSet.getArchive());
			sol.write(solFile);
		}
	}
	
	private static void test121() throws IOException {
		// hack to replace target mfds, fake solution
		File mainDir = new File("/data/kevin/ucerf4/batch_inversions/"
				+ "2022_02_15-coulomb-fm31-ref_branch-seg_model_adjustments-U3_ZENG-Shaw09Mod-DsrUni-SupraB0.8-TotNuclRate-ShawR0_3");
		File inputSolFile = new File(mainDir, "JumpProb/solution.zip");
		File inputRupSetFile = new File(mainDir, "StictEquivJumpProb/rupture_set.zip");
		File destSolFile = new File(mainDir, "StictEquivJumpProb/solution.zip");
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(inputRupSetFile);
		FaultSystemSolution sol = FaultSystemSolution.load(inputSolFile);
		sol.getRupSet().addModule(rupSet.requireModule(InversionTargetMFDs.class));
		sol.write(destSolFile);
	}
	
	private static void test122() throws IOException {
		File inputSol = new File("/tmp/bad_inv/solution.zip");
		InversionConfiguration config = FaultSystemSolution.load(inputSol).requireModule(InversionConfiguration.class);
		InversionConfiguration.writeJSON(config, new File(inputSol.getParentFile(), "config.json"));
	}
	
	private static void test123() throws IOException, InterruptedException, ExecutionException {
		File inputDir = new File("/tmp/avg_debug");
		
		List<FaultSystemSolution> sols = new ArrayList<>();
		for (int i=0; i<10; i++) {
			File solFile = new File(inputDir, "solution_"+i+".zip");
			FaultSystemSolution sol = FaultSystemSolution.load(solFile);
			
			sols.add(sol);
		}
		
		List<List<ConstraintRange>> solRanges = new ArrayList<>();
		for (int i=0; i<sols.size(); i++) {
			FaultSystemSolution sol = sols.get(i);
			InversionMisfits misfits = sol.requireModule(InversionMisfits.class);
			double[] eqMisfits = misfits.getMisfits();
			double[] ineqMisfits = misfits.getInequalityMisfits();
			System.out.println("Solution "+i+" has "+eqMisfits.length+" equality misfits");
			System.out.println("Solution "+i+" has "+ineqMisfits.length+" inequality misfits");
			solRanges.add(misfits.getConstraintRanges());
		}

		System.out.println("Constraint Ranges 0");
		for (ConstraintRange range : solRanges.get(0))
			System.out.println("\t"+range);
		System.out.println("Constraint Ranges 9");
		for (ConstraintRange range : solRanges.get(9))
			System.out.println("\t"+range);
		
////		int expected = getNumRows(FaultSystemRupSetCalc.getCharInversionSectMFD_Constraints(sols.get(0).getRupSet()));
////		System.out.println("Reproduce with RS 0: "+expected+" rows");
////		System.out.println("Reproduce with RS 9: "+getNumRows(
////				FaultSystemRupSetCalc.getCharInversionSectMFD_Constraints(sols.get(9).getRupSet()))+" rows");
//		// do it in parallel to try to catch a synchronization error
//		List<Future<Integer>> futures = new ArrayList<>();
//		ExecutorService exec = Executors.newFixedThreadPool(50);
////		ExecutorService exec = Executors.newFixedThreadPool(1);
//		for (int i=0; i<100; i++) {
//			FaultSystemRupSet rupSet = sols.get(i % sols.size()).getRupSet();
//			futures.add(exec.submit(new Callable<Integer>() {
//
//				@Override
//				public Integer call() throws Exception {
//					return getNumRows(FaultSystemRupSetCalc.getCharInversionSectMFD_Constraints(rupSet));
//				}
//			}));
//		}
//		
//		Map<Integer, Integer> rowCountCounts = new HashMap<>();
//		for (Future<Integer> future : futures) {
//			int numRows = future.get();
//			
//			Integer prevCount = rowCountCounts.get(numRows);
//			if (prevCount == null)
//				prevCount = 0;
//			rowCountCounts.put(numRows, prevCount+1);
//		}
//		for (int numRows : rowCountCounts.keySet()) {
//			System.out.println(numRows+": "+rowCountCounts.get(numRows)+" times");
//		}
//		exec.shutdown();
		
		InversionConfiguration config0 = sols.get(0).requireModule(InversionConfiguration.class);
		InversionConfiguration config9 = sols.get(9).requireModule(InversionConfiguration.class);
		
		U3MFDSubSectNuclInversionConstraint constr0 = null;
		for (InversionConstraint constr : config0.getConstraints())
			if (constr instanceof U3MFDSubSectNuclInversionConstraint)
				constr0 = (U3MFDSubSectNuclInversionConstraint)constr;
		U3MFDSubSectNuclInversionConstraint constr9 = null;
		for (InversionConstraint constr : config9.getConstraints())
			if (constr instanceof U3MFDSubSectNuclInversionConstraint)
				constr9 = (U3MFDSubSectNuclInversionConstraint)constr;
		
		System.out.println("Rows0: "+constr0.getNumRows());
		System.out.println("Rows1: "+constr9.getNumRows());
		
		List<U3SectionMFD_constraint> mfds0 = constr0.getConstraints();
		List<U3SectionMFD_constraint> mfds1 = constr9.getConstraints();
		
		Preconditions.checkState(mfds0.size() == mfds1.size());
		for (int i=0; i<mfds0.size(); i++) {
			FaultSection sect = sols.get(0).getRupSet().getFaultSectionData(i);
			String str = i+". "+sect.getName();
			if (UCERF2_A_FaultMapper.wasUCERF2_TypeAFault(sect.getParentSectionId()))
				str += " (type A)";
			
			U3SectionMFD_constraint mfd0 = mfds0.get(i);
			U3SectionMFD_constraint mfd1 = mfds1.get(i);
			
			if ((mfd0 == null) != (mfd1 == null)) {
				System.out.println("Null mismatch for "+str+": mfd0 ? "+(mfd0 == null)+", mfd1 ? "+(mfd1 == null));
			} else if (mfd0 != null && mfd1 != null) {
				// see if they're the same
				boolean valsSame = true;
				if (mfd0.getNumMags() != mfd1.getNumMags()) {
					System.out.println("Mag count mismatch for "+str+": mfd0="+mfd0.getNumMags()+", mfd1="+mfd1.getNumMags());
				} else {
					for (int j=0; j<mfd0.getNumMags(); j++) {
						double rate0 = mfd0.getRate(j);
						double rate1 = mfd1.getRate(j);
						if ((rate0 > 0) != (rate1 > 0)) {
							System.out.println(">0 mismatch for "+str+", bin "+j+" ("+mfd0.getMag(j)+"): "+rate0+" != "+rate1);
						}
						valsSame = valsSame && (float)rate0 == (float)rate1;
					}
				}
				if (!valsSame)
					System.out.println("Value(s) mismatch for "+str);
			}
		}
		
//		AverageSolutionCreator.buildAverage(sols.toArray(new FaultSystemSolution[0]));
	}
	
	private static int getNumRows(List<U3SectionMFD_constraint> constraints) {
		int numRows = 0;
		for (U3SectionMFD_constraint constraint : constraints)
			if (constraint != null)
				for (int i=0; i<constraint.getNumMags(); i++)
					if (constraint.getRate(i) > 0)
						numRows++;
		return numRows;
	}
	
	private static void test124() throws IOException, InterruptedException, ExecutionException {
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(new File("/tmp/avg_debug/solution_9.zip"));
		
////		int expected = getNumRows(FaultSystemRupSetCalc.getCharInversionSectMFD_Constraints(sols.get(0).getRupSet()));
////		System.out.println("Reproduce with RS 0: "+expected+" rows");
////		System.out.println("Reproduce with RS 9: "+getNumRows(
////				FaultSystemRupSetCalc.getCharInversionSectMFD_Constraints(sols.get(9).getRupSet()))+" rows");
		// do it in parallel to try to catch a synchronization error
		List<Future<Integer>> futures = new ArrayList<>();
		ExecutorService exec = Executors.newFixedThreadPool(16);
//		ExecutorService exec = Executors.newFixedThreadPool(1);
		for (int i=0; i<10000; i++) {
			futures.add(exec.submit(new Callable<Integer>() {

				@Override
				public Integer call() throws Exception {
					try {
						UCERF3InversionConfiguration.getUCERF2MagsAndrates(rupSet, FaultModels.FM3_1);
						return 0;
					} catch (Exception e) {
						e.printStackTrace();
						System.err.flush();
						System.exit(1);
						return -1;
					}
				}
			}));
		}
		
		Map<Integer, Integer> rowCountCounts = new HashMap<>();
		for (Future<Integer> future : futures) {
			int numRows = future.get();
			
			Integer prevCount = rowCountCounts.get(numRows);
			if (prevCount == null)
				prevCount = 0;
			rowCountCounts.put(numRows, prevCount+1);
		}
		for (int numRows : rowCountCounts.keySet()) {
			System.out.println(numRows+": "+rowCountCounts.get(numRows)+" times");
		}
		exec.shutdown();
		
//		AverageSolutionCreator.buildAverage(sols.toArray(new FaultSystemSolution[0]));
	}
	
	private static void test125() {
		GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(6.05, 20, 0.1, 1e16, 0.8);
		System.out.println(gr);
		
		double sectArea = 7d*14d;
		Ellsworth_B_WG02_MagAreaRel ma = new Ellsworth_B_WG02_MagAreaRel();
		
		IncrementalMagFreqDist particMFD = new IncrementalMagFreqDist(gr.getMinX(), gr.size(), gr.getDelta());
		for (int i=0; i<gr.size(); i++) {
			double rupArea = ma.getMedianArea(gr.getX(i));
			particMFD.set(i, gr.getY(i)*rupArea/sectArea);
		}
		
		System.out.println(particMFD);
		
		double totParticRate = particMFD.calcSumOfY_Vals();
		double particRateAbove7p5 = 0d;
		for (int i=0; i<gr.size(); i++)
			if (particMFD.getX(i) >= 7.5)
				particRateAbove7p5 += particMFD.getY(i);
		
		System.out.println("Fract rate >7.5 = "+particRateAbove7p5+" / "+totParticRate+" = "+(particRateAbove7p5/totParticRate));
	}
	
	private static void test126() throws IOException {
//		List<? extends FaultSection> sects = NSHM23_DeformationModels.GEOL_V1p3.buildGeolFullSects(
//				NSHM23_FaultModels.NSHM23_v1p4, "v1p3");
//		
//		GeoJSONFaultReader.writeFaultSections(new File("/tmp/nshm23_dm1p3_sects.geojson"), sects);
	}
	
	private static void test127() throws IOException {
		List<? extends FaultSection> subSects = NSHM23_DeformationModels.GEOLOGIC.build(NSHM23_FaultModels.NSHM23_v1p4);
		
		FaultSectionDataWriter.writeSectionsToFile(subSects, null, new File("/tmp/nshm23_dm1p3_sub_sects.txt"), false);
		GeoJSONFaultReader.writeFaultSections(new File("/tmp/nshm23_dm1p3_sub_sects.geojson"), subSects);
	}
	
	private static void test128() throws IOException {
		double minMag = Double.POSITIVE_INFINITY;
//		for (RSQSimEvent e : Catalogs.BRUCE_2585_1MYR.instance().loader().skipYears(20000).iterable()) {
		for (RSQSimEvent e : Catalogs.BRUCE_5295.instance().loader().skipYears(5000).iterable()) {
			double mag = e.getMagnitude();
			if (mag < minMag && Double.isFinite(mag)) {
				System.out.println("New min mag: "+(float)mag);
				minMag = mag;
			}
		}
		System.exit(0);
	}
	
	private static void test129() throws IOException {
		WC1994_MagLengthRelationship wc94 = new WC1994_MagLengthRelationship();
		
		wc94.setRake(-90);
		System.out.println(wc94.getMedianMag(45.33335));
		System.out.println(wc94.getMedianMag(39.34729));
	}
	
	private static void test130() throws IOException {
		double moment = MagUtils.magToMoment(7.2);
		double len = 64.29441;
		double width = 15d/Math.sin(45d*Math.PI/ 180);
		double area = len*width;
		System.out.println("Area = "+(float)len+" x "+(float)width+" = "+(float)area);
		area *= 1e6; // km^2 -> m^2
		double slip = FaultMomentCalc.getSlip(area, moment);
		System.out.println("Slip: "+(float)slip);
		double slipRate = slip * 1e-7;
		System.out.println("Slip rate: "+(float)slipRate+" m/yr = "+(float)(slipRate*1e3)+" mm/yr");
	}
	
	private static void test131() throws IOException {
		FeatureCollection features = FeatureCollection.read(new File("/tmp/nshm18_test/fault_sections.geojson"));
		List<Feature> modFeatures = new ArrayList<>();
		
		int index = 0;
		for (int i=803; i<features.features.size(); i++) {
			Feature feature = features.features.get(i);
			Feature modFeature = new Feature(index, feature.geometry, feature.properties);
			modFeature.properties.set(GeoJSONFaultSection.FAULT_ID, index);
			index++;
			modFeatures.add(modFeature);
		}
		
		FeatureCollection.write(new FeatureCollection(modFeatures), new File("/tmp/nshm18_test/second_half.geojson"));
	}
	
	private static void test132() throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/markdown/inversions/2022_05_06-nshm23-cluster_specific-u3_reduced/solution.zip"));
		sol.write(new File("/tmp/asdf.zip"));
	}
	
	private static void test133() throws IOException {
		FaultSystemSolution sol1 = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_05_09-nshm23_u3_hybrid_branches-cluster_specific_inversion-shift_seg_1km-FM3_1-CoulombRupSet"
				+ "-DsrUni-NuclMFD-SubB1-ThreshAvg/results_FM3_1_CoulombRupSet_branch_averaged.zip"));
		FaultSystemSolution sol2 = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_05_09-nshm23_u3_hybrid_branches-cluster_specific_inversion-shift_seg_1km-FM3_1-CoulombRupSet"
				+ "-DsrUni-TotNuclRate-SubB1-ThreshAvg/results_FM3_1_CoulombRupSet_branch_averaged.zip"));
		
		IncrementalMagFreqDist mfd1 = sol1.getRupSet().requireModule(InversionTargetMFDs.class).getTotalOnFaultSupraSeisMFD();
		IncrementalMagFreqDist mfd2 = sol2.getRupSet().requireModule(InversionTargetMFDs.class).getTotalOnFaultSupraSeisMFD();
		
		System.out.println("Moment rates: "+mfd1.getTotalMomentRate()+"\t"+mfd2.getTotalMomentRate());
		System.out.println("Incr rates: "+mfd1.getTotalIncrRate()+"\t"+mfd2.getTotalIncrRate());
	}
	
	private static void test134() {
		// From Baken et al (2004), Figure S2
		List<int[]> parfieldYearsList = new ArrayList<>();
		List<String> labels = new ArrayList<>();
		
		labels.add("Basic Sequence");
		parfieldYearsList.add(new int[] { 1857, 1881, 1901, 1922, 1934, 1966, 2004 });
		
		labels.add("Alt 1");
		parfieldYearsList.add(new int[] { 1857, 1877, 1881, 1901, 1922, 1934, 1966, 2004 });
		
		labels.add("Alt 2");
		parfieldYearsList.add(new int[] { 1881, 1901, 1922, 1934, 1966, 2004 });
		
		System.out.println("Bakun et al. (2004) Figure S2:");
		
		DecimalFormat yrDF = new DecimalFormat("0.0");
		DecimalFormat pDF = new DecimalFormat("0.0%");
		
		for (int s=0; s<labels.size(); s++) {
			String label = labels.get(s);
			int[] parkfieldYears = parfieldYearsList.get(s);
			int[] parkfieldRIs = new int[parkfieldYears.length-1];
			for (int i=1; i<parkfieldYears.length; i++)
				parkfieldRIs[i-1] = parkfieldYears[i] - parkfieldYears[i-1];
			
			double[] risDouble = new double[parkfieldRIs.length];
			for (int i=0; i<risDouble.length; i++)
				risDouble[i] = parkfieldRIs[i];
			
			double mean = StatUtils.mean(risDouble);
			double sd = Math.sqrt(StatUtils.variance(risDouble));
			double sdom = sd/Math.sqrt(parkfieldRIs.length);
			double fractSDOM = sdom/mean;
			System.out.println("\n"+label);
			System.out.println("Parkfield event years: "+Joiner.on(",").join(Ints.asList(parkfieldYears)));
			System.out.println("Parkfield event RIs: "+Joiner.on(",").join(Ints.asList(parkfieldRIs)));
			System.out.println("mean="+yrDF.format(mean)+"\tSD="+yrDF.format(sd)+"\tSDOM="+yrDF.format(sdom)+" ("+pDF.format(fractSDOM)+")");
			
			RealDistribution riMeanDist = new NormalDistribution(mean, sdom);
			int trials = 10000000;
			double[] testRates = new double[trials];
			for (int i=0; i<trials; i++)
				testRates[i] = 1d/Math.max(1d, riMeanDist.sample());
			double testSD = Math.sqrt(StatUtils.variance(testRates));
			double rate = 1d/mean;
			System.out.println("Rate SDOM from normal dist: "+(float)testSD+" ("+pDF.format(testSD/rate)+")");
			
			double shape = sdom/mean;
			System.out.println("LN shape: "+shape);
			riMeanDist = new LogNormalDistribution(Math.log(mean), shape);
			testRates = new double[trials];
			for (int i=0; i<trials; i++)
				testRates[i] = 1d/Math.max(1d, riMeanDist.sample());
			testSD = Math.sqrt(StatUtils.variance(testRates));
			rate = 1d/mean;
			System.out.println("Rate SDOM from log-normal dist: "+(float)testSD+" ("+pDF.format(testSD/rate)+")");
		}
	}
	
	private static void test135() {
		double rateM5 = TotalMag5Rate.RATE_7p9.getRateMag5();
		
		GutenbergRichterMagFreqDist u3Target = new GutenbergRichterMagFreqDist(U3InversionTargetMFDs.MIN_MAG,
				U3InversionTargetMFDs.NUM_MAG, U3InversionTargetMFDs.DELTA_MAG);
		double roundedMmaxOnFault = 8.55d;
		u3Target.setAllButTotMoRate(U3InversionTargetMFDs.MIN_MAG, roundedMmaxOnFault, rateM5*1e5, 1.0);
		
		double u3CumRate = u3Target.getCumRateDistWithOffset().getY(5d);
		
		System.out.println("Target rate M>=5: "+rateM5);
		System.out.println("U3 MFD rate M>=5: "+u3CumRate);
		
		GutenbergRichterMagFreqDist testTarget = new GutenbergRichterMagFreqDist(1.0, TotalMag5Rate.RATE_7p9.getRateMag5(), 5.05, 9.95, 50);
		double testCumRate = testTarget.getCumRateDistWithOffset().getY(5d);
		
		System.out.println("Test MFD rate M>=5: "+testCumRate);
		
		testTarget = new GutenbergRichterMagFreqDist(5.05, 9.95, 50);
		testTarget.setAllButTotMoRate(5.05, roundedMmaxOnFault, rateM5, 1d);
		 testCumRate = testTarget.getCumRateDistWithOffset().getY(5d);
		
		System.out.println("Test2 MFD rate M>=5: "+testCumRate);
	}
	
	private static void test136() throws IOException {
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(
//				new File("/data/kevin/markdown/inversions/fm3_1_u3ref_uniform_reproduce_ucerf3.zip"));
				new File("/data/kevin/markdown/inversions/fm3_1_u3ref_uniform_coulomb.zip"));
		
		rupSet = FaultSystemRupSet.buildFromExisting(rupSet)
//				.replaceFaultSections(DeformationModels.MEAN_UCERF3.build(FaultModels.FM3_1))
				.replaceFaultSections(U3_UncertAddDeformationModels.U3_MEAN.build(FaultModels.FM3_1))
				.forScalingRelationship(ScalingRelationships.MEAN_UCERF3)
				.build();
		
		File outputDir = new File("/tmp/test_sect_by_sect");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		SectBySectDetailPlots page = new SectBySectDetailPlots();
		page.writePlot(rupSet, null, "Rupture Set", outputDir);
	}
	
	private static void test137() throws IOException {
		LogicTreeBranch<LogicTreeNode> branch = NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT;
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_05_27-nshm23_branches-NSHM23_v1p4-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvg/"
				+ "results_NSHM23_v1p4_CoulombRupSet_branch_averaged.zip"));
		
//		FaultSystemRupSet rupSet = FaultSystemRupSet.load(
////				new File("/data/kevin/markdown/inversions/fm3_1_u3ref_uniform_reproduce_ucerf3.zip"));
//				new File("/data/kevin/markdown/inversions/fm3_1_u3ref_uniform_coulomb.zip"));
////		LogicTreeBranch<LogicTreeNode> branch = NSHM23_U3_HybridLogicTreeBranch.DEFAULT.copy();
//		
//		branch.setValue(U3_UncertAddDeformationModels.U3_ABM);
//		branch.setValue(ScalingRelationships.ELLSWORTH_B);
//		branch.setValue(SupraSeisBValues.B_0p0);
//		branch.setValue(SegmentationModels.SHAW_R0_3);
//		branch.setValue(DistDependSegShift.TWO_KM);
//		branch.setValue(SegmentationMFD_Adjustment.REL_GR_THRESHOLD_AVG_ITERATIVE);

//		File outputFile = new File("/tmp/nshm_seg_wt_100.zip");
//		InversionConfigurationFactory factory = new NSHM23_InvConfigFactory.ClusterSpecificSegWeight100();
		File outputFile = new File("/tmp/nshm_seg_wt_10.zip");
		InversionConfigurationFactory factory = new NSHM23_InvConfigFactory();
		
		rupSet = factory.updateRuptureSetForBranch(rupSet, branch);
		
		FaultSystemSolution fss = Inversions.run(rupSet, factory, branch, 16, null);
		
		fss.write(outputFile);
	}
	
	private static void test138() throws IOException {
		SolutionLogicTree.load(new File("/data/kevin/ucerf4/batch_inversions/2022_06_10-nshm23_u3_hybrid_branches-FM3_1-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvgIterRelGR-IncludeThruCreep/results.zip"));
	}
	
	private static void test139() throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/tmp/FM3_1_branch_averaged.zip"));
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		CSVFile<String> taperedCSVs = new CSVFile<>(false);
		taperedCSVs.addLine("Rupture Index", "Section 1 Slip (m)", "...");
		
		AveSlipModule aveSlips = rupSet.requireModule(AveSlipModule.class);
		
		Tapered dsr = new SlipAlongRuptureModel.Tapered();
		
		for (int r=0; r<rupSet.getNumRuptures(); r++) {
			List<String> line = new ArrayList<>();
			line.add(r+"");
			double[] slips = dsr.calcSlipOnSectionsForRup(rupSet, aveSlips, r);
			for (double slip : slips)
				line.add((float)slip+"");
			taperedCSVs.addLine(line);
		}
		
		taperedCSVs.writeToFile(new File("/tmp/u3_fm3_1_tapered_dsr.csv"));
	}
	
	private static void test140() throws IOException {
//		FeatureProperties props = new FeatureProperties();
//		Feature feature = new Feature(new Geometry.Point(new Location(34, -118)), props);
//		
//		XY_DataSet xy1 = new DefaultXY_DataSet();
//		xy1.set(3d, 5d);
//		xy1.set(4d, 76d);
//		
//		props.set("TestFunc1", xy1);
//		props.set("Color", Color.BLUE);
//		IncrementalMagFreqDist mfd = new GutenbergRichterMagFreqDist(1d, 5d, 5d, 10d, 10);
//		props.set("MFD", mfd);
//		props.set("Name", "My feature");
//		
//		String json = feature.toJSON();
//		System.out.println(json);
//		
//		Feature feature2 = Feature.fromJSON(json);
		Feature feature2 = Feature.read(new File("/tmp/test.json"));
		
		for (String key : feature2.properties.keySet()) {
			Object obj = feature2.properties.get(key);
			if (obj == null)
				System.out.println(key+": "+null);
			else
				System.out.println(key+": "+obj.getClass());
		}
	}
	
	private static void test141() throws IOException {
		File runDir = new File("/data/kevin/nshm23/batch_inversions/"
				+ "2022_07_21-nshm23_branches-NSHM23_v1p4-CoulombRupSet-DsrUni-TotNuclRate-SubB1-ShawR0_3-Shift2km-"
				+ "ThreshAvgIterRelGR-IncludeThruCreep");
		LogicTree<?> tree = LogicTree.read(new File(runDir, "logic_tree.json"));
		
		BranchAverageSolutionCreator baCreator = new BranchAverageSolutionCreator(new BranchWeightProvider.OriginalWeights());
		
		NSHM23_FaultModels fm = NSHM23_FaultModels.NSHM23_v1p4;
//		NSHM23_DeformationModels dm = NSHM23_DeformationModels.GEOLOGIC;
		NSHM23_DeformationModels dm = NSHM23_DeformationModels.EVANS;
		
		File resultsDir = new File(runDir, "results");
		
		for (LogicTreeBranch<?> branch : tree) {
			if (branch.hasValue(fm) && branch.hasValue(dm)) {
				File branchDir = new File(resultsDir, branch.buildFileName());
				Preconditions.checkState(branchDir.exists(), "Branch dir doesn't exist: %s", branchDir.getName());
				File solFile = new File(branchDir, "solution.zip");
				Preconditions.checkState(branchDir.exists(), "Branch solution doesn't exist: %s", solFile.getAbsolutePath());
				FaultSystemSolution sol = FaultSystemSolution.load(solFile);
				// re-attach modules
				sol.getRupSet().removeModuleInstances(RegionsOfInterest.class);
				fm.attachDefaultModules(sol.getRupSet());
				
				baCreator.addSolution(sol, branch);
			}
		}
		
		FaultSystemSolution baSol = baCreator.build();
		
		baSol.write(new File(runDir, "results_"+fm.getFilePrefix()+"_CoulombRupSet_"+dm.getFilePrefix()+"_branch_averaged.zip"));
	}
	
	private static void test142() throws IOException {
		File runDir = new File("/project/scec_608/kmilner/nshm23/batch_inversions/"
				+ "2022_07_23-nshm23_branches-NSHM23_v1p4-CoulombRupSet-DsrUni-TotNuclRate-SubB1-Shift2km-ThreshAvgIterRelGR-IncludeThruCreep");
		LogicTree<?> tree = LogicTree.read(new File(runDir, "logic_tree.json"));
		
		File resultsDir = new File(runDir, "results");
		
		for (LogicTreeBranch<?> branch : tree) {
			File branchDir = new File(resultsDir, branch.buildFileName());
			Preconditions.checkState(branchDir.exists(), "Branch dir doesn't exist: %s", branchDir.getName());
			File solFile = new File(branchDir, "solution.zip");
			Preconditions.checkState(branchDir.exists(), "Branch solution doesn't exist: %s", solFile.getAbsolutePath());
			FaultSystemSolution sol = FaultSystemSolution.load(solFile);
			
			InversionTargetMFDs targetMFDs = sol.getRupSet().requireModule(InversionTargetMFDs.class);
			
			List<? extends IncrementalMagFreqDist> supraSeisNuclMFDs = targetMFDs.getOnFaultSupraSeisNucleationMFDs();
			for (int i=0; i<supraSeisNuclMFDs.size(); i++) {
				IncrementalMagFreqDist supraMFD = supraSeisNuclMFDs.get(i);
				if (!(supraMFD instanceof UncertainIncrMagFreqDist)) {
					System.out.println(branch);
					System.out.println("Section "+i+" is of type "+(supraMFD == null ? "null" : supraMFD.getClass()));
					System.out.println(supraMFD);
					System.out.flush();
					System.exit(1);
				}
			}
		}
	}
	
	private static void test143() throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/tmp/solution.zip"));
		
		InversionTargetMFDs mfds = sol.getRupSet().requireModule(InversionTargetMFDs.class);
		
		List<? extends IncrementalMagFreqDist> supraMFDs = mfds.getOnFaultSupraSeisNucleationMFDs();
		
		int sectIndex = 3707;
		IncrementalMagFreqDist mfd = supraMFDs.get(sectIndex);
		System.out.println("MFD type: "+mfd.getClass());
		System.out.println(mfd);
		SectSlipRates slipRates = sol.getRupSet().getSectSlipRates();
		System.out.println("Target for "+sol.getRupSet().getFaultSectionData(sectIndex).getSectionName()+": "
				+slipRates.getSlipRate(sectIndex)+" +/- "+slipRates.getSlipRateStdDev(sectIndex));
	}
	
	private static void test144() throws IOException {
		File runDir = new File("/project/scec_608/kmilner/nshm23/batch_inversions/"
				+ "2022_07_25-nshm23_branches-NSHM23_v1p4-CoulombRupSet-DsrUni-TotNuclRate-SubB1-ThreshAvgIterRelGR-IncludeThruCreep");
		LogicTree<?> tree = LogicTree.read(new File(runDir, "logic_tree.json"));
		
		File outputDir = new File(runDir, "partial_bas");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		NSHM23_DeformationModels[] restrictDMs = { NSHM23_DeformationModels.EVANS,
				NSHM23_DeformationModels.GEOLOGIC, NSHM23_DeformationModels.POLLITZ };
		
		OriginalWeights weightProv = new BranchWeightProvider.OriginalWeights();
		Map<NSHM23_SegmentationModels, BranchAverageSolutionCreator> segBACreators = new HashMap<>();
		BranchAverageSolutionCreator fullCreator = new BranchAverageSolutionCreator(weightProv);
		
		NSHM23_FaultModels fm = NSHM23_FaultModels.NSHM23_v1p4;
		
		File resultsDir = new File(runDir, "results");
		
		for (LogicTreeBranch<?> branch : tree) {
			boolean hasDM = false;
			for (NSHM23_DeformationModels dm : restrictDMs) {
				if (branch.hasValue(dm)) {
					hasDM = true;
					break;
				}
			}
			if (branch.hasValue(fm) && hasDM) {
				File branchDir = new File(resultsDir, branch.buildFileName());
				Preconditions.checkState(branchDir.exists(), "Branch dir doesn't exist: %s", branchDir.getName());
				File solFile = new File(branchDir, "solution.zip");
				Preconditions.checkState(branchDir.exists(), "Branch solution doesn't exist: %s", solFile.getAbsolutePath());
				FaultSystemSolution sol = FaultSystemSolution.load(solFile);
				// re-attach modules
				sol.getRupSet().removeModuleInstances(RegionsOfInterest.class);
				fm.attachDefaultModules(sol.getRupSet());
				
				fullCreator.addSolution(sol, branch);
				NSHM23_SegmentationModels segModel = branch.requireValue(NSHM23_SegmentationModels.class);
				BranchAverageSolutionCreator segBACreator = segBACreators.get(segModel);
				if (segBACreator == null) {
					segBACreator = new BranchAverageSolutionCreator(weightProv);
					segBACreators.put(segModel, segBACreator);
				}
				segBACreator.addSolution(sol, branch);
			}
		}
		
		FaultSystemSolution baSol = fullCreator.build();
		
		baSol.write(new File(outputDir, "results_"+fm.getFilePrefix()+"_CoulombRupSet_branch_averaged.zip"));
		
		for (NSHM23_SegmentationModels segModel : segBACreators.keySet()) {
			FaultSystemSolution segSol = segBACreators.get(segModel).build();
			segSol.write(new File(outputDir, "results_"+fm.getFilePrefix()+"_CoulombRupSet_"+segModel.getFilePrefix()+"_branch_averaged.zip"));
		}
	}
	
	private static void test145() throws IOException {
		int sectIndex = 4412;
//		int sectIndex = 4411;
		
		File baDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_07_25-nshm23_branches-NSHM23_v1p4-CoulombRupSet-DsrUni-TotNuclRate-SubB1-ThreshAvgIterRelGR"
				+ "-IncludeThruCreep/node_branch_averaged");
		
		List<FaultSystemSolution> sols = new ArrayList<>();
		List<String> names = new ArrayList<>();
		List<Color> colors = new ArrayList<>();
		
		sols.add(FaultSystemSolution.load(new File(baDir, "SegModel_LowSeg.zip")));
		names.add("Low-Seg");
		colors.add(Color.GREEN.darker());
		
		sols.add(FaultSystemSolution.load(new File(baDir, "SegModel_MidSeg.zip")));
		names.add("Mid-Seg");
		colors.add(Color.BLUE.darker());
		
		sols.add(FaultSystemSolution.load(new File(baDir, "SegModel_HighSeg.zip")));
		names.add("High-Seg");
		colors.add(Color.RED.darker());
		
		List<IncrementalMagFreqDist> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		for (int i=0; i<sols.size(); i++) {
			FaultSystemSolution sol = sols.get(i);
			String name = names.get(i);
			Color color = colors.get(i);
			
			InversionTargetMFDs targets = sol.getRupSet().getModule(InversionTargetMFDs.class);
			
			EvenlyDiscretizedFunc refMFD;
			if (targets != null) {
				IncrementalMagFreqDist target = targets.getOnFaultSupraSeisNucleationMFDs().get(sectIndex);
				refMFD = target;
				target.setName(name+" Target");
				
				funcs.add(target);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, color));
			} else {
				refMFD = SupraSeisBValInversionTargetMFDs.buildRefXValues(sol.getRupSet());
			}
			IncrementalMagFreqDist nucl = sol.calcNucleationMFD_forSect(
					sectIndex, refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
			nucl.setName(name+" Nucl");
			IncrementalMagFreqDist partic = sol.calcParticipationMFD_forSect(
					sectIndex, refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
			partic.setName(name+" Partic");
			
			funcs.add(nucl);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, color));
			
			funcs.add(partic);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, color));
			
			if (sol.hasModule(SolutionSlipRates.class)) {
				SolutionSlipRates solSlips = sol.getModule(SolutionSlipRates.class);
				SectSlipRates targetSlips = sol.getRupSet().getModule(SectSlipRates.class);
				
				double solSlip = solSlips.get(sectIndex)*1e3;
				double targetSlip = targetSlips.getSlipRate(sectIndex)*1e3;
				double targetSD = targetSlips.getSlipRateStdDev(sectIndex)*1e3;
				double z = (solSlip - targetSlip)/targetSD;
				
				System.out.println(name+" Slip Rates");
				System.out.println("\tTarget: "+(float)targetSlip+" +/- "+(float)targetSD);
				System.out.println("\tSolution: "+(float)solSlip+"\tz="+(float)z);
			}
		}
		
		String title = sectIndex+". "+sols.get(0).getRupSet().getFaultSectionData(sectIndex).getSectionName();
		
		for (boolean cumulative : new boolean[] {false, true}) {
			List<DiscretizedFunc> myFuncs = new ArrayList<>();
			if (cumulative) {
				for (IncrementalMagFreqDist func : funcs)
					myFuncs.add(func.getCumRateDistWithOffset());
			} else {
				myFuncs.addAll(funcs);
			}
			PlotSpec spec = new PlotSpec(myFuncs, chars, title, "Magnitude",
					cumulative ? "Cumulative Rate" : "Incremental Rate");
			spec.setLegendVisible(true);
			GraphWindow gw = new GraphWindow(spec);
			gw.setYLog(true);
			gw.setAxisRange(6d, 7.3, 1e-7, cumulative ? 1e-3 : 1e-4);
			gw.setVisible(true);
			gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
		}
	}
	
	private static void test146() throws IOException {
		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
		
		LogicTreeBranch<LogicTreeNode> branch = NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT.copy();
		branch.setValue(RupturePlausibilityModels.AZIMUTHAL_REDUCED);
		
		int threads = 8;
		
		FaultSystemRupSet rupSet = factory.buildRuptureSet(branch, threads);
		
		File rupSetFile = new File("/tmp/rupture_set.zip");
		rupSet.write(rupSetFile);
		
		rupSet = FaultSystemRupSet.load(rupSetFile);
		PaleoseismicConstraintData data = rupSet.requireModule(PaleoseismicConstraintData.class);
		
		System.out.println("Prob model: "+data.getPaleoProbModel().getClass());
		
		InversionConfiguration config = factory.buildInversionConfig(rupSet, branch, threads);
		config = InversionConfiguration.builder(config).avgThreads(4, TimeCompletionCriteria.getInMinutes(1))
				.completion(TimeCompletionCriteria.getInMinutes(10)).build();
		
		FaultSystemSolution sol = Inversions.run(rupSet, config);
		
		sol.write(new File("/tmp/solution.zip"));
		
		ReportPageGen report = new ReportPageGen(rupSet, sol, "Test Solution", new File("/tmp/report"),
				ReportPageGen.getDefaultSolutionPlots(PlotLevel.DEFAULT));
		report.setReplot(true);
		
		report.generatePage();
	}
	
	private static void test147() throws IOException {
		File rupSetFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_07_28-nshm23_branches-NSHM23_v1p4-CoulombRupSet-NSHM23_Avg-DsrUni-TotNuclRate-SubB1-"
				+ "ThreshAvgIterRelGR-IncludeThruCreep/results_NSHM23_v1p4_CoulombRupSet_branch_averaged.zip");
		RupSetScalingRelationship[] scales = NSHM23_ScalingRelationships.values();
		
//		File rupSetFile = new File("/home/kevin/OpenSHA/UCERF4/rup_sets/fm3_1_ucerf3.zip");
//		RupSetScalingRelationship[] scales = ScalingRelationships.values();
		
		FaultSystemRupSet rawRupSet = FaultSystemRupSet.load(rupSetFile);
		
		for (RupSetScalingRelationship scale : scales) {
			if (scale.getNodeWeight(null) == 0d)
				continue;
			FaultSystemRupSet rupSet = FaultSystemRupSet.buildFromExisting(rawRupSet).forScalingRelationship(scale).build();
			
			System.out.println("***********************");
			System.out.println("Scaling relationship: "+scale.getName());
			System.out.println("***********************");
			for (ParkfieldSelectionCriteria criteria : ParkfieldSelectionCriteria.values()) {
				System.out.println("Criteria: "+criteria);
				List<Integer> parkfieldRups = NSHM23_ConstraintBuilder.findParkfieldRups(rupSet, criteria);
				
				System.out.println("Found "+parkfieldRups.size()+" ruptures");
				for (int rupIndex : parkfieldRups) {
					double mag = rupSet.getMagForRup(rupIndex);
					System.out.println("\t"+rupIndex+". M="+(float)mag
							+",\t"+rupSet.getSectionsIndicesForRup(rupIndex).size()+" sects");
				}
				System.out.println();
			}
			System.out.println("***********************");
		}
	}
	
	private static void test148() throws IOException {
		File rupSetFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_07_28-nshm23_branches-NSHM23_v1p4-CoulombRupSet-NSHM23_Avg-DsrUni-TotNuclRate-SubB1-"
				+ "ThreshAvgIterRelGR-IncludeThruCreep/results_NSHM23_v1p4_CoulombRupSet_branch_averaged.zip");
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(rupSetFile);
		
		int printMod = 10;
		
		RupSetMapMaker mapMaker = new RupSetMapMaker(rupSet, RupSetMapMaker.buildBufferedRegion(rupSet.getFaultSectionDataList()));
		mapMaker.setWriteGeoJSON(false);
		mapMaker.setWritePDFs(false);
		
		Stopwatch watch = Stopwatch.createStarted();
		for (int i=0; i<100; i++) {
			mapMaker.plot(new File("/tmp"), "test_map", "Test Map");
			int cnt = i+1;
			if (cnt % printMod == 0) {
				double secs = watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
				double pps = (double)cnt/secs;
				System.out.println("DONE "+cnt+", rate: "+(float)pps+" plots/sec");
			}
		}
		watch.stop();
	}
	
	private static void test149() throws IOException {
		Map<Integer, FaultSection> sects1p4 = NSHM23_FaultModels.NSHM23_v1p4.getFaultSectionIDMap();
		Map<Integer, FaultSection> sects2 = NSHM23_FaultModels.NSHM23_v2.getFaultSectionIDMap();
		
		Map<NSHM23_DeformationModels, List<? extends FaultSection>> dmSects = new HashMap<>();
		for (NSHM23_DeformationModels dm : NSHM23_DeformationModels.values())
			if (dm.getNodeWeight(null) > 0)
				dmSects.put(dm, dm.build(NSHM23_FaultModels.NSHM23_v1p4));
		
		for (int id : new int[] {2922, 1098}) {
			FaultSection sect1p4 = sects1p4.get(id);
			FaultSection sect2 = sects2.get(id);
			
			System.out.println(id+". "+sect1p4.getName());
			System.out.println("Trace, v1.4:");
			for (Location loc : sect1p4.getFaultTrace())
				System.out.println("\t"+(float)loc.getLatitude()+", "+(float)loc.getLongitude());
			System.out.println("Trace, v2:");
			for (Location loc : sect2.getFaultTrace())
				System.out.println("\t"+(float)loc.getLatitude()+", "+(float)loc.getLongitude());
			System.out.println();
			
			FaultSection s1 = sect1p4.clone();
			s1.setSectionId(0);
			FaultSection s2 = sect2.clone();
			s2.setSectionId(1);
			List<FaultSection> sects = List.of(s1, s2);
			RupSetMapMaker mapMaker = new RupSetMapMaker(sects, RupSetMapMaker.buildBufferedRegion(sects));
			
			mapMaker.setWritePDFs(false);
			mapMaker.setWriteGeoJSON(true);
			mapMaker.plot(new File("/tmp"), "sect_"+id, sect2.getName()+" Debug");
			
			sects = List.of(s1);
			mapMaker = new RupSetMapMaker(sects, RupSetMapMaker.buildBufferedRegion(sects));
			mapMaker.setWritePDFs(false);
			mapMaker.setWriteGeoJSON(false);
			mapMaker.plot(new File("/tmp"), "sect_"+id+"_prev", sect2.getName()+", V1.4");
			
			s2.setSectionId(0);
			sects = List.of(s2);
			mapMaker = new RupSetMapMaker(sects, RupSetMapMaker.buildBufferedRegion(sects));
			mapMaker.setWritePDFs(false);
			mapMaker.setWriteGeoJSON(false);
			mapMaker.plot(new File("/tmp"), "sect_"+id+"_new", sect2.getName()+", V2");
		}
		
		for (int id : sects2.keySet()) {
			FaultSection sect1p4 = sects1p4.get(id);
			FaultSection sect2 = sects2.get(id);
			
			if (sect1p4 == null) {
				System.err.println("ID mismatch, v1.4 didn't have "+id+". "+sect2.getName());
			} else {
				boolean bad = false;
				if (!sect2.getFaultTrace().equals(sect1p4.getFaultTrace())) {
					System.err.println("Trace mismatch for "+id+". "+sect2.getName());
					bad = true;
				}
				if (sect2.getAveDip() != sect1p4.getAveDip()) {
					System.err.println("Dip mismatch for "+id+". "+sect2.getName());
					bad = true;
				}
				if (sect2.getAveRake() != sect1p4.getAveRake()) {
					System.err.println("Rake mismatch for "+id+". "+sect2.getName()+": "+sect2.getAveRake()+" != "+sect1p4.getAveRake());
					bad = true;
				}
				if (sect2.getOrigAveUpperDepth() != sect1p4.getOrigAveUpperDepth()) {
					System.err.println("Upper depth mismatch for "+id+". "+sect2.getName());
					bad = true;
				}
				if (sect2.getAveLowerDepth() != sect1p4.getAveLowerDepth()) {
					System.err.println("Lower depth mismatch for "+id+". "+sect2.getName());
					bad = true;
				}
				
				if (bad) {
					System.out.println("Average slip rates:");
					for (NSHM23_DeformationModels dm : NSHM23_DeformationModels.values()) {
						List<? extends FaultSection> dmSubSects = dmSects.get(dm);
						if (dmSubSects != null) {
							int count = 0;
							double slipRate = 0d;
							List<Double> rakes = new ArrayList<>();
							for (FaultSection sect : dmSubSects) {
								if (sect.getParentSectionId() == id) {
									count++;
									slipRate += sect.getOrigAveSlipRate();
									rakes.add(sect.getAveRake());
								}
							}
							
							slipRate /= count;
							double rake = FaultUtils.getInRakeRange(FaultUtils.getAngleAverage(rakes));
							System.out.println("\t"+dm.getShortName()+": slip="+(float)slipRate+", rake="+(float)rake);
						}
					}
				}
			}
		}
	}
	
	private static void test150() throws IOException {
		File rupSetFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_07_28-nshm23_branches-NSHM23_v1p4-CoulombRupSet-NSHM23_Avg-DsrUni-TotNuclRate-SubB1-"
				+ "ThreshAvgIterRelGR-IncludeThruCreep/results_NSHM23_v1p4_CoulombRupSet_branch_averaged.zip");
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(rupSetFile);
		
		int creepingID = FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "San", "Andreas", "Creeping");
		int numThrough = 0;
		for (ClusterRupture cRup : rupSet.requireModule(ClusterRuptures.class)) {
			if (NSHM23_ConstraintBuilder.isRupThroughCreeping(creepingID, cRup))
				numThrough++;
		}
		
		System.out.println("Found "+numThrough+" / "+rupSet.getNumRuptures()+" ruptures through creeping");
	}
	
	private static void test151() throws IOException {
		Gson gson = new GsonBuilder().create();
		
		BufferedReader reader = new BufferedReader(
				new InputStreamReader(NSHM23_FaultModels.class.getResourceAsStream("/data/erf/nshm23/fault_models/v2/named_faults.json")));
		Type type = TypeToken.getParameterized(Map.class, String.class,
				TypeToken.getParameterized(List.class, Integer.class).getType()).getType();
		Map<String, List<Integer>> namedFaults = gson.fromJson(reader, type);
		
		for (String name : namedFaults.keySet())
			System.out.println(name+": "+Joiner.on(",").join(namedFaults.get(name)));
		
		Preconditions.checkState(!namedFaults.isEmpty(), "No named faults found");
		new NamedFaults(null, namedFaults);
	}
	
	private static void test152() throws IOException {
		List<? extends FaultSection> subSects = NSHM23_DeformationModels.GEOLOGIC.build(NSHM23_FaultModels.NSHM23_v2);
		List<SectMappedUncertainDataConstraint> slipData = NSHM23_PaleoDataLoader.loadU3PaleoSlipData(subSects);

		boolean applySlipRateUncertainty = true;
		MinMaxAveTracker percentTrack = new MinMaxAveTracker();
		MinMaxAveTracker percentTrackFromSlip = new MinMaxAveTracker();
		for (SectMappedUncertainDataConstraint constraint : slipData) {
			double targetSlipRate = subSects.get(constraint.sectionIndex).getReducedAveSlipRate()*1e-3;
			double targetSlipRateStdDev = subSects.get(constraint.sectionIndex).getOrigSlipRateStdDev()*1e-3;
			
			BoundedUncertainty slipUncertainty;
			UncertaintyBoundType refType;
			if (constraint.uncertainties[0] instanceof BoundedUncertainty) {
				// estimate slip rate bounds in the same units as the original uncertainty estimate
				slipUncertainty = (BoundedUncertainty)constraint.uncertainties[0];
				refType = slipUncertainty.type;
			} else {
				refType = UncertaintyBoundType.TWO_SIGMA;
				slipUncertainty = constraint.estimateUncertaintyBounds(refType);
			}
			
			System.out.println("Inferring rate constraint from paleo slip constraint on "+constraint.sectionName);
			System.out.println("\tslip="+(float)constraint.bestEstimate+"\tslipUuncert="+slipUncertainty);
			System.out.println("\tslip rate="+(float)targetSlipRate+" +/- "+targetSlipRateStdDev);
			
			// r = s / d
			double meanRate = targetSlipRate / constraint.bestEstimate;
			/*
			 * uncertainty propagation:
			 * 		r +/- deltaR = (s +/- deltaS)/(d +/- deltaD)
			 * simplifies to (see https://www.geol.lsu.edu/jlorenzo/geophysics/uncertainties/Uncertaintiespart2.html):
			 * 		deltaR/r = sqrt((deltaS/s)^2 + (deltaD/d)^2)
			 * 		deltaR = r*sqrt((deltaS/s)^2 + (deltaD/d)^2)
			 */
			double rateSD;
			if (applySlipRateUncertainty) {
				rateSD = meanRate * Math.sqrt(Math.pow(targetSlipRateStdDev/targetSlipRate, 2)
						+ Math.pow(constraint.getPreferredStdDev()/constraint.bestEstimate, 2));
//				rateSD = meanRate * (targetSlipRateStdDev/targetSlipRate + constraint.getPreferredStdDev()/constraint.bestEstimate);
			} else {
				rateSD = meanRate * constraint.getPreferredStdDev()/constraint.bestEstimate;
			}
			Uncertainty rateUncertainty = new Uncertainty(rateSD);
			
			System.out.println("\tRate:\t"+meanRate+" +/- "+rateSD);
			if (applySlipRateUncertainty) {
				double rateSD_without = meanRate * constraint.getPreferredStdDev()/constraint.bestEstimate;
				percentTrackFromSlip.addValue(100d*(rateSD - rateSD_without)/rateSD_without);
			}
			
			// now do it the previous way
//			double lowerRateTarget, upperRateTarget;
//			if (applySlipRateUncertainty) {
//				BoundedUncertainty slipRateUncertainty = refType.estimate(
//						targetSlipRate, targetSlipRateStdDev);
//				lowerRateTarget = Math.max(0d, slipRateUncertainty.lowerBound);
//				upperRateTarget = slipRateUncertainty.upperBound;
////				System.out.println("\tSlip Rate Uncertainties: "+slipRateUncertainty);
//			} else {
//				lowerRateTarget = targetSlipRate;
//				upperRateTarget = targetSlipRate;
//			}
//			rateUncertainty = new Uncertainty(refType.estimateStdDev(meanRate,
//					lowerRateTarget / slipUncertainty.upperBound,
//					upperRateTarget / slipUncertainty.lowerBound));
//			rateUncertainty = new Uncertainty(Math.sqrt(Math.pow((targetSlipRateStdDev+targetSlipRate)/constraint.bestEstimate, 2)
//					+ Math.pow(targetSlipRateStdDev/(constraint.bestEstimate+constraint.getPreferredStdDev()), 2)));
			rateUncertainty = new Uncertainty((targetSlipRateStdDev+targetSlipRate)/(constraint.bestEstimate+constraint.getPreferredStdDev()));
			System.out.println("\tOLD:\t"+meanRate+" +/- "+rateUncertainty.stdDev);
			
			percentTrack.addValue(100d*(rateSD - rateUncertainty.stdDev)/rateUncertainty.stdDev);
		}
		System.out.println("New vs Old % diffs: "+percentTrack);
		if (applySlipRateUncertainty)
			System.out.println("% change from slip rate uncertainty: "+percentTrackFromSlip);
//		double constrMean = 3.15;
////		double constrLower = -1.1;
////		double constrUpper = 7.4;
////		double constrStdDev = 2.125;
//		double constrLower = 2;
//		double constrUpper = 5;
//		double constrStdDev = UncertaintyBoundType.TWO_SIGMA.estimateStdDev(constrMean, constrLower, constrUpper);
//		
//		double targetSlipRate = 0.023757033;
//		double targetSlipRateStdDev = 0.003452321;
//		
//		
//		UncertaintyBoundType refType = UncertaintyBoundType.CONF_95;
//		BoundedUncertainty slipUncertainty = new BoundedUncertainty(refType, constrLower, constrUpper, constrStdDev);
//		
//		System.out.println("Inferring rate constraint from paleo slip constraint");
//		System.out.println("\tslip="+(float)constrMean+"\tslipUuncert="+slipUncertainty);
//		System.out.println("\tslip rate="+(float)targetSlipRate);
//		
//		// r = s / d
//		double meanRate = targetSlipRate / constrMean;
//		/*
//		 * uncertainty propagation:
//		 * 		r +/- deltaR = (s +/- deltaS)/(d +/- deltaD)
//		 * simplifies to (see https://www.geol.lsu.edu/jlorenzo/geophysics/uncertainties/Uncertaintiespart2.html):
//		 * 		deltaR/r = sqrt((deltaS/s)^2 + (deltaD/d)^2)
//		 * 		deltaR = r*sqrt((deltaS/s)^2 + (deltaD/d)^2)
//		 */
//		double rateSD;
//		if (applySlipRateUncertainty) {
//			rateSD = meanRate * Math.sqrt(Math.pow(targetSlipRateStdDev/targetSlipRate, 2)+Math.pow(constrStdDev/constrMean, 2));
//		} else {
//			rateSD = meanRate * constrStdDev/constrMean;
//		}
//		Uncertainty rateUncertainty = new Uncertainty(rateSD);
//		
//		System.out.println("RATE: "+meanRate);
//		
//		System.out.println("NEW SD: "+rateUncertainty);
//		
//		// now do it the previous way
//		double lowerRateTarget, upperRateTarget;
//		if (applySlipRateUncertainty) {
//			BoundedUncertainty slipRateUncertainty = refType.estimate(
//					targetSlipRate, targetSlipRateStdDev);
//			lowerRateTarget = Math.max(0d, slipRateUncertainty.lowerBound);
//			upperRateTarget = slipRateUncertainty.upperBound;
//			System.out.println("\tSlip Rate Uncertainties: "+slipRateUncertainty);
//		} else {
//			lowerRateTarget = targetSlipRate;
//			upperRateTarget = targetSlipRate;
//		}
//		rateUncertainty = new Uncertainty(refType.estimateStdDev(meanRate,
//				lowerRateTarget / slipUncertainty.upperBound,
//				upperRateTarget / slipUncertainty.lowerBound));
//		System.out.println("OLD SD: "+rateUncertainty);
	}
	
	private static void test153() throws IOException {
		File solFile = new File("/data/kevin/nshm23/batch_inversions/"
				+ "2022_08_12-nshm23_branches-NSHM23_v2-CoulombRupSet-NSHM23_Avg-TotNuclRate-SubB1-ThreshAvgIterRelGR/"
				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged.zip");
		FaultSystemSolution sol = FaultSystemSolution.load(solFile);
		
		NSHM23_FaultModels.NSHM23_v2.attachDefaultModules(sol.getRupSet());
		
		sol.write(new File(solFile.getParentFile(), "mod_"+solFile.getName()));
	}
	
	private static void test154() throws IOException {
		// fix bruce's weird CSVs
		File inputFile = new File("/tmp/slipLengthScaling.csv");
//		File inputFile = new File("C:\\Users\\Kevin Milner\\Downloads\\slipLengthScaling.csv");
		CSVFile<String> csv = CSVFile.readFile(inputFile, true);
		
		int numFixed = 0;
		int numFine = 0;
		for (int row=0; row<csv.getNumRows(); row++) {
			for (int col=0; col<csv.getNumCols(); col++) {
				String val = csv.get(row, col);
				if (val.startsWith("b'") && val.endsWith("'")) {
					csv.set(row, col, val.substring(2, val.length()-1));
					numFixed++;
				} else {
					numFine++;
				}
			}
		}
		System.out.println("Fixed "+numFixed+"/"+(numFine+numFixed));
		csv.writeToFile(new File(inputFile.getParentFile(), "mod_"+inputFile.getName()));
	}
	
	private static void test155() throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_08_18-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged.zip"));
		
		FaultSystemRupSet rupSet = sol.getRupSet();
		PaleoseismicConstraintData paleoData = sol.getRupSet().requireModule(PaleoseismicConstraintData.class);
		
		for (SectMappedUncertainDataConstraint data : paleoData.getPaleoSlipConstraints()) {
			if (data.sectionName.contains("Mojave")) {
				System.out.println("Paleo slip for "+data.name+" mapped to "+data.sectionName);
				
				PaleoseismicConstraintData.inferRatesFromSlipConstraints(rupSet.requireModule(SectSlipRates.class), List.of(data), true);
			}
		}
	}
	
	private static void test156() throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_08_18-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged.zip"));
		
		FaultSystemRupSet rupSet = sol.getRupSet();
		NSHM23_FaultModels.NSHM23_v2.attachDefaultModules(rupSet);
		rupSet.getModule(RegionsOfInterest.class);
	}
	
	private static void test157() throws IOException {
		List<? extends FaultSection> allSects = NSHM23_DeformationModels.GEOLOGIC.build(NSHM23_FaultModels.NSHM23_v2);
		System.out.println("State fractional slip std devs (geologic):");
		for (NSHM23_SingleStates state : NSHM23_SingleStates.values()) {
			MinMaxAveTracker track = new MinMaxAveTracker();
			for (FaultSection sect : allSects) {
				Preconditions.checkState(sect instanceof GeoJSONFaultSection);
				GeoJSONFaultSection geoSect = (GeoJSONFaultSection)sect;
				if (state.contains(geoSect))
					track.addValue(sect.getOrigSlipRateStdDev()/sect.getOrigAveSlipRate());
			}
			System.out.println(state+": "+track);
		}
	}
	
	private static void test158() {
		EvenlyDiscretizedFunc refMFD = SupraSeisBValInversionTargetMFDs.buildRefXValues(8.5d);
		double totMoRate = 1e17;
		GutenbergRichterMagFreqDist gr1 = new GutenbergRichterMagFreqDist(
				refMFD.getMinX(), refMFD.size(), refMFD.getDelta(), totMoRate, -0.5);
		GutenbergRichterMagFreqDist gr2 = new GutenbergRichterMagFreqDist(
				refMFD.getMinX(), refMFD.size(), refMFD.getDelta(), totMoRate, 0.9);
		// now average them
		double weight1 = Math.random();
		double weight2 = 1d - weight1;
		
		IncrementalMagFreqDist avg = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		
		for (int i=0; i<avg.size(); i++)
			avg.set(i, gr1.getY(i)*weight1 + gr2.getY(i)*weight2);
		
		System.out.println("GR1 moment: "+gr1.getTotalMomentRate());
		System.out.println("GR2 moment: "+gr2.getTotalMomentRate());
		System.out.println("Avg moment: "+avg.getTotalMomentRate());
	}
	
	private static void test159() {
		for (DeformationModels dm : DeformationModels.values()) {
			if (!dm.isApplicableTo(FaultModels.FM3_1) || dm == DeformationModels.MEAN_UCERF3)
				continue;
			List<? extends FaultSection> subSects = dm.build(FaultModels.FM3_1);
			double min = Double.POSITIVE_INFINITY;
			for (FaultSection sect : subSects)
//				min = Math.min(min, sect.getReducedAveSlipRate());
				min = Math.min(min, sect.getOrigAveSlipRate());
			System.out.println(dm+": min slip rate: "+min);
		}
	}
	
	private static void test160() throws IOException {
		GriddedRegion relm1 = new CaliforniaRegions.RELM_TESTING_GRIDDED(0.1d);
		GriddedRegion relm2 = new GriddedRegion(SeismicityRegions.CONUS_U3_RELM.load(), 0.1d, GriddedRegion.ANCHOR_0_0);
		
		System.out.println("orig has "+relm1.getNodeCount()+" locs");
		System.out.println("new has "+relm2.getNodeCount()+" locs");
	}
	
	private static void test161() throws IOException {
		File dir = new File("/home/kevin/Downloads/tmp/DeclusteringSmoothingGrids/CEUS");
		for (File file : dir.listFiles()) {
			System.out.println(file.getName());
			CSVFile<String> csv = CSVFile.readFile(file, true);
			double sum = 0d;
			for (int row=0; row<csv.getNumRows(); row++)
				sum+= csv.getDouble(row, 2);
			System.out.println(sum);
		}
	}
	
	private static void test162() throws IOException {
		LogicTreeBranch<?> branch = new LogicTreeBranch<>(List.of(NSHM23_LogicTreeBranch.SINGLE_STATES), List.of(NSHM23_SingleStates.NM));
		NSHM23_FaultModels.getDefaultRegion(branch);
	}
	
	private static void test163() throws IOException {
		Region reg = NSHM23_RegionLoader.loadFullConterminousWUS();
		
		double refSpacing = 0.2;
		int refCount = new GriddedRegion(reg, refSpacing, GriddedRegion.ANCHOR_0_0).getNodeCount();
		
		DecimalFormat refDF = new DecimalFormat("0.00");
		
		double[] spacings = { 0.1, 0.2, 0.25, 0.3, 1d/3d, 0.4, 0.5 };
		for (double spacing : spacings) {
			GriddedRegion gridReg = new GriddedRegion(reg, spacing, GriddedRegion.ANCHOR_0_0);
			double ratio = (double)gridReg.getNodeCount()/(double)refCount;
			System.out.println((float)spacing+": "+gridReg.getNodeCount()+" points\t("+refDF.format(ratio)+"x "+(float)refSpacing+")");
		}
	}
	
	private static void test164() throws IOException {
		File mainDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/2022_08_22-nshm23_branches-NSHM23_v2-"
				+ "CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
		File resultsFile = new File(mainDir, "results_gridded_branches.zip");
		
		SolutionLogicTree slt = SolutionLogicTree.load(resultsFile);
		
		for (LogicTreeBranch<?> branch : slt.getLogicTree()) {
			FaultSystemSolution sol = slt.forBranch(branch);
			GridSourceProvider gridProv = sol.getGridSourceProvider();
			System.out.println(gridProv.getName()+": "+gridProv.getClass().getName());
		}
	}
	
	private static void test165() throws IOException {
		NSHM23_FaultModels fm = NSHM23_FaultModels.NSHM23_v2;
		System.out.println(fm.getFaultSections().size()+" sections");
		NSHM23_DeformationModels dm = NSHM23_DeformationModels.GEOLOGIC;
		System.out.println(dm.build(fm).size()+" subsections");
	}
	
	private static void test166() throws IOException {
		EvenlyDiscretizedFunc refMFD = new EvenlyDiscretizedFunc(
				0.05, 120, 0.1);
		System.out.println(refMFD);
		File dir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_08_17-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-SubB1-ThreshAvgIterRelGR-360_samples");
		File resultsFile = new File(dir, "results.zip");
		SolutionLogicTree slt = SolutionLogicTree.load(resultsFile);
		LogicTree<?> tree = slt.getLogicTree();
		tree = tree.sample(10, false, new Random(tree.size()));
		boolean process = false;
		File outputFile = new File(dir, "results_NSHM23_v2_CoulombRupSet_branch_averaged_rebuild.zip");
		BranchAverageSolutionCreator baCreator = new BranchAverageSolutionCreator(tree.getWeightProvider());
		for (LogicTreeBranch<?> branch : tree) {
			FaultSystemSolution sol = slt.forBranch(branch, process);
			if (!process)
				// add ROI, etc
				NSHM23_FaultModels.NSHM23_v2.attachDefaultModules(sol.getRupSet());
			baCreator.addSolution(sol, branch);
		}
		FaultSystemSolution ba = baCreator.build();
		NSHM23_FaultModels.NSHM23_v2.attachDefaultModules(ba.getRupSet());
		ba.write(outputFile);
	}
	
	private static void test167() throws IOException {
		EvenlyDiscretizedFunc refMFD = SupraSeisBValInversionTargetMFDs.buildRefXValues(8d);
		
		int index5 = refMFD.getClosestXIndex(5.01d);
		GutenbergRichterMagFreqDist fullGR = new GutenbergRichterMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		fullGR.setAllButTotMoRate(refMFD.getX(index5), refMFD.getX(refMFD.size()-1), 1d, 1d);
		double totRate = fullGR.calcSumOfY_Vals();
		for (int i=index5; i<refMFD.size(); i++) {
			double rateAbove = 0d;
			for (int j=i; j<fullGR.size(); j++)
				rateAbove += fullGR.getY(j);
			System.out.println("Fraction of rate >= "+(float)refMFD.getX(i)+": "+(float)(rateAbove/totRate));
		}
	}
	
	private static void test168() throws IOException {
		for (FaultSection sect : NSHM23_FaultModels.NSHM23_v2.getFaultSections()) {
			if (sect.getAveDip() < 30) {
				System.out.println(sect.getName()+": dip="+(float)sect.getAveDip()+"; depths=["
						+(float)sect.getOrigAveUpperDepth()+","+(float)sect.getAveLowerDepth()+"]");
			}
		}
	}
	
	private static void test169() throws IOException {
		File jsonFile = new File("C:\\Users\\Kevin Milner\\Downloads\\grid_region.geojson");
		Feature feature = Feature.read(jsonFile);
		for (Geometry geom : ((GeometryCollection)feature.geometry).geometries) {
			if (geom.type == GeoJSON_Type.MultiPoint) {
				System.out.println("MultiPoint has "+((MultiPoint)geom).points.size()+" nodes");
			}
		}
		GriddedRegion gridReg = GriddedRegion.fromFeature(feature);
		System.out.println("Grid reg has "+gridReg.getNodeCount()+" nodes");
	}
	
	private static void test170() throws IOException {
		File outputDir = new File("C:\\Users\\Kevin Milner\\Downloads");
		List<? extends FaultSection> geoSects = DeformationModels.GEOLOGIC.build(FaultModels.FM3_1);
		GeoJSONFaultReader.writeFaultSections(new File(outputDir, "u3_fm3_1_geol_sub_sects.geojson"), geoSects);
		geoSects = NSHM23_DeformationModels.GEOLOGIC.build(NSHM23_FaultModels.NSHM23_v2);
		GeoJSONFaultReader.writeFaultSections(new File(outputDir, "nshm23_geol_sub_sects.geojson"), geoSects);
	}
	
	private static void test171() throws IOException {
		int branchCount = 2250*54;
		int nodeCount = new GriddedRegion(NSHM23_RegionLoader.loadFullConterminousWUS(), 0.2, GriddedRegion.ANCHOR_0_0).getNodeCount();
		long valCount = nodeCount * branchCount;
		System.out.println(branchCount+" branches x "+nodeCount+" nodes = "+valCount+" values");
		long bytes = valCount * 8l;
		long kb = bytes / 1024l;
		long mb = kb / 1024l;
		long gb = mb / 1024l;
		System.out.println("\t"+bytes+" bytes");
		System.out.println("\t"+kb+" kb");
		System.out.println("\t"+mb+" mb");
		System.out.println("\t"+gb+" gb");
	}
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		test166();
	}

}
