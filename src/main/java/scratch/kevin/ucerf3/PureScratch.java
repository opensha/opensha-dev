package scratch.kevin.ucerf3;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.opensha.commons.data.IntegerSampler.ExclusionIntegerSampler;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.IntegerPDF_FunctionSampler;
import org.opensha.commons.data.uncertainty.UncertaintyBoundType;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.modules.OpenSHA_Module;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfigurationFactory;
import org.opensha.sha.earthquake.faultSysSolution.inversion.Inversions;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.InversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.UncertainDataConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.TimeCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupMFDsModule;
import org.opensha.sha.earthquake.faultSysSolution.modules.WaterLevelRates;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.SolMFDPlot;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.GeoJSONFaultReader;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;

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
		
		LogicTreeBranch<?> branch = NSHM23_LogicTreeBranch.DEFAULT;
		
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
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		test120();
	}

}
