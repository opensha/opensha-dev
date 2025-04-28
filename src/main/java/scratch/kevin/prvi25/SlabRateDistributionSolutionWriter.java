package scratch.kevin.prvi25;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.TimeUnit;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter;
import org.opensha.commons.hpc.mpj.MPJExpressShellScriptWriter;
import org.opensha.commons.hpc.mpj.NoMPJSingleNodeShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.USC_CARC_ScriptWriter;
import org.opensha.commons.logicTree.Affects;
import org.opensha.commons.logicTree.DoesNotAffect;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.logicTree.LogicTreeNode.RandomlySampledNode;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.io.archive.ArchiveOutput;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_LogicTreeHazardCalc;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_SingleSolHazardCalc;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.MFDGridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.faultSysSolution.util.SolLogicTreeSampler;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.PRVI25_GridSourceBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateFileLoader.RateType;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_DeclusteringAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SeisSmoothingAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionCaribbeanSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionMuertosSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader.PRVI25_SeismicityRegions;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;

public class SlabRateDistributionSolutionWriter {

	public static void main(String[] args) throws IOException {
		File ratesDir = new File("/home/kevin/OpenSHA/nshm23/prvi/rate_raw_data");
		File invsDir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/");
		FaultSystemSolution subSol = FaultSystemSolution.load(new File(invsDir,
				"2025_01_17-prvi25_subduction_branches/results_PRVI_SUB_FM_LARGE_branch_averaged.zip"));
		
		File fullRandOutDir = new File(invsDir, "2025_04_21-prvi_slab_variability-full_random");
		File threeBranchExactOutDir = new File(invsDir, "2025_04_21-prvi_slab_variability-three_branch-exact");
		File threeBranchAltOutDir = new File(invsDir, "2025_04_21-prvi_slab_variability-three_branch-m1_to_mmax");
		
		Preconditions.checkState(fullRandOutDir.exists() || fullRandOutDir.mkdir());
		Preconditions.checkState(threeBranchExactOutDir.exists() || threeBranchExactOutDir.mkdir());
		Preconditions.checkState(threeBranchAltOutDir.exists() || threeBranchAltOutDir.mkdir());
		
		List<double[]> carPairs = loadRates(new File(ratesDir, "rbpairs-CAR Intraslab-Full-v3.csv"));
		List<double[]> muePairs = loadRates(new File(ratesDir, "rbpairs-MUE Intraslab-Full-v3.csv"));
		int numSamples = 5000;
		int numTotalSamples = numSamples*3*2; // times declustering and smoothing branches
		Random rand = new Random(numSamples);
		
		List<LogicTreeLevel<? extends LogicTreeNode>> fullRandLevels = new ArrayList<>();
		fullRandLevels.add(PRVI25_LogicTreeBranch.SEIS_DECLUSTER);
		fullRandLevels.add(PRVI25_LogicTreeBranch.SEIS_SMOOTH);
		CarRateSamplingLevel carSampler = new CarRateSamplingLevel(carPairs);
		carSampler.buildNodes(rand, numTotalSamples);
		MueRateSamplingLevel mueSampler = new MueRateSamplingLevel(muePairs);
		mueSampler.buildNodes(rand, numTotalSamples);
		fullRandLevels.add(carSampler);
		fullRandLevels.add(mueSampler);
		List<LogicTreeLevel<? extends LogicTreeNode>> threeBranchLevels = new ArrayList<>();
		threeBranchLevels.add(PRVI25_LogicTreeBranch.SEIS_DECLUSTER);
		threeBranchLevels.add(PRVI25_LogicTreeBranch.SEIS_SMOOTH);
		threeBranchLevels.add(PRVI25_LogicTreeBranch.CAR_SEIS_RATE);
		threeBranchLevels.add(PRVI25_LogicTreeBranch.MUE_SEIS_RATE);
		
		Region reg = PRVI25_RegionLoader.loadPRVI_Tight();
		GriddedRegion gridReg = new GriddedRegion(reg, 0.05, GriddedRegion.ANCHOR_0_0);
		System.out.println("Region has "+gridReg.getNodeCount()+" nodes");
		writeHazardScripts(fullRandOutDir, gridReg);
		writeHazardScripts(threeBranchExactOutDir, gridReg);
		writeHazardScripts(threeBranchAltOutDir, gridReg);
//		System.exit(0);
		
		SolutionLogicTree.FileBuilder fullRandBuilder = new SolutionLogicTree.FileBuilder(
//				new ArchiveOutput.AsynchronousZipFileOutput(new File(fullRandOutDir, "results.zip")));
				new ArchiveOutput.ParallelZipFileOutput(new File(fullRandOutDir, "results.zip"), 20));
		fullRandBuilder.setSerializeGridded(true);
		SolutionLogicTree.FileBuilder threeBranchExactBuilder = new SolutionLogicTree.FileBuilder(
				new ArchiveOutput.AsynchronousZipFileOutput(new File(threeBranchExactOutDir, "results.zip")));
		threeBranchExactBuilder.setSerializeGridded(true);
		SolutionLogicTree.FileBuilder threeBranchAltBuilder = new SolutionLogicTree.FileBuilder(
				new ArchiveOutput.AsynchronousZipFileOutput(new File(threeBranchAltOutDir, "results.zip")));
		threeBranchAltBuilder.setSerializeGridded(true);
		
		List<CarSamplingNode> carSamples = carSampler.getNodes();
		List<MueSamplingNode> mueSamples = mueSampler.getNodes();
		
		IncrementalMagFreqDist refMFD = FaultSysTools.initEmptyMFD(PRVI25_GridSourceBuilder.OVERALL_MMIN, PRVI25_GridSourceBuilder.SLAB_MMAX);
		
		subSol.setVerbose(false);
		DecimalFormat pDF = new DecimalFormat("0.0%");
		
		int sampleIndex = 0;
		for (PRVI25_DeclusteringAlgorithms decluster : PRVI25_DeclusteringAlgorithms.values()) {
			if (decluster.getNodeWeight(null) == 0d)
				continue;
			for (PRVI25_SeisSmoothingAlgorithms smooth : PRVI25_SeisSmoothingAlgorithms.values()) {
				if (smooth.getNodeWeight(null) == 0d)
					continue;
				System.out.println("Building 3-branch for "+decluster+", "+smooth);
				for (PRVI25_SubductionCaribbeanSeismicityRate carRate : PRVI25_SubductionCaribbeanSeismicityRate.values()) {
					if (carRate.getNodeWeight(null) == 0d)
						continue;
					for (PRVI25_SubductionMuertosSeismicityRate mueRate : PRVI25_SubductionMuertosSeismicityRate.values()) {
						if (mueRate.getNodeWeight(null) == 0d)
							continue;
						LogicTreeBranch<LogicTreeNode> branch = new LogicTreeBranch<>(threeBranchLevels);
						branch.setValue(decluster);
						branch.setValue(smooth);
						branch.setValue(carRate);
						branch.setValue(mueRate);
						
						PRVI25_SubductionCaribbeanSeismicityRate.TYPE = RateType.EXACT;
						PRVI25_SubductionMuertosSeismicityRate.TYPE = RateType.EXACT;
						GridSourceList exact = PRVI25_GridSourceBuilder.buildSlabGridSourceList(branch);
						
						subSol.setGridSourceProvider(exact);
						threeBranchExactBuilder.solution(subSol, branch);
						
						PRVI25_SubductionCaribbeanSeismicityRate.TYPE = RateType.M1_TO_MMAX;
						PRVI25_SubductionMuertosSeismicityRate.TYPE = RateType.M1_TO_MMAX;
						GridSourceList m1tommax = PRVI25_GridSourceBuilder.buildSlabGridSourceList(branch);
						
						subSol.setGridSourceProvider(m1tommax);
						threeBranchAltBuilder.solution(subSol, branch);
					}
				}
				System.out.println("Building "+numSamples+" for "+decluster+", "+smooth);
				CompletableFuture<Void> prevFuture = null;
				
				Stopwatch totalWatch = Stopwatch.createStarted();
				Stopwatch ioWatch = Stopwatch.createUnstarted();
				for (int s=0; s<numSamples; s++) {
					if (s % 10 == 00) {
						double ioTime = ioWatch.elapsed(TimeUnit.MILLISECONDS) / 1000d;
						double totTime = totalWatch.elapsed(TimeUnit.MILLISECONDS) / 1000d;
						System.out.println("\tBuilding sample "+s+"/"+numSamples+" ("+pDF.format(ioTime/totTime)+" waiting on blocking I/O)");
					}
					CarSamplingNode carSample = carSamples.get(sampleIndex);
					MueSamplingNode mueSample = mueSamples.get(sampleIndex);
					sampleIndex++;
					LogicTreeBranch<LogicTreeNode> branch = new LogicTreeBranch<>(fullRandLevels);
					branch.setValue(decluster);
					branch.setValue(smooth);
					branch.setValue(carSample);
					branch.setValue(mueSample);
					
					GutenbergRichterMagFreqDist carGR = buildGR(carSample.rate, carSample.b, PRVI25_GridSourceBuilder.SLAB_MMAX, 5d, refMFD);
					GutenbergRichterMagFreqDist mueGR = buildGR(mueSample.rate, mueSample.b, PRVI25_GridSourceBuilder.SLAB_MMAX, 5d, refMFD);
					
					if (sampleIndex == 1) {
						System.out.println(carGR);
						System.out.println(mueGR);
					}
					
					GridSourceList carList = PRVI25_GridSourceBuilder.buildSlabGridSourceList(branch, PRVI25_SeismicityRegions.CAR_INTRASLAB, carGR);
					GridSourceList mueList = PRVI25_GridSourceBuilder.buildSlabGridSourceList(branch, PRVI25_SeismicityRegions.MUE_INTRASLAB, mueGR);
					
					GridSourceList combList = GridSourceList.combine(carList, mueList);
					
					if (prevFuture != null) {
						ioWatch.start();
						prevFuture.join();
						ioWatch.stop();
					}
					
					prevFuture = CompletableFuture.runAsync(new Runnable() {
						
						@Override
						public void run() {
							subSol.setGridSourceProvider(combList);
							try {
								fullRandBuilder.solution(subSol, branch);
							} catch (IOException e) {
								e.printStackTrace();
								System.exit(1);
							}
						}
					});
				}
				prevFuture.join();
				
				totalWatch.stop();
			}
		}
		
		fullRandBuilder.close();
		threeBranchExactBuilder.close();
		threeBranchAltBuilder.close();
	}
	
	private static GutenbergRichterMagFreqDist buildGR(double rateAboveM1, double b, double mMax, double m1, IncrementalMagFreqDist refMFD) {
		GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		// this sets shape, min/max
		// subtract a tiny amount from mMax so that if it's exactly at a bin edge, e.g. 7.9, it rounds down, e.g. to 7.85
		gr.setAllButTotCumRate(refMFD.getX(0), refMFD.getX(refMFD.getClosestXIndex(mMax-0.001)), 1e16, b);
		// this scales it to match
		// similarly, add a tiny amount to M1 so that if it's exactly at a bin edge (which it should be as it's determined
		// using cumulative binning), it rounds up to the incremental bin for that cumulative edge
		gr.scaleToCumRate(refMFD.getClosestXIndex(m1+0.001), rateAboveM1);
		return gr;
	}
	
	private static List<double[]> loadRates(File csvFile) throws IOException {
		CSVFile<String> csv = CSVFile.readFile(csvFile, false);
		
		boolean reading = false;
		List<double[]> ret = new ArrayList<>();
		
		MinMaxAveTracker rateTrack = new MinMaxAveTracker();
		MinMaxAveTracker bTrack = new MinMaxAveTracker();
		for (int row=0; row<csv.getNumRows(); row++) {
			if (reading) {
				double b = csv.getDouble(row, 0);
				double rate = csv.getDouble(row, 1);
				rateTrack.addValue(rate);
				bTrack.addValue(b);
				ret.add(new double[] {rate, b});
			} else if (csv.get(row, 0).startsWith("next ") && csv.get(row, 0).contains("lines are b")) {
				reading = true;
			}
		}
		System.out.println("Loaded "+ret.size()+" rate and b-value pairs from "+csvFile.getName());
		System.out.println("\trates: "+rateTrack);
		System.out.println("\tbs: "+bTrack);
		return ret;
	}
	
	@DoesNotAffect(FaultSystemRupSet.SECTS_FILE_NAME)
	@DoesNotAffect(FaultSystemRupSet.RUP_SECTS_FILE_NAME)
	@DoesNotAffect(FaultSystemRupSet.RUP_PROPS_FILE_NAME)
	@DoesNotAffect(FaultSystemSolution.RATES_FILE_NAME)
	@DoesNotAffect(GridSourceProvider.ARCHIVE_GRID_REGION_FILE_NAME)
	@DoesNotAffect(MFDGridSourceProvider.ARCHIVE_MECH_WEIGHT_FILE_NAME)
	@DoesNotAffect(GridSourceList.ARCHIVE_GRID_LOCS_FILE_NAME)
	@Affects(MFDGridSourceProvider.ARCHIVE_SUB_SEIS_FILE_NAME)
	@Affects(MFDGridSourceProvider.ARCHIVE_UNASSOCIATED_FILE_NAME)
	@Affects(GridSourceList.ARCHIVE_GRID_SOURCES_FILE_NAME)
	private static class CarSamplingNode implements RandomlySampledNode {

		private String name;
		private String shortName;
		private String prefix;
		private double weight;
		private long seed;
		private double rate;
		private double b;

		private CarSamplingNode() {
			
		}
		
		private CarSamplingNode(String name, String shortName, String prefix, double weight, long seed, double rate, double b) {
			super();
			this.name = name;
			this.shortName = shortName;
			this.prefix = prefix;
			this.weight = weight;
			this.seed = seed;
			this.rate = rate;
			this.b = b;
		}

		@Override
		public double getNodeWeight(LogicTreeBranch<?> fullBranch) {
			return weight;
		}

		@Override
		public String getFilePrefix() {
			return prefix;
		}

		@Override
		public String getShortName() {
			return shortName;
		}

		@Override
		public String getName() {
			return name;
		}

		@Override
		public long getSeed() {
			return seed;
		}

		@Override
		public void init(String name, String shortName, String prefix, double weight, long seed) {
			this.name = name;
			this.shortName = shortName;
			this.prefix = prefix;
			this.weight = weight;
			this.seed = seed;
		}
		
	}
	
	private static class CarRateSamplingLevel extends LogicTreeLevel.RandomlySampledLevel<CarSamplingNode> {
		
		private List<double[]> samples;
		private List<double[]> randomizedSamples;
		
		private CarRateSamplingLevel() {}

		public CarRateSamplingLevel(List<double[]> samples) {
			this.samples = samples;
		}

		@Override
		public String getShortName() {
			return "CAR-Sampling";
		}

		@Override
		public String getName() {
			return "CAR Rate/b Distribution Sampling";
		}

		@Override
		public CarSamplingNode buildNodeInstance(int index, long seed, double weight) {
			if (randomizedSamples == null) {
				randomizedSamples = new ArrayList<>(samples);
				Collections.shuffle(randomizedSamples, new Random(seed));
			}
			double[] sample = randomizedSamples.get(index);
			return new CarSamplingNode("CAR Sample "+index, "CAR-Sample"+index, "car_sample_"+index, weight, seed,
					sample[0], sample[1]);
		}

		@Override
		public Class<? extends CarSamplingNode> getType() {
			return CarSamplingNode.class;
		}
		
	}
	
	@DoesNotAffect(FaultSystemRupSet.SECTS_FILE_NAME)
	@DoesNotAffect(FaultSystemRupSet.RUP_SECTS_FILE_NAME)
	@DoesNotAffect(FaultSystemRupSet.RUP_PROPS_FILE_NAME)
	@DoesNotAffect(FaultSystemSolution.RATES_FILE_NAME)
	@DoesNotAffect(GridSourceProvider.ARCHIVE_GRID_REGION_FILE_NAME)
	@DoesNotAffect(MFDGridSourceProvider.ARCHIVE_MECH_WEIGHT_FILE_NAME)
	@DoesNotAffect(GridSourceList.ARCHIVE_GRID_LOCS_FILE_NAME)
	@Affects(MFDGridSourceProvider.ARCHIVE_SUB_SEIS_FILE_NAME)
	@Affects(MFDGridSourceProvider.ARCHIVE_UNASSOCIATED_FILE_NAME)
	@Affects(GridSourceList.ARCHIVE_GRID_SOURCES_FILE_NAME)
	private static class MueSamplingNode implements RandomlySampledNode {

		private String name;
		private String shortName;
		private String prefix;
		private double weight;
		private long seed;
		private double rate;
		private double b;

		private MueSamplingNode() {
			
		}
		
		private MueSamplingNode(String name, String shortName, String prefix, double weight, long seed, double rate, double b) {
			super();
			this.name = name;
			this.shortName = shortName;
			this.prefix = prefix;
			this.weight = weight;
			this.seed = seed;
			this.rate = rate;
			this.b = b;
		}

		@Override
		public double getNodeWeight(LogicTreeBranch<?> fullBranch) {
			return weight;
		}

		@Override
		public String getFilePrefix() {
			return prefix;
		}

		@Override
		public String getShortName() {
			return shortName;
		}

		@Override
		public String getName() {
			return name;
		}

		@Override
		public long getSeed() {
			return seed;
		}

		@Override
		public void init(String name, String shortName, String prefix, double weight, long seed) {
			this.name = name;
			this.shortName = shortName;
			this.prefix = prefix;
			this.weight = weight;
			this.seed = seed;
		}
		
	}
	
	private static class MueRateSamplingLevel extends LogicTreeLevel.RandomlySampledLevel<MueSamplingNode> {
		
		private List<double[]> samples;
		private List<double[]> randomizedSamples;
		
		private MueRateSamplingLevel() {}

		public MueRateSamplingLevel(List<double[]> samples) {
			this.samples = samples;
		}

		@Override
		public String getShortName() {
			return "MUE-Sampling";
		}

		@Override
		public String getName() {
			return "MUE Rate/b Distribution Sampling";
		}

		@Override
		public MueSamplingNode buildNodeInstance(int index, long seed, double weight) {
			if (randomizedSamples == null) {
				randomizedSamples = new ArrayList<>(samples);
				Collections.shuffle(randomizedSamples, new Random(seed));
			}
			double[] sample = randomizedSamples.get(index);
			return new MueSamplingNode("MUE Sample "+index, "MUE-Sample"+index, "mue_sample_"+index, weight, seed,
					sample[0], sample[1]);
		}

		@Override
		public Class<? extends MueSamplingNode> getType() {
			return MueSamplingNode.class;
		}
		
	}
	
	private static void writeHazardScripts(File outputDir, GriddedRegion gridReg) throws IOException {
		String dirName = outputDir.getName();
		File localDir = outputDir;
		
		File remoteMainDir = new File("/project/scec_608/kmilner/nshm23/batch_inversions");
		int remoteTotalThreads = 20;
		int remoteTotalMemGB = 50;
		String queue = "scec";
		int nodes = 36;
//		JavaShellScriptWriter mpjWrite = new MPJExpressShellScriptWriter(
//				USC_CARC_ScriptWriter.JAVA_BIN, remoteTotalMemGB*1024, null, USC_CARC_ScriptWriter.MPJ_HOME);
		JavaShellScriptWriter mpjWrite = new FastMPJShellScriptWriter(
				USC_CARC_ScriptWriter.JAVA_BIN, remoteTotalMemGB*1024, null, USC_CARC_ScriptWriter.FMPJ_HOME);
//		JavaShellScriptWriter mpjWrite = new NoMPJSingleNodeShellScriptWriter(USC_CARC_ScriptWriter.JAVA_BIN,
//				remoteTotalMemGB*1024, null); nodes = 1; remoteInversionsPerBundle = 2;
		BatchScriptWriter pbsWrite = new USC_CARC_ScriptWriter();
		
		mpjWrite.setEnvVar("MAIN_DIR", remoteMainDir.getAbsolutePath());
		String mainDirPath = "$MAIN_DIR";
		mpjWrite.setEnvVar("DIR", mainDirPath+"/"+dirName);
		String dirPath = "$DIR";
		
		List<File> classpath = new ArrayList<>();
		classpath.add(new File(dirPath+"/opensha-dev-all.jar"));
		if (mpjWrite instanceof NoMPJSingleNodeShellScriptWriter)
			classpath.add(new File("/project/scec_608/kmilner/git/opensha/lib/mpj-0.38.jar"));
		
		mpjWrite.setClasspath(classpath);
		if (mpjWrite instanceof MPJExpressShellScriptWriter)
			((MPJExpressShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
		else if (mpjWrite instanceof FastMPJShellScriptWriter)
			((FastMPJShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
		
		File gridRegFile = new File(outputDir, "gridded_region.geojson");
		Feature.write(gridReg.toFeature(), gridRegFile);
		
		Double sigmaTrunc = 3d;
		
		String mapScriptName = "batch_hazard_slab_only.slurm";
		
		String argz = "--input-file "+dirPath+"/results.zip";
		argz += " --output-dir "+dirPath+"/results_slab_only";
		argz += " --output-file "+dirPath+"/results_slab_only.zip";
		argz += " --gridded-seis "+IncludeBackgroundOption.ONLY.name();
		argz += " --quick-grid-calc";
		argz += " --region "+dirPath+"/"+gridRegFile.getName();
		argz += " --gmpe "+AttenRelRef.USGS_PRVI_SLAB.name();
		if (sigmaTrunc != null)
			argz += " --gmm-sigma-trunc-one-sided "+sigmaTrunc.floatValue();
		argz += " "+MPJTaskCalculator.argumentBuilder().minDispatch(1).maxDispatch(100).threads(remoteTotalThreads).build();
		List<String> script = mpjWrite.buildScript(MPJ_LogicTreeHazardCalc.class.getName(), argz);
		pbsWrite.writeScript(new File(localDir, mapScriptName), script, 1440, nodes, remoteTotalThreads, queue);
	}

}
