package scratch.kevin.ucerf3.eal;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.zip.ZipFile;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.io.IOUtils;
import org.dom4j.DocumentException;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.util.ClassUtils;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FilePathComparator;
import org.opensha.sha.earthquake.param.BackgroundRupParam;
import org.opensha.sha.earthquake.param.BackgroundRupType;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sra.calc.parallel.MPJ_CondLossCalc;
import org.opensha.sra.gui.portfolioeal.Asset;
import org.opensha.sra.gui.portfolioeal.Portfolio;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import com.google.common.collect.Lists;
import com.google.common.io.Files;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import mpi.MPI;
import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.mean.TrueMeanBuilder;
import scratch.UCERF3.griddedSeismicity.GridSourceProvider;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;
import scratch.UCERF3.logicTree.LogicTreeBranchNode;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.kevin.ucerf3.eal.branches.U3_EAL_GMM_Epistemic;
import scratch.kevin.ucerf3.eal.branches.U3_EAL_GMMs;

import scratch.kevin.ucerf3.eal.branches.U3_EAL_LogicTreeBranch;
import scratch.kevin.ucerf3.eal.branches.U3_EAL_ProbModels;
import scratch.kevin.ucerf3.eal.branches.U3_EAL_Vs30Model;

public class MPJ_UCERF3_EAL_Combiner extends MPJTaskCalculator {
	
	private FaultSystemSolution trueMeanSol;
	private FaultSystemSolutionERF erf;
	private Map<U3LogicTreeBranch, List<Integer>> mappings;
	private CompoundFaultSystemSolution cfss;
	
	private double erfProbsDuration;
	private Map<U3_EAL_ProbModels, ZipFile> probsZipFiles;
	
	private Map<U3_EAL_Vs30Model, File> vs30Dirs;

	private File outputDir;
	private File resultsDir;
	private File resultsFile;
	private File tractCacheDir;

	private List<String> allTracts;
	
	private File tmpDir;
	
	private List<U3_EAL_LogicTreeBranch> branches = new ArrayList<>();
	
	private String[] tractNames;
	
	private LoadingCache<File, double[][]> rupLossesCache;
	private LoadingCache<File, DiscretizedFunc[]> griddedLossesCache;
	private LoadingCache<File, TractLoader> tractCache;
	
	private ExecutorService exec;
	
	private boolean consolidateOnly;
	private CSVFile<String> resultsCSV;
	
	private DiscretizedFunc lecXVals;
	private double[] lecProbLevels;
	private CSVFile<String> lecResultsCSV;
	private File lecResultsFileCSV;
	
	private LossCOV_Model covModel;

	public MPJ_UCERF3_EAL_Combiner(CommandLine cmd, File outputDir) throws IOException, DocumentException {
		super(cmd);
		this.shuffle = false;
		this.outputDir = outputDir;
		
		File trueMeanSolFile = new File(cmd.getOptionValue("true-mean-sol"));
		File compoundSolFile = new File(cmd.getOptionValue("compound-sol"));
		
		consolidateOnly = cmd.hasOption("consolidate-only");
		
		if (rank == 0)
			Preconditions.checkState(outputDir.exists() || outputDir.mkdir(),
					"Output dir doesn't exist and could not be created: %s", outputDir.getAbsolutePath());
		resultsDir = new File(outputDir, "results");
		if (consolidateOnly && rank == 0) {
			debug("Consolidating only");
			Preconditions.checkState(resultsDir.exists(), "Consolidate only but no results dir");
			size = 0;
			for (File file : resultsDir.listFiles())
				if (file.getName().startsWith("results_") && file.getName().endsWith(".csv"))
					size++;
			debug("new size="+size+" for consolidation");
		}
		
		if (rank == 0)
			debug("Loading true mean solution from: "+trueMeanSolFile.getAbsolutePath());
		trueMeanSol = FaultSystemIO.loadSol(trueMeanSolFile);
		// now load in the mappings
		if (rank == 0)
			debug("Loading true mean branch mappings");
		mappings = TrueMeanBuilder.loadRuptureMappings(trueMeanSolFile);
		if (rank == 0)
			debug("Loading CFSS: "+compoundSolFile.getAbsolutePath());
		cfss = new CachedGridSourceCFSS(new ZipFile(compoundSolFile));
		
		if (rank == 0)
			Preconditions.checkState(resultsDir.exists() || resultsDir.mkdir(),
					"Results dir doesn't exist and could not be created: %s", resultsDir.getAbsolutePath());
		
		erfProbsDuration = Double.parseDouble(cmd.getOptionValue("erf-probs-duration"));
		
		File probsZipDir = new File(cmd.getOptionValue("erf-probs-dir"));
		Preconditions.checkState(probsZipDir.exists(), "Probs zip file doesn't exist: %s", probsZipDir.getAbsolutePath());
		
		probsZipFiles = new HashMap<>();
		for (U3_EAL_ProbModels probModel : U3_EAL_ProbModels.values()) {
			File probFile = new File(probsZipDir, "probs_"+(float)+erfProbsDuration+"yr_"+probModel.getShortName()+".zip");
			if (!probFile.exists() && (float)Math.round(erfProbsDuration) == (float)erfProbsDuration)
				probFile = new File(probsZipDir, "probs_"+(int)+erfProbsDuration+"yr_"+probModel.getShortName()+".zip");
			if (probFile.exists())
				probsZipFiles.put(probModel, new ZipFile(probFile));
		}
		Preconditions.checkState(!probsZipFiles.isEmpty(), "No prob zip files with duration=%s found in %s",
				(float)erfProbsDuration, probsZipDir.getAbsolutePath());
		
		vs30Dirs = new HashMap<>();
		if (cmd.hasOption("wills-dir")) {
			File willsDir = new File(cmd.getOptionValue("wills-dir"));
			Preconditions.checkState(willsDir.exists(), "Wills dir doesn't exist: %s", willsDir.getAbsolutePath());
			vs30Dirs.put(U3_EAL_Vs30Model.WILLS_2015, willsDir);
		}
		if (cmd.hasOption("wald-dir")) {
			File waldDir = new File(cmd.getOptionValue("wald-dir"));
			Preconditions.checkState(waldDir.exists(), "Wald dir doesn't exist: %s", waldDir.getAbsolutePath());
			vs30Dirs.put(U3_EAL_Vs30Model.WALD_ALLEN, waldDir);
		}
		Preconditions.checkArgument(!vs30Dirs.isEmpty(), "No Vs30 model directories supplied!");
		tractNames = null;
		if (cmd.hasOption("tract")) {
			String tractStr = cmd.getOptionValue("tract");
			tractNames = tractStr.split(",");
			if (rank == 0)
				debug("Calculating for "+tractNames.length+" tracts");
			Preconditions.checkState(!cmd.hasOption("tract-location") && !cmd.hasOption("tract-radius"),
					"Tract location and radius not supported when --tract option used");
		}
		if (cmd.hasOption("tract-location")) {
			Preconditions.checkArgument(cmd.hasOption("portfolio"), "Must supply --portfolio option with --tract-location");
			String locStr = cmd.getOptionValue("tract-location");
			Preconditions.checkState(locStr.contains(","), "--tract-location format should be lat,lon");
			String[] locSplit = locStr.split(",");
			Preconditions.checkState(locSplit.length == 2, "--tract-location format should be lat,lon");
			double lat = Double.parseDouble(locSplit[0]);
			double lon = Double.parseDouble(locSplit[1]);
			Location tractLoc = new Location(lat, lon);
			
			File portfolioFile = new File(cmd.getOptionValue("portfolio"));
			Portfolio portfolio = Portfolio.createPortfolio(portfolioFile);
			
			HashSet<String> assetNames = new HashSet<>();
			if (cmd.hasOption("tract-radius")) {
				double radius = Double.parseDouble(cmd.getOptionValue("tract-radius"));
				for (Asset asset : portfolio.getAssetList()) {
					double dist = LocationUtils.horzDistanceFast(tractLoc, asset.getLocation());
					if (dist <= radius)
						assetNames.add(MPJ_CondLossCalc.getTractName(asset));
				}
				if (rank == 0)
					debug("found "+assetNames.size()+" tracts within "+(float)radius+" km of "+locStr);
			} else {
				double minDist = Double.POSITIVE_INFINITY;
				Asset closest = null;
				for (Asset asset : portfolio.getAssetList()) {
					double dist = LocationUtils.horzDistanceFast(tractLoc, asset.getLocation());
					if (dist <= minDist) {
						minDist = dist;
						closest = asset;
					}
				}
				assetNames.add(MPJ_CondLossCalc.getTractName(closest));
				if (rank == 0)
					debug("closest tract to "+locStr+" is "+(float)minDist+" km away: "+MPJ_CondLossCalc.getTractName(closest));
			}
			tractNames = assetNames.toArray(new String[0]);
		}
		
		if (cmd.hasOption("tract-branch-eals")) {
			File portfolioFile = new File(cmd.getOptionValue("portfolio"));
			Portfolio portfolio = Portfolio.createPortfolio(portfolioFile);
			
			HashSet<String> tractNamesSet = new HashSet<>();
			for (Asset asset : portfolio.getAssetList())
				tractNamesSet.add(MPJ_CondLossCalc.getTractName(asset));
			allTracts = new ArrayList<>(tractNamesSet);
			Collections.sort(allTracts);
			
			// now see if any are already done
			for (int i=allTracts.size(); --i>=0;) {
				File tractFile = new File(resultsDir, allTracts.get(i)+".csv");
				if (tractFile.exists() && tractFile.length() > 0)
					allTracts.remove(i);
			}
			
			if (rank == 0)
				debug("Calculating for "+allTracts.size()+" tracts");
		}
		
		if (tractNames != null || allTracts != null) {
			Preconditions.checkArgument(cmd.hasOption("background-type"), "Must supply --background-type with census tracts");
			BackgroundRupType rupType = BackgroundRupType.valueOf(cmd.getOptionValue("background-type"));
			erf = new FaultSystemSolutionERF(trueMeanSol);
			erf.setCacheGridSources(true); // otherwise crazy slow
			erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.INCLUDE);
			erf.setParameter(BackgroundRupParam.NAME, rupType);
			erf.updateForecast();
		}
		
		if (rank == 0)
			debug("Building EAL branches");
		
		if (cmd.hasOption("lec-cov"))
			covModel = LossCOV_Model.valueOf(cmd.getOptionValue("lec-cov"));
		
		branches = new ArrayList<>();
		HashSet<U3_EAL_ProbModels> probModels = new HashSet<>();
		HashSet<U3_EAL_GMMs> gmms = new HashSet<>();
		HashSet<U3_EAL_GMM_Epistemic> gmmEpis = new HashSet<>();
		HashSet<U3_EAL_Vs30Model> vs30s = new HashSet<>();
		
		int numTI = mappings.keySet().size();
		for (U3_EAL_ProbModels probModel : probsZipFiles.keySet()) {
			for (U3_EAL_GMMs gmm : U3_EAL_GMMs.values()) {
				for (U3_EAL_GMM_Epistemic gmmEpi : U3_EAL_GMM_Epistemic.values()) {
					for (U3_EAL_Vs30Model vs30 : vs30Dirs.keySet()) {
						File vs30Dir = vs30Dirs.get(vs30);
						String prefix = gmm.getShortName();
						if (gmmEpi != U3_EAL_GMM_Epistemic.NONE)
							prefix += "_"+gmmEpi.getShortName();
						
						File fssFile = new File(vs30Dir, prefix+"_fss_index.bin");
						if (!fssFile.exists())
							fssFile = new File(vs30Dir, prefix+"_fss_index.bin.gz");
						if (!fssFile.exists())
							continue;
						
						File griddedFile = new File(vs30Dir, prefix+"_fss_gridded.bin");
						if (!griddedFile.exists())
							griddedFile = new File(vs30Dir, prefix+"_fss_gridded.bin.gz");
						if (!griddedFile.exists())
							griddedFile = null;
						
						File tractDir = new File(vs30Dir, prefix+"_tract_results");
						if (!tractDir.exists())
							tractDir = null;
						
						probModels.add(probModel);
						gmms.add(gmm);
						gmmEpis.add(gmmEpi);
						vs30s.add(vs30);
						for (U3LogicTreeBranch tiBranch : mappings.keySet())
							branches.add(new U3_EAL_LogicTreeBranch(tiBranch, probModel, gmm, gmmEpi, vs30,
									fssFile, griddedFile, tractDir));
					}
				}
			}
		}
		if (rank == 0) {
			debug(numTI+" TI branches");
			debug("Prob Models: "+Joiner.on(", ").join(probModels));
			debug("GMMs: "+Joiner.on(", ").join(gmms));
			debug("GMM Epis: "+Joiner.on(", ").join(gmmEpis));
			debug("Vs30s: "+Joiner.on(", ").join(vs30s));
			int calcNum = numTI*probModels.size()*gmms.size()*gmmEpis.size()*vs30s.size();
			debug("Calculated number if fully specified: "+calcNum+" (fully specified ? "+(calcNum == branches.size())+")");
		}
		Collections.sort(branches, new ReadOptimizedBranchComparator());
		debug("Built "+branches.size()+" branches");
		Preconditions.checkState(!branches.isEmpty(), "No branches found!");
		
		int maxCacheSize = 50;
		if (tractNames != null)
			maxCacheSize = Integer.min(maxCacheSize, tractNames.length);
		rupLossesCache = CacheBuilder.newBuilder().maximumSize(maxCacheSize).build(new CacheLoader<File, double[][]>() {

			@Override
			public double[][] load(File key) throws Exception {
				if (tmpDir != null)
					key = getCachedFile(key);
				return MPJ_CondLossCalc.loadResults(key);
			}
			
		});
		griddedLossesCache = CacheBuilder.newBuilder().maximumSize(maxCacheSize).build(new CacheLoader<File, DiscretizedFunc[]>() {

			@Override
			public DiscretizedFunc[] load(File key) throws Exception {
				if (tmpDir != null)
					key = getCachedFile(key);
				return MPJ_CondLossCalc.loadGridSourcesFile(key,
						trueMeanSol.getGridSourceProvider().getGriddedRegion());
			}
			
		});
		tractCache = CacheBuilder.newBuilder().maximumSize(10).build(new CacheLoader<File, TractLoader>() {

			@Override
			public TractLoader load(File key) throws Exception {
				return new TractLoader(key);
			}
			
		});
		
		exec = Executors.newFixedThreadPool(getNumThreads());
		
		resultsCSV = new CSVFile<>(true);
		U3_EAL_LogicTreeBranch branch0 = branches.get(0);
		List<String> header = Lists.newArrayList("Index", "Branch Weight", "Total EAL", "Fault EAL", "Gridded EAL");
		for (int i=0; i<branch0.size(); i++)
			header.add(branch0.getValue(i).getBranchLevelName());
		resultsCSV.addLine(header);
		resultsFile = new File(resultsDir, "results_"+rank+".csv");
		
		if (tractNames != null) {
			// do the tract combines first in parallel
			List<File> tractDirs = new ArrayList<>();
			for (U3_EAL_GMMs gmm : U3_EAL_GMMs.values()) {
				for (U3_EAL_GMM_Epistemic gmmEpi : U3_EAL_GMM_Epistemic.values()) {
					for (U3_EAL_Vs30Model vs30 : vs30Dirs.keySet()) {
						File vs30Dir = vs30Dirs.get(vs30);
						String prefix = gmm.getShortName();
						if (gmmEpi != U3_EAL_GMM_Epistemic.NONE)
							prefix += "_"+gmmEpi.getShortName();
						
						File tractDir = new File(vs30Dir, prefix+"_tract_results");
						tractDirs.add(tractDir);
					}
				}
			}
			tractCacheDir = new File(outputDir, "tract_caches");
			if (rank == 0) {
				debug("caching tracts to "+tractCacheDir.getAbsolutePath());
				Preconditions.checkState(tractCacheDir.exists() || tractCacheDir.mkdir(),
						"Tract cache dir doesn't exist and could not be created: %s", outputDir.getAbsolutePath());
			}
			Collections.sort(tractDirs, new FilePathComparator());
			MPI.COMM_WORLD.Barrier();
			for (int i=0; i<tractDirs.size(); i++) {
				if (i % size == rank) {
					try {
						new TractLoader(tractDirs.get(i));
					} catch (ExecutionException e) {
						throw ExceptionUtils.asRuntimeException(e);
					}
				}
			}
			MPI.COMM_WORLD.Barrier();
			if (rank == 0)
				debug("DONE caching tracts!");
		}
		
		if (cmd.hasOption("node-temp-dir")) {
			tmpDir = new File(cmd.getOptionValue("node-temp-dir"));
			Preconditions.checkState(tmpDir.exists() || tmpDir.mkdir());
		}
		
		if (cmd.hasOption("calc-lec")) {
			int num = (int)Math.round((9d)/0.1d)+1;
			EvenlyDiscretizedFunc logXVals = new EvenlyDiscretizedFunc(0d, num, 0.1);
//			System.out.println(logXVals);
			lecXVals = new ArbitrarilyDiscretizedFunc();
			lecXVals.set(0d, 0d);
			for (Point2D pt : logXVals)
				lecXVals.set(Math.pow(10, pt.getX()), 0d);
			if (rank == 0)
				debug("Calculating LECs with "+lecXVals.size()+" values");
			lecProbLevels = new double[] { 0.01, 0.004, 0.0025, 0.0018, 0.0004 };
			for (double level : lecProbLevels)
				header.add("Loss @ "+(float)level);
			resultsCSV = new CSVFile<>(true);
			resultsCSV.addLine(header);
			
			lecResultsCSV = new CSVFile<>(true);
			header = Lists.newArrayList("Index", "Branch Weight");
			for (int i=0; i<branch0.size(); i++)
				header.add(branch0.getValue(i).getBranchLevelName());
			for (Point2D pt : lecXVals)
				header.add((float)pt.getX()+"");
			lecResultsCSV.addLine(header);
			lecResultsFileCSV = new File(resultsDir, "results_lec_"+rank+".csv");
		}
	}
	
	private synchronized File getCachedFile(File file) throws IOException {
		if (tmpDir == null || file.getParent().contains("tract_results"))
			return file;
		String unique = file.getParentFile().getParentFile().getName()
				+"_"+file.getParentFile().getName()+"_"+file.getAbsolutePath().hashCode();
		String destName = unique+"_"+file.getName();
		if (destName.toLowerCase().endsWith(".gz"))
			destName = destName.replaceAll(".gz", "").replaceAll(".GZ", "");
		File tmpFile = new File(tmpDir, destName);
		if (tmpFile.exists())
			return tmpFile;
		// need to write it
		debug("caching "+file.getAbsolutePath()+" to "+tmpFile.getAbsolutePath());
		if (file.getName().toLowerCase().endsWith(".gz")) {
			// unzip it
			InputStream zIn = MPJ_CondLossCalc.getInputStream(file);
			OutputStream fOut = new FileOutputStream(tmpFile); //IOUtils will buffer
			IOUtils.copy(zIn, fOut);
			zIn.close();
			fOut.close();
		} else {
			Files.copy(file, tmpFile);
		}
		return tmpFile;
	}
	
	private class ReadOptimizedBranchComparator implements Comparator<U3_EAL_LogicTreeBranch> {
		
		List<Class<? extends LogicTreeBranchNode<?>>> sortOrderClasses;
		
		public ReadOptimizedBranchComparator() {
			sortOrderClasses = new ArrayList<>();
			
			// sort by vs30/gmm/gmm epi first, as files are stored based on that. this will dispatch jobs for the same
			// binary files together, meaning more cache hits
			sortOrderClasses.add(U3_EAL_Vs30Model.class);
			sortOrderClasses.add(U3_EAL_GMMs.class);
			sortOrderClasses.add(U3_EAL_GMM_Epistemic.class);
			sortOrderClasses.add(U3_EAL_ProbModels.class);
			sortOrderClasses.addAll(U3LogicTreeBranch.getLogicTreeNodeClasses());
		}

		@Override
		public int compare(U3_EAL_LogicTreeBranch b1, U3_EAL_LogicTreeBranch b2) {
			Preconditions.checkState(b1.size() == sortOrderClasses.size());
			Preconditions.checkState(b2.size() == sortOrderClasses.size());
			for (Class<? extends LogicTreeBranchNode<?>> clazz : sortOrderClasses) {
				LogicTreeBranchNode<?> val = b1.getValueUnchecked(clazz);
				LogicTreeBranchNode<?> oval = b2.getValueUnchecked(clazz);
				int cmp = val.getShortName().compareTo(oval.getShortName());
				if (cmp != 0)
					return cmp;
			}
			return 0;
		}
		
	}
	
	private class CachedGridSourceCFSS extends CompoundFaultSystemSolution {
		
		private LoadingCache<U3LogicTreeBranch, GridSourceProvider> gridProvCache;

		public CachedGridSourceCFSS(ZipFile zip) {
			super(zip);
			
			gridProvCache = CacheBuilder.newBuilder().maximumSize(10).build(new CacheLoader<U3LogicTreeBranch, GridSourceProvider>() {

				@Override
				public GridSourceProvider load(U3LogicTreeBranch branch) throws Exception {
					try {
						return CachedGridSourceCFSS.super.loadGridSourceProviderFile(branch);
					} catch (Exception e) {}
					System.out.println("Building gridProv for "+branch.buildFileName());
					return CachedGridSourceCFSS.this.getSolution(branch).getGridSourceProvider();
				}
				
			});
		}

		@Override
		public GridSourceProvider loadGridSourceProviderFile(U3LogicTreeBranch branch)
				throws DocumentException, IOException {
			try {
				return gridProvCache.get(branch);
			} catch (ExecutionException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
		}
		
	}

	@Override
	protected int getNumTasks() {
		if (consolidateOnly)
			return 1;
		if (allTracts != null)
			return allTracts.size();
		return branches.size();
	}

	@Override
	protected void calculateBatch(int[] batch) throws Exception {
		if (consolidateOnly)
			return;
		
		if (allTracts == null) {
			List<Future<CalcTask>> futures = new ArrayList<>();
			
			// batch is branch indexes
			for (int index : batch)
				futures.add(exec.submit(new CalcTask(index, null)));
			
			debug("Waiting on "+futures.size()+" futures");
			for (Future<CalcTask> future : futures) {
				CalcTask task = null;
				try {
					task = future.get();
				} catch (Exception e) {
					abortAndExit(e);
				}
				List<String> line = new ArrayList<>();
				line.add(task.index+"");
				line.add(task.branch.getAprioriBranchWt()+"");
				line.add(task.totalEAL+"");
				line.add(task.faultEAL+"");
				line.add(task.griddedEAL+"");
				for (LogicTreeBranchNode<?> node : task.branch)
					line.add(node.getShortName());
				DiscretizedFunc lec = task.lec;
				if (lecXVals != null && lecProbLevels != null && lecProbLevels.length > 0) {
					for (double rate : lecProbLevels) {
						double loss;
						if (rate > lec.getMaxY())
							loss = lec.getMinX();
						else if (rate < lec.getMinY())
							loss = lec.getMaxX();
						else
							loss = lec.getFirstInterpolatedX(rate);
						line.add(loss+"");
					}
				}
				resultsCSV.addLine(line);
				if (lec != null) {
					line = new ArrayList<>();
					line.add(task.index+"");
					line.add(task.branch.getAprioriBranchWt()+"");
					for (LogicTreeBranchNode<?> node : task.branch)
						line.add(node.getShortName());
					for (Point2D pt : lec)
						line.add(pt.getY()+"");
					lecResultsCSV.addLine(line);
				}
			}
			debug("finished batch and flushing CSV");
			resultsCSV.writeToFile(resultsFile);
			if (lecResultsCSV != null)
				lecResultsCSV.writeToFile(lecResultsFileCSV);
		} else {
			// tract index, calculate for all branches
			for (int tractIndex : batch) {
				String tractName = allTracts.get(tractIndex);
				
				List<Future<CalcTask>> futures = new ArrayList<>();
				
				for (int index=0; index<branches.size(); index++)
					futures.add(exec.submit(new CalcTask(index, tractName)));
				
				CSVFile<String> csv = new CSVFile<>(true);
				csv.addLine(resultsCSV.getLine(0));
				for (int index=0; index<futures.size(); index++) {
					CalcTask task = null;
					try {
						task = futures.get(index).get();
					} catch (Exception e) {
						abortAndExit(e);
					}
					List<String> line = new ArrayList<>();
					line.add(task.index+"");
					line.add(task.branch.getAprioriBranchWt()+"");
					line.add(task.totalEAL+"");
					line.add(task.faultEAL+"");
					line.add(task.griddedEAL+"");
					for (LogicTreeBranchNode<?> node : task.branch)
						line.add(node.getShortName());
					csv.addLine(line);
				}
				csv.writeToFile(new File(resultsDir, tractName+".csv"));
			}
		}
	}
	
	private class TractLoadRunnable implements Runnable {
		
		private File tractFile;
		private double[][] totTractLosses;

		public TractLoadRunnable(File tractFile, double[][] totTractLosses) {
			this.tractFile = tractFile;
			this.totTractLosses = totTractLosses;
		}

		@Override
		public void run() {
			double[][] tractLosses;
			try {
				tractLosses = MPJ_CondLossCalc.loadResults(tractFile);
			} catch (IOException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
//			try {
//				tractLosses = rupLossesCache.get(tractFile);
//			} catch (ExecutionException e) {
//				throw ExceptionUtils.asRuntimeException(e);
//			}
			synchronized (totTractLosses) {
				// add them
				Preconditions.checkState(totTractLosses.length == tractLosses.length);
				for (int i=0; i<tractLosses.length; i++) {
					if (tractLosses[i] != null) {
						if (totTractLosses[i] == null) {
							totTractLosses[i] = Arrays.copyOf(tractLosses[i], tractLosses[i].length);
						} else {
							Preconditions.checkState(tractLosses[i].length == totTractLosses[i].length);
							for (int j=0; j<tractLosses[i].length; j++)
								totTractLosses[i][j] += tractLosses[i][j];
						}
					}
				}
			}
		}
		
	}
	
	private class TractLoader {
		double[][] fssLosses;
		DiscretizedFunc[] griddedLosses;

		public TractLoader(File tractDir) throws ExecutionException, IOException {
			// first sum all tract losses, erf indexed
			double[][] totTractLosses = null;
			
			String unique = tractDir.getParentFile().getName()
					+"_"+tractDir.getName()+"_"+tractDir.getAbsolutePath().hashCode();
			unique += "_"+tractNames.length+"tracts";
			File consolidatedFile = new File(tractCacheDir, unique+".bin");
			synchronized (MPJ_UCERF3_EAL_Combiner.class) {
				if (consolidatedFile.exists())
					totTractLosses = rupLossesCache.get(consolidatedFile);
				else
					debug("Need to cache tract results for "+tractDir.getAbsolutePath());
			}
			
			if (totTractLosses == null) {
				List<Future<?>> futures = new ArrayList<>();
				totTractLosses = new double[erf.getNumSources()][];
				for (String tractName : tractNames) {
					File tractFile = new File(tractDir, tractName+".bin");
					if (!tractFile.exists())
						tractFile = new File(tractDir, tractName+".bin.gz");
					Preconditions.checkState(tractFile.exists(), "Tract file doesn't exist: %s", tractFile.getAbsolutePath());
					
					futures.add(exec.submit(new TractLoadRunnable(tractFile, totTractLosses)));
				}
				
				for (int i=0; i<futures.size(); i++) {
					try {
						futures.get(i).get();
					} catch (InterruptedException e) {
						throw ExceptionUtils.asRuntimeException(e);
					}
					if (i % 100 == 0 || i == futures.size()-1)
						debug("Processed tract "+i+"/"+futures.size());
				}
			}
			
			synchronized (MPJ_UCERF3_EAL_Combiner.class) {
				if (consolidatedFile != null && !consolidatedFile.exists()) {
					debug("Caching tract results from "+tractDir.getAbsolutePath()+" to "+consolidatedFile.getAbsolutePath());
					File writeFile = new File(consolidatedFile.getAbsolutePath()+".tmp."+rank);
					MPJ_CondLossCalc.writeResults(writeFile, totTractLosses);
					Files.move(writeFile, consolidatedFile);
				}
			}
			
			// then convert the summed to FSS and gridded
			fssLosses = MPJ_CondLossCalc.mapResultsToFSS(erf, totTractLosses);
			griddedLosses = MPJ_CondLossCalc.mapResultsToGridded(erf, totTractLosses);
		}
	}
	
	private class CalcTask implements Callable<CalcTask> {
		
		private int index;
		private U3_EAL_LogicTreeBranch branch;
		private String tractName;
		
		private double faultEAL;
		private double griddedEAL;
		private double totalEAL;
		
		private DiscretizedFunc lec;

		public CalcTask(int index, String tractName) {
			this.index = index;
			this.branch = branches.get(index);
			this.tractName = tractName;
		}

		@Override
		public CalcTask call() throws Exception {
			Map<U3LogicTreeBranch, List<Integer>> taskMappings = new HashMap<>();
			taskMappings.put(branch.getTIBranch(), mappings.get(branch.getTIBranch()));
			
			double[][] fssLosses = null;
			DiscretizedFunc[] griddedLosses = null;
			if (tractName != null) {
				File tractDir = branch.getTractDir();
				File tractFile = new File(tractDir, tractName+".bin");
				if (!tractFile.exists())
					tractFile = new File(tractDir, tractName+".bin.gz");
				
				double[][] totTractLosses = rupLossesCache.get(tractFile);
				
				fssLosses = MPJ_CondLossCalc.mapResultsToFSS(erf, totTractLosses);
				griddedLosses = MPJ_CondLossCalc.mapResultsToGridded(erf, totTractLosses);
			} else if (tractNames == null) {
				fssLosses = rupLossesCache.get(branch.getFSSIndexedBinFile());
				griddedLosses = griddedLossesCache.get(branch.getGriddedBinFile());
			} else {
				TractLoader tractLoader = tractCache.get(branch.getTractDir());
				fssLosses = tractLoader.fssLosses;
				griddedLosses = tractLoader.griddedLosses;
			}
			
			ZipFile erfProbsZipFile = probsZipFiles.get(branch.getValue(U3_EAL_ProbModels.class));
			
			UCERF3_EAL_Combiner calc = new UCERF3_EAL_Combiner(cfss, taskMappings, trueMeanSol, fssLosses, griddedLosses,
					erfProbsZipFile, erfProbsDuration, lecXVals, covModel);
			
			faultEAL = calc.getFaultEALs()[0];
			griddedEAL = calc.getGriddedEALs()[0];
			totalEAL = calc.getTotalEALs()[0];
			
			if (lecXVals != null)
				lec = calc.getLECs()[0];
			
			return this;
		}
		
	}
	
	private class BranchNodeVals {
		private double weightTotal;
		private double meanLoss;
		private double meanFault;
		private double meanGridded;
	}

	@Override
	protected void doFinalAssembly() throws Exception {
		exec.shutdown();
		rupLossesCache.invalidateAll();
		griddedLossesCache.invalidateAll();
		if (rank == 0 && allTracts == null) {
			debug("Consolidating CSVs");
			List<List<String>> allLines = new ArrayList<>();
			for (int i=0; i<branches.size(); i++)
				allLines.add(null);
			List<List<String>> allLECLines = null;
			if (lecXVals != null) {
				allLECLines = new ArrayList<>();
				for (int i=0; i<branches.size(); i++)
					allLECLines.add(null);
			}
			int loaded = 0;
			List<Map<LogicTreeBranchNode<?>, BranchNodeVals>> nodeVals = new ArrayList<>();
			for (int i=0; i<branches.get(0).size(); i++)
				nodeVals.add(new HashMap<>());
			double totalWeight = 0d;
			double totalMean = 0d;
			double faultMean = 0d;
			double griddedMean = 0d;
			for (int r=0; r<size; r++) {
				CSVFile<String> csv;
				CSVFile<String> lecCSV = null;
				if (r == 0 && !consolidateOnly) {
					csv = this.resultsCSV;
					lecCSV = this.lecResultsCSV;
				} else {
					File resultsFile = new File(resultsDir, "results_"+r+".csv");
					if (!resultsFile.exists()) {
						debug("No results for rank "+r);
						continue;
					}
					csv = CSVFile.readFile(resultsFile, true);
					if (lecXVals != null) {
						resultsFile = new File(resultsDir, "results_lec_"+r+".csv");
						lecCSV = CSVFile.readFile(resultsFile, true);
					}
				}
				for (int row=1; row<csv.getNumRows(); row++) {
					loaded++;
					int index = Integer.parseInt(csv.get(row, 0));
					Preconditions.checkState(allLines.get(index) == null, "Duplicate found for index "+index+" in rank "+r);
					allLines.set(index, csv.getLine(row));
					if (lecCSV != null)
						allLECLines.set(index, lecCSV.getLine(row));
					U3_EAL_LogicTreeBranch branch = branches.get(index);
					for (int i=0; i<branch.size(); i++) {
//						getBranchLevelName()
						String choice = branch.getValue(i).getShortName();
						String testChoice = csv.get(row, i+5);
						if (!choice.equals(testChoice)) {
							System.err.println("Branch mismatch for rank "+r+", index "+index);
							System.err.println("\tOriginal Branch: "+branch);
							List<String> line = csv.getLine(row);
							System.err.println("\tCSV Branch: "+Joiner.on(", ").join(line.subList(5, line.size())));
							System.err.flush();
							throw new IllegalStateException("Branch mismatch for rank "+r+", index "+index);
						}
					}
					double weight = branch.getAprioriBranchWt();
					double totEAL = Double.parseDouble(csv.get(row, 2));
					double faultEAL = Double.parseDouble(csv.get(row, 3));
					double griddedEAL = Double.parseDouble(csv.get(row, 4));
					totalWeight += weight;
					totalMean += totEAL*weight;
					faultMean += faultEAL*weight;
					griddedMean += griddedEAL*weight;
					for (int i=0; i<branch.size(); i++) {
						BranchNodeVals vals = nodeVals.get(i).get(branch.getValue(i));
						LogicTreeBranchNode<?> node = branch.getValue(i);
						if (vals == null) {
							vals = new BranchNodeVals();
							nodeVals.get(i).put(node, vals);
						}
						vals.weightTotal += weight;
						vals.meanLoss += totEAL*weight;
						vals.meanFault += faultEAL*weight;
						vals.meanGridded += griddedEAL*weight;
					}
				}
			}
			debug("Loaded "+loaded+" branches");
			debug("Total weight: "+totalWeight);
			Preconditions.checkState(loaded == branches.size(), "Did not load all branches. Expected %s, loaded %s", branches.size(), loaded);
			debug("Writing master CSV");
			CSVFile<String> csv = new CSVFile<>(true);
			csv.addLine(resultsCSV.getLine(0)); // header
			csv.addAll(allLines);
			csv.writeToFile(new File(outputDir, "all_branch_results.csv"));
			if (lecXVals != null) {
				csv = new CSVFile<>(true);
				csv.addLine(lecResultsCSV.getLine(0)); // header
				csv.addAll(allLECLines);
				csv.writeToFile(new File(outputDir, "all_branch_lec_results.csv"));
			}
			debug("Writing branch levels CSV");
			csv = new CSVFile<>(true);
			csv.addLine("Branch Level", "Branch Choice", "Total Weight", "Weighted Total Mean EAL",
					"Weighted Fault Mean EAL", "Weighted Gridded Mean EAL");
			Comparator<LogicTreeBranchNode<?>> nodeComparator = new Comparator<LogicTreeBranchNode<?>>() {

				@Override
				public int compare(LogicTreeBranchNode<?> o1, LogicTreeBranchNode<?> o2) {
					return o1.getShortName().compareTo(o2.getShortName());
				}

			};
			
			for (int i=0; i<nodeVals.size(); i++) {
				Map<LogicTreeBranchNode<?>, BranchNodeVals> valsMap = nodeVals.get(i);
				List<LogicTreeBranchNode<?>> choices = new ArrayList<>(valsMap.keySet());
				choices.sort(nodeComparator);
				for (LogicTreeBranchNode<?> choice : choices) {
					List<String> line = new ArrayList<>();
					line.add(choice.getBranchLevelName());
					line.add(choice.getShortName());
					BranchNodeVals vals = valsMap.get(choice);
					// normalize
					vals.meanLoss /= vals.weightTotal;
					vals.meanFault /= vals.weightTotal;
					vals.meanGridded /= vals.weightTotal;
					line.add((float)(vals.weightTotal/totalWeight)+"");
					line.add(vals.meanLoss+"");
					line.add(vals.meanFault+"");
					line.add(vals.meanGridded+"");
					csv.addLine(line);
				}
			}
			// normalize totals
			totalMean /= totalWeight;
			faultMean /= totalWeight;
			griddedMean /= totalWeight;
			csv.addLine("COMPLETE MODEL", "MEAN", "1.0", totalMean+"", faultMean+"", griddedMean+"");
			csv.writeToFile(new File(outputDir, "branch_level_summary.csv"));
		}
	}
	
	public static Options createOptions() {
		Options options = MPJTaskCalculator.createOptions();
		
		Option willsDir = new Option("wills", "wills-dir", true, "Directory containing Wills 2015 results");
		willsDir.setRequired(false);
		options.addOption(willsDir);
		
		Option waldDir = new Option("wald", "wald-dir", true, "Directory containing Wald & Allen results");
		waldDir.setRequired(false);
		options.addOption(waldDir);
		
		Option erfProbsDir = new Option("probs", "erf-probs-dir", true, "ERF probabiltiies directory");
		erfProbsDir.setRequired(true);
		options.addOption(erfProbsDir);
		
		Option erfProbsDuration = new Option("dur", "erf-probs-duration", true, "ERF probabiltiies duration");
		erfProbsDuration.setRequired(true);
		options.addOption(erfProbsDuration);
		
		Option trueMeanSol = new Option("tms", "true-mean-sol", true, "True mean solution file (with mappings)");
		trueMeanSol.setRequired(true);
		options.addOption(trueMeanSol);
		
		Option compoundSol = new Option("cfss", "compound-sol", true, "Compound FSS File");
		compoundSol.setRequired(true);
		options.addOption(compoundSol);
		
		Option consolidateOnly = new Option("co", "consolidate-only", false, "Flag to consolidate only");
		consolidateOnly.setRequired(false);
		options.addOption(consolidateOnly);
		
		Option tract = new Option("tr", "tract", true, "Census tract by name (or comma separated tracts)");
		tract.setRequired(false);
		options.addOption(tract);
		
		Option tractLoc = new Option("trl", "tract-location", true, "Census tract location (lat,lon). "
				+ "Must be used with --portfolio. If --trace-radius supplied, then all tracts within the given "
				+ "radius will be included. Otherwise only the closest");
		tractLoc.setRequired(false);
		options.addOption(tractLoc);
		
		Option tractRadius = new Option("trr", "tract-radius", true, "Census tract radius (km) from tract location");
		tractRadius.setRequired(false);
		options.addOption(tractRadius);
		
		Option tractBranchEALs = new Option("trbe", "tract-branch-eals", false,
				"Flag to write out full branch EALs CSV files for each tract");
		tractBranchEALs.setRequired(false);
		options.addOption(tractBranchEALs);
		
		Option portfolio = new Option("p", "portfolio", true,
				"Portfolio file, used when searching for census tracts with --tract-location");
		portfolio.setRequired(false);
		options.addOption(portfolio);
		
		Option backType = new Option("bgt", "background-type", true,
				"Background rup type, required with --tract or --tract-location to parse tract files");
		backType.setRequired(false);
		options.addOption(backType);
		
		Option tmpDir = new Option("ntd", "node-temp-dir", true,
				"Node-local temp dir to cache input files");
		tmpDir.setRequired(false);
		options.addOption(tmpDir);
		
		Option doLEC = new Option("lec", "calc-lec", false,
				"Flag to calculate LEC values");
		doLEC.setRequired(false);
		options.addOption(doLEC);
		
		Option cov = new Option("cov", "lec-cov", true,
				"Loss COV model for LECs");
		cov.setRequired(false);
		options.addOption(cov);
		
		return options;
	}

	public static void main(String[] args) {
		System.setProperty("java.awt.headless", "true");
		try {
			args = MPJTaskCalculator.initMPJ(args);
			
			Options options = createOptions();
			
			CommandLine cmd = parse(options, args, MPJ_UCERF3_EAL_Combiner.class);
			
			args = cmd.getArgs();
			
			if (args.length != 1) {
				System.err.println("USAGE: "+ClassUtils.getClassNameWithoutPackage(MPJ_UCERF3_EAL_Combiner.class)
						+" <output-dir>");
				abortAndExit(2);
			}
			
			File outputDir = new File(args[0]);
			
			MPJ_UCERF3_EAL_Combiner driver = new MPJ_UCERF3_EAL_Combiner(cmd, outputDir);
			
			driver.run();
			
			finalizeMPJ();
			
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}

}
