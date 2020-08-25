package scratch.kevin.ucerf3.eal.spatialCorr;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.stat.StatUtils;
import org.dom4j.Document;
import org.dom4j.DocumentException;
import org.dom4j.Element;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.util.ClassUtils;
import org.opensha.commons.util.FileNameComparator;
import org.opensha.commons.util.XMLUtils;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.imr.AbstractIMR;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sra.calc.parallel.MPJ_CondLossCalc;
import org.opensha.sra.gui.portfolioeal.Asset;
import org.opensha.sra.gui.portfolioeal.Portfolio;
import org.opensha.sra.vulnerability.VulnerabilityFetcher;

import com.google.common.base.Preconditions;
import com.google.common.collect.Range;
import com.google.common.collect.Table;

import edu.usc.kmilner.mpj.taskDispatch.MPJTaskCalculator;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.kevin.ucerf3.eal.branches.U3_EAL_GMM_Epistemic;
import scratch.kevin.ucerf3.eal.branches.U3_EAL_GMMs;
import scratch.kevin.ucerf3.eal.branches.U3_EAL_LogicTreeBranch;
import scratch.kevin.ucerf3.eal.branches.U3_EAL_Vs30Model;

public class MPJ_SpatiallyCorrelatedLossCalc extends MPJTaskCalculator {
	
	private FaultSystemSolution trueMeanSol;
	private FaultSystemSolutionERF erf;
	
	private Map<U3_EAL_Vs30Model, File> vs30Dirs;

	private File resultsDir;
	
	private List<U3_EAL_LogicTreeBranch> branches = new ArrayList<>();
	
	private DiscretizedFunc lecXVals;
	
	private ExecutorService exec;
	
	private List<Asset> assets;
	private RandomFieldLoader[] fields;
	
	// for preset between-event standard deviations
	private double[] betweenEventStdDevs;
	// for random tau samples
	private int randTaus;
	private NormalDistribution normDist;
	
	public MPJ_SpatiallyCorrelatedLossCalc(CommandLine cmd, File outputDir) throws IOException, DocumentException {
		super(cmd);
		this.shuffle = false;
		
		File trueMeanSolFile = new File(cmd.getOptionValue("true-mean-sol"));
		
		if (rank == 0)
			Preconditions.checkState(outputDir.exists() || outputDir.mkdir(),
					"Output dir doesn't exist and could not be created: %s", outputDir.getAbsolutePath());
		resultsDir = new File(outputDir, "results");
		
		if (rank == 0)
			debug("Loading true mean solution from: "+trueMeanSolFile.getAbsolutePath());
		trueMeanSol = FaultSystemIO.loadSol(trueMeanSolFile);
		
		if (rank == 0)
			Preconditions.checkState(resultsDir.exists() || resultsDir.mkdir(),
					"Results dir doesn't exist and could not be created: %s", resultsDir.getAbsolutePath());
		
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
		
		LogicTreeBranch emptyTIBranch = (LogicTreeBranch)LogicTreeBranch.DEFAULT.clone();
		for (int i=0; i<emptyTIBranch.size(); i++)
			emptyTIBranch.clearValue(i);
		
		branches = new ArrayList<>();
		for (U3_EAL_GMMs gmm : U3_EAL_GMMs.values()) {
			if (gmm == U3_EAL_GMMs.IDRISS_2014)
				// TODO: Idriss
				continue;
			for (U3_EAL_GMM_Epistemic gmmEpi : U3_EAL_GMM_Epistemic.values()) {
				for (U3_EAL_Vs30Model vs30 : vs30Dirs.keySet()) {
					File vs30Dir = vs30Dirs.get(vs30);
					String prefix = gmm.getShortName();
					if (gmmEpi != U3_EAL_GMM_Epistemic.NONE)
						prefix += "_"+gmmEpi.getShortName();
					
					File binFile = new File(vs30Dir, prefix+".bin");
					if (!binFile.exists())
						binFile = new File(vs30Dir, prefix+".bin.gz");
					if (!binFile.exists())
						continue;
					
					branches.add(new U3_EAL_LogicTreeBranch(
							emptyTIBranch, null, gmm, gmmEpi, vs30, binFile, null, null));
				}
			}
		}
		Collections.sort(branches);
		
		int num = (int)Math.round((9d)/0.1d)+1;
		EvenlyDiscretizedFunc logXVals = new EvenlyDiscretizedFunc(0d, num, 0.1);
//		System.out.println(logXVals);
		lecXVals = new ArbitrarilyDiscretizedFunc();
		lecXVals.set(0d, 0d);
		for (Point2D pt : logXVals)
			lecXVals.set(Math.pow(10, pt.getX()), 0d);
		if (rank == 0)
			debug("Have "+lecXVals.size()+" x-values");
		
		File portfolioFile = new File(cmd.getOptionValue("portfolio"));
		Portfolio portfolio = Portfolio.createPortfolio(portfolioFile);
		assets = portfolio.getAssetList();
		
		File vulnFile = new File(cmd.getOptionValue("vuln-file"));
		System.out.println("trying to load vulnerabilities from: "+vulnFile.getAbsolutePath());
		VulnerabilityFetcher.getVulnerabilities(vulnFile);
		System.out.println("DONE loading vulns.");
		
		if (cmd.hasOption("taus")) {
			Preconditions.checkState(!cmd.hasOption("rand-taus"),
					"cannot supply both --taus and --rand-taus options");
			String tauStr = cmd.getOptionValue("taus");
			String[] tauSplit = tauStr.split(",");
			betweenEventStdDevs = new double[tauSplit.length];
			for (int i=0; i<betweenEventStdDevs.length; i++)
				betweenEventStdDevs[i] = Double.parseDouble(tauSplit[i]);
		} else {
			Preconditions.checkState(cmd.hasOption("rand-taus"),
					"must supply --taus or --rand-taus options");
			randTaus = Integer.parseInt(cmd.getOptionValue("rand-taus"));
			Well19937c rng = new Well19937c(System.nanoTime()*(rank+1));
			normDist = new NormalDistribution(rng, 0d, 1d);
		}
		
		File fieldsDir = new File(cmd.getOptionValue("fields-dir"));
		Preconditions.checkArgument(fieldsDir.exists());
		
		File[] fieldFiles = fieldsDir.listFiles();
		Arrays.sort(fieldFiles, new FileNameComparator());
		
		double gridSpacing = Double.parseDouble(cmd.getOptionValue("field-spacing"));
		
		List<RandomFieldLoader> randFields = new ArrayList<>();
		for (File file : fieldFiles)
			if (file.getName().endsWith(".csv"))
				randFields.add(RandomFieldLoader.load(file, gridSpacing));
		if (rank == 0)
			debug("Loaded "+randFields.size()+" random fields");
		fields = randFields.toArray(new RandomFieldLoader[0]);
		
		exec = Executors.newFixedThreadPool(getNumThreads());
	}

	@Override
	protected int getNumTasks() {
		return branches.size();
	}

	@Override
	protected void calculateBatch(int[] batch) throws Exception {
		for (int index : batch) {
			U3_EAL_LogicTreeBranch branch = branches.get(index);
			File binFile = branch.getFSSIndexedBinFile(); // not fss index, but full
			double[][] results = MPJ_CondLossCalc.loadResults(binFile);
			
			List<Future<CalcCallable>> futures = new ArrayList<>();
			
			U3_EAL_GMMs gmmBranch = branch.getValue(U3_EAL_GMMs.class);
			U3_EAL_GMM_Epistemic gmmEpiBranch = branch.getValue(U3_EAL_GMM_Epistemic.class);
			U3_EAL_Vs30Model vs30Branch = branch.getValue(U3_EAL_Vs30Model.class);
			
			String dataPrefix = binFile.getAbsolutePath();
			dataPrefix = dataPrefix.substring(0, dataPrefix.indexOf(".bin"));
			File xmlFile = new File(dataPrefix+".xml");
			Preconditions.checkState(xmlFile.exists(),
					"GMPE xml file doesn't exist: %s", xmlFile.getAbsolutePath());
			
			Document doc = XMLUtils.loadDocument(xmlFile);
			Element root = doc.getRootElement();
			
			if (erf == null) {
				// load the ERF
				erf = (FaultSystemSolutionERF) AbstractERF.fromXMLMetadata(root.element(AbstractERF.XML_METADATA_NAME));
				erf.setSolution(trueMeanSol);
				erf.updateForecast();
			}
			
			Preconditions.checkState(results.length == erf.getNumSources(),
					"Source count mismatch: %s != %s", results.length, erf.getNumSources());
			
			double erfDuration = erf.getTimeSpan().getDuration();
			
			debug("Calculating "+index);
			
			for (int i=0; i<lecXVals.size(); i++) {
				Range<Double> lossRange;
				double lossXVal = lecXVals.getX(i);
				if (i == lecXVals.size()-1)
					lossRange = Range.closedOpen(lossXVal, Double.POSITIVE_INFINITY);
				else
					lossRange = Range.closedOpen(lossXVal, lecXVals.getX(i+1));
				
				ScalarIMR gmpe = (ScalarIMR)AbstractIMR.fromXMLMetadata(
						root.element(AbstractIMR.XML_METADATA_NAME), null);
				
				double maxRate = 0d;
				ProbEqkRupture rup = null;
				double meanLossAtMaxRate = Double.NaN;
				
				for (int sourceID=0; sourceID<results.length; sourceID++) {
					if (results[sourceID] == null)
						continue;
					ProbEqkSource source = erf.getSource(sourceID);
					Preconditions.checkState(results[sourceID].length == source.getNumRuptures(),
							"Rupture count mismatch for source %s: %s != %s",
							sourceID, results[sourceID].length, source.getNumRuptures());
					
					for (int rupID=0; rupID<source.getNumRuptures(); rupID++) {
						if (lossRange.contains(results[sourceID][rupID])) {
							ProbEqkRupture testRup = source.getRupture(rupID);
							double rate = testRup.getMeanAnnualRate(erfDuration);
							if (rate > maxRate) {
								maxRate = rate;
								rup = testRup;
								meanLossAtMaxRate = results[sourceID][rupID];
							}
						}
					}
				}
				
				if (rup != null) {
					double[] betweenEventStdDevs;
					if (this.betweenEventStdDevs == null)
						// randomly sample between-event standard deviations
						betweenEventStdDevs = normDist.sample(randTaus);
					else
						betweenEventStdDevs = this.betweenEventStdDevs;
					futures.add(exec.submit(new CalcCallable(gmpe, rup, lossXVal,
							meanLossAtMaxRate, betweenEventStdDevs)));
				}
			}
			
			debug("Waiting on "+futures.size()+" futures for "+index);
			
			List<CalcCallable> calls = new ArrayList<>();
			for (Future<CalcCallable> future : futures) {
				try {
					calls.add(future.get());
				} catch (Exception e) {
					abortAndExit(e, 1);
				}
			}
			
			if (betweenEventStdDevs == null) {
				// random samples
				CSVFile<String> fullCSV = new CSVFile<>(true);
				List<String> header = new ArrayList<>();
				header.add("GMPE");
				header.add("GMPE Epistempic Branch");
				header.add("Vs30 Model");
				header.add("Loss Bin");
				header.add("Modal Rupture Rate");
				header.add("Modal Rupture Mag");
				header.add("Modal Rupture Mean Loss");
				header.add("Between-Event Index");
				header.add("Between-Event Term");
				for (int f=0; f<fields.length; f++)
					header.add("Loss for Field "+f);
				fullCSV.addLine(header);
				CSVFile<String> summaryCSV = new CSVFile<>(true);
				header = new ArrayList<>();
				header.add("GMPE");
				header.add("GMPE Epistempic Branch");
				header.add("Vs30 Model");
				header.add("Loss Bin");
				header.add("Modal Rupture Rate");
				header.add("Modal Rupture Mag");
				header.add("Modal Rupture Mean Loss");
				header.add("Ln-Mean Calculated Loss");
				header.add("Ln Loss Standard Deviation");
				summaryCSV.addLine(header);
				
				for (CalcCallable call : calls) {
					if (call.rup == null)
						continue;
					
					double[] logLosses = new double[call.betweenEventStdDevs.length*fields.length];
					int lossIndex = 0;
					for (int t=0; t<call.betweenEventStdDevs.length; t++) {
						double between = call.betweenEventStdDevs[t];
						List<String> line = new ArrayList<>();
						line.add(gmmBranch.getShortName());
						line.add(gmmEpiBranch.getShortName());
						line.add(vs30Branch.getShortName());
						line.add((float)call.lossXVal+"");
						line.add(call.rup.getMeanAnnualRate(erfDuration)+"");
						line.add(call.rup.getMag()+"");
						line.add(call.rupMeanLoss+"");
						line.add(t+"");
						line.add((float)between+"");
						for (int i=0; i<fields.length; i++) {
							double val = call.lossTable.get(between, fields[i]);
							Preconditions.checkState(Double.isFinite(val));
							line.add(val+"");
							logLosses[lossIndex++] = Math.log(val);
						}
						fullCSV.addLine(line);
					}
					Preconditions.checkState(logLosses.length == lossIndex);
					List<String> line = new ArrayList<>();
					line.add(gmmBranch.getShortName());
					line.add(gmmEpiBranch.getShortName());
					line.add(vs30Branch.getShortName());
					line.add((float)call.lossXVal+"");
					line.add(call.rup.getMeanAnnualRate(erfDuration)+"");
					line.add(call.rup.getMag()+"");
					line.add(call.rupMeanLoss+"");
					double lnMean = StatUtils.mean(logLosses);
					double std = Math.sqrt(StatUtils.variance(logLosses));
					line.add(lnMean+"");
					line.add(std+"");
					summaryCSV.addLine(line);
				}
				
				String outputPrefix = gmmBranch.encodeChoiceString()+"_"+gmmEpiBranch.encodeChoiceString()
					+"_"+vs30Branch.encodeChoiceString();
				String outputName = outputPrefix+".csv";
				debug("Writing "+outputName);
				fullCSV.writeToFile(new File(resultsDir, outputName));

				outputName = outputPrefix+"_summary.csv";
				debug("Writing "+outputName);
				summaryCSV.writeToFile(new File(resultsDir, outputName));
			} else {
				CSVFile<String> csv = new CSVFile<>(true);
				List<String> header = new ArrayList<>();
				header.add("GMPE");
				header.add("GMPE Epistempic Branch");
				header.add("Vs30 Model");
				header.add("Loss Bin");
				header.add("Modal Rupture Rate");
				header.add("Modal Rupture Mag");
				header.add("Modal Rupture Mean Loss");
				header.add("Modal Rupture Between-Event Term");
				header.add("Field ID");
				header.add("Loss");
				csv.addLine(header);
				
				for (CalcCallable call : calls) {
					if (call.rup == null)
						continue;
					
					for (double between : betweenEventStdDevs) {
						List<String> linePrefix = new ArrayList<>();
						linePrefix.add(gmmBranch.getShortName());
						linePrefix.add(gmmEpiBranch.getShortName());
						linePrefix.add(vs30Branch.getShortName());
						linePrefix.add((float)call.lossXVal+"");
						linePrefix.add(call.rup.getMeanAnnualRate(erfDuration)+"");
						linePrefix.add(call.rup.getMag()+"");
						linePrefix.add(call.rupMeanLoss+"");
						linePrefix.add((float)between+"");
						for (int i=0; i<fields.length; i++) {
							List<String> line = new ArrayList<>(linePrefix);
							line.add(i+"");
							line.add(call.lossTable.get(between, fields[i]).toString());
							csv.addLine(line);
						}
					}
				}
				
				String outputName = gmmBranch.encodeChoiceString()+"_"+gmmEpiBranch.encodeChoiceString()
					+"_"+vs30Branch.encodeChoiceString()+".csv";
				debug("Writing "+outputName);
				csv.writeToFile(new File(resultsDir, outputName));
			}
			
			debug("Done with "+index);
		}
	}
	
	private class CalcCallable implements Callable<CalcCallable> {
		
		private ScalarIMR gmpe;
		private ProbEqkRupture rup;
		private double lossXVal;
		private double rupMeanLoss;
		private double[] betweenEventStdDevs;
		
		private Table<Double, RandomFieldLoader, Double> lossTable;

		public CalcCallable(ScalarIMR gmpe, ProbEqkRupture rup, double lossXVal,
				double rupMeanLoss, double[] betweenEventStdDevs) {
			super();
			this.gmpe = gmpe;
			this.rup = rup;
			this.lossXVal = lossXVal;
			this.rupMeanLoss = rupMeanLoss;
			this.betweenEventStdDevs = betweenEventStdDevs;
		}

		@Override
		public CalcCallable call() throws Exception {
			if (rup == null)
				// no rupture in this bin
				return this;
			
			Location rupCentroid = SpatiallyCorrelatedLossCalc.calcRupCentroid(rup.getRuptureSurface());
			lossTable =	SpatiallyCorrelatedLossCalc.calcSpatiallyCorrelatedLoss(
							gmpe, assets, rup, rupCentroid, betweenEventStdDevs, fields);
			
			return this;
		}
		
	}

	@Override
	protected void doFinalAssembly() throws Exception {
		exec.shutdown();
		
		long totalInterps = 0;
		long totalTilings = 0;
		for (RandomFieldLoader field : fields) {
			totalInterps += field.getTotalNumCalcs();
			totalTilings += field.getNumWrappedCalcs();
		}
		
		double percent = 100d*(double)totalTilings/(double)totalInterps;
		
		debug(totalTilings+"/"+totalInterps+" ("+(float)percent+" %) calculations were wrapped/tiled");
	}
	
	public static Options createOptions() {
		Options options = MPJTaskCalculator.createOptions();
		
		Option vulnOp = new Option("v", "vuln-file", true, "VUL06 file");
		vulnOp.setRequired(true);
		options.addOption(vulnOp);
		
		Option willsDir = new Option("wills", "wills-dir", true, "Directory containing Wills 2015 results");
		willsDir.setRequired(false);
		options.addOption(willsDir);
		
		Option waldDir = new Option("wald", "wald-dir", true, "Directory containing Wald & Allen results");
		waldDir.setRequired(false);
		options.addOption(waldDir);
		
		Option trueMeanSol = new Option("tms", "true-mean-sol", true, "True mean solution file (with mappings)");
		trueMeanSol.setRequired(true);
		options.addOption(trueMeanSol);
		
		Option portfolio = new Option("p", "portfolio", true,
				"Portfolio file");
		portfolio.setRequired(true);
		options.addOption(portfolio);
		
		Option betweenStdDevs = new Option("ts", "taus", true,
				"Between-event standard deviations (comma separated)");
		betweenStdDevs.setRequired(false);
		options.addOption(betweenStdDevs);
		
		Option randBetweenStdDevs = new Option("rts", "rand-taus", true,
				"Number of randomly-sampled between-event standard deviations");
		randBetweenStdDevs.setRequired(false);
		options.addOption(randBetweenStdDevs);
		
		Option fieldsDir = new Option("fd", "fields-dir", true,
				"Directory containing random field CSV files");
		fieldsDir.setRequired(true);
		options.addOption(fieldsDir);
		
		Option fieldSpacing = new Option("fs", "field-spacing", true,
				"Field grid spacing in km");
		fieldSpacing.setRequired(true);
		options.addOption(fieldSpacing);
		
		return options;
	}

	public static void main(String[] args) {
		System.setProperty("java.awt.headless", "true");
		try {
			args = MPJTaskCalculator.initMPJ(args);
			
			Options options = createOptions();
			
			CommandLine cmd = parse(options, args, MPJ_SpatiallyCorrelatedLossCalc.class);
			
			args = cmd.getArgs();
			
			if (args.length != 1) {
				System.err.println("USAGE: "+ClassUtils.getClassNameWithoutPackage(MPJ_SpatiallyCorrelatedLossCalc.class)
						+" <output-dir>");
				abortAndExit(2);
			}
			
			File outputDir = new File(args[0]);
			
			MPJ_SpatiallyCorrelatedLossCalc driver = new MPJ_SpatiallyCorrelatedLossCalc(cmd, outputDir);
			
			driver.run();
			
			finalizeMPJ();
			
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}

}
