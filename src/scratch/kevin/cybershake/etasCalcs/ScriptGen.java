package scratch.kevin.cybershake.etasCalcs;

import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import java.util.List;

import org.dom4j.DocumentException;
import org.opensha.commons.geo.Location;
import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter;
import org.opensha.commons.hpc.mpj.MPJExpressShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.StampedeScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;
import org.opensha.sha.cybershake.etas.ETASModProbConfig.ETAS_CyberShake_Scenarios;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_Params.U3ETAS_ProbabilityModelOptions;
import scratch.UCERF3.utils.FaultSystemIO;
import scratch.UCERF3.erf.ETAS.ETAS_Simulator.TestScenario;
import scratch.kevin.ucerf3.etas.ETAS_BinaryCatalogFilterByMag;
import scratch.kevin.ucerf3.etas.ETAS_BinaryCatalogFilterDependents;
import scratch.kevin.ucerf3.etas.MPJ_ETAS_Simulator;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public class ScriptGen {
	
private static final String args_continue_newline = " \\\n\t";
	
	public static final DateFormat df = new SimpleDateFormat("yyyy_MM_dd");

	public static void main(String[] args) throws IOException, DocumentException {
		File localDir = new File("/home/kevin/OpenSHA/UCERF3/cybershake_etas/sims");
		
		String cacheDirName = "cache_u3rups_u2mapped";
		FaultModels fm = FaultModels.FM3_1;
		String solTypeStr = "u2mapped";
		String fssName = "ucerf2_mapped_sol.zip";
		
		File origSolFile = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip");
		File mappedSolFile = new File(localDir, fssName);
		
		boolean stampede = false;
//		int threads = 1;
//		String queue = null;
//		String pbsNameAdd = null;
		int threads = 8;
		String queue = "scec";
		String pbsNameAdd = "-scec";
		boolean smallTest = false;
		
		boolean writeConsolidate = true;
		boolean bundleConsolidate = true;
		
//		double duration = 10000;
//		int numSims = 100;
//		int hours = 24;
//		int nodes = 50;
		
//		double duration = 1000;
//		int numSims = 5000;
//		int hours = 24;
//		int nodes = 60;
		
//		double duration = 30;
//		int numSims = 5000;
//		int hours = 24;
//		int nodes = 60;
		
//		double duration = 10;
//		int numSims = 25000;
//		int hours = 24;
//		int nodes = 60;
		
		// for scenarios
		double duration = 10; // SIM duration!!
//		int numSims = 25000;
//		int hours = 24;
//		int nodes = 60;
		int numSims = 100000;
		int hours = 24;
//		int nodes = 100;
		int nodes = 24;
		
		// TODO deal with differnt rup IDs for supra
//		TestScenario[] scenarios = { TestScenario.BOMBAY_BEACH_M6 };
		TestScenario[] scenarios = { TestScenario.PARKFIELD };
		boolean includeSpontaneous = false;
		String customCatalog = null;
		long customOT = Long.MIN_VALUE;
		
//		TestScenario[] scenarios = { null };
//		boolean includeSpontaneous = false;
//		String customCatalog = "2016_bombay_swarm.txt";
////		long customOT = 1474920000000l;
//		long customOT = 1474990200000l;
		
		fixRupIndexes(scenarios, origSolFile, mappedSolFile);
		
		boolean griddedOnly = false;
		
//		U3ETAS_ProbabilityModelOptions[] probModels = U3ETAS_ProbabilityModelOptions.values();
//		U3ETAS_ProbabilityModelOptions[] probModels = {U3ETAS_ProbabilityModelOptions.FULL_TD,
//				U3ETAS_ProbabilityModelOptions.NO_ERT};
		U3ETAS_ProbabilityModelOptions[] probModels = {U3ETAS_ProbabilityModelOptions.FULL_TD};
//		double totRateScaleFactor = 1.14;
		double totRateScaleFactor = 1.0;
//		U3ETAS_ProbabilityModelOptions[] probModels = {U3ETAS_ProbabilityModelOptions.NO_ERT};
//		double totRateScaleFactor = 1.0;
//		U3ETAS_ProbabilityModelOptions[] probModels = {U3ETAS_ProbabilityModelOptions.POISSON};
//		boolean[] grCorrs = { false, true };
		boolean[] grCorrs = { false };
//		double[] maxCharFactors = { U3ETAS_MaxCharFactorParam.DEFAULT_VALUE };
//		double[] maxCharFactors = { 10 };
//		boolean applyLongTermRates = false;
		boolean gridSeisCorr = true; // TODO ?
		boolean applySubSeisForSupraNucl = false; // TODO ?
		
		String nameAdd = null;
//		String nameAdd = "100krun1";
//		String nameAdd = "20kmore4";
//		String nameAdd = "scaleMFD1p14";
//		String nameAdd = "newNuclWt";
//		String nameAdd = "4000more";
//		String nameAdd = "mc10-applyGrGridded";
//		String nameAdd = "FelzerParams-mc20";
		boolean nameAddAtEnd = true;
		
		boolean histCatalog = false; // keep false
		if (!includeSpontaneous)
			histCatalog = false;
		int startYear = 2012;
		int mins = hours*60;
		
		if (smallTest && stampede) {
			queue = "development";
			numSims = 200;
			nodes = 10;
			mins = 2*60;
			if (nameAdd == null)
				nameAdd = "";
			else
				nameAdd += "-";
			nameAdd += "quick_test";
		}
		
		String dateStr = df.format(new Date());
//		String dateStr = "2016_11_02";
//		String dateStr = "2016_02_24";
		
		boolean timeIndep = false;
		
		boolean binary = numSims >= 1000 || duration > 200;
		
		int memGigs;
		int ppn;
		if (stampede)
			ppn = 16;
		else
			ppn = 8;
		
		File remoteDir, remoteSolFile, cacheDir;
		JavaShellScriptWriter mpjWrite;
		BatchScriptWriter pbsWrite;
		
		if (stampede) {
			memGigs = 26;
			remoteDir = new File("/work/00950/kevinm/ucerf3/etas_sim/cybershake");
			remoteSolFile = new File(remoteDir, fssName);
			mpjWrite = new FastMPJShellScriptWriter(StampedeScriptWriter.JAVA_BIN, memGigs*1024,
					null, StampedeScriptWriter.FMPJ_HOME);
			((FastMPJShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
			pbsWrite = new StampedeScriptWriter();
			cacheDir = new File(remoteDir, cacheDirName);
		} else {
			if (queue == null)
				memGigs = 9;
			else
				memGigs = 60;
			remoteDir = new File("/home/scec-02/kmilner/ucerf3/etas_sim/cybershake");
			remoteSolFile = new File(remoteDir, fssName);
			boolean fmpj = nodes < 25;
			fmpj = false;
			if (fmpj) {
				mpjWrite = new FastMPJShellScriptWriter(USC_HPCC_ScriptWriter.JAVA_BIN, memGigs*1024,
						null, USC_HPCC_ScriptWriter.FMPJ_HOME);
				((FastMPJShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
			} else {
				mpjWrite = new MPJExpressShellScriptWriter(USC_HPCC_ScriptWriter.JAVA_BIN, memGigs*1024,
						null, USC_HPCC_ScriptWriter.MPJ_HOME);
			}
			pbsWrite = new USC_HPCC_ScriptWriter();
			cacheDir = new File(remoteDir, cacheDirName);
		}
		
		File histCatalogFile = null;
		File rupSurfacesFile = null;
//		if (histCatalog) {
//			histCatalogFile = new File(remoteDir, "ofr2013-1165_EarthquakeCat.txt");
//		}
//		if (histCatalog || customCatalog != null)
//			rupSurfacesFile = new File(remoteDir, "finite_fault_mappings.xml");
		
		mpjWrite.setAutoMemDetect(false);
		
		List<File> classpath = new ArrayList<File>();
		classpath.add(new File(remoteDir, "commons-cli-1.2.jar"));
		
		boolean exactDispatch = numSims / nodes == threads;
		
		for (TestScenario scenario : scenarios) {
			String scenarioName;
			if (scenario == null) {
				if (customCatalog != null && !customCatalog.isEmpty()) {
					scenarioName = customCatalog;
					if (scenarioName.contains("."))
						scenarioName = scenarioName.substring(0, scenarioName.indexOf("."));
				} else {
					scenarioName = "spontaneous";
				}
			} else {
				scenarioName = scenario.name().toLowerCase();
			}
//			if (duration > 1d) {
				if (duration == Math.floor(duration))
					scenarioName += "-"+(int)duration+"yr";
				else
					scenarioName += "-"+(float)duration+"yr";
//			}
			for (U3ETAS_ProbabilityModelOptions probModel : probModels) {
//				for (double maxCharFactor : maxCharFactors) {
				for (boolean grCorr : grCorrs) {
					String grStr;
					if (grCorr)
						grStr = "-grCorr";
					else
						grStr = "";
					String jobName = dateStr+"-"+scenarioName+"-"+solTypeStr;
					if (nameAdd != null && !nameAdd.isEmpty() && !nameAddAtEnd)
						jobName += "-"+nameAdd;
					if (griddedOnly) {
						jobName += "-gridded-only";
					} else {
//						jobName += "-"+probModel.name().toLowerCase()+"-maxChar"+(float)maxCharFactor;
						jobName += "-"+probModel.name().toLowerCase()+grStr;
						if (timeIndep)
							jobName += "-indep";
//						if (applyLongTermRates)
//							jobName += "-applyLTR";
						if (applySubSeisForSupraNucl)
							jobName += "-subSeisSupraNucl";
						if (gridSeisCorr)
							jobName += "-gridSeisCorr";
						if (totRateScaleFactor != 1)
							jobName += "-scale"+(float)totRateScaleFactor;
						if (!includeSpontaneous)
							jobName += "-noSpont";
					}
					
					if (nameAdd != null && !nameAdd.isEmpty() && nameAddAtEnd)
						jobName += "-"+nameAdd;
					
					File localJobDir = new File(localDir, jobName);
					if (!localJobDir.exists())
						localJobDir.mkdir();
					File remoteJobDir = new File(remoteDir, jobName);
					
					System.out.println(jobName);
					
					List<File> subClasspath = Lists.newArrayList(classpath);
					subClasspath.add(new File(remoteJobDir, "OpenSHA_complete.jar"));
					mpjWrite.setClasspath(subClasspath);
					
					String pbsName = jobName;
					if (pbsNameAdd != null)
						pbsName += pbsNameAdd;
					pbsName += ".pbs";
					File pbsFile = new File(localJobDir, pbsName);
					
					String argz;
					String sep;
					if (mpjWrite instanceof MPJExpressShellScriptWriter)
						sep = " ";
					else
						sep = args_continue_newline;
					
					if (exactDispatch) {
						argz = sep+"--min-dispatch "+threads
								+" --max-dispatch "+threads+" --exact-dispatch "+threads;
					} else {
						argz = sep+"--min-dispatch 1 --max-dispatch "+threads*10;
					}
					
					argz += sep+"--threads "+threads
							+sep+"--num "+numSims
							+sep+"--sol-file "+remoteSolFile.getAbsolutePath();
					
					argz += sep+"--duration "+(float)duration;
					if (customOT > Long.MIN_VALUE)
						argz += sep+"--millis "+customOT;
					else
						argz += sep+"--start-year "+startYear;
					
					argz += sep+"--prob-model "+probModel.name();
					
//					argz += " --max-char-factor "+maxCharFactor;
					if (grCorr)
						argz += sep+"--impose-gr";
					
//					argz += sep+"--apply-long-term-rates "+applyLongTermRates;
					argz += sep+"--apply-sub-seis-for-supra-nucl "+applySubSeisForSupraNucl;
					
					if (gridSeisCorr)
						argz += sep+"--grid-seis-correction";
					
					argz += sep+"--tot-rate-scale-factor "+totRateScaleFactor;
					
					if (timeIndep)
						argz += sep+"--indep";
					
					if (binary)
						argz += sep+"--binary";
					
					if (scenario != null) {
						if (scenario.getFSS_Index() >= 0)
							argz += sep+"--trigger-rupture-id "+scenario.getFSS_Index();
						Location loc = scenario.getLocation();
						if (loc != null)
							argz += sep+"--trigger-loc "+(float)loc.getLatitude()
								+","+(float)loc.getLongitude()+","+(float)loc.getDepth();
						if (scenario.getMagnitude() > 0)
							argz += sep+"--trigger-mag "+(float)scenario.getMagnitude();
					}
					if (customCatalog != null && !customCatalog.isEmpty()) {
						File myHistFile = new File(remoteJobDir, customCatalog);
						argz += sep+"--trigger-catalog "+myHistFile.getAbsolutePath();
						if (rupSurfacesFile != null)
							argz += sep+"--rupture-surfaces "+rupSurfacesFile.getAbsolutePath();
					} else {
						if (histCatalogFile != null)
							argz += sep+"--trigger-catalog "+histCatalogFile.getAbsolutePath();
						if (rupSurfacesFile != null)
							argz += sep+"--rupture-surfaces "+rupSurfacesFile.getAbsolutePath();
					}
					if (!includeSpontaneous)
						argz += sep+"--no-spontaneous";
					if (griddedOnly)
						argz += sep+"--gridded-only";
					
					argz += sep+cacheDir.getAbsolutePath()+sep+remoteJobDir.getAbsolutePath();
					
					List<String> script = mpjWrite.buildScript(MPJ_ETAS_Simulator.class.getName(), argz);
					
					List<String> consolidationLines = null;
					if (writeConsolidate) {
						consolidationLines = Lists.newArrayList();
//						consolidationLines.add("# M4 consolidation");
						File resultsDir = new File(remoteJobDir, "results");
						File m4File = new File(remoteJobDir, "results_m4_preserve.bin");
//						consolidationLines.add(mpjWrite.buildCommand(ETAS_CatalogIO.class.getName(),
//								resultsDir.getAbsolutePath()+" "+m4File.getAbsolutePath()+" 4"));
//						consolidationLines.add("");
						File resultsFile = new File(remoteJobDir, "results.bin");
						if (scenario == null) {
							if (!binary) {
								// build binary results.bin file
								consolidationLines.add("# create results.bin binary file");
								consolidationLines.add(mpjWrite.buildCommand(ETAS_CatalogIO.class.getName(),
										resultsDir.getAbsolutePath()+" "+resultsFile.getAbsolutePath()));
							}
						} else {
							// descendents file
							consolidationLines.add("# create descendents binary file");
							File descendentsFile = new File(remoteJobDir, "results_descendents.bin");
							consolidationLines.add(mpjWrite.buildCommand(ETAS_BinaryCatalogFilterDependents.class.getName(),
									resultsFile.getAbsolutePath()+" "+descendentsFile.getAbsolutePath()+" 0"));
						}
						// build m4 file, preserving descendents
						consolidationLines.add("# create results_m4_preserve.bin binary file");
						consolidationLines.add(mpjWrite.buildCommand(ETAS_BinaryCatalogFilterByMag.class.getName(),
								resultsFile.getAbsolutePath()+" "+m4File.getAbsolutePath()+" 4 true"));
						
						if (bundleConsolidate) {
							String lastLine = script.get(script.size()-1);
							if (lastLine.startsWith("exit")) {
								script.remove(script.size()-1);
								script.addAll(consolidationLines);
								script.add("");
								script.add(lastLine);
							} else {
								script.add("");
								script.addAll(consolidationLines);
							}
						}
					}
					
					script = pbsWrite.buildScript(script, mins, nodes, ppn, queue);
					pbsWrite.writeScript(pbsFile, script);
					
					if (writeConsolidate && !bundleConsolidate && stampede) {
						// write consolidation script as well (separately)
						script = Lists.newArrayList();
						
						script.add("#!/bin/bash");
						script.add("");
						script.add("JVM_MEM_MB="+memGigs*1024);
						script.add("");
						script.addAll(consolidationLines);
						
						pbsWrite.writeScript(new File(localJobDir, "consolidate_dev.pbs"), script, 60, 1, 16, "development");
						pbsWrite.writeScript(new File(localJobDir, "consolidate_norm.pbs"), script, 60, 1, 16, "normal");
					}
				}
			}
		}
	}
	
	private static void fixRupIndexes(TestScenario[] scenarios, File origSolFile, File mappedSolFile)
			throws IOException, DocumentException {
		FaultSystemRupSet origRupSet = null;
		FaultSystemRupSet mappedRupSet = null;
		
		for (TestScenario scenario : scenarios) {
			if (scenario == null || scenario.getFSS_Index() < 0)
				continue;
			
			if (origRupSet == null) {
				origRupSet = FaultSystemIO.loadRupSet(origSolFile);
				mappedRupSet = FaultSystemIO.loadRupSet(mappedSolFile);
				
				if (origRupSet.getNumRuptures() == mappedRupSet.getNumRuptures())
					System.out.println("Rupture counts identical, mapping likely not needed");
			}
			
			int index = scenario.getFSS_Index();
			HashSet<Integer> sects = new HashSet<Integer>(origRupSet.getSectionsIndicesForRup(index));
			
			int matchIndex = -1;
			for (int rupIndex : mappedRupSet.getRupturesForSection(sects.iterator().next())) {
				List<Integer> mySects = mappedRupSet.getSectionsIndicesForRup(rupIndex);
				if (mySects.size() != sects.size())
					continue;
				boolean match = true;
				for (int sect : mySects) {
					if (!sects.contains(sect)) {
						match = false;
						break;
					}
				}
				if (match) {
					matchIndex = rupIndex;
					break;
				}
			}
			Preconditions.checkState(matchIndex >= 0, "No match for scenario "+scenario.name()+" in mapped solution!");
			System.out.println("Mapped rupture "+index+" to: "+matchIndex);
			scenario.setFSS_Index(matchIndex);
		}
	}

	public static void main2(String[] args) throws IOException {
		File localDir = new File("/home/kevin/OpenSHA/UCERF3/cybershake_etas/sims");
		
		String cacheDirName = "cache_u3rups_u2mapped";
		FaultModels fm = FaultModels.FM3_1;
		String solTypeStr = "u2mapped";
		String fssName = "ucerf2_mapped_sol.zip";
//		String cacheDirName = "cache_u2rups_u3inverted";
//		FaultModels fm = FaultModels.FM2_1;
//		String solTypeStr = "u3inverted";
//		String fssName = "ucerf2_u3inverted_sol.zip";
		
//		String dateStr = new SimpleDateFormat("yyyy_MM_dd").format(new Date());
//		String dateStr = "2015_03_23";
//		String dateStr = "2015_04_09";
		String dateStr = "2015_06_15";
		
//		Scenarios scenario = Scenarios.LA_HABRA;
//		Scenarios[] scenarios = Scenarios.values();
//		Scenarios[] scenarios = {Scenarios.BOMBAY_BEACH};
//		Scenarios[] scenarios = {Scenarios.BOMBAY_BEACH_M6};
//		Scenarios[] scenarios = {Scenarios.PARKFIELD};
		ETAS_CyberShake_Scenarios[] scenarios = {
//				ETAS_CyberShake_Scenarios.BOMBAY_BEACH_BRAWLEY_FAULT_M6,
//				ETAS_CyberShake_Scenarios.PARKFIELD,
//				ETAS_CyberShake_Scenarios.MOJAVE_S_POINT_M6};
//				ETAS_CyberShake_Scenarios.PARKFIELD};
//				ETAS_CyberShake_Scenarios.BOMBAY_BEACH_M6,
//				ETAS_CyberShake_Scenarios.PARKFIELD,
				ETAS_CyberShake_Scenarios.MOJAVE_S_POINT_M6};
		boolean timeIndep = false;
		int numSims = 50000;
		String nameAdd = "-round2";
//		String nameAdd = "-nospont-round6";
//		String nameAdd = "-round5";
		
		File remoteDir, remoteSolFile;
		FastMPJShellScriptWriter mpjWrite;
		BatchScriptWriter pbsWrite;
		
//		int memGigs = 10;
//		int perNodeMemGB = 32;
//		int mins = 24*60;
//		int nodes = 99;
//		int ppn = 8;
//		String queue = null;
//		int threads = 1;
//		remoteDir = new File("/home/scec-02/kmilner/ucerf3/etas_sim/cybershake");
//		remoteSolFile = new File(remoteDir, fssName);
//		mpjWrite = new FastMPJShellScriptWriter(USC_HPCC_ScriptWriter.JAVA_BIN, memGigs*1024,
//				null, USC_HPCC_ScriptWriter.FMPJ_HOME, false);
//		pbsWrite = new USC_HPCC_ScriptWriter();
//		((USC_HPCC_ScriptWriter)pbsWrite).setPerNodeMemGB(perNodeMemGB);
		
		int memGigs = 25;
		int mins = 12*60;
		int nodes = 128;
		int ppn = 16;
		String queue = null;
		int threads = 8;
		remoteDir = new File("/work/00950/kevinm/ucerf3/etas_sim/cybershake");
		remoteSolFile = new File(remoteDir, fssName);
		mpjWrite = new FastMPJShellScriptWriter(StampedeScriptWriter.JAVA_BIN, memGigs*1024,
				null, StampedeScriptWriter.FMPJ_HOME);
		pbsWrite = new StampedeScriptWriter();
		
		List<File> classpath = new ArrayList<File>();
		classpath.add(new File(remoteDir.getParentFile(), "commons-cli-1.2.jar"));
		
		for (ETAS_CyberShake_Scenarios scenario : scenarios) {
			String jobName = dateStr+"-"+solTypeStr+"-"+scenario.name().toLowerCase();
			if (timeIndep)
				jobName += "-indep";
			if (nameAdd != null)
				jobName += nameAdd;
			
			File localJobDir = new File(localDir, jobName);
			if (!localJobDir.exists())
				localJobDir.mkdir();
			File remoteJobDir = new File(remoteDir, jobName);
			
			List<File> subClasspath = Lists.newArrayList(classpath);
			subClasspath.add(new File(remoteJobDir, "OpenSHA_complete.jar"));
			mpjWrite.setClasspath(subClasspath);
			
			File pbsFile = new File(localJobDir, jobName+".pbs");
			
			String argz = "--min-dispatch 1 --max-dispatch "+threads+" --threads "+threads+" --no-spontaneous"
					+" --num "+numSims+" --sol-file "+remoteSolFile.getAbsolutePath();
//			case BOMBAY_BEACH_CAT:
//				argz += " --trigger-catalog "+(new File(remoteDir, "bombay_catalog.txt")).getAbsolutePath();
//				break;
//			case BOMBAY_BEACH_SINGLE:
//				argz += " --trigger-loc 33.31833333333334,-115.72833333333335,5.8 --trigger-mag 4.8";
//				break;
			Preconditions.checkState(scenario.getTriggerRupIndex(fm) >= 0
					|| (scenario.getTriggerLoc() != null && scenario.getTriggerMag() > 0d));
			if (scenario.getTriggerRupIndex(fm) >= 0) {
				argz += " --trigger-rupture-id "+scenario.getTriggerRupIndex(fm);
			} else {
				Location loc = scenario.getTriggerLoc();
				argz += " --trigger-loc "+loc.getLatitude()+","+loc.getLongitude()+","+loc.getDepth();
			}
			if (scenario.getTriggerMag() > 0d)
				argz += " --trigger-mag "+scenario.getTriggerMag();
			if (timeIndep)
				argz += " --indep";
			argz += " "+new File(remoteDir, cacheDirName).getAbsolutePath()+" "+remoteJobDir.getAbsolutePath();
			
			List<String> script = mpjWrite.buildScript(MPJ_ETAS_Simulator.class.getName(), argz);
			
			script = pbsWrite.buildScript(script, mins, nodes, ppn, queue);
			pbsWrite.writeScript(pbsFile, script);
		}
	}

}
