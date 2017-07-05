package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Map;

import org.opensha.commons.geo.Location;
import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter;
import org.opensha.commons.hpc.mpj.MPJExpressShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.StampedeScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;

import scratch.UCERF3.erf.ETAS.ETAS_Simulator.TestScenario;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_Params.U3ETAS_MaxCharFactorParam;
import scratch.UCERF3.erf.ETAS.ETAS_Params.U3ETAS_ProbabilityModelOptions;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class MPJ_ETAS_SimulatorScriptGen {
	
	private static final String args_continue_newline = " \\\n\t";
	
	public static final DateFormat df = new SimpleDateFormat("yyyy_MM_dd");

	public static void main(String[] args) throws IOException {
		File localDir = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations");
		
		boolean stampede = false;
//		int threads = 2;
//		String queue = null;
//		String pbsNameAdd = null;
		boolean largeSCEC = false;
		int threads;
		String pbsNameAdd;
		if (largeSCEC) {
			threads = 20;
			pbsNameAdd = "-scec-large";
		} else {
			threads = 8;
			pbsNameAdd = "-scec";
		}
		String queue = "scec";
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
		
		// CSEP benchmarks
//		double duration = 1d/365.25;
		double duration = 3d/12d;
		int numSims = 10000;
		int hours = 24;
		int nodes = 10;
		
		// for scenarios
//		double duration = 10; // SIM duration!!
////		int numSims = 25000;
////		int hours = 24;
////		int nodes = 60;
//		int numSims = 10000;
//		int hours = 24;
//////		int nodes = 100;
//		int nodes = 34;
//		if (largeSCEC && !stampede)
//			nodes = 4;
////		int numSims = 1000;
////		int hours = 2;
//////		int nodes = 100;
////		int nodes = 5;
		
//		Scenarios scenario = Scenarios.LA_HABRA;
//		Scenarios[] scenarios = Scenarios.values();
//		Scenarios[] scenarios = {Scenarios.MOJAVE_7};
//		Scenarios[] scenarios = {Scenarios.NAPA};
//		Scenarios[] scenarios = {Scenarios.SPONTANEOUS};
		
//		TestScenario.NORTHRIDGE.updateMag(6.7);
//		TestScenario[] scenarios = {TestScenario.NORTHRIDGE};
//		TestScenario[] scenarios = {TestScenario.MOJAVE_M7};
//		TestScenario[] scenarios = {TestScenario.PARKFIELD};
//		TestScenario[] scenarios = {TestScenario.HAYWIRED_M7};
//		TestScenario[] scenarios = {TestScenario.SAN_JACINTO_0_M4p8};
//		TestScenario[] scenarios = {TestScenario.MOJAVE_M5, TestScenario.MOJAVE_M5p5,
//				TestScenario.MOJAVE_M6pt3_ptSrc, TestScenario.MOJAVE_M6pt3_FSS, TestScenario.MOJAVE_M7};
//		TestScenario[] scenarios = {TestScenario.MOJAVE_M7pt4, TestScenario.MOJAVE_M7pt8};
//		TestScenario[] scenarios = {TestScenario.SAN_JACINTO_0_M5p5, TestScenario.SURPRISE_VALLEY_5p0,
//					TestScenario.SAF_PENINSULA_M5p5, TestScenario.SAF_PENINSULA_M6p3, TestScenario.SAF_PENINSULA_M7};
//		TestScenario[] scenarios = {TestScenario.SURPRISE_VALLEY_5p5, TestScenario.CENTRAL_VALLEY_M5p5};
//		boolean includeSpontaneous = true;
		
//		TestScenario[] scenarios = {TestScenario.SAN_JACINTO_0_M4p8};
//		boolean includeSpontaneous = false;
		
//		TestScenario[] scenarios = {TestScenario.BOMBAY_BEACH_M4pt8};
//		boolean includeSpontaneous = false;
		
//		TestScenario[] scenarios = {TestScenario.MOJAVE_M5p5, TestScenario.MOJAVE_M6pt3_ptSrc,
//				TestScenario.MOJAVE_M6pt3_FSS, TestScenario.MOJAVE_M7};
//		boolean includeSpontaneous = true;
//		TestScenario[] scenarios = {TestScenario.HAYWIRED_M7};
//		TestScenario[] scenarios = { null };
//		boolean includeSpontaneous = false;
//		String customCatalog = null;
//		long customOT = Long.MIN_VALUE;
//		String resetSectsArg = null;
//		boolean griddedOnly = false;
//		boolean customCatIncludeHistSurfaces = false;
		
//		TestScenario[] scenarios = { null };
////		boolean includeSpontaneous = false;
////		String customCatalog = "2016_bombay_swarm.txt";
//		boolean includeSpontaneous = true;
//		String customCatalog = "2016_bombay_swarm_with_hist.txt";
////		long customOT = 1474920000000l;
//		long customOT = 1474990200000l;
		
		// CSEP like benchmark tests
		TestScenario[] scenarios = { null };
		boolean includeSpontaneous = false;
		String customCatalog = "CSEP_1yr_Input_Catalog.txt";
		boolean customCatIncludeHistSurfaces = false;
		long customOT = 1498201200000l;
		boolean griddedOnly = false;
		String resetSectsArg = null;
		
		// USGS 5/17/17 exercise
//		TestScenario[] scenarios = { null };
////		boolean includeSpontaneous = false;
////		String customCatalog = "2016_bombay_swarm.txt";
//		boolean includeSpontaneous = false;
//		String customCatalog = "USGS_Exercise_Catalog.txt";
////		long customOT = 1495039984000l; // make sure this is at least a few seconds after the last event
//		long customOT = 1495051260000l;
//		String resetSectsArg = "1495039984000:2370,2369,2368,2367";
		
//		U3ETAS_ProbabilityModelOptions[] probModels = U3ETAS_ProbabilityModelOptions.values();
//		U3ETAS_ProbabilityModelOptions[] probModels = {U3ETAS_ProbabilityModelOptions.FULL_TD,
//				U3ETAS_ProbabilityModelOptions.NO_ERT};
		U3ETAS_ProbabilityModelOptions[] probModels = {U3ETAS_ProbabilityModelOptions.FULL_TD};
		double totRateScaleFactor = 1.14;
//		U3ETAS_ProbabilityModelOptions[] probModels = {U3ETAS_ProbabilityModelOptions.NO_ERT};
//		double totRateScaleFactor = 1.0;
//		U3ETAS_ProbabilityModelOptions[] probModels = {U3ETAS_ProbabilityModelOptions.POISSON};
//		boolean[] grCorrs = { false, true };
		boolean[] grCorrs = { false };
//		double[] maxCharFactors = { U3ETAS_MaxCharFactorParam.DEFAULT_VALUE };
//		double[] maxCharFactors = { 10 };
//		boolean applyLongTermRates = false;
		boolean gridSeisCorr = true;
		boolean applySubSeisForSupraNucl = true;
		
		if (griddedOnly)
			totRateScaleFactor = 1d;
		
		String nameAdd = null;
//		String nameAdd = "sect-reset-1pm";
//		String nameAdd = "small-speed-test";
//		String nameAdd = "100krun1";
//		String nameAdd = "20kmore4";
//		String nameAdd = "scaleMFD1p14";
//		String nameAdd = "newNuclWt";
//		String nameAdd = "4000more";
//		String nameAdd = "mc10-applyGrGridded";
//		String nameAdd = "FelzerParams-mc20";
		boolean nameAddAtEnd = true;
		
		boolean histCatalog = true;
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
			remoteDir = new File("/work/00950/kevinm/ucerf3/etas_sim");
			remoteSolFile = new File("/work/00950/kevinm/ucerf3/inversion/compound_plots/2013_05_10-ucerf3p3-production-10runs/"
					+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_SpatSeisU3_MEAN_BRANCH_AVG_SOL.zip");
			mpjWrite = new FastMPJShellScriptWriter(StampedeScriptWriter.JAVA_BIN, memGigs*1024,
					null, StampedeScriptWriter.FMPJ_HOME);
			((FastMPJShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
			pbsWrite = new StampedeScriptWriter();
			cacheDir = new File(remoteDir, "cache_fm3p1_ba");
		} else {
			if (queue == null)
				memGigs = 9;
			else if (largeSCEC)
				memGigs = 220;
			else
				memGigs = 60;
			remoteDir = new File("/home/scec-02/kmilner/ucerf3/etas_sim");
			remoteSolFile = new File("/home/scec-02/kmilner/ucerf3/inversion_compound_plots/"
					+ "2013_05_10-ucerf3p3-production-10runs/"
					+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_SpatSeisU3_MEAN_BRANCH_AVG_SOL.zip");
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
			cacheDir = new File(remoteDir, "cache_fm3p1_ba");
		}
		
		File histCatalogFile = null;
		File rupSurfacesFile = null;
		if (histCatalog) {
			histCatalogFile = new File(remoteDir, "ofr2013-1165_EarthquakeCat.txt");
		}
		if ((histCatalog && customCatalog == null) || (customCatalog != null && customCatIncludeHistSurfaces))
			rupSurfacesFile = new File(remoteDir, "finite_fault_mappings.xml");
		
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
				if (scenario.getFSS_Index() >= 0 && scenario.getMagnitude() > 0)
					scenarioName += "-m"+(float)scenario.getMagnitude();
			}
//			if (duration > 1d) {
			if (duration == Math.floor(duration)) {
				scenarioName += "-"+(int)duration+"yr";
			} else {
				if (duration < 1d) {
					// check for months, weeks, days
					double months = duration * 12d;
					double days = duration * 365.25;
					if ((float)months >= 1f) {
						if ((float)months == (float)Math.floor(months))
							scenarioName += "-"+(int)months+"mo";
						else
							scenarioName += "-"+(float)months+"mo";
					} else {
						if ((float)days == (float)Math.floor(days))
							scenarioName += "-"+(int)days+"day";
						else
							scenarioName += "-"+(float)days+"day";
					}
				} else {
					scenarioName += "-"+(float)duration+"yr";
				}
			}
//			}
			for (U3ETAS_ProbabilityModelOptions probModel : probModels) {
//				for (double maxCharFactor : maxCharFactors) {
				for (boolean grCorr : grCorrs) {
					String grStr;
					if (grCorr)
						grStr = "-grCorr";
					else
						grStr = "";
					String jobName = dateStr+"-"+scenarioName;
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
					}
					if (!includeSpontaneous)
						jobName += "-noSpont";
					
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
						argz = sep+"--min-dispatch 1 --max-dispatch "+threads*40;
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
					if (resetSectsArg != null && !resetSectsArg.isEmpty())
						argz += sep+"--reset-sections "+resetSectsArg;
					
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

}
