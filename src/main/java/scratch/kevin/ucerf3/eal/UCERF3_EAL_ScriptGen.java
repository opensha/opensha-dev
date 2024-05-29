package scratch.kevin.ucerf3.eal;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import org.dom4j.Document;
import org.dom4j.Element;
import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.StampedeScriptWriter;
import org.opensha.commons.hpc.pbs.USC_CARC_ScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;
import org.opensha.commons.util.XMLUtils;
import org.opensha.sha.earthquake.param.BackgroundRupParam;
import org.opensha.sha.earthquake.param.BackgroundRupType;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_WrapperFullParam;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_WrapperFullParam.EpistemicOption;
import org.opensha.sra.calc.parallel.MPJ_CondLossCalc;

import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.simulatedAnnealing.hpc.LogicTreePBSWriter;

public class UCERF3_EAL_ScriptGen {

	public static void main(String[] args) throws IOException {
		File writeDir = new File("/home/kevin/OpenSHA/UCERF3/eal");
		if (!writeDir.exists())
			writeDir.mkdir();
		
//		String runSubDirName = "2013_11_04-eal-calc-small-test";
//		String runSubDirName = "2014_01_15-ucerf3-eal-calc-NGA2s-2013";
//		String runSubDirName = "2014_03_19-ucerf3-eal-calc-CB2014-recalc";
//		String runSubDirName = "2014_05_05-ucerf3-eal-calc-wald-vs30";
//		String runSubDirName = "2014_05_16-ucerf3-99percent-wills";
//		String runSubDirName = "2014_05_19-ucerf3-fatality";
//		String runSubDirName = "2014_05_28-ucerf3-99percent-wills-smaller";
//		String runSubDirName = "2014_05_28-ucerf3-fatality-smaller";
//		String runSubDirName = "2016_10_06-ucerf3-90percent-wald-san-bernardino";
//		String runSubDirName = "2016_10_21-ucerf3-90percent-wald-coachella-valley";
//		String runSubDirName = "2016_11_28-ucerf3-90percent-wills-san-bernardino";
//		String runSubDirName = "2016_11_28-ucerf3-90percent-wills-coachella-valley";
//		String runSubDirName = "2017_05_22-stampede-2-test";
//		String runSubDirName = "2017_05_24-ucerf3-ngaw2-cea-proxy-wills2015";
//		String runSubDirName = "2017_05_26-ucerf3-ngaw2-cea-proxy-wald";
//		String runSubDirName = "2020_03_17-ucerf3-ngaw2-cea-proxy-100pct-wills2015";
//		String runSubDirName = "2020_03_17-ucerf3-ngaw2-cea-proxy-100pct-wald2007";
//		String runSubDirName = "2024_05_24-ucerf3-fm31_ba-ngaw2-cea-proxy-100pct-wills2015";
		String runSubDirName = "2024_05_24-ucerf3-fm31_ba-pt_gridded-ngaw2-cea-proxy-100pct-wills2015";
		
//		EpistemicOption ngaEpistemic = EpistemicOption.UPPER;
//		EpistemicOption ngaEpistemic = EpistemicOption.LOWER;
		EpistemicOption ngaEpistemic = null;
		
		writeDir = new File(writeDir, runSubDirName);
		if (!writeDir.exists())
			writeDir.mkdir();
		
		BatchScriptWriter pbsWrite = new USC_CARC_ScriptWriter();
		File remoteMainDir = new File("/project/scec_608/kmilner/ucerf3/eal");
		File remoteSubDir = new File(remoteMainDir, runSubDirName);
		File javaBin = USC_CARC_ScriptWriter.JAVA_BIN;
		File mpjHome = USC_CARC_ScriptWriter.FMPJ_HOME;
		int maxHeapMB = 40000;
		int numThreads = 20;
		int maxDispatch = numThreads*5;
		int nodes = 36;
		String queue = "scec";
		int mins = 12*60;
		int ppn = 20;
		
//		BatchScriptWriter pbsWrite = new StampedeScriptWriter(true);
////		File remoteMainDir = new File("/work/00950/kevinm/ucerf3/eal");
//		File remoteMainDir = new File("/scratch/00950/kevinm/ucerf3/eal");
//		File remoteSubDir = new File(remoteMainDir, runSubDirName);
//		File javaBin = StampedeScriptWriter.JAVA_BIN;
//		File mpjHome = StampedeScriptWriter.FMPJ_HOME;
////		int maxHeapMB = 26000;
////		int numThreads = -1;
//		int maxHeapMB = 72*1024;
//		int numThreads = 68*2;
//		int maxDispatch = numThreads*5;
//		int nodes = 20;
//		String queue = "normal";
//		int mins = 12*60;
		
		boolean gzip = false;
		boolean tractResults = false;
		
//		String meanSolFileName = "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_TRUE_HAZARD_MEAN_SOL_WITH_MAPPING.zip";
//		File meanSolFile = new File(remoteMainDir, meanSolFileName);
		File meanSolFile = new File("/project/scec_608/kmilner/ucerf3/ucerf3-etas-launcher/inputs/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_SpatSeisU3_MEAN_BRANCH_AVG_SOL.zip");
		
//		String vulnFileName = "2012_01_02_VUL06.txt";
		String vulnFileName = "2014_05_16_VUL06.txt"; // updated for 99%
//		String vulnFileName = "2014_05_16b_VUL06.txt"; // fatalities
		File vulnFile = new File(remoteMainDir, vulnFileName);
		
//		String portfolioFileName = "Porter (30 Oct 2013) CEA proxy portfolio.csv"; // 2013 values, Wills Vs30
//		String portfolioFileName = "Porter (05 May 2014) CEA proxy portfolio 2013 values Wald-Allen Vs30.csv";
//		String portfolioFileName = "small_test_port.csv";
//		String portfolioFileName = "Porter (10 May 2014) CA 99pct portfolio 2013 values Wills Vs30.txt";
//		String portfolioFileName = "Porter (16 May 2014) SCEC UCERF3 CA fatality portfolio.txt"; // fatality portfolio, Wills Vs30
//		String portfolioFileName = "Porter-22-May-14-CA-CAS4-90pct-Wills.txt"; // smaller fatality portfolio, Wills Vs30
//		String portfolioFileName = "Porter-22-May-14-CA-ppty-90pct-Wills.txt"; // smaller 90%
//		String portfolioFileName = "Porter-02-Jun-16-CA-ppty-90pct-Wald.txt"; // smaller 90%, Wald version
//		String portfolioFileName = "san_bernardino_Porter-02-Jun-16-CA-ppty-90pct-Wald.txt"; // 90% Wald, san bernardino city
//		String portfolioFileName = "coachella_valley_Porter-02-Jun-16-CA-ppty-90pct-Wald.txt"; // 90% Wald, Coachella Valley (20km circle, Rancho Mirage)
//		String portfolioFileName = "san_bernardino_Porter-22-May-14-CA-ppty-90pct-Wills.txt"; // 90% Wills, san bernardino city
//		String portfolioFileName = "coachella_valley_Porter-22-May-14-CA-ppty-90pct-Wills.txt"; // 90% Wills, Coachella Valley (20km circle, Rancho Mirage)
		// 2017 CEA proxy
//		String portfolioFileName = "Porter-24May2017-CA-RES1-2017-Wills2015.csv"; // Wills 2015
//		String portfolioFileName = "Porter-24May2017-CA-RES1-2017-Wald.csv"; // Wald
		String portfolioFileName = "Porter-09-Feb-2020-CEA-100-pct-procy-portfolio-wills2015.csv"; // From Keith 
//		String portfolioFileName = "Porter-09-Feb-2020-CEA-100-pct-procy-portfolio-wald2007.csv"; // From Keith 
		File portfolioFile = new File(remoteMainDir, portfolioFileName);
		
		FastMPJShellScriptWriter javaWrite = new FastMPJShellScriptWriter(javaBin, maxHeapMB,
				LogicTreePBSWriter.getClasspath(remoteMainDir, remoteSubDir), mpjHome);
		javaWrite.setUseLaunchWrapper(true);
		javaWrite.setInitialHeapSizeMB(maxHeapMB);
		
//		JavaShellScriptWriter javaWrite = new JavaShellScriptWriter(javaBin, maxHeapMB,
//				LogicTreePBSWriter.getClasspath(remoteDir, remoteDir));
		
////		int nodes = 80;
//		int nodes = 20;
//		int ppn = 8;
//		if (numThreads > 0)
//			ppn = numThreads;
//		if (numThreads > 68)
//			ppn = 68;
//		String queue = "normal";
		
		String className = MPJ_CondLossCalc.class.getName();
		
//		AttenRelRef[] imrs = { AttenRelRef.CB_2014, AttenRelRef.CY_2014,
//				AttenRelRef.ASK_2014, AttenRelRef.BSSA_2014, AttenRelRef.IDRISS_2014 };
//		AttenRelRef[] imrs = { AttenRelRef.CB_2014 };
//		AttenRelRef[] imrs = { AttenRelRef.ASK_2014 };
		AttenRelRef[] imrs = { AttenRelRef.NGAWest_2014_AVG_NOIDRISS };
//		AttenRelRef[] imrs = { AttenRelRef.CY_2014,
//				AttenRelRef.ASK_2014, AttenRelRef.BSSA_2014, AttenRelRef.IDRISS_2014 };
		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF();
		erf.setParameter(FaultSystemSolutionERF.FILE_PARAM_NAME, meanSolFile);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.INCLUDE);
//		erf.setParameter(BackgroundRupParam.NAME, BackgroundRupType.CROSSHAIR);
		erf.setParameter(BackgroundRupParam.NAME, BackgroundRupType.POINT);
		
		for (AttenRelRef ref : imrs) {
			String name = ref.name();
			if (ngaEpistemic != null)
				name += "_"+ngaEpistemic.name();
			File localXML = new File(writeDir, name+".xml");
			File remoteXML = new File(remoteSubDir, name+".xml");
			
			Document doc = XMLUtils.createDocumentWithRoot();
			Element root = doc.getRootElement();
			erf.toXMLMetadata(root);
			ScalarIMR imr = ref.instance(null);
			imr.setParamDefaults();
			if (ngaEpistemic != null)
				imr.getParameter(NGAW2_WrapperFullParam.EPISTEMIC_PARAM_NAME).setValue(ngaEpistemic);
			imr.toXMLMetadata(root);
			
			XMLUtils.writeDocumentToFile(localXML, doc);
			
			File remoteOutput = new File(remoteSubDir, name+".bin");
			
			String jobArgs = "";
			if (maxDispatch > 0) {
				if (maxDispatch >= numThreads && numThreads > 0)
					jobArgs += "--min-dispatch "+numThreads+" ";
				jobArgs += "--max-dispatch "+maxDispatch+" ";
			}
			if (numThreads > 0)
				jobArgs += "--threads "+numThreads+" ";
			if (gzip)
				jobArgs += "--gzip ";
			if (tractResults)
				jobArgs += "--tract-results ";
			
			jobArgs += "--vuln-file \""+vulnFile.getAbsolutePath()+"\" \""+portfolioFile.getAbsolutePath()+"\" "
					+remoteXML.getAbsolutePath()+" "+remoteOutput.getAbsolutePath();
			
			File jobFile = new File(writeDir, name+".slurm");
			pbsWrite.writeScript(jobFile, javaWrite.buildScript(className, jobArgs), mins, nodes, ppn, queue);
		}
	}

}
