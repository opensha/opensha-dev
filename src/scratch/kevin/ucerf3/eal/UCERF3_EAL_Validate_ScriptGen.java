package scratch.kevin.ucerf3.eal;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Random;

import org.dom4j.Document;
import org.dom4j.Element;
import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.StampedeScriptWriter;
import org.opensha.commons.util.XMLUtils;
import org.opensha.sha.earthquake.param.BackgroundRupParam;
import org.opensha.sha.earthquake.param.BackgroundRupType;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sra.calc.parallel.MPJ_CondLossCalc;
import org.opensha.sra.calc.parallel.MPJ_EAL_Calc;

import com.google.common.collect.Lists;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.simulatedAnnealing.hpc.LogicTreePBSWriter;
import scratch.UCERF3.utils.FaultSystemIO;

public class UCERF3_EAL_Validate_ScriptGen {

	public static void main(String[] args) throws IOException {
		File writeDir = new File("/home/kevin/OpenSHA/UCERF3/eal");
		if (!writeDir.exists())
			writeDir.mkdir();
		
//		String runSubDirName = "2013_11_04-eal-calc-small-test";
		String runSubDirName = "2013_11_05-ucerf3-eal-calc-CB-2013-validate";
		
		writeDir = new File(writeDir, runSubDirName);
		if (!writeDir.exists())
			writeDir.mkdir();
		
//		BatchScriptWriter pbsWrite = new USC_HPCC_ScriptWriter();
//		File remoteDir = new File("/auto/scec-02/kmilner/ucerf3/curves/MeanUCERF3-curves");
//		File javaBin = USC_HPCC_ScriptWriter.JAVA_BIN;
//		File mpjHome = USC_HPCC_ScriptWriter.FMPJ_HOME;
//		int maxHeapMB = 9000;
		
		BatchScriptWriter pbsWrite = new StampedeScriptWriter();
		File remoteMainDir = new File("/work/00950/kevinm/ucerf3/eal");
		File remoteSubDir = new File(remoteMainDir, runSubDirName);
		File javaBin = StampedeScriptWriter.JAVA_BIN;
		File mpjHome = StampedeScriptWriter.FMPJ_HOME;
		int maxHeapMB = 26000;
		
		File compoundFile = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip");
		CompoundFaultSystemSolution cfss = CompoundFaultSystemSolution.fromZipFile(compoundFile);
		
		int numBranches = 5;
		
//		String meanSolFileName = "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_TRUE_HAZARD_MEAN_SOL_WITH_MAPPING.zip";
//		File meanSolFile = new File(remoteSubDir, meanSolFileName);
		
		String vulnFileName = "2012_01_02_VUL06.txt";
		File vulnFile = new File(remoteSubDir, vulnFileName);
		
		String portfolioFileName = "Porter (30 Oct 2013) CEA proxy portfolio.csv";
//		String portfolioFileName = "small_test_port.csv";
		File portfolioFile = new File(remoteSubDir, portfolioFileName);
		
		FastMPJShellScriptWriter javaWrite = new FastMPJShellScriptWriter(javaBin, maxHeapMB,
				LogicTreePBSWriter.getClasspath(remoteSubDir, remoteSubDir), mpjHome);
		
//		JavaShellScriptWriter javaWrite = new JavaShellScriptWriter(javaBin, maxHeapMB,
//				LogicTreePBSWriter.getClasspath(remoteDir, remoteDir));
		
		int mins = 300;
		int nodes = 20;
		int ppn = 8;
		String queue = null;
		
//		String className = MPJ_CondLossCalc.class.getName();
		String className = MPJ_EAL_Calc.class.getName();
		
//		AttenRelRef[] imrs = { AttenRelRef.CB_2014, AttenRelRef.CY_2014,
//				AttenRelRef.ASK_2014, AttenRelRef.BSSA_2014, AttenRelRef.IDRISS_2014 };
		AttenRelRef[] imrs = { AttenRelRef.CB_2014 };
		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF();
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.INCLUDE);
		erf.setParameter(BackgroundRupParam.NAME, BackgroundRupType.CROSSHAIR);
		
		Random r = new Random();
		List<LogicTreeBranch> branches = Lists.newArrayList(cfss.getBranches());
		
		IncludeBackgroundOption[] bgIncls = { IncludeBackgroundOption.EXCLUDE, IncludeBackgroundOption.ONLY };
		
		for (AttenRelRef ref : imrs) {
			for (int i=0; i<numBranches; i++) {
				LogicTreeBranch branch = branches.get(r.nextInt(branches.size()));
				
				FaultSystemSolution sol = cfss.getSolution(branch);
				
				File localSolFile = new File(writeDir, branch.buildFileName()+"_sol.zip");
				File remoteSolFile = new File(remoteSubDir, branch.buildFileName()+"_sol.zip");
				FaultSystemIO.writeSol(sol, localSolFile);
				
				erf.setParameter(FaultSystemSolutionERF.FILE_PARAM_NAME, remoteSolFile);
				
				for (IncludeBackgroundOption bgIncl : bgIncls) {
					erf.setParameter(IncludeBackgroundParam.NAME, bgIncl);
					
					String name = branch.buildFileName()+"_"+ref.name()+"_bg"+bgIncl.name();
					File localXML = new File(writeDir, name+".xml");
					File remoteXML = new File(remoteSubDir, name+".xml");
					
					Document doc = XMLUtils.createDocumentWithRoot();
					Element root = doc.getRootElement();
					erf.toXMLMetadata(root);
					ScalarIMR imr = ref.instance(null);
					imr.toXMLMetadata(root);
					
					XMLUtils.writeDocumentToFile(localXML, doc);
					
					File remoteOutput = new File(remoteSubDir, name+".txt");
					
					String jobArgs = "--dist-func --vuln-file \""+vulnFile.getAbsolutePath()
							+"\" \""+portfolioFile.getAbsolutePath()+"\" "
							+remoteXML.getAbsolutePath()+" "+remoteOutput.getAbsolutePath();
					
					File jobFile = new File(writeDir, name+".pbs");
					pbsWrite.writeScript(jobFile, javaWrite.buildScript(className, jobArgs), mins, nodes, ppn, queue);
				}
			}
		}
	}

}
