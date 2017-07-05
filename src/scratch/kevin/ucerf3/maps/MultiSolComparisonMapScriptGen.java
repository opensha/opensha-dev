package scratch.kevin.ucerf3.maps;

import java.io.File;
import java.io.IOException;
import java.util.Date;
import java.util.List;

import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.nshmp2.tmp.TestGrid;
import org.opensha.nshmp2.util.SourceIMR;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.io.Files;

import scratch.UCERF3.simulatedAnnealing.hpc.LogicTreePBSWriter;
import scratch.UCERF3.simulatedAnnealing.hpc.LogicTreePBSWriter.RunSites;
import scratch.peter.ucerf3.calc.UC3_CalcMPJ_Map;
import scratch.peter.ucerf3.calc.UC3_CalcMPJ_MapCompound;

public class MultiSolComparisonMapScriptGen {

	public static void main(String[] args) throws IOException {
		File localMainDir = new File("/home/kevin/OpenSHA/UCERF3/biasi_downsample_tests");
		RunSites site = RunSites.HPCC;
		File remoteMainDir = new File("/home/scec-02/kmilner/ucerf3/maps");
		
		String runName = "biasi-downsample-pga";
		
		List<File> solFiles = Lists.newArrayList();
		List<String> solNames = Lists.newArrayList();
		
		solFiles.add(new File(localMainDir,
				"FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate7.9_MMaxOff7.6_NoFix_SpatSeisU3_mean_sol.zip"));
		solNames.add("ref_branch_orig");
		solFiles.add(new File(localMainDir,
				"FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate7.9_MMaxOff7.6_NoFix_SpatSeisU3_VarFilteredRupsDistilledBoth_mean_sol.zip"));
		solNames.add("ref_branch_downsampled_inverted_both");
		solFiles.add(new File(localMainDir,
				"FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate7.9_MMaxOff7.6_NoFix_SpatSeisU3_VarFilteredRupsDistilledEnds_mean_sol.zip"));
		solNames.add("ref_branch_downsampled_inverted_ends");
		solFiles.add(new File(localMainDir,
				"FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate7.9_MMaxOff7.6_NoFix_SpatSeisU3_VarFilteredRupsDistilledStarts_mean_sol.zip"));
		solNames.add("ref_branch_downsampled_inverted_starts");
		
		int nodes = 15;
		int jobMins = 2*60; // TODO
		String resCode = "0.1";
		String imtCode = "GM0P00";
//		String threadsArg = "";
		// trailing space is important
		String threadsArg = "--threads 8 ";
		String imrCode = SourceIMR.WUS_FAULT_14.name();
		String gridCode = TestGrid.CA_RELM.name();
		IncludeBackgroundOption bgInclude = IncludeBackgroundOption.EXCLUDE;
		runName = LogicTreePBSWriter.df.format(new Date())+"-"+runName;
		
		File remoteDir = new File(remoteMainDir, runName);
		File writeDir = new File(localMainDir, runName);
		if (!writeDir.exists())
			writeDir.mkdir();
		
		List<File> classpath = LogicTreePBSWriter.getClasspath(remoteMainDir, remoteDir);
		
		FastMPJShellScriptWriter mpjWrite = new FastMPJShellScriptWriter(site.getJAVA_BIN(), site.getMaxHeapSizeMB(null),
					classpath, site.getMPJ_HOME());
		mpjWrite.setUseLaunchWrapper(true);
		BatchScriptWriter batchWrite = site.forBranch(null);
		
		for (int i=0; i<solFiles.size(); i++) {
			File solFile = solFiles.get(i);
			String solName = solNames.get(i);
			Preconditions.checkState(solFile.exists());
			
			Files.copy(solFile, new File(writeDir, solName+".zip"));
			File remoteSolFile = new File(remoteDir, solName+".zip");
			
			String className = UC3_CalcMPJ_Map.class.getName();
			String classArgs = threadsArg+" --deadlock "+remoteSolFile.getAbsolutePath()+" "+imrCode +" "
					+gridCode+" "+resCode+" "+imtCode+" "+bgInclude.name()+" "+remoteDir.getAbsolutePath()+" false false";
			
			File pbsFile = new File(writeDir, solName+".pbs");
			
			List<String> script = mpjWrite.buildScript(className, classArgs);
			
			batchWrite.writeScript(pbsFile, script, jobMins, nodes, site.getPPN(null), null);
		}
	}

}
