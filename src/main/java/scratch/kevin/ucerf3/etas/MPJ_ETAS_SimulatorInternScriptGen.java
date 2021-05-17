package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.StampedeScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;

import com.google.common.collect.Lists;

public class MPJ_ETAS_SimulatorInternScriptGen {

	public static void main(String[] args) throws IOException {
		File localDir = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/2014_interns");
		
		boolean stampede = false;
		
//		1. 99316- Bartlett 7.33
//		2. 90852- Calaveras 6.92
//		3. 96759- Concord 6.77
//		4. 93119- Diablo 6.72
//		5. 97691- Great Valley (3) 6.97
//		6. 97615- Great Valley (4a) 6.64
//		7. 97639- Great Valley (4b) 6.67
//		8. 100934- Great Valley (5) 6.48*
//		9. 96886- Greenville 7.05
//		10. 106140- Hayward RC 7.08
//		11. 101292- Hayward HN+HS 7.07 
//		12. 98144- Hunting Creek/Berryessa  6.44**
//		13. 119732- North San Andreas 7.49 
//		14. 73385- North San Andreas 7.21
//		15. 64481- North San Andreas 7.15
//		16. 251623- Ortegalita 7.09
//		17. 117665- San Gergorio 7.54
//		18. 93905- West Napa 6.74
//		int[] scenarios = {99316, 90852, 96759, 93119, 97691, 97615, 97639, 100934, 96886, 106140,
//				101292, 98144, 119732, 73385, 64481, 251623, 117665, 93905};
		
//		1. 225699, Newport-Inglewood Fault
//		2. 183888, Sur-Nacimiento Fault 
//		3. 224048, S San Andreas Fault 
//		4. 241894, Masson Hill & Compton Thrust Fault 
//		5. 226852, Malibu Coast & Santa Monica 
		int[] scenarios = {225699, 183888, 224048, 241894, 226852};
		
		boolean timeIndep = false;
		int numSims = 500;
		
		int memGigs;
		int mins = 24*60;
		int nodes = 20;
		int ppn;
		if (stampede)
			ppn = 16;
		else
			ppn = 8;
		String queue = null;
		
		File remoteDir, remoteSolFile;
		FastMPJShellScriptWriter mpjWrite;
		BatchScriptWriter pbsWrite;
		
		if (stampede) {
			memGigs = 26;
			remoteDir = new File("/work/00950/kevinm/ucerf3/etas_sim");
			remoteSolFile = new File("/work/00950/kevinm/ucerf3/inversion/compound_plots/2013_05_10-ucerf3p3-production-10runs/"
					+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip");
			mpjWrite = new FastMPJShellScriptWriter(StampedeScriptWriter.JAVA_BIN, memGigs*1024,
					null, StampedeScriptWriter.FMPJ_HOME);
			pbsWrite = new StampedeScriptWriter();
		} else {
			memGigs = 9;
			remoteDir = new File("/home/scec-02/kmilner/ucerf3/etas_sim");
			remoteSolFile = new File("/home/scec-02/kmilner/ucerf3/inversion_compound_plots/"
					+ "2013_05_10-ucerf3p3-production-10runs/"
					+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip");
			mpjWrite = new FastMPJShellScriptWriter(USC_HPCC_ScriptWriter.JAVA_BIN, memGigs*1024,
					null, USC_HPCC_ScriptWriter.FMPJ_HOME);
			pbsWrite = new USC_HPCC_ScriptWriter();
		}
		File remoteInternDir = new File(remoteDir, "2014_interns");
		
		List<File> classpath = new ArrayList<File>();
		classpath.add(new File(remoteDir, "commons-cli-1.2.jar"));
		classpath.add(new File(remoteInternDir, "OpenSHA_complete.jar"));
		mpjWrite.setClasspath(classpath);
		
		for (int scenario : scenarios) {
			String jobName = new SimpleDateFormat("yyyy_MM_dd").format(new Date())+"-"+scenario;
			if (timeIndep)
				jobName += "-indep";
			
			File localJobDir = new File(localDir, jobName);
			if (!localJobDir.exists())
				localJobDir.mkdir();
			File remoteJobDir = new File(remoteInternDir, jobName);
			
			File pbsFile = new File(localJobDir, jobName+".pbs");
			
			String argz = "--min-dispatch 1 --max-dispatch 1 --num "+numSims+" --sol-file "+remoteSolFile.getAbsolutePath();
			argz += " --trigger-rupture-id "+scenario;
			if (timeIndep)
				argz += " --indep";
			argz += " "+remoteDir.getAbsolutePath()+" "+remoteJobDir.getAbsolutePath();
			
			List<String> script = mpjWrite.buildScript(MPJ_ETAS_Simulator.class.getName(), argz);
			
			script = pbsWrite.buildScript(script, mins, nodes, ppn, queue);
			pbsWrite.writeScript(pbsFile, script);
		}
	}

}
