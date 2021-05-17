package scratch.kevin.surfDistCache;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;
import org.opensha.sha.faultSurface.cache.SurfaceCachingPolicy;
import org.opensha.sha.faultSurface.cache.SurfaceCachingPolicy.CacheTypes;

import com.google.common.collect.Lists;

public class SurfDistCacheTestScriptWriter {
	
	public static void main(String[] args) throws IOException {
		int jobsPerConfig = 3;
		int numSites = 80;
		int[] cacheSizes = { 1, 2, 4, 8, 12 };
		int[] numThreads = { 1, 2, 4, 8, 12, 16 };
		CacheTypes[] forces = { null, CacheTypes.HYBRID };
		
		boolean noExpiration = true;
		SurfDistCacheTests.TestType type = SurfDistCacheTests.TestType.HAZARD_CURVE;
		String jobName = "hazard_5_hybrid_no_exp";
		
		File localDir = new File("/home/kevin/OpenSHA/dist_cache");
		File remoteDir = new File("/auto/scec-02/kmilner/dist_cache_tests");
		
		localDir = new File(localDir, jobName);
		if (!localDir.exists())
			localDir.mkdir();
		
		remoteDir = new File(remoteDir, jobName);
		if (!remoteDir.exists())
			remoteDir.mkdir();
		
		File solFile = new File("/home/scec-02/kmilner/ucerf3/inversion_compound_plots/"
				+ "2013_05_10-ucerf3p3-production-10runs/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip");
		
		int mins = 60;
		int nodes = 1;
		int ppn = 8;
		String queue = "nbns";
		
		List<File> classpath = Lists.newArrayList(new File(remoteDir, "OpenSHA_complete.jar"));
		
//		File javaBin = USC_HPCC_ScriptWriter.JAVA_BIN;
		File javaBin = new File("/usr/usc/jdk/1.7.0_25/bin/java");
		
		JavaShellScriptWriter javaWrite = new JavaShellScriptWriter(javaBin, 7000, classpath);
		USC_HPCC_ScriptWriter writer = new USC_HPCC_ScriptWriter("pe1950");
		
		if (noExpiration) {
			javaWrite.setProperty(SurfaceCachingPolicy.EXP_TIME_PROP, 0+"");
		}
		
		String className = SurfDistCacheTests.class.getName();
		
		for (int threads : numThreads) {
			for (int size : cacheSizes) {
				javaWrite.setProperty(SurfaceCachingPolicy.SIZE_PROP, size+"");
				
				for (CacheTypes force : forces) {
					if (force == null)
						javaWrite.clearProperty(SurfaceCachingPolicy.FORCE_TYPE);
					else
						javaWrite.setProperty(SurfaceCachingPolicy.FORCE_TYPE, force.name()+"");
					String argz = type.name()+" "+solFile.getAbsolutePath()+" "+threads+" "+numSites;
					for (int j=0; j<jobsPerConfig; j++) {
						File outputFile = new File(localDir, "threads_"+threads+"_size_"+size+"_force_"+force+"_run_"+j+".pbs");
						List<String> script = writer.buildScript(javaWrite.buildScript(className, argz), mins, nodes, ppn, queue);
						writer.writeScript(outputFile, script, mins, nodes, ppn, queue);
					}
				}
			}
		}
	}

}
