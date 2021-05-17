package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.dom4j.DocumentException;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;

import com.google.common.collect.Lists;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.utils.FaultSystemIO;

public class U3ERF_MemTest {

	public static void main(String[] args) throws IOException, DocumentException {
		FaultSystemSolution sol = FaultSystemIO.loadSol(new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/"
				+ "InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		
		Runtime runtime = Runtime.getRuntime();
		System.gc();
		double maxMemory = runtime.maxMemory()/1073741824d;
	    double allocatedMemory = runtime.totalMemory()/1073741824d;
	    double freeMemory = runtime.freeMemory()/1073741824d;
	    double usedMemory = allocatedMemory - freeMemory;
		double initialUsed = usedMemory;
		
		System.out.println("Initial used memory: "+initialUsed+" GB");
		
		List<FaultSystemSolutionERF> erfs = Lists.newArrayList();
		
		while (true) {
			System.gc();
			FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
			erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_PREF_BLEND);
			erf.updateForecast();
			
			maxMemory = runtime.maxMemory()/1073741824d;
		    allocatedMemory = runtime.totalMemory()/1073741824d;
		    freeMemory = runtime.freeMemory()/1073741824d;
		    usedMemory = allocatedMemory - freeMemory;
		    
		    erfs.add(erf);
		    double memEach = (usedMemory - initialUsed)/erfs.size();
		    System.out.println("Built "+erfs.size()+" ERFs.\tUsed "
		    		+(float)usedMemory+"/"+(float)allocatedMemory+" GB ("+(float)maxMemory+" tot)\t"+(float)memEach+" GB/ERF");
		}
	}

}
