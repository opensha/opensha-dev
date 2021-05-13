package scratch.kevin;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;

import com.google.common.collect.Lists;

public class MorganSTREC_ScriptGen {

	public static void main(String[] args) throws IOException {
		File strecBin = new File("/home/scec-00/kmilner/anaconda2/bin/getstrec_bulk.py");
		
		File inputDir = new File("/home/scec-00/kmilner/strec/morgan_map_2");
		File outputDir = new File(inputDir, "output");
		File writeDir = new File("/tmp/morgan_strec/pbs");
		
		int num = 109;
		int digits = ((num-1)+"").length();
		
		int mins = 96*60;
		int nodes = 1;
		int ppn = 20;
		String queue = "scec";
		
		int threadsPerJob = 10;
		
		USC_HPCC_ScriptWriter pbsWrite = new USC_HPCC_ScriptWriter();
		
		int jobCount = 0;
		int catCount = 0;
		
		while (catCount < num) {
			List<String> script = Lists.newArrayList();
			script.add("#!/bin/bash");
			script.add("");
			script.add("mkdir -p "+outputDir.getAbsolutePath());
			script.add("");
			for (int i=0; i<threadsPerJob; i++) {
				if (catCount < num) {
					File inputFile = new File(inputDir, "catalog"+catCount+".in");
					File outputFile = new File(outputDir, "catalog"+catCount+".out");
					
					script.add("echo \"input file: "+inputFile.getAbsolutePath()+"\"");
					script.add("echo \"output file: "+outputFile.getAbsolutePath()+"\"");
					script.add("");
					String command = strecBin.getAbsolutePath()+"  --batch-input "+inputFile.getAbsolutePath()
							+" --csv-out > "+outputFile.getAbsolutePath();
					if (threadsPerJob > 1)
						command += " &";
					script.add(command);
					script.add("");
					script.add("");
				}
				
				catCount++;
			}
			
			if (threadsPerJob > 1) {
				script.add("echo \"waiting\"");
				script.add("wait");
				script.add("echo \"DONE!\"");
			}
			
			String runNum = (jobCount++)+"";
			while (runNum.length() < digits)
				runNum = "0"+runNum;
			
			File opbsFile = new File(writeDir, "job"+runNum+".pbs");
			pbsWrite.writeScript(opbsFile, script, mins, nodes, ppn, queue);
		}
		for (int i=0; i<num; i++) {
			
			
			
		}
	}

}
