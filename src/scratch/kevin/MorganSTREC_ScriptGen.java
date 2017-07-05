package scratch.kevin;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;

import com.google.common.collect.Lists;

public class MorganSTREC_ScriptGen {

	public static void main(String[] args) throws IOException {
		File strecBin = new File("/home/scec-00/kmilner/anaconda2/bin/getstrec_bulk.py");
		
		File inputDir = new File("/home/scec-00/kmilner/strec/morgan_map");
		File outputDir = new File(inputDir, "output");
		File writeDir = new File("/tmp/morgan_strec/pbs");
		
		int num = 143;
		int digits = ((num-1)+"").length();
		
		int mins = 48*60;
		int nodes = 1;
		int ppn = 8;
		String queue = "nbns";
		
		USC_HPCC_ScriptWriter pbsWrite = new USC_HPCC_ScriptWriter();
		
		for (int i=0; i<num; i++) {
			File inputFile = new File(inputDir, "catalog"+i+".in");
			File outputFile = new File(outputDir, "catalog"+i+".out");
			
			List<String> script = Lists.newArrayList();
			script.add("#!/bin/bash");
			script.add("");
			script.add("echo \"input file: "+inputFile.getAbsolutePath()+"\"");
			script.add("echo \"output file: "+outputFile.getAbsolutePath()+"\"");
			script.add("");
			script.add("mkdir -p "+outputDir.getAbsolutePath());
			script.add("");
			script.add(strecBin.getAbsolutePath()+"  --batch-input "+inputFile.getAbsolutePath()
				+" --csv-out > "+outputFile.getAbsolutePath());
			
			String runNum = i+"";
			while (runNum.length() < digits)
				runNum = "0"+runNum;
			
			File opbsFile = new File(writeDir, "job"+runNum+".pbs");
			pbsWrite.writeScript(opbsFile, script, mins, nodes, ppn, queue);
		}
	}

}
