package scratch.peter.ucerf3.scripts;

import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.hpc.pbs.BatchScriptWriter;

public class HPCC_ScriptWriter extends BatchScriptWriter {
	
	String nodesAddition;
	
	HPCC_ScriptWriter() {
		this(null);
	}
	
	HPCC_ScriptWriter(String nodesAddition) {
		this.nodesAddition = nodesAddition;
	}

	@Override
	public List<String> getBatchHeader(int hours, int nodes,
			int ppn, String queue) {
		ArrayList<String> pbs = new ArrayList<String>();
		
		if (queue != null && !queue.isEmpty())
			pbs.add("#PBS -q "+queue);
		String dashL = "#PBS -l walltime="+hours+":00:00,nodes="+nodes;
		if (nodesAddition != null && !nodesAddition.isEmpty())
			dashL += ":"+nodesAddition;
		if (ppn > 0)
			dashL += ":ppn="+ppn;
		pbs.add(dashL);
		pbs.add("#PBS -V");
		pbs.add("");
		pbs.add("NEW_NODEFILE=\"/tmp/${USER}-hostfile-${PBS_JOBID}\"");
		pbs.add("echo \"creating PBS_NODEFILE: $NEW_NODEFILE\"");
		pbs.add("cat $PBS_NODEFILE | sort | uniq > $NEW_NODEFILE");
		pbs.add("export PBS_NODEFILE=$NEW_NODEFILE");
		pbs.add("");
		
		return pbs;
	}

}
