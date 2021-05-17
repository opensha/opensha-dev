package scratch.kevin;

import java.awt.Color;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Scanner;

import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.GraphWindow;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.base.Splitter;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class MPJTaskLogNodeFailure {
	
	private static List<List<String>> splitLog(List<String> logFile) {
		List<List<String>> logs = Lists.newArrayList();
		
		List<String> currentLog = null;
		for (String line : logFile) {
			if (line.contains("FastMPJ Runtime: Starting")) {
				if (currentLog != null)
					logs.add(currentLog);
				currentLog = Lists.newArrayList();
			}
			if (currentLog != null)
				currentLog.add(line);
		}
		if (currentLog != null && !currentLog.isEmpty())
			logs.add(currentLog);
		
		return logs;
	}
	
	private static boolean handleLog(List<String> log, Map<Integer, Integer> procFailuresMap, Map<String, Integer> hostFailuresMap) {
		Map<Integer, String> processHostsMap = Maps.newHashMap();
		
		boolean success = false;
		String lastWaitings = null;
		
		for (String line : log) {
			if (line.contains("Process 0 DONE!"))
				success = true;
			
			if (line.contains("Process") && line.contains("[") && line.contains("]")) {
				line = line.substring(line.indexOf("["), line.indexOf("]"));
				// split by spaces
				String[] split = line.split(" ");
				if (split.length != 4)
					continue;
				Integer procNum = Integer.parseInt(split[split.length-1]);
				if (!processHostsMap.containsKey(procNum)) {
					String host = split[1].replaceAll("\\(", "").replace("\\)", "");
					processHostsMap.put(procNum, host);
				}
			} else if (line.contains("not yet. waiting on:")) {
				lastWaitings = line;
			}
		}
		
		for (Integer proc : processHostsMap.keySet()) {
			if (!procFailuresMap.containsKey(proc))
				procFailuresMap.put(proc, 0);
		}
		
		for (String host : processHostsMap.values()) {
			if (!hostFailuresMap.containsKey(host))
				hostFailuresMap.put(host, 0);
		}
		
		if (!success) {
			// find the failures
			if (lastWaitings != null) {
				lastWaitings = lastWaitings.substring(lastWaitings.indexOf("waiting on:")+("waiting on: ").length());
				lastWaitings = lastWaitings.trim();
				if (lastWaitings.contains(" "))
					lastWaitings = lastWaitings.split(" ")[0];
				Iterable<String> procStrs = Splitter.on(",").split(lastWaitings);
				for (String procStr : procStrs) {
					Integer proc = Integer.parseInt(procStr);
					Preconditions.checkNotNull(proc);
					String host = processHostsMap.get(proc);
					procFailuresMap.put(proc, procFailuresMap.get(proc)+1);
					if (host == null)
						System.err.println("WARNING: Host never found for proc "+proc);
					else
						hostFailuresMap.put(host, hostFailuresMap.get(host)+1);
				}
			}
		}
		
		return success;
	}

	/**
	 * @param args
	 * @throws FileNotFoundException 
	 */
	public static void main(String[] args) throws FileNotFoundException {
		Map<Integer, Integer> procFailuresMap = Maps.newHashMap();
		Map<String, Integer> hostFailuresMap = Maps.newHashMap();
		File logDir = new File("/tmp/outs");
		
		int cnt = 0;
		int failures = 0;
		for (File file : logDir.listFiles()) {
			if (!file.getName().contains("pbs.o"))
				continue;
			StringBuilder text = new StringBuilder();
			String NL = System.getProperty("line.separator");
			Scanner scanner = new Scanner(
					new BufferedInputStream(new FileInputStream(file)));
			try {
				while (scanner.hasNextLine()){
					text.append(scanner.nextLine() + NL);
				}
			} finally{
				scanner.close();
			}
			List<String> fullList = Lists.newArrayList(text.toString().split("\n"));
			List<List<String>> logs = splitLog(fullList);
			
			for (List<String> log : logs) {
				boolean success = handleLog(log, procFailuresMap, hostFailuresMap);
				cnt++;
				if (!success)
					failures++;
			}
		}
		System.out.println(failures+"/"+cnt+" failed ("+getPercent(failures, cnt)+")");
		// make plot of process failure rates
		ArbitrarilyDiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
		for (Integer proc : procFailuresMap.keySet())
			func.set((double)proc, (double)procFailuresMap.get(proc));
		
		ArrayList<DiscretizedFunc> funcs = Lists.newArrayList();
		funcs.add(func);
		ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList(
				new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLUE));
		new GraphWindow(funcs, "Process # Failures", chars);
		
		List<String> hosts = Lists.newArrayList(hostFailuresMap.keySet());
		Collections.sort(hosts);
		System.out.println(hosts.size()+" unique hosts");
		ArbitrarilyDiscretizedFunc hostFunc = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<hosts.size(); i++) {
			String host = hosts.get(i);
			hostFunc.set((double)i, (double)hostFailuresMap.get(host));
		}
		
		funcs = Lists.newArrayList();
		funcs.add(hostFunc);
		new GraphWindow(funcs, "Host Failures", chars);
		
		EvenlyDiscretizedFunc hostDist = new EvenlyDiscretizedFunc(0d, (int)hostFunc.getMaxY()+1, 1d);
		for (int i=0; i<hosts.size(); i++) {
			String host = hosts.get(i);
			hostDist.add(hostFailuresMap.get(host), 1d);
		}
		
		funcs = Lists.newArrayList();
		funcs.add(hostDist);
		new GraphWindow(funcs, "Host Failure Distribution", chars);
	}
	
	private static String getPercent(int subVal, int tot) {
		double percent = 100d*(double)subVal / (double) tot;
		return (float)percent+" %";
	}

}
