package scratch.kevin;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

public class JMapDeltaReader {
	
	private static class JMapData { 
		private long instances;
		private long bytes;
		private String name;
		
		public JMapData(long instances, long bytes, String name) {
			super();
			this.instances = instances;
			this.bytes = bytes;
			this.name = name;
		}
	}
	
	private static DecimalFormat pDF = new DecimalFormat("0.00%");
	
	private static class JMapDiff { 
		private JMapData data1;
		private JMapData data2;
		
		private long instDelta;
		private double sizeDelta;
		
		public JMapDiff(JMapData data1, JMapData data2) {
			super();
			this.data1 = data1;
			this.data2 = data2;
			instDelta = data2.instances - data1.instances;
			sizeDelta = byteToMB(data2.bytes - data1.bytes);
		}
		
		public String toString() {
			String str = data1.name+":";
			str += "\n\tInstances: "+data1.instances+" -> "+data2.instances+";\tdelta="
					+(instDelta)+" ("+pDF.format((double)instDelta/(double)data1.instances)+")";
			str += "\n\tSize: "+(float)byteToMB(data1.bytes)+" -> "+(float)byteToMB(data2.bytes)+";\tdelta="
					+(float)sizeDelta+" ("+pDF.format((double)sizeDelta/byteToMB(data1.bytes))+")";
			return str;
		}
	}
	
	private static double byteToMB(long bytes) {
		double kb = (double)bytes / 1024d;
		return kb / 1024d;
	}
	
	private static Map<String, JMapData> load(File file) throws IOException {
		Map<String, JMapData> ret = new HashMap<>();
		for (String line : Files.readLines(file, Charset.defaultCharset())) {
			line = line.trim();
			
			if (line.isBlank())
				continue;
			if (!Character.isDigit(line.charAt(0)) || !line.contains(":"))
				continue;
			
			while (line.contains("  "))
				line = line.replaceAll("  ", " ");
			line = line.trim();
			line = line.replaceAll(" ", "\t");
			System.out.println(line);
			String[] split = line.split("\t");
			Preconditions.checkState(split.length >= 4);
			
			long instances = Long.parseLong(split[1]);
			long bytes = Long.parseLong(split[2]);
			
			String name = split[3];
			for (int i=4; i<split.length; i++)
				name += " "+split[i];
			ret.put(name, new JMapData(instances, bytes, name));
		}
		return ret;
	}
	
	private static double totSize(Collection<JMapData> datas) {
		long bytes = 0l;
		for (JMapData data : datas)
			bytes += data.bytes;
		return byteToMB(bytes);
	}
	
	public static void main(String[] args) throws IOException {
//		Map<String, JMapData> orig = load(new File("/tmp/jmap.txt"));
		Map<String, JMapData> orig = load(new File("/tmp/jmap8.txt"));
		Map<String, JMapData> last = load(new File("/tmp/jmap9.txt"));
		
		List<JMapDiff> diffs = new ArrayList<>();
		for (String key : last.keySet())
			if (orig.containsKey(key))
				diffs.add(new JMapDiff(orig.get(key), last.get(key)));
		
		int num = 10;
		System.out.println("Instance Counts Increases: ");
		diffs.sort(new Comparator<JMapDiff>() {

			@Override
			public int compare(JMapDiff o1, JMapDiff o2) {
				return Long.compare(o2.instDelta, o1.instDelta);
			}
		});
		for (int i=0; i<num && i < diffs.size(); i++)
			System.out.println(i+". "+diffs.get(i));
		System.out.println("Size Increases: ");
		diffs.sort(new Comparator<JMapDiff>() {

			@Override
			public int compare(JMapDiff o1, JMapDiff o2) {
				return Double.compare(o2.sizeDelta, o1.sizeDelta);
			}
		});
		for (int i=0; i<num && i < diffs.size(); i++)
			System.out.println(i+". "+diffs.get(i));
		
		double origMB = totSize(orig.values());
		double lastMB = totSize(last.values());
		
		System.out.println();
		System.out.println("Orig size: "+(int)origMB);
		System.out.println("Final size: "+(int)lastMB);
		System.out.println("Delta: "+(int)(lastMB - origMB)+" ("+pDF.format((lastMB-origMB)/origMB)+")");
	}

}
