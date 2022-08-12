package scratch.kevin.nshm23;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.io.Files;

public class DefModelTraceFlipper {

	public static void main(String[] args) throws IOException {
		File inputDir = new File("/home/kevin/workspace/opensha/src/main/resources/data/erf/nshm23/def_models/geodetic/fm_v1p4");
		File outputDir = new File("/home/kevin/workspace/opensha/src/main/resources/data/erf/nshm23/def_models/geodetic/fm_v2");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		int[] traceFlipIDs = { 2922, 1098 };
		
		for (File file : inputDir.listFiles()) {
			String name = file.getName();
			if (!file.isFile() || !name.endsWith(".txt"))
				continue;
			if (!name.startsWith("EVANS"))
				continue;
			System.out.println("Processing "+name);
			
			List<String> regLines = new ArrayList<>();
			List<String> flipBuffer = null;
			int curFlipID = -1;
			
			for (String line : Files.readLines(file, Charset.defaultCharset())) {
				line = line.trim();
				while (line.contains("  "))
					line = line.replaceAll("  ", " ");
				line = line.replaceAll(" ", "\t");
				while (line.contains("\t\t"))
					line = line.replaceAll("\t\t", "\t");
				Preconditions.checkState(!line.contains(" "));
				if (line.startsWith("#\t"))
					line = "#"+line.substring(2);
				
				int matchFlipID = -1;
				for (int flipID : traceFlipIDs) {
					if (line.startsWith(flipID+"\t")) {
						Preconditions.checkState(matchFlipID < 0);
						matchFlipID = flipID;
					}
				}
				
				if (curFlipID >= 0) {
					// we're working on a flip
					if (matchFlipID >= 0) {
						// continuation
						Preconditions.checkState(matchFlipID == curFlipID);
						flipBuffer.add(line);
					} else {
						// done
						regLines.addAll(flip(flipBuffer, curFlipID));
						flipBuffer = null;
						curFlipID = -1;
						regLines.add(line);
					}
				} else {
					// we weren't working on a flip;
					if (matchFlipID >= 0) {
						// but we are now
						curFlipID = matchFlipID;
						flipBuffer = new ArrayList<>();
						flipBuffer.add(line);
					} else {
						// and we still aren't
						regLines.add(line);
					}
				}
			}
			if (flipBuffer != null)
				regLines.addAll(flip(flipBuffer, curFlipID));
			
			FileWriter fw = new FileWriter(new File(outputDir, name));
			for (String line : regLines)
				fw.write(line+"\n");
			fw.close();
		}
	}
	
	private static List<String> flip(List<String> flipBuffer, int flipID) {
		Collections.reverse(flipBuffer);
		for (int i=0; i<flipBuffer.size(); i++) {
			int origID = flipBuffer.size()-i;
			int newID = i+1;
			String line = flipBuffer.get(i);
			
			String[] split = line.split("\t");
			Preconditions.checkState(split.length == 9);
			
			// overwrite minisection id
			split[1] = newID+"";
			
			// reverse trace locations
			String tmpLat = split[2];
			String tmpLon = split[3];
			
			split[2] = split[4];
			split[3] = split[5];
			split[4] = tmpLat;
			split[5] = tmpLon;
			
			String newLine = Joiner.on("\t").join(split);
			System.out.println("Processing line for "+flipID+", "+origID+" -> "+newID);
			System.out.println("\tPREV: "+line);
			System.out.println("\tNEW:  "+newLine);
			flipBuffer.set(i, newLine);
		}
		return flipBuffer;
	}

}
