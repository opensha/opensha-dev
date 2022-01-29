package scratch.kevin;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.google.common.base.Charsets;
import com.google.common.base.Preconditions;
import com.google.common.io.Files;

public class Apache2HeaderDelete {

	public static void main(String[] args) throws IOException {
//		File parentDir = new File("/tmp/license_strip/opensha-cybershake/src");
		File parentDir = new File("/home/kevin/workspace/opensha-oaf/src");
		int found = processDir(parentDir, true);
		System.out.println("Found "+found+" matches");
	}
	
	private static int processDir(File dir, boolean dryRun) throws IOException {
		int numFound = 0;
		for (File file : dir.listFiles()) {
			String name = file.getName();
			if (name.startsWith("."))
				continue;
			if (file.isDirectory()) {
				numFound += processDir(file, dryRun);
				continue;
			}
			if (name.toLowerCase().endsWith(".java")) {
				System.out.println("Checking: "+file.getAbsolutePath());
				List<String> lines = Files.readLines(file, Charsets.UTF_8);
				if (lines.isEmpty())
					continue;
				
				List<String> extraHeaderLines = new ArrayList<>();
				
				String header = null;
				boolean found = false;
				
				int curBlockStart = -1;
				int curBlockEnd = -1;
				boolean insideBlock = false;
				
				int apacheBlockStart = -1;
				int apacheBlockEnd = -1;
				boolean insideApacheBlock = false;
				
				int numBlocks = 0;
				int packageIndex = -1;
				for (int i=0; i<lines.size(); i++) {
					String line = lines.get(i);
					line = line.trim();
					if (line.isBlank())
						continue;
					if (line.startsWith("package")) {
						packageIndex = i;
						break;
					}
					if (header == null)
						header = "";
					else
						header += "\n";
					header += line;
					
					if (line.contains("OpenSHA.org in partnership with")) {
						found = true;
						Preconditions.checkState(insideBlock);
						insideApacheBlock = true;
						apacheBlockStart = curBlockStart;
					}
					
					if (line.startsWith("/**")) {
						Preconditions.checkState(!insideBlock);
						insideBlock = true;
						numBlocks++;
						curBlockStart = i;
					} else if (line.startsWith("**") && line.endsWith("**/")) {
						insideBlock = false;
						curBlockEnd = i;
						if (insideApacheBlock) {
							insideApacheBlock = false;
							apacheBlockEnd = i;
						}
					} else if (!line.startsWith("*")) {
						// header line
						extraHeaderLines.add(line);
					}
				}
				if (found) {
					int len = 1 + apacheBlockEnd - apacheBlockStart;
					System.out.println("\tFOUND header with "+len+" lines. "+numBlocks+" total blocks, "
							+extraHeaderLines.size()+" extra header lines");
					Preconditions.checkState(numBlocks == 1, "We have "+numBlocks+" blocks, fix manually.");
					Preconditions.checkState(extraHeaderLines.isEmpty(),
							"We have "+extraHeaderLines.size()+" extra header lines, fix manually.");
					Preconditions.checkState(len >= 17 && len <= 19, "Expected 17-19 header lines but have "+len);
					Preconditions.checkState(packageIndex > 0, "Never found package declaration");
					
					if (!dryRun) {
						System.out.println("\tStripping header!");
						List<String> linesMod = lines.subList(packageIndex, lines.size());
						FileWriter fw = new FileWriter(file, Charsets.UTF_8);
						for (String line : linesMod)
							fw.write(line+"\n");
						fw.close();
					}
					numFound++;
				}
			}
		}
		return numFound;
	}

}
