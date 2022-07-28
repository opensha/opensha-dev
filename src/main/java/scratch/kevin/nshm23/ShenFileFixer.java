package scratch.kevin.nshm23;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

public class ShenFileFixer {

	public static void main(String[] args) throws IOException {
		File inputFile = new File("/home/kevin/workspace/opensha/src/main/resources/data/erf/nshm23/def_models/"
//				+ "geodetic/2022_06_27/SHEN_BIRD-include_ghost_corr.txt");
				+ "geodetic/2022_06_27/SHEN_BIRD-no_ghost_corr.txt");
		
		File outputFile = new File(inputFile.getAbsolutePath()+".mod");
		
		FileWriter fw = new FileWriter(outputFile);
		
		for (String line : Files.readLines(inputFile, Charset.defaultCharset())) {
			line = line.trim();
			if (line.isBlank())
				continue;
			boolean header = line.startsWith("#");
			if (header)
				line = line.substring(1).trim();
			line.replaceAll("\t", " ");
			while (line.contains("  "))
				line = line.replaceAll("  ", " ");
			String[] split = line.split(" ");
			Preconditions.checkState(split.length == 9);
			String tmp = split[2];
			split[2] = split[3];
			split[3] = tmp;
			tmp = split[4];
			split[4] = split[5];
			split[5] = tmp;
			line = null;
			for (String str : split) {
				if (line == null)
					line = "";
				else
					line += "\t";
				line += str;
			}
			if (header)
				line = "# "+line;
			fw.write(line+"\n");
		}
		
		fw.close();
	}

}
