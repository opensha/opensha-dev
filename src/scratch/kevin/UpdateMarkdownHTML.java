package scratch.kevin;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.List;

import org.opensha.commons.util.MarkdownUtils;

import com.google.common.io.Files;

public class UpdateMarkdownHTML {
	
	public static void main(String[] args) throws IOException {
		File baseDir;
		if (args.length == 1)
			baseDir = new File(args[0]);
		else
			baseDir = new File("/data-0/kevin/markdown/cybershake-analysis");
		updateCSS(baseDir);
	}
	
	private static void updateCSS(File dir) throws IOException {
		for (File file : dir.listFiles()) {
			if (file.isDirectory())
				updateCSS(file);
			if (file.getName().equals("README.md")) {
				System.out.println("Updating "+file.getAbsolutePath());
				List<String> lines = Files.readLines(file, Charset.defaultCharset());
				MarkdownUtils.writeReadmeAndHTML(lines, dir);
			}
		}
	}

}
