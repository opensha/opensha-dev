package scratch.kevin;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import org.commonmark.Extension;
import org.commonmark.ext.gfm.tables.TablesExtension;
import org.commonmark.node.Node;
import org.commonmark.parser.Parser;
import org.commonmark.renderer.html.HtmlRenderer;

import com.google.common.base.Preconditions;

public class MarkdownUtils {
	
	public static class TableBuilder {
		
		private LinkedList<String> lines;
		
		private TableBuilder() {
			lines = new LinkedList<>();
		}
		
		public TableBuilder addLine(List<String> vals) {
			return addLine(vals.toArray(new String[vals.size()]));
		}
		
		public TableBuilder addLine(String... vals) {
			lines.add(tableLine(vals));
			if (lines.size() == 1)
				lines.add(generateTableDashLine(vals.length));
			
			return this;
		}
		
		public List<String> build() {
			Preconditions.checkState(lines.size() > 1);
			return lines;
		}
	}
	
	public static TableBuilder tableBuilder() {
		return new TableBuilder();
	}
	
	private static String generateTableDashLine(int numVals) {
		Preconditions.checkState(numVals > 1);
		String[] vals = new String[numVals];
		for (int i=0; i<vals.length; i++)
			vals[i] = "-----";
		return tableLine(vals).replaceAll(" ", "");
	}
	
	private static String tableLine(String[] vals) {
		Preconditions.checkState(vals.length > 1);
		StringBuilder line = new StringBuilder().append("| ");
		for (int i=0; i<vals.length; i++) {
			if (i > 0)
				line.append(" | ");
			line.append(vals[i]);
		}
		line.append(" |");
		return line.toString();
	}
	
	public static void writeReadmeAndHTML(List<String> lines, File outputDir) throws IOException {
		File mdFile = new File(outputDir, "README.md");
		// write markdown
		FileWriter fw = new FileWriter(mdFile);
		StringBuilder str = new StringBuilder();
		for (String line : lines) {
			fw.write(line+"\n");
			str.append(line).append("\n");
		}
		fw.close();

		// write html
		File htmlFile = new File(outputDir, "index.html");
		writeHTML(str.toString(), htmlFile);
	}
	
	public static void writeHTML(List<String> lines, File outputFile) throws IOException {
		StringBuilder str = new StringBuilder();
		for (String line : lines)
			str.append(line).append("\n");
		writeHTML(str.toString(), outputFile);
	}
	public static void writeHTML(String markdown, File outputFile) throws IOException {
		List<Extension> extensions = Arrays.asList(TablesExtension.create());
		Parser parser = Parser.builder().extensions(extensions).build();
		Node document = parser.parse(markdown);
		HtmlRenderer renderer = HtmlRenderer.builder().extensions(extensions).build();
		FileWriter fw = new FileWriter(outputFile);
		renderer.render(document, fw);
		fw.close();
	}

}
