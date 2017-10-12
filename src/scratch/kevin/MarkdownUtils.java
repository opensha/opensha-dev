package scratch.kevin;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import org.commonmark.Extension;
import org.commonmark.ext.gfm.tables.TablesExtension;
import org.commonmark.ext.heading.anchor.HeadingAnchorExtension;
import org.commonmark.node.Node;
import org.commonmark.parser.Parser;
import org.commonmark.renderer.html.HtmlRenderer;

import com.google.common.base.CharMatcher;
import com.google.common.base.Preconditions;
import com.google.common.io.Files;

public class MarkdownUtils {
	
	public static class TableBuilder {
		
		private LinkedList<String> lines;
		
		private List<String> curLine;
		
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
		
		public TableBuilder initNewLine() {
			if (curLine != null && !curLine.isEmpty())
				finalizeLine();
			curLine = new ArrayList<>();
			return this;
		}
		
		public TableBuilder addColumn(String val) {
			if (curLine == null)
				initNewLine();
			curLine.add(val);
			return this;
		}
		
		public TableBuilder finalizeLine() {
			Preconditions.checkState(curLine != null && !curLine.isEmpty());
			addLine(curLine);
			curLine = null;
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
	
	public static String getTitle(File markdownPage) throws IOException {
		for (String line : Files.readLines(markdownPage, Charset.defaultCharset())) {
			if (line.startsWith("#"))
				return line.substring(1).trim();
		}
		return null;
	}
	
	private static final CharMatcher ALNUM = CharMatcher.inRange('a', 'z').or(CharMatcher.inRange('A', 'Z'))
			  .or(CharMatcher.inRange('0', '9')).or(CharMatcher.is('-'));
	
	public static List<String> buildTOC(List<String> lines, int minLevel) {
		LinkedList<String> toc = new LinkedList<>();
		
		for (String line : lines) {
			if (line.startsWith("#")) {
				String headerPart = line.substring(0, line.lastIndexOf('#')+1);
				int level = headerPart.length();
				if (level >= minLevel) {
					String tocLine = "";
					while ((level > minLevel)) {
						tocLine += "  ";
						level--;
					}
					String title = line.substring(headerPart.length()).trim();
					String anchor = getAnchorName(title);
					tocLine += "* ["+title+"](#"+anchor+")";
					toc.add(tocLine);
				}
			}
		}
		return toc;
	}
	
	public static String getAnchorName(String heading) {
		while (heading.startsWith("#"))
			heading = heading.substring(1);
		heading = heading.trim();
		return ALNUM.retainFrom(heading.toLowerCase().replaceAll(" ", "-")).toLowerCase();
	}
	
	public static void writeHTML(List<String> lines, File outputFile) throws IOException {
		StringBuilder str = new StringBuilder();
		for (String line : lines)
			str.append(line).append("\n");
		writeHTML(str.toString(), outputFile);
	}
	public static void writeHTML(String markdown, File outputFile) throws IOException {
		List<Extension> extensions = Arrays.asList(TablesExtension.create(), HeadingAnchorExtension.create());
		Parser parser = Parser.builder().extensions(extensions).build();
		Node document = parser.parse(markdown);
		HtmlRenderer renderer = HtmlRenderer.builder().extensions(extensions).build();
		FileWriter fw = new FileWriter(outputFile);
		renderer.render(document, fw);
		fw.close();
	}

}
