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
		
		private List<String[]> lines;
		
		private List<String> curLine;
		
		private TableBuilder() {
			lines = new LinkedList<>();
		}
		
		public TableBuilder addLine(List<String> vals) {
			return addLine(vals.toArray(new String[vals.size()]));
		}
		
		public TableBuilder addLine(String... vals) {
			lines.add(vals);
			
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
		
		public TableBuilder wrap(int maxDataCols, int headerCols) {
			if (curLine != null)
				finalizeLine();
			int curDataCols = lines.get(0).length - headerCols;
			if (curDataCols <= maxDataCols)
				return this;
			int numWraps = (int)Math.ceil((double)curDataCols / (double)maxDataCols);
			int newDataCols = (int)Math.ceil((double)curDataCols/(double)numWraps);
			
			System.out.println("Wrapping data from "+lines.size()+"x"+curDataCols+" to "+(lines.size()*numWraps)+"x"+newDataCols);
			
			List<String[]> newLines = new ArrayList<>(lines.size()*numWraps);
			
			// init new lines with headers if necessary
			for (int i=0; i<lines.size()*numWraps; i++) {
				String[] newLine = new String[headerCols+newDataCols];
				for (int h=0; h<headerCols; h++)
					newLine[h] = lines.get(i % numWraps)[h];
				newLines.add(newLine);
			}
			
			// fill in data
			for (int i=0; i<lines.size(); i++) {
				for (int c=0; c<curDataCols; c++) {
					int row = i + (c/newDataCols)*lines.size();
					int col = headerCols + c%newDataCols;
					newLines.get(row)[col] = lines.get(i)[c];
				}
			}
			
			this.lines = newLines;
			
			return this;
		}
		
		public List<String> build() {
			Preconditions.checkState(lines.size() >= 1);
			
			List<String> strings = new ArrayList<>(lines.size()+1);
			
			for (int i=0; i<lines.size(); i++) {
				strings.add(tableLine(lines.get(i)));
				if (i == 0)
					strings.add(generateTableDashLine(lines.get(i).length));
			}
			
			return strings;
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
