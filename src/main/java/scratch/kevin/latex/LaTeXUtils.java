	package scratch.kevin.latex;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.commons.lang3.exception.ExceptionUtils;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.FileUtils;

import com.google.common.base.Preconditions;
import com.google.common.io.Files;

public class LaTeXUtils {

	private static final String VAR_PREFIX = "dynval";

	/**
	 * Defines a new value command. The command will start with the common prefix, dynval, then the given name.
	 * This is to make them easy to find, and potentially programatically replace when building the final document.
	 * 
	 * The name must only contain letters. The value will be escaped using {@link #escapeLaTeX(String)}.
	 * 
	 * @param name
	 * @param value
	 * @return variable command definition
	 * @throws IllegalStateException if the name contains any spaces, dashes, underscores, numbers, or special characters
	 */
	public static String defineValueCommand(String name, String value) {
		return defineValueCommand(name, value, true);
	}
	
	/**
	 * Defines a new value command. The command will start with the common prefix, dynval, then the given name.
	 * This is to make them easy to find, and potentially programatically replace when building the final document.
	 * 
	 * The name must only contain letters. The value will be escaped using {@link #escapeLaTeX(String)} if escapeValue is true
	 * 
	 * @param name
	 * @param value
	 * @param escapeValue
	 * @return
	 */
	public static String defineValueCommand(String name, String value, boolean escapeValue) {
		String sanitizedName = sanitizeForCommandName(name);

		Preconditions.checkState(sanitizedName.equals(name),
				"Only letters and numbers are allowed in latex commands: %s (sanitized: %s)", name, sanitizedName);

		if (escapeValue)
			value = escapeLaTeX(value);
		
		return "\\newcommand{\\"+VAR_PREFIX+name+"}{"+value+"}";
	}
	
	public static String sanitizeForCommandName(String name) {
		return name.replaceAll("[^a-zA-Z]", "");
	}

	private static final Map<Character, String> LATEX_ESCAPE_MAP = new HashMap<>();

	static {
		LATEX_ESCAPE_MAP.put('&', "\\&");
		LATEX_ESCAPE_MAP.put('%', "\\%");
		LATEX_ESCAPE_MAP.put('$', "\\$");
		LATEX_ESCAPE_MAP.put('#', "\\#");
		LATEX_ESCAPE_MAP.put('_', "\\_");
		LATEX_ESCAPE_MAP.put('{', "\\{");
		LATEX_ESCAPE_MAP.put('}', "\\}");
		LATEX_ESCAPE_MAP.put('~', "\\textasciitilde{}");
		LATEX_ESCAPE_MAP.put('^', "\\textasciicircum{}");
		LATEX_ESCAPE_MAP.put('\\', "\\textbackslash{}");
	}

	public static String escapeLaTeX(String input) {
		if (input == null) {
			return null;
		}
		StringBuilder escaped = new StringBuilder();
		boolean isEscaped = false;

		for (int i = 0; i < input.length(); i++) {
			char c = input.charAt(i);

			// Check if the current character is already escaped
			if (c == '\\' && i + 1 < input.length() && LATEX_ESCAPE_MAP.containsKey(input.charAt(i + 1))) {
				isEscaped = true;
			} else if (LATEX_ESCAPE_MAP.containsKey(c) && !isEscaped) {
				escaped.append(LATEX_ESCAPE_MAP.get(c));
				continue;
			} else {
				isEscaped = false;
			}

			escaped.append(c);
		}
		return escaped.toString();
	}
	
	public static String numberExpFormat(Number number) {
		return numberExpFormat(numberToString(number));
	}
	
	public static String numberExpFormatSigFigs(Number number, int sigFigs) {
		return numberExpFormat(numberToString(DataUtils.roundSigFigs(number.doubleValue(), sigFigs)));
	}
	
	public static String numberExpFormatFixedDecimal(Number number, int decimalPlaces) {
		return numberExpFormat(numberToString(DataUtils.roundFixed(number.doubleValue(), decimalPlaces)));
	}
	
	private static synchronized String numberToString(Number number) {
		String str = number+"";
		if (!str.toLowerCase().contains("e"))
			str = nonExpOptionalDF.format(number);
		return str;
	}

	private static final DecimalFormat nonExpOptionalDF = new DecimalFormat("0.###########");
	
	public static String numberExpFormat(String number) {
		if (number == null || number.isEmpty()) {
			throw new IllegalArgumentException("Input number cannot be null or empty");
		}

		// Check if the input contains 'e' or 'E' for scientific notation
		int index = Math.max(number.indexOf('e'), number.indexOf('E'));
		if (index == -1) {
			// If no 'e' or 'E', return the number as is
			return number;
		}

		// Split the number into base and exponent
		String base = number.substring(0, index);
		String exponent = number.substring(index + 1);

		// Return the formatted LaTeX string
		return String.format("$%s \\times 10^{%s}$", base, exponent);
	}
	
	private static final DecimalFormat groupedIntDF = new DecimalFormat("0");
	static {
		groupedIntDF.setGroupingUsed(true);
		groupedIntDF.setGroupingSize(3);
	}
	
	public static String groupedIntNumber(Number number) {
		return groupedIntDF.format(number);
	}
	
	private static final DecimalFormat percentDF = new DecimalFormat("0.##########%");
	
	public static String numberAsPercent(Number number, int decimalPlaces) {
		return numberExpFormatFixedDecimal(number, decimalPlaces)+LATEX_ESCAPE_MAP.get('%');
	}
	
	public static final String BEGIN_DYNVAL_INCLUDES_HEADER = "% begin dynval includes";
	public static final String END_DYNVAL_INCLUDES_HEADER = "% end dynval includes";
	
	public static List<String> embedDynvalIncludes(List<String> inLines, File refDir) throws IOException {
		boolean processingIncludes = false;
		Map<String, String> dynVals = null;
		
		List<String> ret = new ArrayList<>(inLines.size());
		for (String line : inLines) {
			if (dynVals == null) {
				if (line.trim().equalsIgnoreCase(BEGIN_DYNVAL_INCLUDES_HEADER)) {
					System.out.println("Found '"+BEGIN_DYNVAL_INCLUDES_HEADER+"'");
					processingIncludes = true;
					dynVals = new HashMap<>();
					continue;
				}
			}
			if (processingIncludes) {
				line = line.trim();
				if (!line.isBlank()) {
					if (line.equalsIgnoreCase(END_DYNVAL_INCLUDES_HEADER)) {
						System.out.println("Found '"+END_DYNVAL_INCLUDES_HEADER+"'");
						System.out.println("Processed "+dynVals.size()+" total dynamic values");
						processingIncludes = false;
						continue;
					} else if (line.startsWith("%")) {
						continue;
					} else {
						Preconditions.checkState(line.startsWith("\\include{"),
								"Inside dynval include block but line doesn't start with an \\include{: %s", line);
						String path = line.substring(line.indexOf('{')+1);
						Preconditions.checkState(path.contains("}"), "Line doesn't close '}': %s", line);
						path = path.substring(0, path.indexOf('}'));
						File includeFile = new File(refDir, path);
						if (!includeFile.exists() && !path.toLowerCase().endsWith(".tex"))
							includeFile = new File(refDir, path+".tex");
						System.out.println("\tProcessing "+includeFile.getName());
						Preconditions.checkState(includeFile.exists(), "Include path not found. raw='%s', abs='%s'",
								path, includeFile.getAbsolutePath());
						for (String varLine : Files.readLines(includeFile, Charset.defaultCharset())) {
							varLine = varLine.trim();
							if (varLine.isBlank() || varLine.startsWith("%"))
								continue;
							Preconditions.checkState(varLine.startsWith("\\newcommand{\\dynval"),
									"Unexpected line in dynval include; should start with '\\newcommand{\\dynval': %s",
									varLine);
							try {
								String varName = varLine.substring(varLine.indexOf("dynval"));
								String varValue = varName.substring(varName.indexOf("{")+1);
								varValue = varValue.substring(0, varValue.lastIndexOf("}"));
								varName = varName.substring(0, varName.indexOf("}"));
								System.out.println("\t\t'"+varName+"': '"+varValue+"'");
								dynVals.put(varName, varValue);
							} catch (RuntimeException e) {
								throw new IllegalStateException("Error parsing dynval from: "+varLine, e);
							}
						}
					}
				}
			} else {
				while (line.contains("\\dynval") && !line.startsWith("%")) {
					Preconditions.checkNotNull(dynVals, "dynval command found before include header was found; surround "
							+ "includes with lines stating '%s' and '%s'\n\tLine: %s",
							BEGIN_DYNVAL_INCLUDES_HEADER, END_DYNVAL_INCLUDES_HEADER, line);
					int startIndex = line.indexOf("\\dynval");
					String subString = line.substring(startIndex);
					Preconditions.checkState(subString.contains("{}"), "dynval commands must terminate with {}: %s", line);
					int endIndex = subString.indexOf("{}");
					endIndex += startIndex;
					String varName = line.substring(startIndex+1, endIndex); // +1 skips the \
					Preconditions.checkState(dynVals.containsKey(varName),
							"Variable not found: %s\nfull line: %s", varName, line);
					
					String newLine = "";
					if (startIndex > 0)
						newLine = line.substring(0, startIndex);
					String varValue = dynVals.get(varName);
					String lineStart;
					if (line.length() > 24)
						lineStart = line.substring(0, 20)+" ...";
					else
						lineStart = line;
					
					System.out.println("Processed '"+varName+"' -> '"+varValue+"';\t\t\tLine: "+lineStart);
					newLine += varValue;
					newLine += line.substring(endIndex+2); // +2 here skips the {}
					line = newLine;
				}
				ret.add(line);
			}
		}
		
		return ret;
	}
	
	public static List<String> embedIncludes(List<String> inLines, File refDir) throws IOException {
		return embedIncludes(inLines, refDir, false);
	}
	
	public static List<String> embedIncludes(List<String> inLines, File refDir, boolean skipDynval) throws IOException {
		List<String> ret = new ArrayList<>(inLines.size());
		boolean insideDynval = false;
		for (String line : inLines) {
			if (skipDynval) {
				String trim = line.trim();
				if (insideDynval) {
					if (trim.startsWith(END_DYNVAL_INCLUDES_HEADER))
						insideDynval = false;
					ret.add(line);
					continue;
				} else if (trim.startsWith(BEGIN_DYNVAL_INCLUDES_HEADER)) {
					System.out.println("Won't embed includes from the dynval section");
					insideDynval = true;
					ret.add(line);
					continue;
				}
			}
			if (line.trim().startsWith("\\include{")) {
				line = line.trim();
				String path = line.substring(line.indexOf('{')+1);
				Preconditions.checkState(path.contains("}"), "Line doesn't close '}': %s", line);
				path = path.substring(0, path.indexOf('}'));
				File includeFile = new File(refDir, path);
				if (!includeFile.exists() && !path.toLowerCase().endsWith(".tex"))
					includeFile = new File(refDir, path+".tex");
				System.out.println("\tEmbedding "+includeFile.getName());
				Preconditions.checkState(includeFile.exists(), "Include path not found. raw='%s', abs='%s'",
						path, includeFile.getAbsolutePath());
				ret.addAll(Files.readLines(includeFile, Charset.defaultCharset()));
			} else {
				ret.add(line);
			}
		}
		
		return ret;
	}
	
	private static List<String> stripCommentBlocks(List<String> lines) {
		List<String> ret = new ArrayList<>(lines.size());
		
		boolean insideComment = false;
		
		for (String line : lines) {
			String trimmed = line.trim();
			if (insideComment) {
				if (trimmed.startsWith("\\end{comment}"))
					insideComment = false;
			} else if (trimmed.startsWith("\\begin{comment}")) {
				insideComment = true;
			} else {
				ret.add(line);
			}
		}
		
		return ret;
	}
	
	public static final String BEGIN_EMBEDDED_COMMANDS_HEADER = "% begin embedded commands";
	public static final String END_EMBEDDED_COMMANDS_HEADER = "% end embedded commands";
	
	private static Entry<String, String> processCommand(String line) {
		line = line.trim();
		Preconditions.checkState(line.startsWith("\\newcommand{"));

		// Must be single-line and end cleanly
		if (!line.endsWith("}")) {
			System.err.println("WARNING: skipping multi-line or unterminated \\newcommand: " + line);
			return null;
		}

		// Strip leading \newcommand{
		int pos = "\\newcommand{".length();
		int nameEnd = line.indexOf('}', pos);
		if (nameEnd < 0) {
			System.err.println("WARNING: malformed \\newcommand (missing name '}'): " + line);
			return null;
		}

		String name = line.substring(pos, nameEnd).trim();
		if (name.isEmpty()) {
			System.err.println("WARNING: empty command name in \\newcommand: " + line);
			return null;
		}

		// Normalize to include leading backslash
		if (name.charAt(0) != '\\')
			name = "\\" + name;

		pos = nameEnd + 1;

		// Check for optional [N] argument count (we skip those)
		if (pos < line.length() && line.charAt(pos) == '[') {
			int argEnd = line.indexOf(']', pos + 1);
			if (argEnd < 0) {
				System.err.println("WARNING: malformed \\newcommand argument spec: " + line);
				return null;
			}
			System.err.println("WARNING: skipping \\newcommand with arguments: " + name + " in line: " + line);
			return null;
		}

		// Next must be {replacement}
		if (pos >= line.length() || line.charAt(pos) != '{') {
			System.err.println("WARNING: malformed \\newcommand (missing replacement '{'): " + line);
			return null;
		}

		// Parse replacement, allowing nested braces but staying on one line
		int depth = 0;
		int replStart = pos + 1;
		int replEnd = -1;

		for (int i = pos; i < line.length(); i++) {
			char c = line.charAt(i);
			if (c == '{') {
				depth++;
			} else if (c == '}') {
				depth--;
				if (depth == 0) {
					replEnd = i;
					break;
				}
			}
		}

		if (replEnd < 0) {
			System.err.println("WARNING: unterminated replacement in \\newcommand: " + line);
			return null;
		}

		// Make sure there's nothing extra after the closing brace
		for (int i = replEnd + 1; i < line.length(); i++) {
			if (!Character.isWhitespace(line.charAt(i))) {
				System.err.println("WARNING: extra trailing content in \\newcommand, skipping: " + line);
				return null;
			}
		}

		String replacement = line.substring(replStart, replEnd);
		return new AbstractMap.SimpleEntry<>(name, replacement);
	}
	
	public static final List<String> embedCommands(List<String> inLines, boolean betweenEmbedHeader) {
		List<String> ret = new ArrayList<>(inLines.size());
		
		boolean insideEmbedHeader = false;
		Map<String, String> commands = new LinkedHashMap<>();
		
		int replacements = 0;
		Map<String, Integer> replaceCounts = new HashMap<>(commands.size());
		
		for (int l=0; l<inLines.size(); l++) {
			String line = inLines.get(l);
			String trim = line.trim();
			if (betweenEmbedHeader) {
				if (trim.equalsIgnoreCase(BEGIN_EMBEDDED_COMMANDS_HEADER)) {
					System.out.println("Found '"+BEGIN_EMBEDDED_COMMANDS_HEADER+"'");
					insideEmbedHeader = true;
					continue;
				}
			} else if (trim.startsWith("\\newcommand{")) {
				Entry<String, String> command = processCommand(line);
				if (command != null) {
					System.out.println("Will process command: "+command);
					commands.put(command.getKey(), command.getValue());
				} else {
					ret.add(line);
				}
			}
			if (insideEmbedHeader) {
				if (!trim.isBlank()) {
					if (trim.equalsIgnoreCase(END_EMBEDDED_COMMANDS_HEADER)) {
						System.out.println("Found '"+END_EMBEDDED_COMMANDS_HEADER+"'");
						System.out.println("Processed "+commands.size()+" total embedded commands");
						insideEmbedHeader = false;
						continue;
					} else if (line.startsWith("%")) {
						continue;
					} else {
						if (trim.startsWith("\\newcommand{")) {
							Entry<String, String> command = processCommand(line);
							if (command != null) {
								System.out.println("Will process command: "+command);
								commands.put(command.getKey(), command.getValue());
							}
						}
					}
				}
			} else {
				if (!commands.isEmpty() && line.contains("\\")) {
					// we have a possible command
					
					for (String cmd : commands.keySet()) {
						String repl = commands.get(cmd);

						int from = 0;
						while (true) {
							int idx = line.indexOf(cmd, from);
							if (idx < 0)
								break;

							int after = idx + cmd.length();

							// Require a word-boundary-ish delimiter after the command name, unless it's followed by '{}'
							// This avoids replacing "\foo" inside "\foobar".
							boolean hasEmptyBraces = after + 1 < line.length()
									&& line.charAt(after) == '{'
									&& line.charAt(after + 1) == '}';

							boolean okBoundary = hasEmptyBraces;
							if (!okBoundary) {
								if (after >= line.length()) {
									okBoundary = true;
								} else {
									char c = line.charAt(after);
									okBoundary = !(Character.isLetterOrDigit(c) || c == '_');
								}
							}

							if (!okBoundary) {
								from = after;
								continue;
							}

							// Perform replacement; if "{}" was present, consume it too
							int consumeEnd = hasEmptyBraces ? (after + 2) : after;
							line = line.substring(0, idx) + repl + line.substring(consumeEnd);
							
							replacements++;
							int prevCount = replaceCounts.containsKey(cmd) ? replaceCounts.get(cmd) : 0;
							replaceCounts.put(cmd, prevCount+1);

							// Continue searching after the inserted replacement
							from = idx + repl.length();
						}
					}
				}
				ret.add(line);
			}
		}
		
		System.out.println("Embedded "+replacements+" command replacements");
		for (String command : commands.keySet()) {
			int count = replaceCounts.containsKey(command) ? replaceCounts.get(command) : 0;
			System.out.println("\t"+command+":\t"+count);
		}
		
		return ret;
	}
	
	public static final List<String> embedBibliography(List<String> inLines, File bblFile) throws IOException {
		List<String> bblLines = Files.readLines(bblFile, Charset.defaultCharset());
		
		List<String> ret = new ArrayList<>(inLines.size() + bblLines.size());
		
		boolean found = false;
		for (String line : inLines) {
			if (line.trim().startsWith("\\bibliography")) {
				System.out.println("Found bibliography, replacing: "+line);
				Preconditions.checkState(!found, "Encounted \\bibliography twice");
				found = true;
				ret.addAll(bblLines);
			} else {
				ret.add(line);
			}
		}
		
		Preconditions.checkState(found, "No lines starting with \\bibliography found");
		
		return ret;
	}
	
	public static final String BEGIN_FIGPANEL_HEADER = "% begin figpanel";
	public static final String END_FIGPANEL_HEADER = "% end figpanel";
	
	public static final List<String> embedFigpanelFigures(List<String> inLines, File refDir, File figpanelDir, String textWidth) throws IOException {
		Preconditions.checkState(figpanelDir.exists() || figpanelDir.mkdirs(),
				"Figpanel dir doesn't exist and couldn't be created: %s", figpanelDir.getAbsolutePath());
		List<String> ret = new ArrayList<>(inLines.size());
		
		boolean readingFigure = false;
		boolean readingFigpanel = false;
		List<String> figpanelLines = null;
		String currentWidth = null;
		
		for (int l=0; l<inLines.size(); l++) {
			String line = inLines.get(l);
			String trimmed = line.trim();
			
			if (trimmed.startsWith("\\input{figpanel"))
				continue; // don't include this line in output
			
			if (trimmed.startsWith("\\begin{figure}")) {
				Preconditions.checkState(!readingFigure, "Already reading a figure but encountered: %s", line);
				Preconditions.checkState(!readingFigpanel, "Already reading figpanel, encountered: %s", line);
				readingFigure = true;
			} else if (trimmed.startsWith("\\end{figure}")) {
				Preconditions.checkState(readingFigure, "Not reading a figure but encountered: %s", line);
				Preconditions.checkState(!readingFigpanel, "Ending a figure but still reading figpanel, encountered: %s", line);
				readingFigure = false;
			} else if (trimmed.startsWith(BEGIN_FIGPANEL_HEADER)) {
				Preconditions.checkState(readingFigure, "Not reading a figure but encountered: %s", line);
				Preconditions.checkState(!readingFigpanel, "Already reading figpanel, encountered: %s", line);
				Preconditions.checkState(figpanelLines == null, "Figpanel lines not reset?");
				readingFigpanel = true;
				figpanelLines = new ArrayList<>();
				currentWidth = null;
				if (trimmed.contains("width=")) {
					currentWidth = trimmed.substring(trimmed.indexOf("width=")).trim();
				}
				continue; // don't include this line in output
			} else if (trimmed.startsWith(END_FIGPANEL_HEADER)) {
				readingFigpanel = false;
				Preconditions.checkState(!figpanelLines.isEmpty(), "No figpanel lines encountered");
				
				String label = findCurrentFigureLabel(inLines, l+1);
				System.out.println("Building figpanel figure: "+label);
				if (currentWidth == null) {
					currentWidth = "width=\\textwidth";
					System.out.println("\tno width specified, using full "+currentWidth);
				} else {
					System.out.println("\tusing custom "+currentWidth);
				}
				
				String figPath = buildFigpanelFigure(refDir, figpanelDir, figpanelLines, label, textWidth);
				line = "    \\includegraphics["+currentWidth+"]{"+figPath+"}";
				
				figpanelLines = null;
			} else if (readingFigpanel) {
				figpanelLines.add(line);
				continue; // don't include this line in output
			}
			ret.add(line);
		}
		
		return ret;
	}
	
	private static String findCurrentFigureLabel(List<String> lines, int startIndex) {
		String currentLabel = null;
		for (int i=startIndex; i<lines.size(); i++) {
			String line = lines.get(i);
			String trimmed = line.trim();
			
			if (line.contains("\\label{")) {
				String label = line.substring(line.indexOf("\\label{")+7);
				Preconditions.checkState(label.contains("}"), "Parsing label, not closed with '}': %s", label);
				label = label.substring(0, label.indexOf("}"));
				if (label.contains(":"))
					label = label.substring(label.lastIndexOf(':')+1);
				label = label.trim();
				Preconditions.checkState(!label.isBlank(), "Blank label encountered: %s", label);
				Preconditions.checkState(currentLabel == null,
						"Duplicate label encountered while reading figpanel figure; prev='%s', new='%s'", currentLabel, label);
				currentLabel = label;
			}
			
			if (trimmed.startsWith("\\end{figure}")) {
				Preconditions.checkState(currentLabel != null, "Figure terminated with no label detected: %s", line);
				break;
			}
		}
		Preconditions.checkState(currentLabel != null, "Figure label not found; checked %s lines", (lines.size()-startIndex));
		return currentLabel;
	}
	
	private static String buildFigpanelFigure(File refDir, File figpanelDir, List<String> figpanelLines, String label,
			String textWidth) throws IOException {
		Path refPath = refDir.toPath().toAbsolutePath().normalize();
		Path figpanelPath = figpanelDir.toPath().toAbsolutePath().normalize();

		if (!figpanelPath.startsWith(refPath))
			throw new IllegalArgumentException("figpanelDir is not under refDir");
		
		File texFile = new File(figpanelDir, label+".tex");
		
		// write the figpanel tex
		FileWriter fw = new FileWriter(texFile);
		
		// tiny padding on top, nowhere else
		fw.write("\\documentclass[border={0pt 0pt 0pt 2pt}]{standalone}\n");
//		fw.write("\\documentclass[border=0pt]{standalone}\n");
		fw.write("\n");
		if (textWidth != null && !textWidth.isBlank()) {
			fw.write("\\setlength{\\textwidth}{"+textWidth+"}\n");
			fw.write("\n");
		}
		fw.write("\\input{figpanel}\n");
		fw.write("\n");
		fw.write("\\begin{document}\n");
		fw.write("\n");
		for (String line : figpanelLines) {
			fw.write(line);
			fw.write("\n");
		}
		fw.write("\n");
		fw.write("\\end{document}\n");
		
		fw.close();
		
		return buildLatexPDF(refDir, texFile, figpanelDir, label);
	}
	
	/**
	 * Builds a PDF for the given LaTeX file using latexmk (must be installed natively and on the path)
	 * 
	 * @param refDir root of the latex directory, latexmk will be run from here
	 * @param texFile input .tex file
	 * @param outputDir output directory where the PDF will be written
	 * @param prefix PDF name prefix; output will be outputDir/prefix.pdf
	 * @return relative path from refDir to the written PDF file
	 * @throws IOException
	 */
	public static String buildLatexPDF(File refDir, File texFile, File outputDir, String prefix)
			throws IOException {

		Path ref = refDir.toPath().toAbsolutePath().normalize();
		Path tex = texFile.toPath().toAbsolutePath().normalize();
		Path out = outputDir.toPath().toAbsolutePath().normalize();

		if (!java.nio.file.Files.isDirectory(ref))
			throw new IllegalArgumentException("refDir is not a directory: " + refDir);

		if (!java.nio.file.Files.isRegularFile(tex))
			throw new IllegalArgumentException("texFile is not a file: " + texFile);

		if (!java.nio.file.Files.isDirectory(out))
			throw new IllegalArgumentException("outputDir is not a directory: " + outputDir);

		if (!tex.startsWith(ref))
			throw new IllegalArgumentException("texFile is not under refDir: " + texFile);

		if (!out.startsWith(ref))
			throw new IllegalArgumentException("outputDir is not under refDir: " + outputDir);

		// latexmk will create prefix.pdf in outputDir (relative to refDir)
		String relTex = ref.relativize(tex).toString().replace(File.separatorChar, '/');
		String relOutDir = ref.relativize(out).toString().replace(File.separatorChar, '/');

		// Isolate auxiliary files (and latexmk temp) into a Java temp directory
		Path auxDir = java.nio.file.Files.createTempDirectory("latexmk-aux-");
		// Ensure cleanup even if the build fails (best-effort)
		try {
			String auxDirArg = auxDir.toAbsolutePath().normalize().toString();

			ProcessBuilder pb = new ProcessBuilder(List.of(
					"latexmk",
					"-pdf",
					"-interaction=nonstopmode",
					"-halt-on-error",
					"-outdir=" + relOutDir,
					"-auxdir=" + auxDirArg,
					"-jobname=" + prefix,
					relTex
			));
			pb.directory(refDir);
			pb.redirectErrorStream(true);

			Process p = pb.start();

			ByteArrayOutputStream buf = new ByteArrayOutputStream();
			try (InputStream in = p.getInputStream()) {
				in.transferTo(buf);
			}

			int exit;
			try {
				exit = p.waitFor();
			} catch (InterruptedException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
			if (exit != 0) {
				String log = buf.toString(StandardCharsets.UTF_8);
				throw new IOException("latexmk failed (exit=" + exit + ")\n" + log);
			}

			Path pdf = out.resolve(prefix + ".pdf").toAbsolutePath().normalize();
			if (!java.nio.file.Files.isRegularFile(pdf))
				throw new IOException("latexmk succeeded but PDF not found: " + pdf);

			// Return path from refDir to the written PDF, using '/' for LaTeX friendliness
			return ref.relativize(pdf).toString().replace(File.separatorChar, '/');
		} finally {
			FileUtils.deleteRecursive(auxDir.toFile());
		}
	}
	public static List<String> readTex(File inputFile) throws IOException {
		return Files.readLines(inputFile, Charset.defaultCharset());
	}
	
	public static void writeTex(List<String> lines, File outFile) throws IOException {
		File tmpOutputFile = new File(outFile.getParentFile(), outFile.getName()+".tmp");
		FileWriter fw = new FileWriter(tmpOutputFile, Charset.defaultCharset());
		
		for (String line : lines) {
			fw.write(line);
			fw.write("\n");
		}
		
		fw.close();
		Files.move(tmpOutputFile, outFile);
	}
	
	public static void reorganizeFiguresForSSA(File inputTexFile, File outputDir) throws IOException {
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir(),
				"Output dir doesn't exist and could not be created: %s", outputDir.getAbsolutePath());
		
		File inputDir = inputTexFile.getParentFile();
		File outputTexFile = new File(outputDir, inputTexFile.getName());
		
		List<String> lines = Files.readLines(inputTexFile, Charset.defaultCharset());
		lines = stripCommentBlocks(lines);
		
		File tmpOutputFile = new File(outputTexFile.getParentFile(), outputTexFile.getName()+".bak");
		FileWriter fw = new FileWriter(tmpOutputFile, Charset.defaultCharset());
		
		String[] extensions = { "pdf", "eps", "ps", "png", "jpg", "jpeg", "gif" };
		
		// figure out total figure count
		int totalNumFigures = 0;
		for (String line : lines)
			if (line.trim().startsWith("\\begin{figure}"))
				totalNumFigures++;
		System.out.println("Identified "+totalNumFigures+" figures");
		
		DecimalFormat numDF = new DecimalFormat((totalNumFigures+"").replaceAll(".", "0"));
		
		boolean readingFigure = false;
		boolean readingSubfigure = false;
		boolean subfigHasCaption = false;
		int curFigNum = 1;
		char curSubfigChar = 'a';
		int curSubfigNumExtra = 0;
		
		for (String line : lines) {
			String trimmed = line.trim();
			
			if (readingFigure) {
				if (trimmed.startsWith("\\end{figure}")) {
					System.out.println("DONE with figure "+curFigNum);
					Preconditions.checkState(!readingSubfigure, "done reading figure but still reading subfigure?");
					readingFigure = false;
					curFigNum++;
					curSubfigChar = 'a';
				} else if (trimmed.startsWith("\\begin{subfigure}")) {
					Preconditions.checkState(!readingSubfigure, "already reading a subfigure?");
					readingSubfigure = true;
					subfigHasCaption = false;
				} else if (readingSubfigure && trimmed.startsWith("\\caption")) {
					subfigHasCaption = true;
				} else if (trimmed.startsWith("\\end{subfigure}")) {
					Preconditions.checkState(readingSubfigure, "ending subfigure, but not reading a subfigure?");
					readingSubfigure = false;
				} else if (trimmed.startsWith("\\includegraphics")) {
					int pathStartIndex = line.indexOf('{')+1;
					Preconditions.checkState(pathStartIndex > 0,
							"on an includegraphics line but path {} not found; must be on a single line: %s", trimmed);
					String beforePath = line.substring(0, pathStartIndex);
					String path = line.substring(pathStartIndex);
					int pathEndIndex = path.indexOf('}');
					Preconditions.checkState(pathEndIndex >= 0, "Line doesn't close '}': %s", trimmed);
					String afterPath = path.substring(pathEndIndex);
					path = path.substring(0, pathEndIndex);
					
					File figInFile = new File(inputDir, path);
					for (int i=0; !figInFile.exists() && i<extensions.length; i++) {
						figInFile = new File(inputDir, path+"."+extensions[i]);
					}
					Preconditions.checkState(figInFile.exists(),
							"Figure input file not found: %s. Also tried these extensions: %s", path, extensions);
					String outName = "figure_"+numDF.format(curFigNum);
					if (readingSubfigure) {
						if (subfigHasCaption) {
							outName += "_"+curSubfigChar;
							Preconditions.checkState(curSubfigChar != 'z',
									"Ran out of subfigure letters for figure %s", curFigNum);
							curSubfigChar++;
						} else if (path.contains("cpt")) {
							outName += "_colorscale";
						} else {
							outName += "_extra";
							curSubfigNumExtra++;
							if (curSubfigNumExtra > 1)
								outName += "_"+curSubfigNumExtra;
						}
					}
					String inName = figInFile.getName();
					String inPrefix = inName;
					String outPrefix = outName;
					if (inPrefix.contains(".")) {
						inPrefix = inPrefix.substring(0, inPrefix.lastIndexOf("."));
						String ext = inName.substring(inName.lastIndexOf(".")+1);
						outName += "."+ext;
					}
					
					File figOutFile = new File(outputDir, outName);
					System.out.println("\t"+path+"\t->\t"+outName);
					
					if (outName.endsWith(".pdf"))
						PdfFirstPageCopier.copyFirstPage(figInFile, figOutFile);
					else
						Files.copy(figInFile, figOutFile);
					if (outName.endsWith(".png")) {
						// see if a tif is available because SSA is weird
						File tifFile = new File(figInFile.getParentFile(), inPrefix+".tif");
						if (!tifFile.exists())
							tifFile = new File(figInFile.getParentFile(), inPrefix+".tiff");
						System.out.flush();
						if (tifFile.exists()) {
							System.out.println("\tAlso copying TIF");
							System.err.println("\tWARNING: make sure you converted it so that the metadata thinks it's "
									+ "300 DPI, otherwise they might get confused and send it back: convert figure.png -density 300 figure.tif");
							Files.copy(tifFile, new File(outputDir, outPrefix+".tif"));
						} else {
							System.err.println("\tWARNING: SSA arbitrarily discriminates agains PNGs and won't accept them! "
									+ "Submit TIFs also, and convert with density set to 300: convert figure.png -density 300 figure.tif");
						}
						System.err.flush();
					}
					
					String newLine = beforePath+outName+afterPath;
					fw.write(newLine);
					fw.write("\n");
					continue;
				}
			} else if (trimmed.startsWith("\\begin{figure}")) {
				readingFigure = true;
			}
			fw.write(line);
			fw.write("\n");
		}
		
		fw.close();
		Files.move(tmpOutputFile, outputTexFile);
	}

	public static void main(String[] args) throws IOException {
//		String raw = "Special characters: & % $ # _ { } ~ ^ \\ and \\& already escaped.";
//		String escaped = escapeLaTeX(raw);
//		System.out.println("Original: " + raw);
//		System.out.println("Escaped: " + escaped);
//		
//		System.out.println(numberAsPercent(100d*0.505/1.3, 0));
//		
//		System.out.println(defineValueCommand("SlipRateExample", "33%"));
//		System.out.println(defineValueCommand("SlipRateExample", "\\expnum("+numberExpFormat(3.14)+")", false));
//		System.out.println(defineValueCommand("SlipRateExample", "\\expnum("+numberExpFormat(3e-10)+")", false));
		
		File prviDir = new File("/home/kevin/Documents/papers/2025_PRVI_ERF");
		
		File mainBranch = new File(prviDir, "prvi25-erf-paper");
		
//		File initialBranch = new File(prviDir, "initial-bssa-submission");
//		File refDir = initialBranch;
//		File inputFile = new File(mainBranch, "submission/original/main_for_diff.tex");
//		File outputFile = new File(mainBranch, "submission/original/embedded_for_diff.tex");
		
		String textWidth = "517.79993pt";
		File refDir = mainBranch;
		File figpanelDir = new File(refDir, "Figures/figpanel");
		File inputFile = new File(mainBranch, "main.tex");
		File outputFile = new File(mainBranch, "main_embedded.tex");
		
		List<String> lines = readTex(inputFile);

		lines = embedIncludes(lines, refDir, true); // embed includes first, but skip dynval ones
		lines = embedDynvalIncludes(lines, refDir);
		lines = embedCommands(lines, true);
		lines = embedBibliography(lines, new File(refDir, "refs_compiled.bbl"));
		lines = embedFigpanelFigures(lines, refDir, figpanelDir, textWidth);
		
		writeTex(lines, outputFile);
		
		File ssaOutputDir = new File(new File(mainBranch, "submission"), "final-bssa-submission");
//		File ssaOutputDir = new File(prviDir, "final-bssa-submission-tests");
		Preconditions.checkState(ssaOutputDir.exists() || ssaOutputDir.mkdir(),
				"SSA output dir doesn't exist and could not be created: %s", ssaOutputDir.getAbsolutePath());
		reorganizeFiguresForSSA(outputFile, ssaOutputDir);
		// rename from embedded to just main
		Files.move(new File(ssaOutputDir, outputFile.getName()), new File(ssaOutputDir, inputFile.getName()));
	}

}
