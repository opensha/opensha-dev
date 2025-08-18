	package scratch.kevin.latex;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.util.DataUtils;

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
	
	public static void embedDynvalIncludes(File inputTexFile, File outputTexFile) throws IOException {
		embedDynvalIncludes(inputTexFile, inputTexFile.getParentFile(), outputTexFile);
	}
	
	public static void embedDynvalIncludes(File inputTexFile, File refDir, File outputTexFile) throws IOException {
		boolean processingIncludes = false;
		Map<String, String> dynVals = null;
		
		List<String> lines = Files.readLines(inputTexFile, Charset.defaultCharset());
		File tmpOutputFile = new File(outputTexFile.getParentFile(), outputTexFile.getName()+".bak");
		FileWriter fw = new FileWriter(tmpOutputFile, Charset.defaultCharset());
		for (String line : lines) {
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
					System.out.println("Processed '"+varName+"' -> '"+varValue+"'");
					newLine += varValue;
					newLine += line.substring(endIndex+2); // +2 here skips the {}
					line = newLine;
				}
				fw.write(line);
				fw.write("\n");
			}
		}
		
		fw.close();
		Files.move(tmpOutputFile, outputTexFile);
	}
	
	public static void embedIncludes(File inputTexFile, File outputTexFile) throws IOException {
		embedIncludes(inputTexFile, inputTexFile.getParentFile(), outputTexFile);
	}
	
	public static void embedIncludes(File inputTexFile, File refDir, File outputTexFile) throws IOException {
		List<String> lines = Files.readLines(inputTexFile, Charset.defaultCharset());
		File tmpOutputFile = new File(outputTexFile.getParentFile(), outputTexFile.getName()+".bak");
		FileWriter fw = new FileWriter(tmpOutputFile, Charset.defaultCharset());
		for (String line : lines) {
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
				for (String includeLine : Files.readLines(includeFile, Charset.defaultCharset())) {
					fw.write(includeLine);
					fw.write("\n");
				}
			} else {
				fw.write(line);
				fw.write("\n");
			}
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
		
		File prviDir = new File("/home/kevin/Documents/papers/2024_PRVI_ERF");
		
		File mainBranch = new File(prviDir, "prvi25-erf-paper");
		File initialBranch = new File(prviDir, "initial-bssa-submission");
		
//		File refDir = initialBranch;
//		File inputFile = new File(mainBranch, "submission/original/main_for_diff.tex");
//		File outputFile = new File(mainBranch, "submission/original/embedded_for_diff.tex");
		
		File refDir = mainBranch;
		File inputFile = new File(mainBranch, "main.tex");
		File outputFile = new File(mainBranch, "main_embedded.tex");
		
		embedDynvalIncludes(inputFile, refDir, outputFile);
		embedIncludes(outputFile, refDir, outputFile);
	}

}
