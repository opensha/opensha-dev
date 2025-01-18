	package scratch.kevin.latex;

import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Map;

import org.opensha.commons.util.DataUtils;

import com.google.common.base.Preconditions;

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

	public static void main(String[] args) {
		String raw = "Special characters: & % $ # _ { } ~ ^ \\ and \\& already escaped.";
		String escaped = escapeLaTeX(raw);
		System.out.println("Original: " + raw);
		System.out.println("Escaped: " + escaped);
		
		System.out.println(numberAsPercent(100d*0.505/1.3, 0));
		
		System.out.println(defineValueCommand("SlipRateExample", "33%"));
		System.out.println(defineValueCommand("SlipRateExample", "\\expnum("+numberExpFormat(3.14)+")", false));
		System.out.println(defineValueCommand("SlipRateExample", "\\expnum("+numberExpFormat(3e-10)+")", false));
	}

}
