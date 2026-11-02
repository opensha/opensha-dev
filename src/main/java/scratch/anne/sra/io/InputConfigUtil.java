package scratch.anne.sra.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class InputConfigUtil {
	
    // ------------------- Helper: Load CSV config -------------------
    public static Map<String, String> loadConfig(Path csvPath) throws IOException {
        Map<String, String> map = new HashMap<>();
        try (BufferedReader br = Files.newBufferedReader(csvPath)) {
            br.readLine(); // skip header
            String line;
            while ((line = br.readLine()) != null) {
                String[] tokens = line.split(",", 2);
                if (tokens.length == 2) {
                    map.put(tokens[0].trim(), tokens[1].trim());
                }
            }
        }
        return map;
    }
    
    public static List<Double> parseDoubleList(String value) {

        if (value == null || value.isBlank())
            return List.of();

        // strip surrounding quotes FIRST
        value = value.trim();

        if ((value.startsWith("\"") && value.endsWith("\"")) ||
            (value.startsWith("'") && value.endsWith("'"))) {
            value = value.substring(1, value.length() - 1);
        }

        // remove brackets
        value = value.replace("[", "")
                     .replace("]", "");

        return Arrays.stream(value.split("[,;]"))
                .map(String::trim)
                .filter(s -> !s.isEmpty())
                .map(Double::parseDouble)
                .toList();
    }

}
