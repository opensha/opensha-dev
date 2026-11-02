package scratch.anne.sra.io.readers;

import java.io.IOException;
import java.nio.file.*;
import java.util.*;
import java.util.stream.*;

import com.google.gson.*;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;

import scratch.anne.sra.domain.structural_response.ResponseModelLibrary;
import scratch.anne.sra.domain.structural_response.fragilities.FragilityModel;
import scratch.anne.sra.domain.structural_response.fragilities.definitions.*;
import scratch.anne.sra.util.InputParsingUtil;
import scratch.anne.sra.util.Metadata;
import scratch.anne.sra.util.StringUtil;
import scratch.anne.sra.util.enums.IMT;
import scratch.anne.sra.util.enums.LimitState;

/**
 * Reader for FragilityModel libraries.
 * <p>
 * Supports:
 * <ul>
 *     <li>JSON input (component -> limit states, mixed lognormal/discrete)</li>
 *     <li>CSV input (long-format, mixed fragility types)</li>
 * </ul>
 * Provides convenience constructors that automatically generate metadata.
 */
public class FragilityLibraryReader {

    // ------------------------------------------------------------------------
    // Convenience reader (auto metadata)
    // ------------------------------------------------------------------------
    public static ResponseModelLibrary<FragilityModel> readLibrary(Path path) throws IOException {
        Metadata metadata = new Metadata.Builder()
                .set("librarySource", path.toString())
                .set("description", "Fragility library loaded from file")
                .set("creationInfo", "Generated " + java.time.Instant.now())
                .build();

        return readLibrary(path, metadata);
    }

    // ------------------------------------------------------------------------
    // Convenience reader (explicit metadata fields)
    // ------------------------------------------------------------------------
    public static ResponseModelLibrary<FragilityModel> readLibrary(
            Path path,
            String librarySource,
            String description,
            String creationInfo) throws IOException {

        Metadata metadata = new Metadata.Builder()
                .set("librarySource", librarySource)
                .set("description", description)
                .set("creationInfo", creationInfo)
                .build();

        return readLibrary(path, metadata);
    }

    // ------------------------------------------------------------------------
    // Core reader (single source of truth)
    // ------------------------------------------------------------------------
    public static ResponseModelLibrary<FragilityModel> readLibrary(
            Path path,
            Metadata metadata) throws IOException {

        String filename = path.getFileName().toString().toLowerCase();
        List<FragilityModel> models;

        if (filename.endsWith(".json")) {
            models = readJson(path);
        } else if (filename.endsWith(".csv")) {
            models = readCsv(path);
        } else {
            throw new IllegalArgumentException("Unsupported file format: " + filename);
        }

        return ResponseModelLibrary.of(models, metadata);
    }

    // ------------------------------------------------------------------------
    // JSON parser
    // ------------------------------------------------------------------------
    private static List<FragilityModel> readJson(Path path) throws IOException {
        Gson gson = new Gson();
        String jsonContent = Files.readString(path);

        JsonObject root = gson.fromJson(jsonContent, JsonObject.class);
        List<FragilityModel> allModels = new ArrayList<>();

        for (Map.Entry<String, JsonElement> compEntry : root.entrySet()) {
            String compID = compEntry.getKey();
            JsonObject compNode = compEntry.getValue().getAsJsonObject();

            allModels.add(loadComponentFromJson(compID, compNode));
        }

        return allModels;
    }

    private static FragilityModel loadComponentFromJson(String compID, JsonObject compNode) {
    	//TODO add description and comments option to fragiliy model
        String description = compNode.has("description") ? compNode.get("description").getAsString() : "";
        String comments = compNode.has("comments") ? compNode.get("comments").getAsString() : "";
     
        LimitState primaryLS = null;
        if (compNode.has("primaryLS") && !compNode.get("primaryLS").isJsonNull()) {
            String lsStr = compNode.get("primaryLS").getAsString().trim();
            try {
                primaryLS = LimitState.valueOf(lsStr.toUpperCase());
            } catch (IllegalArgumentException e) {
                throw new IllegalArgumentException("Invalid primary LimitState for component " + compID + ": " + lsStr);
            }
        }
        
        String imtStr = compNode.has("imt") ? compNode.get("imt").getAsString() : null;
        IMT imt;
        Double period = null;
        try {
        	StringUtil.ImtPeriod imtPeriod = StringUtil.stringToImtPeriod(imtStr);
            imt = imtPeriod.imt;
            period = imtPeriod.period;        	
        } catch (IllegalArgumentException e) {
            throw new IllegalArgumentException("Invalid primary LimitState for component " + compID + ": " + imtStr);
        }
        
        JsonObject limitStatesNode = compNode.getAsJsonObject("limitStates");
        List<LimitStateFragility> lsList = new ArrayList<>();

        for (Map.Entry<String, JsonElement> lsEntry : limitStatesNode.entrySet()) {
            String lsName = lsEntry.getKey();
            JsonObject lsNode = lsEntry.getValue().getAsJsonObject();
            
            // Convert the LS string to the enum         
            LimitState ls;
            try {
                ls = LimitState.valueOf(lsName.toUpperCase());
            } catch (IllegalArgumentException e) {
            	throw new IllegalArgumentException("Invalid LimitState for component " + compID + ": " + lsName);
            }

            FragilityDefinition def;
            String type = lsNode.get("type").getAsString().toLowerCase();

            switch (type) {
                case "lognormal":
                    double median = lsNode.get("median").getAsDouble();
                    double beta = lsNode.get("beta").getAsDouble();
                    def = new LognormalFragilityDefinition(median, beta);
                    break;

                case "discrete":
                    def = parseDiscrete(lsNode, compID, lsName);
                    break;

                default:
                    throw new IllegalArgumentException("Unknown fragility type: " + type);
            }

            lsList.add(new LimitStateFragility(ls, def));
        }

        return new FragilityModel(compID, imt, period, lsList, primaryLS);
    }

    private static DiscreteFragilityDefinition parseDiscrete(
            JsonObject dsNode,
            String compID,
            String lsName) {

        Gson gson = new Gson();

        double[] imLevels =
                InputParsingUtil.parseJsonNumericArray(
                        gson.fromJson(
                                dsNode.get("imLevels"),
                                Object.class));

        double[] probability =
                InputParsingUtil.parseJsonNumericArray(
                        gson.fromJson(
                                dsNode.get("probability"),
                                Object.class));

        InputParsingUtil.checkMonotonicIM(
                compID + ":" + lsName,
                imLevels);

        return new DiscreteFragilityDefinition(
                imLevels,
                probability);
    }

    // ------------------------------------------------------------------------
    // CSV parser
    // ------------------------------------------------------------------------
    private static List<FragilityModel> readCsv(Path path) throws IOException {

        List<String[]> rows;
        try (Stream<String> lines = Files.lines(path)) {
            rows = lines.skip(1)
                    .map(line -> Arrays.stream(line.split(",", -1))
                    		.map(InputParsingUtil::cleanCsvString)
                            .toArray(String[]::new))
                    .filter(tokens -> tokens.length >= 9)
                    .collect(Collectors.toList());
        }

        if (rows.isEmpty())
            return Collections.emptyList();

        // ------------------------------------------------------------
        // group by component
        // ------------------------------------------------------------
        Map<String, List<String[]>> byComponent =
                rows.stream().collect(
                        Collectors.groupingBy(r -> r[0],
                                LinkedHashMap::new,
                                Collectors.toList()));

        List<FragilityModel> models = new ArrayList<>();

        for (Map.Entry<String, List<String[]>> compEntry : byComponent.entrySet()) {

            String compID = compEntry.getKey();
            List<String[]> compRows = compEntry.getValue();

            /* --------------------------------------------------------
             * Primary LimitState
             * -------------------------------------------------------- */
            LimitState primaryLS = null;
            
            Set<String> primaryValues = compRows.stream()
                    .map(r -> r[1].trim())
                    .filter(s -> !s.isEmpty())
                    .collect(Collectors.toSet());

            if (primaryValues.size() > 1) {
            	throw new IllegalArgumentException("Multiple primaryLS for component " + compID + ": " + primaryValues);
            }


            if (!primaryValues.isEmpty()) {
                String primaryStr = primaryValues.iterator().next();

                try {
                    primaryLS = LimitState.valueOf(primaryStr.toUpperCase());
                } catch (IllegalArgumentException e) {
                	throw new IllegalArgumentException("Invalid primary LimitState for component "
                            + compID + ": " + primaryStr);
                }
            }

            /* --------------------------------------------------------
             * IMT (single per component)
             * -------------------------------------------------------- */
            
            Set<String> imtValues = compRows.stream()
                    .map(r -> r[2].trim())
                    .filter(s -> !s.isEmpty())
                    .collect(Collectors.toSet());

            
            if (imtValues.size() != 1) {
                throw new IllegalArgumentException(
                        "Component " + compID +
                        " must define exactly one IMT, found: " + imtValues);
            }

            String imtStr = imtValues.iterator().next();

            IMT imt;
            Double period;

            try {
                StringUtil.ImtPeriod imtPeriod =
                        StringUtil.stringToImtPeriod(imtStr);

                imt = imtPeriod.imt;
                period = imtPeriod.period;

            } catch (IllegalArgumentException e) {
                throw new IllegalArgumentException(
                		"Invalid imt for component "
                                + compID + ": " + imtStr);
            }

            /* --------------------------------------------------------
             * group rows by limit state
             * -------------------------------------------------------- */
            Map<String, List<String[]>> byLS =
                    compRows.stream().collect(
                            Collectors.groupingBy(
                                    r -> r[3],
                                    LinkedHashMap::new,
                                    Collectors.toList()));

            List<LimitStateFragility> lsList = new ArrayList<>();

            for (Map.Entry<String, List<String[]>> lsEntry : byLS.entrySet()) {

                String lsName = lsEntry.getKey();
                List<String[]> lsRows = lsEntry.getValue();

                LimitState ls;
                try {
                    ls = LimitState.valueOf(lsName.toUpperCase());
                } catch (IllegalArgumentException e) {
                	throw new IllegalArgumentException("Invalid LimitState for component "
                            + compID + ": " + lsName);
                }

                String type = lsRows.get(0)[4].toLowerCase();
                FragilityDefinition def;

                switch (type) {

                    case "lognormal":
                        double median =
                        		InputParsingUtil.parseDoubleOrThrow(lsRows.get(0)[7]);
                        double beta =
                        		InputParsingUtil.parseDoubleOrThrow(lsRows.get(0)[8]);

                        def = new LognormalFragilityDefinition(
                                median, beta);
                        break;

                    case "discrete":

                        double[] imLevels = lsRows.stream()
                                .mapToDouble(r -> InputParsingUtil.parseDoubleOrThrow(r[5]))
                                .toArray();

                        double[] probabilities = lsRows.stream()
                                .mapToDouble(r -> InputParsingUtil.parseDoubleOrThrow(r[6]))
                                .toArray();

                        def = new DiscreteFragilityDefinition(
                                imLevels, probabilities);
                        break;

                    default:
                        throw new IllegalArgumentException(
                                "Unknown fragility type: " + type);
                }

                lsList.add(new LimitStateFragility(ls, def));
            }

            models.add(
                    new FragilityModel(
                            compID,
                            imt,
                            period,
                            lsList,
                            primaryLS));
        }

        return models;
    }
 
}
