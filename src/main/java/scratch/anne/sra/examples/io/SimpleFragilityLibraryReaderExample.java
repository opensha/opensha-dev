package scratch.anne.sra.examples.io;

import java.nio.file.*;
import java.util.*;

import scratch.anne.sra.domain.structural_response.ResponseModelLibrary;
import scratch.anne.sra.domain.structural_response.SimpleImResponse;
import scratch.anne.sra.domain.structural_response.SimpleImResponseLibrary;
import scratch.anne.sra.domain.structural_response.fragilities.FragilityModel;
import scratch.anne.sra.domain.structural_response.fragilities.SimpleImFragilityLibraryPreparer;
import scratch.anne.sra.domain.structural_response.fragilities.SimpleImFragilityLibraryPreparer.PrepareSpec;
import scratch.anne.sra.io.readers.FragilityLibraryReader;
import scratch.anne.sra.util.AssetKeys.ImKey;
import scratch.anne.sra.util.AssetKeys.ImKey.ImDomain;
import scratch.anne.sra.util.NumericUtil.ImArrayParam;

public class SimpleFragilityLibraryReaderExample {

    public static void main(String[] args) {

        try {

            // =========================================================
            // STEP 1: FILE PATH SETUP (JSON or CSV)
            // =========================================================

            Path resourceFolder = Path.of(
                    "C:/Users/ahulsey/git/opensha-dev/src/main/resources/scratch/anne/sra/mini"
            );

            Path fragLibraryJSON = resourceFolder.resolve("fragilities.json");
            Path fragLibraryCSV  = resourceFolder.resolve("fragilities.csv");

            // Choose input format
            Path filePath = fragLibraryJSON; // or fragLibraryCSV

            // =========================================================
            // STEP 2: READ FRAGILITY MODEL LIBRARY
            // =========================================================

            ResponseModelLibrary<FragilityModel> modelLibrary =
                    FragilityLibraryReader.readLibrary(
                            filePath,
                            "fragility library source",
                            "fragility description",
                            "date, author, workflow v1.0"
                    );

            System.out.println("\n=== RAW FRAGILITY LIBRARY ===");
            System.out.println(modelLibrary);

            // =========================================================
            // STEP 3: DEFINE PORTFOLIO FILTER
            // =========================================================

            Set<String> portfolioFragilities = new HashSet<>(Arrays.asList(
                    "Fragility_A",
                    "Fragility_C"
            ));

            // =========================================================
            // STEP 4: BUILD PREPARATION SPEC
            // =========================================================

            PrepareSpec defaultSpec = PrepareSpec.builder()
                    .build();
            
            PrepareSpec generatedSpec = PrepareSpec.builder()
                    .filterNames(portfolioFragilities)
                    .selectedLS("Collapse")
                    .imParam(
                            ImArrayParam.builder()
                                    .logStep(0.05)
                                    .spread(5.0)
                                    .build()
                    )
                    .build();
            
            PrepareSpec customImSpec = PrepareSpec.builder()
                    .filterNames(portfolioFragilities)
                    .selectedLS("Collapse")
                    .customImValues(new double[] {
                    		0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1
                    })
                    .build();
            
            // select one of the above specs
            PrepareSpec spec = customImSpec;

            // =========================================================
            // STEP 5: PREPARE IM-RESPONSE FRAGILITY LIBRARY
            // =========================================================

            SimpleImResponseLibrary fragilityLib =
                    SimpleImFragilityLibraryPreparer.prepare(modelLibrary, spec);

            System.out.println("\n=== PREPARED FRAGILITY LIBRARY ===");
            System.out.println(fragilityLib.summary());

            // =========================================================
            // STEP 6: PRINT ALL MODELS + RESPONSE TABLES
            // =========================================================

            System.out.println("\n=== ALL FRAGILITIES ===");

            for (SimpleImResponse frag : fragilityLib.getModels()) {

                System.out.println("\n-------------------------");
                System.out.println("Name: " + frag.getName());
                System.out.println("IMT:  " + frag.getImtString());

                System.out.println("\nResponse Table:");
                System.out.println(frag.responseTable());
            }

            // =========================================================
            // STEP 7: IMKEY GROUPING EXAMPLE
            // =========================================================

            SimpleImResponse example = fragilityLib.getModels().iterator().next();

            ImKey key = new ImKey(
                    example.getImtString(),
                    example.getImValues(),
                    ImDomain.LOG_IM
            );

            List<SimpleImResponse> grouped =
                    fragilityLib.getByImKey(key);

            System.out.println("\n=== GROUPED BY IM KEY ===");
            System.out.println("Reference: " + example.getName());

            for (SimpleImResponse f : grouped) {
                System.out.println(" - " + f.getName());
            }

            // =========================================================
            // STEP 8: LOOKUP BY NAME
            // =========================================================

            String testName = portfolioFragilities.iterator().next();

            SimpleImResponse byName = fragilityLib.getByName(testName);

            System.out.println("\n=== LOOKUP BY NAME ===");
            System.out.println("Name: " + testName);
            System.out.println(byName);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
