package scratch.anne.sra.examples.PBR;

import java.nio.file.*;
import java.util.*;

import scratch.anne.sra.domain.structuralResponse.ResponseModelLibrary;
import scratch.anne.sra.domain.structuralResponse.SimpleImResponse;
import scratch.anne.sra.domain.structuralResponse.SimpleImResponseLibrary;
import scratch.anne.sra.domain.structuralResponse.fragilities.FragilityModel;
import scratch.anne.sra.domain.structuralResponse.fragilities.SimpleImFragilityLibraryPreparer;
import scratch.anne.sra.domain.structuralResponse.fragilities.SimpleImFragilityLibraryPreparer.PrepareSpec;
import scratch.anne.sra.io.readers.FragilityLibraryReader;
import scratch.anne.sra.util.AssetKeys.ImKey;
import scratch.anne.sra.util.AssetKeys.ImKey.ImDomain;
import scratch.anne.sra.util.NumericUtil.ImArrayParam;

/**
 * Demonstrates reading the PBR LJB1 fragility (via
 * {@link FragilityLibraryReader#readLibrary}) and preparing it into a
 * {@link SimpleImResponseLibrary} via {@link SimpleImFragilityLibraryPreparer},
 * then inspecting the prepared library the same way as
 * {@code SimpleFragilityLibraryReaderExample}.
 *
 * <p><b>Choice point -- edit before running:</b> STEP 4's {@code spec}
 * selects which {@code PrepareSpec} is used -- {@code defaultSpec} (no
 * filtering/resampling), {@code generatedSpec} (log-spaced IM
 * resampling), or {@code customImSpec} (the exact IM values Devin
 * McPhillips provided for this fragility, recreated explicitly since
 * {@code LJB1_fragility.json} only stores the inferred lognormal fit,
 * not these original points).</p>
 */
public class LBJ1Reader {

    public static void main(String[] args) {

        try {

            // =========================================================
            // STEP 1: FILE PATH SETUP (JSON or CSV)
            // =========================================================

            Path resourceFolder = Path.of(
                    "C:/Users/ahulsey/git/opensha-dev/src/main/resources/scratch/anne/sra/PBR"
            );

            Path fragLibraryJSON = resourceFolder.resolve("LJB1_fragility.json");

            // Choose input format
            Path filePath = fragLibraryJSON;

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
            System.out.println(modelLibrary.getModelNames());

            // =========================================================
            // STEP 3: DEFINE PORTFOLIO FILTER
            // =========================================================

            Set<String> portfolioFragilities = new HashSet<>(Arrays.asList(
                    "LJB1"
            ));

            // =========================================================
            // STEP 4: BUILD PREPARATION SPEC
            // =========================================================

            PrepareSpec defaultSpec = PrepareSpec.builder()
                    .build();
            
            PrepareSpec generatedSpec = PrepareSpec.builder()
                    .filterNames(portfolioFragilities)
                    .selectedLS("Topple")
                    .imParam(
                            ImArrayParam.builder()
                                    .logStep(0.01)
                                    .spread(5.0)
                                    .build()
                    )
                    .build();
            
            // these are the im values given by Devin McPhilips for the example fragility
            // the LJB1_fragility.json only contains the inferred lognormal distribution
            // so this array will recreate the actual values he provided
            PrepareSpec customImSpec = PrepareSpec.builder()
                    .filterNames(portfolioFragilities)
                    .selectedLS("Topple")
                    .customImValues(new double[] {
                    	    0.02, 0.020562691, 0.021141213, 0.021736011, 0.022347544, 0.022976282,
                    	    0.023622709, 0.024287323, 0.024970635, 0.025673173, 0.026395476, 0.0271381,
                    	    0.027901618, 0.028686618, 0.029493702, 0.030323494, 0.031176632, 0.032053772,
                    	    0.03295559, 0.03388278, 0.034836057, 0.035816153, 0.036823824, 0.037859846,
                    	    0.038925015, 0.040020153, 0.041146101, 0.042303728, 0.043493924, 0.044717605,
                    	    0.045975715, 0.04726922, 0.048599118, 0.049966432, 0.051372215, 0.052817548,
                    	    0.054303546, 0.055831351, 0.057402141, 0.059017123, 0.060677543, 0.062384678,
                    	    0.064139842, 0.065944387, 0.067799702, 0.069707216, 0.071668396, 0.073684754,
                    	    0.07575784, 0.077889252, 0.080080631, 0.082333663, 0.084650082, 0.087031674,
                    	    0.08948027, 0.091997756, 0.094586071, 0.097247206, 0.099983212, 0.102796194,
                    	    0.105688317, 0.10866181, 0.11171896, 0.114862121, 0.118093714, 0.121416227,
                    	    0.124832217, 0.128344314, 0.131955222, 0.135667722, 0.139484671, 0.143409008,
                    	    0.147443754, 0.151592017, 0.155856988, 0.160241953, 0.164750287, 0.169385461,
                    	    0.174151043, 0.179050702, 0.184088212, 0.189267449, 0.194592402, 0.20006717,
                    	    0.205695968, 0.211483129, 0.21743311, 0.223550491, 0.229839981, 0.236306423,
                    	    0.242954796, 0.249790218, 0.256817951, 0.264043406, 0.271472146, 0.27910989,
                    	    0.286962519, 0.295036078, 0.303336782, 0.311871023, 0.320645371, 0.329666582,
                    	    0.338941599, 0.348477566, 0.358281822, 0.368361917, 0.37872561, 0.389380881,
                    	    0.400335933, 0.411599201, 0.423179356, 0.435085312, 0.447326238, 0.459911556,
                    	    0.472850957, 0.486154401, 0.499832132, 0.51389468, 0.52835287, 0.543217836,
                    	    0.55850102, 0.57421419, 0.590369442, 0.606979216, 0.624056297, 0.641613835,
                    	    0.659665345, 0.678224727, 0.697306268, 0.71692466, 0.737095006, 0.757832835,
                    	    0.779154114, 0.801075257, 0.823613142, 0.84678512, 0.87060903, 0.895103215,
                    	    0.920286533, 0.946178372, 0.972798666, 1.00016791, 1.028307175, 1.057238125,
                    	    1.086983035, 1.117564803, 1.149006976, 1.181333759, 1.214570042, 1.248741413,
                    	    1.283874179, 1.31999539, 1.357132854, 1.395315163, 1.434571715, 1.474932731,
                    	    1.516429286, 1.559093328, 1.602957703, 1.648056182, 1.694423486, 1.742095313,
                    	    1.791108364, 1.841500376, 1.893310143, 1.946577555, 2.001343621, 2.057650504,
                    	    2.115541556, 2.175061347, 2.236255699, 2.299171726, 2.363857867, 2.430363923,
                    	    2.498741096, 2.569042029, 2.641320847, 2.715633196, 2.792036289, 2.870588948,
                    	    2.95135165, 3.034386574, 3.119757647, 3.207530596, 3.297772996, 3.390554325,
                    	    3.485946014, 3.584021505, 3.684856304, 3.788528044, 3.895116541, 4.004703856,
                    	    4.117374361, 4.233214798, 4.352314352, 4.474764718, 4.600660169, 4.73009763,
                    	    4.863176755, 5.0
                    	})
                    .build();
            
            // select one of the above specs
            PrepareSpec spec = defaultSpec;

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
