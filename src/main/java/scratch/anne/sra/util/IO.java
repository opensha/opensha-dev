package scratch.anne.sra.util;

import java.nio.file.Files;
import java.nio.file.Path;

public class IO {
	
    public static void verifyFileExists(Path path, String description) {
        if (!Files.exists(path)) {
            throw new IllegalArgumentException(description + " does not exist: " + path);
        }
    }

}
