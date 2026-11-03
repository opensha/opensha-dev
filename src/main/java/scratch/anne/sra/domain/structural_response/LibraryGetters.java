package scratch.anne.sra.domain.structural_response;

import java.util.List;
import java.util.Set;

import scratch.anne.sra.util.Metadata;
import scratch.anne.sra.util.enums.IMT;

/**
 * Common read access for a library of {@link NamedResponseModel}s
 * (e.g. {@code ResponseModelLibrary}, {@code SimpleImResponseLibrary}),
 * by name or by intensity measure type.
 */
public interface LibraryGetters<T extends NamedResponseModel> {
    List<T> getModels();				// all models
    Set<String> getModelNames();		// all model names
    T getModelByName(String name);		// single model, throws/null if absent (see implementation)
    
    Set<IMT> getIMTs();         		// distinct IMTs across all models
    List<T> getByIMT(IMT imt); 			// models with the given IMT
    Set<String> getIMTStrings();       // distinct IMT strings (e.g. "SA(1.0)")
    List<T> getByIMTString(String imtString);  // models with the given IMT string

    Metadata getMetadata();             // library metadata
    int size();							// number of models
}
