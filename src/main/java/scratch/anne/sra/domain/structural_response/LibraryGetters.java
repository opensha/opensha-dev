package scratch.anne.sra.domain.structural_response;

import java.util.List;
import java.util.Set;

import scratch.anne.sra.util.Metadata;
import scratch.anne.sra.util.enums.IMT;

public interface LibraryGetters<T extends NamedResponseModel> {
    List<T> getModels();				// all models
    Set<String> getModelNames();		// all model names
    T getModelByName(String name);		// single model
    
    Set<IMT> getIMTs();         
    List<T> getByIMT(IMT imt); 
    Set<String> getIMTStrings();               
    List<T> getByIMTString(String imtString);  

    Metadata getMetadata();             // library metadata
    int size();
}
