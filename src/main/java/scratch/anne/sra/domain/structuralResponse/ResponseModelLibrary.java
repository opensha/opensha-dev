package scratch.anne.sra.domain.structuralResponse;

import java.util.*;

import scratch.anne.sra.util.Metadata;
import scratch.anne.sra.util.enums.IMT;

/**
 * Generic library for named response models (vulnerabilities, fragilities, etc.).
 * Stores models, metadata, and supports efficient lookup by name or IMT.
 */
public class ResponseModelLibrary<T extends NamedResponseModel> implements LibraryGetters<T> {

    private final List<T> models;
    private final Map<String, T> nameMap;

    // IMT lookup
    private final Map<IMT, List<T>> imtMap;
    private final Map<String, List<T>> imtStringMap;

    // Sets of unique IMTs and IMT strings
    private final Set<IMT> imts;
    private final Set<String> imtStrings;

    private final Metadata metadata;

    // ----------------Private Constructor --------------------
    private ResponseModelLibrary(List<T> models, Metadata metadata) {
        this.models = Collections.unmodifiableList(new ArrayList<>(models));
        this.metadata = metadata;

        // Build name → model map
        Map<String, T> nMap = new LinkedHashMap<>();
        for (T m : models) nMap.put(m.getName(), m);
        this.nameMap = Collections.unmodifiableMap(nMap);

        // Build IMT maps
        Map<IMT, List<T>> imtM = new LinkedHashMap<>();
        Map<String, List<T>> imtStrM = new LinkedHashMap<>();
        Set<IMT> imtSet = new LinkedHashSet<>();
        Set<String> imtStrSet = new LinkedHashSet<>();

        for (T m : models) {
            // IMT
            imtSet.add(m.getImt());
            imtM.computeIfAbsent(m.getImt(), k -> new ArrayList<>()).add(m);

            // IMT string
            String s = m.getImtString();
            imtStrSet.add(s);
            imtStrM.computeIfAbsent(s, k -> new ArrayList<>()).add(m);
        }

        this.imtMap = Collections.unmodifiableMap(imtM);
        this.imtStringMap = Collections.unmodifiableMap(imtStrM);
        this.imts = Collections.unmodifiableSet(imtSet);
        this.imtStrings = Collections.unmodifiableSet(imtStrSet);
    }
    
    // ---------------- STATIC FACTORY ----------------
    public static <T extends NamedResponseModel>
    ResponseModelLibrary<T> of(List<T> models, Metadata metadata) {
        return new ResponseModelLibrary<T>(models, metadata);
    }

    // -------------------- Accessors --------------------
    @Override
    public List<T> getModels() { return models; }

    @Override
    public Set<String> getModelNames() { return nameMap.keySet(); }

    @Override
    public T getModelByName(String name) {
        T model = nameMap.get(name);
        if (model == null) throw new NoSuchElementException("No model with name: " + name);
        return model;
    }

    @Override
    public Set<IMT> getIMTs() { return imts; }

    @Override
    public List<T> getByIMT(IMT imt) { 
        return imtMap.getOrDefault(imt, Collections.emptyList()); 
    }

    @Override
    public Set<String> getIMTStrings() { return imtStrings; }

    @Override
    public List<T> getByIMTString(String imtString) { 
        return imtStringMap.getOrDefault(imtString, Collections.emptyList()); 
    }

    @Override
    public Metadata getMetadata() { return metadata; }

    @Override
    public int size() { return models.size(); }

    @Override
    public String toString() {
        return "ResponseModelLibrary{" +
               "nItems=" + size() +
               ", nIMTs=" + imts.size() +
               ", nIMTStrings=" + imtStrings.size() +
               ", metadata=" + metadata +
               '}';
    }

    // -------------------- Builder --------------------
    public static <T extends NamedResponseModel> Builder<T> builder() {
        return new Builder<>();
    }

    public static class Builder<T extends NamedResponseModel> {
        private List<T> models = new ArrayList<>();
        private Metadata metadata = new Metadata.Builder().build();

        public Builder<T> models(Collection<T> models) {
            this.models = new ArrayList<>(models);
            return this;
        }

        public Builder<T> metadata(Metadata metadata) {
            this.metadata = metadata;
            return this;
        }

        public ResponseModelLibrary<T> build() {
            return new ResponseModelLibrary<>(models, metadata);
        }
    }
}