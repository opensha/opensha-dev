package scratch.anne.sra.domain.structuralResponse;

import java.util.*;

import scratch.anne.sra.util.Metadata;
import scratch.anne.sra.util.AssetKeys.ImKey;
import scratch.anne.sra.util.enums.IMT;

/**
 * Library of {@link SimpleImResponse} objects prepared for convolution.
 *
 * <p>Supports efficient lookup by:
 * <ul>
 *   <li>Response name</li>
 *   <li>IMKey (IMT string + IM values array)</li>
 * </ul>
 *
 * <p>This class assumes all responses are already validated and
 * consistent (e.g., via a preparer). It does not modify input objects.
 *
 * <p>All internal maps are built once at construction and exposed as unmodifiable views.
 */
public class SimpleImResponseLibrary implements LibraryGetters<SimpleImResponse> {

    private final Map<String, SimpleImResponse> nameMap;
    private final Map<ImKey, List<SimpleImResponse>> imKeyMap;
    private final Metadata metadata;

    // -----------------------
    // Constructors
    // -----------------------

    public SimpleImResponseLibrary(List<SimpleImResponse> responses, Metadata metadata) {
        this.metadata = metadata;

        Map<String, SimpleImResponse> nameMapTemp = new LinkedHashMap<>();
        Map<ImKey, List<SimpleImResponse>> imKeyMapTemp = new HashMap<>();
        initMaps(responses, nameMapTemp, imKeyMapTemp);

        this.nameMap = Collections.unmodifiableMap(nameMapTemp);
        this.imKeyMap = Collections.unmodifiableMap(imKeyMapTemp);
    }

    public SimpleImResponseLibrary(List<SimpleImResponse> responses) {
        this(responses, null);
    }

    /**
     * Static factory method to create a library.
     */
    public static SimpleImResponseLibrary of(List<SimpleImResponse> responses, Metadata metadata) {
        return new SimpleImResponseLibrary(responses, metadata);
    }

    private void initMaps(
            List<SimpleImResponse> responses,
            Map<String, SimpleImResponse> nameMapTemp,
            Map<ImKey, List<SimpleImResponse>> imKeyMapTemp) {

        for (SimpleImResponse response : responses) {
            // ---- name map ----
            String name = response.getName();
            if (name == null || name.isEmpty())
                throw new IllegalArgumentException("Response must have a non-null name");
            if (nameMapTemp.containsKey(name))
                throw new IllegalArgumentException("Duplicate response name: " + name);
            nameMapTemp.put(name, response);

            // ---- imKey map ----
            ImKey key = response.getImKey();
            imKeyMapTemp.computeIfAbsent(key, k -> new ArrayList<>()).add(response);
        }

        // wrap inner lists
        for (var entry : imKeyMapTemp.entrySet()) {
            entry.setValue(Collections.unmodifiableList(entry.getValue()));
        }
    }

    // -----------------------
    // Accessors
    // -----------------------

    public SimpleImResponse getByName(String name) {
        return nameMap.get(name);
    }

    public List<SimpleImResponse> getByImKey(ImKey key) {
        return imKeyMap.getOrDefault(key, Collections.emptyList());
    }

    public Map<ImKey, List<SimpleImResponse>> getImKeyMap() {
        return imKeyMap;
    }

    public int size() {
        return nameMap.size();
    }

    @Override
    public List<SimpleImResponse> getModels() {
        return Collections.unmodifiableList(new ArrayList<>(nameMap.values()));
    }

    @Override
    public Set<String> getModelNames() {
        return nameMap.keySet();
    }

    @Override
    public SimpleImResponse getModelByName(String name) {
        return getByName(name);
    }

    @Override
    public Set<IMT> getIMTs() {
        Set<IMT> set = new LinkedHashSet<>();
        for (SimpleImResponse v : nameMap.values())
            set.add(v.getImt());
        return set;
    }

    @Override
    public List<SimpleImResponse> getByIMT(IMT imt) {
        List<SimpleImResponse> list = new ArrayList<>();
        for (SimpleImResponse v : nameMap.values()) {
            if (v.getImt() == imt) list.add(v);
        }
        return list;
    }

    @Override
    public Set<String> getIMTStrings() {
        Set<String> set = new LinkedHashSet<>();
        for (SimpleImResponse v : nameMap.values())
            set.add(v.getImtString());
        return set;
    }

    @Override
    public List<SimpleImResponse> getByIMTString(String imtString) {
        List<SimpleImResponse> list = new ArrayList<>();
        for (SimpleImResponse v : nameMap.values()) {
            if (v.getImtString().equals(imtString)) list.add(v);
        }
        return list;
    }

    @Override
    public Metadata getMetadata() {
        return metadata;
    }
    
    @Override
    public String toString() {
        return "SimpleImResponseLibrary[size=" + size()
                + ", imKeys=" + imKeyMap.size()
                + ", metadata=" + (metadata != null ? "present" : "none")
                + "]";
    }
    
    public String summary() {

        StringBuilder sb = new StringBuilder();

        sb.append("Simple IM Response Library Summary\n");
        sb.append("----------------------------------\n");

        sb.append("Total responses : ").append(size()).append("\n");
        sb.append("Unique IM grids : ").append(imKeyMap.size()).append("\n");

        sb.append("IMTs            : ").append(getIMTs()).append("\n");

        if (metadata != null)
            sb.append("Metadata        : ").append(metadata).append("\n");
        else
            sb.append("Metadata        : none\n");

        sb.append("\nResponses by IMKey:\n");

        for (Map.Entry<ImKey, List<SimpleImResponse>> e : imKeyMap.entrySet()) {
            sb.append("  ")
              .append(e.getKey())
              .append(" -> ")
              .append(e.getValue().size())
              .append(" responses\n");
        }

        return sb.toString();
    }
}