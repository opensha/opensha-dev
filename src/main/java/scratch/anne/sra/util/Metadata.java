package scratch.anne.sra.util;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

/**
 * Generic immutable key-value metadata storage with builder support.
 * Can be used anywhere in the system.
 */
public final class Metadata {

    private final Map<String, String> values;

    private Metadata(Map<String, String> values) {
        this.values = Collections.unmodifiableMap(new HashMap<>(values));
    }

    /** Generic access (dictionary-like) */
    public String get(String key) {
        return values.get(key);
    }

    /** Safe access with default */
    public String getOrDefault(String key, String defaultValue) {
        return values.getOrDefault(key, defaultValue);
    }

    public Map<String,String> asMap() {
        return values;
    }

    @Override
    public String toString() {
        return values.toString();
    }

    // ------------------- Builder -------------------

    public static class Builder {
        private final Map<String, String> values = new HashMap<>();

        public Builder set(String key, String value) {
            if (key != null && value != null) {
                values.put(key, value);
            }
            return this;
        }

        public Metadata build() {
            return new Metadata(values);
        }
    }
}
