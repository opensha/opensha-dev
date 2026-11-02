package scratch.anne.sra.util;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Deterministic, order-preserving grouping utilities for portfolio processing.
 *
 * <p>All ordering is derived from input iteration order and preserved using
 * LinkedHashMap and LinkedHashSet.</p>
 */
public final class PortfolioGroupingUtils {

    private PortfolioGroupingUtils() {}

    // ---------------------------------------------------------------------
    // 1. SIMPLE GROUPING (K -> List<T>)
    // ---------------------------------------------------------------------

    /**
     * Groups items by a single deterministic key.
     *
     * Order rule:
     * - keys appear in first-seen order
     * - values preserve asset iteration order
     */
    public static <T, K> Map<K, List<T>> groupBy(
            Collection<T> items,
            Function<T, K> keyFn
    ) {
        Map<K, List<T>> map = new LinkedHashMap<>();

        for (T item : items) {
            K key = keyFn.apply(item);
            map.computeIfAbsent(key, k -> new ArrayList<>()).add(item);
        }

        return map.entrySet().stream()
                .collect(Collectors.toMap(
                        Map.Entry::getKey,
                        e -> List.copyOf(e.getValue()),
                        (a, b) -> a,
                        LinkedHashMap::new
                ));
    }

    // ---------------------------------------------------------------------
    // 2. NESTED INVERSION GROUPING (field -> value -> List<T>)
    // ---------------------------------------------------------------------

    /**
     * Inverts per-asset maps into:
     * field → value → List<T>
     *
     * Order rule:
     * - field order = first appearance in asset stream
     * - value order = first appearance within field
     */
    public static <T> Map<String, Map<String, List<T>>> invertNestedFields(
            Collection<T> items,
            Function<T, Map<String, String>> extractor
    ) {

        Map<String, Map<String, List<T>>> result = new LinkedHashMap<>();

        for (T item : items) {
            Map<String, String> fields = extractor.apply(item);

            if (fields == null) continue;

            for (Map.Entry<String, String> e : fields.entrySet()) {
                if (e.getValue() == null) continue;

                String field = e.getKey();
                String value = e.getValue();

                result
                        .computeIfAbsent(field, f -> new LinkedHashMap<>())
                        .computeIfAbsent(value, v -> new ArrayList<>())
                        .add(item);
            }
        }

        return result.entrySet().stream()
                .collect(Collectors.toMap(
                        Map.Entry::getKey,
                        fieldEntry -> fieldEntry.getValue().entrySet().stream()
                                .collect(Collectors.toMap(
                                        Map.Entry::getKey,
                                        v -> List.copyOf(v.getValue()),
                                        (a, b) -> a,
                                        LinkedHashMap::new
                                )),
                        (a, b) -> a,
                        LinkedHashMap::new
                ));
    }

    // ---------------------------------------------------------------------
    // 3. ORDERED DISTINCT KEY EXTRACTION (single-valued projection)
    // ---------------------------------------------------------------------

    /**
     * Extracts ordered distinct keys from a 1-to-1 mapping:
     *
     * T → K
     *
     * Used for:
     * - IMKeys
     * - SiteKeys
     * - any deterministic single-value projection
     */
    public static <T, K> List<K> orderedDistinctKeys(
            Collection<T> items,
            Function<T, K> keyFn
    ) {
        LinkedHashSet<K> seen = new LinkedHashSet<>();

        for (T item : items) {
            seen.add(keyFn.apply(item));
        }

        return List.copyOf(seen);
    }

    // ---------------------------------------------------------------------
    // 4. ORDERED DISTINCT FLATTENED KEYS (multi-key per item)
    // ---------------------------------------------------------------------

    /**
     * Extracts ordered distinct keys when each item contains multiple keys.
     *
     * Used for:
     * - additional field names
     * - flattened metadata keys
     */
    public static <T> List<String> orderedDistinctFlattenedKeys(
            Collection<T> items,
            Function<T, Set<String>> keySetFn
    ) {
        LinkedHashSet<String> seen = new LinkedHashSet<>();

        for (T item : items) {
            Set<String> keys = keySetFn.apply(item);
            if (keys != null) {
                seen.addAll(keys);
            }
        }

        return List.copyOf(seen);
    }

    // ---------------------------------------------------------------------
    // 5. ORDERED VALUE EXTRACTION (field-specific projection)
    // ---------------------------------------------------------------------

    /**
     * Extracts ordered distinct values for a specific field context.
     *
     * Used for:
     * - additionalFieldValues(field)
     */
    public static <T> List<String> orderedDistinctValues(
            Collection<T> items,
            Function<T, String> valueFn
    ) {
        LinkedHashSet<String> seen = new LinkedHashSet<>();

        for (T item : items) {
            String v = valueFn.apply(item);
            if (v != null) {
                seen.add(v);
            }
        }

        return List.copyOf(seen);
    }

    // ---------------------------------------------------------------------
    // 6. ORDERED KEY MAP HELPERS (for nested structures)
    // ---------------------------------------------------------------------

    public static <K1, K2, T> List<K1> orderedOuterKeys(
            Map<K1, Map<K2, List<T>>> nestedMap
    ) {
        return List.copyOf(nestedMap.keySet());
    }

    public static <K1, K2, T> Map<K1, List<K2>> orderedInnerKeys(
            Map<K1, Map<K2, List<T>>> nestedMap
    ) {
        Map<K1, List<K2>> result = new LinkedHashMap<>();

        for (Map.Entry<K1, Map<K2, List<T>>> e : nestedMap.entrySet()) {
            result.put(e.getKey(), List.copyOf(e.getValue().keySet()));
        }

        return result;
    }
}