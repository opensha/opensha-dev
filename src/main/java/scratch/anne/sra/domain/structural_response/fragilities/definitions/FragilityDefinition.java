package scratch.anne.sra.domain.structural_response.fragilities.definitions;

/**
 * Strategy interface for a fragility curve's shape -- e.g.
 * {@link DiscreteFragilityDefinition} (explicit points) or
 * {@link LognormalFragilityDefinition} (median/beta). Implementations
 * provide no computation here beyond a verbose string representation;
 * see {@link scratch.anne.sra.domain.structural_response.fragilities.SimpleImFragilityLibraryPreparer}
 * for how a definition is turned into IM/response arrays.
 */
public interface FragilityDefinition {
    String toVerboseString();
}
