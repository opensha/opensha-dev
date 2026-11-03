package scratch.anne.sra.domain.structural_response.fragilities;

import java.util.List;

import scratch.anne.sra.domain.structural_response.NamedResponseModel;
import scratch.anne.sra.domain.structural_response.fragilities.definitions.LimitStateFragility;
import scratch.anne.sra.util.StringUtil;
import scratch.anne.sra.util.enums.IMT;
import scratch.anne.sra.util.enums.LimitState;

/**
 * A named fragility response model: one component ID and IMT, with one
 * {@link LimitStateFragility} per limit state it defines (e.g. SLIGHT
 * through COLLAPSE). {@link #getPrimaryLimitState()} picks out the
 * single limit state used when the model is prepared for convolution
 * (see {@link SimpleImFragilityLibraryPreparer}).
 */
public class FragilityModel implements NamedResponseModel {

    private final String name;
    private final IMT imt;
    private final Double period;
    private final String imtString;

    private final List<LimitStateFragility> limitStates;
    private final LimitState primaryLs;

    public FragilityModel(
            String name,
            IMT imt,
            Double period,
            List<LimitStateFragility> limitStates,
            LimitState primaryLs) {

        if (limitStates == null || limitStates.isEmpty())
            throw new IllegalArgumentException("At least one limit state required");

        this.name = name;
        this.imt = imt;
        this.period = period;
        this.imtString = StringUtil.imtToString(imt, period);
        this.limitStates = List.copyOf(limitStates);
        
        if (primaryLs != null && limitStates.stream().noneMatch(ls -> ls.getLimitState() == primaryLs)) {
            throw new IllegalArgumentException("Primary limit state must exist in limitStates list");
        }
        this.primaryLs = primaryLs;

        validate();
    }
    
    //* convenience constructor */
    public FragilityModel(
            String name,
            IMT imt,
            Double period,
            List<LimitStateFragility> limitStates) {
    	this(name, imt, period, limitStates, null);
    }

    private void validate() {
        NamedResponseModel.validateImtPeriod(imt, period);
    }

    public List<LimitStateFragility> getLimitStateFragilities() {
        return limitStates;
    }
    
    // Getter for primary limit state (nullable fallback to last LS)
    public LimitStateFragility getPrimaryLimitState() {
        if (primaryLs != null) {
            return limitStates.stream()
                    .filter(ls -> ls.getLimitState().equals(primaryLs))
                    .findFirst()
                    .orElseThrow(() -> new IllegalStateException("Primary limit state not found"));
        } else {
            return limitStates.get(limitStates.size() - 1); // last limit state as fallback
        }
    }

    // Getter for a specific limit state by name
    public LimitStateFragility getLimitState(String lsName) {
        return limitStates.stream()
                .filter(ls -> ls.getLimitState().name().equals(lsName))
                .findFirst()
                .orElse(null);
    }

    @Override
    public String getName() { return name; }

    @Override
    public IMT getImt() { return imt; }

    @Override
    public Double getPeriod() { return period; }

    @Override
    public String getImtString() { return imtString; }
    
	@Override
	public String toString() {
	  return "FragilityModel{name='" + name + "', imt=" + imtString + "}";
	}
	
	public String toVerboseString() {

	    StringBuilder sb = new StringBuilder();

	    sb.append("Fragility Model: ").append(name).append("\n");
	    sb.append("IMT: ").append(imtString).append("\n");
	    sb.append("Primary LS: ").append(primaryLs).append("\n");
	    sb.append("Limit States:\n");

	    for (LimitStateFragility lsFrag : limitStates) {

	        sb.append("\n");
	        sb.append("  Limit State: ")
	          .append(lsFrag.getLimitState())
	          .append("\n");
	        
	        sb.append(lsFrag.getDefinition().toVerboseString());
	    }
	    return sb.toString();
	}
}