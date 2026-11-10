package scratch.anne.sra.util.enums;
/**
 * Enum defining whether a limit state fragility is defined
 * as:
 *  * <ul>
 *   <li>discrete probability values at specified intensities</li>
 *   <li>a parametric distribution</li>
	 * <ul>
	 *   <li>lognormal distribution</li>
	 *   <li>another type of distribution that could be included here</li>
	 * </ul>  
 * </ul>
 */
public enum FragilityType {
    DISCRETE,
    LOGNORMAL
}
