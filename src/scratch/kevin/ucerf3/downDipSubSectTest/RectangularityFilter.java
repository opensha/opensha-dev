package scratch.kevin.ucerf3.downDipSubSectTest;

import java.util.List;

import org.opensha.commons.util.IDPairing;
import org.opensha.sha.faultSurface.FaultSection;

import scratch.UCERF3.inversion.laughTest.AbstractPlausibilityFilter;
import scratch.UCERF3.inversion.laughTest.PlausibilityResult;

class RectangularityFilter extends AbstractPlausibilityFilter {

	private DownDipSubSectBuilder builder;
	private int minDimension;

	RectangularityFilter(DownDipSubSectBuilder builder, int minDimension) {
		this.builder = builder;
		this.minDimension = minDimension;
	}

	@Override
	public String getShortName() {
		return "Rectangularity";
	}

	@Override
	public String getName() {
		return getShortName();
	}

	@Override
	public PlausibilityResult applyLastSection(List<? extends FaultSection> rupture, List<IDPairing> pairings,
			List<Integer> junctionIndexes) {
		int maxRow = 0;
		int minRow = Integer.MAX_VALUE;
		int maxCol = 0;
		int minCol = Integer.MAX_VALUE;
		int matchCount = 0;
		for (FaultSection sect : rupture) {
			if (sect.getParentSectionId() == builder.getParentID()) {
				int row = builder.getRow(sect);
				int col = builder.getColumn(sect);
				maxRow = Integer.max(maxRow, row);
				minRow = Integer.min(minRow, row);
				maxCol = Integer.max(maxCol, col);
				minCol = Integer.min(minCol, col);
				matchCount++;
			}
		}
		if (matchCount == 0)
			return PlausibilityResult.PASS;
		int rowSpan = 1 + maxRow - minRow;
		int colSpan = 1 + maxCol - minCol;
		if (rowSpan < minDimension || colSpan < minDimension)
			return PlausibilityResult.FAIL_FUTURE_POSSIBLE;
		// if it's rectangular, then the count will be rowSpan x colSpan
		int expectedNum = rowSpan*colSpan;
		if (matchCount == expectedNum)
			return PlausibilityResult.PASS;
		// if multiple columns have holes, it will never pass, bail here
		int minToContinue = (colSpan-1)*rowSpan + 1;
		if (expectedNum < minToContinue)
			return PlausibilityResult.FAIL_HARD_STOP;
		return PlausibilityResult.FAIL_FUTURE_POSSIBLE;
	}

	@Override
	public boolean isApplyJunctionsOnly() {
		// applies within parent fault sections
		return false;
	}

}
