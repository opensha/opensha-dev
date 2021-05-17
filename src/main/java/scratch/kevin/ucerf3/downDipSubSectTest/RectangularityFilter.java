package scratch.kevin.ucerf3.downDipSubSectTest;

import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.FaultSubsectionCluster;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Jump;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.PlausibilityFilter;
import org.opensha.sha.faultSurface.FaultSection;

import scratch.UCERF3.inversion.laughTest.PlausibilityResult;

class RectangularityFilter implements PlausibilityFilter {

	private DownDipSubSectBuilder builder;
	private int minDimension;
	private double maxAspectRatio;

	RectangularityFilter(DownDipSubSectBuilder builder, int minDimension, double maxAspectRatio) {
		this.builder = builder;
		this.minDimension = minDimension;
		this.maxAspectRatio = maxAspectRatio;
	}

	@Override
	public String getShortName() {
		return "Rectangularity";
	}

	@Override
	public String getName() {
		return getShortName();
	}
	
	private PlausibilityResult apply(FaultSubsectionCluster cluster, boolean verbose) {
		if (cluster.parentSectionID != builder.getParentID())
			return PlausibilityResult.PASS;
		int maxRow = 0;
		int minRow = Integer.MAX_VALUE;
		int maxCol = 0;
		int minCol = Integer.MAX_VALUE;
		for (FaultSection sect : cluster.subSects) {
			int row = builder.getRow(sect);
			int col = builder.getColumn(sect);
			maxRow = Integer.max(maxRow, row);
			minRow = Integer.min(minRow, row);
			maxCol = Integer.max(maxCol, col);
			minCol = Integer.min(minCol, col);
		}
		int rowSpan = 1 + maxRow - minRow;
		int colSpan = 1 + maxCol - minCol;
		if (verbose)
			System.out.println(getShortName()+": testing with rowSpan="+rowSpan+" and colSpan="+colSpan);
		if (rowSpan < minDimension || colSpan < minDimension) {
			if (verbose)
				System.out.println(getShortName()+": failing because below min dimension of "+minDimension);
			return PlausibilityResult.FAIL_HARD_STOP;
		}
		double aspectRatio = Math.max((double)rowSpan/(double)colSpan, (double)colSpan/(double)rowSpan);
		if (aspectRatio > maxAspectRatio) {
			if (verbose)
				System.out.println(getShortName()+": failing because of aspect ratio of "+aspectRatio);
			return PlausibilityResult.FAIL_HARD_STOP;
		}
		// if it's rectangular, then the count will be rowSpan x colSpan
		int expectedNum = rowSpan*colSpan;
		if (cluster.subSects.size() == expectedNum) {
			if (verbose)
				System.out.println(getShortName()+": passing with exact match of "+expectedNum+" sects");
			return PlausibilityResult.PASS;
		}
		if (verbose)
			System.out.println(getShortName()+": failing because of hole(s). have "
					+cluster.subSects.size()+", complete would be "+expectedNum);
		return PlausibilityResult.FAIL_HARD_STOP;
	}

	@Override
	public PlausibilityResult apply(ClusterRupture rupture, boolean verbose) {
		PlausibilityResult result = PlausibilityResult.PASS;
		for (FaultSubsectionCluster cluster : rupture.clusters) {
			result = result.logicalAnd(apply(cluster, verbose));
			if (!result.canContinue())
				return result;
		}
		for (ClusterRupture splay : rupture.splays.values()) {
			result = result.logicalAnd(apply(splay, verbose));
			if (!result.canContinue())
				return result;
		}
		return result;
	}

}
