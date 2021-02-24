package scratch.kevin.ucerf3.downDipSubSectTest;

import java.util.ArrayList;
import java.util.List;

import org.opensha.sha.earthquake.faultSysSolution.ruptures.FaultSubsectionCluster;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Jump;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.strategies.ClusterPermutationStrategy;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.primitives.Ints;

public class DownDipTestPermutationStrategy implements ClusterPermutationStrategy {
	
	private DownDipSubSectBuilder downDipBuilder;
	private static final boolean D = false;

	public DownDipTestPermutationStrategy(DownDipSubSectBuilder downDipBuilder) {
		this.downDipBuilder = downDipBuilder;
	}

	@Override
	public List<FaultSubsectionCluster> getPermutations(
			FaultSubsectionCluster fullCluster, FaultSection firstSection) {
		List<FaultSection> clusterSects = fullCluster.subSects;
		int myInd = fullCluster.subSects.indexOf(firstSection);
		Preconditions.checkState(myInd >= 0, "first section not found in cluster");
		List<FaultSection> newSects = new ArrayList<>();
		newSects.add(firstSection);
		
		List<FaultSubsectionCluster> permuations = new ArrayList<>();
		
		if (fullCluster.parentSectionID == downDipBuilder.getParentID()) {
			// this is a down-dip fault section, only build rectangular permutations
			int startCol = downDipBuilder.getColumn(firstSection);
			int startRow = downDipBuilder.getRow(firstSection);
			
			if (D) System.out.println("Building permutations from "+startRow+", "+startCol);
			
			int rows = downDipBuilder.getNumRows();
			int cols = downDipBuilder.getNumCols();
			
			// build down-dip first, starting with single row
			if (D) System.out.println("\tbuilding down-dip");
			for (int endRow=startRow; endRow<rows; endRow++) {
				// build to the right first (including single column)
				for (int endCol=startCol; endCol<cols; endCol++)
					permuations.add(buildRectangularPermutation(
							fullCluster, startRow, startCol, endRow, endCol));
				// build to the left
				for (int endCol=startCol; --endCol>=0;)
					permuations.add(buildRectangularPermutation(
							fullCluster, startRow, startCol, endRow, endCol));
			}
			
			// build up-dip
			if (D) System.out.println("\tbuilding up-dip");
			for (int endRow=startRow; --endRow>=0;) {
				// build to the right first (including single column)
				for (int endCol=startCol; endCol<cols; endCol++)
					permuations.add(buildRectangularPermutation(
							fullCluster, startRow, startCol, endRow, endCol));
				// build to the left
				for (int endCol=startCol; --endCol>=0;)
					permuations.add(buildRectangularPermutation(
							fullCluster, startRow, startCol, endRow, endCol));
			}
		} else {
			// just this section
			permuations.add(buildCopyJumps(fullCluster, newSects));
			
			// build toward the smallest ID
			for (int i=myInd; --i>=0;) {
				FaultSection nextSection = clusterSects.get(i);
				newSects.add(nextSection);
				permuations.add(buildCopyJumps(fullCluster, newSects));
			}
			newSects = new ArrayList<>();
			newSects.add(firstSection);
			// build toward the largest ID
			for (int i=myInd+1; i<clusterSects.size(); i++) {
				FaultSection nextSection = clusterSects.get(i);
				newSects.add(nextSection);
				permuations.add(buildCopyJumps(fullCluster, newSects));
			}
		}
		return permuations;
	}
	
	private static FaultSubsectionCluster buildCopyJumps(FaultSubsectionCluster fullCluster,
			List<FaultSection> subsetSects) {
		FaultSubsectionCluster permutation = new FaultSubsectionCluster(new ArrayList<>(subsetSects));
		for (FaultSection sect : subsetSects)
			for (Jump jump : fullCluster.getConnections(sect))
				permutation.addConnection(new Jump(sect, permutation,
						jump.toSection, jump.toCluster, jump.distance));
		return permutation;
	}
	
	private FaultSubsectionCluster buildRectangularPermutation(FaultSubsectionCluster fullCluster,
			int startRow, int startCol, int endRow, int endCol) {
		List<FaultSection> subsetSects = new ArrayList<>();
		
		// this is a list of all exit points from this rupture (sections from which we can jump to another
		// fault without it being considered a splay jump). we'll define that now as the far edges from the
		// starting point, though other strategies likely exist
		List<FaultSection> exitPoints = new ArrayList<>();
		
		if (D) System.out.println("\t\trow span: "+startRow+" => "+endRow+": "+printIndexes(startRow, endRow));
		if (D) System.out.println("\t\tcol span: "+startCol+" => "+endCol+": "+printIndexes(startCol, endCol));
		
		// build with rows in the outer loop, as this should make azimuth calculations work
		for (int row : indexes(startRow, endRow)) {
			for (int col : indexes(startCol, endCol)) {
//				if (D) System.out.println("row="+row+", col="+col);
				FaultSection sect = downDipBuilder.getSubSect(row, col);
				Preconditions.checkNotNull(sect);
				subsetSects.add(sect);
				if (row == endRow && col == endCol)
					exitPoints.add(sect);
			}
		}
		Preconditions.checkState(subsetSects.get(0).equals(downDipBuilder.getSubSect(startRow, startCol)));
		FaultSubsectionCluster permutation = new FaultSubsectionCluster(subsetSects, exitPoints);
		// add possible jumps out of this permutation
		for (FaultSection sect : subsetSects)
			for (Jump jump : fullCluster.getConnections(sect))
				permutation.addConnection(new Jump(sect, permutation,
						jump.toSection, jump.toCluster, jump.distance));
		return permutation;
	}
	
	private static int[] indexes(int start, int end) {
		if (start <= end) {
			int[] ret = new int[1+end-start];
			for (int i=0; i<ret.length; i++)
				ret[i] = start+i;
			return ret;
		}
		int[] ret = new int[1+start-end];
		for (int i=0; i<ret.length; i++)
			ret[i] = start-i;
		return ret;
	}
	
	private static String printIndexes(int start, int end) {
		int[] indexes = indexes(start, end);
		return Joiner.on(",").join(Ints.asList(indexes));
	}

	@Override
	public String getName() {
		return "Down Dip Test";
	}

}
