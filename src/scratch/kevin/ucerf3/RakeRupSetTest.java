package scratch.kevin.ucerf3;

import java.util.List;

import com.google.common.collect.Lists;

import scratch.UCERF3.inversion.InversionFaultSystemRupSetFactory;
import scratch.UCERF3.inversion.laughTest.LaughTestFilter;
import scratch.UCERF3.logicTree.LogicTreeBranch;

public class RakeRupSetTest {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		LogicTreeBranch branch = LogicTreeBranch.DEFAULT;
		double defaultAseismicityValue = 0.1;
		LaughTestFilter laughTest = LaughTestFilter.getDefault();
		
//		List<Double> cmlRakes = Lists.newArrayList(0d, 30d, 60d, 90d, 120d, 150d, 180d, 210d,
//				240d, 270d, 300d, 330d, 360d, Double.POSITIVE_INFINITY);
		List<Double> cmlRakes = Lists.newArrayList(Double.POSITIVE_INFINITY);
		List<Integer> rupCounts = Lists.newArrayList();
		for (double cmlRake : cmlRakes) {
			laughTest.setMaxCmlRakeChange(cmlRake);
			rupCounts.add(InversionFaultSystemRupSetFactory.forBranch(
					laughTest, defaultAseismicityValue, branch).getNumRuptures());
		}
		
		System.out.println("\nMaxCMLRakeChange\tNum Rups");
		for (int i=0; i<cmlRakes.size(); i++) {
			System.out.println(cmlRakes.get(i)+"\t"+rupCounts.get(i));
		}
	}

}
