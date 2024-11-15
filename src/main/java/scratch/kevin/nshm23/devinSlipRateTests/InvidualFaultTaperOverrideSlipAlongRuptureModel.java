package scratch.kevin.nshm23.devinSlipRateTests;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.modules.SlipAlongRuptureModel;

import com.google.common.base.Preconditions;

public class InvidualFaultTaperOverrideSlipAlongRuptureModel extends SlipAlongRuptureModel.NamedSlipAlongRuptureModel {

	private Set<Integer> parentsForTaper;
	private boolean splitParents;
	private static SlipAlongRuptureModel.Tapered TAPERED = new SlipAlongRuptureModel.Tapered();

	public InvidualFaultTaperOverrideSlipAlongRuptureModel(Set<Integer> parentsForTaper, boolean splitParents) {
		this.parentsForTaper = parentsForTaper;
		this.splitParents = splitParents;
	}

	@Override
	public String getName() {
		return "Individual Fault Taper-Override Dsr (splitParents="+splitParents+")";
	}

	@Override
	public double[] calcSlipOnSectionsForRup(FaultSystemRupSet rupSet, int rthRup, double[] sectArea, double aveSlip) {
		List<List<Integer>> parentBundles = new ArrayList<>();
		List<Boolean> bundleTapers = new ArrayList<>();
		boolean anyTaper = false;
		
		int curParent = -1;
		List<Integer> curBundle = null;
		
		int sectCount = 0;
		for (int sectIndex : rupSet.getSectionsIndicesForRup(rthRup)) {
			int parent = rupSet.getFaultSectionData(sectIndex).getParentSectionId();
			Preconditions.checkState(parent >= 0);
			if (parent != curParent) {
				curParent = parent;
				curBundle = new ArrayList<>();
				parentBundles.add(curBundle);
				boolean bundle = parentsForTaper.contains(parent);
				bundleTapers.add(bundle);
				anyTaper |= bundle;
			}
			curBundle.add(sectIndex);
			sectCount++;
		}
		Preconditions.checkState(sectCount == sectArea.length);
		
		if (!anyTaper)
			// all uniform
			return calcUniformSlipAlong(sectCount, aveSlip);
		
		if (parentBundles.size() > 1 && !splitParents) {
			// re-combine any contiguous tapered bundles
			
			List<List<Integer>> combParentBundles = new ArrayList<>(parentBundles.size());
			List<Boolean> compBundleTapers = new ArrayList<>(parentBundles.size());
			
			combParentBundles.add(parentBundles.get(0));
			compBundleTapers.add(bundleTapers.get(0));
			boolean prevTaper = bundleTapers.get(0);
			for (int i=1; i<parentBundles.size(); i++) {
				List<Integer> bundle = parentBundles.get(i);
				boolean taper = compBundleTapers.get(i);
				if (taper && prevTaper) {
					// add it into the previous one
					combParentBundles.get(combParentBundles.size()-1).addAll(bundle);
				} else {
					// start of a new tapered bundle, or both weren't tapered
					prevTaper = taper;
					combParentBundles.add(bundle);
					compBundleTapers.add(taper);
				}
			}
			
			parentBundles = combParentBundles;
			bundleTapers = compBundleTapers;
		}
		
		double[] ret = new double[sectCount];
		int index = 0;
		for (int b=0; b<parentBundles.size(); b++) {
			List<Integer> bundle = parentBundles.get(b);
			boolean taper = bundleTapers.get(b);
			
			if (taper) {
				double[] subAreas = new double[bundle.size()];
				double areaSum = 0;
				for (int s=0; s<bundle.size(); s++) {
					subAreas[s] = rupSet.getAreaForSection(bundle.get(s));
					areaSum += subAreas[s];
				}
				double[] subSlips = TAPERED.calcSlipOnSectionsForRup(subAreas, areaSum, aveSlip);
				Preconditions.checkState(subSlips.length == bundle.size());
				for (int i=0; i<subSlips.length; i++)
					ret[index++] = subSlips[i];
			} else {
				for (int i=0; i<bundle.size(); i++)
					ret[index++] = aveSlip;
			}
		}
		Preconditions.checkState(index == sectCount);
		
		return ret;
	}

}
