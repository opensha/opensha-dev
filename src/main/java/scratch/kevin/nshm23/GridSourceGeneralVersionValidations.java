package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.EnumMap;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.function.Supplier;

import org.apache.commons.math3.util.Precision;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.util.modules.AverageableModule.AveragingAccumulator;
import org.opensha.commons.util.modules.ModuleContainer;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultGridAssociations;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

public class GridSourceGeneralVersionValidations {

	public static void main(String[] args) throws IOException {
		ModuleContainer.VERBOSE_DEFAULT = false;
		
		SolutionLogicTree slt = SolutionLogicTree.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_02_02-nshm23_branches-WUS_FM_v3/results.zip"));
		LogicTree<?> faultTree = slt.getLogicTree();
		LogicTree<?> gridTree = LogicTree.buildExhaustive(NSHM23_LogicTreeBranch.levelsOffFault, true);
		int numTestBranches = 2;
		AveragingAccumulator<GridSourceProvider> origListAvg = null;
		AveragingAccumulator<GridSourceProvider> newListAvg = null;
		boolean convertToPrecomputed = true;
		boolean averageWithPrev = true;
		
//		final int threads = 1;
		final int threads = 16;
		ExecutorService exec = threads > 1 ? Executors.newFixedThreadPool(threads) : Executors.newSingleThreadExecutor();
		
		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
		for (int i=0; i<numTestBranches; i++) {
			LogicTreeBranch<?> faultBranch = faultTree.getBranch(i);
			double faultWeight = faultBranch.getBranchWeight();
			System.out.println("Processing branch "+i+": "+faultBranch);
			FaultSystemSolution sol = slt.forBranch(faultBranch);
			
			factory.preGridBuildHook(sol, faultBranch);
			
			final int branchIndex = i;
			
			List<CompletableFuture<GridSourceList[]>> futures = new ArrayList<>(gridTree.size());
			for (int g=0; g<gridTree.size(); g++) {
				final int gridBranchIndex = g;
				LogicTreeBranch<?> gridBranch = gridTree.getBranch(g);
				Supplier<GridSourceList[]> supplier = () -> {
					try {
						if (threads < 2)
							System.out.println("Building original for "+branchIndex+"-"+gridBranchIndex+": "+gridBranch);
						GridSourceList origGridList = factory.buildGridSourceProvider(sol, gridBranch);
						if (threads < 2)
							System.out.println("Building updated for "+branchIndex+"-"+gridBranchIndex+": "+gridBranch);
						GridSourceList modGridList = NSHM23_InvConfigFactory.buildUpdatedGridSourceProv(sol, gridBranch);
						System.out.println("Done building for "+branchIndex+"-"+gridBranchIndex+": "+gridBranch);
						if (convertToPrecomputed && origGridList instanceof GridSourceList.DynamicallyBuilt) {
							// do the conversion
							EnumMap<TectonicRegionType, List<List<GriddedRupture>>> trtRupLists = new EnumMap<>(TectonicRegionType.class);
							for (TectonicRegionType trt : origGridList.getTectonicRegionTypes()) {
								List<List<GriddedRupture>> rupLists = new ArrayList<>();
								trtRupLists.put(trt, rupLists);
								for (int gridIndex=0; gridIndex<origGridList.getNumLocations(); gridIndex++)
									rupLists.add(origGridList.getRuptures(trt, gridIndex));
							}
							origGridList = new GridSourceList.Precomputed(origGridList.getGriddedRegion(), trtRupLists);
						}
						return new GridSourceList[] {origGridList, modGridList};
					} catch (IOException e) {
						e.printStackTrace();
						System.exit(1);
						return null;
					}
				};
				futures.add(CompletableFuture.supplyAsync(supplier, exec));
			}
			
			AveragingAccumulator<GridSourceProvider> origSubAvg = null;
			AveragingAccumulator<GridSourceProvider> modSubAvg = null;
			
			AveragingAccumulator<GridSourceProvider> origPrevAvg = null;
			AveragingAccumulator<GridSourceProvider> modPrevAvg = null;
			
			for (int b=0; b<gridTree.size(); b++) {
				LogicTreeBranch<?> gridBranch = gridTree.getBranch(b);
				double gridWeight = gridBranch.getBranchWeight();
				GridSourceList[] provs = futures.get(b).join();
				futures.set(b, null);
				GridSourceList origProv = provs[0];
				GridSourceList modList = provs[1];
				
				System.out.println("Built for "+b+"/"+gridTree.size()+" "+gridBranch+", validating...");
				validate(origProv, modList, false);
				
				if (origSubAvg == null)
					origSubAvg = origProv.averagingAccumulator();
				origSubAvg.process(origProv, gridWeight);
				if (modSubAvg == null)
					modSubAvg = modList.averagingAccumulator();
				modSubAvg.process(modList, gridWeight);
				
				if (averageWithPrev) {
					if (b > 0) {
						System.out.println("Averaging with prior and comparing");
						origPrevAvg.process(origProv, gridWeight);
						GridSourceList mfdWithPrev = (GridSourceList)origPrevAvg.getAverage();
						modPrevAvg.process(modList, gridWeight);
						GridSourceList listWithPrev = (GridSourceList)modPrevAvg.getAverage();
						validate(mfdWithPrev, listWithPrev, false);
					}
					
					if (b < gridTree.size()-1) {
						origPrevAvg = origProv.averagingAccumulator();
						origPrevAvg.process(origProv, gridWeight);
						modPrevAvg = modList.averagingAccumulator();
						modPrevAvg.process(modList, gridWeight);
					}
				}
			}
			
			GridSourceList avgOfOrig = (GridSourceList)origSubAvg.getAverage();
			GridSourceList avgOfMod = (GridSourceList)modSubAvg.getAverage();
			System.out.println("Built all for "+faultBranch+", validating averages");
			validate(avgOfOrig, avgOfMod, false);
			if (origListAvg == null)
				origListAvg = avgOfOrig.averagingAccumulator();
			origListAvg.process(avgOfOrig, faultWeight);
			if (newListAvg == null)
				newListAvg = avgOfMod.averagingAccumulator();
			newListAvg.process(avgOfMod, faultWeight);
		}
		
		exec.shutdown();
		
		System.out.println("Built all, validating averages");
		GridSourceList avgOfOrig = (GridSourceList)origListAvg.getAverage();
		GridSourceList avgOfMod = (GridSourceList)newListAvg.getAverage();
		validate(avgOfOrig, avgOfMod, true);
		
		System.out.println("DONE");
		System.exit(0);
	}
	
	private static void validate(GridSourceList origList, GridSourceList modList, boolean verbose) {
		for (int i=0; i<origList.getNumLocations(); i++) {
			IncrementalMagFreqDist mfdSubSeis1 = origList.getMFD_SubSeisOnFault(i);
			IncrementalMagFreqDist mfdSubSeis2 = modList.getMFD_SubSeisOnFault(null, i);
			IncrementalMagFreqDist mfdUnassoc1 = origList.getMFD_Unassociated(i);
			IncrementalMagFreqDist mfdUnassoc2 = modList.getMFD_Unassociated(null, i);
			
			try {
				assertEquals(mfdUnassoc1, mfdUnassoc2, i, false);
				assertEquals(mfdSubSeis1, mfdSubSeis2, i, true);
				if (verbose) {
					if (mfdSubSeis1 == null && mfdSubSeis2 == null && mfdUnassoc1 == null && mfdUnassoc2 == null)
						System.out.println("All null for "+i);
					else
						System.out.println("Validated "+i+"/"+modList.getNumLocations()
								+"; assoc="+(mfdSubSeis1 != null && mfdSubSeis1.calcSumOfY_Vals() > 0d || mfdSubSeis2 != null && mfdSubSeis2.calcSumOfY_Vals() > 0d)
								+", unassoc="+(mfdUnassoc1 != null && mfdUnassoc1.calcSumOfY_Vals() > 0d || mfdUnassoc2 != null && mfdUnassoc2.calcSumOfY_Vals() > 0d));
				}
			} catch (Exception e) {
				System.out.println("Failed for Grid Node "+i+"/"+modList.getNumLocations());
				IncrementalMagFreqDist refMFD = mfdUnassoc1;
				if (refMFD == null)
					refMFD = mfdSubSeis1;
				else if (mfdSubSeis1 != null && mfdSubSeis1.size() > refMFD.size())
					refMFD = mfdSubSeis1;
				if (refMFD != null) {
					System.out.println("MFDs\tUnAssoc1\tUnAssoc2\tSubSeis1\tSubSeis2\tSumEach1\tSumEach2");
					for (int x=0; x<refMFD.size(); x++) {
						double unassoc1 = 0d;
						if (mfdUnassoc1 != null && mfdUnassoc1.size() > x)
							unassoc1 = mfdUnassoc1.getY(x);
						double unassoc2 = 0d;
						if (mfdUnassoc2 != null && mfdUnassoc2.size() > x)
							unassoc2 = mfdUnassoc2.getY(x);
						double subSeis1 = 0d;
						if (mfdSubSeis1 != null && mfdSubSeis1.size() > x)
							subSeis1 = mfdSubSeis1.getY(x);
						double subSeis2 = 0d;
						if (mfdSubSeis2 != null && mfdSubSeis2.size() > x)
							subSeis2 = mfdSubSeis2.getY(x);
						
						System.out.println((float)refMFD.getX(x)+"\t"+(float)unassoc1+"\t"+(float)unassoc2
								+"\t"+(float)subSeis1+"\t"+(float)subSeis2
								+"\t"+(float)(unassoc1+subSeis1)+"\t"+(float)(unassoc2+subSeis2));
					}
				}
				System.out.flush();
				e.printStackTrace();
				System.exit(1);
			}

			// now verify the ruptures themselves
			try {
				assertEquals(origList, modList, i);
			} catch (Exception e) {
				System.out.println("Failed for Grid Node "+i+"/"+modList.getNumLocations());
				System.out.flush();
				e.printStackTrace();
				System.exit(1);
			}
		}
		System.out.println("Validated!");
	}
	
	private static List<GriddedRupture> getSorted(List<GriddedRupture> rups) {
		if (rups == null)
			return new ArrayList<>();
		List<GriddedRupture> ret = new ArrayList<>(rups);
		ret.sort(gridRupComp);
		for (int i=ret.size(); --i>=0;)
			if (ret.get(i).properties.magnitude < 2.5d)
				ret.remove(i);
		return ret;
	}
	
	private static Comparator<GriddedRupture> gridRupComp = new Comparator<GridSourceList.GriddedRupture>() {
		
		@Override
		public int compare(GriddedRupture o1, GriddedRupture o2) {
			// magnitude first
			int cmp = Float.compare((float)o1.properties.magnitude, (float)o2.properties.magnitude);
			if (cmp != 0)
				return cmp;
			// rake
			cmp = Float.compare((float)o1.properties.rake, (float)o2.properties.rake);
			if (cmp != 0)
				return cmp;
			return o1.compareTo(o2);
		}
	};
	
	private static boolean mfdEquals(double v1, double v2) {
		return equalsWithTol(v1, v2, 1e-3, 0d, 1e-10);
	}
	
	private static boolean rateEquals(double v1, double v2) {
		return equalsWithTol(v1, v2, 1e-3, 1e-10, 1e-15);
	}
	
	private static boolean assocEquals(double v1, double v2) {
		return equalsWithTol(v1, v2, 1e-3, 1e-6, 1e-10);
	}
	
	private static boolean equalsWithTol(double v1, double v2, double relativeTol, double absTol,
			double ignoreBelowThresh) {
		if (v1 < ignoreBelowThresh && v2 < ignoreBelowThresh)
			return true;
		if (relativeTol > 0 && !Precision.equalsWithRelativeTolerance(v1, v2, relativeTol))
			return false;
		if (absTol > 0 && !Precision.equals(v1, v2, absTol))
			return false;
		return true;
	}
	
	private static void assertEquals(IncrementalMagFreqDist mfd1, IncrementalMagFreqDist mfd2, int gridIndex, boolean assoc) {
		if (mfd1 == null || mfd1.calcSumOfY_Vals() == 0d) {
			Preconditions.checkState(mfd2 == null || mfd2.calcSumOfY_Vals() == 0d,
					"GridIndex=%s, assoc=%s: MFDList is null, SourceList is %s", gridIndex, assoc, mfd2);
		} else {
			double delta = Math.max(mfd1.getDelta(), mfd2.getDelta());
			
			for (int i=0; i<mfd1.size(); i++) {
				double x1 = mfd1.getX(i);
				double y1 = mfd1.getY(i);
				int i2 = mfd2.getClosestXIndex(x1);
				double x2 = mfd2.getX(i2);
				double y2;
				if ((float)Math.abs(x2 - x1) >= (float)0.5*delta) {
					Preconditions.checkState(x1 > x2, "x1=%s, x2=%s, delta=%s", x1, x2, delta);
					x2 = x1;
					y2 = 0d;
				} else {
					y2 = mfd2.getY(i2);
				}
				Preconditions.checkState(mfdEquals(y1, y2) || y1 < 1e-10 && y2 < 1e-10,
						"Mismatch at GridIndex=%s, assoc=%s:\tMFDProv[%s]=%s\tSourceList[%s]=%s",
						gridIndex, assoc, (float)x1, (float)y1, (float)x2, (float)y2);
			}
			for (int i=mfd1.size(); i<mfd2.size(); i++) {
				double x1 = mfd2.getX(i);
				double y1 = mfd2.getY(i);
				int i2 = mfd1.getClosestXIndex(x1);
				double x2 = mfd1.getX(i2);
				double y2;
				if ((float)Math.abs(x2 - x1) >= (float)0.5*mfd2.getDelta()) {
					Preconditions.checkState(x1 > x2);
					x2 = x1;
					y2 = 0d;
				} else {
					y2 = mfd1.getY(i2);
				}
				Preconditions.checkState(mfdEquals(y1, y2) || y1 < 1e-10 && y2 < 1e-10,
						"Mismatch at GridIndex=%s, assoc=%s:\tMFDProv[%s]=%s\tSourceList[%s]=%s",
						gridIndex, assoc, (float)x2, (float)y2, (float)x1, (float)y1);
			}
		}
	}
	
	private static void assertEquals(GridSourceList gridList1, GridSourceList gridList2, int gridIndex) {
		List<GriddedRupture> origRups = getSorted(gridList1.getRuptures(null, gridIndex));
		List<GriddedRupture> modRups = getSorted(gridList2.getRuptures(null, gridIndex));
		
		double totRate1 = origRups.stream().mapToDouble(R->R.rate).sum();
		double totRate2 = modRups.stream().mapToDouble(R->R.rate).sum();
		
		Preconditions.checkState(origRups.size() == modRups.size(),
				"Original list has %s rups (totRate=%s), mod has %s (totRate=%s); gridIndex=%s",
				origRups.size(), totRate1, modRups.size(), totRate2, gridIndex);
		
		for (int i=0; i<origRups.size(); i++) {
			GriddedRupture rup1 = origRups.get(i);
			GriddedRupture rup2 = modRups.get(i);
			Preconditions.checkState(rup1.properties.equals(rup2.properties),
					"Props mismatch for gridIndex=%s, rupIndex=%s:\n\t%s\n\t%s\n\tTotal rates: orig=%s, mod=%s, delta=%s",
					gridIndex, i, rup1.properties, rup2.properties, totRate1, totRate2, totRate2-totRate1);
			Preconditions.checkState(rateEquals(rup1.rate, rup2.rate),
					"Rate mismatch for gridIndex=%s, rupIndex=%s: %s != %s\n\t%s\n\t%s\n\tTotal rates: orig=%s, mod=%s, delta=%s",
					gridIndex, i, (float)rup1.rate, (float)rup2.rate, rup1.properties, rup2.properties, totRate1, totRate2, totRate2-totRate1);
			int[] assoc1 = rup1.associatedSections;
			int[] assoc2 = rup2.associatedSections;
			if (assoc1 == null || assoc1.length == 0) {
				Preconditions.checkState(assoc2 == null || assoc2.length == 0,
						"Original is unassociated, mod is associated. gridIndex=%s, rupIndex=%s\n\tprops=%s",
						gridIndex, i, rup1.properties);
			} else {
				Preconditions.checkState(assoc1.length == assoc2.length,
						"Association count mismatch for gridIndex=%s, rupIndex=%s, count1=%s, count2=%s, \n\tprops=%s",
						gridIndex, i, rup1.properties);
				for (int i1=0; i1<assoc1.length; i1++) {
					int sectID = assoc1[i1];
					double fract = rup1.associatedSectionFracts[i1];
					boolean found = false;
					for (int i2=0; i2<assoc2.length; i2++) {
						if (sectID == assoc2[i2]) {
							found = true;
							double fract2 = rup2.associatedSectionFracts[i2];
//							Preconditions.checkState((float)fract == (float)fract2,
							Preconditions.checkState(assocEquals(fract, fract2),
									"Association fraction mismatch for gridIndex=%s, rupIndex=%s, sect=%s: %s != %s",
									gridIndex, i, sectID, (float)fract, (float)fract2);
						}
					}
					Preconditions.checkState(found, "No association found for %s in updated gridIndex=%s, rupIndex=%s",
							sectID, gridIndex, i);
				}
			}
		}
	}

}
