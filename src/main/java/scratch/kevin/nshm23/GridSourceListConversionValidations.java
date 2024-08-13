package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CompletableFuture;
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
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.gridded.NSHM23_SingleRegionGridSourceProvider;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.faultSysSolution.modules.MFDGridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

public class GridSourceListConversionValidations {

	public static void main(String[] args) throws IOException {
//		FaultSystemSolution sol1 = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
//				+ "2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged_gridded.zip"));
//		FaultSystemSolution sol2 = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
//				+ "2024_02_02-nshm23_branches-WUS_FM_v3-gridded_rebuild/results_WUS_FM_v3_branch_averaged_gridded.zip"));
//		
//		MFDGridSourceProvider mfdProv = sol1.requireModule(MFDGridSourceProvider.class);
//		FaultGridAssociations assoc = sol1.getRupSet().requireModule(FaultGridAssociations.class);
//		GridSourceList sourceList = sol2.requireModule(GridSourceList.class);
//		
//		validate(mfdProv, assoc, sourceList);
		
		ModuleContainer.VERBOSE_DEFAULT = false;
		
		SolutionLogicTree slt = SolutionLogicTree.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_02_02-nshm23_branches-WUS_FM_v3/results.zip"));
		LogicTree<?> faultTree = slt.getLogicTree();
		LogicTree<?> gridTree = LogicTree.buildExhaustive(NSHM23_LogicTreeBranch.levelsOffFault, true);
		int numTestBranches = 2;
		AveragingAccumulator<GridSourceProvider> mfdAvg = null;
		AveragingAccumulator<GridSourceProvider> listAvg = null;
		AveragingAccumulator<FaultGridAssociations> assocAvg = null;
		boolean convertToPrecomputed = true;
		boolean averageWithPrev = true;
		for (int i=0; i<numTestBranches; i++) {
			LogicTreeBranch<?> faultBranch = faultTree.getBranch(i);
			double faultWeight = faultBranch.getBranchWeight();
			System.out.println("Processing branch "+i+": "+faultBranch);
			FaultSystemSolution sol = slt.forBranch(faultBranch);
			
			new NSHM23_InvConfigFactory().preGridBuildHook(sol, faultBranch);
			
			List<CompletableFuture<GridSourceProvider[]>> futures = new ArrayList<>(gridTree.size());
			for (LogicTreeBranch<?> gridBranch : gridTree) {
				futures.add(CompletableFuture.supplyAsync(new Supplier<GridSourceProvider[]>() {

					@Override
					public GridSourceProvider[] get() {
						try {
							NSHM23_SingleRegionGridSourceProvider mfdProv = (NSHM23_SingleRegionGridSourceProvider)
									NSHM23_InvConfigFactory.buildGridSourceProv(sol, gridBranch);
							GridSourceList gridList = mfdProv.convertToGridSourceList(0d);
							if (convertToPrecomputed && gridList instanceof GridSourceList.DynamicallyBuilt) {
								// do the conversion
								EnumMap<TectonicRegionType, List<List<GriddedRupture>>> trtRupLists = new EnumMap<>(TectonicRegionType.class);
								for (TectonicRegionType trt : gridList.getTectonicRegionTypes()) {
									List<List<GriddedRupture>> rupLists = new ArrayList<>();
									trtRupLists.put(trt, rupLists);
									for (int gridIndex=0; gridIndex<gridList.getNumLocations(); gridIndex++)
										rupLists.add(gridList.getRuptures(trt, gridIndex));
								}
								gridList = new GridSourceList.Precomputed(gridList.getGriddedRegion(), trtRupLists);
							}
							return new GridSourceProvider[] {mfdProv, gridList};
						} catch (IOException e) {
							e.printStackTrace();
							System.exit(1);
							return null;
						}
					}
				}));
			}
			
			AveragingAccumulator<GridSourceProvider> mfdSubAvg = null;
			AveragingAccumulator<GridSourceProvider> listSubAvg = null;
			FaultGridAssociations assoc = sol.getRupSet().requireModule(FaultGridAssociations.class);
			
			AveragingAccumulator<GridSourceProvider> mfdPrevAvg = null;
			AveragingAccumulator<GridSourceProvider> listPrevAvg = null;
			
			for (int b=0; b<gridTree.size(); b++) {
				LogicTreeBranch<?> gridBranch = gridTree.getBranch(b);
				double gridWeight = gridBranch.getBranchWeight();
				GridSourceProvider[] provs = futures.get(b).join();
				MFDGridSourceProvider mfdProv = (MFDGridSourceProvider)provs[0];
				GridSourceList gridList = (GridSourceList)provs[1];
				
				System.out.println("Built for "+b+"/"+gridTree.size()+" "+gridBranch+", validating...");
				validate(mfdProv, assoc, gridList, false);
				
				if (mfdSubAvg == null)
					mfdSubAvg = mfdProv.averagingAccumulator();
				mfdSubAvg.process(mfdProv, gridWeight);
				if (listSubAvg == null)
					listSubAvg = gridList.averagingAccumulator();
				listSubAvg.process(gridList, gridWeight);
				
				if (averageWithPrev) {
					if (b > 0) {
						System.out.println("Averaging with prior and comparing");
						mfdPrevAvg.process(mfdProv, gridWeight);
						MFDGridSourceProvider mfdWithPrev = (MFDGridSourceProvider)mfdPrevAvg.getAverage();
						listPrevAvg.process(gridList, gridWeight);
						GridSourceList listWithPrev = (GridSourceList)listPrevAvg.getAverage();
						validate(mfdWithPrev, assoc, listWithPrev, false);
					}
					
					if (b < gridTree.size()-1) {
						mfdPrevAvg = mfdProv.averagingAccumulator();
						mfdPrevAvg.process(mfdProv, gridWeight);
						listPrevAvg = gridList.averagingAccumulator();
						listPrevAvg.process(gridList, gridWeight);
					}
				}
			}
			
			if (assocAvg == null)
				assocAvg = assoc.averagingAccumulator();
			assocAvg.process(assoc, faultWeight);
			GridSourceProvider avgOfMFD = mfdSubAvg.getAverage();
			GridSourceProvider avgOfList = listSubAvg.getAverage();
			System.out.println("Built all for "+faultBranch+", validating averages");
			validate((MFDGridSourceProvider)avgOfMFD, assoc, (GridSourceList)avgOfList, false);
			if (mfdAvg == null)
				mfdAvg = avgOfMFD.averagingAccumulator();
			mfdAvg.process(avgOfMFD, faultWeight);
			if (listAvg == null)
				listAvg = avgOfList.averagingAccumulator();
			listAvg.process(avgOfList, faultWeight);
		}
		
		System.out.println("Built all, validating averages");
		FaultGridAssociations avgAssoc = assocAvg.getAverage();
		GridSourceProvider avgOfMFD = mfdAvg.getAverage();
		GridSourceProvider avgOfList = listAvg.getAverage();
		validate((MFDGridSourceProvider)avgOfMFD, avgAssoc, (GridSourceList)avgOfList, true);
	}
	
	private static void validate(MFDGridSourceProvider mfdProv, FaultGridAssociations assoc, GridSourceList sourceList, boolean verbose) {
		// filter the MFD provider to be the same minimum magnitude
		double sourceListMinMag = Double.POSITIVE_INFINITY;
		for (int i=0; i<mfdProv.getNumLocations(); i++)
			for (GriddedRupture rup : sourceList.getRuptures(mfdProv.getTectonicRegionType(i), i))
				sourceListMinMag = Math.min(sourceListMinMag, rup.properties.magnitude);
		mfdProv = (MFDGridSourceProvider) mfdProv.getAboveMinMag((float)sourceListMinMag);
		
		for (int i=0; i<mfdProv.getNumLocations(); i++) {
			TectonicRegionType trt = mfdProv.getTectonicRegionType(i);
			IncrementalMagFreqDist mfdSubSeis1 = mfdProv.getMFD_SubSeisOnFault(i);
			IncrementalMagFreqDist mfdSubSeis2 = sourceList.getMFD_SubSeisOnFault(trt, i);
			IncrementalMagFreqDist mfdUnassoc1 = mfdProv.getMFD_Unassociated(i);
			IncrementalMagFreqDist mfdUnassoc2 = sourceList.getMFD_Unassociated(trt, i);
			
			try {
				assertEquals(mfdUnassoc1, mfdUnassoc2, i, false);
				assertEquals(mfdSubSeis1, mfdSubSeis2, i, true);
				validateAssoc(mfdSubSeis1, mfdUnassoc1, assoc, i);
				validateAssoc(mfdSubSeis2, mfdUnassoc2, assoc, i);
				if (verbose) {
					if (mfdSubSeis1 == null && mfdSubSeis2 == null && mfdUnassoc1 == null && mfdUnassoc2 == null)
						System.out.println("All null for "+i);
					else
						System.out.println("Validated "+i+"/"+mfdProv.getNumLocations()
								+"; assoc="+(mfdSubSeis1 != null && mfdSubSeis1.calcSumOfY_Vals() > 0d || mfdSubSeis2 != null && mfdSubSeis2.calcSumOfY_Vals() > 0d)
								+", unassoc="+(mfdUnassoc1 != null && mfdUnassoc1.calcSumOfY_Vals() > 0d || mfdUnassoc2 != null && mfdUnassoc2.calcSumOfY_Vals() > 0d));
				}
			} catch (Exception e) {
				System.out.println("Failed for Grid Node "+i+"/"+mfdProv.getNumLocations());
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
			
		}
		System.out.println("Validated!");
	}
	
	private static final double RELATIVE_TOL = 1e-3;
	
	private static void assertEquals(IncrementalMagFreqDist mfd1, IncrementalMagFreqDist mfd2, int gridIndex, boolean assoc) {
		if (mfd1 == null) {
			Preconditions.checkState(mfd2 == null || mfd2.calcSumOfY_Vals() == 0d,
					"GridIndex=%s, assoc=%s: MFDList is null, SourceList is %s", gridIndex, assoc, mfd2);
		} else {
			for (int i=0; i<mfd1.size(); i++) {
				double x1 = mfd1.getX(i);
				double y1 = mfd1.getY(i);
				int i2 = mfd2.getClosestXIndex(x1);
				double x2 = mfd2.getX(i2);
				double y2;
				if ((float)Math.abs(x2 - x1) >= (float)0.5*mfd1.getDelta()) {
					Preconditions.checkState(x1 > x2);
					x2 = x1;
					y2 = 0d;
				} else {
					y2 = mfd2.getY(i2);
				}
				Preconditions.checkState(Precision.equalsWithRelativeTolerance(y1, y2, RELATIVE_TOL) || y1 < 1e-10 && y2 < 1e-10,
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
				Preconditions.checkState(Precision.equalsWithRelativeTolerance(y1, y2, RELATIVE_TOL) || y1 < 1e-10 && y2 < 1e-10,
						"Mismatch at GridIndex=%s, assoc=%s:\tMFDProv[%s]=%s\tSourceList[%s]=%s",
						gridIndex, assoc, (float)x2, (float)y2, (float)x1, (float)y1);
			}
		}
	}
	
	private static void validateAssoc(IncrementalMagFreqDist subSeisMFD, IncrementalMagFreqDist unAssocMFD, FaultGridAssociations assoc, int nodeIndex) {
		Map<Integer, Double> scaledSectFracts = assoc.getScaledSectFracsOnNode(nodeIndex);
		double fractAssociated = 0d;
		if (scaledSectFracts != null)
			for (double fract : scaledSectFracts.values())
				fractAssociated += fract;
		if (fractAssociated == 0d) {
			Preconditions.checkState(subSeisMFD == null || subSeisMFD.calcSumOfY_Vals() == 0d,
					"FaultGridAssociations says fractAssociated=%s, but we have a SubSeisMFD: %s", fractAssociated, subSeisMFD);
		} else if (fractAssociated == 1d) {
			Preconditions.checkState(unAssocMFD == null || unAssocMFD.calcSumOfY_Vals() == 0d,
					"FaultGridAssociations says fractAssociated=%s, but we have an UnAssocMFD: %s", fractAssociated, subSeisMFD);
		} else {
			// these don't match, but not sure they actually are supposed to...(they don't for the original MFDGridSourceProvider when averaged)
//			Preconditions.checkNotNull(unAssocMFD, "FaultGridAssociations says fractAssociated=%s, but we don't have an UnAssocMFD!", fractAssociated);
//			Preconditions.checkNotNull(subSeisMFD, "FaultGridAssociations says fractAssociated=%s, but we don't have a SubSeisMFD!", subSeisMFD);
//			double sumUnassoc = unAssocMFD.calcSumOfY_Vals();
//			double sumSubSeis = subSeisMFD.calcSumOfY_Vals();
//			double calcSumFract = sumSubSeis/(sumSubSeis+sumUnassoc);
//			double mminUnassoc = unAssocMFD.getY(0);
//			double mminSubSeis = subSeisMFD.getY(0);
//			double calcMinMagFract = mminSubSeis/(mminSubSeis+mminUnassoc);
//			Preconditions.checkState(Precision.equalsWithRelativeTolerance(calcSumFract, fractAssociated, 1e-3)
//					|| calcSumFract < 1e-10 && fractAssociated < 1e-10,
//					"FaultGridAssociations says fractAssoc=%s, MFD sumRatio=%s, mminRatio=%s",
//					(float)fractAssociated, (float)calcSumFract, (float)calcMinMagFract);
		}
	}

}
