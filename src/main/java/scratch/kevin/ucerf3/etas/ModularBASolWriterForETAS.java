package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.CompletableFuture;

import org.opensha.commons.logicTree.BranchWeightProvider;
import org.opensha.commons.util.modules.AverageableModule.AveragingAccumulator;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchAveragingOrder;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchParentSectParticMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchRegionalMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchSectBVals;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchSectNuclMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchSectParticMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultGridAssociations;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.InversionTargetMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupMFDsModule;
import org.opensha.sha.earthquake.faultSysSolution.modules.SubSeismoOnFaultMFDs;
import org.opensha.sha.earthquake.faultSysSolution.util.BranchAverageSolutionCreator;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

import scratch.UCERF3.U3CompoundFaultSystemSolution;
import scratch.UCERF3.U3FaultSystemSolutionFetcher;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.enumTreeBranches.MaxMagOffFault;
import scratch.UCERF3.enumTreeBranches.MomentRateFixes;
import scratch.UCERF3.enumTreeBranches.SpatialSeisPDF;
import scratch.UCERF3.griddedSeismicity.AbstractGridSourceProvider.Precomputed;
import scratch.UCERF3.griddedSeismicity.FaultPolyMgr;
import scratch.UCERF3.griddedSeismicity.UCERF3_GridSourceGenerator;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.inversion.U3InversionTargetMFDs;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;

public class ModularBASolWriterForETAS {

	public static void main(String[] args) throws IOException {
		FaultModels fm = FaultModels.FM3_2;
		SpatialSeisPDF spatSeisPDF = SpatialSeisPDF.UCERF3;
		
		// build BA polys
		File origSolFile = new File("/home/kevin/OpenSHA/UCERF3/rup_sets/modular/"
				+fm.encodeChoiceString()+"_SpatSeisU3_branch_averaged_full_modules.zip");
		FaultSystemSolution origSol = FaultSystemSolution.load(origSolFile);
		FaultSystemRupSet origRupSet = origSol.getRupSet();
		FaultPolyMgr polyMGR = FaultPolyMgr.create(origRupSet.getFaultSectionDataList(), U3InversionTargetMFDs.FAULT_BUFFER);
		
		File inputDir = new File("/home/kevin/OpenSHA/UCERF3/rup_sets/orig");
		File comoundFile = new File(inputDir, "full_ucerf3_compound_sol.zip");
		
		U3FaultSystemSolutionFetcher cfss = U3CompoundFaultSystemSolution.fromZipFile(comoundFile);
		
		int branchCount = 0;
		double totWeight = 0;
		BranchAverageSolutionCreator baBuilder = new BranchAverageSolutionCreator(new BranchWeightProvider.CurrentWeights());
		CompletableFuture<Void> avgFuture = null;
		for (U3LogicTreeBranch branch : cfss.getBranches()) {
			if (!branch.hasValue(spatSeisPDF) || !branch.hasValue(fm))
				continue;
			double weight = branch.getBranchWeight();
			totWeight += weight;
			System.out.println("***********************");
			System.out.println("Processing branch "+branchCount+" with weight="+(float)weight+": "+branch);
			
			InversionFaultSystemSolution sol = cfss.getSolution(branch);
			FaultSystemRupSet rupSet = sol.getRupSet();
			
			// replace poly associations
			rupSet.removeModuleInstances(FaultGridAssociations.class);
			rupSet.addModule(polyMGR);
			
			// build new grid source prov
			UCERF3_GridSourceGenerator gridSources = new UCERF3_GridSourceGenerator(sol, spatSeisPDF,
					MomentRateFixes.NONE,
					rupSet.requireModule(InversionTargetMFDs.class),
					sol.requireModule(SubSeismoOnFaultMFDs.class),
					branch.requireValue(MaxMagOffFault.class).getMaxMagOffFault(),
					polyMGR);
			
			sol.removeModuleInstances(GridSourceProvider.class);
			sol.addModule(gridSources);
			
			if (avgFuture != null)
				avgFuture.join();
			avgFuture = CompletableFuture.runAsync(new Runnable() {
				
				@Override
				public void run() {
					baBuilder.addSolution(sol, branch);
				}
			});
			System.out.println("DONE branch "+branchCount);
			System.out.println("***********************");
			branchCount++;
			
//			if (branchCount == 10)
//				break;
		}
		if (avgFuture != null)
			avgFuture.join();
		System.out.println("Processed "+branchCount+" branches with total weight: "+(float)totWeight);
		
		FaultSystemSolution avgSol = baBuilder.build();
		IncrementalMagFreqDist totMFD = gridSourceTotMFD(avgSol.getGridSourceProvider());
		System.out.println("Total Average Gridded Seis MFD:\n"+totMFD);
		avgSol.removeModuleInstances(RupMFDsModule.class);
		avgSol.removeModuleInstances(BranchAveragingOrder.class);
		avgSol.removeModuleInstances(BranchRegionalMFDs.class);
		avgSol.removeModuleInstances(BranchSectNuclMFDs.class);
		avgSol.removeModuleInstances(BranchSectParticMFDs.class);
		avgSol.removeModuleInstances(BranchParentSectParticMFDs.class);
		avgSol.removeModuleInstances(BranchSectBVals.class);
		File outputFile = new File("/tmp/"+fm.encodeChoiceString()+"_"+spatSeisPDF.encodeChoiceString()+"_sol_ba_modular.zip");
		avgSol.write(outputFile);
//		File origSolFile = new File("/home/kevin/OpenSHA/UCERF3/rup_sets/modular/FM3_1_SpatSeisU3_branch_averaged_full_modules.zip");
//		File outputFile = new File("/tmp/sol_modular.zip");
//		
//		FaultSystemSolution origSol = FaultSystemSolution.load(origSolFile);
//		
//		FaultSystemRupSet rupSet = origSol.getRupSet();
//		
//		// this won't have PolygonFaultAssociations because polygons vary by branch
//		
//		// build BA polys
//		FaultPolyMgr polyMGR = FaultPolyMgr.create(rupSet.getFaultSectionDataList(), U3InversionTargetMFDs.FAULT_BUFFER);
//		
//		// build one for each Mmax
//		AveragingAccumulator<GridSourceProvider> accumulator = null;
//		for (MaxMagOffFault mMaxOff : MaxMagOffFault.values()) {
//			double weight = mMaxOff.getRelativeWeight(InversionModels.CHAR_CONSTRAINED);
//			if (weight == 0d)
//				continue;
//			System.out.println("Building for mMaxOff="+mMaxOff);
//			UCERF3_GridSourceGenerator gridSources = new UCERF3_GridSourceGenerator(origSol, spatSeisPDF,
//					MomentRateFixes.NONE,
//					rupSet.requireModule(InversionTargetMFDs.class),
//					origSol.requireModule(SubSeismoOnFaultMFDs.class),
//					mMaxOff.getMaxMagOffFault(),
//					polyMGR);
//			IncrementalMagFreqDist totMFD = gridSourceTotMFD(gridSources);
//			System.out.println("Total gridded MFD:\n"+totMFD);
//			if (accumulator == null)
//				accumulator = new Precomputed(gridSources).averagingAccumulator();
//			accumulator.process(gridSources, weight);
//		}
//		GridSourceProvider avgGridProv = accumulator.getAverage();
//		System.out.println("Built average grid source prov of type: "+avgGridProv.getClass().getName());
//		IncrementalMagFreqDist totMFD = gridSourceTotMFD(avgGridProv);
//		System.out.println("Total gridded MFD:\n"+totMFD);
//		
//		rupSet.removeModuleInstances(FaultGridAssociations.class);
//		rupSet.addModule(polyMGR);
//		origSol.removeModuleInstances(GridSourceProvider.class);
//		origSol.addModule(avgGridProv);
//		origSol.write(outputFile);
	}

//	public static void main(String[] args) throws IOException {
//		SpatialSeisPDF spatSeisPDF = SpatialSeisPDF.UCERF3;
//		File origSolFile = new File("/home/kevin/OpenSHA/UCERF3/rup_sets/modular/FM3_1_SpatSeisU3_branch_averaged_full_modules.zip");
//		File outputFile = new File("/tmp/sol_modular.zip");
//		
//		FaultSystemSolution origSol = FaultSystemSolution.load(origSolFile);
//		
//		FaultSystemRupSet rupSet = origSol.getRupSet();
//		
//		// this won't have PolygonFaultAssociations because polygons vary by branch
//		
//		// build BA polys
//		FaultPolyMgr polyMGR = FaultPolyMgr.create(rupSet.getFaultSectionDataList(), U3InversionTargetMFDs.FAULT_BUFFER);
//		
//		// build one for each Mmax
//		AveragingAccumulator<GridSourceProvider> accumulator = null;
//		for (MaxMagOffFault mMaxOff : MaxMagOffFault.values()) {
//			double weight = mMaxOff.getRelativeWeight(InversionModels.CHAR_CONSTRAINED);
//			if (weight == 0d)
//				continue;
//			System.out.println("Building for mMaxOff="+mMaxOff);
//			UCERF3_GridSourceGenerator gridSources = new UCERF3_GridSourceGenerator(origSol, spatSeisPDF,
//					MomentRateFixes.NONE,
//					rupSet.requireModule(InversionTargetMFDs.class),
//					origSol.requireModule(SubSeismoOnFaultMFDs.class),
//					mMaxOff.getMaxMagOffFault(),
//					polyMGR);
//			IncrementalMagFreqDist totMFD = gridSourceTotMFD(gridSources);
//			System.out.println("Total gridded MFD:\n"+totMFD);
//			if (accumulator == null)
//				accumulator = new Precomputed(gridSources).averagingAccumulator();
//			accumulator.process(gridSources, weight);
//		}
//		GridSourceProvider avgGridProv = accumulator.getAverage();
//		System.out.println("Built average grid source prov of type: "+avgGridProv.getClass().getName());
//		IncrementalMagFreqDist totMFD = gridSourceTotMFD(avgGridProv);
//		System.out.println("Total gridded MFD:\n"+totMFD);
//		
//		rupSet.removeModuleInstances(FaultGridAssociations.class);
//		rupSet.addModule(polyMGR);
//		origSol.removeModuleInstances(GridSourceProvider.class);
//		origSol.addModule(avgGridProv);
//		origSol.write(outputFile);
//	}
	
	private static IncrementalMagFreqDist gridSourceTotMFD(GridSourceProvider gridSources) {
		IncrementalMagFreqDist totMFD = null;
		for (int i=0; i<gridSources.size(); i++) {
			IncrementalMagFreqDist mfd = gridSources.getMFD(i);
			if (mfd == null)
				continue;
			if (totMFD == null) {
				totMFD = new IncrementalMagFreqDist(mfd.getMinX(), mfd.size(), mfd.getDelta());
			} else {
				Preconditions.checkState(mfd.getMinX() == totMFD.getMinX());
				Preconditions.checkState(mfd.getDelta() == totMFD.getDelta());
				if (mfd.size() > totMFD.size()) {
					// grow it
					IncrementalMagFreqDist tmp = new IncrementalMagFreqDist(mfd.getMinX(), mfd.size(), mfd.getDelta());
					for (int j=0; j<totMFD.size(); j++)
						totMFD.set(j, totMFD.getY(j));
				}
			}
			for (int j=0; j<mfd.size(); j++)
				totMFD.add(j, mfd.getY(j));
		}
		return totMFD;
	}

}
