package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CompletableFuture;

import org.apache.commons.lang3.exception.ExceptionUtils;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.util.modules.AverageableModule.AveragingAccumulator;
import org.opensha.commons.util.modules.OpenSHA_Module;
import org.opensha.sha.earthquake.PointSource;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultGridAssociations;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.MFDGridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.ModelRegion;
import org.opensha.sha.earthquake.faultSysSolution.modules.RegionsOfInterest;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.faultSysSolution.util.MaxMagOffFaultBranchNode;
import org.opensha.sha.earthquake.faultSysSolution.util.SolModuleStripper;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.gridded.NSHM23_GridFocalMechs;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.gridded.NSHM23_SingleRegionGridSourceProvider.NSHM23_WUS_FiniteRuptureConverter;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeclusteringAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_RegionalSeismicity;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SeisSmoothingAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.SeismicityRegions;
import org.opensha.sha.earthquake.util.GriddedSeismicitySettings;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

public class WriteGridWithoutCarveout {

	public static void main(String[] args) throws IOException {
		FaultSystemSolution baSol = FaultSystemSolution.load(
				new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
						+ "2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged_gridded.zip"));
		final FaultSystemSolution origSol = baSol;
		
		File outputDir = new File("/home/kevin/markdown/inversions/2026_01_15-nshm23-ba-with-no-fault-grid-prov");
		
		LogicTree<?> gridTree = LogicTree.buildExhaustive(NSHM23_LogicTreeBranch.levelsOffFault, true);
		
		SeismicityRegions region = SeismicityRegions.CONUS_WEST;
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(
				Math.max(7.95, baSol.getRupSet().getMaxMag()));
		
		FaultGridAssociations assoc = baSol.getRupSet().requireModule(FaultGridAssociations.class);
		GriddedRegion gridReg = baSol.getGridSourceProvider().getGriddedRegion();
		
		EnumMap<TectonicRegionType, Region> trtRegions = NSHM23_InvConfigFactory.getTRT_Regions();
		TectonicRegionType[] trts = new TectonicRegionType[gridReg.getNodeCount()];
		for (int i=0; i<trts.length; i++)
			trts[i] = TectonicRegionType.ACTIVE_SHALLOW;
		for (int i=0; i<trts.length; i++) {
			Location loc = gridReg.getLocation(i);
			for (TectonicRegionType trt : trtRegions.keySet()) {
				if (trtRegions.get(trt).contains(loc)) {
					trts[i] = trt;
					break;
				}
			}
		}
		
		// focal mechanisms
		double[] fractStrikeSlip = NSHM23_GridFocalMechs.getFractStrikeSlip(region, gridReg);
		double[] fractReverse = NSHM23_GridFocalMechs.getFractReverse(region, gridReg);
		double[] fractNormal = NSHM23_GridFocalMechs.getFractNormal(region, gridReg);
		
		List<CompletableFuture<GridSourceList>> treeFutures = new ArrayList<>();
		
		for (int b=0; b<gridTree.size(); b++) {
			LogicTreeBranch<?> branch = gridTree.getBranch(b);
//			System.out.println("Processing branch "+b+"/"+gridTree.size()+": "+branch);
			
			NSHM23_RegionalSeismicity seisBranch = branch.requireValue(NSHM23_RegionalSeismicity.class);
			NSHM23_DeclusteringAlgorithms declusteringAlg = branch.requireValue(NSHM23_DeclusteringAlgorithms.class);
			NSHM23_SeisSmoothingAlgorithms seisSmooth = branch.requireValue(NSHM23_SeisSmoothingAlgorithms.class);
			double maxMagOff = branch.requireValue(MaxMagOffFaultBranchNode.class).getMaxMagOffFault();
			
			// total G-R up to Mmax
			IncrementalMagFreqDist totalGR = seisBranch.build(region, refMFD, maxMagOff);
			
			// spatial seismicity PDF
			double[] pdf = seisSmooth.load(region, declusteringAlg);
			
			treeFutures.add(CompletableFuture.supplyAsync(() -> {
				NoFaultGridProv gridProv = new NoFaultGridProv(gridReg, trts, totalGR, pdf, fractStrikeSlip, fractReverse, fractNormal);
				return GridSourceList.convert(gridProv, assoc, new NSHM23_WUS_FiniteRuptureConverter());
			}));
		}
		
		GridSourceProvider avgProv = buildAvgGridProv(gridTree, treeFutures);
		
		List<OpenSHA_Module> rsModulesToKeep = new ArrayList<>();
		rsModulesToKeep.add(baSol.getRupSet().getModule(RegionsOfInterest.class));
		rsModulesToKeep.add(baSol.getRupSet().getModule(ModelRegion.class));
		rsModulesToKeep.add(new FaultGridAssociations.Precomputed(assoc));
		System.out.println("Simplifying solution");
		baSol = SolModuleStripper.stripModules(baSol, 0d, true, true);
		for (OpenSHA_Module module : rsModulesToKeep)
			baSol.getRupSet().addModule(module);
		
		baSol.setGridSourceProvider(avgProv);
		
		System.out.println("Writing solution");
		baSol.write(new File(outputDir, "solution_gridded_pure_gr.zip"));
		
		System.out.println("Now loading and rebuilding with original grid-branch SLT");
		
		SolutionLogicTree slt = SolutionLogicTree.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_02_02-nshm23_branches-WUS_FM_v3/results_gridded_branches.zip"));
		
		List<LogicTreeLevel<? extends LogicTreeNode>> gridFileLevels = new ArrayList<>();
		gridFileLevels.add(NSHM23_LogicTreeBranch.FM);
		gridFileLevels.addAll(NSHM23_LogicTreeBranch.levelsOffFault);
		treeFutures = new ArrayList<>();
		// first load the original ones
		for (int b=0; b<gridTree.size(); b++) {
			LogicTreeBranch<?> branch = gridTree.getBranch(b);
//			System.out.println("Processing branch "+b+"/"+gridTree.size()+": "+branch);
			
			LogicTreeBranch<LogicTreeNode> gridFileBranch = new LogicTreeBranch<>(gridFileLevels);
			gridFileBranch.setValue(NSHM23_FaultModels.WUS_FM_v3);
			for (LogicTreeNode node : branch)
				gridFileBranch.setValue(node);
			
			treeFutures.add(CompletableFuture.supplyAsync(() -> {
				try {
					return (GridSourceList)slt.loadGridProvForBranch(gridFileBranch);
				} catch (IOException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
			}));
		}
		
		avgProv = buildAvgGridProv(gridTree, treeFutures);
		List<GridSourceList> origGridMFDs = new ArrayList<>();
		for (CompletableFuture<GridSourceList> future : treeFutures)
			origGridMFDs.add(future.join());
		baSol.setGridSourceProvider(avgProv);
		
		System.out.println("Writing solution");
		baSol.write(new File(outputDir, "solution_orig_rebuild.zip"));
		
		System.out.println("Now rebuilding with overall rate-balancing");
		
		treeFutures = new ArrayList<>();
		for (int b=0; b<gridTree.size(); b++) {
			LogicTreeBranch<?> branch = gridTree.getBranch(b);
//			System.out.println("Processing branch "+b+"/"+gridTree.size()+": "+branch);

			NSHM23_DeclusteringAlgorithms declusteringAlg = branch.requireValue(NSHM23_DeclusteringAlgorithms.class);
			NSHM23_SeisSmoothingAlgorithms seisSmooth = branch.requireValue(NSHM23_SeisSmoothingAlgorithms.class);
			
			// spatial seismicity PDF
			double[] pdf = seisSmooth.load(region, declusteringAlg);
			
			GridSourceList origGridProv = origGridMFDs.get(b);
			Preconditions.checkState(gridReg.getNodeCount() == origGridProv.getNumLocations());
			
			treeFutures.add(CompletableFuture.supplyAsync(() -> {
				
				EvenlyDiscretizedFunc origRef = origGridProv.getRefMFD();
				SummedMagFreqDist totalGR = new SummedMagFreqDist(origRef.getMinX(), origRef.getMaxX(), origRef.size());
				for (int i=0; i<gridReg.getNodeCount(); i++) {
					IncrementalMagFreqDist gridMFD = origGridProv.getMFD(i);
					Preconditions.checkState(gridMFD.size() <= origRef.size());
					Preconditions.checkState((float)gridMFD.getMinX() == (float)origRef.getMinX());
					for (int j=0; j<gridMFD.size(); j++)
						totalGR.add(j, gridMFD.getY(j));
				}
				
				NoFaultGridProv gridProv = new NoFaultGridProv(gridReg, trts, totalGR, pdf, fractStrikeSlip, fractReverse, fractNormal);
				return GridSourceList.convert(gridProv, assoc, new NSHM23_WUS_FiniteRuptureConverter());
			}));
		}
		
		avgProv = buildAvgGridProv(gridTree, treeFutures);
		baSol.setGridSourceProvider(avgProv);
		
		System.out.println("Writing solution");
		baSol.write(new File(outputDir, "solution_global_rate_balance.zip"));
		
		System.out.println("Now rebuilding using the regular NSHM23 methodology but the BA solution");
		
		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
		LogicTreeBranch<LogicTreeNode> avgFullBranch = NSHM23_LogicTreeBranch.AVERAGE_COMBINED;
		treeFutures = new ArrayList<>();
		for (int b=0; b<gridTree.size(); b++) {
			LogicTreeBranch<?> gridBranch = gridTree.getBranch(b);
			
			LogicTreeBranch<LogicTreeNode> fullBranch = avgFullBranch.copy();
			for (LogicTreeNode node : gridBranch)
				fullBranch.setValue(node);
			
			if (b == 0)
				factory.preGridBuildHook(origSol, fullBranch);
			
			treeFutures.add(CompletableFuture.supplyAsync(() -> {
				try {
					return factory.buildGridSourceProvider(origSol, fullBranch);
				} catch (IOException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
			}));
		}

		avgProv = buildAvgGridProv(gridTree, treeFutures);
		baSol.setGridSourceProvider(avgProv);
		
		System.out.println("Writing solution");
		baSol.write(new File(outputDir, "solution_normal_method.zip"));
	}
	
	private static GridSourceList buildAvgGridProv(LogicTree<?> gridTree, List<CompletableFuture<GridSourceList>> treeFutures) {
		Preconditions.checkState(treeFutures.size() == gridTree.size());
		AveragingAccumulator<GridSourceProvider> averager = null;
		for (int b=0; b<gridTree.size(); b++) {
			LogicTreeBranch<?> branch = gridTree.getBranch(b);
			GridSourceList gridList = treeFutures.get(b).join();
			if (averager == null)
				averager = gridList.averagingAccumulator();
			averager.process(gridList, gridTree.getBranchWeight(branch));
			System.out.println("Finished branch "+b+"/"+gridTree.size()+": "+branch);
		}
		System.out.println("Building average");
		return (GridSourceList) averager.getAverage();
	}
	
	private static class NoFaultGridProv extends MFDGridSourceProvider.Abstract {
		
		private GriddedRegion gridReg;
		private TectonicRegionType[] trts;
		private IncrementalMagFreqDist totalGR;
		private double[] pdf;
		private double[] fractStrikeSlip;
		private double[] fractReverse;
		private double[] fractNormal;

		public NoFaultGridProv(GriddedRegion gridReg, TectonicRegionType[] trts, IncrementalMagFreqDist totalGR,
				double[] pdf, double[] fractStrikeSlip, double[] fractReverse, double[] fractNormal) {
			this.gridReg = gridReg;
			this.trts = trts;
			this.totalGR = totalGR;
			this.pdf = pdf;
			this.fractStrikeSlip = fractStrikeSlip;
			this.fractReverse = fractReverse;
			this.fractNormal = fractNormal;
		}

		@Override
		public TectonicRegionType getTectonicRegionType(int gridIndex) {
			return trts[gridIndex];
		}

		@Override
		public MFDGridSourceProvider newInstance(Map<Integer, IncrementalMagFreqDist> nodeSubSeisMFDs,
				Map<Integer, IncrementalMagFreqDist> nodeUnassociatedMFDs, double[] fracStrikeSlip, double[] fracNormal,
				double[] fracReverse, TectonicRegionType[] trts) {
			throw new IllegalStateException();
		}

		@Override
		public IncrementalMagFreqDist getMFD_Unassociated(int gridIndex) {
			IncrementalMagFreqDist mfd = this.totalGR.deepClone();
			mfd.scale(pdf[gridIndex]);
			return mfd;
		}

		@Override
		public IncrementalMagFreqDist getMFD_SubSeisOnFault(int gridIndex) {
			return null;
		}

		@Override
		public GriddedRegion getGriddedRegion() {
			return gridReg;
		}

		@Override
		public double getFracStrikeSlip(int gridIndex) {
			return fractStrikeSlip[gridIndex];
		}

		@Override
		public double getFracReverse(int gridIndex) {
			return fractReverse[gridIndex];
		}

		@Override
		public double getFracNormal(int gridIndex) {
			return fractNormal[gridIndex];
		}

		@Override
		public void scaleAll(double[] valuesArray) {
			throw new IllegalStateException();
		}

		@Override
		public String getName() {
			return "No-Faults MFD";
		}

		@Override
		protected PointSource buildSource(int gridIndex, IncrementalMagFreqDist mfd, double duration,
				GriddedSeismicitySettings gridSourceSettings) {
			throw new IllegalStateException();
		}
		
	}

}
