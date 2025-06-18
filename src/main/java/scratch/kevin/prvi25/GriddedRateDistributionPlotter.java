package scratch.kevin.prvi25;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.commons.lang3.exception.ExceptionUtils;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertaintyBoundType;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.util.io.archive.ArchiveInput;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;

import net.mahdilamb.colormap.Colors;

public class GriddedRateDistributionPlotter {

	public static void main(String[] args) throws IOException {
		File invsDir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/");
		
		File randDir = new File(invsDir, "2025_05_13-prvi_gridded_variability-full_random");
		File threeBranchDir = new File(invsDir, "2025_05_13-prvi_gridded_variability-three_branch");
		
		SolutionLogicTree randSLT = SolutionLogicTree.load(
				new ArchiveInput.PreloadingZipFileInput(new File(randDir, "results.zip")));
		SolutionLogicTree threeBranchSLT = SolutionLogicTree.load(
				new ArchiveInput.PreloadingZipFileInput(new File(threeBranchDir, "results.zip")));
		
//		randSLT = new SolutionLogicTree.SubsetSolutionLogicTree(randSLT, randSLT.getLogicTree().sample(500, true));
//		threeBranchSLT = new SolutionLogicTree.SubsetSolutionLogicTree(threeBranchSLT, threeBranchSLT.getLogicTree().sample(500, true));
		
		File outputDir = new File(randDir, "mfds");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(5.01, 8.55);

		System.out.println("Loading MFDs for "+randSLT.getLogicTree().size()+" fully random branches");
		Table<LogicTreeBranch<?>, TectonicRegionType, IncrementalMagFreqDist> randMFDs = loadBranchMFDs(randSLT, refMFD);
		System.out.println("Loading MFDs for "+threeBranchSLT.getLogicTree().size()+" three-branch branches");
		Table<LogicTreeBranch<?>, TectonicRegionType, IncrementalMagFreqDist> threeBranchMFDs = loadBranchMFDs(threeBranchSLT, refMFD);
		
		List<TectonicRegionType> trts = new ArrayList<>();
		trts.add(null);
		trts.addAll(randMFDs.columnKeySet());
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		Range xRange = new Range(5d, 8.5d);
		Range yRange = new Range(1e-6d, 1e1d);
		
		int transAlpha = 60;
		Color randColor = Colors.tab_orange;
		Color threeBranchColor = Colors.tab_blue;
		Color randAlphaColor = new Color(randColor.getRed(), randColor.getGreen(), randColor.getBlue(), transAlpha);
		Color threeBranchAlphaColor = new Color(threeBranchColor.getRed(), threeBranchColor.getGreen(), threeBranchColor.getBlue(), transAlpha);
		Color indvColor = new Color(0, 0, 0, 20);
		PlotCurveCharacterstics indvChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, indvColor);
		
		for (TectonicRegionType trt : trts) {
			Map<LogicTreeBranch<?>, IncrementalMagFreqDist> trtRandMFDs = trt == null ? sumTRTs(refMFD, randMFDs) : randMFDs.column(trt);
			IncrementalMagFreqDist randMean = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
			UncertainBoundedIncrMagFreqDist[] boundsRand = calcBounds(trtRandMFDs, randMean, randSLT.getLogicTree());
			Map<LogicTreeBranch<?>, IncrementalMagFreqDist> trtThreeBranchMFDs = trt == null ? sumTRTs(refMFD, threeBranchMFDs) : threeBranchMFDs.column(trt);
			IncrementalMagFreqDist threeBranchMean = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
			UncertainBoundedIncrMagFreqDist[] boundsThreeBranch = calcBounds(trtThreeBranchMFDs, threeBranchMean, threeBranchSLT.getLogicTree());
			
			for (int i=0; i<3; i++) {
				List<DiscretizedFunc> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				if (i == 0 || i == 2) {
					randMean.setName("Full Random");
					funcs.add(randMean);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, randColor));
					
					boundsRand[0].setName("95% and 68% bounds");
					funcs.add(boundsRand[0]);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, randAlphaColor));
					
					boundsRand[1].setName(null);
					funcs.add(boundsRand[1]);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, randAlphaColor));
				}
				
				if (i == 1 || i == 2) {
					threeBranchMean.setName("Three Branch");
					funcs.add(threeBranchMean);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, threeBranchColor));
					
					boundsThreeBranch[0].setName("95% and 68% bounds");
					funcs.add(boundsThreeBranch[0]);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, threeBranchAlphaColor));
					
					boundsThreeBranch[1].setName(null);
					funcs.add(boundsThreeBranch[1]);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, threeBranchAlphaColor));
				}
				
				if (i == 2) {
					// add means again on top
					IncrementalMagFreqDist meanCopy = new IncrementalMagFreqDist(randMean);
					meanCopy.setName(null);
					funcs.add(meanCopy);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, randColor));

					meanCopy = new IncrementalMagFreqDist(threeBranchMean);
					meanCopy.setName(null);
					funcs.add(meanCopy);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, threeBranchColor));
				}
				
				String title = trt == null ? "All Regimes" : trt.toString();
				String prefix;
				if (i == 0)
					prefix = "full_rand";
				else if (i == 1)
					prefix = "three_branch";
				else
					prefix = "combined";
				if (trt != null)
					prefix += "_"+trt.name();
				PlotSpec plot = new PlotSpec(funcs, chars, title, "Magnitude", "Incremental Rate (1/yr)");
				plot.setLegendInset(RectangleAnchor.TOP_RIGHT);
				
				gp.drawGraphPanel(plot, false, true, xRange, yRange);
				
				PlotUtils.writePlots(outputDir, prefix, gp, 800, 800, true, true, false);
				
				if (i < 2) {
					// idividual
					funcs.clear();
					chars.clear();
					
					Collection<IncrementalMagFreqDist> mfds = i == 0 ? trtRandMFDs.values() : trtThreeBranchMFDs.values();
					
					for (IncrementalMagFreqDist mfd : mfds) {
						if (funcs.isEmpty())
							mfd.setName("Individual branches");
						else
							mfd.setName(null);
						funcs.add(mfd);
						chars.add(indvChar);
					}
					
					if (i == 0) {
						funcs.add(randMean);
						chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, randColor));
					} else {
						funcs.add(threeBranchMean);
						chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, threeBranchColor));
					}
					
					plot = new PlotSpec(funcs, chars, title, "Magnitude", "Incremental Rate (1/yr)");
					plot.setLegendInset(RectangleAnchor.TOP_RIGHT);
					
					gp.drawGraphPanel(plot, false, true, xRange, yRange);
					
					PlotUtils.writePlots(outputDir, prefix+"_indv", gp, 800, 800, true, true, false);
				}
			}
		}
	}
	
	private static Table<LogicTreeBranch<?>, TectonicRegionType, IncrementalMagFreqDist> loadBranchMFDs(
			SolutionLogicTree slt, EvenlyDiscretizedFunc refMFD) throws IOException {
		Table<LogicTreeBranch<?>, TectonicRegionType, IncrementalMagFreqDist> ret = HashBasedTable.create();
		
		slt.setVerbose(false);
		
		double minMag = refMFD.getMinX()-0.51*refMFD.getDelta();
		
		LinkedList<CompletableFuture<Void>> futures = new LinkedList<>();
		int maxParallel = 5;
		final AtomicInteger count = new AtomicInteger(0);
		int size = slt.getLogicTree().size();
		for (LogicTreeBranch<?> branch : slt.getLogicTree()) {
			if (futures.size() == maxParallel)
				futures.removeFirst().join();
			futures.add(CompletableFuture.runAsync(new Runnable() {
				
				@Override
				public void run() {
					GridSourceProvider prov;
					try {
						prov = slt.loadGridProvForBranch(branch);
					} catch (IOException e) {
						throw ExceptionUtils.asRuntimeException(e);
					}
					Preconditions.checkState(prov instanceof GridSourceList);
					GridSourceList gridList = (GridSourceList)prov;
					Map<TectonicRegionType, IncrementalMagFreqDist> map = new HashMap<>();
					for (TectonicRegionType trt : gridList.getTectonicRegionTypes()) {
						IncrementalMagFreqDist mfd = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
						for (int l=0; l<gridList.getNumLocations(); l++)
							for (GriddedRupture rup : gridList.getRuptures(trt, l))
								if (rup.properties.magnitude > minMag)
									mfd.add(mfd.getClosestXIndex(rup.properties.magnitude), rup.rate);
						map.put(trt, mfd);
					}
					synchronized (ret) {
						for (TectonicRegionType trt : map.keySet())
							ret.put(branch, trt, map.get(trt));
						int myCount = count.incrementAndGet();
						System.out.print(".");
						if (myCount % 100 == 0)
							System.out.println(" "+count+"/"+size);
					}
				}
			}));
		}
		
		for (CompletableFuture<Void> future : futures)
			future.join();
		System.out.println();
		
		return ret;
	}
	
	public static Map<LogicTreeBranch<?>, IncrementalMagFreqDist> sumTRTs(EvenlyDiscretizedFunc refMFD,
			Table<LogicTreeBranch<?>, TectonicRegionType, IncrementalMagFreqDist> table) {
		Map<LogicTreeBranch<?>, IncrementalMagFreqDist> ret = new HashMap<>();
		for (LogicTreeBranch<?> branch : table.rowKeySet()) {
			SummedMagFreqDist sum = new SummedMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
			for (IncrementalMagFreqDist mfd : table.row(branch).values())
				sum.addIncrementalMagFreqDist(mfd);
			ret.put(branch, sum);
		}
		return ret;
	}
	
	public static UncertainBoundedIncrMagFreqDist[] calcBounds(Map<LogicTreeBranch<?>, IncrementalMagFreqDist> branchMFDs,
			IncrementalMagFreqDist mean, LogicTree<?> tree) {
		double[] weights = new double[tree.size()];
		double[] means = new double[mean.size()];
		LightFixedXFunc[] normCDFs = new LightFixedXFunc[mean.size()];
		int maxNonzeroIndex = 0;
		for (int i=0; i<mean.size(); i++) {
			double[] values = new double[weights.length];
			boolean anyNonzero = false;
			double weightedSum = 0d;
			double sumWeights = 0d;
			for (int b=0; b<weights.length; b++) {
				LogicTreeBranch<?> branch = tree.getBranch(b);
				if (i == 0)
					weights[b] = tree.getBranchWeight(b);
				values[b] = branchMFDs.get(branch).getY(i);
				anyNonzero |= values[b] > 0d;
				weightedSum += values[b]*weights[b];
				sumWeights += weights[b];
			}
			if (anyNonzero) {
				maxNonzeroIndex = i;
				normCDFs[i] = ArbDiscrEmpiricalDistFunc.calcQuickNormCDF(values, weights);
				means[i] = weightedSum/sumWeights;
			}
		}
		Preconditions.checkState(maxNonzeroIndex > 0, "No nonzero MFDs?");
		
		IncrementalMagFreqDist p2p5 = new IncrementalMagFreqDist(mean.getMinX(), maxNonzeroIndex+1, mean.getDelta());
		IncrementalMagFreqDist p16 = new IncrementalMagFreqDist(p2p5);
		IncrementalMagFreqDist p50 = new IncrementalMagFreqDist(p2p5);
		IncrementalMagFreqDist p84 = new IncrementalMagFreqDist(p2p5);
		IncrementalMagFreqDist p97p5 = new IncrementalMagFreqDist(p2p5);
		for (int i=0; i<p2p5.size(); i++) {
			Preconditions.checkNotNull(normCDFs[i], "No normCDF for %s. %s, but mean=%s", i, mean.getX(i), means[i]);
			mean.set(i, means[i]);
			p2p5.set(i, ArbDiscrEmpiricalDistFunc.calcFractileFromNormCDF(normCDFs[i], 0.025));
			p16.set(i, ArbDiscrEmpiricalDistFunc.calcFractileFromNormCDF(normCDFs[i], 0.16));
			p50.set(i, ArbDiscrEmpiricalDistFunc.calcFractileFromNormCDF(normCDFs[i], 0.5));
			p84.set(i, ArbDiscrEmpiricalDistFunc.calcFractileFromNormCDF(normCDFs[i], 0.84));
			p97p5.set(i, ArbDiscrEmpiricalDistFunc.calcFractileFromNormCDF(normCDFs[i], 0.975));
		}
		
		UncertainBoundedIncrMagFreqDist[] ret = {
				new UncertainBoundedIncrMagFreqDist(p50, p2p5, p97p5, UncertaintyBoundType.CONF_95),
				new UncertainBoundedIncrMagFreqDist(p50, p16, p84, UncertaintyBoundType.CONF_68)
		};
		return ret;
	}

}
