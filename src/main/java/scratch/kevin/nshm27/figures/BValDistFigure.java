package scratch.kevin.nshm27.figures;

import static scratch.kevin.nshm27.figures.NSHM27_PaperPaths.*;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;

import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.commons.statistics.distribution.ContinuousDistribution;
import org.apache.commons.statistics.distribution.UniformContinuousDistribution;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeFigureWriter;
import org.opensha.commons.logicTree.LogicTreeLevel.SamplingMethod;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.util.modules.ModuleContainer;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SectionSupraSeisBValues;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SectionSupraSeisBValues.DistributionSamplingLevel;
import org.opensha.sha.util.TectonicRegionType;

import gov.usgs.earthquake.nshmp.erf.nshm27.NSHM27_InvConfigFactory;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_InterfaceHingedBValue;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_InterfaceObsSeisDMAdjustment;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_LogicTree;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_SeisClassificationMethod;
import gov.usgs.earthquake.nshmp.erf.nshm27.logicTree.NSHM27_SeisRateModel;
import gov.usgs.earthquake.nshmp.erf.nshm27.util.NSHM27_RegionLoader.NSHM27_SeismicityRegions;
import gov.usgs.earthquake.nshmp.erf.seismicity.SeismicityRateFileLoader.PureGR;
import net.mahdilamb.colormap.Colors;

public class BValDistFigure {

	public static void main(String[] args) throws IOException {
		System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", FaultSysTools.defaultNumThreads()+"");
//		int numSamples = 100;
//		int numSamples = 1000;
//		int numSamples = 2000;
//		int numSamples = 5000;
		int numSamples = 10000;
//		int numSamples = 50000;
//		double binWidth = 0.025;
		ModuleContainer.VERBOSE_DEFAULT = false;
		double binWidth = 0.05;
		Range bRange = new Range(-0.6, 1.6);
		Range yRange = new Range(0, 1.4);
		
		Color hingedColor = Colors.tab_blue;
		hingedColor = new Color(hingedColor.getRed(), hingedColor.getGreen(), hingedColor.getBlue(), 200);
		Color hingedBelowOverlayColor = new Color(255, 255, 255, 160);
		Color extrapColor = Colors.tab_green;
		extrapColor = new Color(extrapColor.getRed(), extrapColor.getGreen(), extrapColor.getBlue(), 200);
		Color combColor = Colors.tab_lightred;
		
		double extrapolateWeight = NSHM27_InterfaceObsSeisDMAdjustment.EXTRAPOLATE.getNodeWeight();
		double otherWeight = 1d-extrapolateWeight;
		double hingeWeight = NSHM27_LogicTree.INTERFACE_B_HINGED_WEIGHT;
		double distWeight = 1d - hingeWeight;
		hingeWeight *= otherWeight;
		distWeight *= otherWeight;
		
		NSHM27_LogicTree.INTERFACE_B_HINGED_WEIGHT = 1d;
		
		EvenlyDiscretizedFunc bValDiscr = HistogramFunction.getEncompassingHistogram(-1, 2, binWidth);
		
		NSHM27_InvConfigFactory factory = new NSHM27_InvConfigFactory();
		
		for (NSHM27_SeismicityRegions seisReg : NSHM27_SeismicityRegions.values()) {
			LogicTree<LogicTreeNode> tree = NSHM27_LogicTree.buildLogicTree(seisReg,
					TectonicRegionType.SUBDUCTION_INTERFACE, numSamples, true,
					SamplingMethod.PAIRWISE_OPTIMIZED_LATIN_HYPERCUBE);
			tree.write(new File("/tmp/test_tree_"+seisReg.name()+".json"));
			ContinuousDistribution dist = NSHM27_LogicTree.getInterfaceBDist(seisReg);
			DistributionSamplingLevel distLevel = new DistributionSamplingLevel("Dist", "Dist", dist);
			String distLabel = LogicTreeFigureWriter.getDistString(distLevel);
			
			List<XY_DataSet> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			FaultSystemRupSet genericRS = factory.buildGenericRupSet(tree.getBranch(0), FaultSysTools.defaultNumThreads());
			
			List<LogicTreeBranch<LogicTreeNode>> hingedBranches = new ArrayList<>();
			List<LogicTreeBranch<LogicTreeNode>> extrapolatedBranches = new ArrayList<>();
			for (LogicTreeBranch<LogicTreeNode> branch : tree) {
				if (branch.hasValue(NSHM27_InterfaceObsSeisDMAdjustment.EXTRAPOLATE))
					extrapolatedBranches.add(branch);
				else
					hingedBranches.add(branch);
			}
			
			// hinge samples
			List<CompletableFuture<Double>> bFutures = new ArrayList<>(numSamples);
			for (LogicTreeBranch<LogicTreeNode> branch : hingedBranches) {
				bFutures.add(CompletableFuture.supplyAsync(()->{
					FaultSystemRupSet rs;
					try {
						rs = factory.updateRuptureSetForBranch(genericRS, branch);
					} catch (IOException e) {
						throw ExceptionUtils.asRuntimeException(e);
					}
					return NSHM27_InterfaceHingedBValue.calcRawInterfaceHingedBValue(rs, branch);
				}));
			}
			
			for (int i=0; i<bFutures.size(); i++) {
				LogicTreeBranch<LogicTreeNode> branch = hingedBranches.get(i);
				try {
					bFutures.get(i).get();
				} catch (Exception e) {
					System.err.println("Failed on "+branch);
					e.printStackTrace();
					System.err.flush();
					System.exit(1);
				}
			}
			
			EvenlyDiscretizedFunc hingedHist = new EvenlyDiscretizedFunc(
					bValDiscr.getMinX(), bValDiscr.getMaxX(), bValDiscr.size());

			double densityEach = hingeWeight/(bFutures.size()*hingedHist.getDelta());
			for (CompletableFuture<Double> future : bFutures) {
				double b = future.join();
				int index = hingedHist.getClosestXIndex(b);
				hingedHist.add(index, densityEach);
			}
			
			for (LogicTreeBranch<LogicTreeNode> branch : hingedBranches) {
				bFutures.add(CompletableFuture.supplyAsync(()->{
					FaultSystemRupSet rs;
					try {
						rs = factory.updateRuptureSetForBranch(genericRS, branch);
					} catch (IOException e) {
						throw ExceptionUtils.asRuntimeException(e);
					}
					return NSHM27_InterfaceHingedBValue.calcRawInterfaceHingedBValue(rs, branch);
				}));
			}
			
			EvenlyDiscretizedFunc extrapHist = new EvenlyDiscretizedFunc(
					bValDiscr.getMinX(), bValDiscr.getMaxX(), bValDiscr.size());

			densityEach = extrapolateWeight/(extrapolatedBranches.size()*hingedHist.getDelta());
			for (LogicTreeBranch<LogicTreeNode> branch : extrapolatedBranches) {
				double b = ((PureGR)branch.requireValue(NSHM27_SeisRateModel.class).getRateRecord(
						seisReg, branch.requireValue(NSHM27_SeisClassificationMethod.class), TectonicRegionType.SUBDUCTION_INTERFACE)).b;
				int index = extrapHist.getClosestXIndex(b);
				extrapHist.add(index, densityEach);
			}

			EvenlyDiscretizedFunc totalHist = hingedHist.deepClone();

			// The hinge calculation is truncated at b = 0.
			int zeroBin = totalHist.getClosestXIndex(0d);
			for (int i=0; i<zeroBin; i++) {
				totalHist.add(zeroBin, totalHist.getY(i));
				totalHist.set(i, 0d);
			}
			
			// add extrap
			for (int i=0; i<totalHist.size(); i++)
				totalHist.add(i, extrapHist.getY(i));

			for (int i=0; i<totalHist.size(); i++) {
				double x = totalHist.getX(i);
				double start = x - 0.5*totalHist.getDelta();
				double end = x + 0.5*totalHist.getDelta();

				double distProb = dist.cumulativeProbability(end)
						- dist.cumulativeProbability(start);

				totalHist.add(i, distWeight*distProb/totalHist.getDelta());
			}
			
			XY_DataSet distFunc;
			if (dist instanceof UniformContinuousDistribution) {
				distFunc = new DefaultXY_DataSet();
				distFunc.set(dist.getSupportLowerBound(), 0d);
				double density = dist.density(0.5*(dist.getSupportUpperBound() + dist.getSupportLowerBound()));
				distFunc.set(dist.getSupportLowerBound(), distWeight*density);
				distFunc.set(dist.getSupportUpperBound(), distWeight*density);
				distFunc.set(dist.getSupportUpperBound(), 0d);
			} else {
				distFunc = new EvenlyDiscretizedFunc(bRange.getLowerBound(), bRange.getUpperBound(), 1000);
				for (int i=0; i<distFunc.size(); i++)
					distFunc.set(i, distWeight*dist.density(distFunc.getX(i)));
			}
			
			distFunc.setName(distLabel);
			funcs.add(distFunc);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
			
			EvenlyDiscretizedFunc hingedBelow = hingedHist.deepClone();
			for (int i=hingedBelow.getClosestXIndex(0.000001); i<hingedBelow.size(); i++)
				hingedBelow.set(i, 0d);
			hingedBelow.setName(null);
			funcs.add(hingedBelow);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, hingedBelowOverlayColor));
			
			hingedHist.setName("Hinged");
			funcs.add(hingedHist);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, hingedColor));
			
			extrapHist.setName("Extrapolated");
			funcs.add(extrapHist);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, extrapColor));
			
			totalHist.setName("Combined");
			funcs.add(totalHist);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, combColor));
			
			PlotSpec plot = new PlotSpec(funcs, chars, seisReg.getTitleCaseAcronym(), "Interface b-value", "Density");
//			plot.setLegendInset(true);
			plot.setLegendInset(RectangleAnchor.TOP_LEFT);
//			plot.setLegendVisible(true);
			
			System.out.println("Hinged fract: "+hingedHist.calcSumOfY_Vals()/totalHist.calcSumOfY_Vals());
			System.out.println("Extrap fract: "+extrapHist.calcSumOfY_Vals()/totalHist.calcSumOfY_Vals());
			System.out.println("Dist fract: "+(1 - (hingedHist.calcSumOfY_Vals()+extrapHist.calcSumOfY_Vals())/totalHist.calcSumOfY_Vals()));
//			System.exit(0);
			
			HeadlessGraphPanel gp = PlotUtils.initPrintHeadless();
			
			gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
			
			gp.drawGraphPanel(plot, false, false, bRange, yRange);
			
			PlotUtils.writePrintPlots(FIGURES_DIR, "interface_b_"+seisReg.name(), gp,
					PlotUtils.DEFAULT_USABLE_PAGE_WIDTH*2d/3d, 3d, 150, true, true, false);
		}
	}

}
