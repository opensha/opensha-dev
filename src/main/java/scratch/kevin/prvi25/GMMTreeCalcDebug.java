package scratch.kevin.prvi25;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CompletableFuture;
import java.util.function.Supplier;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.params.filters.SourceFilter;
import org.opensha.sha.calc.params.filters.SourceFilterManager;
import org.opensha.sha.calc.params.filters.SourceFilters;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.PointSource.PoissonPointSource;
import org.opensha.sha.earthquake.PointSource.PoissonPointSourceData;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.erf.BaseFaultSystemSolutionERF;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_LogicTreeHazardCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysHazardCalcSettings;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_GMM_GenericEpistemicModel;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_GMM_SlabEpistemicModel;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_GMM_SlabSigmaModel;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_LogicTree;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionSlabGMMs;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.nshmp.NSHMP_GMM_Wrapper;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import gov.usgs.earthquake.nshmp.gmm.Gmm;
import gov.usgs.earthquake.nshmp.gmm.GmmInput;
import gov.usgs.earthquake.nshmp.gmm.GroundMotion;
import gov.usgs.earthquake.nshmp.gmm.GroundMotionModel;
import gov.usgs.earthquake.nshmp.gmm.Imt;

public class GMMTreeCalcDebug {
	
	public static void main(String[] args) throws IOException {
		simpleGMMTest(Gmm.PRVI_2025_ACTIVE_CRUST, Imt.PGA);
//		System.exit(0);
//		Site site = new Site(new Location(18.4653, -66.1167));
		Site site = new Site(new Location(18, -68));
		
//		ScalarIMR avgGMM = AttenRelRef.USGS_PRVI_ACTIVE.get();
//		List<LogicTreeLevel<? extends LogicTreeNode>> gmmLevels = PRVI25_LogicTree.levelsCrustalGMM;
//		TectonicRegionType trt = TectonicRegionType.ACTIVE_SHALLOW;
		
//		ScalarIMR avgGMM = AttenRelRef.USGS_PRVI_INTERFACE.get();
//		List<LogicTreeLevel<? extends LogicTreeNode>> gmmLevels = PRVI25_LogicTree.levelsInterfaceGMM;
//		TectonicRegionType trt = TectonicRegionType.SUBDUCTION_INTERFACE;

		NSHMP_GMM_Wrapper avgGMM = (NSHMP_GMM_Wrapper)AttenRelRef.USGS_PRVI_SLAB.get();
		LogicTreeNode[] required = null;
		List<LogicTreeLevel<? extends LogicTreeNode>> gmmLevels = PRVI25_LogicTree.levelsSlabGMM;
//		NSHMP_GMM_Wrapper avgGMM = (NSHMP_GMM_Wrapper)AttenRelRef.USGS_PRVI_SLAB.get();
//		PRVI25_GMM_SlabEpistemicModel.EPI_OFF.setParams(avgGMM);
//		LogicTreeNode[] required = { PRVI25_GMM_SlabEpistemicModel.EPI_OFF };
//		List<LogicTreeLevel<? extends LogicTreeNode>> gmmLevels = PRVI25_LogicTree.levelsSlabGMM;
//		NSHMP_GMM_Wrapper avgGMM = (NSHMP_GMM_Wrapper)PRVI25_SubductionSlabGMMs.AS_PROVIDED.getSupplier().get();
//		LogicTreeNode[] required = { PRVI25_SubductionSlabGMMs.AS_PROVIDED };
//		NSHMP_GMM_Wrapper avgGMM = (NSHMP_GMM_Wrapper)PRVI25_SubductionSlabGMMs.DATA_ADJUSTED.getSupplier().get();
//		LogicTreeNode[] required = { PRVI25_SubductionSlabGMMs.DATA_ADJUSTED };
//		List<LogicTreeLevel<? extends LogicTreeNode>> gmmLevels = PRVI25_LogicTree.levelsSlabGMM;
		TectonicRegionType trt = TectonicRegionType.SUBDUCTION_SLAB;
		
		FaultSystemSolution sol = FaultSystemSolution.load(
				new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
						+ "2025_08_01-prvi25_crustal_subduction_combined_branches/combined_branch_averaged_solution.zip"));
		
		BaseFaultSystemSolutionERF erf = new BaseFaultSystemSolutionERF();
		erf.setSolution(sol);
		erf.getTimeSpan().setDuration(1d);
		erf.updateForecast();
		
		TRTWrappedERF trtERF = new TRTWrappedERF(erf, trt, true);
		
		site.addParameterList(avgGMM.getSiteParams());
		
		avgGMM.setIntensityMeasure(PGA_Param.NAME);
		
		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(PGA_Param.NAME);
		ArbitrarilyDiscretizedFunc newXVals = new ArbitrarilyDiscretizedFunc();
		newXVals.set(1e-8, 1d);
		newXVals.set(1e-7, 1d);
		newXVals.set(1e-6, 1d);
		newXVals.set(1e-5, 1d);
		for (Point2D pt : xVals)
			newXVals.set(pt);
		xVals = newXVals;
		
		CompletableFuture<DiscretizedFunc> avgCurveFuture = CompletableFuture.supplyAsync(new CurveCalc(site, trtERF, avgGMM, xVals));
		
		LogicTree<LogicTreeNode> gmmTree = LogicTree.buildExhaustive(gmmLevels, true, required);
		avgGMM.setCurrentGmmInput(defaultInput);
		System.out.println("Test rupture input: "+defaultInput+"; IMT="+avgGMM.getCurrentIMT());
		System.out.println("Average tree:\n"+avgGMM.getGroundMotionTree());
		
		double sumWeight = 0d;
		for (LogicTreeBranch<?> branch : gmmTree) {
			double weight = gmmTree.getBranchWeight(branch);
			sumWeight += weight;
		}
		
		List<Double> treeWeights = new ArrayList<>();
		List<CompletableFuture<DiscretizedFunc>> treeCurveFutures = new ArrayList<>();
		
		for (LogicTreeBranch<?> branch : gmmTree) {
			double weight = gmmTree.getBranchWeight(branch);
			if ((float)sumWeight != 1f)
				weight /= sumWeight;
			Map<TectonicRegionType, ? extends Supplier<ScalarIMR>> gmmSupplier = FaultSysHazardCalcSettings.getGMM_Suppliers(branch, null, false);
			NSHMP_GMM_Wrapper gmm = (NSHMP_GMM_Wrapper)gmmSupplier.get(trt).get();
//			NSHMP_GMM_Wrapper gmm = new NSHMP_GMM_Wrapper.Single(branch.requireValue(PRVI25_SubductionSlabGMMs.class).getGMM(), false);
//			branch.requireValue(PRVI25_GMM_SlabEpistemicModel.class).getModel().setParams(gmm);
//			branch.requireValue(PRVI25_GMM_SlabSigmaModel.class).getModel().setParams(gmm);
			gmm.setIntensityMeasure(avgGMM.getIntensityMeasure().getName());
			System.out.println("Branch "+branch+" with weight="+(float)weight+" has GMM: "+gmm.getName());
//			gmm.setSite(site);
//			gmm.setEqkRupture(testRup);
			gmm.setCurrentGmmInput(defaultInput);
//			System.out.println("Test rupture input: "+gmm.getCurrentGmmInput()+"; IMT="+gmm.getCurrentIMT());
			System.out.println( gmm.getGroundMotionTree());
			
			CompletableFuture<DiscretizedFunc> curveFuture = CompletableFuture.supplyAsync(new CurveCalc(site, trtERF, gmm, xVals));
			treeWeights.add(weight);
			treeCurveFutures.add(curveFuture);
		}
		
		DiscretizedFunc calcAvgCurve = xVals.deepClone();
		calcAvgCurve.scale(0d);
		
		for (int i=0; i<treeCurveFutures.size(); i++) {
			DiscretizedFunc func = treeCurveFutures.get(i).join();
//			System.out.println("Curve:\n"+func);
			double weight = treeWeights.get(i);
			for (int j = 0; j < func.size(); j++) {
				double y = func.getY(j);
				Preconditions.checkState(y >= 0 && (RATE_CURVES || y <= 1), "Bad y=%s", y);
				calcAvgCurve.set(j, calcAvgCurve.getY(j) + y * weight);
			}
		}
		
		DiscretizedFunc avgCurve = avgCurveFuture.join();
		
//		System.out.println("Average curve:\n"+avgCurve);
//		System.out.println("Calculated average curve:\n"+calcAvgCurve);
		
		for (int i=0; i<avgCurve.size(); i++) {
			double x = avgCurve.getX(i);
			double y1 = avgCurve.getY(i);
			double y2 = calcAvgCurve.getY(i);
			double diff = y2 - y1;
			double pDiff = 100d * diff / y1;
			System.out.println("X: "+(float)x+"\tAvg: "+(float)y1+"\tCalc Avg: "+(float)y2+"\tDiff: "+(float)diff+" ("+(float)pDiff+" %)");
		}
		
		ReturnPeriods rp = ReturnPeriods.TWO_IN_50;
		double rpProb = RATE_CURVES ? 1d/rp.returnPeriod : rp.oneYearProb;
		System.out.println("Return period: "+rp+" at y="+rpProb);
		double rp1 = avgCurve.getFirstInterpolatedX_inLogXLogYDomain(rpProb);
		double rp2 = calcAvgCurve.getFirstInterpolatedX_inLogXLogYDomain(rpProb);
		double diff = rp2 - rp1;
		double pDiff = 100d * diff / rp1;
		System.out.println("Avg 2in50: "+(float)rp1+"\tCalc Avg 2in50: "+(float)rp2+"\tDiff: "+(float)diff+" ("+(float)pDiff+" %)");
	}
	
	static class TRTWrappedERF extends AbstractERF {
		
		private List<ProbEqkSource> sources;
		private TectonicRegionType trt;

		public TRTWrappedERF(ERF erf, TectonicRegionType trt, Boolean poisson) {
			this.trt = trt;
			sources = new ArrayList<>();
			for (ProbEqkSource source : erf) {
				if (source.getTectonicRegionType() == trt) {
					if (poisson != null && poisson != source.isSourcePoissonian())
						sources.add(new PoissonOverrideSource(source, poisson));
					else
						sources.add(source);
				}
			}
		}

		@Override
		public int getNumSources() {
			return sources.size();
		}

		@Override
		public ProbEqkSource getSource(int idx) {
			return sources.get(idx);
		}

		@Override
		public void updateForecast() {
			
		}

		@Override
		public String getName() {
			return trt.toString();
		}
		
	}
	
	private static class PoissonOverrideSource extends ProbEqkSource {
		
		private ProbEqkSource orig;

		public PoissonOverrideSource(ProbEqkSource orig, boolean poisson) {
			this.orig = orig;
			this.isPoissonian = poisson;
		}

		@Override
		public RuptureSurface getSourceSurface() {
			return orig.getSourceSurface();
		}
		
		@Override
		public LocationList getAllSourceLocs() {
			return orig.getAllSourceLocs();
		}
		
		@Override
		public ProbEqkRupture getRupture(int nRupture) {
			return orig.getRupture(nRupture);
		}
		
		@Override
		public int getNumRuptures() {
			return orig.getNumRuptures();
		}
		
		@Override
		public double getMinDistance(Site site) {
			return orig.getNumRuptures();
		}
		
	}
	
	private static class CurveCalc implements Supplier<DiscretizedFunc> {
		
		private Site site;
		private AbstractERF erf;
		private ScalarIMR gmpe;
		private DiscretizedFunc xVals;

		public CurveCalc(Site site, AbstractERF erf, ScalarIMR gmpe, DiscretizedFunc xVals) {
			this.site = site;
			this.erf = erf;
			this.gmpe = gmpe;
			this.xVals = xVals;
		}

		@Override
		public DiscretizedFunc get() {
			double[] logXValsArray = new double[xVals.size()];
			for (int i=0; i<xVals.size(); i++)
				logXValsArray[i] = Math.log(xVals.getX(i));
			LightFixedXFunc logXVals = new LightFixedXFunc(logXValsArray, new double[xVals.size()]);
			HazardCurveCalculator calc = new HazardCurveCalculator(new SourceFilterManager(SourceFilters.TRT_DIST_CUTOFFS));
			
			gmpe.setIntensityMeasure(PGA_Param.NAME);
			
			if (RATE_CURVES)
				rateBasedHazardCalc(site, erf, gmpe, logXVals, calc.getSourceFilters());
			else
				calc.getHazardCurve(logXVals, site, gmpe, erf);
			
			DiscretizedFunc ret = xVals.deepClone();
			for (int i=0; i<xVals.size(); i++)
				ret.set(i, logXVals.getY(i));
			
			return ret;
		}
		
	}
	
	private static final boolean RATE_CURVES = false;
	
	private static final GmmInput defaultInput = GmmInput.builder().withDefaults().mag(7d).rJB(10d).rRup(10d).rX(1d).dip(90d).width(15d)
			.zTor(0d).zHyp(5d).rake(0d).vs30(760).z1p0(0.2).z2p5(1d).build();
	
	private static void simpleGMMTest(Gmm avgGMM, Imt imt) {
		System.out.println("GMM input: "+defaultInput);
		GroundMotionModel avgModel = avgGMM.instance(imt);
		System.out.println("Average GMM: "+avgGMM);
		gov.usgs.earthquake.nshmp.tree.LogicTree<GroundMotion> avgTree = avgModel.calc(defaultInput);
		System.out.println("Average tree:\n"+avgTree);
	}
	
	private static void rateBasedHazardCalc(Site site, AbstractERF erf, ScalarIMR gmpe, LightFixedXFunc curve, Collection<SourceFilter> filters) {
		curve.scale(0d);
		gmpe.setSite(site);
		
		LightFixedXFunc exceedProbs = new LightFixedXFunc(curve.getXVals(), new double[curve.size()]);
		
		for (ProbEqkSource source : erf) {
			if (HazardCurveCalculator.canSkipSource(filters, source, site))
				continue;
			for (int rupID=0; rupID<source.getNumRuptures(); rupID++) {
				ProbEqkRupture rup = source.getRupture(rupID);
				if (HazardCurveCalculator.canSkipRupture(filters, rup, site))
					continue;
				
				double rate;
				if (source instanceof PoissonPointSource) {
					PoissonPointSourceData data = ((PoissonPointSource)source).getData();
					if (data.getNumRuptures() == source.getNumRuptures())
						rate = data.getRate(rupID);
					else
						rate = rup.getMeanAnnualRate(1d);
				} else {
					rate = rup.getMeanAnnualRate(1d);
				}
				
				gmpe.setEqkRupture(rup);
				
				gmpe.getExceedProbabilities(exceedProbs);
				
				for (int i=0; i<curve.size(); i++)
					curve.set(i, curve.getY(i) + rate*exceedProbs.getY(i));
			}
		}
	}

}
