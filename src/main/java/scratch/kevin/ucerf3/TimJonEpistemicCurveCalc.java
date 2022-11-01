package scratch.kevin.ucerf3;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.param.impl.DoubleParameter;
import org.opensha.commons.param.impl.StringParameter;
import org.opensha.commons.param.impl.WarningDoubleParameter;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.param.BPTAveragingTypeOptions;
import org.opensha.sha.earthquake.param.BPTAveragingTypeParam;
import org.opensha.sha.earthquake.param.HistoricOpenIntervalParam;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.imr.param.SiteParams.Vs30_TypeParam;

import com.google.common.base.Preconditions;

import scratch.UCERF3.erf.FaultSystemSolutionERF;

public class TimJonEpistemicCurveCalc {

	public static void main(String[] args) throws IOException {
		/*
		 * Curve calculation requested by Tim O'Donnell and Jon Stewart, via e-mail, 8/19/2022
		 * Subject: Continuing Epistemic Uncertainty in PSHA Research
		 */
		
		AttenRelRef[] gmpeRefs = {
				AttenRelRef.ASK_2014,
				AttenRelRef.BSSA_2014,
				AttenRelRef.CB_2014,
				AttenRelRef.CY_2014
		};
		double[] periods = { 0.2, 1d };
		boolean[] vsInferreds = {false, true};
		
//		String sitePrefix = "Davis";
//		Location siteLoc = new Location(38.3210, -121.4530);
		
		String sitePrefix = "Berkeley";
		Location siteLoc = new Location(37.5216, -122.1527);
		
		double[] vs30Vals = {
				// current set
				179,
				212,
				238,
				266,
				296,
				333,
				394
				// previous set
//				139.3,
//				171.8,
//				199.4,
//				228.2,
//				261.2,
//				303.2,
//				373.9
		};
		
		double[] percentiles = {
				0.1,
				1.3,
				2.5,
				5,
				9.25,
				16,
				33,
				50,
				67,
				84,
				90.75,
				95,
				97.5,
				98.7,
				99.9
		};
		
		// "USGS-2002 SA 0.3, 0.4, 0.5 and 1.0 Values" from the control panel
		DiscretizedFunc xVals = new ArbitrarilyDiscretizedFunc();
		xVals.set(.0025,1);
		xVals.set(.00375,1);
		xVals.set(.00563 ,1);
		xVals.set(.00844,1);
		xVals.set(.0127,1);
		xVals.set(.0190,1);
		xVals.set(.0285,1);
		xVals.set(.0427,1);
		xVals.set(.0641,1);
		xVals.set(.0961,1);
		xVals.set(.144,1);
		xVals.set(.216,1);
		xVals.set(.324,1);
		xVals.set(.487,1);
		xVals.set(.730,1);
		xVals.set(1.09,1);
		xVals.set(1.64,1);
		xVals.set(2.46,1);
		xVals.set(3.69,1);
		xVals.set(5.54,1);
		
		SolutionLogicTree slt = SolutionLogicTree.load(new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2021_11_30-u3_branches-orig_calcs-5h/results.zip"));
		LogicTree<?> tree = slt.getLogicTree();
//		tree = tree.sample(16, false);
		
		File outputDir = new File("/home/kevin/OpenSHA/UCERF3/2022_10-tim-jon-calcs");
		
		int threads = 28;
		ExecutorService exec = Executors.newFixedThreadPool(threads);
		
		DecimalFormat oDF = new DecimalFormat("0.##");
		
		for (AttenRelRef gmpeRef : gmpeRefs) {
			List<Future<Map<CalcKey, DiscretizedFunc>>> futures = new ArrayList<>();
			List<Double> weights = new ArrayList<>();
			for (LogicTreeBranch<?> branch : tree) {
				futures.add(exec.submit(new CurveCalcCallable(slt, branch, siteLoc, gmpeRef, periods, vs30Vals, vsInferreds, xVals)));
				weights.add(tree.getBranchWeight(branch));
			}
			
			Map<CalcKey, List<DiscretizedFunc>> curvesMap = new HashMap<>();
			
			for (int i=0; i<futures.size(); i++) {
				Future<Map<CalcKey, DiscretizedFunc>> future = futures.get(i);
				try {
					Map<CalcKey, DiscretizedFunc> curves = future.get();
					for (CalcKey key : curves.keySet()) {
						if (!curvesMap.containsKey(key))
							curvesMap.put(key, new ArrayList<>());
						curvesMap.get(key).add(curves.get(key));
					}
				} catch (InterruptedException | ExecutionException e) {
					throw ExceptionUtils.asRuntimeException(e);
				}
				System.out.println("DONE branch "+(i+1)+"/"+futures.size()+"\t"+gmpeRef.name());
			}
			
			for (double vs30 : vs30Vals) {
				for (boolean vsInf : vsInferreds) {
					for (double period : periods) {
						CSVFile<String> csv = new CSVFile<>(true);
						List<String> header = new ArrayList<>();
						header.add("Branch Index");
						header.add("Branch Weight");
						for (Point2D pt : xVals)
							header.add((float)pt.getX()+"");
						csv.addLine(header);
						
						List<DiscretizedFunc> curves = curvesMap.get(new CalcKey(period, vs30, vsInf));
						
						Preconditions.checkState(curves.size() == slt.getLogicTree().size());
						
						ArbDiscrEmpiricalDistFunc[] dists = new ArbDiscrEmpiricalDistFunc[xVals.size()];
						for (int i=0; i<dists.length; i++)
							dists[i] = new ArbDiscrEmpiricalDistFunc();
						
						for (int b=0; b<curves.size(); b++) {
							double weight = weights.get(b);
							DiscretizedFunc curve = curves.get(b);
							
							List<String> line = new ArrayList<>();
							line.add(b+"");
							line.add(weight+"");
							for (int i=0; i<curve.size(); i++) {
								line.add(curve.getY(i)+"");
								dists[i].set(curve.getY(i), weight);
							}
							csv.addLine(line);
						}
						
						List<String> meanLine = new ArrayList<>();
						meanLine.add("MEAN"); meanLine.add("");
						List<List<String>> percentileLines = new ArrayList<>();
						for (double p : percentiles) {
							List<String> percentileLine = new ArrayList<>();
							percentileLine.add("p"+oDF.format(p));
							percentileLine.add("");
							percentileLines.add(percentileLine);
						}
						
						for (int i=0; i<dists.length; i++) {
							meanLine.add(dists[i].getMean()+"");
							for (int p=0; p<percentiles.length; p++)
								percentileLines.get(p).add(dists[i].getInterpolatedFractile(percentiles[p]/100d)+"");
						}
						csv.addLine(meanLine);
						for (List<String> percentileLine : percentileLines)
							csv.addLine(percentileLine);
						
						String prefix = "curves";
						
						prefix += "_"+sitePrefix;
						
						prefix += "_sa_"+oDF.format(period)+"s";
						
						prefix += "_"+gmpeRef.name();
						
						prefix += "_vs30_"+oDF.format(vs30);
						if (vsInf)
							prefix += "_inferred";
						else
							prefix += "_measured";
						
						csv.writeToFile(new File(outputDir, prefix+".csv"));
					}
				}
			}
		}
		
		exec.shutdown();
	}
	
	private static class CalcKey {
		public final double period;
		public final double vs30;
		public final boolean vsInf;
		
		public CalcKey(double period, double vs30, boolean vsInf) {
			super();
			this.period = period;
			this.vs30 = vs30;
			this.vsInf = vsInf;
		}

		@Override
		public int hashCode() {
			return Objects.hash(period, vs30, vsInf);
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			CalcKey other = (CalcKey) obj;
			return Double.doubleToLongBits(period) == Double.doubleToLongBits(other.period)
					&& Double.doubleToLongBits(vs30) == Double.doubleToLongBits(other.vs30) && vsInf == other.vsInf;
		}
	}
	
	private static class CurveCalcCallable implements Callable<Map<CalcKey, DiscretizedFunc>> {
		
		private SolutionLogicTree slt;
		private LogicTreeBranch<?> branch;
		private Location siteLoc;
		private AttenRelRef gmpeRef;
		private double[] periods;
		private double[] vs30s;
		private boolean[] vsInfs;
		private DiscretizedFunc xVals;

		public CurveCalcCallable(SolutionLogicTree slt, LogicTreeBranch<?> branch, Location siteLoc,
				AttenRelRef gmpeRef, double[] periods, double[] vs30s, boolean[] vsInfs, DiscretizedFunc xVals) {
			this.slt = slt;
			this.branch = branch;
			this.siteLoc = siteLoc;
			this.gmpeRef = gmpeRef;
			this.periods = periods;
			this.vs30s = vs30s;
			this.vsInfs = vsInfs;
			this.xVals = xVals;
		}

		@Override
		public Map<CalcKey, DiscretizedFunc> call() throws Exception {
			Map<CalcKey, DiscretizedFunc> ret = new HashMap<>();
			
			ScalarIMR gmpe = gmpeRef.instance(null);
			
			gmpe.setParamDefaults();
			gmpe.setIntensityMeasure(SA_Param.NAME);
			
			FaultSystemSolution sol = slt.forBranch(branch);
			FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
			
			erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.INCLUDE);
			erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_PREF_BLEND);
			erf.setParameter(HistoricOpenIntervalParam.NAME, 146d);
			erf.setParameter(BPTAveragingTypeParam.NAME, BPTAveragingTypeOptions.AVE_RI_AVE_NORM_TIME_SINCE);
			erf.getTimeSpan().setDuration(1d);
			erf.getTimeSpan().setStartTime(2021);
			erf.updateForecast();
			
			HazardCurveCalculator calc = new HazardCurveCalculator();
			
			for (double vs30 : vs30s) {
				for (boolean vsInf : vsInfs) {
					Site site = new Site(siteLoc);
					boolean vsSet = false;
					boolean vsInfSet = false;
					for (Parameter<?> param : gmpe.getSiteParams()) {
						param = (Parameter<?>) param.clone();
						if (param.getName().equals(Vs30_Param.NAME)) {
							((WarningDoubleParameter)param).setValueIgnoreWarning(vs30);
							vsSet = true;
						}
						if (param.getName().equals(Vs30_TypeParam.NAME)) {
							((StringParameter)param).setValue(
									vsInf ? Vs30_TypeParam.VS30_TYPE_INFERRED : Vs30_TypeParam.VS30_TYPE_MEASURED);
							vsInfSet = true;
						}
						site.addParameter(param);
					}
					Preconditions.checkState(vsSet);
					Preconditions.checkState(vsInfSet);
					for (double period : periods) {
						SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), period);
						
						DiscretizedFunc logCurve = new ArbitrarilyDiscretizedFunc();
						for (Point2D pt : xVals)
							logCurve.set(Math.log(pt.getX()), 0d);
						
						calc.getHazardCurve(logCurve, site, gmpe, erf);
						
						DiscretizedFunc curve = new ArbitrarilyDiscretizedFunc();
						for (int i=0; i<xVals.size(); i++)
							curve.set(xVals.getX(i), logCurve.getY(i));
						ret.put(new CalcKey(period, vs30, vsInf), curve);
					}
				}
			}
			
			return ret;
		}
		
	}

}
