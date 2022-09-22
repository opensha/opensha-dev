package scratch.kevin.ucerf3;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
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

import scratch.UCERF3.erf.FaultSystemSolutionERF;

public class TimothyJonEpistemicCurveCalc {

	public static void main(String[] args) throws IOException {
		/*
		 * Curve calculation requested by Tim O'Donnell and Jon Stewart, via e-mail, 8/19/2022
		 * Subject: Continuing Epistemic Uncertainty in PSHA Research
		 */
		
		AttenRelRef gmpeRef = AttenRelRef.BSSA_2014;
		double period = 0.2d;
		
		Location siteLoc = new Location(38.3210, -121.4530);
		double[] vs30Vals = {
				139.3,
				171.8,
				199.4,
				228.2,
				261.2,
				303.2,
				373.9
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
		
		List<Site> sites = new ArrayList<>();
		ScalarIMR gmpe = gmpeRef.instance(null);
		gmpe.setParamDefaults();
		for (double vs30 : vs30Vals) {
			Site site = new Site(siteLoc);
			for (Parameter<?> param : gmpe.getSiteParams()) {
				param = (Parameter<?>) param.clone();
				if (param.getName().equals(Vs30_Param.NAME))
					((WarningDoubleParameter)param).setValueIgnoreWarning(vs30);
				site.addParameter(param);
			}
			sites.add(site);
		}
		
		SolutionLogicTree slt = SolutionLogicTree.load(new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2021_11_30-u3_branches-orig_calcs-5h/results.zip"));
		LogicTree<?> tree = slt.getLogicTree();
//		tree = tree.sample(16, false);
		
		File outputDir = new File("/home/kevin/OpenSHA/UCERF3/2022_09-tim-jon-calcs");
		
		int threads = 16;
		ExecutorService exec = Executors.newFixedThreadPool(threads);
		
		List<Future<DiscretizedFunc[]>> futures = new ArrayList<>();
		List<Double> weights = new ArrayList<>();
		for (LogicTreeBranch<?> branch : tree) {
			futures.add(exec.submit(new CurveCalcCallable(slt, branch, sites, gmpeRef, period, xVals)));
			weights.add(tree.getBranchWeight(branch));
		}
		
		List<DiscretizedFunc[]> curves = new ArrayList<>();
		
		for (Future<DiscretizedFunc[]> future : futures) {
			try {
				curves.add(future.get());
			} catch (InterruptedException | ExecutionException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
			System.out.println("DONE branch "+curves.size()+"/"+futures.size());
		}
		
		exec.shutdown();
		
		for (int s=0; s<sites.size(); s++) {
			double vs30 = vs30Vals[s];
			CSVFile<String> csv = new CSVFile<>(true);
			List<String> header = new ArrayList<>();
			header.add("Branch Index");
			header.add("Branch Weight");
			for (Point2D pt : xVals)
				header.add((float)pt.getX()+"");
			csv.addLine(header);
			
			ArbDiscrEmpiricalDistFunc[] dists = new ArbDiscrEmpiricalDistFunc[xVals.size()];
			for (int i=0; i<dists.length; i++)
				dists[i] = new ArbDiscrEmpiricalDistFunc();
			
			for (int b=0; b<curves.size(); b++) {
				DiscretizedFunc[] branchCurves = curves.get(b);
				double weight = weights.get(b);
				DiscretizedFunc curve = branchCurves[s];
				
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
			List<String> medianLine = new ArrayList<>();
			medianLine.add("MEDIAN"); medianLine.add("");
			
			for (int i=0; i<dists.length; i++) {
				meanLine.add(dists[i].getMean()+"");
				medianLine.add(dists[i].getInterpolatedFractile(00.5d)+"");
			}
			csv.addLine(meanLine);
			csv.addLine(medianLine);
			
			csv.writeToFile(new File(outputDir, "curves_vs30_"+(float)vs30+".csv"));
		}
	}
	
	private static class CurveCalcCallable implements Callable<DiscretizedFunc[]> {
		
		private SolutionLogicTree slt;
		private LogicTreeBranch<?> branch;
		private List<Site> sites;
		private AttenRelRef gmpeRef;
		private double period;
		private DiscretizedFunc xVals;

		public CurveCalcCallable(SolutionLogicTree slt, LogicTreeBranch<?> branch, List<Site> sites,
				AttenRelRef gmpeRef, double period, DiscretizedFunc xVals) {
			this.slt = slt;
			this.branch = branch;
			this.sites = sites;
			this.gmpeRef = gmpeRef;
			this.period = period;
			this.xVals = xVals;
		}

		@Override
		public DiscretizedFunc[] call() throws Exception {
			ScalarIMR gmpe = gmpeRef.instance(null);
			
			gmpe.setParamDefaults();
			gmpe.setIntensityMeasure(SA_Param.NAME);
			SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), period);
			
			FaultSystemSolution sol = slt.forBranch(branch);
			FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
			
			erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.INCLUDE);
			erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_PREF_BLEND);
			erf.setParameter(HistoricOpenIntervalParam.NAME, 146d);
			erf.setParameter(BPTAveragingTypeParam.NAME, BPTAveragingTypeOptions.AVE_RI_AVE_NORM_TIME_SINCE);
			erf.getTimeSpan().setDuration(1d);
			erf.getTimeSpan().setStartTime(2021);
			erf.updateForecast();
			
			DiscretizedFunc[] ret = new DiscretizedFunc[sites.size()];
			
			HazardCurveCalculator calc = new HazardCurveCalculator();
			
			for (int s=0; s<sites.size(); s++) {
				Site site = sites.get(s);
				
				DiscretizedFunc logCurve = new ArbitrarilyDiscretizedFunc();
				for (Point2D pt : xVals)
					logCurve.set(Math.log(pt.getX()), 0d);
				
				calc.getHazardCurve(logCurve, site, gmpe, erf);
				
				DiscretizedFunc curve = new ArbitrarilyDiscretizedFunc();
				for (int i=0; i<xVals.size(); i++)
					curve.set(xVals.getX(i), logCurve.getY(i));
				ret[s] = curve;
			}
			
			return ret;
		}
		
	}

}
