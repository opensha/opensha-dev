package scratch.kevin.ucerf3;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.dom4j.DocumentException;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;

import com.google.common.collect.Lists;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.utils.FaultSystemIO;

public class TimeDepNoDataTest {

	public static void main(String[] args) throws IOException, DocumentException {
		FaultSystemSolution sol = FaultSystemIO.loadSol(new File("/home/kevin/workspace/OpenSHA/"
				+ "dev/scratch/UCERF3/data/scratch/InversionSolutions/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		double duration = 1;
		FaultSystemSolutionERF depERF = new FaultSystemSolutionERF(sol);
		depERF.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_BPT);
		depERF.getTimeSpan().setDuration(duration);
		depERF.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
		System.out.println("Updating time dep forecast");
		depERF.updateForecast();
		System.out.println("Done updating time dep forecast");
		
		FaultSystemSolutionERF indepERF = new FaultSystemSolutionERF(sol);
		indepERF.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
		indepERF.getTimeSpan().setDuration(duration);
		System.out.println("Updating time indep forecast");
		indepERF.updateForecast();
		System.out.println("Done updating time indep forecast");
		
		double[] depRates = new double[depERF.getNumSources()];
		double[] indepRates = new double[depERF.getNumSources()];
		
		for (int sourceID=0; sourceID<depERF.getNumSources(); sourceID++) {
			depRates[sourceID] = depERF.getSource(sourceID).computeTotalProb();
			indepRates[sourceID] = indepERF.getSource(sourceID).computeTotalProb();
		}
		
		DefaultXY_DataSet ratioData = new DefaultXY_DataSet(indepRates, depRates);
		
		double max = StatUtils.max(depRates);
		max = Math.max(max, StatUtils.max(indepRates));
		double min = StatUtils.min(depRates);
		min = Math.min(min, StatUtils.min(indepRates));
		min = Math.max(1e-10, min);
		ArbitrarilyDiscretizedFunc eventRatio = new ArbitrarilyDiscretizedFunc();
		eventRatio.set(0d, 0d);
		eventRatio.set(min, min);
		eventRatio.set(max, max);
		
		List<XY_DataSet> elems = Lists.newArrayList();
		elems.add(ratioData);
		elems.add(eventRatio);
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 1f, Color.BLACK));
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		
		GraphWindow gw = new GraphWindow(elems, "Time Dep vs Indep with no last event data", chars);
		gw.setX_AxisLabel("Indep Probs");
		gw.setY_AxisLabel("Dep Probs");
		Range range = new Range(min, max);
		gw.setAxisRange(range, range);
		gw.setSize(1000, 800);
		gw.setXLog(true);
		gw.setYLog(true);
	}

}
