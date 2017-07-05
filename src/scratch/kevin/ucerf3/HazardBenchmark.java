package scratch.kevin.ucerf3;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.TimeUnit;

import org.dom4j.DocumentException;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.param.Parameter;
import org.opensha.nshmp2.imr.NSHMP08_WUS;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.hazardMap.HazardCurveSetCalculator;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.gcim.ui.infoTools.IMT_Info;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;

import com.google.common.base.Stopwatch;
import com.google.common.collect.Lists;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.utils.FaultSystemIO;

public class HazardBenchmark {
	
	private static void makeSubSectArticulationsHist(FaultSystemSolution fss) {
		HistogramFunction hist = new HistogramFunction(0d, 10, 1d);
		
		for (FaultSectionPrefData sect : fss.getRupSet().getFaultSectionDataList()) {
			int num = sect.getFaultTrace().size()-2;
			hist.add((double)num, 1d);
		}
		
		ArrayList<DiscretizedFunc> funcs = Lists.newArrayList();
		funcs.add(hist);
		ArrayList<PlotCurveCharacterstics> chars = Lists.newArrayList(
				new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));
		
		new GraphWindow(funcs, "Sub Sect Bends Histogram", chars);
	}

	/**
	 * @param args
	 * @throws DocumentException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException, DocumentException {
		Site site = new Site(new Location(34.055, -118.2467));
		ScalarIMR imr = new NSHMP08_WUS();
		imr.setParamDefaults();
		for (Parameter<?> param : imr.getSiteParams())
			site.addParameter(param);
		
		imr.setIntensityMeasure(PGA_Param.NAME);
		
		InversionFaultSystemSolution fss = FaultSystemIO.loadInvSol(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/" +
						"scratch/InversionSolutions/FM3_1_ZENG_Shaw09Mod_DsrTap_" +
						"CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_mean_sol.zip"));
		
		makeSubSectArticulationsHist(fss);
		
		ERF erf = new FaultSystemSolutionERF(fss);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.ONLY);
		erf.updateForecast();
		
		ArbitrarilyDiscretizedFunc xValues = IMT_Info.getUSGS_PGA_Function();
		ArbitrarilyDiscretizedFunc curve = HazardCurveSetCalculator.getLogFunction(xValues);
		
		HazardCurveCalculator calc = new HazardCurveCalculator();
		
		Stopwatch watch = Stopwatch.createStarted();
		calc.getHazardCurve(curve, site, imr, erf);
		watch.stop();
		System.out.println("Took "+(watch.elapsed(TimeUnit.MILLISECONDS) / 1000d)+" secs to calculate");
		watch.reset();
		watch.start();
		calc.getHazardCurve(curve, site, imr, erf);
		watch.stop();
		System.out.println("Took "+(watch.elapsed(TimeUnit.MILLISECONDS) / 1000d)+" secs to calculate again");
		watch.reset();
		watch.start();
		calc.getHazardCurve(curve, site, imr, erf);
		watch.stop();
		System.out.println("Took "+(watch.elapsed(TimeUnit.MILLISECONDS) / 1000d)+" secs to calculate again");
	}

}
