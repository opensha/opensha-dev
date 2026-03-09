package scratch.ned.nshm23;

import java.awt.Color;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.jfree.data.Range;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.calc.ERF_Calculator;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships_StableContinental;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.AnalysisRegions;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.SeismicityRegions;
import org.opensha.sha.faultSurface.CompoundSurface;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import gov.usgs.earthquake.nshmp.model.NshmErf;
import gov.usgs.earthquake.nshmp.model.NshmSurface;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.ned.FSS_Inversion2019.PlottingUtils;

public class MiscPlots {	
	
	
	// NSHM23_ScalingRelationships 
	
	
	
	public static void makeSlipLengthPlot(double downDipWidth, int maxLength, boolean saveFiles) {
		
		
		// test
		double lengthKm2 = 600;
		double length2 = lengthKm2*1e3;
		double area2 = length2*downDipWidth*1e3;
		double rake2 = Double.NaN;
		NSHM23_ScalingRelationships test = NSHM23_ScalingRelationships.LOGA_C4p2_SQRT_LEN;
		NSHM23_ScalingRelationships test2 = NSHM23_ScalingRelationships.LOGA_C4p1_SQRT_LEN;
		System.out.println(test.getAveSlip(area2, length2, downDipWidth*1e3, downDipWidth*1e3, rake2));
		System.out.println(test2.getAveSlip(area2, length2, downDipWidth*1e3, downDipWidth*1e3, rake2));


		

		
		ArbitrarilyDiscretizedFunc ellA_func = new ArbitrarilyDiscretizedFunc();
		ellA_func.setName("LogA+4.1, From Moment");
		ArbitrarilyDiscretizedFunc sh09_funcMod = new ArbitrarilyDiscretizedFunc();
		sh09_funcMod.setName("Width Limited, From Moment");
		ArbitrarilyDiscretizedFunc ellB_func = new ArbitrarilyDiscretizedFunc();
		ellB_func.setName("LogA+4.2, From Moment");
		ArbitrarilyDiscretizedFunc ellC_func = new ArbitrarilyDiscretizedFunc();
		ellC_func.setName("LogA+4.3, From Moment");
		ArbitrarilyDiscretizedFunc hb_func = new ArbitrarilyDiscretizedFunc();
		hb_func.setName("HANKS_BAKUN_08, From Moment");
		ArbitrarilyDiscretizedFunc ellA_sqrtL_func = new ArbitrarilyDiscretizedFunc();
		ellA_sqrtL_func.setName("LogA+4.1, Sqrt Length");
		ArbitrarilyDiscretizedFunc ellB_sqrtL_func = new ArbitrarilyDiscretizedFunc();
		ellB_sqrtL_func.setName("LogA+4.2, Sqrt Length");
		ArbitrarilyDiscretizedFunc sh12_csd_func = new ArbitrarilyDiscretizedFunc();
		sh12_csd_func.setName("Width Limited, Const Stress Drop");
		
		
		NSHM23_ScalingRelationships ellA = NSHM23_ScalingRelationships.LOGA_C4p1;
		NSHM23_ScalingRelationships sh09_Mod = NSHM23_ScalingRelationships.WIDTH_LIMITED;
		NSHM23_ScalingRelationships ellB = NSHM23_ScalingRelationships.LOGA_C4p2;
		NSHM23_ScalingRelationships ellC = NSHM23_ScalingRelationships.LOGA_C4p3;
		ScalingRelationships hb = ScalingRelationships.HANKS_BAKUN_08;
		NSHM23_ScalingRelationships ellA_sqrtL = NSHM23_ScalingRelationships.LOGA_C4p1_SQRT_LEN;
		NSHM23_ScalingRelationships sh12_csd = NSHM23_ScalingRelationships.WIDTH_LIMITED_CSD;
		NSHM23_ScalingRelationships ellB_sqrtL = NSHM23_ScalingRelationships.LOGA_C4p2_SQRT_LEN;

		
		// log10 area from 1 to 5
    	for(int i=1; i<=maxLength; i++) {
    		double lengthKm = (double)i;
    		double length = lengthKm*1e3;
    		double area = length*downDipWidth*1e3;
    		double rake = Double.NaN;
    		ellA_func.set(lengthKm,ellA.getAveSlip(area, length, downDipWidth*1e3, downDipWidth*1e3, rake));
    		sh09_funcMod.set(lengthKm,sh09_Mod.getAveSlip(area, length, downDipWidth*1e3, downDipWidth*1e3, rake));
    		ellB_func.set(lengthKm,ellB.getAveSlip(area, length, downDipWidth*1e3, downDipWidth*1e3, rake));
    		ellC_func.set(lengthKm,ellC.getAveSlip(area, length, downDipWidth*1e3, downDipWidth*1e3, rake));
    		hb_func.set(lengthKm,hb.getAveSlip(area, length, downDipWidth*1e3, downDipWidth*1e3, rake));
    		ellA_sqrtL_func.set(lengthKm,ellA_sqrtL.getAveSlip(area, length, downDipWidth*1e3, downDipWidth*1e3, rake));
    		ellB_sqrtL_func.set(lengthKm,ellB_sqrtL.getAveSlip(area, length, downDipWidth*1e3, downDipWidth*1e3, rake));
    		sh12_csd_func.set(lengthKm,sh12_csd.getAveSlip(area, length, downDipWidth*1e3, downDipWidth*1e3, rake));
    	}
    	
    	ArrayList<ArbitrarilyDiscretizedFunc> funcs = new ArrayList<ArbitrarilyDiscretizedFunc>();
    	funcs.add(sh09_funcMod);
    	funcs.add(ellA_func);
    	funcs.add(ellB_func);
    	funcs.add(hb_func);
    	funcs.add(sh12_csd_func);
//    	funcs.add(ellA_sqrtL_func);
    	funcs.add(ellC_func);
    	funcs.add(ellB_sqrtL_func);
    	
    	ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.GREEN));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.MAGENTA));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.GRAY));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.CYAN)); // this is same as the one above

    	
		GraphWindow graph = new GraphWindow(funcs, "Slip-Length Relationships; DDW="+downDipWidth+" km", plotChars); 
		graph.setX_AxisLabel("Length (km)");
		graph.setY_AxisLabel("Slip (m)");
		graph.setPlotLabelFontSize(18);
		graph.setAxisLabelFontSize(18);
		graph.setTickLabelFontSize(16);
		
		if(saveFiles) {
			try {
				graph.saveAsPDF("slipLengthScaling2023_Plot.pdf");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

	}

	
	
	/**
	 * This assumes no aseismicity
	 * @param saveFiles
	 */
	public static void makeMagAreaPlot(boolean saveFiles) {
		
		double downDipWidth=11;	// orig down-dip width equals reduced
		
		ArbitrarilyDiscretizedFunc sh09mod_func = new ArbitrarilyDiscretizedFunc();
		sh09mod_func.setName("Width Limited"+downDipWidth);
		ArbitrarilyDiscretizedFunc ellA_func = new ArbitrarilyDiscretizedFunc();
		ellA_func.setName("LogA+4.1");
		ArbitrarilyDiscretizedFunc ellB_func = new ArbitrarilyDiscretizedFunc();
		ellB_func.setName("LogA+4.2");
		ArbitrarilyDiscretizedFunc ellC_func = new ArbitrarilyDiscretizedFunc();
		ellB_func.setName("LogA+4.3");
		ArbitrarilyDiscretizedFunc hb_func = new ArbitrarilyDiscretizedFunc();
		hb_func.setName("HANKS_BAKUN_08");
		
		NSHM23_ScalingRelationships sh09mod = NSHM23_ScalingRelationships.WIDTH_LIMITED;
		NSHM23_ScalingRelationships ellA = NSHM23_ScalingRelationships.LOGA_C4p1;
		NSHM23_ScalingRelationships ellB = NSHM23_ScalingRelationships.LOGA_C4p2;
		NSHM23_ScalingRelationships ellC = NSHM23_ScalingRelationships.LOGA_C4p3;
		ScalingRelationships hb = ScalingRelationships.HANKS_BAKUN_08;
		
		// log10 area from 1 to 5
    	for(int i=50; i<=20000; i+=10) {
    		double area = (double)i;
    		double rake = Double.NaN;
     		sh09mod_func.set(area,sh09mod.getMag(area*1e6,1e3*area/downDipWidth, downDipWidth*1e3, downDipWidth*1e3, rake));
    		ellA_func.set(area,ellA.getMag(area*1e6,1e3*area/downDipWidth, downDipWidth*1e3, downDipWidth*1e3, rake));
    		ellB_func.set(area,ellB.getMag(area*1e6,1e3*area/downDipWidth, downDipWidth*1e3, downDipWidth*1e3, rake));
    		ellC_func.set(area,ellC.getMag(area*1e6,1e3*area/downDipWidth, downDipWidth*1e3, downDipWidth*1e3, rake));
    		hb_func.set(area,hb.getMag(area*1e6,1e3*area/downDipWidth, downDipWidth*1e3, downDipWidth*1e3, rake));
    	}
    	
    	ArrayList<ArbitrarilyDiscretizedFunc> funcs = new ArrayList<ArbitrarilyDiscretizedFunc>();
    	funcs.add(sh09mod_func);
    	funcs.add(ellA_func);
    	funcs.add(ellB_func);
    	funcs.add(ellC_func);
    	funcs.add(hb_func);
    	
    	ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.GRAY));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.GREEN));

    	
		GraphWindow graph = new GraphWindow(funcs, "Mag-Area Relationships",plotChars); 
		graph.setX_AxisLabel("Area (km-sq)");
		graph.setY_AxisLabel("Magnitude");
		graph.setXLog(true);
		graph.setX_AxisRange(50, 2e4);
		graph.setY_AxisRange(5, 9);
		graph.setPlotLabelFontSize(18);
		graph.setAxisLabelFontSize(18);
		graph.setTickLabelFontSize(16);
		
		if(saveFiles) {
			try {
				graph.saveAsPDF("magAreaScaling2023_Plot.pdf");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	
	
	/**
	 * This assumes no aseismicity
	 * @param saveFiles
	 */
	public static void makeMagAreaPlot_StableContinental(boolean saveFiles) {
		
		double downDipWidth=15;	// orig down-dip width equals reduced
		
		ArbitrarilyDiscretizedFunc widthLimited_func = new ArbitrarilyDiscretizedFunc();
		widthLimited_func.setName("Width Limited"+downDipWidth);
		ArbitrarilyDiscretizedFunc logA_4pt2_func = new ArbitrarilyDiscretizedFunc();
		logA_4pt2_func.setName("LogA+4.2");
		ArbitrarilyDiscretizedFunc logA_4pt3_func = new ArbitrarilyDiscretizedFunc();
		logA_4pt3_func.setName("LogA+4.3");
		ArbitrarilyDiscretizedFunc wcLegth_func = new ArbitrarilyDiscretizedFunc();
		wcLegth_func.setName("Wells & Coppersmith Length");
		
		
		NSHM23_ScalingRelationships_StableContinental widthLimited = NSHM23_ScalingRelationships_StableContinental.WIDTH_LIMITED;
		NSHM23_ScalingRelationships_StableContinental logA_4pt2 = NSHM23_ScalingRelationships_StableContinental.LOGA_C4p2;
		NSHM23_ScalingRelationships_StableContinental logA_4pt3 = NSHM23_ScalingRelationships_StableContinental.LOGA_C4p3;
		WC1994_MagLengthRelationship wc = new WC1994_MagLengthRelationship();
		
		String infoString = "AreaSource\tlength\tarea\tlogA+4.2\tlogA+4.3\twidthLimited\twc94-Length\n";
		
		// test
		double[] lengthArray = {68, 21, 84, 53, 35, 121};
		String[] names = {"Cental VA Regional", "Cental VA Local", "Crowleys Ridge (South)", "Crowleys Ridge (West)", "Joyner Ridge", "Saline River"};
		for(int i=0;i<6; i++) {
			double area = lengthArray[i]*downDipWidth;
			double rake = Double.NaN;
			double mag1 = logA_4pt2.getMag(area*1e6,1e3*area/downDipWidth, downDipWidth*1e3, downDipWidth*1e3, rake);
			double mag2 = logA_4pt3.getMag(area*1e6,1e3*area/downDipWidth, downDipWidth*1e3, downDipWidth*1e3, rake);
			double mag3 = widthLimited.getMag(area*1e6,1e3*area/downDipWidth, downDipWidth*1e3, downDipWidth*1e3, rake);
			double mag4 = wc.getMedianMag(area/downDipWidth);
			infoString += names[i]+"\t"+lengthArray[i]+"\t"+(float)area+"\t"+(float)mag1+"\t"+(float)mag2+"\t"+(float)mag3+"\t"+(float)mag4+"\n";
		}
		System.out.println(infoString);			


		
		// log10 area from 1 to 5
    	for(int i=50; i<=20000; i+=10) {
    		double area = (double)i;
    		double rake = Double.NaN;
     		widthLimited_func.set(area,widthLimited.getMag(area*1e6,1e3*area/downDipWidth, downDipWidth*1e3, downDipWidth*1e3, rake));
    		logA_4pt2_func.set(area,logA_4pt2.getMag(area*1e6,1e3*area/downDipWidth, downDipWidth*1e3, downDipWidth*1e3, rake));
    		logA_4pt3_func.set(area,logA_4pt3.getMag(area*1e6,1e3*area/downDipWidth, downDipWidth*1e3, downDipWidth*1e3, rake));
    		wc.getMedianMag(area/downDipWidth);
    		wcLegth_func.set(area,wc.getMedianMag(area/downDipWidth));
    	}
    	
    	ArrayList<ArbitrarilyDiscretizedFunc> funcs = new ArrayList<ArbitrarilyDiscretizedFunc>();
    	funcs.add(widthLimited_func);
    	funcs.add(logA_4pt2_func);
    	funcs.add(logA_4pt3_func);
    	funcs.add(wcLegth_func);
    	
    	ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.GREEN));

    	
		GraphWindow graph = new GraphWindow(funcs, "Stable Continental Mag-Area Relationships (Shaw, 2023)",plotChars); 
		graph.setX_AxisLabel("Area (km-sq)");
		graph.setY_AxisLabel("Magnitude");
		graph.setXLog(true);
		graph.setX_AxisRange(50, 2e4);
		graph.setY_AxisRange(5, 9);
		graph.setPlotLabelFontSize(18);
		graph.setAxisLabelFontSize(18);
		graph.setTickLabelFontSize(16);
		
		if(saveFiles) {
			try {
				graph.saveAsPDF("magAreaScaling2023_stableContinental_Plot.pdf");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	
	/**
	 * This assumes no aseismicity
	 * @param saveFiles
	 */
	public static void makeRegMFD_Plots(boolean saveFiles) {
		
//		ceus_full.in
//		mean rate for M>= 5 0.433 with b= 0.94
//		mean - 2 *sigma rate for M>= 5 0.374 with b= 0.98
//		mean + 2 *sigma rate for M>= 5 0.508 with b= 0.9
//		ceus_gk.in
//		mean rate for M>= 5 0.579 with b= 0.8
//		mean - 2 *sigma rate for M>= 5 0.499 with b= 0.84
//		mean + 2 *sigma rate for M>= 5 0.681 with b= 0.76
//		ceus_nn.in
//		mean rate for M>= 5 0.226 with b= 0.98
//		mean - 2 *sigma rate for M>= 5 0.186 with b= 1.04
//		mean + 2 *sigma rate for M>= 5 0.282 with b= 0.92
//		ceus_r85.in
//		mean rate for M>= 5 0.409 with b= 0.94
//		mean - 2 *sigma rate for M>= 5 0.353 with b= 0.98
//		mean + 2 *sigma rate for M>= 5 0.48 with b= 0.9
//		deepn_full.in deeps_full.in
//		mean rate for M>= 5 0.116 with b= 0.71
//		mean - 2 *sigma rate for M>= 5 0.0915 with b= 0.87
//		mean + 2 *sigma rate for M>= 5 0.272 with b= 0.55
//		deepn_gk.in deeps_gk.in
//		mean rate for M>= 5 0.101 with b= 0.7
//		mean - 2 *sigma rate for M>= 5 0.0775 with b= 0.86
//		mean + 2 *sigma rate for M>= 5 0.256 with b= 0.54
//		deepn_nn.in deeps_nn.in
//		mean rate for M>= 5 0.11 with b= 0.81
//		mean - 2 *sigma rate for M>= 5 0.0806 with b= 1.01
//		mean + 2 *sigma rate for M>= 5 0.314 with b= 0.61
//		deepn_r85.in deeps_r85.in
//		mean rate for M>= 5 0.114 with b= 0.71
//		mean - 2 *sigma rate for M>= 5 0.0889 with b= 0.85
//		mean + 2 *sigma rate for M>= 5 0.268 with b= 0.57
//		imw_full.in    pnw_full.in    ucerf3_full.in other_full.in 
//		mean rate for M>= 5 12 with b= 0.84
//		mean - 2 *sigma rate for M>= 5 11.3 with b= 0.88
//		mean + 2 *sigma rate for M>= 5 13.2 with b= 0.8
//		imw_gk.in    pnw_gk.in    ucerf3_gk.in other_gk.in 
//		mean rate for M>= 5 7.18 with b= 0.71
//		mean - 2 *sigma rate for M>= 5 6.7 with b= 0.73
//		mean + 2 *sigma rate for M>= 5 8 with b= 0.69
//		imw_nn.in    pnw_nn.in    ucerf3_nn.in other_nn.in 
//		mean rate for M>= 5 8.02 with b= 0.8
//		mean - 2 *sigma rate for M>= 5 7.45 with b= 0.84
//		mean + 2 *sigma rate for M>= 5 9.01 with b= 0.76
//		imw_r85.in    pnw_r85.in    ucerf3_r85.in other_r85.in 
//		mean rate for M>= 5 9.74 with b= 0.82
//		mean - 2 *sigma rate for M>= 5 9.13 with b= 0.86
//		mean + 2 *sigma rate for M>= 5 10.8 with b= 0.78
//		ucerf3_full.in
//		mean rate for M>= 5 8.25 with b= 0.84
//		mean - 2 *sigma rate for M>= 5 7.55 with b= 0.88
//		mean + 2 *sigma rate for M>= 5 9.13 with b= 0.8
//		ucerf3_gk.in
//		mean rate for M>= 5 4.85 with b= 0.71
//		mean - 2 *sigma rate for M>= 5 4.35 with b= 0.73
//		mean + 2 *sigma rate for M>= 5 5.47 with b= 0.69
//		ucerf3_nn.in
//		mean rate for M>= 5 5.4 with b= 0.8
//		mean - 2 *sigma rate for M>= 5 4.79 with b= 0.84
//		mean + 2 *sigma rate for M>= 5 6.17 with b= 0.76
//		ucerf3_r85.in
//		mean rate for M>= 5 6.44 with b= 0.82
//		mean - 2 *sigma rate for M>= 5 5.79 with b= 0.86
//		mean + 2 *sigma rate for M>= 5 7.25 with b= 0.78
		
		GutenbergRichterMagFreqDist gr_WUS = new GutenbergRichterMagFreqDist(5.05, 8.95, 40);
		GutenbergRichterMagFreqDist gr_WUS_low = new GutenbergRichterMagFreqDist(5.05, 8.95, 40);
		GutenbergRichterMagFreqDist gr_WUS_high = new GutenbergRichterMagFreqDist(5.05, 8.95, 40);
		GutenbergRichterMagFreqDist gr_WUS_GKdecl = new GutenbergRichterMagFreqDist(5.05, 8.95, 40);
		GutenbergRichterMagFreqDist gr_WUS_NN = new GutenbergRichterMagFreqDist(5.05, 8.95, 40);
		GutenbergRichterMagFreqDist gr_WUS_Reas = new GutenbergRichterMagFreqDist(5.05, 8.95, 40);


		gr_WUS.setAllButTotMoRate(5.05, 8.45, 12.0, 0.84);
		gr_WUS_low.setAllButTotMoRate(5.05, 8.45, 11.3, 0.88);
		gr_WUS_high.setAllButTotMoRate(5.05, 8.45, 13.2, 0.80);
		gr_WUS_GKdecl.setAllButTotMoRate(5.05, 8.45, 7.18, 0.71);
		gr_WUS_GKdecl.setName("GK");
		gr_WUS_NN.setAllButTotMoRate(5.05, 8.45, 8.02, 0.8);
		gr_WUS_NN.setName("NN");
		gr_WUS_Reas.setAllButTotMoRate(5.05, 8.45, 9.74, 0.82);
		gr_WUS_NN.setName("Reas");

		ArrayList<EvenlyDiscretizedFunc> funcs = new ArrayList<EvenlyDiscretizedFunc>();
//    	funcs.add(gr_WUS);
//    	funcs.add(gr_WUS_low);
//    	funcs.add(gr_WUS_high);
    	funcs.add(gr_WUS.getCumRateDistWithOffset());
    	funcs.add(gr_WUS_low.getCumRateDistWithOffset());
    	funcs.add(gr_WUS_high.getCumRateDistWithOffset());
    	funcs.add(gr_WUS_GKdecl.getCumRateDistWithOffset());
    	funcs.add(gr_WUS_NN.getCumRateDistWithOffset());
//    	funcs.add(gr_WUS_Reas.getCumRateDistWithOffset());
    	
    	ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, null, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, null, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.GREEN));

    	
		GraphWindow graph = new GraphWindow(funcs, "WUS MFD",plotChars); 
		graph.setX_AxisLabel("Magnitude");
		graph.setY_AxisLabel("Cumulative Rate (per yr)");
		graph.setYLog(true);
//		graph.setX_AxisRange(50, 2e4);
		graph.setY_AxisRange(0.00001,20);
		graph.setPlotLabelFontSize(18);
		graph.setAxisLabelFontSize(18);
		graph.setTickLabelFontSize(16);
		
		if(saveFiles) {
			try {
				graph.saveAsPDF("WUS_MFD_Plot.pdf");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		GutenbergRichterMagFreqDist gr_CEUS = new GutenbergRichterMagFreqDist(5.05, 8.95, 40);
		GutenbergRichterMagFreqDist gr_CEUS_low = new GutenbergRichterMagFreqDist(5.05, 8.95, 40);
		GutenbergRichterMagFreqDist gr_CEUS_high = new GutenbergRichterMagFreqDist(5.05, 8.95, 40);
		GutenbergRichterMagFreqDist gr_CEUS_GKdecl = new GutenbergRichterMagFreqDist(5.05, 8.95, 40);
		GutenbergRichterMagFreqDist gr_CEUS_NN = new GutenbergRichterMagFreqDist(5.05, 8.95, 40);
		GutenbergRichterMagFreqDist gr_CEUS_Reas = new GutenbergRichterMagFreqDist(5.05, 8.95, 40);

		
		gr_CEUS.setAllButTotMoRate(5.05, 8.45, 0.43,0.94);
		gr_CEUS_low.setAllButTotMoRate(5.05, 8.45, 0.37,0.98);
		gr_CEUS_high.setAllButTotMoRate(5.05, 8.45, 0.51,0.90);
		gr_CEUS_GKdecl.setAllButTotMoRate(5.05, 8.45, 0.579, 0.8);
		gr_CEUS_GKdecl.setName("GK");
		gr_CEUS_NN.setAllButTotMoRate(5.05, 8.45, 0.226, 0.98);
		gr_CEUS_NN.setName("NN");
		gr_CEUS_Reas.setAllButTotMoRate(5.05, 8.45, 0.409, 0.94);
		gr_CEUS_Reas.setName("Reas");

    	
    	ArrayList<EvenlyDiscretizedFunc> funcsCEUS = new ArrayList<EvenlyDiscretizedFunc>();
//    	funcsCEUS.add(gr_CEUS);
//    	funcsCEUS.add(gr_CEUS_low);
//    	funcsCEUS.add(gr_CEUS_high);
    	funcsCEUS.add(gr_CEUS.getCumRateDistWithOffset());
    	funcsCEUS.add(gr_CEUS_low.getCumRateDistWithOffset());
    	funcsCEUS.add(gr_CEUS_high.getCumRateDistWithOffset());
    	funcsCEUS.add(gr_CEUS_GKdecl.getCumRateDistWithOffset());
    	funcsCEUS.add(gr_CEUS_NN.getCumRateDistWithOffset());
//    	funcsCEUS.add(gr_CEUS_Reas.getCumRateDistWithOffset());
    	
 		GraphWindow graphCEUS = new GraphWindow(funcsCEUS, "CEUS MFD",plotChars); 
		graphCEUS.setX_AxisLabel("Magnitude");
		graphCEUS.setY_AxisLabel("Cumulative Rate (per yr)");
		graphCEUS.setYLog(true);
//		graphCEUS.setX_AxisRange(50, 2e4);
		graphCEUS.setY_AxisRange(0.00001,20);
		graphCEUS.setPlotLabelFontSize(18);
		graphCEUS.setAxisLabelFontSize(18);
		graphCEUS.setTickLabelFontSize(16);
		
		if(saveFiles) {
			try {
				graphCEUS.saveAsPDF("CEUS_MFD_Plot.pdf");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

	}
	
	public static void listSrcClassesForERF(NshmErf erf) {
		
		ArrayList<String> stringSrcTypesList = new ArrayList<String>();
		ArrayList<String> stringSurfTypesList = new ArrayList<String>();
		int totNum=0;

		for (int s = 0; s < erf.getNumSources(); ++s) {
			ProbEqkSource source = erf.getSource(s);
			if(!stringSrcTypesList.contains(source.getClass().toString()))
				stringSrcTypesList.add(source.getClass().toString());
			
			for (int r = 0; r < source.getNumRuptures(); ++r) {
				ProbEqkRupture rupture = source.getRupture(r);
				RuptureSurface surf = rupture.getRuptureSurface();
				if(!stringSurfTypesList.contains(surf.getClass().toString()))
					stringSurfTypesList.add(surf.getClass().toString());
				int numPts = surf.getEvenlyDiscritizedListOfLocsOnSurface().size();
				if(numPts==1)
					totNum+=1;
			}
		}
		System.out.println(stringSrcTypesList);
		System.out.println("\n"+stringSurfTypesList);
		System.out.println("\ntotNum="+totNum);
		return;
	}

	
	
	public static void OLDmakeCEUS_ModelMFD_Plots(boolean saveFiles) {
		
//		File outputDir = new File("/tmp/wrapper_tests");
//		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		Region region=null;
		try {
//			region = AnalysisRegions.CONUS_EAST.load();
			region = SeismicityRegions.CONUS_WEST.load();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	
		// Full MFDs
		IncludeBackgroundOption bgOption = IncludeBackgroundOption.INCLUDE;
		
		Set<TectonicRegionType> trts = EnumSet.of(TectonicRegionType.ACTIVE_SHALLOW,TectonicRegionType.STABLE_SHALLOW);
		Path erfPath = Path.of("/Users/field/nshm-haz_data/nshm-conus-5.3.0");
		NshmErf erf = new NshmErf(erfPath, trts, bgOption);
		erf.getTimeSpan().setDuration(1.0);
		erf.updateForecast();
		SummedMagFreqDist mfd2018_both = ERF_Calculator.getMagFreqDistInRegion(erf, region, 5.05,40,0.1, true);
		mfd2018_both.setName("2018 CEUS MFD full model");
		mfd2018_both.setInfo("rate >= m5 = "+(float)mfd2018_both.getTotalIncrRate());
		
		
		erfPath = Path.of("/Users/field/nshm-haz_data/nshm-conus-6.a.6");
		erf = new NshmErf(erfPath, trts, bgOption);
		erf.getTimeSpan().setDuration(1.0);
		erf.updateForecast();
		SummedMagFreqDist mfd2023_both = ERF_Calculator.getMagFreqDistInRegion(erf, region, 5.05,40,0.1, true);
		mfd2023_both.setName("2023 CEUS MFD full model");
		mfd2023_both.setInfo("rate >= m5 = "+(float)mfd2023_both.getTotalIncrRate());
		
		
		// Background MFDs
		bgOption = IncludeBackgroundOption.ONLY;
		trts = EnumSet.of(TectonicRegionType.ACTIVE_SHALLOW,TectonicRegionType.STABLE_SHALLOW);
		erfPath = Path.of("/Users/field/nshm-haz_data/nshm-conus-5.3.0");
		erf = new NshmErf(erfPath, trts, bgOption);
		erf.getTimeSpan().setDuration(1.0);
		erf.updateForecast();
		SummedMagFreqDist mfd2018_bg = ERF_Calculator.getMagFreqDistInRegion(erf, region, 5.05,40,0.1, true);
		mfd2018_bg.setName("2018 CEUS MFD background only");
		mfd2018_bg.setInfo("rate >= m5 = "+(float)mfd2018_bg.getTotalIncrRate());
		
		
		erfPath = Path.of("/Users/field/nshm-haz_data/nshm-conus-6.a.6");
		erf = new NshmErf(erfPath, trts, bgOption);
		erf.getTimeSpan().setDuration(1.0);
		erf.updateForecast();
		SummedMagFreqDist mfd2023_bg = ERF_Calculator.getMagFreqDistInRegion(erf, region, 5.05,40,0.1, true);
		mfd2023_bg.setName("2023 CEUS MFD background only");
		mfd2023_bg.setInfo("rate >= m5 = "+(float)mfd2023_bg.getTotalIncrRate());

		
		// Fault MFDs
		bgOption = IncludeBackgroundOption.EXCLUDE;
		trts = EnumSet.of(TectonicRegionType.ACTIVE_SHALLOW,TectonicRegionType.STABLE_SHALLOW);
		erfPath = Path.of("/Users/field/nshm-haz_data/nshm-conus-5.3.0");
		erf = new NshmErf(erfPath, trts, bgOption);
		erf.getTimeSpan().setDuration(1.0);
		erf.updateForecast();
		SummedMagFreqDist mfd2018_faults = ERF_Calculator.getMagFreqDistInRegion(erf, region, 5.05,40,0.1, true);
		mfd2018_faults.setName("2018 CEUS MFD faults only");
		mfd2018_faults.setInfo("rate >= m5 = "+(float)mfd2018_faults.getTotalIncrRate());
		
		
		erfPath = Path.of("/Users/field/nshm-haz_data/nshm-conus-6.a.6");
		erf = new NshmErf(erfPath, trts, bgOption);
		erf.getTimeSpan().setDuration(1.0);
		erf.updateForecast();
		SummedMagFreqDist mfd2023_faults = ERF_Calculator.getMagFreqDistInRegion(erf, region, 5.05,40,0.1, true);
		mfd2023_faults.setName("2023 CEUS MFD faults only");
		mfd2023_faults.setInfo("rate >= m5 = "+(float)mfd2023_faults.getTotalIncrRate());

		
//		// TEMP
		GutenbergRichterMagFreqDist gr_CEUS = new GutenbergRichterMagFreqDist(5.05, 8.95, 40);
//		GutenbergRichterMagFreqDist gr_CEUS_low = new GutenbergRichterMagFreqDist(5.05, 8.95, 40);
//		GutenbergRichterMagFreqDist gr_CEUS_high = new GutenbergRichterMagFreqDist(5.05, 8.95, 40);
//		gr_CEUS.setAllButTotMoRate(5.05, 8.95, 0.43,0.94);
		gr_CEUS.setAllButTotMoRate(5.05, 8.95, mfd2023_both.getTotalIncrRate(),0.94);
//		gr_CEUS_low.setAllButTotMoRate(5.05, 8.45, 0.37,0.98);
//		gr_CEUS_high.setAllButTotMoRate(5.05, 8.45, 0.51,0.90);


		ArrayList<EvenlyDiscretizedFunc> funcs = new ArrayList<EvenlyDiscretizedFunc>();
    	funcs.add(mfd2018_both);
    	funcs.add(mfd2023_both);
    	funcs.add(mfd2018_bg);
    	funcs.add(mfd2023_bg);
    	funcs.add(mfd2018_faults);
    	funcs.add(mfd2023_faults);
    	funcs.add(gr_CEUS);
    	
		ArrayList<EvenlyDiscretizedFunc> funcsCumulative = new ArrayList<EvenlyDiscretizedFunc>();
		funcsCumulative.add(mfd2018_both.getCumRateDistWithOffset());
		funcsCumulative.add(mfd2023_both.getCumRateDistWithOffset());
		funcsCumulative.add(mfd2018_bg.getCumRateDistWithOffset());
		funcsCumulative.add(mfd2023_bg.getCumRateDistWithOffset());
		funcsCumulative.add(mfd2018_faults.getCumRateDistWithOffset());
		funcsCumulative.add(mfd2023_faults.getCumRateDistWithOffset());
		funcsCumulative.add(gr_CEUS.getCumRateDistWithOffset());

     	
    	ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, null, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, null, 1f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, null, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, null, 1f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, null, 1f, Color.BLACK));

    	
		GraphWindow graph = new GraphWindow(funcs, "CEUS Model MFDs",plotChars); 
		graph.setX_AxisLabel("Magnitude");
		graph.setY_AxisLabel("Incremental Rate (per yr)");
		graph.setYLog(true);
		graph.setX_AxisRange(5, 8.5);
		graph.setY_AxisRange(1e-6,1);
		graph.setPlotLabelFontSize(18);
		graph.setAxisLabelFontSize(18);
		graph.setTickLabelFontSize(16);
		
		if(saveFiles) {
			try {
				graph.saveAsPDF("CEUS_Model_IncrMFDs_Plot.pdf");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		GraphWindow graphCum = new GraphWindow(funcsCumulative, "CEUS Model MFDs",plotChars); 
		graphCum.setX_AxisLabel("Magnitude");
		graphCum.setY_AxisLabel("Cumulative Rate (per yr)");
		graphCum.setYLog(true);
		graphCum.setX_AxisRange(5, 8.5);
		graphCum.setY_AxisRange(1e-6,1);
		graphCum.setPlotLabelFontSize(18);
		graphCum.setAxisLabelFontSize(18);
		graphCum.setTickLabelFontSize(16);
		
		if(saveFiles) {
			try {
				graphCum.saveAsPDF("CEUS_Model_CumMFDs_Plot.pdf");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

	}
	
	
	
	public static void makeModelMFD_Plots(ArrayList<Region> regionList, boolean saveFiles) {
		
//		File outputDir = new File("/tmp/wrapper_tests");
//		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
			
		// Full MFDs
		IncludeBackgroundOption bgOption = IncludeBackgroundOption.INCLUDE;
		
		Set<TectonicRegionType> trts = EnumSet.of(TectonicRegionType.ACTIVE_SHALLOW,TectonicRegionType.STABLE_SHALLOW);
		Path erfPath = Path.of("/Users/field/nshm-haz_data/nshm-conus-5.3.0");
		NshmErf erf = new NshmErf(erfPath, trts, bgOption);
		erf.getTimeSpan().setDuration(1.0);
		erf.updateForecast();
		System.out.println("Done with ERF1");
		ArrayList<SummedMagFreqDist> mfd2018_bothList= new ArrayList<SummedMagFreqDist>();
		for(Region region:regionList) {
			System.out.println("\tworking on "+region.getName());
			SummedMagFreqDist mfd2018_both = ERF_Calculator.getMagFreqDistInRegion(erf, region, 5.05,40,0.1, true);
			mfd2018_both.setName("2018 full model MFD, "+region.getName());
			mfd2018_both.setInfo("rate >= m5 = "+(float)mfd2018_both.getTotalIncrRate());
			mfd2018_bothList.add(mfd2018_both);
		}
		
		
		erfPath = Path.of("/Users/field/nshm-haz_data/nshm-conus-6.a.6");
		erf = new NshmErf(erfPath, trts, bgOption);
		erf.getTimeSpan().setDuration(1.0);
		erf.updateForecast();
		System.out.println("Done with ERF2");
		ArrayList<SummedMagFreqDist> mfd2023_bothList= new ArrayList<SummedMagFreqDist>();
		for(Region region:regionList) {
			System.out.println("\tworking on "+region.getName());
			SummedMagFreqDist mfd2023_both = ERF_Calculator.getMagFreqDistInRegion(erf, region, 5.05,40,0.1, true);
			mfd2023_both.setName("2023 full model MFD, "+region.getName());
			mfd2023_both.setInfo("rate >= m5 = "+(float)mfd2023_both.getTotalIncrRate());
			mfd2023_bothList.add(mfd2023_both);
		}
		
		
		// Background MFDs
		bgOption = IncludeBackgroundOption.ONLY;
		trts = EnumSet.of(TectonicRegionType.ACTIVE_SHALLOW,TectonicRegionType.STABLE_SHALLOW);
		erfPath = Path.of("/Users/field/nshm-haz_data/nshm-conus-5.3.0");
		erf = new NshmErf(erfPath, trts, bgOption);
		erf.getTimeSpan().setDuration(1.0);
		erf.updateForecast();
		System.out.println("Done with ERF3");
		ArrayList<SummedMagFreqDist> mfd2018_bgList= new ArrayList<SummedMagFreqDist>();
		for(Region region:regionList) {
			System.out.println("\tworking on "+region.getName());
			SummedMagFreqDist mfd2018_bg = ERF_Calculator.getMagFreqDistInRegion(erf, region, 5.05,40,0.1, true);
			mfd2018_bg.setName("2018 background only MFD, "+region.getName());
			mfd2018_bg.setInfo("rate >= m5 = "+(float)mfd2018_bg.getTotalIncrRate());	
			mfd2018_bgList.add(mfd2018_bg);
		}

		
		
		erfPath = Path.of("/Users/field/nshm-haz_data/nshm-conus-6.a.6");
		erf = new NshmErf(erfPath, trts, bgOption);
		erf.getTimeSpan().setDuration(1.0);
		erf.updateForecast();
		System.out.println("Done with ERF4");
		ArrayList<SummedMagFreqDist> mfd2023_bgList= new ArrayList<SummedMagFreqDist>();
		for(Region region:regionList) {
			System.out.println("\tworking on "+region.getName());
			SummedMagFreqDist mfd2023_bg = ERF_Calculator.getMagFreqDistInRegion(erf, region, 5.05,40,0.1, true);
			mfd2023_bg.setName("2023 background only MFD, "+region.getName());
			mfd2023_bg.setInfo("rate >= m5 = "+(float)mfd2023_bg.getTotalIncrRate());	
			mfd2023_bgList.add(mfd2023_bg);
		}

		
		// Fault MFDs
		bgOption = IncludeBackgroundOption.EXCLUDE;
		trts = EnumSet.of(TectonicRegionType.ACTIVE_SHALLOW,TectonicRegionType.STABLE_SHALLOW);
		erfPath = Path.of("/Users/field/nshm-haz_data/nshm-conus-5.3.0");
		erf = new NshmErf(erfPath, trts, bgOption);
		erf.getTimeSpan().setDuration(1.0);
		erf.updateForecast();
		System.out.println("Done with ERF5");
		ArrayList<SummedMagFreqDist> mfd2018_faultsList= new ArrayList<SummedMagFreqDist>();
		for(Region region:regionList) {
			System.out.println("\tworking on "+region.getName());
			SummedMagFreqDist mfd2018_faults = ERF_Calculator.getMagFreqDistInRegion(erf, region, 5.05,40,0.1, true);
			mfd2018_faults.setName("2018 faults only MFD, "+region.getName());
			mfd2018_faults.setInfo("rate >= m5 = "+(float)mfd2018_faults.getTotalIncrRate());
			mfd2018_faultsList.add(mfd2018_faults);
		}

		
		
		erfPath = Path.of("/Users/field/nshm-haz_data/nshm-conus-6.a.6");
		erf = new NshmErf(erfPath, trts, bgOption);
		erf.getTimeSpan().setDuration(1.0);
		erf.updateForecast();
		System.out.println("Done with ERF6");
		ArrayList<SummedMagFreqDist> mfd2023_faultsList= new ArrayList<SummedMagFreqDist>();
		for(Region region:regionList) {
			System.out.println("\tworking on "+region.getName());
			SummedMagFreqDist mfd2023_faults = ERF_Calculator.getMagFreqDistInRegion(erf, region, 5.05,40,0.1, true);
			mfd2023_faults.setName("2023 faults only MFD, "+region.getName());
			mfd2023_faults.setInfo("rate >= m5 = "+(float)mfd2023_faults.getTotalIncrRate());
			mfd2023_faultsList.add(mfd2023_faults);
		}

    	ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, null, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, null, 1f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, null, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, null, 1f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, null, 1f, Color.BLACK));

		
//		// TEMP
//		GutenbergRichterMagFreqDist gr_CEUS = new GutenbergRichterMagFreqDist(5.05, 8.95, 40);
//		GutenbergRichterMagFreqDist gr_CEUS_low = new GutenbergRichterMagFreqDist(5.05, 8.95, 40);
//		GutenbergRichterMagFreqDist gr_CEUS_high = new GutenbergRichterMagFreqDist(5.05, 8.95, 40);
//		gr_CEUS.setAllButTotMoRate(5.05, 8.95, 0.43,0.94);
//		gr_CEUS.setAllButTotMoRate(5.05, 8.95, mfd2023_both.getTotalIncrRate(),0.94);
//		gr_CEUS_low.setAllButTotMoRate(5.05, 8.45, 0.37,0.98);
//		gr_CEUS_high.setAllButTotMoRate(5.05, 8.45, 0.51,0.90);

		for(int i=0;i<regionList.size();i++) {

			Region region = regionList.get(i);



			ArrayList<EvenlyDiscretizedFunc> funcs = new ArrayList<EvenlyDiscretizedFunc>();
			funcs.add(mfd2018_bothList.get(i));
			funcs.add(mfd2023_bothList.get(i));
			funcs.add(mfd2018_bgList.get(i));
			funcs.add(mfd2023_bgList.get(i));
			funcs.add(mfd2018_faultsList.get(i));
			funcs.add(mfd2023_faultsList.get(i));
			//  	funcs.add(gr_CEUS);

			ArrayList<EvenlyDiscretizedFunc> funcsCumulative = new ArrayList<EvenlyDiscretizedFunc>();
			funcsCumulative.add(mfd2018_bothList.get(i).getCumRateDistWithOffset());
			funcsCumulative.add(mfd2023_bothList.get(i).getCumRateDistWithOffset());
			funcsCumulative.add(mfd2018_bgList.get(i).getCumRateDistWithOffset());
			funcsCumulative.add(mfd2023_bgList.get(i).getCumRateDistWithOffset());
			funcsCumulative.add(mfd2018_faultsList.get(i).getCumRateDistWithOffset());
			funcsCumulative.add(mfd2023_faultsList.get(i).getCumRateDistWithOffset());
			//		funcsCumulative.add(gr_CEUS.getCumRateDistWithOffset());
			
			GraphWindow graph = new GraphWindow(funcs, region.getName(),plotChars); 
			graph.setX_AxisLabel("Magnitude");
			graph.setY_AxisLabel("Incremental Rate (per yr)");
			graph.setYLog(true);
			graph.setX_AxisRange(5, 8.5);
			graph.setY_AxisRange(1e-6,20);
			graph.setPlotLabelFontSize(18);
			graph.setAxisLabelFontSize(18);
			graph.setTickLabelFontSize(16);
			
			if(saveFiles) {
				try {
					graph.saveAsPDF(region.getName().replace(" ", "_")+"_MFD_Plot.pdf");
					graph.saveAsTXT(region.getName().replace(" ", "_")+"_MFD_Plot.txt");
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			
			GraphWindow graphCum = new GraphWindow(funcsCumulative, region.getName(),plotChars); 
			graphCum.setX_AxisLabel("Magnitude");
			graphCum.setY_AxisLabel("Cumulative Rate (per yr)");
			graphCum.setYLog(true);
			graphCum.setX_AxisRange(5, 8.5);
			graphCum.setY_AxisRange(1e-6,20);
			graphCum.setPlotLabelFontSize(18);
			graphCum.setAxisLabelFontSize(18);
			graphCum.setTickLabelFontSize(16);
			
			if(saveFiles) {
				try {
					graphCum.saveAsPDF(region.getName().replace(" ", "_")+"_CumMFD_Plot.pdf");
					graphCum.saveAsTXT(region.getName().replace(" ", "_")+"_CumMFD_Plot.txt");
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}

    	

	}
	
	
	/**
	 * Need this because Peter's surface does not work in our calculator
	 * @param erf
	 * @param region
	 * @param minMag
	 * @param numMag
	 * @param deltaMag
	 * @param preserveRates
	 * @return
	 */
	public static SummedMagFreqDist OLDgetMagFreqDistInRegion(NshmErf erf, Region region,
			double minMag,int numMag,double deltaMag, boolean preserveRates) {
		
		String tempString="";

		SummedMagFreqDist magFreqDist = new SummedMagFreqDist(minMag, numMag, deltaMag);
		double duration = erf.getTimeSpan().getDuration();
		ArrayList<String> stringList = new ArrayList<String>();

		for (int s = 0; s < erf.getNumSources(); ++s) {
			ProbEqkSource source = erf.getSource(s);
			if(!stringList.contains(source.getClass().toString()))
				stringList.add(source.getClass().toString());

			for (int r = 0; r < source.getNumRuptures(); ++r) {
				ProbEqkRupture rupture = source.getRupture(r);
				RuptureSurface surf = rupture.getRuptureSurface();
				double mag = rupture.getMag();
				double equivRate = rupture.getMeanAnnualRate(duration);

				if(surf instanceof NshmSurface) {
					Location loc;
					try {
						loc = ((NshmSurface)rupture.getRuptureSurface()).centroid();
						if(region.contains(loc)) {
							magFreqDist.addResampledMagRate(mag, equivRate, preserveRates);
							tempString+=loc.getLongitude()+"\t"+loc.getLatitude()+"\n";
//							System.out.println(loc.getLongitude()+"\t"+loc.getLatitude());
						}
					} catch (Exception e) {
//						e.printStackTrace();
//						System.out.println(rupture.getRuptureSurface().getClass().toString());
					}
				}
				else {	// it's a CompoundSurface
					CompoundSurface compSurf = (CompoundSurface) surf;
					List<? extends RuptureSurface> surfList = compSurf.getSurfaceList();
					double ptRate = equivRate/surfList.size();
					for(RuptureSurface surface: surfList) {
						if(region.contains(((NshmSurface)surface).centroid()))
							magFreqDist.addResampledMagRate(mag, ptRate, preserveRates);
					}
				}
			}
		}
		System.out.println(stringList);
		
//		File file = new File("/Users/field/junkTest.txt");
//		DataOutputStream out;
//		try {
//			out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file)));
//			out.writeUTF(tempString);
////			for (double val : array) {
////				out.writeDouble(val);
////			}
//			out.close();
//		} catch (Exception e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}

		return magFreqDist;
	}

	
	
	
	public static double[] getSlipRates(NSHM23_DeformationModels dm) {
		double[] valArray=null;
		try {
			List<? extends FaultSection> subSects =dm.build(NSHM23_FaultModels.WUS_FM_v1p4);
			valArray = new double[subSects.size()];
			for(int i=0;i<valArray.length;i++)
				valArray[i] = subSects.get(i).getOrigAveSlipRate();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return valArray;
	}
	
	
	public static HashMap<String,Double> getParSectSlipRates(NSHM23_DeformationModels dm) {
		double[] subSectSR=null;
		HashMap<String,Double> mapAveParSlipRates = new HashMap<String,Double>();
		try {
			List<? extends FaultSection> subSects =dm.build(NSHM23_FaultModels.WUS_FM_v1p4);
			subSectSR = new double[subSects.size()];
			for(int i=0;i<subSectSR.length;i++)
				subSectSR[i] = subSects.get(i).getOrigAveSlipRate();

			// get parent section average slip rates
			HashMap<String,ArrayList<Double>> mapSlipRatesForPar = new HashMap<String,ArrayList<Double>>();
			for(int i=0;i<subSectSR.length;i++) {
				String parSectName=subSects.get(i).getParentSectionName();
				if(!mapSlipRatesForPar.keySet().contains(parSectName)) {
					mapSlipRatesForPar.put(parSectName, new ArrayList<Double>());
				}
				mapSlipRatesForPar.get(parSectName).add(subSectSR[i]);
			}
			
			Mean meanCalc = new Mean();
			for(String parFltName : mapSlipRatesForPar.keySet()) {
				double aveForParFlt = meanCalc.evaluate(mapSlipRatesForPar.get(parFltName).stream().mapToDouble(d -> d).toArray());
				mapAveParSlipRates.put(parFltName, aveForParFlt);
			}

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return mapAveParSlipRates;
	}

	


	
	
	

	
	
	public static void defModelAnalysis(double wtArray[], String filenameSuffix) {
		
		double plotMin=1e-3;
		double plotMax=1e2;
		Range plotRange = new Range(plotMin,plotMax);
		DefaultXY_DataSet equalFunc = new DefaultXY_DataSet();
		equalFunc.set(plotMin,plotMin);
		equalFunc.set(plotMax,plotMax);

		
//		List<? extends FaultSection> subSects=null;
//		try {
//			subSects = NSHM23_DeformationModels.GEOLOGIC.build(NSHM23_FaultModels.NSHM23_v1p4);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//		System.out.println("Num subsections = "+subSects.size());

		
//		for(NSHM23_DeformationModels dm:NSHM23_DeformationModels.values()) {
//			System.out.println(dm.getShortName()+"\t"+dm.getNodeWeight(null));
//		}
//		System.exit(0);
		
		
		NSHM23_DeformationModels[] dmArray = {
				NSHM23_DeformationModels.EVANS,
				NSHM23_DeformationModels.POLLITZ,
				NSHM23_DeformationModels.GEOLOGIC,
				NSHM23_DeformationModels.ZENG,
				NSHM23_DeformationModels.SHEN_BIRD };
		
		
		PlotSymbol plotSymbol = PlotSymbol.FILLED_CIRCLE;
    	ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(plotSymbol, 1f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(plotSymbol, 1f, Color.GREEN));
		plotChars.add(new PlotCurveCharacterstics(plotSymbol, 1f, Color.CYAN));
		plotChars.add(new PlotCurveCharacterstics(plotSymbol, 1f, Color.GRAY));
		plotChars.add(new PlotCurveCharacterstics(plotSymbol, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(plotSymbol, 1f, Color.BLACK));
 

		
		String[] shortNameArray = new String[dmArray.length];
		ArrayList<double[]> doubleArrayList = new ArrayList<double[]>();
		
		// make list of parent names
		HashMap<String,Double> tempMap = getParSectSlipRates(NSHM23_DeformationModels.GEOLOGIC);
		ArrayList<String> parSectNameList = new ArrayList<String>();
		for(String parName:tempMap.keySet())
			parSectNameList.add(parName);
		
		int numParSect = parSectNameList.size();
		
		int i=0;
		for(NSHM23_DeformationModels dm:dmArray) {
			shortNameArray[i] = dm.getShortName();
			HashMap<String,Double> slipRateMap = getParSectSlipRates(dm);
			double[] slipRateArray = new double[slipRateMap.size()];
			for(int j=0;j<slipRateMap.size();j++)
				slipRateArray[j] = slipRateMap.get(parSectNameList.get(j));
			doubleArrayList.add(slipRateArray);
			i++;
		}
		
		
		double[] medianSlipRates = new double[numParSect];
		double[] meanSlipRates = new double[numParSect];
		Median median = new Median();
		for(i=0;i<numParSect;i++) {
			double[] valArray = new double[doubleArrayList.size()];
			for(int j=0;j<doubleArrayList.size();j++)
				valArray[j]=doubleArrayList.get(j)[i];
			medianSlipRates[i] = median.evaluate(valArray);
//			if(i==0 ) {
//				for(double val : valArray)
//					System.out.println(val);
//				System.out.println("median: "+medianSlipRates[i]);
//			}
			
			for(int j=0;j<valArray.length;j++) meanSlipRates[i] += wtArray[j]*valArray[j];
		}
		

		ArrayList<XY_DataSet> dmFuncs = new ArrayList<XY_DataSet>();
		for(i=0;i<dmArray.length;i++) {
			DefaultXY_DataSet func = new DefaultXY_DataSet(medianSlipRates,doubleArrayList.get(i));
			func.setName(shortNameArray[i]);
	   		int numBelowPlotMin=0;
    		for(int j=0;j<func.size();j++) {
    			if(func.getY(j)<plotMin) 
    				numBelowPlotMin+=1;;
    		}
    		func.setInfo("numBelowPlotMin = "+numBelowPlotMin+"\tfraction: "+(double)numBelowPlotMin/(double)func.size());
			dmFuncs.add(func);
		}
    	
		DefaultXY_DataSet meanFunc = new DefaultXY_DataSet(medianSlipRates,meanSlipRates); 
		meanFunc.setName("meanFunc");
//    	funcs.add(meanFunc);
		
		// plot each separately
		for (i=0;i<dmFuncs.size();i++) {
			GraphWindow graph = new GraphWindow(dmFuncs.get(i), dmFuncs.get(i).getName(), plotChars.get(i)); 
			graph.setX_AxisLabel("Median Slip Rate");
			graph.setY_AxisLabel("Model Slip Rate");
			graph.setYLog(true);
			graph.setXLog(true);
			graph.setX_AxisRange(plotMin, plotMax);
			graph.setY_AxisRange(plotMin, plotMax);
			graph.setPlotLabelFontSize(18);
			graph.setAxisLabelFontSize(18);
			graph.setTickLabelFontSize(16);
		}
    	
//		GraphWindow graph = new GraphWindow(dmFuncs, "", plotChars); 
//		graph.setX_AxisLabel("Median Slip Rate");
//		graph.setY_AxisLabel("Model Slip Rate");
//		graph.setYLog(true);
//		graph.setXLog(true);
//		graph.setX_AxisRange(plotMin, plotMax);
//		graph.setY_AxisRange(plotMin, plotMax);
//		graph.setPlotLabelFontSize(18);
//		graph.setAxisLabelFontSize(18);
//		graph.setTickLabelFontSize(16);
//		if(filenameSuffix != null) {
//			try {
//				graph.saveAsPDF("SlipRateScatter_"+filenameSuffix+".pdf");
//			} catch (IOException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//		}
		
		if(filenameSuffix != null)
			PlottingUtils.writeAndOrPlotFuncs(
				dmFuncs, 
				plotChars, 
				"",
				"Median Slip Rate (cm/yr)",
				"Model Slip Rate (cm/yr)",
				plotRange,
				plotRange,
				true,
				true,
				3.5,
				3d,
				"SlipRateScatterPlot_"+filenameSuffix, 
				true);
				
//		ArrayList<XY_DataSet> funcs, 
//		ArrayList<PlotCurveCharacterstics> plotChars, 
//		String plotName,
//		String xAxisLabel,
//		String yAxisLabel,
//		Range xAxisRange,
//		Range yAxisRange,
//		boolean logX,
//		boolean logY,
//		double widthInches,
//		double heightInches,
//		String fileNamePrefix, 
//		boolean popupWindow) {

		
		
		ArrayList<XY_DataSet> meanFuncs = new ArrayList<XY_DataSet>();
		meanFuncs.add(meanFunc);
		meanFuncs.add(equalFunc);
    	ArrayList<PlotCurveCharacterstics> plotCharsMean = new ArrayList<PlotCurveCharacterstics>();
    	plotCharsMean.add(new PlotCurveCharacterstics(plotSymbol, 1f, Color.BLACK));
    	plotCharsMean.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.BLACK));
    	if(filenameSuffix != null)
    		PlottingUtils.writeAndOrPlotFuncs(
				meanFuncs, 
				plotCharsMean, 
				"",
				"Median Slip Rate (cm/yr)",
				"Mean Slip Rate (cm/yr)",
				plotRange,
				plotRange,
				true,
				true,
				3.5,
				3d,
				"MeanSlipRateScatter_"+filenameSuffix, 
				true);

    	
//    	GraphWindow graphMean = new GraphWindow(meanFuncs, "", plotCharsMean); 
//		graphMean.setX_AxisLabel("Median Slip Rate");
//		graphMean.setY_AxisLabel("Mean Slip Rate");
//		graphMean.setYLog(true);
//		graphMean.setXLog(true);
//		graphMean.setX_AxisRange(plotMin, plotMax);
//		graphMean.setY_AxisRange(plotMin, plotMax);
//		graphMean.setPlotLabelFontSize(18);
//		graphMean.setAxisLabelFontSize(18);
//		graphMean.setTickLabelFontSize(16);
//		if(filenameSuffix != null) {
//			try {
//				graphMean.saveAsPDF("MeanSlipRateScatter_"+filenameSuffix+".pdf");
//			} catch (IOException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//		}


		
		
		double thresh = 2d;
		ArrayList<DefaultXY_DataSet> outlierFuncs = new ArrayList<DefaultXY_DataSet>();
		for(i=0;i<dmArray.length;i++) {
			DefaultXY_DataSet outlierFunc = new DefaultXY_DataSet();
			outlierFunc.setName(shortNameArray[i]+" Outlier Func");
			outlierFuncs.add(outlierFunc);
		}
		DefaultXY_DataSet meanVsMedianAboveThresh = new DefaultXY_DataSet(); 
		int numBelow=0, numAbove=0;
		System.out.println("index\tmedian\tmean\tmeanVsmedian\toutVal\toutValVsMedian\tmodel\tname");
		for(i=0;i<medianSlipRates.length;i++) {
			double ratio = meanSlipRates[i]/medianSlipRates[i];
			if(ratio>thresh) {
				numAbove+=1;
				meanVsMedianAboveThresh.set(medianSlipRates[i],meanSlipRates[i]);
				double max= Double.NEGATIVE_INFINITY;
				int index = -1;
				for(int j=0;j<doubleArrayList.size();j++) {
					if(doubleArrayList.get(j)[i]*wtArray[j]>max) {
						max= doubleArrayList.get(j)[i]*wtArray[j];
						index=j;
					}
				}
				outlierFuncs.get(index).set(medianSlipRates[i],doubleArrayList.get(index)[i]);
				double outlierRatio = doubleArrayList.get(index)[i]/medianSlipRates[i];
				double meanVsMedian = meanSlipRates[i]/medianSlipRates[i];
				System.out.println(i+"\t"+(float)medianSlipRates[i]+"\t"+(float)meanSlipRates[i]+"\t"+
						+(float)meanVsMedian+"\t"+
						(float)doubleArrayList.get(index)[i]+"\t"+(float)outlierRatio+"\t"+
						shortNameArray[index]+"\t"+parSectNameList.get(i));
			}
			if(ratio<1.0/thresh) numBelow+=1;
		}
		
		// fraction above thresh:
		double fracAbove = (double)meanVsMedianAboveThresh.size()/(double)medianSlipRates.length;
		System.out.println("\nFraction with mean/median > "+thresh+ " is "+(float)fracAbove+"\n");

		// check for empty funcs to avoid plotting problems
		double sum=0;
		for(i=0;i<outlierFuncs.size();i++) {
			DefaultXY_DataSet func = outlierFuncs.get(i);
			double frac = (double)func.size()/(double)meanVsMedianAboveThresh.size();
			System.out.println("Frac of mean/median>"+thresh+ " cases from "+shortNameArray[i]+" = "+(float)frac);
			sum += frac;
			if(func.size()==0)
				func.set(plotMin/10d,plotMin/10d); // add on off plot value
		}
		System.out.println("\t sum="+(float)sum+"\n");

		
		System.out.println("numBelowThresh="+numBelow+"\tfrac="+(double)numBelow/(double)medianSlipRates.length);
		System.out.println("numAboveThresh="+numAbove+"\tfrac="+(double)numAbove/(double)medianSlipRates.length);
		
		
		outlierFuncs.add(meanVsMedianAboveThresh);
		outlierFuncs.add(equalFunc);
		GraphWindow graphMeanThresh = new GraphWindow(outlierFuncs, "Means Outside Threshold", plotChars); 
		graphMeanThresh.setX_AxisLabel("Median Slip Rate");
		graphMeanThresh.setY_AxisLabel("Mean & Outlier Slip Rate");
		graphMeanThresh.setYLog(true);
		graphMeanThresh.setXLog(true);
		graphMeanThresh.setX_AxisRange(plotMin, plotMax);
		graphMeanThresh.setY_AxisRange(plotMin, plotMax);
		graphMeanThresh.setPlotLabelFontSize(18);
		graphMeanThresh.setAxisLabelFontSize(18);
		graphMeanThresh.setTickLabelFontSize(16);
		
		
		HistogramFunction valVsMeanHist = new HistogramFunction(-3, 61, 0.1);
		for(i=0;i<medianSlipRates.length;i++) {
			for(double[] slipRates: doubleArrayList) {
				double logRatio = Math.log10(slipRates[i]/meanSlipRates[i]);
				if(logRatio<valVsMeanHist.getMinX()) logRatio=valVsMeanHist.getMinX();
				else if(logRatio>valVsMeanHist.getMaxX())logRatio=valVsMeanHist.getMaxX();
				valVsMeanHist.add(logRatio, 1.0);
			}
		}
		valVsMeanHist.normalizeBySumOfY_Vals();
		double lowVal = valVsMeanHist.getY(0);
		valVsMeanHist.set(0, 0.0);
//		System.out.println(valVsMeanHist);
		String stdDev = Double.toString(valVsMeanHist.computeStdDev());
		valVsMeanHist.set(0, lowVal);
		valVsMeanHist.setInfo(stdDev);
		GraphWindow valVsMeanHistWindow = new GraphWindow(valVsMeanHist, "Val vs Mean Hist", new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 2f, null, 1f, Color.BLUE)); 
		valVsMeanHistWindow.setX_AxisLabel("Val/Mean");
		valVsMeanHistWindow.setY_AxisLabel("Density");
		valVsMeanHistWindow.setPlotLabelFontSize(18);
		valVsMeanHistWindow.setAxisLabelFontSize(18);
		valVsMeanHistWindow.setTickLabelFontSize(16);

		
		HistogramFunction valVsMedianHist = new HistogramFunction(-3, 61, 0.1);
		for(i=0;i<medianSlipRates.length;i++) {
			for(double[] slipRates: doubleArrayList) {
				double logRatio = Math.log10(slipRates[i]/medianSlipRates[i]);
				if(logRatio<valVsMedianHist.getMinX()) logRatio=valVsMedianHist.getMinX();
				else if(logRatio>valVsMedianHist.getMaxX())logRatio=valVsMedianHist.getMaxX();
				valVsMedianHist.add(logRatio, 1.0);
			}
		}
		valVsMedianHist.normalizeBySumOfY_Vals();
		lowVal = valVsMedianHist.getY(0);
		valVsMedianHist.set(0, 0.0);
		stdDev = Double.toString(valVsMeanHist.computeStdDev());
		valVsMedianHist.set(0, lowVal);
		valVsMedianHist.setInfo(stdDev);
		GraphWindow valVsMedianHistWindow = new GraphWindow(valVsMedianHist, "Val vs Median Hist", new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 2f, null, 1f, Color.BLUE)); 
		valVsMedianHistWindow.setX_AxisLabel("Val/Median");
		valVsMedianHistWindow.setY_AxisLabel("Density");
//		valVsMedianHistWindow.setYLog(true);
//		valVsMedianHistWindow.setXLog(true);
//		valVsMedianHistWindow.setX_AxisRange(1e-3, 1e2);
//		valVsMedianHistWindow.setY_AxisRange(1e-8, 1e3);
		valVsMedianHistWindow.setPlotLabelFontSize(18);
		valVsMedianHistWindow.setAxisLabelFontSize(18);
		valVsMedianHistWindow.setTickLabelFontSize(16);

		
		
		HistogramFunction meanVsMedianHist = new HistogramFunction(-100, 201, 1.0);
		for(i=0;i<medianSlipRates.length;i++)
			meanVsMedianHist.add(meanSlipRates[i]/medianSlipRates[i], 1.0);
		meanVsMedianHist.normalizeBySumOfY_Vals();
		GraphWindow hist = new GraphWindow(meanVsMedianHist, "", new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 2f, null, 1f, Color.BLUE)); 
		hist.setX_AxisLabel("Mean/Median");
		hist.setY_AxisLabel("Density");
//		hist.setYLog(true);
//		hist.setXLog(true);
//		hist.setX_AxisRange(1e-3, 1e2);
//		hist.setY_AxisRange(1e-8, 1e3);
		hist.setPlotLabelFontSize(18);
		hist.setAxisLabelFontSize(18);
		hist.setTickLabelFontSize(16);
		
		GraphWindow histCum = new GraphWindow(meanVsMedianHist.getCumulativeDistFunction(), "", new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.BLUE)); 
		histCum.setX_AxisLabel("Mean/Median");
		histCum.setY_AxisLabel("Cum Density");
//		histCum.setX_AxisRange(1e-3, 1e2);
//		histCum.setY_AxisRange(1e-8, 1e3);
		histCum.setPlotLabelFontSize(18);
		histCum.setAxisLabelFontSize(18);
		histCum.setTickLabelFontSize(16);


		
		
		
		
		

	}



	public static void main(String[] args) {
		

		
		
//		double wtArray[] = {0.1,0.2,0.2,0.25,0.25};
//		defModelAnalysis(wtArray, "OrigWts");
//
//		double wtArray2[] = {0.02,0.08,0.26,0.32,0.32};
//		defModelAnalysis(wtArray2, "RevisedWts");

//		double wtArray[] = {0.0,0.1,0.26,0.32,0.32};

//		test();
		
//		// Model MFD plots
//		ArrayList<Region> regionList = new ArrayList<Region>();
//		try {
//			regionList.add(SeismicityRegions.CONUS_WEST.load());
//			regionList.add(SeismicityRegions.CONUS_EAST.load());
//			regionList.add(AnalysisRegions.CONUS_EAST.load());
//			regionList.add(AnalysisRegions.CONUS_IMW.load());
//			regionList.add(AnalysisRegions.CONUS_PNW.load());
//			regionList.add(AnalysisRegions.CONUS_U3_RELM.load());
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//		makeModelMFD_Plots(regionList, true);
		
//		makeRegMFD_Plots(true);
		
		makeMagAreaPlot(true);
		makeSlipLengthPlot(11, 1000, true);
		
//		makeMagAreaPlot_StableContinental(true);
//		
//		PoissonDistribution pd;
//		GammaDistribution gd;
		
		
//		// write out region names
//		Region region;
//		try {
//			region = SeismicityRegions.CONUS_WEST.load();
//			System.out.println(region.getName());
//			region = SeismicityRegions.CONUS_EAST.load();
//			System.out.println(region.getName());
//			region = AnalysisRegions.CONUS_EAST.load();
//			System.out.println(region.getName());
//			region = AnalysisRegions.CONUS_IMW.load();
//			System.out.println(region.getName());
//			region = AnalysisRegions.CONUS_PNW.load();
//			System.out.println(region.getName());
//			region = AnalysisRegions.CONUS_U3_RELM.load();
//			System.out.println(region.getName());
//
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//		System.exit(0);

		

	}

}
