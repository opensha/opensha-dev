package scratch.ned.nshm23;

import java.awt.Color;
import java.io.IOException;
import java.util.ArrayList;

import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;

import scratch.UCERF3.enumTreeBranches.ScalingRelationships;

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
    	funcs.add(ellA_sqrtL_func);
    	funcs.add(ellB_sqrtL_func);
    	
    	ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.GREEN));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.MAGENTA));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.CYAN));
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
		ArbitrarilyDiscretizedFunc hb_func = new ArbitrarilyDiscretizedFunc();
		hb_func.setName("HANKS_BAKUN_08");
		
		NSHM23_ScalingRelationships sh09mod = NSHM23_ScalingRelationships.WIDTH_LIMITED;
		NSHM23_ScalingRelationships ellA = NSHM23_ScalingRelationships.LOGA_C4p1;
		NSHM23_ScalingRelationships ellB = NSHM23_ScalingRelationships.LOGA_C4p2;
		ScalingRelationships hb = ScalingRelationships.HANKS_BAKUN_08;
		
		// log10 area from 1 to 5
    	for(int i=50; i<=20000; i+=10) {
    		double area = (double)i;
    		double rake = Double.NaN;
     		sh09mod_func.set(area,sh09mod.getMag(area*1e6,1e3*area/downDipWidth, downDipWidth*1e3, downDipWidth*1e3, rake));
    		ellA_func.set(area,ellA.getMag(area*1e6,1e3*area/downDipWidth, downDipWidth*1e3, downDipWidth*1e3, rake));
    		ellB_func.set(area,ellB.getMag(area*1e6,1e3*area/downDipWidth, downDipWidth*1e3, downDipWidth*1e3, rake));
    		hb_func.set(area,hb.getMag(area*1e6,1e3*area/downDipWidth, downDipWidth*1e3, downDipWidth*1e3, rake));
    	}
    	
    	ArrayList<ArbitrarilyDiscretizedFunc> funcs = new ArrayList<ArbitrarilyDiscretizedFunc>();
    	funcs.add(sh09mod_func);
    	funcs.add(ellA_func);
    	funcs.add(ellB_func);
    	funcs.add(hb_func);
    	
    	ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.RED));
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


	public static void main(String[] args) {
		
		makeMagAreaPlot(true);
		makeSlipLengthPlot(11, 1000, true);

	}

}
