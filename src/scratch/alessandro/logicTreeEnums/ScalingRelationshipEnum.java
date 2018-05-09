/**
 * 
 */
package scratch.alessandro.logicTreeEnums;

import java.awt.Color;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.Ellsworth_B_WG02_MagAreaRel;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.HanksBakun2002_MagAreaRel;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.Shaw_2009_ModifiedMagAreaRel;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.eq.MagUtils;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.GraphWindow;

import com.google.common.collect.Lists;

/**
 * This was adopted from the class in the UCERF3 project, although this one does not extend LogicTreeBranchNode<ScalingRelationships>, 
 * but does keep the methods of the latter (except getRelativeWeight() because it was hard coded anyway). 
 * 
 * @author field
 *
 */
public enum ScalingRelationshipEnum {
		
	// Wells and Coppersmith (1994, BSSA v.84, p.974-1002) from surface rupture length (SRL) for All slip types
	WC94_SRL_ALL("Wells & Coppersmith (1994) SRL, All", "WC94_SRL_All") {
		
		public double getAveSlip(double area, double length, double origWidth) {
			double mag = getMag(area, length, origWidth);
			double moment = MagUtils.magToMoment(mag);
			return FaultMomentCalc.getSlip(area, moment);
		}
		
		public double getMag(double area, double length, double origWidth) {
			return  5.08 +1.16*Math.log10(length*1e-3);
		}		
	},

	
	// Wells and Coppersmith (1994, BSSA v.84, p.974-1002) from rupture area (RA) for All slip types
	WC94_RA_ALL("Wells & Coppersmith (1994) RA, All", "WC94_RA_All") {
		
		public double getAveSlip(double area, double length, double origWidth) {
			double mag = getMag(area, length, origWidth);
			double moment = MagUtils.magToMoment(mag);
			return FaultMomentCalc.getSlip(area, moment);
		}
		
		public double getMag(double area, double length, double origWidth) {
			return  4.07 +0.98*Math.log10(area*1e-6);
		}		
	},

	
	// Stirling (2002, BSSA v.92, p.812-830) from surface rupture length (SRL) for All slip types
	STIRLING_02_SRL_ALL("Stirling et al. (2002) SRL, All", "S02_SRL_All") {
		
		public double getAveSlip(double area, double length, double origWidth) {
			double mag = getMag(area, length, origWidth);
			double moment = MagUtils.magToMoment(mag);
			return FaultMomentCalc.getSlip(area, moment);
		}
		
		public double getMag(double area, double length, double origWidth) {
			return  5.88 +0.80*Math.log10(length*1e-3);
		}		
	},

	
	// Stirling (2002, BSSA v.92, p.812-830) from rupture area (RA) for All slip types
	STIRLING_02_RA_ALL("Stirling et al. (2002) RA, All", "S02_RA_All") {
		
		public double getAveSlip(double area, double length, double origWidth) {
			double mag = getMag(area, length, origWidth);
			double moment = MagUtils.magToMoment(mag);
			return FaultMomentCalc.getSlip(area, moment);
		}
		
		public double getMag(double area, double length, double origWidth) {
			return  5.09 +0.73*Math.log10(area*1e-6);
		}		
	},
	
	// Wesnousky (2008, BSSA v.98, p.1609-1632) from surface rupture length (SRL) for All slip types
	WESNOUSKY_08_SRL_ALL("Wesnousky (2008) SRL, All", "W08_SRL_All") {
		
		public double getAveSlip(double area, double length, double origWidth) {
			double mag = getMag(area, length, origWidth);
			double moment = MagUtils.magToMoment(mag);
			return FaultMomentCalc.getSlip(area, moment);
		}
		
		public double getMag(double area, double length, double origWidth) {
			return  5.30 +1.02*Math.log10(length*1e-3);
		}		
	},

	// Thingbaijam et al., (2017, BSSA v.107, p.2225-2246)  from surface rupture length (SRL) for normal (N) faulting
	THINGBAIJAM_17_SRL_N("Thingbaijam et al., (2017) SRL, N", "T17_SRL_N") {
		
		public double getAveSlip(double area, double length, double origWidth) {
			double mag = getMag(area, length, origWidth);
			double moment = MagUtils.magToMoment(mag);
			return FaultMomentCalc.getSlip(area, moment);
		}
		
		public double getMag(double area, double length, double origWidth) {
			return (Math.log10(length*1e-3)+1.722)/0.485;
		}		
	},


	AVE_UCERF2("Average UCERF2", "AveU2") {
		
		public double getAveSlip(double area, double length, double origWidth) {
			double areaKm = area/1e6;
			double mag = (ellB_magArea.getMedianMag(areaKm) + hb_magArea.getMedianMag(areaKm))/2;
			double moment = MagUtils.magToMoment(mag);
			return FaultMomentCalc.getSlip(area, moment);
		}
		
		public double getMag(double area, double length, double origWidth) {
			double areaKm = area/1e6;
			return (ellB_magArea.getMedianMag(areaKm) + hb_magArea.getMedianMag(areaKm))/2;
		}		
		
	},
	
	
	SHAW_2009_MOD("Shaw (2009) Modified", "Shaw09Mod") {
		
		public double getAveSlip(double area, double length, double origWidth) {
			double mag = getMag(area, length, origWidth);
			double moment = MagUtils.magToMoment(mag);
			return FaultMomentCalc.getSlip(area, moment);
		}
		
		public double getMag(double area, double length, double origWidth) {
			return sh09_ModMagArea.getWidthDepMedianMag(area*1e-6, origWidth*1e-3);
		}		
		
	},

	HANKS_BAKUN_08("Hanks & Bakun (2008)", "HB08") {
		
		public double getAveSlip(double area, double length, double origWidth) {
			double mag = hb_magArea.getMedianMag(area/1e6);
			double moment = MagUtils.magToMoment(mag);
			return FaultMomentCalc.getSlip(area, moment);
		}
		
		public double getMag(double area, double length, double origWidth) {
			return hb_magArea.getMedianMag(area*1e-6);
		}		
		
	},
	
	

	ELLSWORTH_B("Ellsworth B", "EllB") {
		
		public double getAveSlip(double area, double length, double origWidth) {
			double mag = ellB_magArea.getMedianMag(area/1e6);
			double moment = MagUtils.magToMoment(mag);
			return FaultMomentCalc.getSlip(area, moment);
		}
		
		public double getMag(double area, double length, double origWidth) {
			return ellB_magArea.getMedianMag(area*1e-6);
		}		
		
	},
	
	
	ELLB_SQRT_LENGTH("EllB M(A) & Shaw12 Sqrt Length D(L)", "EllBsqrtLen") {
		
		public double getAveSlip(double area, double length, double origWidth) {
			double impliedAseis = 1.0 - (area/length)/origWidth;
//System.out.println("impliedAseis="+impliedAseis);
			if(impliedAseis>=0.2) {
				double moment = MagUtils.magToMoment(getMag(area, length, origWidth));
				return FaultMomentCalc.getSlip(area, moment);
			}
			double c6 = 5.69e-5;
			double xi = 1.25;
			double w = 15e3;  // units of m
//			double w = xi*area/length;  // units of m
			return c6*Math.sqrt(length*w);
		}
		
		public double getMag(double area, double length, double origWidth) {
			return ellB_magArea.getMedianMag(area*1e-6);
		}		
		
	},

	SHAW_CONST_STRESS_DROP("Shaw09 M(A) & Shaw12 Const Stress Drop D(L)", "ShConStrDrp") {
		
		public double getAveSlip(double area, double length, double origWidth) {
			double impliedAseis = 1.0 - (area/length)/origWidth;
			if(impliedAseis>=0.2) {
				double moment = MagUtils.magToMoment(getMag(area, length, origWidth));
				return FaultMomentCalc.getSlip(area, moment);
			}
			double stressDrop = 4.54;  // MPa
			double xi = 1.25;
			double w = 15e3; // unit of meters
//			double w = xi*area/length; // unit of meters
			double temp = 1.0/(7.0/(3.0*length) + 1.0/(2.0*w))*1e6;
			return stressDrop*temp/FaultMomentCalc.SHEAR_MODULUS;
		}
		
		public double getMag(double area, double length, double origWidth) {
			return sh09_ModMagArea.getWidthDepMedianMag(area*1e-6, origWidth*1e-3);
		}		
		
	},
	
	MEAN_UCERF3("Mean UCERF3 Scaling Relationship", "MeanU3Scale") {
				
		public double getAveSlip(double area, double length, double origWidth) {
			double aveSlip = ScalingRelationshipEnum.SHAW_2009_MOD.getAveSlip(area, length, origWidth) *0.2;
			aveSlip += ScalingRelationshipEnum.SHAW_CONST_STRESS_DROP.getAveSlip(area, length, origWidth) *0.2;
			aveSlip += ScalingRelationshipEnum.ELLB_SQRT_LENGTH.getAveSlip(area, length, origWidth) *0.2;
			aveSlip += ScalingRelationshipEnum.ELLSWORTH_B.getAveSlip(area, length, origWidth) *0.2;
			aveSlip += ScalingRelationshipEnum.HANKS_BAKUN_08.getAveSlip(area, length, origWidth) *0.2;
			return aveSlip;
		}

		@Override
		public double getMag(double area, double length, double origWidth) {
			double aveMag = ScalingRelationshipEnum.SHAW_2009_MOD.getMag(area, length, origWidth) *0.2;
			aveMag += ScalingRelationshipEnum.SHAW_CONST_STRESS_DROP.getAveSlip(area, length, origWidth) *0.2;
			aveMag += ScalingRelationshipEnum.ELLB_SQRT_LENGTH.getAveSlip(area, length, origWidth) *0.2;
			aveMag += ScalingRelationshipEnum.ELLSWORTH_B.getAveSlip(area, length, origWidth) *0.2;
			aveMag += ScalingRelationshipEnum.HANKS_BAKUN_08.getAveSlip(area, length, origWidth) *0.2;
			return aveMag;
		}
	};
	
	private static Ellsworth_B_WG02_MagAreaRel ellB_magArea = new Ellsworth_B_WG02_MagAreaRel();
	private static HanksBakun2002_MagAreaRel hb_magArea = new HanksBakun2002_MagAreaRel();
	private static Shaw_2009_ModifiedMagAreaRel sh09_ModMagArea = new Shaw_2009_ModifiedMagAreaRel();

	private String name, shortName;
	
	private ScalingRelationshipEnum(String name, String shortName) {
		this.name = name;
		this.shortName = shortName;
	}
	
	
	
	/**
	 * This returns the slip (m) for the given rupture area (m-sq) or rupture length (m)
	 * @param area (m)
	 * @param length (m)
	  * @param origWidth (m) - the original down-dip width (before reducing by aseismicity factor)
	 * @return
	 */
	 public abstract double getAveSlip(double area, double length, double origWidth);
	 
	 /**
	  * This returns the magnitude for the given rupture area (m-sq) and width (m)
	 * @param area (m)
	 * @param length (m)
	  * @param origWidth (m) - the original down-dip width (before reducing by aseismicity factor)
	  * @return
	  */
	 public abstract double getMag(double area, double length, double origWidth);
	 

	public String getName() {
		return name;
	}
	
	public String getShortName() {
		return shortName;
	}
	
	public String getBranchLevelName() {
		return "Scaling Relationship";
	}
	
	public static void makeSlipLengthPlot(double downDipWidth, int maxLength, boolean saveFiles) {
		
		ArbitrarilyDiscretizedFunc u2_func = new ArbitrarilyDiscretizedFunc();
		u2_func.setName("AVE_UCERF2");
		ArbitrarilyDiscretizedFunc sh09_funcMod = new ArbitrarilyDiscretizedFunc();
		sh09_funcMod.setName("SHAW_2009_MOD");
		ArbitrarilyDiscretizedFunc ellB_func = new ArbitrarilyDiscretizedFunc();
		ellB_func.setName("ELLSWORTH_B");
		ArbitrarilyDiscretizedFunc hb_func = new ArbitrarilyDiscretizedFunc();
		hb_func.setName("HANKS_BAKUN_08");
		ArbitrarilyDiscretizedFunc sh12_sqrtL_func = new ArbitrarilyDiscretizedFunc();
		sh12_sqrtL_func.setName("ELLB_SQRT_LENGTH");
		ArbitrarilyDiscretizedFunc sh12_csd_func = new ArbitrarilyDiscretizedFunc();
		sh12_csd_func.setName("SHAW_CONST_STRESS_DROP");
		
		
		ScalingRelationshipEnum u2 = ScalingRelationshipEnum.AVE_UCERF2;
		ScalingRelationshipEnum sh09_Mod = ScalingRelationshipEnum.SHAW_2009_MOD;
		ScalingRelationshipEnum ellB = ScalingRelationshipEnum.ELLSWORTH_B;
		ScalingRelationshipEnum hb = ScalingRelationshipEnum.HANKS_BAKUN_08;
		ScalingRelationshipEnum sh12_sqrtL = ScalingRelationshipEnum.ELLB_SQRT_LENGTH;
		ScalingRelationshipEnum sh12_csd = ScalingRelationshipEnum.SHAW_CONST_STRESS_DROP;
		
		
		// log10 area from 1 to 5
    	for(int i=1; i<=maxLength; i++) {
    		double lengthKm = (double)i;
    		double length = lengthKm*1e3;
    		double area = length*downDipWidth*1e3;
    		u2_func.set(lengthKm,u2.getAveSlip(area, length, downDipWidth*1e3));
    		sh09_funcMod.set(lengthKm,sh09_Mod.getAveSlip(area, length, downDipWidth*1e3));
    		ellB_func.set(lengthKm,ellB.getAveSlip(area, length, downDipWidth*1e3));
    		hb_func.set(lengthKm,hb.getAveSlip(area, length, downDipWidth*1e3));
    		sh12_sqrtL_func.set(lengthKm,sh12_sqrtL.getAveSlip(area, length, downDipWidth*1e3));
    		sh12_csd_func.set(lengthKm,sh12_csd.getAveSlip(area, length, downDipWidth*1e3));
    	}
    	
    	ArrayList<ArbitrarilyDiscretizedFunc> funcs = new ArrayList<ArbitrarilyDiscretizedFunc>();
    	funcs.add(sh09_funcMod);
    	funcs.add(ellB_func);
    	funcs.add(hb_func);
    	funcs.add(sh12_sqrtL_func);
    	funcs.add(sh12_csd_func);
    	
    	ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.GREEN));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.MAGENTA));

    	
		GraphWindow graph = new GraphWindow(funcs, "Slip-Length Relationships; DDW="+downDipWidth+" km", plotChars); 
		graph.setX_AxisLabel("Length (km)");
		graph.setY_AxisLabel("Slip (m)");
		graph.setPlotLabelFontSize(18);
		graph.setAxisLabelFontSize(18);
		graph.setTickLabelFontSize(16);
		
		if(saveFiles) {
			try {
				graph.saveAsPDF("slipLengthScalingPlot.pdf");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

	}
	
	
	public static void makeSlipMagPlot(double downDipWidth, int maxLength, boolean saveFiles) {
		
		ArbitrarilyDiscretizedFunc heckerEq9b = new ArbitrarilyDiscretizedFunc();
		heckerEq9b.setName("Hecker Eq 9b");
		ArbitrarilyDiscretizedFunc sh09_funcMod = new ArbitrarilyDiscretizedFunc();
		sh09_funcMod.setName("SHAW_2009_MOD");
		ArbitrarilyDiscretizedFunc ellB_func = new ArbitrarilyDiscretizedFunc();
		ellB_func.setName("ELLSWORTH_B");
		ArbitrarilyDiscretizedFunc hb_func = new ArbitrarilyDiscretizedFunc();
		hb_func.setName("HANKS_BAKUN_08");
		ArbitrarilyDiscretizedFunc sh12_sqrtL_func = new ArbitrarilyDiscretizedFunc();
		sh12_sqrtL_func.setName("ELLB_SQRT_LENGTH");
		ArbitrarilyDiscretizedFunc sh12_csd_func = new ArbitrarilyDiscretizedFunc();
		sh12_csd_func.setName("SHAW_CONST_STRESS_DROP");
		
		
		ScalingRelationshipEnum sh09_Mod = ScalingRelationshipEnum.SHAW_2009_MOD;
		ScalingRelationshipEnum ellB = ScalingRelationshipEnum.ELLSWORTH_B;
		ScalingRelationshipEnum hb = ScalingRelationshipEnum.HANKS_BAKUN_08;
		ScalingRelationshipEnum sh12_sqrtL = ScalingRelationshipEnum.ELLB_SQRT_LENGTH;
		ScalingRelationshipEnum sh12_csd = ScalingRelationshipEnum.SHAW_CONST_STRESS_DROP;
		
		
		// log10 area from 1 to 5
    	for(int i=(int)downDipWidth; i<=maxLength; i++) {
    		double lengthKm = (double)i;
    		double length = lengthKm*1e3;
    		double area = length*downDipWidth*1e3;
    		sh09_funcMod.set(sh09_Mod.getMag(area, length, downDipWidth*1e3),sh09_Mod.getAveSlip(area, length, downDipWidth*1e3));
    		ellB_func.set(ellB.getMag(area, length, downDipWidth*1e3),ellB.getAveSlip(area, length, downDipWidth*1e3));
    		hb_func.set(hb.getMag(area, length, downDipWidth*1e3),hb.getAveSlip(area, length, downDipWidth*1e3));
    		sh12_sqrtL_func.set(ellB.getMag(area, length, downDipWidth*1e3),sh12_sqrtL.getAveSlip(area, length, downDipWidth*1e3));
    		sh12_csd_func.set(sh09_Mod.getMag(area, length, downDipWidth*1e3),sh12_csd.getAveSlip(area, length, downDipWidth*1e3));
    		
    		double heckerMag = hb.getMag(area, length, downDipWidth*1e3);
//    		WC1994_MagAreaRelationship wc94 = new WC1994_MagAreaRelationship();
//    		double heckerMag = wc94.getMedianMag(area*1e-6);
    		double heckerSlip = Math.pow(10.0,0.41*heckerMag-2.79);
    		heckerEq9b.set(heckerMag,heckerSlip);
    		
//    		System.out.println((float)area*1e-6+"\t"+(float)length*1e-3+"\t"+(float)(downDipWidth)+
//    				"\t"+(float)sh09_Mod.getMag(area, downDipWidth*1e3)+
//    				"\t"+(float)ellB.getMag(area, downDipWidth*1e3)+
//    				"\t"+(float)hb.getMag(area, downDipWidth*1e3));

    	}
    	
    	ArrayList<ArbitrarilyDiscretizedFunc> funcs = new ArrayList<ArbitrarilyDiscretizedFunc>();
    	funcs.add(sh09_funcMod);
    	funcs.add(ellB_func);
    	funcs.add(hb_func);
    	funcs.add(sh12_sqrtL_func);
    	funcs.add(sh12_csd_func);
    	funcs.add(heckerEq9b);
    	
    	ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.GREEN));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.MAGENTA));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.CYAN));

    	
		GraphWindow graph = new GraphWindow(funcs, "Implied Slip vs Mag Relationships; DDW="+downDipWidth+" km", plotChars); 
		graph.setX_AxisLabel("Magnitude");
		graph.setY_AxisLabel("Slip (m)");
		graph.setX_AxisRange(6.0, 8.5);
		graph.setY_AxisRange(0.0, 20.0);
		graph.setPlotLabelFontSize(18);
		graph.setAxisLabelFontSize(18);
		graph.setTickLabelFontSize(16);
		
		if(saveFiles) {
			try {
				graph.saveAsPDF("slipMagScalingPlot.pdf");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

	}

	
	
	/**
	 * This tests the magnitudes and implied slip amounts for creeping-section faults
	 * assuming a length and DDW
	 */
	public static void testCreepingSectionSlips() {
		double lengthKm = 150;
		double origWidthKm = 11;
		double widthKm = 1.2;
		double areaKm = lengthKm*widthKm;
		
		ArrayList<ScalingRelationshipEnum> magAreaList = new ArrayList<ScalingRelationshipEnum>();
		magAreaList.add(ScalingRelationshipEnum.ELLSWORTH_B);
		magAreaList.add(ScalingRelationshipEnum.HANKS_BAKUN_08);
		magAreaList.add(ScalingRelationshipEnum.SHAW_2009_MOD);
		
		ArrayList<ScalingRelationshipEnum> aveSlipForRupModelsList= new ArrayList<ScalingRelationshipEnum>();
		aveSlipForRupModelsList.add(ScalingRelationshipEnum.ELLSWORTH_B);
		aveSlipForRupModelsList.add(ScalingRelationshipEnum.HANKS_BAKUN_08);
		aveSlipForRupModelsList.add(ScalingRelationshipEnum.SHAW_2009_MOD);
		aveSlipForRupModelsList.add(ScalingRelationshipEnum.ELLB_SQRT_LENGTH);
		aveSlipForRupModelsList.add(ScalingRelationshipEnum.SHAW_CONST_STRESS_DROP);
		
		
		String result = "CREEPING SECTION Mag and AveSlip (assuming length=150, origDDW=11, and DDW=1.2 km):\n";
		
		for(ScalingRelationshipEnum scale : ScalingRelationshipEnum.values()) {
			double lengthKM = areaKm/origWidthKm;
			double mag = scale.getMag(areaKm*1e6, lengthKM*1e3, origWidthKm*1e3);
			double slip = scale.getAveSlip(areaKm*1e6, lengthKm*1e3, origWidthKm*1e3);
			mag = Math.round(mag*100)/100.;
			slip = Math.round(slip*100)/100.;
			result += (float)mag+"\t"+(float)slip+"\tfor\t"+scale.getShortName()+"\n";
		}

		
		System.out.println(result);

	}
	
	
	/**
	 * This assumes no aseismicity
	 * @param saveFiles
	 */
	public static void makeMagAreaPlot(boolean saveFiles) {
		
		double downDipWidth=11;	// orig down-dip width equals reduced
		
		ArbitrarilyDiscretizedFunc sh09mod_func = new ArbitrarilyDiscretizedFunc();
		sh09mod_func.setName("SHAW_2009_Mod; downDipWidth="+downDipWidth);
		ArbitrarilyDiscretizedFunc ellB_func = new ArbitrarilyDiscretizedFunc();
		ellB_func.setName("ELLSWORTH_B");
		ArbitrarilyDiscretizedFunc hb_func = new ArbitrarilyDiscretizedFunc();
		hb_func.setName("HANKS_BAKUN_08");
		
		ScalingRelationshipEnum sh09mod = ScalingRelationshipEnum.SHAW_2009_MOD;
		ScalingRelationshipEnum ellB = ScalingRelationshipEnum.ELLSWORTH_B;
		ScalingRelationshipEnum hb = ScalingRelationshipEnum.HANKS_BAKUN_08;
		
		// log10 area from 1 to 5
		for(int i=50; i<=20000; i+=10) {
			double area = (double)i;
			double length = area/downDipWidth;
			sh09mod_func.set(area,sh09mod.getMag(area*1e6,length*1e3,downDipWidth*1e3));
			ellB_func.set(area,ellB.getMag(area*1e6,length*1e3,downDipWidth*1e3));
			hb_func.set(area,hb.getMag(area*1e6,length*1e3,downDipWidth*1e3));
		}
    	
    	ArrayList<ArbitrarilyDiscretizedFunc> funcs = new ArrayList<ArbitrarilyDiscretizedFunc>();
    	funcs.add(sh09mod_func);
    	funcs.add(ellB_func);
    	funcs.add(hb_func);
    	
    	ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, null, 1f, Color.BLUE));
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
				graph.saveAsPDF("magAreaScalingPlot.pdf");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}


	
	//public 
	public static void main(String[] args) throws IOException {
		
		 makeSlipLengthPlot(11, 1000, true);
		 makeMagAreaPlot(true);

	}


}
