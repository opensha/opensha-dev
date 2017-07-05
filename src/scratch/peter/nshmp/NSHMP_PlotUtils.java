/**
 * 
 */
package scratch.peter.nshmp;

import static org.opensha.commons.mapping.gmt.GMT_MapGenerator.*;
import static org.opensha.nshmp2.tmp.TestGrid.*;
import static org.opensha.nshmp2.util.Period.*;
import static scratch.peter.curves.ProbOfExceed.*;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.data.xyz.GeoDataSetMath;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.GMT_MapGenerator;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.mapping.gmt.gui.GMT_MapGuiBean;
import org.opensha.commons.mapping.gmt.gui.ImageViewerWindow;
import org.opensha.commons.param.impl.CPTParameter;
import org.opensha.commons.util.FileUtils;
import org.opensha.nshmp2.tmp.TestGrid;
import org.opensha.nshmp2.util.Period;
import org.opensha.sra.rtgm.RTGM;
import org.opensha.sra.rtgm.RTGM.Frequency;

import scratch.peter.curves.ProbOfExceed;

import com.google.common.io.Files;

/**
 *  Class of static methods for plotting USGS NSHMP and OpenSHA hazard data.
 *
 *
 * @author Peter Powers
 * @version $Id:$
 */
public class NSHMP_PlotUtils {

	private final static String DL_DIR = "/Users/pmpowers/Documents/OpenSHA/NSHMPdev2/hazfigs/";
//	private final static String DL_DIR = "/Volumes/Sorrento/hazfigs/";
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
//		makeNationalPE("HazardTool", "DataR3", 0.05, GM0P20, PE2IN50);
//		makeNationalPE("DataR1", "DataR3", 0.05, GM0P20, PE2IN50);
//		makeNationalPE("FortranUpdate", "HazardTool", 0.05, GM1P00, PE2IN50);

//		makeRegionalPE(LOS_ANGELES, "trunk", "FortranUpdate", 0.1, GM0P00, PE2IN50);
//		makeRegionalPE(LOS_ANGELES, "trunk", "FortranUpdate", 0.1, GM0P20, PE2IN50);
//		makeRegionalPE(LOS_ANGELES, "trunk", "FortranUpdate", 0.1, GM1P00, PE2IN50);

//		makeRegionalPE(SAN_FRANCISCO, "trunk", "FortranUpdate", 0.05, GM0P00, PE2IN50);
//		makeRegionalPE(SAN_FRANCISCO, "trunk", "FortranUpdate", 0.05, GM0P20, PE2IN50);
//		makeRegionalPE(SAN_FRANCISCO, "trunk", "FortranUpdate", 0.05, GM1P00, PE2IN50);
//
//		makeRegionalPE(SALT_LAKE_CITY, "trunk", "FortranUpdate", 0.05, GM0P00, PE2IN50);
//		makeRegionalPE(SALT_LAKE_CITY, "trunk", "FortranUpdate", 0.05, GM0P20, PE2IN50);
//		makeRegionalPE(SALT_LAKE_CITY, "trunk", "FortranUpdate", 0.05, GM1P00, PE2IN50);
//
//		makeRegionalPE(SEATTLE, "trunk", "FortranUpdate", 0.05, GM0P00, PE2IN50);
//		makeRegionalPE(SEATTLE, "trunk", "FortranUpdate", 0.05, GM0P20, PE2IN50);
//		makeRegionalPE(SEATTLE, "trunk", "FortranUpdate", 0.05, GM1P00, PE2IN50);
//
//		makeRegionalPE(MEMPHIS, "trunk", "FortranUpdate", 0.05, GM0P00, PE2IN50);
//		makeRegionalPE(MEMPHIS, "trunk", "FortranUpdate", 0.05, GM0P20, PE2IN50);
//		makeRegionalPE(MEMPHIS, "trunk", "FortranUpdate", 0.05, GM1P00, PE2IN50);
		
//		makeRegionalPE_Ratio("hpc", "nshmp_extestopt", GRID_TEST, "FortranUpdateEXT", 0.02, GM0P00, PE2IN50);
//		makeRegionalPE_Ratio("hpc", "nshmp_extest", GRID_TEST, "FortranUpdateEXT", 0.02, GM0P00, PE2IN50);
		

//		makeRegionalPE_RatioSHA("hpc", "muc2-1", "muc2-2", CA_RELM, 0.1, GM0P00, PE2IN50);
//		makeRegionalPE_RatioSHA("hpc", "muc2-1", "muc2-3", CA_RELM, 0.1, GM0P00, PE2IN50);
//		makeRegionalPE_RatioSHA("hpc", "muc2-2", "muc2-3", CA_RELM, 0.1, GM0P00, PE2IN50);
		
//		makeRegionalPE_Ratio("hpc", "nshmp_us", NATIONAL_POLY, "FortranUpdate", 0.1, GM0P00, PE2IN50);
//		tmpRatio("hpc", "muc2_bg", "hpc", "nshmp_ca_bg", CA_RELM, CA_RELM, 0.1, GM0P00, PE2IN50);
//		tmpScaleRatio("hpc", "nshmp_ca_bg", CA_RELM, 0.1, GM0P00, PE2IN50);
		
//		makeRegionalPE_RatioSHA("hpc", "nshmp_ca_nobgopt", "~muc2-1", CA_RELM, 0.1, GM0P00, PE2IN50);
//		makeRegionalPE_RatioSHA("hpc", "muc2_pt", "nshmp_ca_pt", CA_RELM, 0.02, GM0P00, PE2IN50);
//		makeRegionalPE_RatioSHA("hpc", "muc2_pt_EXT", "nshmp_ca_pt_EXT", CA_RELM, 0.02, GM0P00, PE2IN50);
//		makeRegionalPE_RatioSHA("hpc", "nshmp_ca_pt_noEXT-WUS", "muc2_pt_noEXT-WUS", CA_RELM, 0.1, GM0P00, PE2IN50);
//		makeRegionalPE_RatioSHA("hpc", "nshmp_ca_fix", "muc2_fix", CA_RELM, 0.1, GM0P00, PE2IN50);
//		makeRegionalPE_RatioSHA("hpc", "nshmp_ca", "nshmp_ca_epi", CA_RELM, 0.1, GM0P00, PE2IN50);
//		makeRegionalPE_RatioSHA("hpc", "nshmp_ca", "muc2", CA_RELM, 0.1, GM0P00, PE2IN50);
//		makeRegionalPE_RatioSHA("hpc", "muc2_fm2p1", "nshmp_ca", CA_RELM, 0.1, GM0P00, PE2IN50);

//		makeRegionalPE_RatioSHA("hpc", "ca_nshmp_ca_noepi", "ca_nshmp_ca", CA_RELM, 0.1, GM0P00, PE2IN50);
//		makeRegionalPE_RatioSHA("hpc", "mean_uc2_noepi", "ca_nshmp_ca_noepi", CA_RELM, 0.1, GM0P00, PE2IN50);
//		makeRegionalPE_RatioSHA("hpc", "muc2-1", "ca_nshmp", CA_RELM, 0.05, GM0P00, PE2IN50);
//		makeRegionalPE_RatioSHA("hpc", "muc2_fm2p1_noepi", "mean_uc2_noepi", CA_RELM, 0.1, GM0P00, PE2IN50);
//		makeRegionalPE_RatioSHA("hpc", "mmuc2_fm2p1_noepi", "muc2_fm2p1_noepi", CA_RELM, 0.1, GM0P00, PE2IN50);
//		makeRegionalPE_RatioSHA("hpc", "uc2_fm2p1_fss_noepi", "mmuc2_fm2p1_noepi", CA_RELM, 0.1, GM0P00, PE2IN50);
//		makeRegionalPE_RatioSHA("hpc", "uc2_fm2p1_fss_noepi", "ca_nshmp_ca_noepi", CA_RELM, 0.1, GM0P00, PE2IN50);

//		makeNationalRTGM("HazardTool", "DesignMap", 0.05, Frequency.SA_0P20);
//		makeNationalRTGM("HazardTool", "DesignMap", 0.05, Frequency.SA_1P00);
//		makeNationalRTGM("FortranLatest", "HazardTool", 0.05, Frequency.SA_0P20);
//		makeNationalRTGM("FortranLatest", "HazardTool", 0.05, Frequency.SA_1P00);
		
//		makeRegionalRTGM(LOS_ANGELES, "trunk", "FortranUpdate", 0.1, Frequency.SA_0P20);
//		makeRegionalRTGM(LOS_ANGELES, "trunk", "FortranUpdate", 0.1, Frequency.SA_1P00);
//		makeRegionalRTGM(SAN_FRANCISCO, "trunk", "FortranUpdate", 0.1, Frequency.SA_0P20);
//		makeRegionalRTGM(SAN_FRANCISCO, "trunk", "FortranUpdate", 0.1, Frequency.SA_1P00);
//		makeRegionalRTGM(SALT_LAKE_CITY, "trunk", "FortranUpdate", 0.1, Frequency.SA_0P20);
//		makeRegionalRTGM(SALT_LAKE_CITY, "trunk", "FortranUpdate", 0.1, Frequency.SA_1P00);
//		makeRegionalRTGM(SEATTLE, "trunk", "FortranUpdate", 0.1, Frequency.SA_0P20);
//		makeRegionalRTGM(SEATTLE, "trunk", "FortranUpdate", 0.1, Frequency.SA_1P00);
//		makeRegionalRTGM(MEMPHIS, "trunk", "FortranUpdate", 0.1, Frequency.SA_0P20);
//		makeRegionalRTGM(MEMPHIS, "trunk", "FortranUpdate", 0.1, Frequency.SA_1P00);
		
//		makeRegionalMap_SHA(MEMPHIS_BIG, "trunk", GM0P00, PE2IN50);
//		makeRegionalMap_SHA(LOS_ANGELES_BIG, "trunk", GM0P00, PE2IN50);
		
		//makeRegionalMap_NSHMP(MEMPHIS_BIG, "FortranUpdate", GM0P00, PE2IN50);
//		makeRegionalMap_NSHMP(CA_RELM, "FortranUpdate", GM0P00, PE2IN50);
//		makeRegionalMap_NSHMP(CA_RELM, "FortranUpdateGridTest", GM0P00, PE2IN50);
//		makeRegionalMap_SHA("test", "stripe_test", STRIPE, GM0P00, PE2IN50);
		
//		makeRegionalMap_SHA("hpc", "muc2_bg", CA_RELM, GM0P00, PE2IN50);
//		makeRegionalMap_SHA("hpc", "nshmp_ca_bg", CA_RELM, GM0P00, PE2IN50);
//		tmpMap("hpc", "ca_nshmp", CA_RELM, GM0P00, PE2IN50);
//		makeRegionalMap_SHA("hpc", "ca_nshmp", CA_RELM, GM0P00, PE2IN50);
//		makeRegionalMap_SHA("hpc", "ca_nshmp_ca", CA_RELM, GM0P00, PE2IN50);
//		makeRegionalMap_SHA("hpc", "nshmp_ca_epi", CA_RELM, GM0P00, PE2IN50);
//		makeRegionalMap_SHA("hpc", "mean_uc2_noepi_cab", CA_RELM, GM0P00, PE2IN50);
//		makeRegionalMap_SHA("hpc", "muc2_fm2p1_noepi", CA_RELM, GM0P00, PE2IN50);
//		makeRegionalMap_SHA("hpc", "mmuc2_fm2p1_noepi", CA_RELM, GM0P00, PE2IN50);
//		makeRegionalMap_SHA("hpc", "uc2_fm2p1_fss_noepi", CA_RELM, GM0P00, PE2IN50);
		
//		makeRegionalMap_SHA("hpc", "uc3_ref_mean", CA_RELM, GM0P00, PE2IN50);
//		for (int i=0; i<20; i++) {
//			makeRegionalMap_SHA("hpc", "uc3_ref_" + i, CA_RELM, GM0P00, PE2IN50);
//			String title = "UC3 RefVar"+ i + " over Ref 2%50 PGA";
//			makeRegionalPE_RatioSHA("hpc", "uc3_ref_" + i, "uc3_ref_mean", CA_RELM, 0.1, GM0P00, PE2IN50, title);
//		}

//		for (int i=0; i<20; i++) {
//			GeoDataSet ref = getRatioData("hpc", "uc3_ref_19", "uc3_ref_mean", CA_RELM, GM0P00, PE2IN50);
//			String title = "UC3 RefVar"+ i + " over Ref 2%50 PGA";
//			ppTmp(ref, "hpc", "uc3_ref_" + i, "uc3_ref_mean", CA_RELM, 0.1, GM0P00, PE2IN50, title);
//		}

//		makeRegionalMap_SHA("hpc", "us_test", NATIONAL_POLY, GM0P00, PE2IN50);
//		makeRegionalMap_SHA("hpc", "us_test", NATIONAL_POLY, GM1P00, PE2IN50);
//		makeRegionalMap_SHA("hpc", "us_test", NATIONAL_POLY, GM0P20, PE2IN50);
		
//		makeRegionalMap_SHA("hpc", "uc3_ref_mean", CA_RELM, GM0P00, PE2IN50);
//		makeRegionalMap_SHA("hpc", "uc3_ref_mean", CA_RELM, GM0P20, PE2IN50);
//		makeRegionalMap_SHA("hpc", "uc3_ref_mean", CA_RELM, GM1P00, PE2IN50);
		
		// =====================
		
//		String title = "OpenSHA/Fortan 2%50 1Hz w/epi";
//		makeRegionalPE_Ratio("hpc", "nshmp_us_epi", NATIONAL_POLY, "FortranUpdate", 0.1, GM1P00, PE2IN50, title);
		
//		File out = new File("tmp/EXTmap.ch.2p50.no-opt.txt");
//		NSHMP_GeoDataUtils.writePE("hpc", "nshmp_extest", GRID_TEST, GM0P00, PE2IN50,out);
		
//		String title = "All Src SHA/Fort 2%50 PGA w/ epi";
//		makeRegionalPE_Ratio("hpc", "nshmp_ca_allsrc_epi", CA_RELM, "FortranUpdate", 0.1, GM0P00, PE2IN50, title);

//		String title = "CA_SHA/All_SHA 2%50 PGA w/ epi";
//		makeRegionalPE_RatioSHA("hpc", "nshmp_ca_epi", "nshmp_ca_allsrc_epi", CA_RELM, 0.1, GM0P00, PE2IN50, title);

//		String title = "All Src SHA/SHA 2%50 PGA w/o epi";
//		makeRegionalPE_RatioSHA("hpc", "nshmp_ca", "nshmp_ca_epi", CA_RELM, 0.1, GM0P00, PE2IN50, title);

//		String title = "CA Src mUC2/SHA-NSHM 2%50 PGA";
//		makeRegionalPE_RatioSHA("hpc", "muc2-originalbg", "nshmp_ca", CA_RELM, 0.1, GM0P00, PE2IN50, title);

//		String title = "mUC2/SHA-NSHM 2%50 PGA w/o opt";
//		makeRegionalPE_RatioSHA("hpc", "muc2", "nshmp_ca_nobgopt", CA_RELM, 0.1, GM0P00, PE2IN50, title);
		
//		String title = "mUC2fm2p1/mUC2 2%50 PGA";
//		makeRegionalPE_RatioSHA("hpc", "muc2_fm2p1", "muc2", CA_RELM, 0.1, GM0P00, PE2IN50, title);

//		String title = "mod-mUC2fm2p1/mUC2fm2p1 2%50 PGA";
//		makeRegionalPE_RatioSHA("hpc", "mmuc2_fm2p1", "muc2_fm2p1", CA_RELM, 0.1, GM0P00, PE2IN50, title);
		
//		String title = "fss_uc2+Munc/mod-mUC2 2%50 PGA";
//		makeRegionalPE_RatioSHA("hpc", "fss_uc2_fm2p1+Munc", "mmuc2_fm2p1", CA_RELM, 0.1, GM0P00, PE2IN50, title);

//		title = "fss_uc2-Munc/mod-mUC2 2%50 PG";
//		makeRegionalPE_RatioSHA("hpc", "fss_uc2_fm2p1-Munc", "mmuc2_fm2p1", CA_RELM, 0.1, GM0P00, PE2IN50, title);

		//		String title = "EXmap.ch SHA no-opt/opt 2%50 PGA";
//		makeRegionalPE_RatioSHA("hpc", "nshmp_extest", "nshmp_extestopt", GRID_TEST, 0.02, GM0P00, PE2IN50, title);

//		String title = "UC3/UC2 2%50 PGA";
//		makeRegionalPE_RatioSHA("hpc", "uc3_ref_mean", "nshmp_ca_nobgopt", CA_RELM, 0.2, GM0P00, PE2IN50, title);

		
		// ================
//		makeRegionalMap_SHA("hpc", "uc3_ref_bgUC2", CA_RELM, GM0P00, PE2IN50);
//		makeRegionalMap_SHA("hpc", "uc3_ref_bgUC3", CA_RELM, GM0P00, PE2IN50);
		makeRegionalMap_SHA("hpc", "nshmp_ca_pt", CA_RELM, GM0P00, PE2IN50);

//		String title = "UC3bg/UC2bg 2%50 PGA";
//		makeRegionalPE_RatioSHA("hpc", "uc3_ref_mean", "nshmp_ca_nobgopt", CA_RELM, 0.2, GM0P00, PE2IN50, title);

	}
	


	/*
	 * Make a regional hazard map
	 */
	private static void makeRegionalMap_SHA(String shaDir, String runGroup, TestGrid tg, Period p, ProbOfExceed pe) {
		String name = " " + runGroup + " " + pe + " " + p;
		GeoDataSet xyz = NSHMP_GeoDataUtils.getPE_SHA(shaDir, runGroup, tg, p, pe);
		double[] minmax = getRange(p);
		makeMapPlot(xyz, tg.bounds(), minmax[0], minmax[1], name, getCPT(p));
	}
	
	private static void tmpMap(String shaDir, String runGroup, TestGrid tg, Period p, ProbOfExceed pe) {
		String name = " " + runGroup + " " + pe + " " + p;
		GeoDataSet xyz = NSHMP_GeoDataUtils.getPE_SHA(shaDir, runGroup, tg, p, pe);
		double[] minmax = getRange(p);
		makeMapPlot(xyz, tg.bounds(), 0.15, 0.4, name, getCPT(p));
	}

	
	/*
	 * Make a regional hazard map
	 */
	private static void makeRegionalMap_NSHMP(TestGrid tg, String nshmDir, Period p, ProbOfExceed pe) {
		String name = tg + " " + nshmDir + " " + pe + " " + p;
		GeoDataSet xyz = NSHMP_GeoDataUtils.getPE_NSHMP(tg, nshmDir, p, pe);
		double[] minmax = getRange(p);
		makeMapPlot(xyz, tg.bounds(), minmax[0], minmax[1], name, getCPT(p));
	}

	/*
	 * Make a national scale prob. exceedance map with the supplied datasets.
	 */
	private static void makeNationalPE_Ratio(String over, String under,
			double maxScale, Period p, ProbOfExceed pe) {
		String name = over + " over " + under + " " + pe + " " + p;
		GriddedRegion gr = NSHMP_UtilsDev.getNSHMP_Region(0.1);
		GeoDataSet xyz = NSHMP_GeoDataUtils.getPE_Ratio(over, under,
			p, pe, gr);
		NSHMP_GeoDataUtils.minusOne(xyz);
		makeRatioPlot(xyz, NATIONAL.bounds(), -maxScale, maxScale, name, name);
	}
	
	/*
	 * Make a regional over national plot
	 */
	private static void makeRegionalPE_Ratio(String shaDir, String dir, TestGrid over, String under,
			double maxScale, Period p, ProbOfExceed pe, String title) {
		String name = dir + " over " + under + " " + pe + " " + p;
		GeoDataSet xyz = NSHMP_GeoDataUtils.getPE_Ratio_SHAoNSHMP(shaDir, dir, over, under,
			p, pe);
//		NSHMP_GeoDataUtils.minusOne(xyz);
		makeRatioPlot(xyz, over.bounds(), 1-maxScale, 1+maxScale, name, title);
	}

	private static GeoDataSet getRatioData(String shaDir, String dir1, String dir2, TestGrid grid,
			Period p, ProbOfExceed pe) {
		GeoDataSet xyz = NSHMP_GeoDataUtils.getPE_Ratio_SHA(shaDir, dir1, dir2, grid,
			p, pe);
		return xyz;
	}
	
	/* tmp dividing uc3 ratio maps M>5 fix */
	private static void ppTmp(GeoDataSet fix, String shaDir, String dir1, String dir2, TestGrid grid,
			double maxScale, Period p, ProbOfExceed pe, String title) {
		String name =  dir1 + " over " + dir2 + " " + pe + " " + p;
		GeoDataSet xyz = NSHMP_GeoDataUtils.getPE_Ratio_SHA(shaDir, dir1, dir2, grid,
			p, pe);
		GeoDataSet ratioFix = GeoDataSetMath.divide(xyz, fix);
//		NSHMP_GeoDataUtils.minusOne(xyz);
		makeRatioPlot(ratioFix, grid.bounds(), 1-maxScale, 1+maxScale, name, title);
	}
	
	
	/*
	 * Make a sha over sha plot
	 */
	private static void makeRegionalPE_RatioSHA(String shaDir, String dir1, String dir2, TestGrid grid,
			double maxScale, Period p, ProbOfExceed pe, String title) {
		String name =  dir1 + " over " + dir2 + " " + pe + " " + p;
		GeoDataSet xyz = NSHMP_GeoDataUtils.getPE_Ratio_SHA(shaDir, dir1, dir2, grid,
			p, pe);
//		NSHMP_GeoDataUtils.minusOne(xyz);
		makeRatioPlot(xyz, grid.bounds(), 1-maxScale, 1+maxScale, name, title);
	}

	private static void tmpRatio(String shaDir1, String dir1, String shaDir2, String dir2, TestGrid grid1, TestGrid grid2,
			double maxScale, Period p, ProbOfExceed pe) {
		String name =  dir1 + " over " + dir2 + " " + pe + " " + p;
		GeoDataSet xyz = NSHMP_GeoDataUtils.tmpPEratio(shaDir1, dir1, shaDir2, dir2, grid1, grid2, p, pe);
		NSHMP_GeoDataUtils.minusOne(xyz);
		makeRatioPlot(xyz, grid1.bounds(), -maxScale, maxScale, name, name);
	}
	
	private static void tmpScaleRatio(String shaDir1, String dir1, TestGrid grid1,
			double maxScale, Period p, ProbOfExceed pe) {
		String name =  "scaled " + dir1 + " over " + dir1 + " " + pe + " " + p;
		GeoDataSet xyz = NSHMP_GeoDataUtils.tmpScalePEratio(shaDir1, dir1, grid1, p, pe);
		NSHMP_GeoDataUtils.minusOne(xyz);
		makeRatioPlot(xyz, grid1.bounds(), -maxScale, maxScale, name, name);
	}

	/*
	 * Make national scale RTGM comparison map
	 */
	private static void makeNationalRTGM(String over, String under,
			double maxScale, Frequency f) {
		String name = "RTGM " + over + " over " + under + " " + f;
		GriddedRegion gr = NSHMP_UtilsDev.getNSHMP_Region(0.1);
		GeoDataSet xyz = NSHMP_GeoDataUtils.getRTGM_Ratio(over, under, f, gr);
		NSHMP_GeoDataUtils.minusOne(xyz);
		makeRatioPlot(xyz, NATIONAL.bounds(), -maxScale, maxScale, name, name);
	}
	
	/*
	 * Make regional scale RTGM comparison map
	 *
	 */
	private static void makeRegionalRTGM(TestGrid over, String shaDir, String under,
			double maxScale, Frequency f) {
		String name = "RTGM " + over + " SHA" + shaDir + " over " + under + " " + f;
		GeoDataSet xyz = NSHMP_GeoDataUtils.getRTGM_Ratio(over, shaDir, under, f);
		NSHMP_GeoDataUtils.minusOne(xyz);
		makeRatioPlot(xyz, over.bounds(), -maxScale, maxScale, name, name);
	}
	

	private static void makeRatioPlot(GeoDataSet xyz, double[] bounds,
			double scaleMin, double scaleMax, String name, String title) {
		GMT_MapGenerator mapGen = create(bounds);
		mapGen.setParameter(COLOR_SCALE_MIN_PARAM_NAME, scaleMin);
		mapGen.setParameter(COLOR_SCALE_MAX_PARAM_NAME, scaleMax);
		CPTParameter cptParam = (CPTParameter) mapGen.getAdjustableParamsList()
				.getParameter(CPT_PARAM_NAME);
			cptParam.setValue(GMT_CPT_Files.GMT_POLAR.getFileName());
		mapGen.setParameter(LOG_PLOT_NAME, false);
		try {
			GMT_Map map = mapGen.getGMTMapSpecification(xyz);
			map.setCustomLabel(title);
			makeMap(map, mapGen, "No metadata", DL_DIR + name + File.separator);
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}
	
	public static GMT_CPT_Files getCPT(Period p) {
		return (p == GM0P20) ? GMT_CPT_Files.NSHMP_5hz : GMT_CPT_Files.NSHMP_1hz;
//		return (p == GM1P00) ? GMT_CPT_Files.NSHMP_1hz : GMT_CPT_Files.NSHMP_5hz;
	}
	
	public static double[] getRange(Period p) {
		return (p == GM0P20) ? new double[] {0.0, 3.0} : new double[] {0.0, 1.0};
//		return (p == GM1P00) ? new double[] {0.0, 1.0} : new double[] {0.0, 3.0};
	}
	
	private static void makeMapPlot(GeoDataSet xyz, double[] bounds,
			double scaleMin, double scaleMax, String name, GMT_CPT_Files cpt) {
		GMT_MapGenerator mapGen = create(bounds);
		mapGen.setParameter(COLOR_SCALE_MIN_PARAM_NAME, scaleMin);
		mapGen.setParameter(COLOR_SCALE_MAX_PARAM_NAME, scaleMax);
		CPTParameter cptParam = (CPTParameter) mapGen.getAdjustableParamsList()
				.getParameter(CPT_PARAM_NAME);
		cptParam.setValue(cpt.getFileName());
		mapGen.setParameter(LOG_PLOT_NAME, false);
		try {
			GMT_Map map = mapGen.getGMTMapSpecification(xyz);
			map.setCustomLabel(name);
			makeMap(map, mapGen, "No metadata", DL_DIR + name + File.separator);
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}

//	public static void makePolarPlot() {
//		GMT_MapGenerator map = getMapGenNational();
//		
//		//override default scale
//		map.setParameter(COLOR_SCALE_MIN_PARAM_NAME, -0.02);
//		map.setParameter(COLOR_SCALE_MAX_PARAM_NAME, 0.02);
////		map.setParameter(COLOR_SCALE_MIN_PARAM_NAME, 0.9);
////		map.setParameter(COLOR_SCALE_MAX_PARAM_NAME, 1.1);
//		
//		CPTParameter cptParam = (CPTParameter) map.getAdjustableParamsList()
//			.getParameter(CPT_PARAM_NAME);
//		cptParam.setValue(GMT_CPT_Files.GMT_POLAR.getFileName());
//
//		map.setParameter(LOG_PLOT_NAME, false);
//		
////		GeoDataSet xyz = SmoothSeismicitySpatialPDF_Fetcher.getUCERF3_PDF();
//		
////		GeoDataSet xyz = CurveComparisons.test();
//		GeoDataSet xyz = CurveComparisons.getFortranOverR3();
//		for (int i=0;i<xyz.size();i++) {
//			xyz.set(i, 1 - xyz.get(i));
//		}
//		
//		try {
//			makeMap(xyz, "2p50", "No metadata", "PP_test2", map);
//		} catch (IOException ioe) {
//			ioe.printStackTrace();
//		}
//	
//	}
	

		
	/**
	 * Make a friggin map.
	 * 
	 * @param xyz
	 * @param gmt_MapGenerator
	 * @param scaleLabel
	 * @param metadata
	 * @param dirName
	 * @param dlDir
	 * @throws IOException
	 */
	public static void makeMap(GMT_Map map,
			GMT_MapGenerator gmt_MapGenerator,
			String metadata, String dlDir) throws IOException {
		try {
//			if(makeMapOnServer) {
				String url = gmt_MapGenerator.makeMapUsingServlet(map, metadata, null);
				metadata += GMT_MapGuiBean.getClickHereHTML(gmt_MapGenerator.getGMTFilesWebAddress());
				File zipFile = new File(dlDir, "allFiles.zip");
				Files.createParentDirs(zipFile);
				// construct zip URL
				String zipURL = url.substring(0, url.lastIndexOf('/')+1)+"allFiles.zip";
				FileUtils.downloadURL(zipURL, zipFile);
				FileUtils.unzipFile(zipFile, new File(dlDir));
				new ImageViewerWindow(url,metadata, true);
//			} else {
//				gmt_MapGenerator.makeMapLocally(xyz, scaleLabel, metadata, dlDir);
//			}
		} catch (GMT_MapException e) {
			e.printStackTrace();
		}
	}


	
	/*
	 ************************************************************************
	 * 
	 * Basic map generator creation
	 * 
	 ************************************************************************
	 */
		
	private final static boolean makeMapOnServer = true;
	private final static double spacing = 0.1;
	private final static String scaling = COLOR_SCALE_MODE_MANUALLY;
	private final static double minScale = -1.0;
	private final static double maxScale = 1.0;
	private final static String topoRes = TOPO_RESOLUTION_NONE;
	private final static String showHwy = SHOW_HIWYS_NONE;
	private final static String showCoast = COAST_DRAW;
	private final static double imgWidth = 6.5;
	private final static boolean smooth = false;
	private final static boolean bgBlack = false;
	
	/**
	 * MapGenerator creator.
	 * @param bounds
	 * @return a map generator
	 */
	public static GMT_MapGenerator create(double[] bounds) {
		GMT_MapGenerator map = new GMT_MapGenerator();
		CPTParameter cptParam = (CPTParameter) map.getAdjustableParamsList()
			.getParameter(CPT_PARAM_NAME);
		cptParam.setValue(GMT_CPT_Files.GMT_OCEAN2.getFileName());
		map.setParameter(GRID_SPACING_PARAM_NAME, spacing);
		map.setParameter(COLOR_SCALE_MODE_NAME, scaling);
		map.setParameter(COLOR_SCALE_MIN_PARAM_NAME, minScale);
		map.setParameter(COLOR_SCALE_MAX_PARAM_NAME, maxScale);
		map.setParameter(TOPO_RESOLUTION_PARAM_NAME, topoRes);
		map.setParameter(SHOW_HIWYS_PARAM_NAME, showHwy);
		map.setParameter(COAST_PARAM_NAME, showCoast);
		map.setParameter(IMAGE_WIDTH_NAME, imgWidth);
		map.setParameter(GMT_SMOOTHING_PARAM_NAME, smooth);
		map.setParameter(BLACK_BACKGROUND_PARAM_NAME, bgBlack);
		initBounds(map, bounds);
		return map;
	}
	
	private static void initBounds(GMT_MapGenerator map, double[] bounds) {
		map.setParameter(MIN_LAT_PARAM_NAME, bounds[0]);
		map.setParameter(MIN_LON_PARAM_NAME, bounds[2]);
		map.setParameter(MAX_LAT_PARAM_NAME, bounds[1]);
		map.setParameter(MAX_LON_PARAM_NAME, bounds[3]);
	}
	
}

//	// DataR3 over HazardTool datas
//private static void makeR1overR3(double maxScale, Period p) {
//	String name = "DataR1 over DataR3 " + PE2IN50 + " " + p;
//	GriddedRegion gr = NSHMP_UtilsDev.getNSHMP_Region(0.1);
//	GeoDataSet xyz = NSHMP_GeoDataUtils.getPE_Ratio("DataR1", "DataR3",
//		p, PE2IN50, gr);
//	for (int i = 0; i < xyz.size(); i++) {
//		xyz.set(i, 1 - xyz.get(i));
//	}
//	makeRatioPlot(xyz, MapExtent.NATIONAL, -maxScale, maxScale, name);
//}
//
//// DataR3 over HazardTool datas
//private static void makeHToverR3(double maxScale, Period p) {
//	String name = "HazardTool over DataR3 " + PE2IN50 + " " + p;
//	GriddedRegion gr = NSHMP_UtilsDev.getNSHMP_Region(0.1);
//	GeoDataSet xyz = NSHMP_GeoDataUtils.getPE_Ratio("HazardTool", "DataR3",
//		p, PE2IN50, gr);
//	for (int i=0;i<xyz.size();i++) {
//		xyz.set(i, 1 - xyz.get(i));
//	}
//	makeRatioPlot(xyz, MapExtent.NATIONAL, -maxScale, maxScale, name);
//}
//
//// FortranLatest over R3 data
//private static void makeFoverHT(double maxScale, Period p) {
//	String name = "Fortran over HazardTool " + PE2IN50 + " " + p;
//	GriddedRegion gr = NSHMP_UtilsDev.getNSHMP_Region(0.1);
//	GeoDataSet xyz = NSHMP_GeoDataUtils.getPE_Ratio("FortranLatest", "HazardTool",
//		p, PE2IN50, gr);
//	for (int i=0;i<xyz.size();i++) {
//		xyz.set(i, 1 - xyz.get(i));
//	}
//	makeRatioPlot(xyz, MapExtent.NATIONAL, -maxScale, maxScale, name);
//}


