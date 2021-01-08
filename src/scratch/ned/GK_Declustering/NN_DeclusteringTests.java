package scratch.ned.GK_Declustering;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.jfree.data.Range;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.earthquake.observedEarthquake.parsers.UCERF3_CatalogParser;
import org.opensha.sha.earthquake.observedEarthquake.parsers.USGS_NSHMP_CatalogParser;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import scratch.ned.FSS_Inversion2019.PlottingUtils;

public class NN_DeclusteringTests {

	
	
	
	public static void main(String[] args) throws IOException {
		
		long startMillis = System.currentTimeMillis();
		
//		File file = new File("/Users/field/MiscDocs/MuellerGK_FortranCode/wmm.c2");
//		ObsEqkRupList origRupList = USGS_NSHMP_CatalogParser.loadCatalog(file);
		
		File file = new File("/Users/field/MiscDocs/Hauksson2012_Catalog/hs_1981_2011_06_comb_K2_A.cat_so_SCSN_v01.txt");
		ObsEqkRupList origRupList = Hauksson2012_CatalogParser.loadCatalog(file);
		origRupList.sortByOriginTime();
		
//		IncrementalMagFreqDist mfd = makeMFD(origRupList);
//		GutenbergRichterMagFreqDist mfd_grFit = getBestFitGR(mfd);
//		mfd.setName("Original Catalog");
//		mfd.setInfo("Best fit b-value = "+(float)mfd_grFit.get_bValue());
//		mfd.setName("Input Catalog MFD");

		System.out.println("numEvents = "+origRupList.size());
		
		NearestNeighborDeclustering nnDeclustering = new NearestNeighborDeclustering(origRupList, 3.0);

//		ObsEqkRupList myDeclusteredRupList = NearestNeighborDeclustering.getDeclusteredCatalog(origRupList);
//		IncrementalMagFreqDist mfdMyDeclustered = makeMFD(myDeclusteredRupList);
//		GutenbergRichterMagFreqDist mfdMyDeclustered_grFit = getBestFitGR(mfdMyDeclustered);
//		mfdMyDeclustered.setName("OpenSHA Declustered Catalog");
//		mfdMyDeclustered.setInfo("Best fit b-value = "+(float)mfdMyDeclustered_grFit.get_bValue());

		ArrayList<XY_DataSet> funcs = new  ArrayList<XY_DataSet>();
		funcs.add(nnDeclustering.getNNDistHistogram());
		funcs.addAll(nnDeclustering.getFitGaussianFunctions());
//		funcs.add(mfdDeclustered);
//		funcs.add(mfdMyDeclustered);
//		funcs.add(mfd_grFit);
//		funcs.add(mfdMyDeclustered_grFit);
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 2f, Color.BLUE));
//		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 4f, Color.BLACK));
//		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GREEN));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.MAGENTA));
//		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLUE));
//		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.RED));

		Range xRange = null;
		Range yRange = null;
		PlottingUtils.writeAndOrPlotFuncs(funcs, plotChars, null, "NN Dist", "Num", xRange, yRange, 
				false, false, 3.5, 3.0, null, true);
		
		double runtime = (double)(System.currentTimeMillis()-startMillis)/(1000.*60.);
		System.out.println("Runtime was "+runtime+" minutes");
		
//		ArrayList<XY_DataSet> funcs_cum = new  ArrayList<XY_DataSet>();
//		funcs_cum.add(mfd.getCumRateDistWithOffset());
//		funcs_cum.add(mfdDeclustered.getCumRateDistWithOffset());
//		funcs_cum.add(mfdMyDeclustered.getCumRateDistWithOffset());
////		funcs_cum.add(mfd_grFit.getCumRateDistWithOffset());
////		funcs_cum.add(mfdMyDeclustered_grFit.getCumRateDistWithOffset());
//		yRange = new Range(1e-1, 1e6);
//		PlottingUtils.writeAndOrPlotFuncs(funcs_cum, plotChars, null, "Magnitude", "Cumulative Num", xRange, yRange, 
//				false, true, 3.5, 3.0, null, true);
//
//		
//		// plot fraction that are main shocks
//		IncrementalMagFreqDist fractMainShocks = new IncrementalMagFreqDist(1.05, 80, 0.1);
//		for(int i=0;i<fractMainShocks.size();i++) {
//			fractMainShocks.set(i,mfdMyDeclustered.getY(i)/mfd.getY(i));
//		}
//		IncrementalMagFreqDist fractMainShocksFit = new IncrementalMagFreqDist(mfd_grFit.getMinX(), mfd_grFit.size(), mfd_grFit.getDelta());
//		for(int i=0;i<fractMainShocksFit.size();i++) {
//			double ratio = mfdMyDeclustered_grFit.getY(i)/mfd_grFit.getY(i);
//			if(ratio>1) ratio=1;
//			fractMainShocksFit.set(i,ratio);
//		}
//
//		ArrayList<XY_DataSet> funcs2 = new  ArrayList<XY_DataSet>();
//		funcs2.add(fractMainShocks);
//		funcs2.add(fractMainShocksFit);
//		double fract = mfdMyDeclustered.getCumRate(5.05)/mfd.getCumRate(5.05);
//		fractMainShocks.setInfo("fract main shock M≥5 is "+fract);
//		yRange = new Range(0, 1.2);
//		PlottingUtils.writeAndOrPlotFuncs(funcs2, plotChars, null, "Magnitude", "Fract Main Shock", xRange, yRange, 
//				false, false, 3.5, 3.0, null, true);


		
	}

}
