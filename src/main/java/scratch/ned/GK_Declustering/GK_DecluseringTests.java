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
import org.opensha.sha.earthquake.observedEarthquake.Declustering.GardnerKnopoffDeclustering;
import org.opensha.sha.earthquake.observedEarthquake.parsers.UCERF3_CatalogParser;
import org.opensha.sha.earthquake.observedEarthquake.parsers.USGS_NSHMP_CatalogParser;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import scratch.ned.FSS_Inversion2019.PlottingUtils;

public class GK_DecluseringTests {

	
	private static IncrementalMagFreqDist makeMFD(ObsEqkRupList rupList) {
		IncrementalMagFreqDist mfd = new IncrementalMagFreqDist(1.05, 80, 0.1);
		mfd.setTolerance(0.1);
		for (ObsEqkRupture rup : rupList) {
			if(rup.getMag()>1)
				mfd.add(rup.getMag(), 1.0);
		}
		return mfd;
	}
	
	public static GutenbergRichterMagFreqDist getBestFitGR(IncrementalMagFreqDist mfd) {
		double mMin=2.55;
		double mMax = 7.95;
		int num = 55;
		GutenbergRichterMagFreqDist gr_bestFit = null;
		double totRate = mfd.getCumRate(mMin)-mfd.getCumRate(mMax);
		double last_sum = Double.MAX_VALUE;
		double best_b = -10;
		for(double b=0.5; b<1.5; b+=0.01) {
			GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(b, totRate, mMin, mMax, num);
			double sum = 0;
			for(int i=0;i<gr.size();i++) {
				double mag = gr.getX(i);
				if(mfd.getY(mag)>0)
					sum += Math.pow(Math.log10(gr.getY(mag))-Math.log10(mfd.getY(mag)), 2);
			}
			if(sum<last_sum) {
				gr_bestFit = gr;
				best_b = b;
				last_sum=sum;
			}
		}
		System.out.println("best-fit b = "+best_b);
		
		return gr_bestFit;
	}
	

	public static GutenbergRichterMagFreqDist getBestFitGR_alt(IncrementalMagFreqDist mfd) {
		// for fitting b:
		double mMin=2.55;
		double mMax = 5.95;
		int num = 55;
		double b = (Math.log10(mfd.getY(mMin))-Math.log10(mfd.getY(mMax)))/(mMax-mMin);
		GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(b, mfd.getCumRate(mMin), mMin, 7.95, num);
		System.out.println("best-fit b = "+b);
		
		return gr;
	}

	
	
	public static void main(String[] args) throws IOException {
		

		File file = new File("/Users/field/MiscDocs/MuellerGK_FortranCode/wmm.c2");
		ObsEqkRupList origRupList = USGS_NSHMP_CatalogParser.loadCatalog(file);
		IncrementalMagFreqDist mfd = makeMFD(origRupList);
		GutenbergRichterMagFreqDist mfd_grFit = getBestFitGR(mfd);
		mfd.setName("Original Catalog");
		mfd.setInfo("Best fit b-value = "+(float)mfd_grFit.get_bValue());

		System.out.println("numEvents = "+origRupList.size());

		File fileDeclustered = new File("/Users/field/MiscDocs/MuellerGK_FortranCode/wmm.c3");
		ObsEqkRupList declusteredRupList = USGS_NSHMP_CatalogParser.loadCatalog(fileDeclustered);
		IncrementalMagFreqDist mfdDeclustered = makeMFD(declusteredRupList);
		mfdDeclustered.setName("Chuck's Declustered Catalog");
		
		ObsEqkRupList myDeclusteredRupList = GardnerKnopoffDeclustering.getDeclusteredCatalog(origRupList);
		IncrementalMagFreqDist mfdMyDeclustered = makeMFD(myDeclusteredRupList);
		GutenbergRichterMagFreqDist mfdMyDeclustered_grFit = getBestFitGR(mfdMyDeclustered);
		mfdMyDeclustered.setName("OpenSHA Declustered Catalog");
		mfdMyDeclustered.setInfo("Best fit b-value = "+(float)mfdMyDeclustered_grFit.get_bValue());

		
		mfd.setName("Input Catalog MFD");
		ArrayList<XY_DataSet> funcs = new  ArrayList<XY_DataSet>();
		funcs.add(mfd);
		funcs.add(mfdDeclustered);
		funcs.add(mfdMyDeclustered);
//		funcs.add(mfd_grFit);
//		funcs.add(mfdMyDeclustered_grFit);
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 4f, Color.BLACK));
//		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GREEN));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 2f, Color.RED));

		Range xRange = new Range(2, 8);
		Range yRange = new Range(1e-1, 1e5);
		PlottingUtils.writeAndOrPlotFuncs(funcs, plotChars, null, "Magnitude", "Num", xRange, yRange, 
				false, true, 3.5, 3.0, null, true);
		
		
		ArrayList<XY_DataSet> funcs_cum = new  ArrayList<XY_DataSet>();
		funcs_cum.add(mfd.getCumRateDistWithOffset());
		funcs_cum.add(mfdDeclustered.getCumRateDistWithOffset());
		funcs_cum.add(mfdMyDeclustered.getCumRateDistWithOffset());
//		funcs_cum.add(mfd_grFit.getCumRateDistWithOffset());
//		funcs_cum.add(mfdMyDeclustered_grFit.getCumRateDistWithOffset());
		yRange = new Range(1e-1, 1e6);
		PlottingUtils.writeAndOrPlotFuncs(funcs_cum, plotChars, null, "Magnitude", "Cumulative Num", xRange, yRange, 
				false, true, 3.5, 3.0, null, true);

		
		// plot fraction that are main shocks
		IncrementalMagFreqDist fractMainShocks = new IncrementalMagFreqDist(1.05, 80, 0.1);
		for(int i=0;i<fractMainShocks.size();i++) {
			fractMainShocks.set(i,mfdMyDeclustered.getY(i)/mfd.getY(i));
		}
		IncrementalMagFreqDist fractMainShocksFit = new IncrementalMagFreqDist(mfd_grFit.getMinX(), mfd_grFit.size(), mfd_grFit.getDelta());
		for(int i=0;i<fractMainShocksFit.size();i++) {
			double ratio = mfdMyDeclustered_grFit.getY(i)/mfd_grFit.getY(i);
			if(ratio>1) ratio=1;
			fractMainShocksFit.set(i,ratio);
		}

		ArrayList<XY_DataSet> funcs2 = new  ArrayList<XY_DataSet>();
		funcs2.add(fractMainShocks);
		funcs2.add(fractMainShocksFit);
		double fract = mfdMyDeclustered.getCumRate(5.05)/mfd.getCumRate(5.05);
		fractMainShocks.setInfo("fract main shock M≥5 is "+fract);
		yRange = new Range(0, 1.2);
		PlottingUtils.writeAndOrPlotFuncs(funcs2, plotChars, null, "Magnitude", "Fract Main Shock", xRange, yRange, 
				false, false, 3.5, 3.0, null, true);


		
	}

}
