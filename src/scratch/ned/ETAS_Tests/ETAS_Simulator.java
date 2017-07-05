package scratch.ned.ETAS_Tests;

import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import org.opensha.commons.data.function.AbstractXY_DataSet;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.PlotColorAndLineTypeSelectorControlPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.mapping.gmt.GMT_MapGenerator;
import org.opensha.commons.mapping.gmt.gui.GMT_MapGuiBean;
import org.opensha.commons.mapping.gmt.gui.ImageViewerWindow;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.calc.ERF_Calculator;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.UCERF2;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;
import org.opensha.sha.faultSurface.AbstractEvenlyGriddedSurface;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.StirlingGriddedSurface;
import org.opensha.sha.gui.infoTools.CalcProgressBar;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.sha.imr.param.PropagationEffectParams.DistanceRupParameter;
import org.opensha.sha.magdist.ArbIncrementalMagFreqDist;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import scratch.UCERF3.erf.ETAS.ETAS_Utils;
import scratch.UCERF3.erf.ETAS.IntegerPDF_FunctionSampler;
import scratch.UCERF3.utils.ModUCERF2.ModMeanUCERF2;

public class ETAS_Simulator {
	
	String dirToSaveData;
	String infoForOutputFile;
	
	ArrayList<EqksInGeoBlock> blockList;
	GriddedRegion griddedRegion;
	AbstractERF erf;
	ETAS_Utils etasUtils;
//	double distDecay=1.4;
	double distDecay;
	double minDist;
//	double minDist=2.0;
	double tMin=0;		//days
	double tMax=360;	//days
	boolean useAdaptiveBlocks=true;
	boolean includeBlockRates=true;
	int sourceID_ToIgnore = -1;
	
	ProbEqkRupture mainShock;
	
	ArrayList<PrimaryAftershock> primaryAftershockList;
	ETAS_PrimaryEventSamplerJUNK2 etas_FirstGenSampler, etas_sampler;
	EqkRupture parentRup;
	double expectedNumPrimaryAftershocks, expectedNum;
	SummedMagFreqDist totalExpectedNumMagDist;
	ArrayList<PrimaryAftershock> allAftershocks;

	
	
	public ETAS_Simulator(AbstractERF erf, GriddedRegion griddedRegion, int sourceID_ToIgnore) {
		this.erf = erf;
		this.griddedRegion = griddedRegion;
		this.sourceID_ToIgnore = sourceID_ToIgnore;
		etasUtils = new ETAS_Utils();
		makeAllEqksInGeoBlocks();
				
	}
	
	/**
	 * This returns a list of randomly sampled primary aftershocks
	 * @param parentRup
	 * @return list of PrimaryAftershock objects
	 */
	public ArrayList<PrimaryAftershock> getPrimaryAftershocksList(PrimaryAftershock parentRup) {
				
		double originTime = parentRup.getOriginTime();
		
		int parentID = parentRup.getID();
		
				
		// compute the number of primary aftershocks:
		expectedNum = etasUtils.getDefaultExpectedNumEvents(parentRup.getMag(), originTime, tMax);
		int numAftershocks = etasUtils.getPoissonRandomNumber(expectedNum);	//
//		System.out.println("\tMag="+(float)parentRup.getMag()+"\tOriginTime="+(float)originTime+"\tExpNum="+
//				(float)expectedNum+"\tSampledNum = "+numAftershocks);
		
		if(numAftershocks == 0)
			return new ArrayList<PrimaryAftershock>();
		
		// Make the ETAS sampler for the given main shock:
		etas_sampler = new ETAS_PrimaryEventSamplerJUNK2(parentRup,blockList, erf, distDecay,minDist, useAdaptiveBlocks, includeBlockRates);
		
		// Write spatial probability dist data to file:
//		ETAS_sampler.writeRelBlockProbToFile();

		// Now make the list of aftershocks:
		ArrayList<PrimaryAftershock> primaryAftershocks = new ArrayList<PrimaryAftershock>();
		for (int i = 0; i < numAftershocks; i++) {
			PrimaryAftershock aftershock = etas_sampler.samplePrimaryAftershock(etasUtils.getDefaultRandomTimeOfEvent(originTime, tMax));
			aftershock.setParentID(parentID);
			double dist = LocationUtils.distanceToSurfFast(aftershock.getHypocenterLocation(), parentRup.getRuptureSurface());
			aftershock.setDistanceToParent(dist);
			primaryAftershocks.add(aftershock);
		}
		// save info if this is a mainshock
		if(parentRup.getGeneration() == 0) {
			etas_FirstGenSampler = etas_sampler;
			expectedNumPrimaryAftershocks = expectedNum;
			primaryAftershockList = primaryAftershocks;
		}
		
		return primaryAftershocks;
	}
	
	public void plotERF_MagFreqDists() {
		
		// This makes an MFD of the original, total ERF (for entire region)
		ArbIncrementalMagFreqDist erfMagDist = new ArbIncrementalMagFreqDist(2.05, 8.95, 70);
		double duration = erf.getTimeSpan().getDuration();
		for(int s=0;s<erf.getNumSources();s++) {
			ProbEqkSource src = erf.getSource(s);
			for(int r=0;r<src.getNumRuptures();r++) {
				ProbEqkRupture rup = src.getRupture(r);
				erfMagDist.addResampledMagRate(rup.getMag(), rup.getMeanAnnualRate(duration), true);
			}			
		}
		erfMagDist.setName("Total Mag-Prob Dist for "+erf.getName());
		erfMagDist.setInfo(" ");
		EvenlyDiscretizedFunc erfCumDist = erfMagDist.getCumRateDistWithOffset();
		erfCumDist.setName("Total Cum Mag-Freq Dist for "+erf.getName());
		erfMagDist.scaleToCumRate(2.05, 1); // normalize to mag-prob dist

		// Plot cum MFDs for ERF (plus Karen's UCERF2 obs range)
		ArrayList funcs3 = new ArrayList();
		funcs3.addAll(UCERF2.getObsCumMFD(true));
		funcs3.add(erfCumDist);
		GraphWindow sr_graph3 = new GraphWindow(funcs3, "Cum MFD for ERF and Karen's Obs for CA"); 
		sr_graph3.setX_AxisLabel("Mag");
		sr_graph3.setY_AxisLabel("Cumulative Rate");
		sr_graph3.setYLog(true);
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 5f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 5f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 5f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		sr_graph3.setPlotChars(plotChars);
		sr_graph3.setY_AxisRange(1e-6, sr_graph3.getY_AxisRange().getUpperBound());
		
	}

		
	
	/**
	 * This generates two mag-freq dist plots for the primary aftershock sequence.
	 * @param srcIndex - include if you want to show that for specific source (leave null otherwise)
	 */
	public String plotMagFreqDists(Integer srcIndex, String info, boolean savePDF_Files) {
		
		double days = tMax-tMin;
		ArrayList<EvenlyDiscretizedFunc> magProbDists = new ArrayList<EvenlyDiscretizedFunc>();
		ArrayList<EvenlyDiscretizedFunc> expNumDists = new ArrayList<EvenlyDiscretizedFunc>();

		// FIRST PLOT MAG-PROB DISTS
		
		// get the expected mag-prob dist for the primary aftershock sequence:
		ArbIncrementalMagFreqDist expMagProbDist = etas_FirstGenSampler.getMagProbDist();
		expMagProbDist.setName("Expected Mag-Prob Dist for first-generation events");
		expMagProbDist.setInfo(" ");
		// get expected cum MFD for sequence
		EvenlyDiscretizedFunc expCumDist = expMagProbDist.getCumRateDist();
		expCumDist.scale(expectedNumPrimaryAftershocks);
		expCumDist.setName("Expected num 1st-generation >=M in "+(float)days+" days");
		expCumDist.setInfo(" ");
		magProbDists.add(expMagProbDist);
		expNumDists.add(expCumDist);


		// get the observed mag-prob dist for the sampled set of events 
		ArbIncrementalMagFreqDist obsMagProbDist = new ArbIncrementalMagFreqDist(2.05,8.95, 70);
		for (PrimaryAftershock event : primaryAftershockList)
			obsMagProbDist.addResampledMagRate(event.getMag(), 1.0, true);
		EvenlyDiscretizedFunc obsCumDist = obsMagProbDist.getCumRateDist();
		obsMagProbDist.scaleToCumRate(2.05, 1);	// normalize to 1.0
		obsMagProbDist.setName("Sampled Mag-Prob Dist for first-generation events");
		obsMagProbDist.setInfo(" ");
		obsCumDist.setName("Sampled num 1st-generation >=M in "+(float)days+" days");
		obsCumDist.setInfo(" ");
		magProbDists.add(obsMagProbDist);
		expNumDists.add(obsCumDist);
		
		// MFDs for the specified source
		EvenlyDiscretizedFunc expCumDistForSource = new EvenlyDiscretizedFunc(0.0,10.0,10); // bogus function just so one exists
		if(srcIndex != null) {
			ArbIncrementalMagFreqDist expMagProbDistForSource = etas_FirstGenSampler.getMagProbDistForSource(srcIndex);
			expCumDistForSource = expMagProbDistForSource.getCumRateDist();
			expCumDistForSource.scale(expectedNumPrimaryAftershocks);
			expCumDistForSource.setName("Expected num 1st-generation >=M in "+(float)days+" days for source: "+erf.getSource(srcIndex).getName());
			expCumDistForSource.setInfo(" ");
			expMagProbDistForSource.setName("Expected Mag-Prob Dist for first-generation events for source: "+erf.getSource(srcIndex).getName());
			expMagProbDistForSource.setInfo(" ");
			magProbDists.add(expMagProbDistForSource);
			expNumDists.add(expCumDistForSource);
		}
		
		// now make the cumulative distributions for the entire sequence
		EvenlyDiscretizedFunc totalExpCumDist = totalExpectedNumMagDist.getCumRateDist();
		totalExpCumDist.setName("Expected total num aftershocks >=M in "+(float)days+" days (all generations)");
		totalExpCumDist.setInfo(" ");
		totalExpectedNumMagDist.scaleToCumRate(2.05, 1.0);	// a permanent change!!!!
		totalExpectedNumMagDist.setName("Expected Mag-Prob Dist for all events (all generations)");
		totalExpectedNumMagDist.setInfo(" ");
		magProbDists.add(totalExpectedNumMagDist);
		expNumDists.add(totalExpCumDist);


		// get the observed mag-prob dist for all sampled events 
		ArbIncrementalMagFreqDist totObsMagProbDist = new ArbIncrementalMagFreqDist(2.05,8.95, 70);
		for (PrimaryAftershock event : allAftershocks)
			totObsMagProbDist.addResampledMagRate(event.getMag(), 1.0, true);
		EvenlyDiscretizedFunc totObsCumDist = totObsMagProbDist.getCumRateDist();
		totObsMagProbDist.scaleToCumRate(2.05, 1);	// normalize to 1.0
		totObsMagProbDist.setName("Sampled Mag-Prob Dist for all events (all generations)");
		totObsMagProbDist.setInfo(" ");
		totObsCumDist.setName("Sampled num for all event >=M in "+(float)days+" days (all generations)");
		totObsCumDist.setInfo(" ");
		magProbDists.add(totObsMagProbDist);
		expNumDists.add(totObsCumDist);

		
		// GR dist for comparison
		GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(1.0, 1.0,2.55,8.95, 65);
		double scaleBy = expectedNumPrimaryAftershocks/gr.getY(2.55);
		gr.scale(scaleBy);
		String name = gr.getName();
		gr.setName(name+" (for comparison)");
		expNumDists.add(gr);
		
		
		// Plot these MFDs
		GraphWindow magProbDistsGraph = new GraphWindow(magProbDists, "Mag-Prob Distributions for "+info+" Aftershocks"); 
		magProbDistsGraph.setX_AxisLabel("Mag");
		magProbDistsGraph.setY_AxisLabel("Probability");
		magProbDistsGraph.setY_AxisRange(1e-8, magProbDistsGraph.getY_AxisRange().getUpperBound());
		magProbDistsGraph.setYLog(true);
		
	
		GraphWindow expNumDistGraph = new GraphWindow(expNumDists, "Mag-Num Distributions for Aftershocks for "+(float)days+" days following "+info); 
		expNumDistGraph.setX_AxisLabel("Mag");
		expNumDistGraph.setY_AxisLabel("Expected Num");
		expNumDistGraph.setY_AxisRange(1e-6, expNumDistGraph.getY_AxisRange().getUpperBound());
		expNumDistGraph.setYLog(true);
		
		
		// plot the contributing sources for primary events
		ArrayList<ArbIncrementalMagFreqDist> mfdList = etas_FirstGenSampler.getMagProbDistForAllSources(10);
		SummedMagFreqDist sumMFD = new  SummedMagFreqDist(2.05, 8.95, 70);
		for(ArbIncrementalMagFreqDist mfd:mfdList)
			sumMFD.addResampledMagFreqDist(mfd, true);
		sumMFD.setName("Total Incremental MFD");
		sumMFD.setInfo(" ");
		ArrayList<EvenlyDiscretizedFunc> mfdListFinal = new ArrayList<EvenlyDiscretizedFunc>();
		mfdListFinal.add(sumMFD);
		mfdListFinal.add(sumMFD.getCumRateDistWithOffset());
		mfdListFinal.get(1).setName("Total Cum MFD");
		mfdListFinal.addAll(mfdList);
		GraphWindow contributingPrimarySrcsGraph = new GraphWindow(mfdListFinal, "Mag-Prob Dists for Contributing Sources to Primary Events for "+info); 
		contributingPrimarySrcsGraph.setX_AxisLabel("Mag");
		contributingPrimarySrcsGraph.setY_AxisLabel("Probability");
		contributingPrimarySrcsGraph.setY_AxisRange(1e-6, contributingPrimarySrcsGraph.getY_AxisRange().getUpperBound());
		contributingPrimarySrcsGraph.setYLog(true);
		
		// plot above as cumulative distributions
		ArrayList<EvenlyDiscretizedFunc> cumMFD_List = new ArrayList<EvenlyDiscretizedFunc>();
		cumMFD_List.add(sumMFD.getCumRateDistWithOffset());
		for(ArbIncrementalMagFreqDist mfd:mfdList)
			cumMFD_List.add(mfd.getCumRateDistWithOffset());
		GraphWindow contrCumPrimarySrcsGraph = new GraphWindow(cumMFD_List, "Cumulative Mag-Prob Dists for Contributing Sources to Primary Events for "+info); 
		contrCumPrimarySrcsGraph.setX_AxisLabel("Mag");
		contrCumPrimarySrcsGraph.setY_AxisLabel("Probability");
		contrCumPrimarySrcsGraph.setY_AxisRange(1e-6, contrCumPrimarySrcsGraph.getY_AxisRange().getUpperBound());
		contrCumPrimarySrcsGraph.setYLog(true);
		
		if(savePDF_Files) {
			try {
				magProbDistsGraph.saveAsPDF(dirToSaveData+"magProbDistsGraph.pdf");
				expNumDistGraph.saveAsPDF(dirToSaveData+"magNumDistsGraph.pdf");
				contributingPrimarySrcsGraph.saveAsPDF(dirToSaveData+"contrPrimarySrcsMPDsGraph.pdf");
				contrCumPrimarySrcsGraph.saveAsPDF(dirToSaveData+"contrPrimarySrcsCumMPDsGraph.pdf");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		// make the infoText
		String srcName=erf.getSource(srcIndex).getName();
		String infoText = "Expected num primary aftershocks above mags:\n\n\tMag\tTotal\t"+srcName+"\t(% from "+srcName+")\n";
		double cumRateSrc, cumRate, ratio;
		expCumDist.setTolerance(0.0001);
		expCumDistForSource.setTolerance(0.0001);
		for(int i=0; i<3; i++) {
			double mag = i*0.5+6.5;
			cumRate = expCumDist.getY(mag+0.05);
			if(srcIndex != null) {
				cumRateSrc = expCumDistForSource.getY(mag+0.05);
				ratio = 100*cumRateSrc/cumRate;				
			}
			else {
				cumRateSrc = Double.NaN;
				ratio = Double.NaN;				
			}
			infoText += "\t"+mag+"\t"+(float)cumRate+"\t"+(float)cumRateSrc+"\t"+(float)ratio+"\n";
		}
		// add row for the main shock mag
		double mag = mainShock.getMag();
		cumRate = expCumDist.getY(mag);
		if(srcIndex != null) {
			cumRateSrc = expCumDistForSource.getY(mag);
			ratio = 100*cumRateSrc/cumRate;				
		}
		else {
			cumRateSrc = Double.NaN;
			ratio = Double.NaN;				
		}
		infoText += "\t"+mag+"\t"+(float)cumRate+"\t"+(float)cumRateSrc+"\t"+(float)ratio+"\n";


		return infoText;
	}
	
	
	/**
	 * 
	 * @param info
	 */
	public void plotNumVsTime(String info, boolean savePDF_File, ProbEqkRupture mainShock) {
		
		double delta = 1.0; // days

		// make the target function & change it to a PDF
		EvenlyDiscretizedFunc targetFunc = etasUtils.getDefaultNumWithTimeFunc(mainShock.getMag(), tMin, tMax, delta);
		targetFunc.setName("Expected Number for First-generation Aftershocks");
		
		int numPts = (int) Math.round(tMax-tMin);
		EvenlyDiscretizedFunc allEvents = new EvenlyDiscretizedFunc(tMin+0.5,numPts,delta);
		EvenlyDiscretizedFunc firstGenEvents= new EvenlyDiscretizedFunc(tMin+0.5,numPts,delta);
		allEvents.setTolerance(2.0);
		firstGenEvents.setTolerance(2.0);
		for (PrimaryAftershock event : allAftershocks) {
			double time = event.getOriginTime();
			allEvents.add(time, 1.0);
			if(event.getGeneration() == 1)
				firstGenEvents.add(time, 1.0);
		}
		allEvents.setName("All aftershocks");
		firstGenEvents.setName("First-generation aftershocks");
		allEvents.setInfo(" ");
		firstGenEvents.setInfo(" ");
		
		ArrayList funcs = new ArrayList();
		funcs.add(allEvents);
		funcs.add(firstGenEvents);
		funcs.add(targetFunc);
		
		GraphWindow graph = new GraphWindow(funcs, "Num aftershocks per day for "+info); 
		graph.setX_AxisLabel("Days (since main shock)");
		graph.setY_AxisLabel("Num Events");
		graph.setX_AxisRange(0.4, 360);
		graph.setY_AxisRange(0.1, graph.getY_AxisRange().getUpperBound());
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		graph.setPlotChars(plotChars);
		graph.setYLog(true);
		graph.setXLog(true);
		if(savePDF_File)
		try {
			graph.saveAsPDF(dirToSaveData+"numAshocksVsTime.pdf");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	private DefaultXY_DataSet getEpicenterLocsXY_DataSet(double magLow, double magHigh, int generation) {
		DefaultXY_DataSet epicenterLocs = new DefaultXY_DataSet();
		for (PrimaryAftershock event : allAftershocks) {
			if(event.getMag()>=magLow && event.getMag()<magHigh && event.getGeneration()==generation)
				epicenterLocs.set(event.getHypocenterLocation().getLongitude(), event.getHypocenterLocation().getLatitude());
		}
		epicenterLocs.setName("Generation "+generation+" Aftershock Epicenters for "+magLow+"<=Mag<"+magHigh);
		return epicenterLocs;
	}
	
	
	
	public void plotEpicenterMap(String info, boolean savePDF_File, ProbEqkRupture mainShock) {
		
		ArrayList<AbstractXY_DataSet> funcs = new ArrayList<AbstractXY_DataSet>();
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();

		// M<5
		DefaultXY_DataSet epLocsGen1_Mlt5 = getEpicenterLocsXY_DataSet(2.0, 5.0, 1);
		if(epLocsGen1_Mlt5.size()>0) {
			funcs.add(epLocsGen1_Mlt5);
			epLocsGen1_Mlt5.setInfo("(circles)");
			plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE, 1f, Color.BLACK));
		}
		DefaultXY_DataSet epLocsGen2_Mlt5 = getEpicenterLocsXY_DataSet(2.0, 5.0, 2);
		if(epLocsGen2_Mlt5.size()>0) {
			funcs.add(epLocsGen2_Mlt5);
			epLocsGen2_Mlt5.setInfo("(circles)");
			plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE, 1f, Color.BLUE));
		}
		DefaultXY_DataSet epLocsGen3_Mlt5 = getEpicenterLocsXY_DataSet(2.0, 5.0, 3);
		if(epLocsGen3_Mlt5.size()>0) {
			funcs.add(epLocsGen3_Mlt5);
			epLocsGen3_Mlt5.setInfo("(circles)");
			plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE, 1f, Color.GREEN));
		}
		DefaultXY_DataSet epLocsGen4_Mlt5 = getEpicenterLocsXY_DataSet(2.0, 5.0, 4);
		if(epLocsGen4_Mlt5.size()>0) {
			funcs.add(epLocsGen4_Mlt5);
			epLocsGen4_Mlt5.setInfo("(circles)");
			plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE, 1f, Color.RED));
		}
		DefaultXY_DataSet epLocsGen5_Mlt5 = getEpicenterLocsXY_DataSet(2.0, 5.0, 5);
		if(epLocsGen5_Mlt5.size()>0) {
			funcs.add(epLocsGen5_Mlt5);
			epLocsGen5_Mlt5.setInfo("(circles)");
			plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE, 1f, Color.ORANGE));
		}
		DefaultXY_DataSet epLocsGen6_Mlt5 = getEpicenterLocsXY_DataSet(2.0, 5.0, 6);
		if(epLocsGen6_Mlt5.size()>0) {
			funcs.add(epLocsGen6_Mlt5);
			epLocsGen6_Mlt5.setInfo("(circles)");
			plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CIRCLE, 1f, Color.YELLOW));
		}


		// 5.0<=M<6.5
		DefaultXY_DataSet epLocsGen1_Mgt5lt65 = getEpicenterLocsXY_DataSet(5.0, 6.5, 1);
		if(epLocsGen1_Mgt5lt65.size()>0) {
			funcs.add(epLocsGen1_Mgt5lt65);
			epLocsGen1_Mgt5lt65.setInfo("(triangles)");
			plotChars.add(new PlotCurveCharacterstics(PlotSymbol.TRIANGLE, 4f, Color.BLACK));
		}
		DefaultXY_DataSet epLocsGen2_Mgt5lt65 = getEpicenterLocsXY_DataSet(5.0, 6.5, 2);
		if(epLocsGen2_Mgt5lt65.size()>0) {
			funcs.add(epLocsGen2_Mgt5lt65);
			epLocsGen2_Mgt5lt65.setInfo("(triangles)");
			plotChars.add(new PlotCurveCharacterstics(PlotSymbol.TRIANGLE, 4f, Color.BLUE));
		}
		DefaultXY_DataSet epLocsGen3_Mgt5lt65 = getEpicenterLocsXY_DataSet(5.0, 6.5, 3);
		if(epLocsGen3_Mgt5lt65.size()>0) {
			funcs.add(epLocsGen3_Mgt5lt65);
			epLocsGen3_Mgt5lt65.setInfo("(triangles)");
			plotChars.add(new PlotCurveCharacterstics(PlotSymbol.TRIANGLE, 4f, Color.GREEN));
		}
		DefaultXY_DataSet epLocsGen4_Mgt5lt65 = getEpicenterLocsXY_DataSet(5.0, 6.5, 4);
		if(epLocsGen4_Mgt5lt65.size()>0) {
			funcs.add(epLocsGen4_Mgt5lt65);
			epLocsGen4_Mgt5lt65.setInfo("(triangles)");
			plotChars.add(new PlotCurveCharacterstics(PlotSymbol.TRIANGLE, 4f, Color.RED));
		}
		DefaultXY_DataSet epLocsGen5_Mgt5lt65 = getEpicenterLocsXY_DataSet(5.0, 6.5, 5);
		if(epLocsGen5_Mgt5lt65.size()>0) {
			funcs.add(epLocsGen5_Mgt5lt65);
			epLocsGen5_Mgt5lt65.setInfo("(triangles)");
			plotChars.add(new PlotCurveCharacterstics(PlotSymbol.TRIANGLE, 4f, Color.ORANGE));
		}
		DefaultXY_DataSet epLocsGen6_Mgt5lt65 = getEpicenterLocsXY_DataSet(5.0, 6.5, 6);
		if(epLocsGen6_Mgt5lt65.size()>0) {
			funcs.add(epLocsGen6_Mgt5lt65);
			epLocsGen6_Mgt5lt65.setInfo("(triangles)");
			plotChars.add(new PlotCurveCharacterstics(PlotSymbol.TRIANGLE, 4f, Color.YELLOW));
		}


		// 6.5<=M<9.0
		DefaultXY_DataSet epLocsGen1_Mgt65 = getEpicenterLocsXY_DataSet(6.5, 9.0, 1);
		if(epLocsGen1_Mgt65.size()>0) {
			funcs.add(epLocsGen1_Mgt65);
			epLocsGen1_Mgt65.setInfo("(squares)");
			plotChars.add(new PlotCurveCharacterstics(PlotSymbol.SQUARE, 8f, Color.LIGHT_GRAY));
		}
		DefaultXY_DataSet epLocsGen2_Mgt65 = getEpicenterLocsXY_DataSet(6.5, 9.0, 2);
		if(epLocsGen2_Mgt65.size()>0) {
			funcs.add(epLocsGen2_Mgt65);
			epLocsGen2_Mgt65.setInfo("(squares)");
			plotChars.add(new PlotCurveCharacterstics(PlotSymbol.SQUARE, 8f, Color.BLUE));
		}
		DefaultXY_DataSet epLocsGen3_Mgt65 = getEpicenterLocsXY_DataSet(6.5, 9.0, 3);
		if(epLocsGen3_Mgt65.size()>0) {
			funcs.add(epLocsGen3_Mgt65);
			epLocsGen3_Mgt65.setInfo("(squares)");
			plotChars.add(new PlotCurveCharacterstics(PlotSymbol.SQUARE, 8f, Color.GREEN));
		}
		DefaultXY_DataSet epLocsGen4_Mgt65 = getEpicenterLocsXY_DataSet(6.5, 9.0, 4);
		if(epLocsGen4_Mgt65.size()>0) {
			funcs.add(epLocsGen4_Mgt65);
			epLocsGen4_Mgt65.setInfo("(squares)");
			plotChars.add(new PlotCurveCharacterstics(PlotSymbol.SQUARE, 8f, Color.RED));
		}
		DefaultXY_DataSet epLocsGen5_Mgt65 = getEpicenterLocsXY_DataSet(6.5, 9.0, 5);
		if(epLocsGen5_Mgt65.size()>0) {
			funcs.add(epLocsGen5_Mgt65);
			epLocsGen5_Mgt65.setInfo("(squares)");
			plotChars.add(new PlotCurveCharacterstics(PlotSymbol.SQUARE, 8f, Color.ORANGE));
		}
		DefaultXY_DataSet epLocsGen6_Mgt65 = getEpicenterLocsXY_DataSet(6.5, 9.0, 6);
		if(epLocsGen6_Mgt65.size()>0) {
			funcs.add(epLocsGen6_Mgt65);
			epLocsGen6_Mgt65.setInfo("(squares)");
			plotChars.add(new PlotCurveCharacterstics(PlotSymbol.SQUARE, 8f, Color.YELLOW));
		}
		
		double minLat=90, maxLat=-90,minLon=360,maxLon=-360;
		for(AbstractXY_DataSet func:funcs) {
			System.out.println(func.getMinX()+"\t"+func.getMaxX()+"\t"+func.getMinY()+"\t"+func.getMaxY());
			if(func.getMaxX()>maxLon) maxLon = func.getMaxX();
			if(func.getMinX()<minLon) minLon = func.getMinX();
			if(func.getMaxY()>maxLat) maxLat = func.getMaxY();
			if(func.getMinY()<minLat) minLat = func.getMinY();
		}
		
		System.out.println("latDada\t"+minLat+"\t"+maxLat+"\t"+minLon+"\t"+maxLon+"\t");
		
		FaultTrace trace = mainShock.getRuptureSurface().getEvenlyDiscritizedUpperEdge();
		ArbitrarilyDiscretizedFunc traceFunc = new ArbitrarilyDiscretizedFunc();
		traceFunc.setName("Main Shock Trace");
		for(Location loc:trace)
			traceFunc.set(loc.getLongitude(), loc.getLatitude());
		funcs.add(traceFunc);
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.MAGENTA));
		
		GraphWindow graph = new GraphWindow(funcs, "Aftershock Epicenters for "+info); 
		graph.setX_AxisLabel("Longitude");
		graph.setY_AxisLabel("Latitude");
		double deltaLat = maxLat-minLat;
		double deltaLon = maxLon-minLon;
		double aveLat = (minLat+maxLat)/2;
		double scaleFactor = 1.57/Math.cos(aveLat*Math.PI/180);	// this is what deltaLon/deltaLat should equal
		if(deltaLat > deltaLon/scaleFactor) {	// expand lon range
			double newLonMax = minLon + deltaLat*scaleFactor;
			graph.setX_AxisRange(minLon, newLonMax);
			graph.setY_AxisRange(minLat, maxLat);
		}
		else { // expand lat range
			double newMaxLat = minLat + deltaLon/scaleFactor;
			graph.setX_AxisRange(minLon, maxLon);
			graph.setY_AxisRange(minLat, newMaxLat);
		}
		graph.setPlotChars(plotChars);
		if(savePDF_File)
		try {
			graph.saveAsPDF(dirToSaveData+"epicenterMap.pdf");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	
	/**
	 * This writes aftershock data to a file
	 */
	public void writeAftershockDataToFile(ArrayList<PrimaryAftershock> events, String filePathAndName) {
		try{
			FileWriter fw1 = new FileWriter(filePathAndName);
			fw1.write("id\ttime\tlat\tlon\tdepth\tmag\tParentID\tGeneration\tdistToMainShock\tdistToParent\tERF_srcID\tERF_rupID\tERF_SrcName\n");
			for(PrimaryAftershock event: events) {
				Location hLoc = event.getHypocenterLocation();
				double dist = LocationUtils.distanceToSurfFast(hLoc, mainShock.getRuptureSurface());
				fw1.write(event.getID()+"\t"+(float)event.getOriginTime()+"\t"
						+(float)hLoc.getLatitude()+"\t"
						+(float)hLoc.getLongitude()+"\t"
						+(float)hLoc.getDepth()+"\t"
						+(float)event.getMag()+"\t"
						+event.getParentID()+"\t"
						+event.getGeneration()+"\t"
						+(float)dist+"\t"
						+(float)event.getDistanceToParent()+"\t"
						+event.getERF_SourceIndex()+"\t"
						+event.getERF_RupIndex()+"\t"
						+erf.getSource(event.getERF_SourceIndex()).getName()+"\n");
			}
			fw1.close();
		}catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	
	
	/**
	 * This creates an EqksInGeoBlock for the given ERF at each point in the GriddedRegion region
	 */
	public void makeAllEqksInGeoBlocks() {

		double calcStartTime=System.currentTimeMillis();
		System.out.println("Starting to make blocks");

		blockList = new ArrayList<EqksInGeoBlock>();
		for(Location loc: griddedRegion) {
			EqksInGeoBlock block = new EqksInGeoBlock(loc,griddedRegion.getSpacing(),0,16);
			blockList.add(block);
		}
		System.out.println("Number of Blocks: "+blockList.size()+" should be("+griddedRegion.getNodeCount()+")");

		double forecastDuration = erf.getTimeSpan().getDuration();
		double rateUnAssigned = 0;
		int numSrc = erf.getNumSources();
		for(int s=0;s<numSrc;s++) {
			ProbEqkSource src = erf.getSource(s);
			
			// skip the source if it's to be ignored (blunt elastic rebound application)
			if(s==sourceID_ToIgnore) {
				System.out.println("Ignoring source "+s+" ("+src.getName()+") in the sampling of events");
				continue;
			}

			int numRups = src.getNumRuptures();
			for(int r=0; r<numRups;r++) {
				ProbEqkRupture rup = src.getRupture(r);
				ArbDiscrEmpiricalDistFunc numInEachNode = new ArbDiscrEmpiricalDistFunc(); // node on x-axis and num on y-axis
				LocationList locsOnRupSurf = rup.getRuptureSurface().getEvenlyDiscritizedListOfLocsOnSurface();
				double rate = rup.getMeanAnnualRate(forecastDuration);
				int numUnAssigned=0;
				for(Location loc: locsOnRupSurf) {
					int nodeIndex = griddedRegion.indexForLocation(loc);
					if(nodeIndex != -1)
						numInEachNode.set((double)nodeIndex,1.0);
					else
						numUnAssigned +=1;
				}
				int numNodes = numInEachNode.size();
				if(numNodes>0) {
					for(int i=0;i<numNodes;i++) {
						int nodeIndex = (int)Math.round(numInEachNode.getX(i));
						double fracInside = numInEachNode.getY(i)/locsOnRupSurf.size();
						double nodeRate = rate*fracInside;	// fraction of rate in node
						blockList.get(nodeIndex).processRate(nodeRate, fracInside, s, r, rup.getMag());
					}
				}
				float fracUnassigned = (float)numUnAssigned/(float)locsOnRupSurf.size();
				if(numUnAssigned>0) System.out.println(fracUnassigned+" of rup "+r+" were unassigned for source "+s+" ("+erf.getSource(s).getName()+")");
				rateUnAssigned += rate*fracUnassigned;
			}
		}

		System.out.println("rateUnAssigned = "+rateUnAssigned);

		double runtime = (System.currentTimeMillis()-calcStartTime)/1000;
		System.out.println("Making blocks took "+runtime+" seconds");

		// This checks to make sure total rate in all blocks (plus rate unassigned) is equal the the total ERF rate
		System.out.println("TESTING RESULT");
		double testRate1=0;
		for(EqksInGeoBlock block: blockList) {
			testRate1+=block.getTotalRateInside();
		}
		testRate1+=rateUnAssigned;
		double testRate2=0;
		for(int s=0;s<numSrc;s++) {
			ProbEqkSource src = erf.getSource(s);
			int numRups = src.getNumRuptures();
			for(int r=0; r<numRups;r++) {
				testRate2 += src.getRupture(r).getMeanAnnualRate(forecastDuration);
			}
		}
		System.out.println("\tRate1="+(float)testRate1+" should equal Rate2="+(float)testRate2+";\tratio="+(float)(testRate1/testRate2));
	//	System.out.println("\tRate2="+testRate2);
		
	}


	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		String rootDir = "/Users/field/workspace/OpenSHA/dev/scratch/ned/ETAS_Tests/computedData/";
		
		// make the gridded region
		CaliforniaRegions.RELM_GRIDDED griddedRegion = new CaliforniaRegions.RELM_GRIDDED();
		
		// write the expected number of primary aftershocks for different mags
		ETAS_Utils utils = new ETAS_Utils();
		System.out.println("Num primary aftershocks for mag:");
		for(int i=2;i<9;i++) {
			System.out.println("\t"+i+"\t"+utils.getDefaultExpectedNumEvents((double)i, 0, 365));
		}
		
		// for keeping track of runtime:
		long startRunTime=System.currentTimeMillis();

		// Create the UCERF2 instance
		System.out.println("Starting ERF instantiation");
		double forecastDuration = 1.0;	// years
		ModMeanUCERF2 meanUCERF2 = new ModMeanUCERF2();
		meanUCERF2.setParameter(UCERF2.RUP_OFFSET_PARAM_NAME, new Double(10.0));
		meanUCERF2.getParameter(UCERF2.PROB_MODEL_PARAM_NAME).setValue(UCERF2.PROB_MODEL_POISSON);
//		meanUCERF2.setParameter(UCERF2.BACK_SEIS_NAME, UCERF2.BACK_SEIS_ONLY);
		meanUCERF2.setParameter(UCERF2.BACK_SEIS_NAME, UCERF2.BACK_SEIS_INCLUDE);
//		meanUCERF2.setParameter(UCERF2.BACK_SEIS_NAME, UCERF2.BACK_SEIS_EXCLUDE);
		meanUCERF2.setParameter(UCERF2.BACK_SEIS_RUP_NAME, UCERF2.BACK_SEIS_RUP_POINT);
		meanUCERF2.getTimeSpan().setDuration(forecastDuration);
		meanUCERF2.updateForecast();
		double runtime = (System.currentTimeMillis()-startRunTime)/1000;
		System.out.println("ERF instantiation took "+runtime+" seconds");
		
		double cumRate = Math.round(10*ERF_Calculator.getTotalMFD_ForERF(meanUCERF2, 2.05,8.95, 70, true).getCumRate(5.05))/10.0;
		// get a slightly different result using ERF_Calculator.getMagFreqDistInRegion(meanUCERF2, CaliforniaRegions.RELM_GRIDDED(),5.05,35,0.1, true), but with rounding it's the same
		System.out.println("\nCumRate >= M5 for MeanUCERF2_ETAS: "+cumRate+"  (should be 7.5)\n");
		
//		meanUCERF2.plotMFD_InRegionNearLanders(10.0);
		
//		ProbEqkSource source = meanUCERF2.getSource(1708);
//		for(int r=0; r<source.getNumRuptures();r++)
//			System.out.println("test 1708\t"+source.getName()+"\t"+(float)source.getRupture(r).getMag()+"\t"+(float)source.getRupture(r).getMeanAnnualRate(forecastDuration));

				
		
		
		/*
		// print info about sources
		for(int s=0;s<meanUCERF2.getNumSources();s++) {
			System.out.println(s+"\t"+meanUCERF2.getSource(s).getName()+"\t #rups="+
				meanUCERF2.getSource(s).getNumRuptures()+"\t"+meanUCERF2.getSource(s).getRupture(0).getRuptureSurface().getAveDip());
			if(s==223) {
				ProbEqkSource src = meanUCERF2.getSource(s);
				for(int r=0; r<src.getNumRuptures(); r++ )
					System.out.println("\tLanders rup "+r+"; M="+src.getRupture(r).getMag());
			}
		}
		*/
		

		// Create the ETAS simulator
//		int srcID_ToIgnore = -1; // none
        int srcID_ToIgnore = 195; // landers
//		int srcID_ToIgnore = 223; // Northridge
		// SSAF is 42 to 98
		// Landers is 195
		// Northridge is 223
		ETAS_Simulator etasSimulator = new ETAS_Simulator(meanUCERF2,griddedRegion,srcID_ToIgnore);
		
//		etasSimulator.plotERF_MagFreqDists();
		
		ProbEqkRupture mainShock;
		
		// NOW RUN TESTS FOR VARIOUS MAIN SHOCKS
/*	
		//  Point source at edge of sub-blocks
		mainShock = new ProbEqkRupture();
		mainShock.setMag(5);
		mainShock.setPointSurface(new Location(34,-118,8));
		etasSimulator.runTests(mainShock,"M5 Pt Src at 34,-118,8 (sub-block edge)",null);
		
		// Point source in midpoint of sub block
		mainShock.setPointSurface(new Location(34.00625,-118.00625,7));
		etasSimulator.runTests(mainShock,"M5 Pt Src at 34.00625,-118.00625,7 (sub-block mid)",null);

		// point source in center of sub block
		mainShock.setPointSurface(new Location(34.0125,	-118.0125,	6.0));
		etasSimulator.runTests(mainShock,"M5 Pt Src at 34.0125,-118.0125,6.0 (sub-block center)",null);

		// 68	S. San Andreas;CH+CC+BB+NM+SM+NSB+SSB+BG+CO	 #rups=8
		mainShock = meanUCERF2.getSource(68).getRupture(4);
		etasSimulator.runTests(mainShock,"SSAF Wall-to-wall Rupture; M="+mainShock.getMag(), null);

		
		FaultTrace trace = new FaultTrace("Test");
		trace.add(new Location (34,-117));
		trace.add(new Location (35,-117));
		StirlingGriddedSurface rupSurf = new StirlingGriddedSurface(trace, 90.0, 0.0,16.0, 1.0);
		mainShock = new ProbEqkRupture(7.0,0,1,rupSurf,null);
		double distDecay = 1.7;
		double minDist = 0.3;
		String info = "Test Event (M="+mainShock.getMag()+"); distDecay="+distDecay;
		etasSimulator.runTempTest(mainShock,info, 195, rootDir+"TestM7/",distDecay,minDist);
*/		

		// 195	Landers Rupture	 #rups=46 (rup 41 for M 7.25)
		mainShock = meanUCERF2.getSource(195).getRupture(41);
		double distDecay = 1.7;
		double minDist = 0.3;
		String info = "Landers Rupture (M="+mainShock.getMag()+"); distDecay="+distDecay;
		etasSimulator.runTests(mainShock,info, 195, rootDir+"Landers_decay1pt7_withOutSrc5/",distDecay,minDist);
//		etasSimulator.runTempTest(mainShock,info, 195, rootDir+"TestLanders/",distDecay,minDist);
/*
				
		// 223	Northridge	 #rups=13	(rup 8 for M 6.75)
		mainShock = meanUCERF2.getSource(223).getRupture(8);
		double distDecay = 1.7;
		double minDist = 0.3;
		String info = "Northridge Rupture (M="+mainShock.getMag()+"); distDecay="+distDecay;
		etasSimulator.runTests(mainShock,info, 223, rootDir+"Northridge_decay1pt7_withSrc/",distDecay,minDist);
		
		

		// 236	Pitas Point (Lower, West)	 #rups=19	13.0 (shallowest dipping rupture I could find)
		mainShock = meanUCERF2.getSource(236).getRupture(10);
		String info = "Pitas Point (shallowest dipping source); M="+mainShock.getMag()+"; AveDip="+mainShock.getRuptureSurface().getAveDip();
		etasSimulator.runTests(mainShock,info,null);
*/

		runtime = (System.currentTimeMillis()-startRunTime)/1000;
		System.out.println("Test Run took "+runtime+" seconds");
		
	}
	
	
	public ArrayList<PrimaryAftershock> getAllAftershocks(PrimaryAftershock mainShock) {
		int generation = 0;
		totalExpectedNumMagDist = new SummedMagFreqDist(2.05, 8.95, 70);
		ArrayList<PrimaryAftershock> mainShocksToProcess = new ArrayList<PrimaryAftershock>();
		ArrayList<PrimaryAftershock> allEvents = new ArrayList<PrimaryAftershock>();
		mainShocksToProcess.add(mainShock);
//		allEvents.add(mainShock);
		int numToProcess = mainShocksToProcess.size();
		
		while(numToProcess > 0) {
			generation +=1;
			System.out.println("WORKING ON GENERATION: "+generation+ "  (numToProcess="+numToProcess+")");
			ArrayList<PrimaryAftershock> aftershockList = new ArrayList<PrimaryAftershock>();
			CalcProgressBar progressBar = new CalcProgressBar("ETAS_Simulator", "Progress on Generation "+generation);
			progressBar.displayProgressBar();
			for(int i=0;i<numToProcess; i++) {
				progressBar.updateProgress(i,numToProcess);
//				System.out.print(numToProcess-i);
				aftershockList.addAll(getPrimaryAftershocksList(mainShocksToProcess.get(i)));
				// add to the total expected MFD
				ArbIncrementalMagFreqDist expMFD = etas_sampler.getMagProbDist();
				expMFD.scale(expectedNum);
				totalExpectedNumMagDist.addIncrementalMagFreqDist(expMFD);
			}	
			
			// set the IDs & generation
			int firstID = allEvents.size();
			for(int i=0;i<aftershockList.size(); i++) {
				PrimaryAftershock aShock = aftershockList.get(i);
				aShock.setID(i+firstID);
				aShock.setGeneration(generation);
			}
			allEvents.addAll(aftershockList);
			mainShocksToProcess = aftershockList;
//			numToProcess = 0;
			numToProcess = mainShocksToProcess.size();
			progressBar.dispose();

		}
		System.out.println("Total Num Aftershocks="+allEvents.size());

		return allEvents;

	}
	
	
	
	public void runTempTest(ProbEqkRupture mainShock, String info, Integer srcIndex, String dirName, double distDecay, double minDist) {
		this.distDecay = distDecay;
		this.minDist = minDist;
		this.mainShock = mainShock;
		includeBlockRates=false;
		
		// Make the ETAS sampler for the given main shock:
//		ETAS_PrimaryEventSamplerTest sampler = new ETAS_PrimaryEventSamplerTest(mainShock,blockList, erf, distDecay,minDist, useAdaptiveBlocks, includeBlockRates);
		
		System.out.println("Starting old way");
		ETAS_PrimaryEventSamplerJUNK2 sampler = new ETAS_PrimaryEventSamplerJUNK2(mainShock,blockList, erf, distDecay,minDist, useAdaptiveBlocks, includeBlockRates);
		sampler.plotBlockProbMap(info+" - old way", true, "testOldWay");
		sampler.plotDistDecayTestFuncs("OldWay", "testOldWay/distDecayPlot");
		
		System.out.println("Starting new way");
		ETAS_PrimaryEventSamplerTest sampler2 = new ETAS_PrimaryEventSamplerTest(mainShock,blockList, erf, distDecay,minDist, useAdaptiveBlocks, includeBlockRates);
		sampler2.plotBlockProbMap(info+" - new way", true, "testNewWay");
		sampler.plotDistDecayTestFuncs("NewWay", "testNewWay/distDecayPlot");

	
		IntegerPDF_FunctionSampler probsOld = sampler.getRandomBlockSampler();
		IntegerPDF_FunctionSampler probsNew = sampler2.getRandomBlockSampler();

		System.out.println("should be 1: "+probsOld.calcSumOfY_Vals());
		System.out.println("should be 1: "+probsNew.calcSumOfY_Vals());

		// PLOT Ratio
		
		ArrayList<EqksInGeoBlock> revisedBlockList = sampler.getRevisedBlockList();
		GMT_MapGenerator mapGen = new GMT_MapGenerator();
		mapGen.setParameter(GMT_MapGenerator.GMT_SMOOTHING_PARAM_NAME, false);
		mapGen.setParameter(GMT_MapGenerator.TOPO_RESOLUTION_PARAM_NAME, GMT_MapGenerator.TOPO_RESOLUTION_NONE);
		mapGen.setParameter(GMT_MapGenerator.MIN_LAT_PARAM_NAME,31.5);		// -R-125.4/-113.0/31.5/43.0
		mapGen.setParameter(GMT_MapGenerator.MAX_LAT_PARAM_NAME,43.0);
		mapGen.setParameter(GMT_MapGenerator.MIN_LON_PARAM_NAME,-125.4);
		mapGen.setParameter(GMT_MapGenerator.MAX_LON_PARAM_NAME,-113.0);
		mapGen.setParameter(GMT_MapGenerator.LOG_PLOT_NAME,false);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MODE_NAME,GMT_MapGenerator.COLOR_SCALE_MODE_MANUALLY);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME,0.0);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME,4.0);

		CaliforniaRegions.RELM_GRIDDED griddedRegion = new CaliforniaRegions.RELM_GRIDDED();
		GriddedGeoDataSet xyzDataSet = new GriddedGeoDataSet(griddedRegion, true);
		
		double[] newProbs = new double[griddedRegion.getNodeCount()];
		double[] oldProbs = new double[griddedRegion.getNodeCount()];
		
		for(int i=0; i<revisedBlockList.size();i++) {
			EqksInGeoBlock block = revisedBlockList.get(i);
			int locIndex = griddedRegion.indexForLocation(block.getBlockCenterLoc());
			newProbs[locIndex] += probsNew.getY(i);
			oldProbs[locIndex] += probsOld.getY(i);
		}
		
		// set ratio
		for(int i=0; i<xyzDataSet.size();i++) xyzDataSet.set(i, newProbs[i]/oldProbs[i]);

//		System.out.println("Min & Max Z: "+xyzDataSet.getMinZ()+"\t"+xyzDataSet.getMaxZ());
		String metadata = "";
		
		try {
			String name = mapGen.makeMapLocally(xyzDataSet, "Ratio", metadata, "testRatio");
//			metadata += GMT_MapGuiBean.getClickHereHTML(mapGen.getGMTFilesWebAddress());
//			ImageViewerWindow imgView = new ImageViewerWindow(name,metadata, true);
//			System.out.println("GMT Plot Filename: "+name);
		} catch (GMT_MapException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (RuntimeException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		/*		*/	

	}

	
	/**
	 * This runs some tests for the given main shock
	 * @param mainShock
	 * @param info - plotting label
	 * @param srcIndex - the index of a source to compare results with
	 */
	public void runTests(ProbEqkRupture mainShock, String info, Integer srcIndex, String dirName, double distDecay, double minDist) {
		
		this.distDecay = distDecay;
		this.minDist = minDist;
		
		dirToSaveData = dirName;
		
		this.mainShock = mainShock;
				
		// make the directory
		File file1 = new File(dirName);
		file1.mkdirs();

		// convert mainShock to PrimaryAftershock so it has the origin time
		PrimaryAftershock mainShockConverted = new PrimaryAftershock(mainShock);
		mainShockConverted.setOriginTime(tMin);
		mainShockConverted.setID(-1);
		mainShockConverted.setGeneration(0);
		
		infoForOutputFile = new String();
		infoForOutputFile += "Run for "+info+"\n\n";
		infoForOutputFile += "sourceID_ToIgnore = "+sourceID_ToIgnore;
		if(sourceID_ToIgnore != -1)
			infoForOutputFile += " ("+erf.getSource(sourceID_ToIgnore).getName()+")\n\n";
		else
			infoForOutputFile += " (null)\n\n";

		infoForOutputFile += "ERF Parameter Settings for "+erf.getName()+": "+erf.getAdjustableParameterList().toString()+"\n\n";
		
		infoForOutputFile += "ETAS Param Settings:\n"+
			"\n\tdistDecay="+distDecay+
			"\n\tminDist="+minDist+
			"\n\ttMin="+tMin+
			"\n\ttMax="+tMax;
		for(String str: etasUtils.getDefaultParametersAsStrings()) {
			infoForOutputFile += "\n\t"+str;
		}
		infoForOutputFile += "\n\n";

		// get the aftershocks
		long startRunTime=System.currentTimeMillis();
		allAftershocks = getAllAftershocks(mainShockConverted);
		double runtime = (System.currentTimeMillis()-startRunTime)/1000;
		infoForOutputFile += "Generating all aftershocks took "+runtime+" seconds ("+(float)(runtime/60)+" minutes)\n\n";

		
		System.out.println("Num aftershocks = "+allAftershocks.size());
		infoForOutputFile += "Num aftershocks = "+allAftershocks.size()+"\n\n";

		// get the number of aftershocks in each generation
		int[] numEventsInEachGeneration = new int[10];
		for(PrimaryAftershock aShock:allAftershocks) {
			int gen = aShock.getGeneration();
			numEventsInEachGeneration[gen] += 1;
		}
		String numInGenString = "Num events in each generation:\n"+
				"\n\t1st\t"+numEventsInEachGeneration[1]+
				"\n\t2nd\t"+numEventsInEachGeneration[2]+
				"\n\t3rd\t"+numEventsInEachGeneration[3]+
				"\n\t4th\t"+numEventsInEachGeneration[4]+
				"\n\t5th\t"+numEventsInEachGeneration[5]+
				"\n\t6th\t"+numEventsInEachGeneration[6]+
				"\n\t7th\t"+numEventsInEachGeneration[7]+
				"\n\t8th\t"+numEventsInEachGeneration[8]+
				"\n\t9th\t"+numEventsInEachGeneration[9];
		System.out.println(numInGenString);
		infoForOutputFile += numInGenString+"\n\n";
		
		
		// write out any occurrences of srcIndex aftershocks:
		infoForOutputFile += erf.getSource(srcIndex).getName()+" (srcId="+srcIndex+") aftershocks:\n\n";
		ArrayList<PrimaryAftershock> srcID_aftershocks = new ArrayList<PrimaryAftershock>();
		for(PrimaryAftershock aShock:allAftershocks)
			if(aShock.getERF_SourceIndex() == srcIndex)
				srcID_aftershocks.add(aShock);
		if(srcID_aftershocks.size()==0)
			infoForOutputFile += "\tNone\n";
		else {
			infoForOutputFile += "\tid\tmag\tgen\tparID\trupID\tsrcID\toriginTime\tdistToParent\n";
			for(PrimaryAftershock aShock:srcID_aftershocks) {
				infoForOutputFile += "\t"+aShock.getID()+
				"\t"+aShock.getMag()+
				"\t"+aShock.getGeneration()+
				"\t"+aShock.getParentID()+
				"\t"+aShock.getERF_RupIndex()+
				"\t"+aShock.getERF_SourceIndex()+
				"\t"+(float)aShock.getOriginTime()+
				"\t"+(float)aShock.getDistanceToParent()+"\n";
			}
		}
		infoForOutputFile += "\n";
		
		// write out event larger that M 6.5
		infoForOutputFile += "Aftershocks larger than M 6.5:\n\n";
		ArrayList<PrimaryAftershock> largeAftershocks = new ArrayList<PrimaryAftershock>();
		for(PrimaryAftershock aShock:allAftershocks)
			if(aShock.getMag()>6.5)
				largeAftershocks.add(aShock);
		if(largeAftershocks.size()==0)
			infoForOutputFile += "\tNone\n";
		else {
			infoForOutputFile += "\tid\tmag\tgen\tparID\trupID\tsrcID\toriginTime\tdistToParent\tsrcName\n";
			for(PrimaryAftershock aShock:largeAftershocks) {
				infoForOutputFile += "\t"+aShock.getID()+
				"\t"+(float)aShock.getMag()+
				"\t"+aShock.getGeneration()+
				"\t"+aShock.getParentID()+
				"\t"+aShock.getERF_RupIndex()+
				"\t"+aShock.getERF_SourceIndex()+
				"\t"+(float)aShock.getOriginTime()+
				"\t"+(float)aShock.getDistanceToParent()+
				"\t"+erf.getSource(aShock.getERF_SourceIndex()).getName()+"\n";
			}
		}
		infoForOutputFile += "\n";

		
		String infoText = plotMagFreqDists(srcIndex, info, true);

		infoForOutputFile += infoText;	// this writes the total num expected above M 6.5, plus that from the srcIndex
		
		writeAftershockDataToFile(allAftershocks, dirToSaveData+"eventList.txt");
				
		plotDistDecayForAshocks(info, true, mainShock);
		
		plotEpicenterMap(info, true, mainShock);
		
		String mapMetadata = etas_FirstGenSampler.plotBlockProbMap(info, true,"test");
		infoForOutputFile += "\n"+mapMetadata+".  The allFiles.zip should have been"+
				" moved to this dir, uncompressed, and the name changed from allFiles to samplingProbMaps.\n\n";
		
		plotNumVsTime(info, true, mainShock);
		
		
		//write info file
		try{
			FileWriter fw = new FileWriter(dirToSaveData+"INFO.txt");
			fw.write(infoForOutputFile);
			fw.close();
		}catch(Exception e) {
			e.printStackTrace();
		}
		
		try {
			ETAS_Utils.writeEQCatFile(new File(dirToSaveData+"catalog.sc"), allAftershocks);
		} catch (IOException e) {
			e.printStackTrace();
		}
			
	}
	
	
	
	

	
	
	
	
	public void plotDistDecayForAshocks(String info, boolean savePDF_File, ProbEqkRupture mainShock) {
		
		double delta = 10;
		ArrayList<EvenlyDiscretizedFunc> distDecayFuncs = etas_FirstGenSampler.getDistDecayTestFuncs(delta);
		distDecayFuncs.get(0).setName("Approx Expected Distance Decay for Primary Aftershocks of "+info);
		distDecayFuncs.get(0).setInfo("Diff from theoretical mostly due to no events to sample outside RELM region, but also spatially variable a-values");
		distDecayFuncs.get(1).setName("Theoretical Distance Decay");
		distDecayFuncs.get(1).setInfo("(dist+minDist)^-distDecay, where minDist="+minDist+" and distDecay="+distDecay+", and where finite discretization accounted for");
		EvenlyDiscretizedFunc tempFunc = distDecayFuncs.get(0);
		EvenlyDiscretizedFunc obsPrimaryDistHist = new EvenlyDiscretizedFunc(delta/2, tempFunc.size(), tempFunc.getDelta());
		obsPrimaryDistHist.setTolerance(tempFunc.getTolerance());
		EvenlyDiscretizedFunc obsAllDistHist = new EvenlyDiscretizedFunc(delta/2, tempFunc.size(), tempFunc.getDelta());
		obsAllDistHist.setTolerance(tempFunc.getTolerance());
		double totAllNum = 0, totPrimaryNum = 0;
		for (PrimaryAftershock event : allAftershocks) {
			if(event.getGeneration()==1) {
				obsPrimaryDistHist.add(event.getDistanceToParent(), 1.0);
				totPrimaryNum += 1;	
				obsAllDistHist.add(event.getDistanceToParent(), 1.0);
				totAllNum += 1;
			}
			// get distance to mainshock for the rest of the ruptures
			else {
				double dist = LocationUtils.distanceToSurfFast(event.getHypocenterLocation(), mainShock.getRuptureSurface());
				obsAllDistHist.add(dist, 1.0);
				totAllNum += 1;
			}
		}
		for(int i=0; i<obsPrimaryDistHist.size();i++) {
			obsPrimaryDistHist.set(i, obsPrimaryDistHist.getY(i)/totPrimaryNum);		// convert to PDF
			obsAllDistHist.set(i, obsAllDistHist.getY(i)/totAllNum);					// convert to PDF
		}
		obsPrimaryDistHist.setName("Sampled Distance-Decay Histogram for Primary Aftershocks of "+info);
		obsPrimaryDistHist.setInfo("(filled circles)");
		distDecayFuncs.add(obsPrimaryDistHist);
		obsAllDistHist.setName("Sampled Distance-Decay Histogram for All Aftershocks of "+info);
		obsAllDistHist.setInfo("(Crosses, and these are distances to the main shock, not to the parent)");
		distDecayFuncs.add(obsAllDistHist);

		GraphWindow graph = new GraphWindow(distDecayFuncs, "Distance Decay for Aftershocks of "+info); 
		graph.setX_AxisLabel("Distance (km)");
		graph.setY_AxisLabel("Fraction of Aftershocks");
		graph.setX_AxisRange(0.4, 1200);
		graph.setY_AxisRange(1e-6, 1);
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.FILLED_CIRCLE, 3f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, Color.GREEN));
		graph.setPlotChars(plotChars);
		graph.setYLog(true);
		graph.setXLog(true);
		if(savePDF_File)
		try {
			graph.saveAsPDF(dirToSaveData+"primaryAshocksDistDecay.pdf");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

}
