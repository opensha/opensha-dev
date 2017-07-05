package scratch.ned.ETAS_ERF.testModels;

import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.calc.magScalingRelations.MagLengthRelationship;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.mapping.gmt.GMT_MapGenerator;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.mapping.gmt.gui.GMT_MapGuiBean;
import org.opensha.commons.mapping.gmt.gui.ImageViewerWindow;
import org.opensha.commons.param.impl.CPTParameter;
import org.opensha.commons.util.FileUtils;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.calc.ERF_Calculator;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.earthquake.param.AleatoryMagAreaStdDevParam;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.griddedSeis.NSHMP_GridSourceGenerator;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.griddedSeis.Point2Vert_FaultPoisSource;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.sha.magdist.ArbIncrementalMagFreqDist;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import java.awt.Toolkit;
import java.io.File;
import java.io.IOException;

import scratch.UCERF3.analysis.GMT_CA_Maps;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.FaultSystemSolutionPoissonERF;
import scratch.UCERF3.utils.UCERF3_DataUtils;
import scratch.UCERF3.utils.ModUCERF2.NSHMP_GridSourceGeneratorMod2;


/**
 * Note that this does not yet include C zones (fixed strike sources)
 * @author field
 *
 */
public class TestModel1_ERF extends FaultSystemSolutionERF {
	
	final static boolean D=true;

	GriddedRegion griddedRegion;
	
	double minGridLat=35.0;
	double maxGridLat=37.0;
	double minGridLon=240-360;
	double maxGridLon=244-360;
	double gridSpacing=0.05;
//	double gridSpacing=0.1;
	
	ArbIncrementalMagFreqDist  onFaultPointMFD, offFaultPointMFD;
	
	ArrayList<Integer> locIndicesOnFault;
	
	WC1994_MagLengthRelationship magLengthRel;
	
	public TestModel1_ERF() {
		super(new TestModel1_FSS());
		
//		System.out.println("SIMULATION_MODE="+SIMULATION_MODE);
		
		griddedRegion = new GriddedRegion(new Location(minGridLat,minGridLon), new Location(maxGridLat,maxGridLon), gridSpacing, gridSpacing, null);
		
		numOtherSources = griddedRegion.getNumLocations();
//		numOtherSources=0;
		System.out.println("numOtherSources="+numOtherSources);
//		System.out.println(griddedRegion.getLocation(0));
//		System.out.println(griddedRegion.getLocation(1));
		
		magLengthRel = new WC1994_MagLengthRelationship();
		
		// get the point sources that are on the fault
		LocationList pointLocsOnFault = ((TestModel1_FSS)this.getSolution()).getPointLocsOnFault();
		locIndicesOnFault = new ArrayList<Integer>();
		for(Location loc : pointLocsOnFault){
			int index = griddedRegion.indexForLocation(loc);
			locIndicesOnFault.add(index);
//			Location testLoc = griddedRegion.getLocation(index);
//			System.out.println(loc+"\t"+testLoc);
		}
		

		// the following is the target MFD (down to M 2.5)
		GutenbergRichterMagFreqDist  targetFaultGR = ((TestModel1_FSS)this.getSolution()).getTargetFaultGR();
		// the following is the MFD for the fault (seismogenic and larger)
		ArbIncrementalMagFreqDist faultGR = ((TestModel1_FSS)this.getSolution()).getFaultGR();
		
		double offFaultSeisReductionFactor = 1;
//		double offFaultSeisReductionFactor = 10;
		int numPtsOnFault = pointLocsOnFault.size();
		if(D) System.out.println("numPtsOnFault+"+numPtsOnFault);
		onFaultPointMFD  = new ArbIncrementalMagFreqDist(2.55, 6.05, 36);
		for(int i=0; i<onFaultPointMFD.size();i++)
			onFaultPointMFD.set(i,targetFaultGR.getY(i)/numPtsOnFault);
		// make point off fault 1/10th rate of those on the faults
		offFaultPointMFD = new ArbIncrementalMagFreqDist(2.55, 7.55, 51);
		for(int i=0; i<offFaultPointMFD.size();i++)
			offFaultPointMFD.set(i,targetFaultGR.getY(i)/(numPtsOnFault*offFaultSeisReductionFactor));
		
		
		if(D) {
			ArrayList<EvenlyDiscretizedFunc> funcs = new ArrayList<EvenlyDiscretizedFunc>();
			funcs.add(faultGR);
			faultGR.setName("faultGR");
			faultGR.setInfo(" ");		
			funcs.add(faultGR.getCumRateDistWithOffset());
			funcs.add(targetFaultGR);
			targetFaultGR.setName("targetFaultGR");
			targetFaultGR.setInfo(" ");		
			funcs.add(targetFaultGR.getCumRateDistWithOffset());
			funcs.add(onFaultPointMFD);
			onFaultPointMFD.setName("onFaultPointMFD");
			onFaultPointMFD.setInfo(" ");		
			funcs.add(onFaultPointMFD.getCumRateDistWithOffset());
			funcs.add(offFaultPointMFD);
			offFaultPointMFD.setName("offFaultPointMFD");
			offFaultPointMFD.setInfo(" ");		
			funcs.add(offFaultPointMFD.getCumRateDistWithOffset());
			GraphWindow graph = new GraphWindow(funcs, ""); 	
			graph.setYLog(true);
		}


	}
	
	public GriddedRegion getGriddedRegion() {
		return griddedRegion;
	}
	
	public ArbIncrementalMagFreqDist getOnFaultMFD() {
		return onFaultPointMFD;
	}
	
	@Override
	protected ProbEqkSource getOtherSource(int iSource) {
		int regionIndex = iSource;
		ArbIncrementalMagFreqDist mfd;
		if(locIndicesOnFault.contains(regionIndex))
			mfd = onFaultPointMFD;
		else
			mfd = offFaultPointMFD;
		double magCutOff =8.0; // all rups below are treated a point sources
		boolean isCrossHair=true;
		return new Point2Vert_FaultPoisSource(griddedRegion.getLocation(regionIndex), mfd,
				magLengthRel,timeSpan.getDuration(), magCutOff ,1.0, 0.0,0.0, isCrossHair);
	}


	
	private void makeNucleationMap(double minMag, double maxMag) throws IOException {
		
		String scaleLabel="TestMod1 Nucl";
		String metadata=""; 
		String dirName="TestMod1 Nucl"; 
		
		GMT_MapGenerator gmt_MapGenerator = GMT_CA_Maps.getDefaultGMT_MapGenerator();
		gmt_MapGenerator.setParameter(GMT_MapGenerator.MIN_LAT_PARAM_NAME, minGridLat);
		gmt_MapGenerator.setParameter(GMT_MapGenerator.MIN_LON_PARAM_NAME, minGridLon);
		gmt_MapGenerator.setParameter(GMT_MapGenerator.MAX_LAT_PARAM_NAME, maxGridLat);
		gmt_MapGenerator.setParameter(GMT_MapGenerator.MAX_LON_PARAM_NAME, maxGridLon);
		gmt_MapGenerator.setParameter(GMT_MapGenerator.GRID_SPACING_PARAM_NAME, gridSpacing);
		gmt_MapGenerator.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME, -6.0);
		gmt_MapGenerator.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME, 2.0);
		
		// must set this parameter this way because the setValue(CPT) method takes a CPT object, and it must be the
		// exact same object as in the constraint (same instance); the setValue(String) method was added for convenience
		// but it won't succeed for the isAllowed(value) call.
		CPTParameter cptParam = (CPTParameter )gmt_MapGenerator.getAdjustableParamsList().getParameter(GMT_MapGenerator.CPT_PARAM_NAME);
		cptParam.setValue(GMT_CPT_Files.MAX_SPECTRUM.getFileName());

		
		File GMT_DIR = new File(UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, "GMT");
		
		GriddedGeoDataSet geoDataSet = ERF_Calculator.getNucleationRatesInRegion(this, griddedRegion, minMag, maxMag);

		try {
				if (!GMT_DIR.exists())
					GMT_DIR.mkdir();
				String url = gmt_MapGenerator.makeMapUsingServlet(geoDataSet, scaleLabel, metadata, "TestModel1_ERF_"+minMag);
				metadata += GMT_MapGuiBean.getClickHereHTML(gmt_MapGenerator.getGMTFilesWebAddress());
				File downloadDir = new File(GMT_DIR, dirName);
				if (!downloadDir.exists())
					downloadDir.mkdir();
				File zipFile = new File(downloadDir, "allFiles.zip");
				// construct zip URL
				String zipURL = url.substring(0, url.lastIndexOf('/')+1)+"allFiles.zip";
				FileUtils.downloadURL(zipURL, zipFile);
				FileUtils.unzipFile(zipFile, downloadDir);
				
				ImageViewerWindow imgView = new ImageViewerWindow(url,metadata, true);
		} catch (GMT_MapException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	
	/**
	 * I commented out the guts of this code because it uses the old sampler (ETAS_PrimaryEventSamplerAlt), which I wanted
	 * to remove from the repository, and I didn't have time to update the implementation here.
	 * @param griddedRegion
	 * @param mainShock
	 * @param rthRup
	 * @param includeEqkRates
	 * @param magThresh
	 */
	public void calcSelfTriggeringPob(GriddedRegion griddedRegion, ObsEqkRupture mainShock, int rthRup, boolean includeEqkRates, 
			double magThresh) {
		
//		if(!SIMULATION_MODE)
//			throw new RuntimeException("This method can only be run if SIMULATION_MODE = true");
//
//		System.out.println("Updating forecast (twice)");
//		// get the total rate over the duration of the forecast
//		updateForecast();	// do this to get annual rate over the entire forecast (used to sample spontaneous events)
//		double origTotRate = totalRate;	// this include ER time dependence, but diff shouldn't be noticeable.
//		System.out.println("origTotRate="+origTotRate);
//		
//		// set to yearly probabilities for simulation forecast (in case input was not a 1-year forecast)
//		timeSpan.setDuration(1.0);	// annualize
//		updateForecast();
//
//		
//		Region regionForRates = new Region(griddedRegion.getBorder(),BorderType.MERCATOR_LINEAR);
//		// first make array of rates for each source
//		double sourceRates[] = new double[getNumSources()];
//		double duration = getTimeSpan().getDuration();
//		for(int s=0;s<getNumSources();s++)
//			sourceRates[s] = getSource(s).computeTotalEquivMeanAnnualRate(duration);
//
//		System.out.println("making ETAS_PrimaryEventSamplerAlt");
//		ETAS_PrimaryEventSamplerAlt etas_PrimEventSampler = new ETAS_PrimaryEventSamplerAlt(regionForRates, this, sourceRates, 
//				gridSpacing, null, includeEqkRates);
//		
//		
//		
////		System.out.println("Plotting expected MFD");
////		GraphWindow graph = new GraphWindow(etas_PrimEventSampler.getExpectedMFD(mainShock), "Expected MFD for MainShock"); 
//
//		// Plot the map
//		System.out.println("Making sampler map");
//		etas_PrimEventSampler.plotSamplerMap(etas_PrimEventSampler.getAveSamplerForRupture(mainShock), "testMod1", "testMod1_Map");
//
//		
////		etas_PrimEventSampler.testRates();
//		
//		System.out.println("getting TriggerProbOfEachSource");
//		double[] srcTrigProb = etas_PrimEventSampler.getTriggerProbOfEachSource(mainShock);
//		
//		System.out.println("getting rupsThatOverlap");
//		List<Integer> rupsThatOverlap = ((TestModel1_FSS)invSol).getRupsThatOverlapGivenRup(rthRup, 11);
//		
//		System.out.println("rupsThatOverlap.size()="+rupsThatOverlap.size());
//		
////		int indexOfMagThresh = offFaultPointMFD.getClosestXIndex(magThresh);
////		System.out.println(magThresh+"\t"+indexOfMagThresh);
////		System.out.println(offFaultPointMFD);
//		
//		
//		// compute expected and overlapping MFDs
//		SummedMagFreqDist magDist = new SummedMagFreqDist(2.05, 8.95, 70);
//		SummedMagFreqDist rupsThatOverlapMFD = new SummedMagFreqDist(2.05, 8.95, 70);
//		double testFault =0, testFaultSame =0, testOffFault=0;
//		for(int s=0; s<srcTrigProb.length;s++) {
//			SummedMagFreqDist srcMFD = ERF_Calculator.getTotalMFD_ForSource(this.getSource(s), 1.0, 2.05, 8.95, 70, true);
//			srcMFD.normalizeByTotalRate();
//			srcMFD.scale(srcTrigProb[s]);
//			if(!Double.isNaN(srcMFD.getTotalIncrRate())) { // not sure why this is needed
//				magDist.addIncrementalMagFreqDist(srcMFD);
//				// test
//				if(s<numNonZeroFaultSystemSources) {
//					testFault += srcMFD.getCumRate(magThresh);
//					if(rupsThatOverlap.contains(s)) {
//						testFaultSame += srcMFD.getCumRate(magThresh);
//					}
//				}
//				else {
//					testOffFault += srcMFD.getCumRate(magThresh);
//				}
//				
//				// end test
//				if(rupsThatOverlap.contains(s)) {
//					rupsThatOverlapMFD.addIncrementalMagFreqDist(srcMFD);
//				}
//			}
//		}
//		System.out.println("\ttestFault="+testFault+"\ttestFaultSame="+testFaultSame+"\ttestOffFault="+testOffFault
//				+"\ttestTotal="+(testOffFault+testFault));
//		ArrayList<EvenlyDiscretizedFunc> funcs = new ArrayList<EvenlyDiscretizedFunc>();
//		funcs.add(magDist.getCumRateDistWithOffset());
//		funcs.add(rupsThatOverlapMFD.getCumRateDistWithOffset());
//		GraphWindow graph2 = new GraphWindow(funcs, "Expected MFD for MainShock"); 
//
//
//		
//
//		System.out.println("computing probabilities");
//		double ptSrcProbAboveMagThresh = offFaultPointMFD.getCumRate(magThresh)/offFaultPointMFD.getTotalIncrRate();
////		GraphWindow graph3 = new GraphWindow(offFaultPointMFD.getCumRateDistWithOffset(), "Expected MFD for MainShock"); 
//
////		System.out.println("offFaultPointMFD.getTotalIncrRate()="+offFaultPointMFD.getTotalIncrRate());
////		System.out.println("offFaultPointMFD.getCumRate(indexOfMagThresh)="+offFaultPointMFD.getCumRate(magThresh));
//		System.out.println("magThresh="+magThresh);
//		System.out.println("ptSrcProbAboveMagThresh="+ptSrcProbAboveMagThresh);
//
//		double sameEventProb = 0;
//		double totalLargeEventProb = 0;
//		double totalProb = 0;
//		double totalPtSrcLargeEventProb=0;
//		double totalLargeFaultProb=0;
//		for(int s=0; s<srcTrigProb.length; s++) {
//			totalProb += srcTrigProb[s];
//			if(s<numNonZeroFaultSystemSources) {
//				if(getSource(s).getNumRuptures() != 1)  throw new RuntimeException("Problem");	// check to make sure there is only one rupture
//				if(getSource(s).getRupture(0).getMag() >= (magThresh-0.05)) {	// 0.05 is half the bin width
//					totalLargeEventProb += srcTrigProb[s];
//					totalLargeFaultProb += srcTrigProb[s];
//					if(rupsThatOverlap.contains(s)) {
//						sameEventProb += srcTrigProb[s];
//					}
//				}
//			}
//			else {
//				// these are all gridded sources
//				if(!locIndicesOnFault.contains(s-numNonZeroFaultSystemSources)){	// make sure this is not a fault location (that would be double counting)
//					totalLargeEventProb += srcTrigProb[s]*ptSrcProbAboveMagThresh;
//					totalPtSrcLargeEventProb += srcTrigProb[s]*ptSrcProbAboveMagThresh;				
//				}
//			}
//		}
//		
//		System.out.println("totalProb="+(float)totalProb+"\t(should be ~1.0)");
//		System.out.println("totalLargeEventProb="+(float)totalLargeEventProb);
//		System.out.println("totalPtSrcLargeEventProb="+(float)totalPtSrcLargeEventProb);
//		System.out.println("sameEventProb="+(float)sameEventProb);
//		System.out.println("totalLargeFaultProb="+(float)totalLargeFaultProb);
//		System.out.println("(sameEventProb/totalLargeEventProb)="+(float)(sameEventProb/totalLargeEventProb));

	}

	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		TestModel1_ERF erf = new TestModel1_ERF();
		erf.getParameter(AleatoryMagAreaStdDevParam.NAME).setValue(0.0);
//		erf.aleatoryMagAreaStdDevParam.setValue(0.0);
//		erf.bpt_AperiodicityParam.setValue(0.2);
		erf.getTimeSpan().setStartTimeInMillis(0);
		erf.getTimeSpan().setDuration(1);
		long runtime = System.currentTimeMillis();
		
		// update forecast
		erf.updateForecast();
		
		// print the nucleation rate map
//		try {
//			erf.makeNucleationMap(2.5,10.0);
//			erf.makeNucleationMap(6.5,10.0);
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
		
		int sthSrc = 892;	// same as source index
		ProbEqkRupture mainshock = erf.getNthRupture(sthSrc);		
		ObsEqkRupture obsMainShock = new ObsEqkRupture();
		obsMainShock.setAveRake(mainshock.getAveRake());
		obsMainShock.setMag(mainshock.getMag());
		obsMainShock.setRuptureSurface(mainshock.getRuptureSurface());
//		obsMainShock.setPointSurface(mainshock.getRuptureSurface().getFirstLocOnUpperEdge());
		obsMainShock.setOriginTime(0);	// occurs at 1970
		System.out.println("main shock: nthRup="+sthSrc+"; mag="+obsMainShock.getMag()+
				"; src name: " +erf.getSource(sthSrc).getName());

		// this applies elastic rebound reduction of probability
//		erf.setRuptureOccurrence(sthSrc, 0);
		erf.calcSelfTriggeringPob(erf.getGriddedRegion(), obsMainShock, sthSrc, false, 6.15);
		
		
		// this is for test simulations
//		ArrayList<ObsEqkRupture> obsEqkRuptureList = new ArrayList<ObsEqkRupture>();
//		obsEqkRuptureList.add(obsMainShock);
//		erf.setRuptureOccurrence(sthSrc, 0);
//		erf.testETAS_Simulation(erf.getGriddedRegion(), obsEqkRuptureList,true, true, false,0.05);

//		erf.testER_Simulation();
//		runtime -= System.currentTimeMillis();
//		System.out.println("simulation took "+(double)runtime/(1000.0*60.0)+" minutes");
	}
}
