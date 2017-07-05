package scratch.ned.ETAS_ERF.tests;

import java.util.ArrayList;

import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.param.AleatoryMagAreaStdDevParam;

import scratch.UCERF3.erf.FaultSystemSolutionPoissonERF;
import scratch.ned.ETAS_ERF.EqksInGeoBlock;
import scratch.ned.ETAS_ERF.EqksInGeoBlockUtils;
import scratch.ned.ETAS_ERF.testModels.TestModel1_ERF;

public class Test_ETAS_PrimaryEventSampler {

	/**
	 * 
	 * CODE COMMENTED OUT DUE TO ERRORS AFTER CHANGING TestModel1_ERF to extend FaultSystemSolutionERF
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
////		CaliforniaRegions.RELM_GRIDDED griddedRegion = new CaliforniaRegions.RELM_GRIDDED();
//		
//		
//		double distDecay = 2;
//		double minDist = 0.3;
//		boolean useAdaptiveBlocks= true; 
//		boolean includeBlockRates = false;
//		double maxBlockDepth=24;
//
//		
//		
//		TestModel1_ERF erf = new TestModel1_ERF();
//		erf.getParameter(AleatoryMagAreaStdDevParam.NAME).setValue(0.0);
//		System.out.println("Running updateForecast()");
//		erf.updateForecast();
//		
//		// SET THE MAINSHOCK
//		ProbEqkRupture mainShock;
////		int nthRup = 892;	// same as source index
////		mainShock = erf.getNthRupture(nthRup);	
//		
//		//  Point source 
//		mainShock = new ProbEqkRupture();
//		mainShock.setMag(5);
////		mainShock.setPointSurface(new Location(36,242-360,8));  	// loc at edge of sub-block
//		mainShock.setPointSurface(new Location(36+0.05/12,(242-360)+0.05/12,9));  	// loc at middle of sub-block
//
//
//		
//		GriddedRegion griddedRegion = erf.getGriddedRegion();  
//		
//		// make the EqksInGeoBlock lists
//		System.out.println("Making initial EqksInGeoBlock lists");
//		ArrayList<EqksInGeoBlock> blockList = EqksInGeoBlockUtils.makeAllEqksInGeoBlocks(erf, griddedRegion, maxBlockDepth);
//		ArrayList<ArrayList<EqksInGeoBlock>> subBlockList1 =  new ArrayList<ArrayList<EqksInGeoBlock>>(); // intermediate level of subdivision
//		ArrayList<ArrayList<EqksInGeoBlock>> subBlockList2 =  new ArrayList<ArrayList<EqksInGeoBlock>>(); // highest level of subdivision
//		// populate the sub-block lists
//		for(int i=0;i<blockList.size();i++) {
//			subBlockList1.add(null);
//			subBlockList2.add(null);
//		}
//
//		long startTime = System.currentTimeMillis();
//		System.out.println("Making ETAS_PrimaryEventSampler");
//		ETAS_PrimaryEventSampler sampler = new ETAS_PrimaryEventSampler(mainShock, blockList, subBlockList1, subBlockList2, erf, distDecay, minDist, 
//				useAdaptiveBlocks, includeBlockRates);
//		System.out.println("... that took "+ ((System.currentTimeMillis()-startTime)/(1000))+" seconds");
//
//		
//		System.out.println("Making Dist Decay Plot");
//		sampler.plotDistDecayTestFuncs("test", null,2.0);

	}

}
