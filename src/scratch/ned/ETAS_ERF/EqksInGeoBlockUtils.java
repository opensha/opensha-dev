package scratch.ned.ETAS_ERF;

import java.util.ArrayList;

import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.earthquake.AbstractERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;

import scratch.UCERF3.erf.FaultSystemSolutionPoissonERF;

public class EqksInGeoBlockUtils {

	/**
	 * @param args
	 */
	public static void main(String[] args) {

		int numSlice =6;
		double minLat = 36.0;
		double maxLat = 36.0 + 0.1/numSlice;
		double minLon = 241-360;
		double maxLon = 241-360 + 0.1/numSlice; 
		double minDepth = 8.0;
		double maxDepth = 8.0 + 24/(numSlice/2);
		
		double distDecay = 2;
		double minDist = 0.3;
		int numDiscr = 100;
		
//		Location loc = new Location(minLat-0.04,minLon-0.04,minDepth-4);	// at corner
		Location loc = new Location(minLat-0.03,minLon-0.03,minDepth-3);	// 8.9 km from center

		EqksInGeoBlock block = new EqksInGeoBlock(minLat, maxLat, minLon, maxLon, minDepth, maxDepth);
		
		double dist = LocationUtils.linearDistanceFast(loc, block.getBlockCenterLoc());
		double wt1 = Math.pow(dist+minDist, -distDecay);
		double dist2 = EqksInGeoBlockUtils.getEquivDistForBlock(block, loc, distDecay, minDist, numDiscr);
		double wt2 =  Math.pow(dist2+minDist, -distDecay);
		System.out.println("Bias: "+wt1/wt2);
		dist2 = EqksInGeoBlockUtils.getEquivDistForBlockFast(block, loc, distDecay, minDist, numDiscr);
		wt2 =  Math.pow(dist2+minDist, -distDecay);
		System.out.println("BiasFast: "+wt1/wt2);

	}
	
	/**
	 * This creates an EqksInGeoBlock for the given ERF at each point in the GriddedRegion region
	 */
	public static ArrayList<EqksInGeoBlock> makeAllEqksInGeoBlocks(FaultSystemSolutionPoissonERF erf, GriddedRegion griddedRegion, double maxDepth) {

		double calcStartTime=System.currentTimeMillis();
		System.out.println("Starting to make blocks");

		ArrayList<EqksInGeoBlock> blockList = new ArrayList<EqksInGeoBlock>();
		for(Location loc: griddedRegion) {
			EqksInGeoBlock block = new EqksInGeoBlock(loc,griddedRegion.getSpacing(),0,maxDepth);
			blockList.add(block);
		}
		System.out.println("Number of Blocks: "+blockList.size()+" should be("+griddedRegion.getNodeCount()+")");

		double forecastDuration = erf.getTimeSpan().getDuration();
		double rateUnAssigned = 0;
		int numSrc = erf.getNumSources();
		double maxRupDepth = 0;
		for(int s=0;s<numSrc;s++) {
			ProbEqkSource src = erf.getSource(s);
			int numRups = src.getNumRuptures();
			for(int r=0; r<numRups;r++) {
				ProbEqkRupture rup = src.getRupture(r);
				ArbDiscrEmpiricalDistFunc numInEachNode = new ArbDiscrEmpiricalDistFunc(); // node on x-axis and num on y-axis
				LocationList locsOnRupSurf = rup.getRuptureSurface().getEvenlyDiscritizedListOfLocsOnSurface();
				double rate = rup.getMeanAnnualRate(forecastDuration);
				int numUnAssigned=0;
				for(Location loc: locsOnRupSurf) {
					if(maxRupDepth<loc.getDepth()) maxRupDepth = loc.getDepth();
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
						double fracInside = numInEachNode.getY(i)/(double)locsOnRupSurf.size();
						double nodeRate = rate*fracInside;	// fraction of rate in node
						int nthRup = erf.getIndexN_ForSrcAndRupIndices(s, r);
						blockList.get(nodeIndex).processRate(nodeRate, fracInside, nthRup, rup.getMag());
					}
				}
				float fracUnassigned = (float)numUnAssigned/(float)locsOnRupSurf.size();
				if(numUnAssigned>0) 
					System.out.println(fracUnassigned+" (fraction) of rup "+r+" were unassigned for source "+s+" ("+erf.getSource(s).getName()+")");
				rateUnAssigned += rate*fracUnassigned;
			}
		}
		
		if(maxRupDepth > maxDepth) {
			throw new RuntimeException("ruptures go deeper than the the given maxDepth:\tmaxRupDepth="+maxRupDepth);
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
		
//		System.out.println("LONGER TEST OF BLOCKS");
//		
//		System.out.println("Testing block 630");
//		blockList.get(630).testThisBlock(erf);
//		
//		
//		int index=0;
//		for(EqksInGeoBlock block: blockList) {
//			System.out.println(index);
//			block.testThisBlock(erf);
//			index +=1;
//		}

		
		return blockList;
		
	}
	
	/**
	 * This checks whether the total rate in sub-block lists is equal to the rate in the main block
	 * @param erf
	 * @param blockList
	 * @param subBlockList1
	 * @param subBlockList2
	 */
	public static void testSubBlockListRates(AbstractERF erf, ArrayList<EqksInGeoBlock> blockList, 
			ArrayList<ArrayList<EqksInGeoBlock>> subBlockList1, ArrayList<ArrayList<EqksInGeoBlock>> subBlockList2) {
		
//		System.out.println("Testing sub-block rates:\t");
		
		boolean blockList1_Tested = false;
		boolean blockList2_Tested = false;
		double totalRateOverAllBlocks =0;
		
		// test blockList1
		for(int b=0; b<blockList.size(); b++) {
			double blockRate=blockList.get(b).getTotalRateInside();
			totalRateOverAllBlocks += blockRate;
			ArrayList<EqksInGeoBlock> blocks1 = subBlockList1.get(b);
			if(blocks1 != null) {
				blockList1_Tested = true;
				double totRate=0;
				for(EqksInGeoBlock blk:blocks1) {
					totRate += blk.getTotalRateInside();
				}
				double ratio = totRate/blockRate;
				if(Math.abs(totRate) < 1e-15 && Math.abs(blockRate) < 1e-15 )
					ratio = 1;
				if(ratio <0.999 || ratio > 1.001) {
					throw new RuntimeException("Descrepancy for block "+b+" of blockList1;\tratio="+ratio+
							"\ttotRate="+totRate+"\tblockRate="+blockRate);
				}

				
			}
			
			ArrayList<EqksInGeoBlock> blocks2 = subBlockList2.get(b);
			if(blocks2 != null) {
				blockList2_Tested = true;
				double totRate=0;
				for(EqksInGeoBlock blk:blocks2) {
					totRate += blk.getTotalRateInside();
				}
				double ratio = totRate/blockRate;
				if(Math.abs(totRate) < 1e-15 && Math.abs(blockRate) < 1e-15 )
					ratio = 1;
				if(ratio <0.999 || ratio > 1.001) {
					throw new RuntimeException("Descrepancy for block "+b+" of blockList2;\tratio="+ratio+
							"\ttotRate="+totRate+"\tblockRate="+blockRate);
				}
			}
		}
//		System.out.println("blockList1 was testable = "+blockList1_Tested+
//				";\tblockList2 was testable = "+blockList2_Tested+
//				"\ttotalRateOverAllBlocks="+totalRateOverAllBlocks);
	}
	
	
	
	public static double getEquivDistForBlock(EqksInGeoBlock block, Location loc, double distDecay, double minDist, int numDiscr) {
		long startTime = System.currentTimeMillis();

		double totSum = 0;
		double deltaLat = (block.maxLat-block.minLat)/numDiscr;
		double deltaLon = (block.maxLon-block.minLon)/numDiscr;
		double deltaDepth = (block.maxDepth-block.minDepth)/numDiscr;
		for(int iLat = 0; iLat < numDiscr; iLat++) {
			double lat = block.minLat + iLat*deltaLat + deltaLat/2;
			for(int iLon = 0; iLon < numDiscr; iLon++) {
				double lon = block.minLon + iLon*deltaLon + deltaLon/2;
				for(int iDep = 0; iDep < numDiscr; iDep++) {
					double depth = block.minDepth + iDep*deltaDepth + deltaDepth/2;
					Location loc2 = new Location(lat,lon,depth);
					double dist = LocationUtils.linearDistanceFast(loc, loc2);
					totSum += Math.pow(dist+minDist, -distDecay);
				}
			}
		}

		totSum /= (double)(numDiscr*numDiscr*numDiscr);
		
		System.out.println("runTime="+ (System.currentTimeMillis()-startTime));

		return Math.pow(totSum,-1.0/distDecay)-minDist;
	}
	
	
	public static double getEquivDistForBlockFast(EqksInGeoBlock block, Location loc, double distDecay, double minDist, int numDiscr) {
//		long startTime = System.currentTimeMillis();
		double totSum = 0;
		double deltaLat = (block.maxLat-block.minLat)/numDiscr;
		double deltaLon = (block.maxLon-block.minLon)/numDiscr;
		double deltaDepth = (block.maxDepth-block.minDepth)/numDiscr;
		for(int iLat = 0; iLat < numDiscr; iLat++) {
			double distLat = (loc.getLatitude() - (block.minLat + iLat*deltaLat + deltaLat/2))*111.0;
			for(int iLon = 0; iLon < numDiscr; iLon++) {
				double distLon = (loc.getLongitude() - (block.minLon + iLon*deltaLon + deltaLon/2)) * 111.0 * Math.cos(loc.getLatitude()*Math.PI/180);
				for(int iDep = 0; iDep < numDiscr; iDep++) {
					double distDepth = loc.getDepth() - (block.minDepth + iDep*deltaDepth + deltaDepth/2);
					double dist = Math.sqrt(distLat*distLat+distLon*distLon+distDepth*distDepth);
					totSum += Math.pow(dist+minDist, -distDecay);
				}
			}
		}

		totSum /= (double)(numDiscr*numDiscr*numDiscr);

//		System.out.println("runTime="+ (System.currentTimeMillis()-startTime));

		return Math.pow(totSum,-1.0/distDecay)-minDist;

	}

	

}
