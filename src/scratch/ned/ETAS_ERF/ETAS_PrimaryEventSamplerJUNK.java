package scratch.ned.ETAS_ERF;

import java.awt.Color;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.xyz.ArbDiscrGeoDataSet;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
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
import org.opensha.sha.faultSurface.AbstractEvenlyGriddedSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.sha.magdist.ArbIncrementalMagFreqDist;

import scratch.UCERF3.erf.FaultSystemSolutionPoissonERF;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.ETAS_Utils;
import scratch.UCERF3.erf.ETAS.IntegerPDF_FunctionSampler;

/**
 * This class allows one to randomly sample which block a primary aftershock nucleates in 
 * based on the original block rates and distance decay.
 * @author field
 *
 */
public class ETAS_PrimaryEventSamplerJUNK {
	
	protected final static boolean D = true;  // for debugging
	
	ArrayList<EqksInGeoBlock> origBlockList;				// the original list of blocks given
	ArrayList<EqksInGeoBlock> revisedBlockList;	// a revised list (those used here)
	ArrayList<Double> revisedBlockDistances;	// distances to parentRup for the revised blocks
	EqkRupture parentRup;						// the main shock
	FaultSystemSolutionPoissonERF erf;							// the ERF from which primary aftershocks will be sampled
	double[] relBlockProb;						// the relative probability that each block will nucleate an event (REDUNDANT WITH randomBlockSampler BELOW)
	int numBlocks;
	IntegerPDF_FunctionSampler randomBlockSampler;	// used to randomly sample blocks
	final static double ADAPT_BLOCK_DIST1 = 8;		// first distance threshold for sampling blocks more densely
	final static double ADAPT_BLOCK_DIST2 = 16;		// second distance threshold for sampling blocks more densely
	double distDecay, minDist;
	double minBlockDist=Double.MAX_VALUE;
	int minDistIndex=-1;
	boolean includeBlockRates;
	int numToProcess;
	
	
	
	/**
	 * 
	 * @param parentRup - the main shock
	 * @param blockList - the list of blocks to sample from
	 * @param erf - the ERF from which primary aftershocks are sampled from
	 * @param distDecay - distance decay parameter
	 * @param minDist - value that prevents singularity at zero distance
	 * @param useAdaptiveBlocks - indicate whether blocks close to parentRup should be sub-sampled
	 * @param includeBlockRates - whether or not to use rates inside of blocks to modify spatial probabilities
	 */
	public ETAS_PrimaryEventSamplerJUNK(EqkRupture parentRup, ArrayList<EqksInGeoBlock> blockList, 
			ArrayList<ArrayList<EqksInGeoBlock>> subBlockList1, ArrayList<ArrayList<EqksInGeoBlock>> subBlockList2,
			FaultSystemSolutionPoissonERF erf, double distDecay, double minDist, 
			boolean useAdaptiveBlocks, boolean includeBlockRates) {
		
		this.parentRup=parentRup;
		origBlockList = blockList;
		this.erf = erf;
		this.distDecay = distDecay;
		this.minDist = minDist;
		this.includeBlockRates = includeBlockRates;
		
		RuptureSurface rupSurf = parentRup.getRuptureSurface();
		
		revisedBlockList = new ArrayList<EqksInGeoBlock>();  // revised is for replacing blocks with sub-blocks close in
		ArrayList<EqksInGeoBlock> subBlocks = null;
		revisedBlockDistances = new ArrayList<Double>();

		if(D) System.out.print("computing block distances & sub-blocks");
		int counter=0, counterThresh = origBlockList.size()/20, counterIncr=counterThresh;
		minBlockDist=Double.MAX_VALUE;
		minDistIndex=-1;
		// compute distance from rupture to center of each block (and subdivide close blocks if specified)
		for(int b=0; b<origBlockList.size(); b++) {
			EqksInGeoBlock origBlock = origBlockList.get(b);
			counter++;
			// write progress
			if(counter>counterThresh) {
				double perc = 100*counter/origBlockList.size();
				if(D) System.out.print(", "+(int)perc);
				counterThresh += counterIncr;
			}
			double dist;
			if(rupSurf.isPointSurface())
				dist = EqksInGeoBlockUtils.getEquivDistForBlockFast(origBlock, parentRup.getRuptureSurface().getFirstLocOnUpperEdge(), distDecay, minDist, 10);	
			else
				dist = LocationUtils.distanceToSurfFast(origBlock.getBlockCenterLoc(), rupSurf);

			if(dist>ADAPT_BLOCK_DIST2  || !useAdaptiveBlocks) {
				revisedBlockDistances.add(dist);
				revisedBlockList.add(origBlock);
				if(dist<minBlockDist) {
					minBlockDist=dist;
					minDistIndex = revisedBlockDistances.size()-1;
				}
			}
			else {  // apply adaptive block sizes
				if (dist > ADAPT_BLOCK_DIST1) {
					subBlocks = subBlockList1.get(b);
					if(subBlocks == null) {
						subBlocks = origBlock.getSubBlocks(3, 6, erf);
						if(parentRup.getMag()>5)
							subBlockList1.set(b,subBlocks);
					}
				}
				else {
					subBlocks = subBlockList2.get(b);
					if(subBlocks == null) {
						subBlocks = origBlock.getSubBlocks(6, 12, erf);
						if(parentRup.getMag()>5)
							subBlockList2.set(b,subBlocks);
					}
//System.out.println("getAveBlockSize()="+subBlocks.get(0).getAveBlockSize());
//System.exit(0);
				}
				double testRate1=0;
				for(EqksInGeoBlock subBlock:subBlocks) {
					testRate1 += subBlock.getTotalRateInside();
//					double dist2 = LocationUtils.distanceToSurfFast(subBlock.getBlockCenterLoc(), rupSurf);
					
					double dist2;
					if(rupSurf.isPointSurface())
						dist2 = EqksInGeoBlockUtils.getEquivDistForBlockFast(subBlock, parentRup.getRuptureSurface().getFirstLocOnUpperEdge(), distDecay, minDist, 10);	
					else
						dist2 = LocationUtils.distanceToSurfFast(subBlock.getBlockCenterLoc(), rupSurf);

					revisedBlockDistances.add(dist2);
					revisedBlockList.add(subBlock);
					if(dist2<minBlockDist) {
						minBlockDist=dist2;
						minDistIndex = revisedBlockDistances.size()-1;
					}
				}
				// Test that block rates are the same
				double testRate2=origBlock.getTotalRateInside();
				double ratio = testRate1/testRate2;
				if(Math.abs(testRate1) < 1e-15 && Math.abs(testRate2) < 1e-15 )
					ratio =1;
				if(ratio<0.9999 || ratio>1.0001) {
					throw new RuntimeException("PROBLEM: ratio="+ratio+";\ttestRate1="+testRate1+"\ttestRate2="+testRate2);
				}
			}
		}
		numBlocks = revisedBlockList.size();
		if(D) System.out.print(" ");
		
//		getSpatialProbDensOfBlocks();
		
		upDataRandomBlockSampler();
	}
	
	/**
	 * This returns the parent rupture (main shocks)
	 * @return
	 */
	public EqkRupture getParentRup() {
		return parentRup;
	}
	
	
	/**
	 * This updates the randomBlockSampler
	 */
	public void upDataRandomBlockSampler() {

		// compute relative probability of each block
		relBlockProb = new double[numBlocks];
		
		double rupLength = parentRup.getRuptureSurface().getAveLength();
		
		if(D) {
			System.out.println("\nnum revised blocks="+numBlocks);
			System.out.println("minBlockDist ="+minBlockDist+" for block index "+minDistIndex);
			if(parentRup.getRuptureSurface().isPointSurface()) {
				Location ptLoc = parentRup.getRuptureSurface().getFirstLocOnUpperEdge();
				Location blLoc = this.revisedBlockList.get(minDistIndex).getBlockCenterLoc();
				System.out.println("\tdLat="+(float)(ptLoc.getLatitude()-blLoc.getLatitude())+
						"\tdLon="+(float)(ptLoc.getLongitude()-blLoc.getLongitude())+
						"\tdDep="+(float)(ptLoc.getDepth()-blLoc.getDepth()));

			}
			System.out.println("computing relative block probabilities");
			System.out.println("rupLength="+rupLength);
		}
		
		double total=0;
		if(rupLength>5.0) {	// a non point source rupture
			for(int i=0; i<numBlocks;i++) {
				double vol = revisedBlockList.get(i).getBlockVolume();
				double dist = revisedBlockDistances.get(i);
				double blockWt = vol/(2*Math.PI*(dist+minDist)+2*rupLength); // fraction of the annulus
				relBlockProb[i] = Math.pow(revisedBlockDistances.get(i)+minDist, -distDecay)*blockWt;
				if(includeBlockRates) relBlockProb[i] *= revisedBlockList.get(i).getTotalRateInside()/vol;
				total += relBlockProb[i];
			}
		}
		else {	// Point-source case
			double closestBlockVal=0;
//			double check = 1.0;
			double check = revisedBlockList.get(minDistIndex).getAveBlockSize()/2;
			if(minBlockDist>check) // check whether nucleation point is closer to the center or edge of closest block; value of "1" depends on smallest block size!
				minDistIndex=-1;  // don't do special case below
			for(int i=0; i<numBlocks;i++) {
				double vol = revisedBlockList.get(i).getBlockVolume();
				double dist = revisedBlockDistances.get(i);
				double blockWt;
				if(i == minDistIndex) {  // this will calculate the fraction directly
					double radius = Math.pow(vol/(4*Math.PI),0.33333);	// 4*PI*r^3;
					closestBlockVal=ETAS_Utils.getDecayFractionInsideDistance(distDecay, minDist, radius);
					relBlockProb[i] = closestBlockVal;
					if(D) System.out.println("Calculated wt of closest block directly; it equals: "+closestBlockVal);
//					double altWt = getBlockWeight(parentRup.getRuptureSurface().getFirstLocOnUpperEdge(), revisedBlockList.get(minDistIndex), distDecay, minDist, 10);
//					if(D) System.out.println("Calculated wt from getBlockWeight: "+altWt);
				}
				else {
//					if(dist <8.0)
//						blockWt = vol/(4*Math.PI*(dist*dist));
//					else
						blockWt = vol/(2*Math.PI*(dist));
					relBlockProb[i] = Math.pow(dist+minDist, -distDecay)*blockWt; 
					total += relBlockProb[i];
				}
			}
			// normalize
			for(int i=0; i<numBlocks;i++)
				if(i != minDistIndex)
					relBlockProb[i] *= (1-closestBlockVal)/total;	// normalize so total is 1.0
			total=1;
			// add block rates
			if(includeBlockRates) {
				total = 0;
				for(int i=0; i<numBlocks;i++) {
					double vol = revisedBlockList.get(i).getBlockVolume();
					relBlockProb[i] *= revisedBlockList.get(i).getTotalRateInside()/vol;
					total += relBlockProb[i];
				}
			}
		}
		
		// final normalization and to create the random block sampler
		randomBlockSampler = new IntegerPDF_FunctionSampler(numBlocks);
		for(int i=0; i<numBlocks;i++) {
			relBlockProb[i] /= total;
			randomBlockSampler.set(i,relBlockProb[i]);
		}
		if(D) {
			System.out.println("check that randomBlockSampler values sum to 1.0:\t"+randomBlockSampler.calcSumOfY_Vals());
			System.out.println("Done computing relative block probabilities");
		}
	}
	
	
	/**
	 * This updates the randomBlockSampler
	 */
	private void getSpatialProbDensOfBlocks() {

		// compute relative probability of each block
		relBlockProb = new double[numBlocks];
		
		ArbDiscrEmpiricalDistFunc distanceDist = new ArbDiscrEmpiricalDistFunc();

		for(int i=0; i<numBlocks;i++) {
			double aveBlockSize = revisedBlockList.get(i).getAveBlockSize();
			double dist = revisedBlockDistances.get(i);
			double blockWt;
			distanceDist.set(dist, 1.0);
			distanceDist.set(dist+aveBlockSize/3, 1.0);
			distanceDist.set(dist-aveBlockSize/3, 1.0);
		}
		
		ArrayList funcs = new ArrayList();
		funcs.add(distanceDist);
		funcs.add(distanceDist.getCumDist());
		GraphWindow graph = new GraphWindow(funcs, "test"); 
		graph.setAxisRange(1, 1200, 1e-6, 1);
		graph.setYLog(true);
//		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
//		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
//		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.BLUE));
//		graph.setPlottingFeatures(plotChars);
//		graph.setX_AxisLabel("Distance (km)");
//		graph.setY_AxisLabel("Probability");

			
	}

	
	
	/**
	 * This returns a randomly sampled block for a primary aftershock
	 * @return
	 */
	public int sampleRandomBlockIndex() {
		return randomBlockSampler.getRandomInt();
	}
	
	
	/** This computes the expected mag-freq dist for the primary aftershock sequence
	 * 
	 * @return
	 */
	public ArbIncrementalMagFreqDist getMagProbDist() {
		ArbIncrementalMagFreqDist magDist = new ArbIncrementalMagFreqDist(2.05, 8.95, 70);
		for(int i=0; i<revisedBlockList.size();i++) {
			ArbIncrementalMagFreqDist blockMagProbDist = revisedBlockList.get(i).getMagProbDist();
			double blockProb=randomBlockSampler.getY(i);
			for(int j=0; j<blockMagProbDist.size(); j++)
				magDist.addResampledMagRate(blockMagProbDist.getX(j), blockProb*blockMagProbDist.getY(j), true);
		}
		// now normalize
		double total = magDist.calcSumOfY_Vals();
		magDist.scaleToCumRate(2.05, 1.0/total);
		return magDist;
	}
	
	
	/** This computes the expected mag-freq dist for a particular source in the primary 
	 * aftershock sequence
	 * 
	 * @return
	 */
	public ArbIncrementalMagFreqDist getMagProbDistForSource(int srcIndex) {
		ArbIncrementalMagFreqDist magDist = new ArbIncrementalMagFreqDist(2.05, 8.95, 70);
		for(int i=0; i<revisedBlockList.size();i++) {
			ArrayList<Double> magList = revisedBlockList.get(i).getMagList();
			ArrayList<Double> rateInsideList = revisedBlockList.get(i).getRateInsideList();
			ArrayList<Integer> rupIndexN_List = revisedBlockList.get(i).getRupIndexN_List();
			double blockProb=randomBlockSampler.getY(i);
			ArbIncrementalMagFreqDist blockMagDist = new ArbIncrementalMagFreqDist(2.05, 8.95, 70);
			double total=0;
			for(int j=0; j<magList.size(); j++) {
				total += rateInsideList.get(j);
				if(erf.getSrcIndexForNthRup(rupIndexN_List.get(j)) == srcIndex)
					blockMagDist.addResampledMagRate(magList.get(j), rateInsideList.get(j), true);
			}
			blockMagDist.scale(blockProb/total);
			for(int j=0; j<blockMagDist.size(); j++) magDist.add(blockMagDist.getX(j), blockMagDist.getY(j));
		}
		return magDist;
	}

	
	/** This returns the mag-prob-dists of the top numMFD sources that contribute to the
	 * primary aftershock sequence
	 * 
	 * @param numMFDs - the number of MFDs desired
	 * @return
	 */
	public ArrayList<ArbIncrementalMagFreqDist> getMagProbDistForAllSources(int numMFDs) {
		HashMap<String,ArbIncrementalMagFreqDist> mfd_map = new HashMap<String,ArbIncrementalMagFreqDist>();
		for(int i=0; i<revisedBlockList.size();i++) {
			double totRateInsideBlock = revisedBlockList.get(i).getTotalRateInside();
			ArrayList<Double> magList = revisedBlockList.get(i).getMagList();
			ArrayList<Double> rateInsideList = revisedBlockList.get(i).getRateInsideList();
			ArrayList<Integer> rupIndexN_List = revisedBlockList.get(i).getRupIndexN_List();
			double blockProb=randomBlockSampler.getY(i);
			ArbIncrementalMagFreqDist  mfd;
			for(int j=0; j<magList.size(); j++) {
				int srcIndex =erf.getSrcIndexForNthRup(rupIndexN_List.get(j));
				String srcName = this.erf.getSource(srcIndex).getName();
				if(mfd_map.containsKey(srcName))
					mfd = mfd_map.get(srcName);
				else {
					mfd = new  ArbIncrementalMagFreqDist(2.05, 8.95, 70);
					mfd_map.put(srcName, mfd);
				}
				mfd.addResampledMagRate(magList.get(j), blockProb*rateInsideList.get(j)/totRateInsideBlock, true);
			}
		}
		//Check results
		double totalRate=0;
		for(String sName: mfd_map.keySet())
			totalRate += mfd_map.get(sName).getTotalIncrRate();
		System.out.println("getMagProbDistForAllSources TEST - This should be 1.0: "+totalRate);

		
		ArrayList<ArbIncrementalMagFreqDist> mfdList = new ArrayList<ArbIncrementalMagFreqDist>();
		while(mfdList.size()<numMFDs) {
			double max = 0;
			String targetName = null;
			for(String sName: mfd_map.keySet()) {
				double totRate = mfd_map.get(sName).getTotalIncrRate();
				if(totRate>max) {
					max = totRate;
					targetName=sName;
				}
			}
			ArbIncrementalMagFreqDist mfd = mfd_map.get(targetName);
			mfd.setName(targetName);
			mfd.setInfo(" ");
			mfdList.add(mfd);
			mfd_map.remove(targetName);
		}
		
		return mfdList;
	}

	
	/**
	 * This samples a random primary aftershock
	 * @return
	 */
	public ETAS_EqkRupture samplePrimaryAftershock(long originTime) {
		int blockIndex = sampleRandomBlockIndex();
		EqksInGeoBlock block = revisedBlockList.get(blockIndex);
		int nthRup = block.getRandomRuptureIndexN();  // this returns the nth index of the randomly chosen rupture
		ProbEqkRupture rup = (ProbEqkRupture)erf.getNthRupture(nthRup);
/*
		if(rup.getMag()>7.8) {
			System.out.println("BIG MAG:\t"+rup.getMag()+"\t"+blockIndex+"\tiSrc="+indices[0]+"\tiRup="+indices[1]);
			block.writeResults();			
		}
*/		
		ETAS_EqkRupture aftershock = new ETAS_EqkRupture(rup);
		aftershock.setNthERF_Index(nthRup);
		// assign a hypocenter
		block.setRandomHypocenterLoc(aftershock);
		aftershock.setOriginTime(originTime);
		return aftershock;
	}
	
	
	/**
	 * This samples a random block
	 * @return
	 */
	public EqksInGeoBlock sampleRandomBlock() {
		int blockIndex = sampleRandomBlockIndex();
		return revisedBlockList.get(blockIndex);
	}

	
	/**
	 * This returns the revised block list
	 * @return
	 */
	public ArrayList<EqksInGeoBlock> getRevisedBlockList () { return revisedBlockList; }
	
	
	public double getDistForRevisedBlock(int blockIndex) {
		return revisedBlockDistances.get(blockIndex);
	}
	
	
	/**
	 * This returns the random block sampler (where the relative rate of each block is on the y-axis, 
	 * and the block index is both the function index and x-axis value)
	 * @return
	 */
	public IntegerPDF_FunctionSampler getRandomBlockSampler() {return randomBlockSampler;}
	
	
	/**
	 * Check and write results
	 */
	public void writeRelBlockProbToFile() {
		
		System.out.println("starting random block sampling test");
		
		/**/
		EvenlyDiscretizedFunc testFunc = new EvenlyDiscretizedFunc(0.0, numBlocks, 1.0);
		int numSamples=100000000;
		int thresh = 10000000, threshIncr=10000000;
		for(int i=0;i<numSamples;i++) {
			testFunc.add(randomBlockSampler.getRandomInt(),1.0);
			if(i==thresh){
				System.out.println("\t"+((float)thresh/(float)numSamples));
				thresh+=threshIncr;
			}
		}
		for(int i=0;i<testFunc.size();i++) testFunc.set(i,testFunc.getY(i)/numSamples);
		
		// make orig dist for plotting comparison
		EvenlyDiscretizedFunc origFunc = new EvenlyDiscretizedFunc(0.0, numBlocks, 1.0);
		for(int i=0;i<numBlocks;i++)  origFunc.set(i,relBlockProb[i]);
		
		try{
			FileWriter fw1 = new FileWriter("/Users/field/workspace/OpenSHA/dev/scratch/ned/ETAS_ERF/relBlockProbs.txt");
			fw1.write("index\tlat\tlon\tdepth\trelProb\trandProbs\tvol\tblockDist\n");
			for(int i=0;i<numBlocks;i++) {
				Location loc = revisedBlockList.get(i).getBlockCenterLoc();
				double vol = revisedBlockList.get(i).getBlockVolume(); // normalize rates by vol for adaptive block spacing
				fw1.write(i+"\t"+(float)loc.getLatitude()+"\t"+(float)loc.getLongitude()+"\t"+(float)loc.getDepth()+"\t"+(float)relBlockProb[i]+"\t"+
						testFunc.getY(i)+"\t"+vol+"\t"+revisedBlockDistances.get(i)+"\n");
			}
			fw1.close();
		}catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	/**
	 * 
	 * @param label - plot label
	 * @param local - whether GMT map is made locally or on server
	 * @param dirName
	 * @return
	 */
	public String plotBlockProbMap(String label, boolean local, String dirName) {
		
		GMT_MapGenerator mapGen = new GMT_MapGenerator();
		mapGen.setParameter(GMT_MapGenerator.GMT_SMOOTHING_PARAM_NAME, false);
		mapGen.setParameter(GMT_MapGenerator.TOPO_RESOLUTION_PARAM_NAME, GMT_MapGenerator.TOPO_RESOLUTION_NONE);
		mapGen.setParameter(GMT_MapGenerator.MIN_LAT_PARAM_NAME,31.5);		// -R-125.4/-113.0/31.5/43.0
		mapGen.setParameter(GMT_MapGenerator.MAX_LAT_PARAM_NAME,43.0);
		mapGen.setParameter(GMT_MapGenerator.MIN_LON_PARAM_NAME,-125.4);
		mapGen.setParameter(GMT_MapGenerator.MAX_LON_PARAM_NAME,-113.0);
		mapGen.setParameter(GMT_MapGenerator.LOG_PLOT_NAME,true);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MODE_NAME,GMT_MapGenerator.COLOR_SCALE_MODE_MANUALLY);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME,-6.);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME,-1.);


		CaliforniaRegions.RELM_GRIDDED griddedRegion = new CaliforniaRegions.RELM_GRIDDED();
		GriddedGeoDataSet xyzDataSet = new GriddedGeoDataSet(griddedRegion, true);
		
		// initialize values to zero
		for(int i=0; i<xyzDataSet.size();i++) xyzDataSet.set(i, 0);
		
		for(int i=0; i<revisedBlockList.size();i++) {
			EqksInGeoBlock block = revisedBlockList.get(i);
			int locIndex = griddedRegion.indexForLocation(block.getBlockCenterLoc());
			double totProbInCell = relBlockProb[i]+xyzDataSet.get(locIndex);
			xyzDataSet.set(locIndex, totProbInCell);
		}
//		System.out.println("Min & Max Z: "+xyzDataSet.getMinZ()+"\t"+xyzDataSet.getMaxZ());
		String metadata = "no metadata";
		
		try {
			String name;
			if(local)
				name = mapGen.makeMapLocally(xyzDataSet, "Prob from "+label, metadata, dirName);
			else {
				name = mapGen.makeMapUsingServlet(xyzDataSet, "Prob from "+label, metadata, dirName);
				metadata += GMT_MapGuiBean.getClickHereHTML(mapGen.getGMTFilesWebAddress());
				ImageViewerWindow imgView = new ImageViewerWindow(name,metadata, true);				
			}

//			System.out.println("GMT Plot Filename: "+name);
		} catch (GMT_MapException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (RuntimeException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return "For Block Prob Map: "+mapGen.getGMTFilesWebAddress()+" (deleted at midnight)";
	}
	
	
	/**
	 * This returns a function with a simulated PDF of the distance decay (from 10000000 samples) plus
	 * the expected distance decay ((dist+minDist)^-distDecay) for the given deltaDist (where the latter, 
	 * to close approximation, takes into account the finite deltaDist effects).
	 * @param deltaDist
	 * @return
	 */
	public ArrayList<EvenlyDiscretizedFunc> getDistDecayTestFuncs(double deltaDist) {
		ArrayList<EvenlyDiscretizedFunc> funcs = new ArrayList<EvenlyDiscretizedFunc>();
		int num = (int)(2200.0/deltaDist)+1;

		EvenlyDiscretizedFunc distHist = new EvenlyDiscretizedFunc(deltaDist/2, num, deltaDist);
		distHist.setTolerance(5);
		EvenlyDiscretizedFunc target = new EvenlyDiscretizedFunc(0.5, 2000, 1.0);
		EvenlyDiscretizedFunc targetHist = new EvenlyDiscretizedFunc(deltaDist/2, num, deltaDist);
		targetHist.setTolerance(5);
		
		double sum=0;
		for(int i=0;i<10000000;i++) {
			double dist = revisedBlockDistances.get(sampleRandomBlockIndex());
			distHist.add(dist,1.0);
			sum+=1;
		}
		for(int i=0;i<distHist.size();i++) distHist.set(i,distHist.getY(i)/sum);
		
		for(int i=0; i<target.size();i++) target.set(i,Math.pow(target.getX(i)+minDist, -distDecay));
		double sum2 = target.calcSumOfY_Vals();
		for(int i=0; i<target.size();i++) target.set(i,target.getY(i)/sum2);
		
		for(int i=0; i<target.size();i++) targetHist.add(target.getX(i), target.getY(i));
		
		distHist.setName("Simulated Distance Decay for Primary Aftershocks");
		targetHist.setName("Target Distance Decay for Primary Aftershocks");
		funcs.add(distHist);
		funcs.add(targetHist);
		
		return funcs;
	}
	
	/**
	 * This compares a simulated distance decay histogram with the expected histogram
	 * @param plotTitle
	 * @param pdfFileNameAndPath - give null is you don't want to save the file
	 */
	public void plotDistDecayTestFuncs(String plotTitle, String pdfFileNameAndPath, double deltaDist) {
		
		ArrayList funcs = getDistDecayTestFuncs(deltaDist);
		GraphWindow graph = new GraphWindow(funcs, plotTitle); 
		graph.setAxisRange(1, 1200, 1e-6, 1);
		graph.setYLog(true);
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 2f, Color.BLUE));
		graph.setPlotChars(plotChars);
		graph.setX_AxisLabel("Distance (km)");
		graph.setY_AxisLabel("Probability");
		
		if(pdfFileNameAndPath != null) {
			try {
				graph.saveAsPDF(pdfFileNameAndPath);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
	}
	
	
			
	
	public double getMinDist() {
		return minDist;
	}
	
	public double getDistDecay() {
		return  distDecay;
	}


	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		System.out.println(ETAS_Utils.getDecayFractionInsideDistance(2,0.3, 100));
//		double r = Math.pow(25.6*0.75/Math.PI, 0.33333);
//		System.out.println(r);
	}

}
