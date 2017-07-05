package scratch.ned.ETAS_ERF.sandbox;

import java.util.ArrayList;
import java.util.Hashtable;

import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.mapping.gmt.GMT_MapGenerator;
import org.opensha.commons.mapping.gmt.gui.GMT_MapGuiBean;
import org.opensha.commons.mapping.gmt.gui.ImageViewerWindow;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.calc.ERF_Calculator;
import org.opensha.sha.gui.infoTools.CalcProgressBar;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.sha.magdist.SummedMagFreqDist;

import scratch.UCERF3.erf.FaultSystemSolutionPoissonERF;
import scratch.UCERF3.erf.ETAS.IntegerPDF_FunctionSampler;
import scratch.UCERF3.erf.UCERF2_Mapped.UCERF2_FM2pt1_FaultSysSolTimeDepERF;

public class ERF_RatesAtPointsInSpace {
	
	// these define the points in space
	int numDepths, numRegLocs, numPoints;
	double maxDepth, depthDiscr;
	GriddedRegion region;
	
	// this is for each point in space
	double[] latForPoint, lonForPoint, depthForPoint;
	
	// this will hold the rate of each ERF source
	double sourceRates[];
	
	// this stores the rates of erf ruptures that go unassigned (outside the region here)
	double rateUnassigned;
	double totRate;
	
	ArrayList<double[]> fractionSrcAtPointList;
	ArrayList<int[]> srcAtPointList;
	
	IntegerPDF_FunctionSampler pointSampler;
	
	

	/**
	 * TO DO
	 * 
	 * @param griddedRegion
	 * @param erf
	 * @param maxDepth
	 * @param depthDiscr
	 * @param pointSrcDiscr - the grid spacing of off-fault/background events
	 */
	public ERF_RatesAtPointsInSpace(GriddedRegion griddedRegion, FaultSystemSolutionPoissonERF erf, double sourceRates[],
			double maxDepth, double depthDiscr, double pointSrcDiscr, String oututFileNameWithPath) {
		
		this.maxDepth=maxDepth;
		this.depthDiscr=depthDiscr;
		numDepths = (int)Math.round(maxDepth/depthDiscr);
		
		this.region = griddedRegion;
		numRegLocs = griddedRegion.getNumLocations();
		numPoints = numRegLocs*numDepths;
		
		this.sourceRates = sourceRates;

		double regSpacing = griddedRegion.getLatSpacing();
		if(griddedRegion.getLonSpacing() != regSpacing)
			throw new RuntimeException("griddedRegion.getLonSpacing() must equal griddedRegion.getLatSpacing()");
		
		latForPoint = new double[numPoints];
		lonForPoint = new double[numPoints];
		depthForPoint = new double[numPoints];
		for(int i=0;i<this.numPoints;i++) {
			int iDep = (int)Math.floor((double)i/(double)numRegLocs);	// depth index
			int iReg = i - iDep*numRegLocs;								// region index
			Location loc = region.getLocation(iReg);
			latForPoint[i] = loc.getLatitude();
			lonForPoint[i] = loc.getLongitude();
			depthForPoint[i] = this.getDepth(iDep);
		}

		System.out.println("Making sourcesAtPointList & fractionsAtPointList");
		ArrayList<ArrayList<Integer>> sourcesAtPointList = new ArrayList<ArrayList<Integer>>();
		ArrayList<ArrayList<Double>> fractionsAtPointList = new ArrayList<ArrayList<Double>>();
		for(int i=0; i<numPoints;i++) {
			sourcesAtPointList.add(new ArrayList<Integer>());
			fractionsAtPointList.add(new ArrayList<Double>());
		}

		int numPtSrcSubPts = (int)Math.round(pointSrcDiscr/regSpacing);
		double extra;	// this is needed to get the right intervals withing each point source
		if (numPtSrcSubPts % 2 == 0) {	// if even
			extra=0;
		}
		else							// if odd
			extra = regSpacing/2.0;
	

		rateUnassigned=0;

		int totNumSrc = erf.getNumSources();
		if(totNumSrc != sourceRates.length)
			throw new RuntimeException("Problem with number of sources");
		
		int numFltSystRups = erf.getNumFaultSystemSources();

		
		CalcProgressBar progressBar = new CalcProgressBar("Events to process", "junk");
		progressBar.displayProgressBar();
		progressBar.showProgress(true);
	
		System.out.println("Starting loop to populate fractionSrcAtPointList & srcAtPointList");
		for(int s=0;s<totNumSrc;s++) {
			ProbEqkSource src = erf.getSource(s);
			progressBar.updateProgress(s, totNumSrc);

			// If it's not a point sources:
			if(s<numFltSystRups) {
				Hashtable<Integer,Double> fractAtPointTable = new Hashtable<Integer,Double>(); // int is ptIndex and double is fraction there
				LocationList locsOnRupSurf = src.getRupture(0).getRuptureSurface().getEvenlyDiscritizedListOfLocsOnSurface();
				int numLocs = locsOnRupSurf.size();
				for(Location loc: locsOnRupSurf) {
					int regIndex = griddedRegion.indexForLocation(loc);
					int depIndex = getDepthIndex(loc.getDepth());
					if(regIndex != -1) {
						int ptIndex = depIndex*numRegLocs+regIndex;
						if(fractAtPointTable.containsKey(ptIndex)) {
							double newFrac = fractAtPointTable.get(ptIndex) + 1.0/(double)numLocs;
							fractAtPointTable.put(ptIndex,newFrac);
						}
						else {
							fractAtPointTable.put(ptIndex,1.0/(double)numLocs);
						}
					}
					else {
						rateUnassigned += sourceRates[s]/numLocs;
					}
				}	
				// now assign this hashTable
				for(Integer ptIndex : fractAtPointTable.keySet()) {
					double fract = fractAtPointTable.get(ptIndex);
					sourcesAtPointList.get(ptIndex).add(s);
					fractionsAtPointList.get(ptIndex).add(fract);
				}
			}
			else {	// It's a point source
				for(ProbEqkRupture rup: src)
					if(!rup.getRuptureSurface().isPointSurface())	// make sure they're all point surfaces
						throw new RuntimeException("All ruptures for source must have point surfaces here");

				Location centerLoc = src.getRupture(0).getRuptureSurface().getFirstLocOnUpperEdge();
				double numPts = numPtSrcSubPts*numPtSrcSubPts*numDepths;
				double ptRate = sourceRates[s]/numPts;
				double ptFrac = 1.0/numPts;
				// distribution this among the locations within the space represented by the point source
				for(int iLat=0; iLat<numPtSrcSubPts;iLat++) {
					double lat = centerLoc.getLatitude()-pointSrcDiscr/2 + iLat*regSpacing+extra;
					for(int iLon=0; iLon<numPtSrcSubPts;iLon++) {
						double lon = centerLoc.getLongitude()-pointSrcDiscr/2 + iLon*regSpacing+extra;
						int regIndex = griddedRegion.indexForLocation(new Location(lat,lon));
						if(regIndex != -1){
							for(int iDep =0; iDep<numDepths; iDep++) {
								int ptIndex = iDep*numRegLocs+regIndex;
								sourcesAtPointList.get(ptIndex).add(s);
								fractionsAtPointList.get(ptIndex).add(ptFrac);
							}
						}
						else {
							rateUnassigned += ptRate*numDepths;
//							System.out.println("1\t"+centerLoc.getLatitude()+"\t"+centerLoc.getLongitude()+"\t"+centerLoc.getDepth());
						}
					}
				}
			}
		}
		progressBar.showProgress(false);
		System.out.println("rateUnassigned="+rateUnassigned);
		
		System.out.println("Converting list types");
		fractionSrcAtPointList = new ArrayList<double[]>();
		srcAtPointList = new ArrayList<int[]> ();
		for(int i=0;i<numPoints;i++) {
			ArrayList<Integer> sourceList = sourcesAtPointList.get(i);
			ArrayList<Double> fractList = fractionsAtPointList.get(i);
			int[] sourceArray = new int[sourceList.size()];
			double[] fractArray = new double[fractList.size()];
			for(int j=0;j<sourceArray.length;j++) {
				sourceArray[j] = sourceList.get(j);
				fractArray[j] = fractList.get(j);
			}
//if(i==500 || i==1000000)
//	for(int z=0;z<sourceArray.length;z++)
//		System.out.println("\t"+sourceArray[z]+"\n"+fractArray[z]+"\t"+(sourceArray[z]>erf.getNumFaultSystemSources()));
			srcAtPointList.add(sourceArray);
			fractionSrcAtPointList.add(fractArray);
		}
		System.out.println("Done converting list types");
		
		// write results to file
//		if(oututFileNameWithPath != null)
//			writeEqksAtPointArrayToFile(oututFileNameWithPath);
		
	}
	
	
	public void declareRateChange() {
		pointSampler = null;
	}
	
	
	
	/**
	 * This method
	 */
	public IntegerPDF_FunctionSampler getPointSampler() {
		if(pointSampler == null) {
			pointSampler = new IntegerPDF_FunctionSampler(numPoints);
			for(int i=0;i<numPoints;i++) {
				int[] sources = srcAtPointList.get(i);
				double[] fract = fractionSrcAtPointList.get(i);
				double totRate=0;
				for(int j=0; j<sources.length;j++) {
					totRate += sourceRates[sources[j]]*fract[j];
				}
				pointSampler.set(i,totRate);
			}
		}
		return pointSampler;
	}
	
	
	public IntegerPDF_FunctionSampler getPointSamplerWithDistDecay(EqkRupture mainshock, ETAS_LocationWeightCalculator etasLocWtCalc) {
		getPointSampler();	// this makes sure it updated
		IntegerPDF_FunctionSampler sampler = new IntegerPDF_FunctionSampler(numDepths*numRegLocs);
		LocationList locList = mainshock.getRuptureSurface().getEvenlyDiscritizedListOfLocsOnSurface();
		for(int index=0; index<this.numPoints; index++) {
			double ptWt = 0;
			// loop over locations on rupture surface
			for(Location loc : locList) {
				double relLat = Math.abs(loc.getLatitude()-latForPoint[index]);
				double relLon = Math.abs(loc.getLongitude()-lonForPoint[index]);
				double relDep = Math.abs(loc.getDepth()-depthForPoint[index]);
				ptWt += etasLocWtCalc.getProbAtPoint(relLat, relLon, relDep, loc.getDepth());
			}
			sampler.set(index,ptWt*pointSampler.getY(index));

		}
		return sampler;
	}

	/**
	 * This sampler ignores the long-term rates
	 * @param mainshock
	 * @param etasLocWtCalc
	 * @return
	 */
	public IntegerPDF_FunctionSampler getPointSamplerWithOnlyDistDecay(EqkRupture mainshock, ETAS_LocationWeightCalculator etasLocWtCalc) {
		IntegerPDF_FunctionSampler sampler = new IntegerPDF_FunctionSampler(numDepths*numRegLocs);
		LocationList locList = mainshock.getRuptureSurface().getEvenlyDiscritizedListOfLocsOnSurface();
		for(int index=0; index<this.numPoints; index++) {
			double ptWt = 0;
			// loop over locations on rupture surface
			for(Location loc : locList) {
				double relLat = Math.abs(loc.getLatitude()-latForPoint[index]);
				double relLon = Math.abs(loc.getLongitude()-lonForPoint[index]);
				double relDep = Math.abs(loc.getDepth()-depthForPoint[index]);
				ptWt += etasLocWtCalc.getProbAtPoint(relLat, relLon, relDep, loc.getDepth());
			}
			sampler.set(index,ptWt);

		}
		return sampler;
	}
	
	public int getRandomSourceAtPoint(int ptIndex) {
		int[] sources = srcAtPointList.get(ptIndex);
		double[] fracts = this.fractionSrcAtPointList.get(ptIndex);
		IntegerPDF_FunctionSampler sampler = new IntegerPDF_FunctionSampler(sources.length);
		for(int s=0; s<sources.length;s++) 
			sampler.set(s,sourceRates[sources[s]]*fracts[s]);
		return sources[sampler.getRandomInt()];
	}

	
	/**
	 * Region index is first element, and depth index is second
	 * @param index
	 * @return
	 */
	private int[] getRegAndDepIndicesForSamplerIndex(int index) {
		
		int[] indices = new int[2];
		indices[1] = (int)Math.floor((double)index/(double)numRegLocs);	// depth index
		if(indices[1] >= this.numDepths )
			System.out.println("PROBLEM: "+index+"\t"+numRegLocs+"\t"+indices[1]+"\t"+numDepths);
		indices[0] = index - indices[1]*numRegLocs;						// region index
		return indices;
	}
	
	public Location getLocationForSamplerIndex(int index) {
		int[] regAndDepIndex = getRegAndDepIndicesForSamplerIndex(index);
		Location regLoc = region.getLocation(regAndDepIndex[0]);
		return new Location(regLoc.getLatitude(),regLoc.getLongitude(),getDepth(regAndDepIndex[1]));
	}
	
	public int getSamplerIndexForLocation(Location loc) {
		int iReg = region.indexForLocation(loc);
		int iDep = getDepthIndex(loc.getDepth());
		return getSamplerIndexForRegAndDepIndices(iReg,iDep);
	}

	private int getSamplerIndexForRegAndDepIndices(int iReg,int iDep) {
		return iDep*numRegLocs+iReg;
	}


	public void testRates(FaultSystemSolutionPoissonERF erf) {
		
		System.out.println("Testing total rate");
		getPointSampler();
		
		totRate=this.pointSampler.calcSumOfY_Vals();
		totRate+=rateUnassigned;
		
		double testRate2=0;
		double duration = erf.getTimeSpan().getDuration();
		for(int s=0;s<erf.getNumSources();s++) {
			ProbEqkSource src = erf.getSource(s);
			int numRups = src.getNumRuptures();
			for(int r=0; r<numRups;r++) {
				testRate2 += src.getRupture(r).getMeanAnnualRate(duration);
			}
		}
		System.out.println("\ttotRateTest="+(float)totRate+" should equal Rate2="+(float)testRate2+";\tratio="+(float)(totRate/testRate2));
	}
	
	
	public void testMagFreqDist(FaultSystemSolutionPoissonERF erf) {
		
		System.out.println("Running testMagFreqDist()");
		SummedMagFreqDist magDist = new SummedMagFreqDist(2.05, 8.95, 70);
//		getPointSampler();	// make sure it exisits
		double duration = erf.getTimeSpan().getDuration();
		for(int i=0; i<numPoints;i++) {
			int[] sources = srcAtPointList.get(i);
			double[] fracts = fractionSrcAtPointList.get(i);
			for(int s=0; s<sources.length;s++) {
				SummedMagFreqDist mfd = ERF_Calculator.getTotalMFD_ForSource(erf.getSource(sources[s]), duration, 2.05, 8.95, 70, true);
				mfd.scale(fracts[s]);
				magDist.addIncrementalMagFreqDist(mfd);
				if(s>erf.getNumFaultSystemSources()) {
					System.out.println("source "+s+"\n"+mfd);
					System.exit(0);
				}
			}
		}
		magDist.setName("MFD from EqksAtPoint list");
		ArrayList<EvenlyDiscretizedFunc> magDistList = new ArrayList<EvenlyDiscretizedFunc>();
		magDistList.add(magDist);
		magDistList.add(magDist.getCumRateDistWithOffset());
		
		SummedMagFreqDist erfMFD = ERF_Calculator.getTotalMFD_ForERF(erf, 2.05, 8.95, 70, true);
		erfMFD.setName("MFD from ERF");
		magDistList.add(erfMFD);
		magDistList.add(erfMFD.getCumRateDistWithOffset());

		// Plot these MFDs
		GraphWindow magDistsGraph = new GraphWindow(magDistList, "Mag-Freq Distributions"); 
		magDistsGraph.setX_AxisLabel("Mag");
		magDistsGraph.setY_AxisLabel("Rate");
		magDistsGraph.setY_AxisRange(1e-6, magDistsGraph.getY_AxisRange().getUpperBound());
		magDistsGraph.setYLog(true);

	}
	
	
	private int getDepthIndex(double depth) {
		return (int)Math.round((depth-depthDiscr/2.0)/depthDiscr);
	}
	
	private double getDepth(int depthIndex) {
		return (double)depthIndex*depthDiscr + depthDiscr/2;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		System.out.println("Instantiating ERF");
		UCERF2_FM2pt1_FaultSysSolTimeDepERF erf = new UCERF2_FM2pt1_FaultSysSolTimeDepERF();
		erf.updateForecast();
		
		double sourceRates[] = new double[erf.getNumSources()];
		double duration = erf.getTimeSpan().getDuration();
		for(int s=0;s<erf.getNumSources();s++)
			sourceRates[s] = erf.getSource(s).computeTotalEquivMeanAnnualRate(duration);

//		CaliforniaRegions.RELM_GRIDDED griddedRegion = new CaliforniaRegions.RELM_GRIDDED();
		
		System.out.println("Instantiating Region");
		GriddedRegion gridedRegion = new GriddedRegion(new CaliforniaRegions.RELM_TESTING(), 0.02, GriddedRegion.ANCHOR_0_0);
		
		long startTime = System.currentTimeMillis();
		System.out.println("Instantiating ERF_RatesAtPointsInSpace");
		
		String testFileName = "/Users/field/workspace/OpenSHA/dev/scratch/ned/ETAS_ERF/testBinaryFile";

//		ERF_RatesAtPointsInSpace erf_RatesAtPointsInSpace = new ERF_RatesAtPointsInSpace(gridedRegion, erf, sourceRates, 24d,2d,0.1,testFileName);
		System.out.println("Instantiating took "+(System.currentTimeMillis()-startTime)/1000+" sec");

		// TESTS
//		erf_RatesAtPointsInSpace.testRates(erf);
//		erf_RatesAtPointsInSpace.testMagFreqDist(erf);
//		erf_RatesAtPointsInSpace.plotRatesMap("test", true, "/Users/field/workspace/OpenSHA/dev/scratch/ned/ETAS_ERF/mapTest");
//		erf_RatesAtPointsInSpace.plotOrigERF_RatesMap("orig test", true, "/Users/field/workspace/OpenSHA/dev/scratch/ned/ETAS_ERF/mapOrigTest", erf);
//		erf_RatesAtPointsInSpace.plotRandomSampleRatesMap("random test", true, "/Users/field/workspace/OpenSHA/dev/scratch/ned/ETAS_ERF/mapRandomTest", erf,10000);

		

	}
	
	
	
	/**
	 * This plots the rates in space for the RELM region (rates are summed inside each spatial bin).
	 * 
	 * Compare this to what produced by the plotOrigERF_RatesMap(*) method (should be same)
	 * 
	 * @param label - plot label
	 * @param local - whether GMT map is made locally or on server
	 * @param dirName
	 * @return
	 */
	public String plotRatesMap(IntegerPDF_FunctionSampler pointSamplerForMap, String label, boolean local, String dirName) {
		
		GMT_MapGenerator mapGen = new GMT_MapGenerator();
		mapGen.setParameter(GMT_MapGenerator.GMT_SMOOTHING_PARAM_NAME, false);
		mapGen.setParameter(GMT_MapGenerator.TOPO_RESOLUTION_PARAM_NAME, GMT_MapGenerator.TOPO_RESOLUTION_NONE);
		mapGen.setParameter(GMT_MapGenerator.MIN_LAT_PARAM_NAME,31.5);		// -R-125.4/-113.0/31.5/43.0
		mapGen.setParameter(GMT_MapGenerator.MAX_LAT_PARAM_NAME,43.0);
		mapGen.setParameter(GMT_MapGenerator.MIN_LON_PARAM_NAME,-125.4);
		mapGen.setParameter(GMT_MapGenerator.MAX_LON_PARAM_NAME,-113.0);
		mapGen.setParameter(GMT_MapGenerator.LOG_PLOT_NAME,true);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MODE_NAME,GMT_MapGenerator.COLOR_SCALE_MODE_MANUALLY);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME,-3.5);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME,1.5);


		CaliforniaRegions.RELM_GRIDDED mapGriddedRegion = new CaliforniaRegions.RELM_GRIDDED();
		GriddedGeoDataSet xyzDataSet = new GriddedGeoDataSet(mapGriddedRegion, true);
		
		// initialize values to zero
		for(int i=0; i<xyzDataSet.size();i++) xyzDataSet.set(i, 0);
		
		getPointSampler();
		for(int i=0;i<numPoints;i++) {
			Location loc = getLocationForSamplerIndex(i);
			int mapLocIndex = mapGriddedRegion.indexForLocation(loc);
			if(mapLocIndex>=0) {
				double oldRate = xyzDataSet.get(mapLocIndex);
				xyzDataSet.set(mapLocIndex, pointSamplerForMap.getY(i)+oldRate);					
			}

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
	 * 
	 * @param label - plot label
	 * @param local - whether GMT map is made locally or on server
	 * @param dirName
	 * @return
	 */
	public String plotOrigERF_RatesMap(String label, boolean local, String dirName, FaultSystemSolutionPoissonERF erf) {
		
		GMT_MapGenerator mapGen = new GMT_MapGenerator();
		mapGen.setParameter(GMT_MapGenerator.GMT_SMOOTHING_PARAM_NAME, false);
		mapGen.setParameter(GMT_MapGenerator.TOPO_RESOLUTION_PARAM_NAME, GMT_MapGenerator.TOPO_RESOLUTION_NONE);
		mapGen.setParameter(GMT_MapGenerator.MIN_LAT_PARAM_NAME,31.5);		// -R-125.4/-113.0/31.5/43.0
		mapGen.setParameter(GMT_MapGenerator.MAX_LAT_PARAM_NAME,43.0);
		mapGen.setParameter(GMT_MapGenerator.MIN_LON_PARAM_NAME,-125.4);
		mapGen.setParameter(GMT_MapGenerator.MAX_LON_PARAM_NAME,-113.0);
		mapGen.setParameter(GMT_MapGenerator.LOG_PLOT_NAME,true);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MODE_NAME,GMT_MapGenerator.COLOR_SCALE_MODE_MANUALLY);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME,-3.5);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME,1.5);


		CaliforniaRegions.RELM_GRIDDED mapGriddedRegion = new CaliforniaRegions.RELM_GRIDDED();
		GriddedGeoDataSet xyzDataSet = new GriddedGeoDataSet(mapGriddedRegion, true);
		
		// initialize values to zero
		for(int i=0; i<xyzDataSet.size();i++) xyzDataSet.set(i, 0);
		
		double duration = erf.getTimeSpan().getDuration();
		for(ProbEqkSource src : erf) {
			for(ProbEqkRupture rup : src) {
				LocationList locList = rup.getRuptureSurface().getEvenlyDiscritizedListOfLocsOnSurface();
				double ptRate = rup.getMeanAnnualRate(duration)/locList.size();
				for(Location loc:locList) {
					int locIndex = mapGriddedRegion.indexForLocation(loc);
					if(locIndex>=0) {
						double oldRate = xyzDataSet.get(locIndex);
						xyzDataSet.set(locIndex, ptRate+oldRate);					
					}
				}
			}
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
	 * 
	 * @param label - plot label
	 * @param local - whether GMT map is made locally or on server
	 * @param dirName
	 * @return
	 */
	public String plotRandomSampleRatesMap(String label, boolean local, String dirName, FaultSystemSolutionPoissonERF erf, int numYrs) {
		
		GMT_MapGenerator mapGen = new GMT_MapGenerator();
		mapGen.setParameter(GMT_MapGenerator.GMT_SMOOTHING_PARAM_NAME, false);
		mapGen.setParameter(GMT_MapGenerator.TOPO_RESOLUTION_PARAM_NAME, GMT_MapGenerator.TOPO_RESOLUTION_NONE);
		mapGen.setParameter(GMT_MapGenerator.MIN_LAT_PARAM_NAME,31.5);		// -R-125.4/-113.0/31.5/43.0
		mapGen.setParameter(GMT_MapGenerator.MAX_LAT_PARAM_NAME,43.0);
		mapGen.setParameter(GMT_MapGenerator.MIN_LON_PARAM_NAME,-125.4);
		mapGen.setParameter(GMT_MapGenerator.MAX_LON_PARAM_NAME,-113.0);
		mapGen.setParameter(GMT_MapGenerator.LOG_PLOT_NAME,true);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MODE_NAME,GMT_MapGenerator.COLOR_SCALE_MODE_MANUALLY);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME,-3.5);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME,1.5);


		CaliforniaRegions.RELM_GRIDDED mapGriddedRegion = new CaliforniaRegions.RELM_GRIDDED();
		GriddedGeoDataSet xyzDataSet = new GriddedGeoDataSet(mapGriddedRegion, true);
		
		// initialize values to zero
		for(int i=0; i<xyzDataSet.size();i++) xyzDataSet.set(i, 0);
		
		// get numYrs yrs worth of samples
		int numSamples = numYrs*(int)totRate;
		System.out.println("num random samples for map test = "+numSamples);
		// do this to make sure it exists
		getPointSampler();
		
		for(int i=0;i<numSamples;i++) {
			int indexFromSampler = pointSampler.getRandomInt();
			int[] regAndDepIndex = getRegAndDepIndicesForSamplerIndex(indexFromSampler);
			int indexForMap = mapGriddedRegion.indexForLocation(region.locationForIndex(regAndDepIndex[0]));	// ignoring depth
			double oldNum = xyzDataSet.get(indexForMap)*numYrs;
			xyzDataSet.set(indexForMap, (1.0+oldNum)/(double)numYrs);
			
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


}
